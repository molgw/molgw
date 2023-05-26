!=========================================================================
! This file is part of MOLGW.
! Author: Ivan Maliyov
!
! This file contains
! the main program to transform TDDFT time data into frequency spectra
!
!=========================================================================
#include "molgw.h"
program spectrum
  use m_definitions
  use m_warning
  implicit none

  integer                    :: ntau,ntau_read,itau,idir,iomega,nomega
  integer                    :: file_transforms, file_dipolar_spectra
  integer                    :: file_dipole_time, file_excitation, io
  integer                    :: file_dipole_damped
  integer                    :: irow,nrow,num_fields
  real(dp)                   :: time_cur, time_min, time_sim
  real(dp)                   :: damp_factor
  real(dp)                   :: omega_factor
  real(dp)                   :: div_factor
  real(dp)                   :: real_dipole(3)
  real(dp)                   :: real_excit_field_dir
  real(dp),allocatable       :: m_times(:)
  complex(dp)                :: average
  complex(dp),allocatable    :: dipole_time_ref(:,:)
  complex(dp),allocatable    :: dipole_time_damped(:,:)
  complex(dp),allocatable    :: trans_m_excit_field_dir(:)
  complex(dp),allocatable    :: m_excit_field_dir(:)
  complex(dp),allocatable    :: trans_dipole_time(:,:)
  logical                    :: file_exists
  logical                    :: output_damped_
  logical                    :: output_transforms_
  character(len=100)         :: name_dipole_time, name_excitation
  character(len=100)         :: name_dipolar_spectra, name_transforms='transforms.dat'
  character(len=100)         :: name_dipole_damped='dipole_damped.dat'
  character(len=500)         :: cur_string
  type(C_PTR)                :: plan
  !=====
#if defined(HAVE_FFTW)
#include "fftw3.f03"


  call get_spectrum_arguments( damp_factor,name_dipole_time,name_excitation,name_dipolar_spectra,output_damped_,output_transforms_)

  inquire(file=name_dipole_time,exist=file_exists)
  if(.NOT. file_exists) then
    call die('No dipole_time file found')
  endif

  inquire(file=name_excitation,exist=file_exists)
  if(.NOT. file_exists) then
    call die('No excitation file found')
  endif

  inquire(file=name_dipolar_spectra,exist=file_exists)
  if(.NOT. file_exists) then
    call issue_warning('Dipolar spectra file will be rewritten')
  endif

  ! Get number of lines nrow and lines without "#" ntau in dipole_time file
  ntau = 0
  nrow = 0
  open(newunit=file_dipole_time, file = name_dipole_time)
  do
    read(file_dipole_time,'(A)',iostat=io)cur_string
    if ( io /= 0 ) exit
    nrow = nrow + 1
    cur_string = ADJUSTL(cur_string)
    if( cur_string(1:1) /= "#" ) then
      num_fields = get_number_of_elements(cur_string)
      if( num_fields /= 4 ) call die('Dipole file must contain 4 columns: time dx dy dz')
      ntau = ntau + 1
    endif
  enddo
  close(file_dipole_time)

  allocate(m_times(ntau))
  allocate(dipole_time_ref(ntau,3))
  allocate(dipole_time_damped(ntau,3))
  allocate(m_excit_field_dir(ntau))
  allocate(trans_m_excit_field_dir(ntau))
  allocate(trans_dipole_time(ntau,3))

  ! Fill the dipole_time array
  itau = 1
  open(newunit=file_dipole_time, file = name_dipole_time)
  do irow=1,nrow
    read(file_dipole_time,'(A)')cur_string
    cur_string = ADJUSTL(cur_string)
    if( cur_string(1:1) /= "#" ) then
      read(cur_string,*)m_times(itau),real_dipole(:)
      dipole_time_ref(itau,:) = real_dipole(:) / au_debye
      itau = itau + 1
    endif
  enddo
  close(file_dipole_time)

  do idir=1,3
    dipole_time_ref(:,idir) = dipole_time_ref(:,idir) - calc_average(dipole_time_ref(:,idir))
  enddo

  time_sim=m_times(ntau)

  ! Get number of lines nrow and lines without "#" ntau_read in dipole_time file
  ntau_read = 0
  nrow = 0
  open(newunit=file_excitation, file = name_excitation)
  do
    read(file_excitation,'(A)',iostat=io)cur_string
    if ( io /= 0 ) exit
    nrow = nrow + 1
    cur_string = ADJUSTL(cur_string)
    if( cur_string(1:1) /= "#" ) then
      num_fields = get_number_of_elements(cur_string)
      if( num_fields /= 2 ) call die('Excitation file must contain 2 columns: time E_excit_dir')
      ntau_read = ntau_read + 1
    endif
  enddo
  close(file_excitation)

  if( ntau /= ntau_read ) then
    call die('Different number of lines in the dipole_time file and excitation file')
  endif

  ! Fill the m_excit_field array
  itau = 1
  open(newunit=file_excitation, file = name_excitation)
  do irow=1,nrow
    read(file_excitation,'(A)') cur_string
    cur_string = ADJUSTL(cur_string)
    if( cur_string(1:1) /= "#" ) then
      read(cur_string,*) time_cur,real_excit_field_dir
      if( time_cur-m_times(itau) > 1.0e-7_dp ) then
        call die('Time values in the dipole_time file and excitation file are not the same')
      endif
      m_excit_field_dir(itau) = real_excit_field_dir
      itau = itau + 1
    endif
  enddo
  close(file_excitation)

  ! Apply damping to the dipole
  if(output_damped_) then
    open(newunit=file_dipole_damped, file = name_dipole_damped)
    write(file_dipole_damped,*) "# time(au)                 Dipole_damped_x(D)" // &
                                "        Dipole_damped_y(D)        Dipole_damped_z(D)   EXP(-m_times(itau) / damp_factor )"
  endif
  do itau=1,ntau
    dipole_time_damped(itau,:)=dipole_time_ref(itau,:)*EXP(-m_times(itau) / damp_factor )
    if(output_damped_) then
      write(file_dipole_damped,*) m_times(itau), REAL(dipole_time_damped(itau,:),dp)*au_debye, EXP(-m_times(itau) / damp_factor )
    endif
  enddo
  if(output_damped_) close(file_dipole_damped)

  ! Fourier transform of dipole_time_damped
  do idir=1,3
    plan = fftw_plan_dft_1d(ntau,dipole_time_damped(1:ntau,idir),trans_dipole_time(1:ntau,idir),FFTW_FORWARD,FFTW_ESTIMATE)
    call fftw_execute_dft(plan,dipole_time_damped(1:ntau,idir),trans_dipole_time(1:ntau,idir))
    call fftw_destroy_plan(plan)
  enddo

  trans_dipole_time =  trans_dipole_time / ntau

  ! Fourier transform of m_excit_field_dir
  plan = fftw_plan_dft_1d(ntau,m_excit_field_dir(1:ntau),trans_m_excit_field_dir(1:ntau),FFTW_FORWARD,FFTW_ESTIMATE)
  call fftw_execute_dft(plan,m_excit_field_dir(1:ntau),trans_m_excit_field_dir(1:ntau))
  call fftw_destroy_plan(plan)

  trans_m_excit_field_dir(:)=trans_m_excit_field_dir(:) / ntau

  ! Write absorption spectra in the dipolar_spectra file
  if(output_transforms_) then
    open(newunit=file_transforms,file=name_transforms)
    write(file_transforms,*) "# omega(eV), |E_excit_dir(omega)|, real(E_excit_dir(omega))," // &
                             " aimag(E_excit_dir(omega)), |d_x(omega)|, real(d_x(omega)), aimag(d_x(omega))"
  endif
  open(newunit=file_dipolar_spectra,file=name_dipolar_spectra)
  write(file_dipolar_spectra,*) "# omega(eV), transform(dipole_x)/|transform(E_dir)|," // &
                                " transform(dipole_y)/|transform(E_dir)|, transform(dipole_z)/|transform(E_dir)|"

  nomega = ntau
  do iomega=1,nomega/2 ! here iomega correspons to frequency
    omega_factor = 4.0_dp * pi * 2 * pi * (iomega-1) / time_sim / c_speedlight
    if(2*pi * iomega / time_sim * Ha_eV < 500.0_dp) then
      write(file_dipolar_spectra,"(x,es16.8E3)",advance='no')  2 * pi * iomega / time_sim * Ha_eV

      do idir=1,3
        if(ABS(trans_m_excit_field_dir(iomega))>1.0e-15) then
          write(file_dipolar_spectra,"(3(3x,es16.8E3))",advance='no') &
                  AIMAG((trans_dipole_time(iomega,idir))/trans_m_excit_field_dir(iomega)) * omega_factor
        else
          write(file_dipolar_spectra,"(3(17x,f5.2))",advance='no') -1.0_dp
        endif
        if(idir==3) write(file_dipolar_spectra,*)
      enddo
      if(output_transforms_) then
        write(file_transforms,*) pi * iomega / time_sim * Ha_eV, ABS(trans_m_excit_field_dir(iomega)), &
                                 REAL(trans_m_excit_field_dir(iomega),dp), &
                                 AIMAG(trans_m_excit_field_dir(iomega)), ABS(trans_dipole_time(iomega,1)), &
                                 REAL(trans_dipole_time(iomega,1),dp), AIMAG(trans_dipole_time(iomega,1))
      endif
    endif
  enddo



  close(file_transforms)
  close(file_dipolar_spectra)
#else
  call die("tddft: calculate_propagation; FFTW is not present")
#endif

contains


!=========================================================================
function calc_average(matrix) result(aver)
  complex(dp)  :: matrix(:)
  complex(dp)  :: aver
  !=====
  integer :: mdim,imat
  !=====

  mdim = SIZE( matrix, DIM=1)

  aver = (0.0_dp,0.0_dp)
  do imat=1,mdim
    aver = aver + matrix(imat)
  enddo
  aver = aver / mdim

end function calc_average

end program spectrum


!====================================================================================
subroutine get_spectrum_arguments(damp_factor,name_dipole_time,name_excitation, &
                                  name_dipolar_spectra,output_damped_,output_transforms_)
  use m_definitions
  use m_warning
  implicit none
  real(dp),          intent(inout)        :: damp_factor
  character(len=100),intent(inout)        :: name_dipole_time, name_excitation
  character(len=100),intent(inout)        :: name_dipolar_spectra
  logical,           intent(inout)        :: output_damped_
  logical,           intent(inout)        :: output_transforms_
  !=====
  integer             :: iarg
  character(len=20)   :: arg_cur
  character(len=20)   :: ch_damp_factor
  !=====

  damp_factor=            500.0_dp
  name_dipole_time=      'dipole_time.dat'
  name_excitation=       'excitation_time.dat'
  name_dipolar_spectra=  'dipolar_spectra.dat'
  output_damped_=.false.
  output_transforms_=.false.

  iarg = 1
  do while (iarg<=COMMAND_ARGUMENT_COUNT())
    call GET_COMMAND_ARGUMENT(iarg,arg_cur)
    select case (arg_cur)
    case ('-id','--input-dipole')
      call GET_COMMAND_ARGUMENT(iarg+1,name_dipole_time)
      iarg = iarg + 1
    case ('-ie','--input-excitation')
      call GET_COMMAND_ARGUMENT(iarg+1,name_excitation)
      iarg = iarg + 1
    case ('-o','--output')
      call GET_COMMAND_ARGUMENT(iarg+1,name_dipolar_spectra)
      iarg = iarg + 1
    case('-d','--damping')
      call GET_COMMAND_ARGUMENT(iarg+1,ch_damp_factor)
      read(ch_damp_factor,*)damp_factor
      iarg = iarg + 1
    case('-od','--output-damped')
      output_damped_=.true.
    case('-ot','--output-transforms')
      output_transforms_=.true.
    case('--help')
      write(stdout,'(1x,a)') "Calculation of dipolar spectrum after RT-TDDFT simulation"
      write(stdout,'(1x,a)') "-id, --input-dipole  file        indicate an input dipole time file"
      write(stdout,'(1x,a)') "                                (if different from dipole_time.dat)"
      write(stdout,'(1x,a)') "-ie, --input-excitation file    indicate an input excitation time file"
      write(stdout,'(1x,a)') "                                (if different from excitation_time.dat)"
      write(stdout,'(1x,a)') "-o , --output file               indicate an ouput spectrum file"
      write(stdout,'(1x,a)') "                                (if different from dipolar_spectra.dat)"
      write(stdout,'(1x,a,a)') "-d , --damping                   damping coefficient in au of time", &
                             " for the Fourier transform of dipolar moment"
      write(stdout,'(1x,a)') "                               (value by default if 500 au)"
      write(stdout,'(1x,a)') "-od, --output-damped            write damped dipole time dependence in to the dipole_damped.dat file"
      write(stdout,'(1x,a)') "-ot, --output-transfoms         write some Fourier transforms in file transforms.dat"
      stop
    case default
      write(stdout,*) "No instructions for: ", arg_cur
      call die("Unknown command line option")
    end select
    iarg = iarg + 1
  enddo

  write(stdout,'(1x,70("="))')
  write(stdout,'(12x,a)') 'MOLGW post processing'
  write(stdout,'(1x,70("="))')
  write(stdout,'(/,10x,a,/)') ' === Input parameters  ==='
  write(stdout,'(1x,a,2x,es16.8)')   'Dipole damping factor       :',damp_factor
  write(stdout,'(1x,a,2x,a50)')      'Input dipole time file      :',name_dipole_time
  write(stdout,'(1x,a,2x,a50)')      'Input excitaiton field file :',name_excitation
  write(stdout,'(1x,a,2x,a50)')      'Output dipolar spectra file :',name_dipolar_spectra
  write(stdout,'(1x,a,2x,l1)')       'Output damped dipole file   :',output_damped_
  write(stdout,'(1x,a,2x,l1)')       'Output Fourier transforms   :',output_transforms_
  write(stdout,*)


end subroutine get_spectrum_arguments


!=========================================================================
