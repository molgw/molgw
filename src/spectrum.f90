program spectrum
 use m_definitions
 use m_warning
 use m_tools
 use ISO_C_BINDING
 implicit none
!=====
 integer                    :: ntau, ntau_read,itau,idir,iomega,nomega
 integer                    :: file_transforms, file_dipolar_spectra 
 integer                    :: file_dipole_time, file_excitation, io
 integer                    :: file_dipole_damped
 integer                    :: irow,nrow,num_fields
 real(dp)                   :: time_cur, time_min, time_sim
 real(dp)                   :: damp_factor
 real(dp)                   :: omega_factor
 real(dp)                   :: average
 real(dp)                   :: div_factor
 real(dp)                   :: real_dipole(3)
 real(dp)                   :: real_excit_field(3)
 real(dp),allocatable       :: m_times(:)
 complex(dp),allocatable    :: dipole_time_ref(:,:)
 complex(dp),allocatable    :: dipole_time_damped(:,:)
 complex(dp),allocatable    :: trans_m_excit_field(:,:)
 complex(dp),allocatable    :: m_excit_field(:,:)
 complex(dp),allocatable    :: trans_dipole_time(:,:)
 logical                    :: file_exists
 logical                    :: output_damped_
 character(len=100)         :: name_dipole_time, name_excitation
 character(len=100)         :: name_dipolar_spectra, name_transforms='transforms.dat'
 character(len=100)         :: name_dipole_damped='dipole_damped.dat'
 character(len=500)         :: cur_string
 type(C_PTR)                :: plan

#ifdef HAVE_FFTW
#include "fftw3.f03"

! ntau=int((time_sim-time_min)/time_step)+1

 call get_spectrum_arguments( damp_factor,name_dipole_time,name_excitation,name_dipolar_spectra,output_damped_)

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
 open (newunit=file_dipole_time, file = name_dipole_time)
 do
   read(file_dipole_time,'(A)',iostat=io)cur_string
   if (io/=0) exit
   nrow = nrow + 1
   cur_string=ADJUSTL(cur_string)
   if(cur_string(1:1)/="#") then
     call get_number_of_elements(cur_string,num_fields)
     if(num_fields/=4) call die('Dipole file must contain 4 columns: time dx dy dz')
     ntau = ntau + 1
   end if
 end do
 close (file_dipole_time)

 allocate(m_times(ntau))
 allocate(dipole_time_ref(ntau,3))
 allocate(dipole_time_damped(ntau,3))
 allocate(m_excit_field(ntau,3))
 allocate(trans_m_excit_field(ntau,3))
 allocate(trans_dipole_time(ntau,3))
  
 ! Fill the dipole_time array
 itau=1
 open (newunit=file_dipole_time, file = name_dipole_time)
 do irow=1,nrow
   read(file_dipole_time,'(A)')cur_string
   cur_string=ADJUSTL(cur_string)
   if(cur_string(1:1)/="#") then
     read(cur_string,*)m_times(itau),real_dipole(:)
     dipole_time_ref(itau,:)=real_dipole(:) / au_debye
     itau=itau+1
   end if
 end do
 close (file_dipole_time)

 time_sim=m_times(ntau) 

 ! Get number of lines nrow and lines without "#" ntau_read in dipole_time file
 ntau_read = 0
 nrow = 0
 open (newunit=file_excitation, file = name_excitation)
 do
   read(file_excitation,'(A)',iostat=io)cur_string
   if (io/=0) exit
   nrow = nrow + 1
   cur_string=ADJUSTL(cur_string)
   if(cur_string(1:1)/="#") then
     call get_number_of_elements(cur_string,num_fields)
     if(num_fields/=4) call die('Excitation file must contain 4 columns: time Ex Ey Ez')
     ntau_read = ntau_read + 1
   end if
 end do
 close (file_excitation)

 if(ntau/=ntau_read) then
   call die('Different number of lines in the dipole_time file and excitation file')
 end if

 ! Fill the m_excit_field array
 itau=1
 open (newunit=file_excitation, file = name_excitation)
 do irow=1,nrow
   read(file_excitation,'(A)')cur_string
   cur_string=ADJUSTL(cur_string)
   if(cur_string(1:1)/="#") then
     read(cur_string,*)time_cur,real_excit_field(:)
     if(time_cur-m_times(itau)>1.0e-7_dp) then
       call die('Time values in the dipole_time file and excitation file are not the same')
     end if
     m_excit_field(itau,:)=real_excit_field(:)
     itau=itau+1
   end if
 end do
 close (file_excitation)

 ! Apply damping to the dipole
 if(output_damped_) then
   open (newunit=file_dipole_damped, file = name_dipole_damped)
   write(file_dipole_damped,*) "# time(au)                 Dipole_damped_x(D)        Dipole_damped_y(D)        Dipole_damped_z(D)   EXP(-m_times(itau) / damp_factor )"
 end if
 do itau=1,ntau
   dipole_time_damped(itau,:)=dipole_time_ref(itau,:)*EXP(-m_times(itau) / damp_factor )
   if(output_damped_) then
     write(file_dipole_damped,*) m_times(itau), REAL(dipole_time_damped(itau,:))*au_debye, EXP(-m_times(itau) / damp_factor )
   end if
 end do
 if(output_damped_) close(file_dipole_damped)

 ! Fourier transform of dipole_time_damped
 do idir=1,3
   plan = fftw_plan_dft_1d(ntau,dipole_time_damped(1:ntau,idir),trans_dipole_time(1:ntau,idir),FFTW_FORWARD,FFTW_ESTIMATE)
   call fftw_execute_dft(plan,dipole_time_damped(1:ntau,idir),trans_dipole_time(1:ntau,idir))
   call fftw_destroy_plan(plan)

 end do

 trans_dipole_time =  trans_dipole_time / ntau

 ! Fourier transform of m_excit_field
 do idir=1,3
   plan = fftw_plan_dft_1d(ntau,m_excit_field(1:ntau,idir),trans_m_excit_field(1:ntau,idir),FFTW_FORWARD,FFTW_ESTIMATE)
   call fftw_execute_dft(plan,m_excit_field(1:ntau,idir),trans_m_excit_field(1:ntau,idir))
   call fftw_destroy_plan(plan)
 end do
 
 trans_m_excit_field(:,:)=trans_m_excit_field(:,:) / ntau

 ! Write absorption spectra in the dipolar_spectra file
 open(newunit=file_transforms,file=name_transforms)
 open(newunit=file_dipolar_spectra,file=name_dipolar_spectra)
 write(file_dipolar_spectra,*) "# omega(eV), Average, xx, xy, xz, yx, yy, yz, zx, zy, zz"
 write(file_transforms,*) "# omega(eV), |E_x(omega)|, real(E_x(omega)), aimag(E_x(omega)), |d_x(omega)|, real(d_x(omega)), aimag(d_x(omega))"

 nomega=ntau
 do iomega=1,nomega/2 ! here iomega correspons to frequency
   omega_factor = 4.0_dp * pi * 2 * pi * (iomega-1) / time_sim / c_speedlight
   if(2*pi * iomega / time_sim * Ha_eV < 500.0_dp) then
     write(file_dipolar_spectra,"(x,es16.8E3)",advance='no')  2 * pi * iomega / time_sim * Ha_eV

     average=0.0_dp
     div_factor=0
     do idir=1,3
       if(ABS(trans_m_excit_field(iomega,idir))>1.0e-15) then
         div_factor=div_factor+1
         average=average + AIMAG((trans_dipole_time(iomega,idir)) / (trans_m_excit_field(iomega,idir)))
       end if
     end do
     average=average / div_factor * omega_factor
     write(file_dipolar_spectra,"(3x,es16.8E3)",advance='no') average

     do idir=1,3
       if(ABS(trans_m_excit_field(iomega,idir))>1.0e-15) then
         write(file_dipolar_spectra,"(3(3x,es16.8E3))",advance='no') AIMAG((trans_dipole_time(iomega,:))/(trans_m_excit_field(iomega,idir))) * omega_factor
       else
         write(file_dipolar_spectra,"(3(17x,f5.2))",advance='no') -1.0_dp, -1.0_dp, -1.0_dp
       end if
       if(idir==3) write(file_dipolar_spectra,*)
     end do
     write(file_transforms,*) pi * iomega / time_sim * Ha_eV, abs(trans_m_excit_field(iomega,1)), REAL(trans_m_excit_field(iomega,1)), AIMAG(trans_m_excit_field(iomega,1)), abs(trans_dipole_time(iomega,1)), REAL(trans_dipole_time(iomega,1)), AIMAG(trans_dipole_time(iomega,1))
   end if
 end do

 close(file_transforms)
 close(file_dipolar_spectra)
#else
 call issue_warning("tddft: calculate_propagation; fftw is not present")
#endif

end program spectrum

subroutine get_spectrum_arguments( damp_factor,name_dipole_time,name_excitation,name_dipolar_spectra,output_damped_)

 use m_definitions
 implicit none
 integer                    :: iarg
 real(dp)                   :: damp_factor
 character(len=20)          :: arg_cur
 character(len=100)         :: name_dipole_time, name_excitation
 character(len=100)         :: name_dipolar_spectra 
 character(len=20)          :: ch_damp_factor 
 logical                    :: output_damped_
!=====

 damp_factor=            500.0_dp
 name_dipole_time=      'dipole_time.dat'
 name_excitation=       'excitation_time.dat'
 name_dipolar_spectra=  'dipolar_spectra.dat'
 output_damped_=.false.

 do iarg=1,COMMAND_ARGUMENT_COUNT()
   call GET_COMMAND_ARGUMENT(iarg,arg_cur)
   select case (arg_cur)
   case ('-id','--input-dipole')
     call GET_COMMAND_ARGUMENT(iarg+1,name_dipole_time)
   case ('-ie','--input-excitation')
     call GET_COMMAND_ARGUMENT(iarg+1,name_excitation)
   case ('-o','--output')
     call GET_COMMAND_ARGUMENT(iarg+1,name_dipolar_spectra)
   case('-d','--damping')
     call GET_COMMAND_ARGUMENT(iarg+1,ch_damp_factor)
     read(ch_damp_factor,*)damp_factor
   case('-od','--output-damped')
     output_damped_=.true.
   case('--help')
     write(stdout,'(1x,a)') "Calculation of dipolar spectrum after RT-TDDFT simulation"
     write(stdout,'(1x,a)') "-id, --input-dipole  file        indicate an input dipole time file"
     write(stdout,'(1x,a)') "                                (if different from dipole_time.dat)"
     write(stdout,'(1x,a)') "-ie, --input-excitation file    indicate an input excitation time file"
     write(stdout,'(1x,a)') "                                (if different from excitation_time.dat)"
     write(stdout,'(1x,a)') "-o , --output file               indicate an ouput spectrum file"
     write(stdout,'(1x,a)') "                                (if different from dipolar_spectra.dat)"
     write(stdout,'(1x,a)') "-d , --damping                   damping coefficient in au of time for the Fourier transform of dipolar moment"
     write(stdout,'(1x,a)') "                               (value by default if 500 au)"
     write(stdout,'(1x,a)') "-od, --output-damped            write damped dipole time dependence in to the dipole_damped.dat file" 
     stop
   end select  
 end do

 write(stdout,'(1x,70("="))')
 write(stdout,'(12x,a)') 'This is the part of MOLGW posprocessing'
 write(stdout,'(1x,70("="))')
 write(stdout,'(/,10x,a,/)') ' === Input parameters  ==='
 write(stdout,'(1x,a,2x,es16.8)')   'Dipole damping factor       :',damp_factor
 write(stdout,'(1x,a,2x,a50)')      'Input dipole time file      :',name_dipole_time
 write(stdout,'(1x,a,2x,a50)')      'Input excitaiton field file :',name_excitation
 write(stdout,'(1x,a,2x,a50)')      'Output dipolar spectra file :',name_dipolar_spectra
 write(stdout,'(1x,a,2x,l1)')       'Output damped dipole file   :',output_damped_
 write(stdout,*)
end subroutine get_spectrum_arguments

