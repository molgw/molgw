!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! propagation of the wavefunction in time
!
!=========================================================================
!module m_tddft_propagator
! use m_definitions
! use m_basis_set 
! use m_scf_loop
! use m_memory
! use m_hamiltonian
!contains


!=========================================================================
subroutine calculate_propagation(nstate,              &  
                                 basis,               &  
                                 occupation,          &   
                                 s_matrix,            & 
                                 s_matrix_sqrt_inv,   &
                                 c_matrix,            &
                                 hamiltonian_kinetic, & 
                                 hamiltonian_nucleus)
 use ISO_C_BINDING
 use m_definitions
 use m_basis_set
 use m_scf_loop
 use m_memory
 use m_hamiltonian
 use m_inputparam
 use m_dft_grid
 use m_tools
 use m_scf
 use m_timing
 implicit none

!=====
 type(basis_set),intent(in)      :: basis
 integer,intent(in)              :: nstate 
 real(dp),intent(in)             :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: c_matrix(basis%nbf,nstate,nspin) 
 real(dp),intent(in)             :: occupation(nstate,nspin)
 real(dp),intent(in)             :: hamiltonian_kinetic(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: hamiltonian_nucleus(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: s_matrix_sqrt_inv(basis%nbf,nstate)
!=====
 integer                    :: itau, jtau, ntau, idir, info, ispin, ibf
 integer                    :: file_time_data, file_excit_field, file_check_matrix
 integer                    :: file_dipolar_spectra, file_dipole_time, file_energies 
 integer                    :: file_check_dipole
 integer                    :: file_str_dipolar_spectra, file_transform_excit,file_eexcit
 real(dp)                   :: t_min, time_cur, excit_field(3) 
 real(dp),allocatable       :: s_matrix_inv(:,:)
 real(dp),allocatable       :: dipole_basis(:,:,:)
 real(dp),allocatable       :: energies_inst(:)
 real(dp)                   :: dipole(3)
 complex(dp),allocatable    :: m_excit_field(:,:)
 complex(dp),allocatable    :: dipole_time(:,:)
 complex(dp),allocatable    :: trans_dipole_time(:,:)
 complex(dp),allocatable    :: trans_m_excit_field(:,:)
 complex(dp), allocatable   :: l_matrix_cmplx(:,:,:) ! Follow the notation of M.A.L.Marques, C.A.Ullrich et al, 
 complex(dp),allocatable    :: b_matrix_cmplx(:,:,:) ! TDDFT Book, Springer (2006), p205
 complex(dp),allocatable    :: unity_matrix_cmplx(:,:)
 complex(dp),allocatable    :: hamiltonian_fock_cmplx(:,:,:)
 complex(dp),allocatable    :: propagator_eigen(:,:)
 complex(dp),allocatable    :: c_matrix_cmplx(:,:,:)
 complex(dp),allocatable    :: p_matrix_cmplx(:,:,:)
 complex(dp),allocatable    :: a_matrix_cmplx(:,:)
 complex(dp),allocatable    :: check_matrix(:,:)
 complex(dp),allocatable    :: h_small_cmplx(:,:)
 logical                    :: calc_excit_
 type(C_PTR)                :: plan
!===== variables for the calc_p_matrix_error
 complex(dp),allocatable    :: p_matrix_time_test_cmplx(:,:,:,:)
 complex(dp),allocatable    :: p_matrix_time_ref_cmplx(:,:,:,:)
 real(dp),allocatable       :: m_time_steps(:)
 integer                    :: istep, iprop_type
 character(len=3),allocatable  :: m_prop_types(:)
 integer                    :: file_p_matrix_error
 character(len=50)          :: name_p_matrix_error
 integer                    :: iwrite, nwrite, mod_write 
 integer                    :: n_time_steps,n_prop_types 
#ifdef HAVE_FFTW
#include "fftw3.f03"
#endif
 
 mod_write = INT( write_step / time_step ) ! write each write_step, assumed that time_step <= write_step

 call init_dft_grid(grid_level)
 
 t_min=0.0_dp
 allocate(s_matrix_inv(basis%nbf,basis%nbf))
 allocate(hamiltonian_fock_cmplx(basis%nbf,basis%nbf,nspin))
 allocate(c_matrix_cmplx(basis%nbf,nstate,nspin))
 allocate(p_matrix_cmplx(basis%nbf,basis%nbf,nspin))
 allocate(unity_matrix_cmplx(basis%nbf,basis%nbf))
 allocate(dipole_basis(basis%nbf,basis%nbf,3))
 allocate(check_matrix(basis%nbf,basis%nbf))
 
 select case ( prop_type )
 case('CN')
   allocate(l_matrix_cmplx(basis%nbf,basis%nbf,nspin))
   allocate(b_matrix_cmplx(basis%nbf,basis%nbf,nspin))
 case ('IEH') 
   allocate(h_small_cmplx(nstate,nstate))
   allocate(energies_inst(nstate))
   allocate(a_matrix_cmplx(basis%nbf,nstate))
   allocate(propagator_eigen(basis%nbf,basis%nbf))
   propagator_eigen(:,:) = ( 0.0_dp , 0.0_dp ) 
 case default
   call die('Invalid choice for the propagation algorithm. Change prop_type value in the input file')
 end select

 c_matrix_cmplx(:,:,:)=c_matrix(:,:,:) 

 call setup_density_matrix_cmplx(basis%nbf,nstate,c_matrix_cmplx,occupation,p_matrix_cmplx)

 call invert(basis%nbf,s_matrix,s_matrix_inv)  
 call fill_unity(unity_matrix_cmplx,basis%nbf)
 call calculate_dipole_basis_cmplx(basis,dipole_basis)

 ntau=int((time_sim-t_min)/time_step)+1
 
 allocate(dipole_time(2*ntau,3))
 allocate(trans_dipole_time(2*ntau,3))
 allocate(m_excit_field(2*ntau,3))
 allocate(trans_m_excit_field(2*ntau,3))

 if(calc_p_matrix_error_) &
          allocate(p_matrix_time_ref_cmplx(basis%nbf,basis%nbf,nspin,INT(time_sim/write_step)+1))

 dipole_time(:,:)= ( 0.0_dp , 0.0_dp )
 m_excit_field(:,:)= ( 0.0_dp , 0.0_dp )

 open(newunit=file_time_data,file="time_data.dat")
 open(newunit=file_excit_field,file="excitation.dat")
 open(newunit=file_dipole_time,file="dipole_time.dat")

 write(file_time_data,*) "# time_cur enuc  ekin  ehart  eexx_hyb  exc   e_total eexcit trace"

 if(print_tddft_matrices_) then
   open(newunit=file_check_matrix,file="check_matrix.dat")
   do idir=1,3
     call print_square_2d_matrix_real("dipole_basis", dipole_basis(:,:,idir) , basis%nbf, file_check_matrix, 3 )
   end do 
 end if


!********start time loop*************
 call start_clock(timing_tddft_loop)
 time_cur=t_min
 do itau=1,ntau
   en%excit=0.0_dp

   if(print_tddft_matrices_) then
     write(file_check_matrix,*) "========================="
     write(file_check_matrix,"(A,F9.4)") "time_cur = ", time_cur
   end if

   !--Hamiltonian - Hartree Exchange Correlation---
   call calculate_hamiltonian_hxc_ri_cmplx(basis,                    &
                                           nstate,                   &
                                           basis%nbf,                &
                                           basis%nbf,                &
                                           basis%nbf,                &
                                           nstate,                   &      
                                           occupation,               &       
                                           c_matrix_cmplx,           &          
                                           p_matrix_cmplx,           &          
                                           file_check_matrix,        &
                                           hamiltonian_fock_cmplx)                   
   do ispin=1, nspin                                                  
     !------
     !--Hamiltonian - Static part--
     hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_fock_cmplx(:,:,ispin) + hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
!      hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
     !--Hamiltonian - Excitation--
     excit_field=0.0_dp
     calc_excit_ = .false.
     calc_excit_ = calc_excit_ .OR. ( excit_type == 'GAU' )
     calc_excit_ = calc_excit_ .OR. ( excit_type == 'HSW'  .AND. abs(time_cur - excit_time0 - excit_omega/2.0_dp)<=excit_omega/2.0_dp )  
     calc_excit_ = calc_excit_ .OR. ( excit_type == 'STEP' .AND. abs(time_cur - excit_time0 - excit_omega/2.0_dp)<=excit_omega/2.0_dp )
     calc_excit_ = calc_excit_ .OR. ( excit_type == 'DEL'  .AND. abs(time_cur - excit_time0)<=time_step  ) 

     if ( calc_excit_ ) then
       call calculate_excit_field(time_cur,excit_field)
       m_excit_field(itau,:)=excit_field(:)
       do idir=1,3
         hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_fock_cmplx(:,:,ispin) - dipole_basis(:,:,idir) * excit_field(idir)
         en%excit=en%excit+real(SUM(dipole_basis(:,:,idir)*excit_field(idir)*p_matrix_cmplx(:,:,ispin)),dp)
       end do     
     end if
    

    if ( print_tddft_matrices_ ) then                                   
      write(file_check_matrix,"(A,I2)") "ispin = ", ispin
      check_matrix = MATMUL(MATMUL(s_matrix(:,:) ,c_matrix_cmplx(:,:,ispin) ) ,TRANSPOSE(CONJG(c_matrix_cmplx(:,:,ispin)))  )
      call print_square_2d_matrix_cmplx("S*C*C**H = ",check_matrix,basis%nbf,file_check_matrix,2)
      call print_square_2d_matrix_cmplx("c_matrix_cmplx = ",c_matrix_cmplx(:,:,ispin),basis%nbf,file_check_matrix,2)
      call print_square_2d_matrix_cmplx("p_matrix_cmplx = ",p_matrix_cmplx(:,:,ispin),basis%nbf,file_check_matrix,2)
      call print_square_2d_matrix_cmplx("Hamiltonian_cmplx = ",hamiltonian_fock_cmplx(:,:,ispin),basis%nbf,file_check_matrix,2)
    end if

    !-------------------------------
    select case (prop_type)

    ! Crank-Nicholson propagator
    case('CN')  
      l_matrix_cmplx(:,:,ispin) = unity_matrix_cmplx + im * time_step / 2.0_dp * matmul( s_matrix_inv,hamiltonian_fock_cmplx(:,:,ispin)  )
      b_matrix_cmplx(:,:,ispin) = unity_matrix_cmplx - im * time_step / 2.0_dp * matmul( s_matrix_inv,hamiltonian_fock_cmplx(:,:,ispin)  )
      call invert(basis%nbf , l_matrix_cmplx(:,:,ispin))

      c_matrix_cmplx(:,:,ispin) = matmul( l_matrix_cmplx(:,:,ispin),matmul( b_matrix_cmplx(:,:,ispin),c_matrix_cmplx(:,:,ispin) ) )

    ! Instanteneous eigenvalues Hamiltonian propagator
    case('IEH') 
      h_small_cmplx(:,:) = MATMUL( TRANSPOSE(s_matrix_sqrt_inv(:,:)) , &
                            MATMUL( hamiltonian_fock_cmplx(:,:,ispin) , s_matrix_sqrt_inv(:,:) ) )
      call diagonalize(nstate,h_small_cmplx(:,:),energies_inst(:),a_matrix_cmplx)
      a_matrix_cmplx(:,1:nstate) = MATMUL( s_matrix_sqrt_inv(:,:) , a_matrix_cmplx(:,:) )
    
      do ibf=1,basis%nbf
        propagator_eigen(ibf,ibf) = exp(-im*time_step*energies_inst(ibf))
      end do   
      c_matrix_cmplx(:,:,ispin) = matmul( matmul( matmul( matmul( a_matrix_cmplx(:,:), propagator_eigen(:,:)  ) , &
              conjg(transpose(a_matrix_cmplx(:,:)))  ), s_matrix(:,:) ), c_matrix_cmplx(:,:,ispin) )

      if ( print_tddft_matrices_ ) then
        ! test of diagonalization
        call print_square_2d_matrix_cmplx("propagator_eigen = ",propagator_eigen(:,:),basis%nbf,file_check_matrix,8)
        call print_square_2d_matrix_cmplx("a_matrx_cmplx = ",a_matrix_cmplx(:,:),basis%nbf,file_check_matrix,4)
        call print_square_2d_matrix_cmplx("h_small_cmplx = ",h_small_cmplx(:,:),basis%nbf,file_check_matrix,4)
        call print_square_2d_matrix_cmplx("test of diagonalization = ", &
                 matmul(matmul(transpose(conjg(a_matrix_cmplx)),hamiltonian_fock_cmplx(:,:,ispin)),a_matrix_cmplx),basis%nbf,file_check_matrix,2)
        write(file_check_matrix,*) "energies: ", energies_inst(:)
      end if

    case default
      call die('Invalid choice for the propagation algorithm. Change prop_type value in the input file')
    end select
     
   end do !spin loop

   call setup_density_matrix_cmplx(basis%nbf,nstate,c_matrix_cmplx,occupation,p_matrix_cmplx)

   en%kin = real(SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
   en%nuc = real(SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
   en%tot = en%nuc + en%kin + en%nuc_nuc + en%hart + en%exx_hyb + en%xc

   call static_dipole_cmplx(nstate,basis,occupation,c_matrix_cmplx,dipole(:))
   dipole_time(itau,:)=dipole(:)

   if(mod(itau-1,mod_write)==0) then
     if( print_cube_rho_tddft_ ) call plot_cube_wfn_cmplx(nstate,basis,occupation,c_matrix_cmplx,itau)
     if( calc_p_matrix_error_ ) p_matrix_time_ref_cmplx(:,:,:,itau/mod_write+1)=p_matrix_cmplx(:,:,:)
     write(file_excit_field,*) time_cur, excit_field
     write(file_dipole_time,*) time_cur, dipole(:) * au_debye
     write(file_time_data,"(F9.4,7(2x,es16.8E3),2x,2(2x,F7.2))") &
        time_cur, en%tot, en%nuc, en%kin, en%hart, en%exx_hyb, en%xc, en%excit, matrix_trace_cmplx(matmul(p_matrix_cmplx(:,:,1),s_matrix(:,:)))
   end if
   time_cur=t_min+itau*time_step
 end do
 call stop_clock(timing_tddft_loop)
!********end time loop*******************

 close(file_time_data)
 close(file_excit_field)
 close(file_dipole_time)

 if(print_tddft_matrices_) close(file_check_matrix) 


#ifdef HAVE_FFTW
 call start_clock(timing_tddft_fourier)
 !---Fourier Transform of dipole_time---
 do idir=1,3
   plan = fftw_plan_dft_1d(ntau,dipole_time(1:ntau,idir),trans_dipole_time(1:ntau,idir),FFTW_FORWARD,FFTW_ESTIMATE)
   call fftw_execute_dft(plan,dipole_time(1:ntau,idir),trans_dipole_time(1:ntau,idir))
   call fftw_destroy_plan(plan)
   
 end do


 do idir=1,3
   plan = fftw_plan_dft_1d(ntau,trans_dipole_time(1:ntau,idir),dipole_time(1:ntau,idir),FFTW_FORWARD,FFTW_ESTIMATE)
   call fftw_execute_dft(plan,trans_dipole_time(1:ntau,idir),dipole_time(1:ntau,idir))
   call fftw_destroy_plan(plan)

 end do

 open(newunit=file_transform_excit,file="transform_excit.dat")
 open(newunit=file_str_dipolar_spectra,file="str_dipolar_spectra.dat")
 open(newunit=file_dipolar_spectra,file="dipolar_spectra.dat")
 do idir=1,3
   plan = fftw_plan_dft_1d(2*ntau,m_excit_field(:,idir),trans_m_excit_field(:,idir),FFTW_FORWARD,FFTW_ESTIMATE)
   call fftw_execute_dft(plan,m_excit_field(:,idir),trans_m_excit_field(:,idir))
   call fftw_destroy_plan(plan)
 end do

 trans_m_excit_field(:,:)=trans_m_excit_field(:,:) / ( 2.0_dp * ntau )

 write(file_dipolar_spectra,*) "#abs(trans_dipole_time(i))/abs(trans_m_excit_field(j))"
 write(file_dipolar_spectra,*) "# omega(eV), Average, xx, xy, xz, yx, yy, yz, zx, zy, zz"
 do itau=1,ntau ! here itau correspons to frequency
    if(2*pi * itau / time_sim * Ha_eV < 1000.0_dp) then
  !   write(file_dipolar_spectra,*) pi*itau / time_sim * Ha_eV , abs(trans_dipole_time(itau,1))/abs(trans_m_excit_field(itau,1)), &
  !   aimag(trans_dipole_time(itau,1)/trans_m_excit_field(itau,1))
     write(file_dipolar_spectra,"(x,es16.8E3)",advance='no') 2*pi * itau / time_sim * Ha_eV
     write(file_str_dipolar_spectra,"(x,es16.8E3)",advance='no') 2*pi * itau / time_sim * Ha_eV
  
     write(file_dipolar_spectra,"(3x,es16.8E3)",advance='no') &
             SUM( abs(trans_dipole_time(itau,:)) / abs(trans_m_excit_field(itau,:)) ) / 3.0_dp
     write(file_str_dipolar_spectra,"(3x,es16.8E3)",advance='no') &
             SUM( aimag((trans_dipole_time(itau,:)) / (trans_m_excit_field(itau,:))) ) / 3.0_dp
  
     do idir=1,3
       if( excit_dir(idir) > 1.0E-10_dp ) then
  !       4.0_dp * pi * REAL(omega(iomega),dp) / c_speedlight 
         write(file_dipolar_spectra,"(3(3x,es16.8E3))",advance='no') abs(trans_dipole_time(itau,:))/abs(trans_m_excit_field(itau,idir))
         write(file_str_dipolar_spectra,"(3(3x,es16.8E3))",advance='no') aimag((trans_dipole_time(itau,:))/(trans_m_excit_field(itau,idir)))
       else
         write(file_dipolar_spectra,"(3(es16.8E3))",advance='no') -1.0_dp, -1.0_dp, -1.0_dp
         write(file_str_dipolar_spectra,"(3(es16.8E3))",advance='no') 0.0_dp, 0.0_dp, 0.0_dp
       end if
       if(idir==3) write(file_dipolar_spectra,*)
       if(idir==3) write(file_str_dipolar_spectra,*)
     end do
     write(file_transform_excit,*) pi * itau / time_sim * Ha_eV, abs(trans_m_excit_field(itau,:))
   end if
 end do

 close(file_transform_excit)
 close(file_dipolar_spectra)
 close(file_str_dipolar_spectra)
 call stop_clock(timing_tddft_fourier)
#else
 call issue_warning("tddft: calculate_propagation; fftw is not present") 
#endif
  
 call destroy_dft_grid()

 
 if( calc_p_matrix_error_ ) then

   allocate(p_matrix_time_test_cmplx(basis%nbf,basis%nbf,nspin,INT(time_sim/write_step)+1))
   call get_number_of_elements(error_prop_types,n_prop_types)
   allocate(m_prop_types(n_prop_types))
   read(error_prop_types,*)m_prop_types(:) 
  
   call get_number_of_elements(error_time_steps,n_time_steps)
   allocate(m_time_steps(n_time_steps))
   read(error_time_steps,*)m_time_steps(:)

   do istep=1, size(m_time_steps)
     do iprop_type=1, size(m_prop_types)
       call calculate_p_matrix_time_cmplx(nstate,                           &
                                          basis,                            &
                                          occupation,                       &
                                          s_matrix,                         &
                                          s_matrix_inv,                     &
                                          s_matrix_sqrt_inv,                &
                                          c_matrix,                         &
                                          hamiltonian_kinetic,              &
                                          hamiltonian_nucleus,              &
                                          m_prop_types(iprop_type),         &
                                          m_time_steps(istep),              &
                                          p_matrix_time_test_cmplx)
       write(name_p_matrix_error,'(3A,F6.4,A)') "p_matrix_error_", TRIM(m_prop_types(iprop_type)), "_step_", m_time_steps(istep), ".dat" 
       open(newunit=file_p_matrix_error,file=name_p_matrix_error)
       time_cur=t_min
       do itau=1, INT(time_sim/write_step)+1 
         write(file_p_matrix_error,*) time_cur, SUM(ABS(p_matrix_time_ref_cmplx(:,:,:,itau) - p_matrix_time_test_cmplx(:,:,:,itau))) / (basis%nbf)**2
         time_cur=t_min+itau*write_step
       end do
       close(file_p_matrix_error)     
     end do
   end do

 deallocate(p_matrix_time_test_cmplx)
 end if

end subroutine calculate_propagation


!==========================================
subroutine calculate_excit_field(time_cur,excit_field)
 use m_definitions
 use m_inputparam
 implicit none
 real(dp),intent(in)      :: time_cur ! time in au
 real(dp),intent(inout)   :: excit_field(3) ! electric field in 3 dimensions
!=====
 real(dp)    :: norm_dir
!=====

 norm_dir=NORM2(excit_dir(:))

 select case(excit_type)
 case('GAU') !Gaussian electic field
   excit_field(:) = excit_kappa * exp( -( time_cur-excit_time0 )**2 / 2.0_dp / excit_omega**2 ) * &
                  & excit_dir(:) / norm_dir
 case('HSW') !Hann sine window
   excit_field(:) = excit_kappa * sin( pi / excit_omega * ( time_cur - excit_time0  ) )**2 * excit_dir(:) / norm_dir
 case('DEL') ! Delta excitation
   excit_field(:) = excit_kappa * excit_dir(:) / norm_dir
 case('STEP') ! Step excitation
   excit_field(:) = excit_kappa * excit_dir(:) / norm_dir
 case default
    call die('Invalid choice for the excitation type. Change excit_type value in the input file')
 end select
 
end subroutine calculate_excit_field


subroutine get_number_of_elements(string,num)
 implicit none
 character(len=100),intent(in)  ::  string
 integer,intent(inout)          :: num
!===
 integer   :: i,pos

 pos=1
 num=0

 do
   i=verify(string(pos:),' ')    !-- Find next non-blank 
   if (i==0) exit                !-- No word found
   num=num+1                     !-- Found something
   pos=pos+i-1                   !-- Move to start of the word 
   i=scan(string(pos:),' ')      !-- Find next blank 
   if (i==0) exit                !-- No blank found
   pos=pos+i-1                   !-- Move to the blank
 end do

end subroutine get_number_of_elements



!=======================================
subroutine fill_unity(unity_matrix_cmplx,M)
 use m_definitions
 implicit none
 integer, intent(in) :: M
 complex(dp) :: unity_matrix_cmplx(M,M)
 integer :: i,j
 do i=1,M
   do j=1,M
     if (i == j) then
       unity_matrix_cmplx(i,j)=1.0_dp
     else
      unity_matrix_cmplx(i,j)=0.0_dp
     end if
   end do
 end do
end subroutine

!=========================================================================
!end module m_tddft_propagator
!=========================================================================


!=======================================
subroutine print_square_2d_matrix_cmplx (desc,matrix_cmplx,size_m,write_unit,prec)
 use m_definitions
 implicit none
 integer, intent(in)      :: prec ! precision
 integer, intent(in)      :: size_m, write_unit
 complex(dp),intent(in)  :: matrix_cmplx(size_m,size_m)
 character(*),intent(in)  :: desc
!=====
 character(100)  :: write_format1, write_format2, write_form
 integer            :: ivar,beg

 beg=4
! beg=3
 write(write_format1,*) '(',size_m," ('( ',F", prec+beg, ".", prec,"' ,',F", prec+beg, ".",prec,",' )  ') " ,')' ! (  1.01 ,  -0.03)  (  0.04 ,  0.10) 
 write(write_format2,*) '(',size_m," (F", prec+beg, ".", prec,"' +  i',F", prec+beg, ".",prec,",'  ') " ,')'   ! 1.01 +  i  -0.03    0.03 +  i  0.10
 write(write_unit,*) desc
 do ivar=1,size_m
   write(write_unit,write_format1) matrix_cmplx(ivar,:)          
 end do

end subroutine print_square_2d_matrix_cmplx


!=======================================
subroutine print_square_2d_matrix_real (desc,matrix_real,size_m,write_unit,prec)
 use m_definitions
 implicit none
 integer, intent(in)      :: prec ! precision
 integer, intent(in)      :: size_m, write_unit
 real(dp),intent(in)      :: matrix_real(size_m,size_m)
 character(*),intent(in)  :: desc
!=====
 character(100)  :: write_format1
 integer            :: ivar

 write(write_format1,*) '(',size_m," (F", prec+4, ".", prec,') ' ,')' 
 write(write_unit,*) desc
 do ivar=1,size_m
   write(write_unit,write_format1) matrix_real(ivar,:)          
 end do
end subroutine print_square_2d_matrix_real


!=======================================
subroutine calculate_p_matrix_time_cmplx(nstate,                           &
                                         basis,                            &
                                         occupation,                       &
                                         s_matrix,                         &
                                         s_matrix_inv,                     &
                                         s_matrix_sqrt_inv,                &
                                         c_matrix,                         &
                                         hamiltonian_kinetic,              &
                                         hamiltonian_nucleus,              &
                                         prop_type_cur,                    &
                                         time_step_cur,                    &
                                         p_matrix_time_cmplx)
 use m_definitions
 use m_basis_set
 use m_scf_loop
 use m_memory
 use m_hamiltonian
 use m_inputparam
 use m_dft_grid
 use m_tools
 use m_scf
 implicit none

!=====
 type(basis_set),intent(in)      :: basis
 integer,intent(in)              :: nstate
 real(dp),intent(in)             :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: c_matrix(basis%nbf,nstate,nspin) 
 real(dp),intent(in)             :: occupation(nstate,nspin)
 real(dp),intent(in)             :: hamiltonian_kinetic(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: hamiltonian_nucleus(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: s_matrix_sqrt_inv(basis%nbf,nstate)
 real(dp),intent(in)             :: s_matrix_inv(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: time_step_cur
 character(len=3),intent(in)     :: prop_type_cur
 complex(dp),intent(inout)  :: p_matrix_time_cmplx(basis%nbf,basis%nbf,nspin,INT(time_sim/write_step)+1)
!=====
 integer                    :: itau, jtau, idir, info, ispin, ibf, ntau, mod_write
 real(dp)                   :: time_cur, excit_field(3),t_min
 real(dp),allocatable       :: dipole_basis(:,:,:)
 real(dp),allocatable       :: energies_inst(:)
 real(dp)                   :: dipole(3)
 complex(dp),allocatable    :: c_matrix_cmplx(:,:,:)
 complex(dp),allocatable    :: dipole_time(:,:)
 complex(dp),allocatable    :: m_excit_field(:,:)
 complex(dp), allocatable   :: l_matrix_cmplx(:,:,:) ! Follow the notation of M.A.L.Marques, C.A.Ullrich et al, 
 complex(dp),allocatable    :: b_matrix_cmplx(:,:,:) ! TDDFT Book, Springer (2006), p205
 complex(dp),allocatable    :: unity_matrix_cmplx(:,:)
 complex(dp),allocatable    :: hamiltonian_fock_cmplx(:,:,:)
 complex(dp),allocatable    :: propagator_eigen(:,:)
 complex(dp),allocatable    :: p_matrix_cmplx(:,:,:)
 complex(dp),allocatable    :: a_matrix_cmplx(:,:)
 complex(dp),allocatable    :: check_matrix(:,:)
 complex(dp),allocatable    :: h_small_cmplx(:,:)
 logical                    :: calc_excit_
!=====

 ntau=int((time_sim-t_min)/time_step_cur)+1
 mod_write = INT( write_step / time_step_cur )

 call init_dft_grid(grid_level)

 t_min=0.0_dp
 allocate(hamiltonian_fock_cmplx(basis%nbf,basis%nbf,nspin))
 allocate(p_matrix_cmplx(basis%nbf,basis%nbf,nspin))
 allocate(unity_matrix_cmplx(basis%nbf,basis%nbf))
 allocate(dipole_basis(basis%nbf,basis%nbf,3))
 allocate(check_matrix(basis%nbf,basis%nbf))
 allocate(c_matrix_cmplx(basis%nbf,nstate,nspin))
 allocate(m_excit_field(2*ntau,3))

 select case ( prop_type_cur )
 case('CN')
   allocate(l_matrix_cmplx(basis%nbf,basis%nbf,nspin))
   allocate(b_matrix_cmplx(basis%nbf,basis%nbf,nspin))
 case ('IEH') 
   allocate(h_small_cmplx(nstate,nstate))
   allocate(energies_inst(nstate))
   allocate(a_matrix_cmplx(basis%nbf,nstate))
   allocate(propagator_eigen(basis%nbf,basis%nbf))
   propagator_eigen(:,:) = ( 0.0_dp , 0.0_dp )
 case default
   call die('Invalid choice for the propagation algorithm. Change prop_type value in the input file')
 end select

 c_matrix_cmplx(:,:,:)=c_matrix(:,:,:)

 call setup_density_matrix_cmplx(basis%nbf,nstate,c_matrix_cmplx,occupation,p_matrix_cmplx)

 call fill_unity(unity_matrix_cmplx,basis%nbf)
 call calculate_dipole_basis_cmplx(basis,dipole_basis)

 
 allocate(dipole_time(2*ntau,3))

 dipole_time(:,:)= ( 0.0_dp , 0.0_dp )
 m_excit_field(:,:)= ( 0.0_dp , 0.0_dp )

!********start time loop*************
 time_cur=t_min
 do itau=1,ntau
   en%excit=0.0_dp
   !--Hamiltonian - Hartree Exchange Correlation---
   call calculate_hamiltonian_hxc_ri_cmplx(basis,                    &
                                           nstate,                   &
                                           basis%nbf,                &
                                           basis%nbf,                &
                                           basis%nbf,                &
                                           nstate,                   &      
                                           occupation,               &       
                                           c_matrix_cmplx,           &          
                                           p_matrix_cmplx,           &          
                                           0,                        &
                                           hamiltonian_fock_cmplx)                   
   do ispin=1, nspin                                                  
     !------
     !--Hamiltonian - Static part--
     hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_fock_cmplx(:,:,ispin) + hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
!      hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
     !--Hamiltonian - Excitation--
     excit_field=0.0_dp
     calc_excit_ = .false.
     calc_excit_ = calc_excit_ .OR. ( excit_type == 'GAU' )
     calc_excit_ = calc_excit_ .OR. ( excit_type == 'HSW'  .AND. abs(time_cur - excit_time0 - excit_omega/2.0_dp)<=excit_omega/2.0_dp )  
     calc_excit_ = calc_excit_ .OR. ( excit_type == 'STEP' .AND. abs(time_cur - excit_time0 - excit_omega/2.0_dp)<=excit_omega/2.0_dp )
     calc_excit_ = calc_excit_ .OR. ( excit_type == 'DEL'  .AND. abs(time_cur - excit_time0)<=time_step_cur ) 

     if ( calc_excit_ ) then
       call calculate_excit_field(time_cur,excit_field)
       m_excit_field(itau,:)=excit_field(:)
       do idir=1,3
         hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_fock_cmplx(:,:,ispin) - dipole_basis(:,:,idir) * excit_field(idir)
         en%excit=en%excit+real(SUM(dipole_basis(:,:,idir)*excit_field(idir)*p_matrix_cmplx(:,:,ispin)),dp)
       end do     
     end if
    
    !-------------------------------
    select case (prop_type_cur)

    ! Crank-Nicholson propagator
    case('CN')  
      l_matrix_cmplx(:,:,ispin) = unity_matrix_cmplx + im * time_step_cur / 2.0_dp * matmul( s_matrix_inv,hamiltonian_fock_cmplx(:,:,ispin)  )
      b_matrix_cmplx(:,:,ispin) = unity_matrix_cmplx - im * time_step_cur / 2.0_dp * matmul( s_matrix_inv,hamiltonian_fock_cmplx(:,:,ispin)  )
      call invert(basis%nbf , l_matrix_cmplx(:,:,ispin))

      c_matrix_cmplx(:,:,ispin) = matmul( l_matrix_cmplx(:,:,ispin),matmul( b_matrix_cmplx(:,:,ispin),c_matrix_cmplx(:,:,ispin) ) )

    ! Instanteneous eigenvalues Hamiltonian propagator
    case('IEH') 
      h_small_cmplx(:,:) = MATMUL( TRANSPOSE(s_matrix_sqrt_inv(:,:)) , &
                            MATMUL( hamiltonian_fock_cmplx(:,:,ispin) , s_matrix_sqrt_inv(:,:) ) )
      call diagonalize(nstate,h_small_cmplx(:,:),energies_inst(:),a_matrix_cmplx)
      a_matrix_cmplx(:,1:nstate) = MATMUL( s_matrix_sqrt_inv(:,:) , a_matrix_cmplx(:,:) )
    
      do ibf=1,basis%nbf
        propagator_eigen(ibf,ibf) = exp(-im*time_step_cur*energies_inst(ibf))
      end do   
      c_matrix_cmplx(:,:,ispin) = matmul( matmul( matmul( matmul( a_matrix_cmplx(:,:), propagator_eigen(:,:)  ) , &
              conjg(transpose(a_matrix_cmplx(:,:)))  ), s_matrix(:,:) ), c_matrix_cmplx(:,:,ispin) )
    case default
      call die('Invalid choice for the propagation algorithm. Change prop_type_cur value in the input file')
    end select
     
   end do !spin loop

   call setup_density_matrix_cmplx(basis%nbf,nstate,c_matrix_cmplx,occupation,p_matrix_cmplx)
   if(mod(itau-1,mod_write)==0) p_matrix_time_cmplx(:,:,:,itau/mod_write+1)=p_matrix_cmplx(:,:,:); write(stdout,*) "error ", time_step_cur, time_cur

   en%kin = real(SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
   en%nuc = real(SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
   en%tot = en%nuc + en%kin + en%nuc_nuc + en%hart + en%exx_hyb + en%xc

   call static_dipole_cmplx(nstate,basis,occupation,c_matrix_cmplx,dipole(:))
   dipole_time(itau,:)=dipole(:)
   time_cur=t_min+itau*time_step_cur
 end do
   !********end time loop*******************

 call destroy_dft_grid()

end subroutine calculate_p_matrix_time_cmplx

