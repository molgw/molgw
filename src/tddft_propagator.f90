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
                                 energy,              &
                                 s_matrix,            & 
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
 implicit none

!=====
 type(basis_set),intent(in)      :: basis
 integer                         :: nstate 
 real(dp),intent(in)             :: energy(nstate,nspin)
 real(dp),intent(in)             :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: c_matrix(basis%nbf,nstate,nspin) 
 real(dp),intent(in)             :: occupation(nstate,nspin)
 real(dp),intent(in)             :: hamiltonian_kinetic(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!=====
 integer                    :: itau, jtau, ntau, idir, info, ispin
 integer                    :: file_time_data, file_excit_field, file_check_matrix
 integer                    :: file_dipolar_spectra, file_dipole_time 
 real(dp)                   :: t_min, time_cur, excit_field(3), eexcit
 real(dp),allocatable       :: s_matrix_inv(:,:)
 real(dp),allocatable       :: dipole_basis(:,:,:)
 real(dp)                   :: dipole(3)
 complex(dpc),allocatable   :: dipole_time(:,:)
 complex(dpc),allocatable   :: trans_dipole_time(:,:)
 complex(dpc), allocatable  :: l_matrix(:,:,:) ! Follow the notation of M.A.L.Marques, C.A.Ullrich et al, 
 complex(dpc),allocatable   :: b_matrix(:,:,:) ! TDDFT Book, Springer (2006), p205
 complex(dpc),allocatable   :: m_unity(:,:)
 complex(dpc),allocatable   :: hamiltonian_fock_cmplx(:,:,:)
 complex(dpc),allocatable   :: hamiltonian_eigen_cmplx(:,:,:)
 complex(dpc),allocatable   :: c_matrix_cmplx(:,:,:)
 complex(dpc),allocatable   :: p_matrix_cmplx(:,:,:)
 complex(dpc),allocatable   :: a_matrix_cmplx(:,:,:)
 complex(dpc),allocatable   :: check_matrix(:,:)
 character(len=36)          :: check_format
 type(C_PTR)                :: plan
#ifdef HAVE_FFTW
#include "fftw3.f03"
#endif

 call init_dft_grid(grid_level)
 
 t_min=0.0_dp
write(stdout,*) "ONLY ONCE", excit_dir
 allocate(s_matrix_inv(basis%nbf,basis%nbf))
 allocate(hamiltonian_fock_cmplx(basis%nbf,basis%nbf,nspin))
 allocate(c_matrix_cmplx(basis%nbf,nstate,nspin))
 allocate(p_matrix_cmplx(basis%nbf,basis%nbf,nspin))
 allocate(a_matrix_cmplx(basis%nbf,basis%nbf,nspin))
 allocate(l_matrix(basis%nbf,basis%nbf,nspin))
 allocate(b_matrix(basis%nbf,basis%nbf,nspin))
 allocate(m_unity(basis%nbf,basis%nbf))
 allocate(dipole_basis(basis%nbf,basis%nbf,3))
 allocate(check_matrix(basis%nbf,basis%nbf))


 ! c_matrix_cmplx(:,:,:)=c_matrix(:,:,:)
 c_matrix_cmplx(:,:,:)=c_matrix(:,:,:) 

 call setup_density_matrix_cmplx(basis%nbf,nstate,c_matrix_cmplx,occupation,p_matrix_cmplx)



 call invert(basis%nbf,s_matrix,s_matrix_inv)  
 call fill_unity(m_unity,basis%nbf)
 call calculate_dipole_basis_cmplx(basis,dipole_basis)

 ntau=int((time_sim-t_min)/time_step)
 
 allocate(dipole_time(2*ntau,3))
 allocate(trans_dipole_time(2*ntau,3))
 
 dipole_time(:,:)= ( 0.0_dp , 0.0_dp )
 
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
 do itau=1,ntau
   time_cur=t_min+itau*time_step
   eexcit=0.0_dp
   write(file_time_data,"(F9.4)",advance='no') time_cur
   write(stdout,"(A,F9.4)") "time_cur = ", time_cur

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
                                           hamiltonian_fock_cmplx,   &                   
                                           hamiltonian_kinetic,      &                
                                           hamiltonian_nucleus,      &                
                                           file_time_data,           &            
                                           file_check_matrix)         
   do ispin=1, nspin                                                  
     !------
     !--Hamiltonian - Static part--
     hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_fock_cmplx(:,:,ispin) + hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)

     !--Hamiltonian - Excitation--
     select case(excit_type)
     case('GEF')
       call calculate_excit_field(time_cur,excit_field)
       write(file_excit_field,*) time_cur, excit_field
       do idir=1,3
         hamiltonian_fock_cmplx(:,:,ispin) = hamiltonian_fock_cmplx(:,:,ispin) - &
                                   & dipole_basis(:,:,idir) * excit_field(idir)
         eexcit=eexcit+real(SUM(dipole_basis(:,:,idir)*excit_field(idir)*SUM(p_matrix_cmplx(:,:,:),DIM=3)),dp)                  
       end do     
      case default
       call die('Invalid choice for the excitation type. Change excit_type value in the input file')
     end select

  if ( print_tddft_matrices_ ) then                                   
     !--ChEcK c_MaTrIx--
    write(file_check_matrix,"(A,I2)") "ispin = ", ispin
    check_matrix = MATMUL(MATMUL(s_matrix(:,:) ,c_matrix_cmplx(:,:,ispin) ) ,TRANSPOSE(CONJG(c_matrix_cmplx(:,:,ispin)))  )
    call print_square_2d_matrix_cmplx("S*C*C**H = ",check_matrix,basis%nbf,file_check_matrix,4)
    call print_square_2d_matrix_cmplx("c_matrx_cmplx = ",c_matrix_cmplx(:,:,ispin),basis%nbf,file_check_matrix,2)
    call print_square_2d_matrix_cmplx("p_matrix_cmplx = ",p_matrix_cmplx(:,:,ispin),basis%nbf,file_check_matrix,2)
    call print_square_2d_matrix_cmplx("Hamiltonian_cmplx = ",hamiltonian_fock_cmplx(:,:,ispin),basis%nbf,file_check_matrix,4)

  end if

    !-------------------------------
    select case (prop_type)
    case('CN') ! Crank-Nocholson propagator
      l_matrix(:,:,ispin) = m_unity + im * time_step / 2.0_dp * matmul( s_matrix_inv,hamiltonian_fock_cmplx(:,:,ispin)  )
      b_matrix(:,:,ispin) = m_unity - im * time_step / 2.0_dp * matmul( s_matrix_inv,hamiltonian_fock_cmplx(:,:,ispin)  )
      call invert(basis%nbf , l_matrix(:,:,ispin))
      c_matrix_cmplx(:,:,ispin) = matmul( l_matrix(:,:,ispin),matmul( b_matrix(:,:,ispin),c_matrix_cmplx(:,:,ispin) ) )

    case('IEH') ! Instanteneous eigenvalues Hamiltonian propagator
   !   call diagonalize_cdp(basis%nbf,hamiltonian_fock_cmplx(:,:,ispin),hamiltonian_eigen(:,:,ispin),a_matrix(:,:,ispin))

    case default
      call die('Invalid choice for the propagation algorithm. Change prop_type value in the input file')
    end select
     
   end do !spin loop

  call setup_density_matrix_cmplx(basis%nbf,nstate,c_matrix_cmplx,occupation,p_matrix_cmplx)
  write(file_time_data,"(F9.4,'   ',2(2x,F7.2))") eexcit, matrix_trace_cmplx(matmul(p_matrix_cmplx(:,:,1),s_matrix(:,:)))
  call static_dipole_cmplx(nstate,basis,occupation,c_matrix_cmplx,dipole(:))
  dipole_time(itau,:)=dipole(:)
  write(file_dipole_time,*) time_cur, dipole(:) * au_debye
 end do
   !********end time loop*******************

 close(file_time_data)
 close(file_excit_field)
 close(file_dipole_time)


 if(print_tddft_matrices_) close(file_check_matrix) 


#ifdef HAVE_FFTW
 !---Fourier Transform of dipole_time---
 do idir=1,3
   plan = fftw_plan_dft_1d(2*ntau,dipole_time(:,idir),trans_dipole_time(:,idir),FFTW_FORWARD,FFTW_ESTIMATE)
   call fftw_execute_dft(plan,dipole_time(:,idir),trans_dipole_time(:,idir))
   call fftw_destroy_plan(plan)
   
 end do

 trans_dipole_time(:,:)=trans_dipole_time(:,:) / ( 2.0_dp * ntau )
 
 open(newunit=file_dipolar_spectra,file="dipolar_spectra.dat")
 do itau=1,2*ntau
   write(file_dipolar_spectra,*) pi * itau / time_sim * Ha_eV, (abs(trans_dipole_time(itau,:)))**2 
 end do


#else
 call issue_warning("tddft: calculate_propagation; fftw is not present") 
#endif
  
 call destroy_dft_grid()



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
 
 excit_field(:) = excit_kappa * exp( -( time_cur-excit_time0 )**2 / 2.0_dp / excit_omega**2 ) * &
                  & excit_dir(:) / norm_dir

end subroutine calculate_excit_field


!=======================================
subroutine fill_unity(m_unity,M)
 use m_definitions
 implicit none
 integer, intent(in) :: M
 complex(dpc) :: m_unity(M,M)
 integer :: i,j
 do i=1,M
   do j=1,M
     if (i == j) then
       m_unity(i,j)=1.0_dp
     else
      m_unity(i,j)=0.0_dp
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
 complex(dpc),intent(in)  :: matrix_cmplx(size_m,size_m)
 character(*),intent(in)  :: desc
!=====
 character(100)  :: write_format1, write_format2
 integer            :: ivar

 write(write_format1,*) '(',size_m," ('( ',F", prec+4, ".", prec,"' ,',F", prec+4, ".",prec,",' )  ') " ,')' ! (  1.01 ,  -0.03)  (  0.04 ,  0.10) 
 write(write_format2,*) '(',size_m," (F", prec+4, ".", prec,"' +  i',F", prec+4, ".",prec,",'  ') " ,')'   ! 1.01 +  i  -0.03    0.03 +  i  0.10
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

