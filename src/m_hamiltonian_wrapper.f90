!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to evaluate the Kohn-Sham Hamiltonian
! with no distribution of the memory
!
!=========================================================================
module m_hamiltonian_wrapper
 use m_definitions
 use m_timing
 use m_mpi
 use m_scalapack
 use m_warning
 use m_inputparam,only: has_auxil_basis,incore_
 use m_basis_set
 use m_hamiltonian_twobodies
 use m_hamiltonian_cmplx
 use m_hamiltonian_tools
 use m_scf

 interface calculate_exchange
   module procedure calculate_exchange_real
 end interface

contains


!=========================================================================
subroutine calculate_hartree(basis,p_matrix,hhartree,eh)
 implicit none
 type(basis_set),intent(in)    :: basis
 real(dp),intent(in)           :: p_matrix(:,:,:)
 real(dp),intent(out)          :: hhartree(:,:)
 real(dp),intent(out),optional :: eh
!=====
 real(dp) :: ehartree
!=====


 !
 if( .NOT. has_auxil_basis ) then
   if( incore_ ) then
     call setup_hartree(p_matrix,hhartree,ehartree)
   else
     call setup_hartree_oneshell(basis,p_matrix,hhartree,ehartree)
   endif
 else
   call setup_hartree_ri(p_matrix,hhartree,ehartree)
 endif

 if( PRESENT(eh) ) eh = ehartree

end subroutine calculate_hartree


!=========================================================================
subroutine calculate_exchange_real(basis,p_matrix,hexx,ex,occupation,c_matrix)
 implicit none
 type(basis_set),intent(in)    :: basis
 real(dp),intent(in)           :: p_matrix(:,:,:)
 real(dp),intent(out)          :: hexx(:,:,:)
 real(dp),intent(out),optional :: ex
 real(dp),intent(in),optional  :: occupation(:,:)
 real(dp),intent(in),optional  :: c_matrix(:,:,:)
!=====
 real(dp),allocatable :: c_matrix_tmp(:,:,:)
 real(dp),allocatable :: occupation_tmp(:,:)
 real(dp)             :: eexx
!=====


 if( .NOT. has_auxil_basis ) then
   if( incore_ ) then
     call setup_exchange(p_matrix,hexx,eexx)
   else
     call issue_warning('no out-of-core exchange implemented yet')
     hexx(:,:,:) = 0.0_dp
     eexx = 0.0_dp
   endif
 else
   if( PRESENT(occupation) .AND. PRESENT(c_matrix) ) then
     call setup_exchange_ri(occupation,c_matrix,p_matrix,hexx,eexx)
   else
     !
     ! c_matrix is not provided, then calculate it from the square-root of P
     allocate(c_matrix_tmp,MOLD=p_matrix)
     allocate(occupation_tmp(SIZE(p_matrix,DIM=1),nspin))
     call get_c_matrix_from_p_matrix(p_matrix,c_matrix_tmp,occupation_tmp)
     call setup_exchange_ri(occupation_tmp,c_matrix_tmp,p_matrix,hexx,eexx)
     deallocate(c_matrix_tmp)
     deallocate(occupation_tmp)

   endif
 endif

 if( PRESENT(ex) ) ex = eexx

end subroutine calculate_exchange_real


!=========================================================================
subroutine calculate_exchange_lr(basis,p_matrix,hexx,ex,occupation,c_matrix)
 implicit none
 type(basis_set),intent(in)    :: basis
 real(dp),intent(in)           :: p_matrix(:,:,:)
 real(dp),intent(out)          :: hexx(:,:,:)
 real(dp),intent(out),optional :: ex
 real(dp),intent(in),optional  :: occupation(:,:)
 real(dp),intent(in),optional  :: c_matrix(:,:,:)
!=====
 real(dp),allocatable :: c_matrix_tmp(:,:,:)
 real(dp),allocatable :: occupation_tmp(:,:)
 real(dp)             :: eexx
!=====


 if( .NOT. has_auxil_basis ) then
   call setup_exchange_longrange(p_matrix,hexx,eexx)
 else
   if( PRESENT(occupation) .AND. PRESENT(c_matrix) ) then
     call setup_exchange_longrange_ri(occupation,c_matrix,p_matrix,hexx,eexx)
   else
     !
     ! c_matrix is not provided, then calculate it from the square-root of P
     allocate(c_matrix_tmp,MOLD=p_matrix)
     allocate(occupation_tmp(SIZE(p_matrix,DIM=1),nspin))
     call get_c_matrix_from_p_matrix(p_matrix,c_matrix_tmp,occupation_tmp)
     call setup_exchange_longrange_ri(occupation_tmp,c_matrix_tmp,p_matrix,hexx,eexx)
     deallocate(c_matrix_tmp)
     deallocate(occupation_tmp)

   endif
 endif

 if( PRESENT(ex) ) ex = eexx

end subroutine calculate_exchange_lr


!=========================================================================
subroutine calculate_hamiltonian_hxc(basis,nstate,occupation,c_matrix,p_matrix,hamiltonian_hxc,ehxc)
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(:,:,:)
 real(dp),intent(in)        :: p_matrix(:,:,:)
 real(dp),intent(out)       :: hamiltonian_hxc(:,:,:)
 real(dp),intent(out)       :: ehxc
!=====
 integer              :: ispin
 real(dp),allocatable :: hamiltonian_tmp(:,:)
 real(dp),allocatable :: hamiltonian_spin_tmp(:,:,:)
 real(dp)             :: ehart,exc,eexx,eexx_hyb
!=====


 allocate(hamiltonian_tmp,MOLD=hamiltonian_hxc(:,:,1))
 allocate(hamiltonian_spin_tmp,MOLD=hamiltonian_hxc(:,:,:))

 !
 ! Hartree contribution to the Hamiltonian
 !
 call calculate_hartree(basis,p_matrix,hamiltonian_tmp,eh=ehart)

 do ispin=1,nspin
   hamiltonian_hxc(:,:,ispin) = hamiltonian_tmp(:,:)
 enddo


 !
 !  XC part of the Hamiltonian
 !

 !
 ! DFT XC potential is added here
 !
 if( calc_type%is_dft ) then
   hamiltonian_spin_tmp(:,:,:) = 0.0_dp

   call dft_exc_vxc_batch(BATCH_SIZE,basis,occupation,c_matrix,hamiltonian_spin_tmp,exc)

   hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_tmp(:,:,:)
 endif


 !
 ! LR Exchange contribution to the Hamiltonian
 !
 if(calc_type%need_exchange_lr) then
   hamiltonian_spin_tmp(:,:,:) = 0.0_dp

   call calculate_exchange_lr(basis,p_matrix,hamiltonian_spin_tmp,ex=eexx,occupation=occupation,c_matrix=c_matrix)
   ! Rescale with alpha_hybrid_lr for range-separated hybrid functionals
   eexx_hyb = alpha_hybrid_lr * eexx
   hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_tmp(:,:,:) * alpha_hybrid_lr

 endif


 !
 ! Exchange contribution to the Hamiltonian
 !
 if( calc_type%need_exchange ) then
   hamiltonian_spin_tmp(:,:,:) = 0.0_dp

   call calculate_exchange(basis,p_matrix,hamiltonian_spin_tmp,ex=eexx,occupation=occupation,c_matrix=c_matrix)
   ! Rescale with alpha_hybrid for hybrid functionals
   eexx_hyb = eexx_hyb + alpha_hybrid * eexx
   hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_tmp(:,:,:) * alpha_hybrid

 endif

 ehxc = ehart + eexx_hyb + exc


end subroutine calculate_hamiltonian_hxc


!=========================================================================
subroutine calculate_hamiltonian_hxc_ri_cmplx(basis,                  &
                                              occupation,             &
                                              c_matrix_cmplx,         &
                                              p_matrix_cmplx,         &
                                              hamiltonian_hxc_cmplx,  &
                                              en_tddft)
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(:,:)
 complex(dp),intent(in)     :: c_matrix_cmplx(:,:,:)
 complex(dp),intent(in)     :: p_matrix_cmplx(:,:,:)
 complex(dp),intent(out)    :: hamiltonian_hxc_cmplx(:,:,:)
 type(energy_contributions),intent(inout) :: en_tddft
!=====
 integer                    :: nstate
 integer                    :: ispin
 real(dp),allocatable       :: hamiltonian_tmp(:,:,:)
!=====

 en_tddft%hart    = 0.0_dp
 en_tddft%xc      = 0.0_dp
 en_tddft%exx     = 0.0_dp
 en_tddft%exx_hyb = 0.0_dp

 nstate = SIZE(occupation,DIM=1)

 ! Initialize real arrays


 hamiltonian_hxc_cmplx = ( 0.0_dp , 0.0_dp )

 !
 ! Exchange contribution to the Hamiltonian
 !
 if( calc_type%need_exchange ) then
   call setup_exchange_ri_cmplx(occupation,c_matrix_cmplx,p_matrix_cmplx,hamiltonian_hxc_cmplx,en_tddft%exx)

   ! Rescale with alpha_hybrid for hybrid functionals
   en_tddft%exx_hyb = alpha_hybrid * en_tddft%exx
   hamiltonian_hxc_cmplx(:,:,:) = hamiltonian_hxc_cmplx(:,:,:) * alpha_hybrid
 endif


   !
   ! Hartree contribution to the Hamiltonian
   ! Hartree contribution is real and depends only on real(p_matrix) but we pass the full p_matrix_cmplx any way
   !
   call setup_hartree_ri(p_matrix_cmplx,hamiltonian_tmp(:,:,1),en_tddft%hart)

 do ispin=1,nspin
   hamiltonian_hxc_cmplx(:,:,ispin) = hamiltonian_hxc_cmplx(:,:,ispin) + hamiltonian_tmp(:,:,1)
 enddo

 !
 !  XC part of the Hamiltonian
 !

 !
 ! DFT XC potential is added here
 !
 if( calc_type%is_dft ) then
   call dft_exc_vxc_batch(BATCH_SIZE,basis,occupation,c_matrix_cmplx,hamiltonian_tmp,en_tddft%xc)

   hamiltonian_hxc_cmplx(:,:,:) = hamiltonian_hxc_cmplx(:,:,:) + hamiltonian_tmp(:,:,:)
 endif

! write(file_time_data,"(6(x,e16.10,2x),'    ')",advance='no') enuc,ekin,ehart, eexx_hyb,exc, enuc+ekin+ehart+eexx_hyb+exc
 !
 ! LR Exchange contribution to the Hamiltonian
 !
 ! if(calc_type%need_exchange_lr) then
 !   hamiltonian_spin_tmp(:,:,:) = 0.0_dp
 !
 !     call setup_exchange_longrange_ri(basis%nbf,nstate,occupation,c_matrix,p_matrix,hamiltonian_spin_tmp,eexx)
 !
 !   ! Rescale with alpha_hybrid_lr for range-separated hybrid functionals
 !   eexx_hyb = alpha_hybrid_lr * eexx
 !   hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_tmp(:,:,:) * alpha_hybrid_lr
 ! endif


end subroutine  calculate_hamiltonian_hxc_ri_cmplx


!=========================================================================
end module m_hamiltonian_wrapper
!=========================================================================
