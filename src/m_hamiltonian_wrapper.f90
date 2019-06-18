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
 use m_hamiltonian_buffer
 use m_hamiltonian_tools

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


end module m_hamiltonian_wrapper


!=========================================================================
