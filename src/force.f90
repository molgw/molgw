!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the calculation of the forces (requires LIBINT with gradients)
!
!=========================================================================
subroutine calculate_force(basis,nstate,occupation,energy,c_matrix)
 use m_definitions
 use m_warning
 use m_timing
 use m_atoms
 use m_inputparam
 use m_basis_set
 use m_hamiltonian_libint
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(inout)     :: energy(nstate,nspin)
 real(dp),intent(inout)     :: c_matrix(basis%nbf,nstate,nspin)
!=====
 integer                 :: istate,iatom,ispin
 real(dp),allocatable    :: grad_tmp(:,:,:,:)
!=====

#ifndef HAVE_LIBINT_ONEBODY
 call issue_warning('calculate_force: impossible to calculate gradient if LIBINT does have the gradients')
 return
#endif


!   allocate(grad_tmp(basis%nbf,basis%nbf,1,1))
!   call setup_overlap_grad_libint(print_matrix_,basis,grad_tmp(:,:,1,1))
!   call setup_kinetic_grad_libint(print_matrix_,basis,grad_tmp(:,:,1,1))
!   deallocate(grad_tmp)

   allocate(grad_tmp(basis%nbf,basis%nbf,natom,3))
   call setup_nucleus_grad_libint(print_matrix_,basis,grad_tmp)
   write(stdout,'(/,1x,a)') ' ====== Forces ====== '
   write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'
   call nucleus_nucleus_force()
   do iatom=1,natom
     do ispin=1,nspin
       do istate=1,nstate
         force(1,iatom) = force(1,iatom)  &
              + occupation(istate,ispin) * DOT_PRODUCT( c_matrix(:,istate,ispin) , &
                   MATMUL( grad_tmp(:,:,iatom,1) , c_matrix(:,istate,ispin) ) )
         force(2,iatom) = force(2,iatom)  &
              + occupation(istate,ispin) * DOT_PRODUCT( c_matrix(:,istate,ispin) , &
                   MATMUL( grad_tmp(:,:,iatom,2) , c_matrix(:,istate,ispin) ) )
         force(3,iatom) = force(3,iatom)  &
              + occupation(istate,ispin) * DOT_PRODUCT( c_matrix(:,istate,ispin) , &
                   MATMUL( grad_tmp(:,:,iatom,3) , c_matrix(:,istate,ispin) ) )

       enddo
     enddo
     write(*,'(1x,a,i4,a,2x,3(2x,f16.8))') 'atom ',iatom,':',force(:,iatom)
   enddo
   write(stdout,'(1x,a,/)') ' ==================== '
   deallocate(grad_tmp)


end subroutine calculate_force


!=========================================================================
