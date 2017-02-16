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
 use m_eri
 use m_eri_calculate
 use m_hamiltonian
 use m_hamiltonian_libint
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(inout)     :: energy(nstate,nspin)
 real(dp),intent(inout)     :: c_matrix(basis%nbf,nstate,nspin)
!=====
 integer                 :: ijshellpair,klshellpair
 integer                 :: istate,iatom,ispin
 integer                 :: ibf,jbf,kbf,lbf
 integer                 :: ishell,jshell,kshell,lshell
 integer                 :: ni,nj,nk,nl
 real(dp),allocatable    :: grad_tmp(:,:,:,:)
 real(dp),allocatable    :: q_matrix(:,:)
 real(dp),allocatable    :: grad_onebody(:,:,:)
 real(dp),allocatable    :: p_matrix(:,:,:)
 real(dp),allocatable    :: shell_gradA(:,:,:,:,:)
 real(dp),allocatable    :: shell_gradB(:,:,:,:,:)
 real(dp),allocatable    :: shell_gradC(:,:,:,:,:)
 real(dp),allocatable    :: shell_gradD(:,:,:,:,:)
!=====

#ifndef HAVE_LIBINT_ONEBODY
 call issue_warning('calculate_force: impossible to calculate gradient if LIBINT does have the gradients')
 return
#endif

 write(stdout,'(/,1x,a)') 'Calculate the forces'

 allocate(p_matrix(basis%nbf,basis%nbf,nspin))
 allocate(q_matrix(basis%nbf,basis%nbf))
 call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)
 call setup_energy_density_matrix(basis%nbf,nstate,c_matrix,occupation,energy,q_matrix)

 call dump_out_matrix(.TRUE.,'=== P-matrix ===',basis%nbf,nspin,p_matrix)
 call dump_out_matrix(.TRUE.,'=== Q-matrix ===',basis%nbf,nspin,q_matrix)


 write(stdout,'(/,1x,a)') ' ====== Pulay Forces ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'

 allocate(grad_onebody(basis%nbf,basis%nbf,3))
 call setup_overlap_grad_libint(print_matrix_,basis,grad_onebody)

 force(:,:) = 0.0_dp
 do iatom=1,natom
   do jbf=1,basis%nbf
     if( basis%bff(jbf)%iatom == iatom ) then
       force(1,iatom) = force(1,iatom) - DOT_PRODUCT( q_matrix(jbf,:) , grad_onebody(:,jbf,1) )
       force(2,iatom) = force(2,iatom) - DOT_PRODUCT( q_matrix(jbf,:) , grad_onebody(:,jbf,2) )
       force(3,iatom) = force(3,iatom) - DOT_PRODUCT( q_matrix(jbf,:) , grad_onebody(:,jbf,3) )
     endif
   enddo
   do ibf=1,basis%nbf
     if( basis%bff(ibf)%iatom == iatom ) then
       force(1,iatom) = force(1,iatom) - DOT_PRODUCT( q_matrix(ibf,:) , grad_onebody(:,ibf,1) )
       force(2,iatom) = force(2,iatom) - DOT_PRODUCT( q_matrix(ibf,:) , grad_onebody(:,ibf,2) )
       force(3,iatom) = force(3,iatom) - DOT_PRODUCT( q_matrix(ibf,:) , grad_onebody(:,ibf,3) )
     endif
   enddo
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force(:,iatom)
 enddo

 call issue_warning('force set to zero here1')
 force(:,:) = 0.0_dp

 write(stdout,'(/,1x,a)') ' ====== Pulay Forces ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'

 call setup_kinetic_grad_libint(print_matrix_,basis,grad_onebody)
 do iatom=1,natom
   do ibf=1,basis%nbf
     if( basis%bff(ibf)%iatom == iatom ) then
       force(1,iatom) = force(1,iatom) + 2.0_dp * DOT_PRODUCT( SUM( p_matrix(ibf,:,:),DIM=2) , grad_onebody(ibf,:,1) )
       force(2,iatom) = force(2,iatom) + 2.0_dp * DOT_PRODUCT( SUM( p_matrix(ibf,:,:),DIM=2) , grad_onebody(ibf,:,2) )
       force(3,iatom) = force(3,iatom) + 2.0_dp * DOT_PRODUCT( SUM( p_matrix(ibf,:,:),DIM=2) , grad_onebody(ibf,:,3) )
     endif
   enddo
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force(:,iatom)
 enddo

 deallocate(grad_onebody)

 call issue_warning('force set to zero here2')
 force(:,:) = 0.0_dp

 write(stdout,'(/,1x,a)') ' ====== Pulay Forces ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'
 allocate(grad_tmp(basis%nbf,basis%nbf,natom+1,3))
 call setup_nucleus_grad_libint(print_matrix_,basis,grad_tmp)
 do iatom=1,natom
   do ibf=1,basis%nbf
     if( basis%bff(ibf)%iatom == iatom ) then
       force(1,iatom) = force(1,iatom) + 2.0_dp * DOT_PRODUCT( SUM( p_matrix(ibf,:,:),DIM=2) , grad_tmp(ibf,:,natom+1,1) )
       force(2,iatom) = force(2,iatom) + 2.0_dp * DOT_PRODUCT( SUM( p_matrix(ibf,:,:),DIM=2) , grad_tmp(ibf,:,natom+1,2) )
       force(3,iatom) = force(3,iatom) + 2.0_dp * DOT_PRODUCT( SUM( p_matrix(ibf,:,:),DIM=2) , grad_tmp(ibf,:,natom+1,3) )
     endif
   enddo
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force(:,iatom)
 enddo
 deallocate(grad_tmp)

 call issue_warning('force set to zero here3')
 force(:,:) = 0.0_dp

 write(stdout,'(/,1x,a)') ' ====== Pulay Forces 2body ints ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'

 do iatom=1,natom
   write(*,*) ' ============ iatom', iatom
   do ijshellpair=1,nshellpair
     ishell = index_shellpair(1,ijshellpair)
     jshell = index_shellpair(2,ijshellpair)
!     if( basis%bff(ibf)%iatom /= iatom ) cycle

     ni = number_basis_function_am( basis%gaussian_type , shell(ishell)%am )
     nj = number_basis_function_am( basis%gaussian_type , shell(jshell)%am )

     do klshellpair=1,nshellpair
       kshell = index_shellpair(1,klshellpair)
       lshell = index_shellpair(2,klshellpair)
       nk = number_basis_function_am( basis%gaussian_type , shell(kshell)%am )
       nl = number_basis_function_am( basis%gaussian_type , shell(lshell)%am )

       allocate(shell_gradA(ni,nj,nk,nl,3))
       allocate(shell_gradB(ni,nj,nk,nl,3))
       allocate(shell_gradC(ni,nj,nk,nl,3))
       allocate(shell_gradD(ni,nj,nk,nl,3))
       call calculate_eri_4center_shell_grad(basis,0.0_dp,ijshellpair,klshellpair,&
                                             shell_gradA,shell_gradB,shell_gradC,shell_gradD)
       write(*,'(a,4(i3,1x),4(1x,f20.12))') ' === ',ishell,jshell,kshell,lshell,shell_gradA(:,:,:,:,1), &
                                                                                shell_gradB(:,:,:,:,1), &
                                                                                shell_gradC(:,:,:,:,1), &
                                                                                shell_gradD(:,:,:,:,1)

       deallocate(shell_gradA,shell_gradB,shell_gradC,shell_gradD)
     enddo
   enddo
 enddo



 call issue_warning('force set to zero here4')
 force(:,:) = 0.0_dp

 allocate(grad_tmp(basis%nbf,basis%nbf,natom+1,3))
 call setup_nucleus_grad_libint(print_matrix_,basis,grad_tmp)
 write(stdout,'(/,1x,a)') ' ====== Hellman Feynman Forces ====== '
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


 deallocate(p_matrix)
 deallocate(q_matrix)

end subroutine calculate_force


!=========================================================================
