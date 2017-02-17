!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the calculation of the forces (requires LIBINT with gradients)
!
!=========================================================================
subroutine calculate_force(basis,nstate,occupation,energy,c_matrix,hkin,hnuc)
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
 real(dp),intent(in)        :: energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)        :: hkin(basis%nbf,basis%nbf)
 real(dp),intent(in)        :: hnuc(basis%nbf,basis%nbf)
!=====
 integer                 :: ijshellpair,klshellpair
 integer                 :: istate,iatom,ispin
 integer                 :: ibf,jbf,kbf,lbf
 integer                 :: ishell,jshell,kshell,lshell
 integer                 :: ni,nj,nk,nl
 real(dp),allocatable    :: grad_tmp(:,:,:,:)
 real(dp),allocatable    :: q_matrix(:,:)
 real(dp),allocatable    :: psp_matrix(:,:,:,:,:)
 real(dp),allocatable    :: grad_onebody(:,:,:)
 real(dp),allocatable    :: p_matrix(:,:,:)
 real(dp),allocatable    :: shell_gradA(:,:,:,:,:)
 real(dp),allocatable    :: shell_gradB(:,:,:,:,:)
 real(dp),allocatable    :: shell_gradC(:,:,:,:,:)
 real(dp),allocatable    :: shell_gradD(:,:,:,:,:)
!=====

#ifndef HAVE_LIBINT_ONEBODY
 call issue_warning('calculate_force: impossible to calculate gradient if LIBINT does support the gradients')
 return
#endif

 write(stdout,'(/,1x,a)') 'Calculate the forces'

 allocate(p_matrix(basis%nbf,basis%nbf,nspin))
 allocate(q_matrix(basis%nbf,basis%nbf))
 allocate(psp_matrix(basis%nbf,basis%nbf,nspin,natom,3))
 call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)
 call setup_energy_density_matrix(basis%nbf,nstate,c_matrix,occupation,energy,q_matrix)


 call dump_out_matrix(.TRUE.,'=== P-matrix ===',basis%nbf,nspin,p_matrix)
 call dump_out_matrix(.TRUE.,'=== Q-matrix ===',basis%nbf,nspin,q_matrix)


 write(stdout,'(/,1x,a)') ' ====== Pulay Forces ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'



 allocate(grad_onebody(basis%nbf,basis%nbf,3))
 call setup_overlap_grad_libint(print_matrix_,basis,grad_onebody)

 psp_matrix(:,:,:,:,:) = 0.0_dp
 do iatom=1,natom
   do ibf=1,basis%nbf
     do jbf=1,basis%nbf
       if( basis%bff(ibf)%iatom == iatom ) &
         psp_matrix(ibf,jbf,1,iatom,:) = psp_matrix(ibf,jbf,1,iatom,:) + grad_onebody(ibf,jbf,:)  * 2.0_dp
     enddo
   enddo
   psp_matrix(:,:,nspin,iatom,:) = psp_matrix(:,:,1,iatom,:)
   do ispin=1,nspin
     psp_matrix(:,:,ispin,iatom,1) = -MATMUL( p_matrix(:,:,ispin) , MATMUL( psp_matrix(:,:,ispin,iatom,1) , p_matrix(:,:,ispin) ) ) / spin_fact
     psp_matrix(:,:,ispin,iatom,2) = -MATMUL( p_matrix(:,:,ispin) , MATMUL( psp_matrix(:,:,ispin,iatom,2) , p_matrix(:,:,ispin) ) ) / spin_fact
     psp_matrix(:,:,ispin,iatom,3) = -MATMUL( p_matrix(:,:,ispin) , MATMUL( psp_matrix(:,:,ispin,iatom,3) , p_matrix(:,:,ispin) ) ) / spin_fact
   enddo
 enddo
 call dump_out_matrix(.TRUE.,'=== PSP-matrix Atom 1 ===',basis%nbf,nspin,psp_matrix(:,:,1,1,1))
 call dump_out_matrix(.TRUE.,'=== PSP-matrix Atom 2 ===',basis%nbf,nspin,psp_matrix(:,:,1,2,1))


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
 deallocate(grad_onebody)

 !
 ! Kinetic energy is working !
 ! 
 call issue_warning('force set to zero here: before kinetic energy')
 force(:,:) = 0.0_dp

 write(stdout,'(/,1x,a,f26.16,/)') 'Kinetic energy',SUM( hkin(:,:) * SUM(p_matrix(:,:,:),DIM=3) )


 allocate(grad_onebody(basis%nbf,basis%nbf,3))
 call setup_kinetic_grad_libint(print_matrix_,basis,grad_onebody)

! call dump_out_matrix(.TRUE.,'=== gradX Kinetic ===',basis%nbf,nspin,grad_onebody(:,:,1))

 write(stdout,'(/,1x,a)') ' ====== Kinetic Energy Forces ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'
 do iatom=1,natom
   do ibf=1,basis%nbf
     if( basis%bff(ibf)%iatom == iatom ) then
       force(:,iatom) = force(:,iatom) + 2.0_dp * MATMUL( SUM( p_matrix(ibf,:,:),DIM=2) , grad_onebody(ibf,:,:) )
     endif
   enddo
 enddo

 deallocate(grad_onebody)

 do iatom=1,natom
   do ispin=1,nspin
     do jbf=1,basis%nbf
       do ibf=1,basis%nbf
         force(:,iatom) = force(:,iatom) + psp_matrix(ibf,jbf,ispin,iatom,:) * hkin(ibf,jbf)
       enddo
     enddo
   enddo
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force(:,iatom)
 enddo

 !
 ! Kinetic energy force is done !
 ! 


 !
 ! Nucleus energy is working !
 ! 
 call issue_warning('force set to zero here: before nucleus energy')
 force(:,:) = 0.0_dp

 write(stdout,'(/,1x,a,f26.16,/)') 'Nucleus energy',SUM( hnuc(:,:) * SUM(p_matrix(:,:,:),DIM=3) )

 allocate(grad_tmp(basis%nbf,basis%nbf,natom+2,3))
 call setup_nucleus_grad_libint(print_matrix_,basis,grad_tmp)


 write(stdout,'(/,1x,a)') ' ====== Nucleus Energy Forces ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'
 do iatom=1,natom
   do ibf=1,basis%nbf
     if( basis%bff(ibf)%iatom == iatom ) then
       force(:,iatom) = force(1,iatom) + 2.0_dp * MATMUL( SUM( p_matrix(ibf,:,:),DIM=2) , grad_tmp(ibf,:,natom+1,:) )
     endif
   enddo
 enddo

 do iatom=1,natom
   force(1,iatom) = force(1,iatom) - SUM( SUM( p_matrix(:,:,:),DIM=3 ) * grad_tmp(:,:,iatom,1) )
   force(2,iatom) = force(2,iatom) - SUM( SUM( p_matrix(:,:,:),DIM=3 ) * grad_tmp(:,:,iatom,2) )
   force(3,iatom) = force(3,iatom) - SUM( SUM( p_matrix(:,:,:),DIM=3 ) * grad_tmp(:,:,iatom,3) )
 enddo

 do iatom=1,natom
   do ispin=1,nspin
     do jbf=1,basis%nbf
       do ibf=1,basis%nbf
         force(:,iatom) = force(:,iatom) + psp_matrix(ibf,jbf,ispin,iatom,:) * hnuc(ibf,jbf)
       enddo
     enddo
   enddo
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force(:,iatom)
 enddo


 deallocate(grad_tmp)
 !
 ! Nucleus energy force section finished
 ! 


 !
 ! Hartree to be tested
 ! Fock to be tested 
 ! 
 call issue_warning('force set to zero here3')
 force(:,:) = 0.0_dp

 write(stdout,'(/,1x,a)') ' ====== Pulay Forces 2body ints ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'

 do ijshellpair=1,nshellpair
   ishell = index_shellpair(1,ijshellpair)
   jshell = index_shellpair(2,ijshellpair)

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
!     write(*,'(a,4(i3,1x),4(1x,f20.12))') ' === ',ishell,jshell,kshell,lshell,shell_gradA(:,:,:,:,1), &
!                                                                              shell_gradB(:,:,:,:,1), &
!                                                                              shell_gradC(:,:,:,:,1), &
!                                                                              shell_gradD(:,:,:,:,1)
#if 0
     do iatom=1,natom
       ibf = ishell
       jbf = jshell
       kbf = kshell
       lbf = lshell
       if( basis%bff(ibf)%iatom == iatom ) then
         write(*,*) 'HARTREE',iatom
         force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradA(1,1,1,1,1)
!         write(*,*) 'EXCHANGE',iatom
!         force(1,iatom) = force(1,iatom) -          p_matrix(ibf,kbf,1) * p_matrix(jbf,lbf,1) * shell_gradA(1,1,1,1,1)

         if( kshell /= lshell ) then
           ibf = ishell
           jbf = jshell
           kbf = lshell
           lbf = kshell
           write(*,*) 'HARTREE',iatom
          force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradA(1,1,1,1,1)
!           write(*,*) 'EXCHANGE',iatom
!           force(1,iatom) = force(1,iatom) -          p_matrix(ibf,kbf,1) * p_matrix(jbf,lbf,1) * shell_gradA(1,1,1,1,1)
         endif

       endif

       if( ishell /= jshell ) then
         ibf = jshell
         jbf = ishell
         kbf = kshell
         lbf = lshell
         if( basis%bff(ibf)%iatom == iatom ) then
           write(*,*) 'HARTREE',iatom
           force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradB(1,1,1,1,1)
         endif
         if( basis%bff(kbf)%iatom == iatom ) then
!           write(*,*) 'EXCHANGE',iatom
!           force(1,iatom) = force(1,iatom) -          p_matrix(ibf,kbf,1) * p_matrix(jbf,lbf,1) * shell_gradC(1,1,1,1,1)
         endif

         if( kshell /= lshell ) then
           ibf = jshell
           jbf = ishell
           kbf = lshell
           lbf = kshell
           if( basis%bff(ibf)%iatom == iatom ) then
             write(*,*) 'HARTREE',iatom
             force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradB(1,1,1,1,1)
           endif
           if( basis%bff(kbf)%iatom == iatom ) then
!             write(*,*) 'EXCHANGE',iatom
!             force(1,iatom) = force(1,iatom) -          p_matrix(ibf,kbf,1) * p_matrix(jbf,lbf,1) * shell_gradD(1,1,1,1,1)
           endif

         endif
       endif

     enddo
#else
     do iatom=1,natom
       ibf = ishell
       jbf = jshell
       kbf = kshell
       lbf = lshell
       if( basis%bff(ibf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradA(1,1,1,1,1)
       if( basis%bff(jbf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradB(1,1,1,1,1)
       if( basis%bff(kbf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradC(1,1,1,1,1)
       if( basis%bff(lbf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradD(1,1,1,1,1)
       if( ishell /= jshell ) then
         ibf = jshell
         jbf = ishell
         kbf = kshell
         lbf = lshell
         if( basis%bff(ibf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradB(1,1,1,1,1)
         if( basis%bff(jbf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradA(1,1,1,1,1)
         if( basis%bff(kbf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradC(1,1,1,1,1)
         if( basis%bff(lbf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradD(1,1,1,1,1)
       endif
       if( kshell /= lshell ) then
         ibf = ishell
         jbf = jshell
         kbf = lshell
         lbf = kshell
         if( basis%bff(ibf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradA(1,1,1,1,1)
         if( basis%bff(jbf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradB(1,1,1,1,1)
         if( basis%bff(kbf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradD(1,1,1,1,1)
         if( basis%bff(lbf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradC(1,1,1,1,1)
       endif
       if( ishell /= jshell .AND. kshell /= lshell ) then
         ibf = jshell
         jbf = ishell
         kbf = lshell
         lbf = kshell
         if( basis%bff(ibf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradB(1,1,1,1,1)
         if( basis%bff(jbf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradA(1,1,1,1,1)
         if( basis%bff(kbf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradD(1,1,1,1,1)
         if( basis%bff(lbf)%iatom == iatom ) force(1,iatom) = force(1,iatom) + 2.0_dp * p_matrix(ibf,jbf,1) * p_matrix(kbf,lbf,1) * shell_gradC(1,1,1,1,1)
       endif

     enddo
#endif


     deallocate(shell_gradA,shell_gradB,shell_gradC,shell_gradD)
   enddo
 enddo

 do iatom=1,natom
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force(:,iatom)
 enddo


 call issue_warning('force set to zero here4')
 force(:,:) = 0.0_dp

 allocate(grad_tmp(basis%nbf,basis%nbf,natom+2,3))
 call setup_nucleus_grad_libint(print_matrix_,basis,grad_tmp)
 write(stdout,'(/,1x,a)') ' ====== Hellman Feynman Forces ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'
 call nucleus_nucleus_force()
 do iatom=1,natom
   force(1,iatom) = force(1,iatom) - SUM( SUM( p_matrix(:,:,:),DIM=3 ) * grad_tmp(:,:,iatom,1) )
   force(2,iatom) = force(2,iatom) - SUM( SUM( p_matrix(:,:,:),DIM=3 ) * grad_tmp(:,:,iatom,2) )
   force(3,iatom) = force(3,iatom) - SUM( SUM( p_matrix(:,:,:),DIM=3 ) * grad_tmp(:,:,iatom,3) )
   write(*,'(1x,a,i4,a,2x,3(2x,f16.8))') 'atom ',iatom,':',force(:,iatom)
 enddo
 write(stdout,'(1x,a,/)') ' ==================== '
 deallocate(grad_tmp)


 deallocate(p_matrix)
 deallocate(q_matrix)

end subroutine calculate_force


!=========================================================================
