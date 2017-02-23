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
 integer                 :: istate,ispin
 integer                 :: iatom,jatom,katom,latom
 integer                 :: ibf,jbf,kbf,lbf
 integer                 :: ishell,jshell,kshell,lshell
 integer                 :: ni,nj,nk,nl
 real(dp)                :: fact
 real(dp),allocatable    :: grad_tmp(:,:,:,:)
 real(dp),allocatable    :: r_matrix(:,:)
#ifdef PSP
 real(dp),allocatable    :: psp_matrix(:,:,:,:,:)
#endif
 real(dp),allocatable    :: grad_onebody(:,:,:)
 real(dp),allocatable    :: p_matrix(:,:,:)
 real(dp),allocatable    :: shell_gradA(:,:,:,:,:)
 real(dp),allocatable    :: shell_gradB(:,:,:,:,:)
 real(dp),allocatable    :: shell_gradC(:,:,:,:,:)
 real(dp),allocatable    :: shell_gradD(:,:,:,:,:)
#if 0
 !debug FBFB
 real(dp),allocatable    :: shellABCD(:,:,:,:)
 real(dp),allocatable    :: shellABCD1(:,:,:,:)
 real(dp),allocatable    :: shellABCD2(:,:,:,:)
 real(dp),allocatable    :: shellABCD3(:,:,:,:)
 real(dp),allocatable    :: shellABCD4(:,:,:,:)
 real(dp),parameter      :: dx=1.0e-5_dp
 integer :: amtot
#endif
!=====

#ifndef HAVE_LIBINT_ONEBODY
 call issue_warning('calculate_force: impossible to calculate gradient if LIBINT does support the gradients')
 return
#endif

 write(stdout,'(/,1x,a)') 'Calculate the forces'

 allocate(p_matrix(basis%nbf,basis%nbf,nspin))
 allocate(r_matrix(basis%nbf,basis%nbf))
#ifdef PSP
 allocate(psp_matrix(basis%nbf,basis%nbf,nspin,natom,3))
#endif
 call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)
 call setup_energy_density_matrix(basis%nbf,nstate,c_matrix,occupation,energy,r_matrix)

 call dump_out_matrix(.TRUE.,'=== P-matrix ===',basis%nbf,nspin,p_matrix)


 allocate(grad_onebody(basis%nbf,basis%nbf,3))
 call setup_overlap_grad_libint(print_matrix_,basis,grad_onebody)

#ifdef PSP
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
#endif

 write(stdout,'(/,1x,a)') ' ====== Overlap Forces ====== '
 force_ovp(:,:) = 0.0_dp
 do ibf=1,basis%nbf
   iatom = basis%bff(ibf)%iatom
   force_ovp(:,iatom) = force_ovp(:,iatom) + 2.0_dp * MATMUL( r_matrix(ibf,:) , grad_onebody(:,ibf,:) )
 enddo
 do iatom=1,natom
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force_ovp(:,iatom)
 enddo
 deallocate(grad_onebody)

 !
 ! Kinetic energy is working !
 ! 
 force_kin(:,:) = 0.0_dp

 write(stdout,'(/,1x,a,f26.16,/)') 'Kinetic energy',SUM( hkin(:,:) * SUM(p_matrix(:,:,:),DIM=3) )


 allocate(grad_onebody(basis%nbf,basis%nbf,3))
 call setup_kinetic_grad_libint(print_matrix_,basis,grad_onebody)

 write(stdout,'(/,1x,a)') ' ====== Kinetic Energy Forces ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'
 do ibf=1,basis%nbf
   iatom = basis%bff(ibf)%iatom
   force_kin(:,iatom) = force_kin(:,iatom) + 2.0_dp * MATMUL( SUM( p_matrix(ibf,:,:),DIM=2) , grad_onebody(ibf,:,:) )
 enddo

 deallocate(grad_onebody)

#ifdef PSP
 do iatom=1,natom
   do ispin=1,nspin
     do jbf=1,basis%nbf
       do ibf=1,basis%nbf
         force_kin(:,iatom) = force_kin(:,iatom) + psp_matrix(ibf,jbf,ispin,iatom,:) * hkin(ibf,jbf)
       enddo
     enddo
   enddo
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force_kin(:,iatom)
 enddo
#endif

 !
 ! Kinetic energy force is done !
 ! 


 !
 ! Nucleus energy is working !
 ! 
 force_nuc(:,:) = 0.0_dp
 force_hl(:,:) = 0.0_dp

 write(stdout,'(/,1x,a,f26.16,/)') 'Nucleus energy',SUM( hnuc(:,:) * SUM(p_matrix(:,:,:),DIM=3) )

 allocate(grad_tmp(basis%nbf,basis%nbf,natom+1,3))
 call setup_nucleus_grad_libint(print_matrix_,basis,grad_tmp)

 write(stdout,'(/,1x,a)') ' ====== Nucleus Energy Forces ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'
 do ibf=1,basis%nbf
   iatom = basis%bff(ibf)%iatom
   force_nuc(:,iatom) = force_nuc(:,iatom) + 2.0_dp * MATMUL( SUM( p_matrix(ibf,:,:),DIM=2) , grad_tmp(ibf,:,natom+1,:) )
 enddo

 do iatom=1,natom
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'NUC atom ',iatom,':',force_nuc(:,iatom)
 enddo

 do iatom=1,natom
   force_hl(1,iatom) = force_hl(1,iatom) + SUM( SUM( p_matrix(:,:,:),DIM=3 ) * grad_tmp(:,:,iatom,1) )
   force_hl(2,iatom) = force_hl(2,iatom) + SUM( SUM( p_matrix(:,:,:),DIM=3 ) * grad_tmp(:,:,iatom,2) )
   force_hl(3,iatom) = force_hl(3,iatom) + SUM( SUM( p_matrix(:,:,:),DIM=3 ) * grad_tmp(:,:,iatom,3) )
 enddo
 force_nuc(:,:) = force_nuc(:,:) + force_hl(:,:)

#ifdef PSP
 do iatom=1,natom
   do ispin=1,nspin
     do jbf=1,basis%nbf
       do ibf=1,basis%nbf
         force_nuc(:,iatom) = force_nuc(:,iatom) + psp_matrix(ibf,jbf,ispin,iatom,:) * hnuc(ibf,jbf)
       enddo
     enddo
   enddo
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force_nuc(:,iatom)
 enddo
#endif


 deallocate(grad_tmp)
 !
 ! Nucleus energy force section finished
 ! 


 !
 ! Hartree ok!
 ! Fock ok!
 ! 
 force_har(:,:) = 0.0_dp
 force_exx(:,:) = 0.0_dp

 write(stdout,'(/,1x,a)') ' ====== Pulay Forces 2body ints ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'

#ifdef PSP
 do iatom=1,natom
   do ispin=1,nspin
     do lbf=1,basis%nbf
       do kbf=1,basis%nbf
         do jbf=1,basis%nbf
           do ibf=1,basis%nbf
             ! Hartree
             force_har(:,iatom) = force_har(:,iatom) + psp_matrix(ibf,jbf,ispin,iatom,:) * eri(ibf,jbf,kbf,lbf) * p_matrix(kbf,lbf,ispin)
             ! Exchange
             force_exx(:,iatom) = force_exx(:,iatom) &
                     - psp_matrix(ibf,kbf,ispin,iatom,:) * eri(ibf,jbf,kbf,lbf) * p_matrix(jbf,lbf,ispin) / spin_fact
           enddo
         enddo
       enddo
     enddo
   enddo
 enddo
#endif


#if 0
 ! FBFB section to be removed
 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)
   nk = number_basis_function_am( basis%gaussian_type , shell(kshell)%am )
   nl = number_basis_function_am( basis%gaussian_type , shell(lshell)%am )
   do ijshellpair=1,nshellpair
     ishell = index_shellpair(1,ijshellpair)
     jshell = index_shellpair(2,ijshellpair)
     ni = number_basis_function_am( basis%gaussian_type , shell(ishell)%am )
     nj = number_basis_function_am( basis%gaussian_type , shell(jshell)%am )
     if( shell(ishell)%am + shell(jshell)%am > shell(kshell)%am + shell(lshell)%am ) cycle

     if( shell(ishell)%iatom == shell(jshell)%iatom ) cycle
     if( shell(ishell)%iatom == shell(kshell)%iatom ) cycle
     if( shell(ishell)%iatom == shell(lshell)%iatom ) cycle
     if( shell(jshell)%iatom == shell(kshell)%iatom ) cycle
     if( shell(jshell)%iatom == shell(lshell)%iatom ) cycle
     if( shell(kshell)%iatom == shell(lshell)%iatom ) cycle
     amtot = shell(ishell)%am &
             +shell(jshell)%am &
             +shell(kshell)%am &
             +shell(lshell)%am

     write(*,*) ijshellpair, klshellpair
     call calculate_eri_4center_shell(basis,0.0_dp,ijshellpair,klshellpair,shellABCD)

     shell(ishell)%x0(1) = shell(ishell)%x0(1) + dx
     call calculate_eri_4center_shell(basis,0.0_dp,ijshellpair,klshellpair,shellABCD1)
     shell(ishell)%x0(1) = shell(ishell)%x0(1) - dx

     shell(jshell)%x0(1) = shell(jshell)%x0(1) + dx
     call calculate_eri_4center_shell(basis,0.0_dp,ijshellpair,klshellpair,shellABCD2)
     shell(jshell)%x0(1) = shell(jshell)%x0(1) - dx

     shell(kshell)%x0(1) = shell(kshell)%x0(1) + dx
     call calculate_eri_4center_shell(basis,0.0_dp,ijshellpair,klshellpair,shellABCD3)
     shell(kshell)%x0(1) = shell(kshell)%x0(1) - dx

     shell(lshell)%x0(1) = shell(lshell)%x0(1) + dx
     call calculate_eri_4center_shell(basis,0.0_dp,ijshellpair,klshellpair,shellABCD4)
     shell(lshell)%x0(1) = shell(lshell)%x0(1) - dx




     allocate(shell_gradA(ni,nj,nk,nl,3))
     allocate(shell_gradB(ni,nj,nk,nl,3))
     allocate(shell_gradC(ni,nj,nk,nl,3))
     allocate(shell_gradD(ni,nj,nk,nl,3))
     call calculate_eri_4center_shell_grad(basis,0.0_dp,ijshellpair,klshellpair,&
                                           shell_gradA,shell_gradB,shell_gradC,shell_gradD)
!     write(201,*) shell_gradA(:,:,:,:,1)
!     write(202,*) shell_gradB(:,:,:,:,1)
!     write(203,*) shell_gradC(:,:,:,:,1)
!     write(204,*) shell_gradD(:,:,:,:,1)

!     write(301,*) ( shellABCD1(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradA(:,:,:,:,1)
!     write(302,*) ( shellABCD2(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradB(:,:,:,:,1)
!     write(303,*) ( shellABCD3(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradC(:,:,:,:,1)
!     write(304,*) ( shellABCD4(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradD(:,:,:,:,1)

     write(*,*) 'FBFB',MAXVAL( ABS( ( shellABCD1(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradA(:,:,:,:,1) ) )
     write(*,*) 'FBFB',MAXVAL( ABS( ( shellABCD2(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradB(:,:,:,:,1) ) )
     write(*,*) 'FBFB',MAXVAL( ABS( ( shellABCD3(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradC(:,:,:,:,1) ) )
     write(*,*) 'FBFB',MAXVAL( ABS( ( shellABCD4(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradD(:,:,:,:,1) ) )
     if( MAXVAL( ABS( ( shellABCD1(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradA(:,:,:,:,1) ) ) > 0.0001 ) then
       write(*,*) '========= Problem with gradAx',amtot
       write(*,*) shell(ishell)%iatom,shell(ishell)%am
       write(*,*) shell(jshell)%iatom,shell(jshell)%am
       write(*,*) shell(kshell)%iatom,shell(kshell)%am
       write(*,*) shell(lshell)%iatom,shell(lshell)%am
       write(*,*) '============================='
     endif

     if( MAXVAL( ABS( ( shellABCD2(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradB(:,:,:,:,1) ) ) > 0.0001 ) then
       write(*,*) '========= Problem with gradBx',amtot
       write(*,*) shell(ishell)%iatom,shell(ishell)%am
       write(*,*) shell(jshell)%iatom,shell(jshell)%am
       write(*,*) shell(kshell)%iatom,shell(kshell)%am
       write(*,*) shell(lshell)%iatom,shell(lshell)%am
       write(*,*) '============================='
     endif

     if( MAXVAL( ABS( ( shellABCD3(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradC(:,:,:,:,1) ) ) > 0.0001 ) then
       write(90,*) shell(ishell)%am,shell(jshell)%am,shell(kshell)%am,shell(lshell)%am, MAXVAL( ABS( ( shellABCD3(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradC(:,:,:,:,1) ) )
!       if( shell(jshell)%am+shell(lshell)%am == 0 ) stop 'WEIRDO'
     else
       write(91,*) shell(ishell)%am,shell(jshell)%am,shell(kshell)%am,shell(lshell)%am, MAXVAL( ABS( ( shellABCD3(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradC(:,:,:,:,1) ) )
!       if( shell(jshell)%am+shell(lshell)%am /= 0 ) stop 'WEIRDODO'
     endif
     if( MAXVAL( ABS( ( shellABCD3(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradC(:,:,:,:,1) ) ) > 0.0001 ) then
       write(*,*) '========= Problem with gradCx',amtot
       write(*,*) shell(ishell)%iatom,shell(ishell)%am
       write(*,*) shell(jshell)%iatom,shell(jshell)%am
       write(*,*) shell(kshell)%iatom,shell(kshell)%am
       write(*,*) shell(lshell)%iatom,shell(lshell)%am
       write(*,*) '============================='
!       write(*,*) shell_gradC(:,:,:,:,1)
!       write(*,*) shell_gradD(:,:,:,:,1)
!       write(*,*) ( shellABCD3(:,:,:,:)-shellABCD(:,:,:,:) ) / dx
!       write(*,*) ( shellABCD4(:,:,:,:)-shellABCD(:,:,:,:) ) / dx
     else
       if( amtot >= 2 ) then
         write(*,*) '======= SURPRISED C',amtot
         write(*,*) shell(ishell)%iatom,shell(ishell)%am
         write(*,*) shell(jshell)%iatom,shell(jshell)%am
         write(*,*) shell(kshell)%iatom,shell(kshell)%am
         write(*,*) shell(lshell)%iatom,shell(lshell)%am
         write(*,*) '============================='
       endif
     endif

     if( MAXVAL( ABS( ( shellABCD4(:,:,:,:)-shellABCD(:,:,:,:) ) / dx - shell_gradD(:,:,:,:,1) ) ) > 0.0001 ) then
       write(*,*) '========= Problem with gradDx',amtot
       write(*,*) shell(ishell)%iatom,shell(ishell)%am
       write(*,*) shell(jshell)%iatom,shell(jshell)%am
       write(*,*) shell(kshell)%iatom,shell(kshell)%am
       write(*,*) shell(lshell)%iatom,shell(lshell)%am
       write(*,*) '============================='
     else
       if( amtot >= 2 ) then
         write(*,*) 'SURPRISED D',amtot
       endif
     endif




     deallocate(shell_gradA,shell_gradB,shell_gradC,shell_gradD)
     deallocate(shellABCD)
     deallocate(shellABCD1)
     deallocate(shellABCD2)
     deallocate(shellABCD3)
     deallocate(shellABCD4)

   enddo
 enddo
 call die('enough')
#endif

 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)
   nk = number_basis_function_am( basis%gaussian_type , shell(kshell)%am )
   nl = number_basis_function_am( basis%gaussian_type , shell(lshell)%am )

   do ijshellpair=1,nshellpair
     ishell = index_shellpair(1,ijshellpair)
     jshell = index_shellpair(2,ijshellpair)
     ni = number_basis_function_am( basis%gaussian_type , shell(ishell)%am )
     nj = number_basis_function_am( basis%gaussian_type , shell(jshell)%am )


     if( shell(ishell)%am + shell(jshell)%am > shell(kshell)%am + shell(lshell)%am ) cycle
     allocate(shell_gradA(ni,nj,nk,nl,3))
     allocate(shell_gradB(ni,nj,nk,nl,3))
     allocate(shell_gradC(ni,nj,nk,nl,3))
     allocate(shell_gradD(ni,nj,nk,nl,3))
     call calculate_eri_4center_shell_grad(basis,0.0_dp,ijshellpair,klshellpair,&
                                           shell_gradA,shell_gradB,shell_gradC,shell_gradD)
!     if( ishell == 1 .AND. jshell == 1 .AND. kshell == 3 .AND. lshell == 4 ) then
!       write(*,*) 'FBFB1 3444',shell_gradA(:,:,:,:,1)
!       write(*,*) 'FBFB2 3444',shell_gradB(:,:,:,:,1)
!       write(*,*) 'FBFB3 3444',shell_gradC(:,:,:,:,1)
!       write(*,*) 'FBFB4 3444',shell_gradD(:,:,:,:,1)
!       write(*,*) eri(1,1,3,6),eri(1,1,4,6),eri(1,1,5,6)
!     endif
     !
     ! Hartree
     !
     if( kshell /= lshell ) then
       fact = -4.0_dp
     else
       fact = -2.0_dp
     endif

     iatom = shell(ishell)%iatom
     do lbf=shell(lshell)%istart,shell(lshell)%iend
       do kbf=shell(kshell)%istart,shell(kshell)%iend
         do jbf=shell(jshell)%istart,shell(jshell)%iend
           do ibf=shell(ishell)%istart,shell(ishell)%iend
             force_har(:,iatom) = force_har(:,iatom)  &
                                   + fact * SUM( p_matrix(ibf,jbf,:) * p_matrix(kbf,lbf,:) ) &
                                      * shell_gradA(ibf-shell(ishell)%istart+1,       &
                                                    jbf-shell(jshell)%istart+1,       &
                                                    kbf-shell(kshell)%istart+1,       &
                                                    lbf-shell(lshell)%istart+1,:)
           enddo
         enddo
       enddo
     enddo

     if( ishell /= jshell ) then
       iatom = shell(jshell)%iatom
       do lbf=shell(lshell)%istart,shell(lshell)%iend
         do kbf=shell(kshell)%istart,shell(kshell)%iend
           do jbf=shell(jshell)%istart,shell(jshell)%iend
             do ibf=shell(ishell)%istart,shell(ishell)%iend
               force_har(:,iatom) = force_har(:,iatom)  &
                                     + fact * SUM( p_matrix(jbf,ibf,:) * p_matrix(kbf,lbf,:) ) &
                                        * shell_gradB(ibf-shell(ishell)%istart+1,       &
                                                      jbf-shell(jshell)%istart+1,       &
                                                      kbf-shell(kshell)%istart+1,       &
                                                      lbf-shell(lshell)%istart+1,:)
             enddo
           enddo
         enddo
       enddo
     endif

     !
     ! When the opposite is not calculated by LIBINT:
     if( shell(ishell)%am + shell(jshell)%am < shell(kshell)%am + shell(lshell)%am ) then
       if( ishell /= jshell ) then
         fact = -4.0_dp
       else
         fact = -2.0_dp
       endif

       iatom = shell(kshell)%iatom
       do lbf=shell(lshell)%istart,shell(lshell)%iend
         do kbf=shell(kshell)%istart,shell(kshell)%iend
           do jbf=shell(jshell)%istart,shell(jshell)%iend
             do ibf=shell(ishell)%istart,shell(ishell)%iend
               force_har(:,iatom) = force_har(:,iatom)  &
                                     + fact * SUM( p_matrix(kbf,lbf,:) * p_matrix(ibf,jbf,:) ) &
                                        * shell_gradC(ibf-shell(ishell)%istart+1,       &
                                                      jbf-shell(jshell)%istart+1,       &
                                                      kbf-shell(kshell)%istart+1,       &
                                                      lbf-shell(lshell)%istart+1,:)
             enddo
           enddo
         enddo
       enddo

       if( kshell /= lshell ) then
         iatom = shell(lshell)%iatom
         do lbf=shell(lshell)%istart,shell(lshell)%iend
           do kbf=shell(kshell)%istart,shell(kshell)%iend
             do jbf=shell(jshell)%istart,shell(jshell)%iend
               do ibf=shell(ishell)%istart,shell(ishell)%iend
                 force_har(:,iatom) = force_har(:,iatom)  &
                                       + fact * SUM( p_matrix(lbf,kbf,:) * p_matrix(ibf,jbf,:) ) &
                                          * shell_gradD(ibf-shell(ishell)%istart+1,       &
                                                        jbf-shell(jshell)%istart+1,       &
                                                        kbf-shell(kshell)%istart+1,       &
                                                        lbf-shell(lshell)%istart+1,:)
               enddo
             enddo
           enddo
         enddo
       endif


     endif

     !
     ! Exchange
     !
     iatom = shell(ishell)%iatom
     do lbf=shell(lshell)%istart,shell(lshell)%iend
       do kbf=shell(kshell)%istart,shell(kshell)%iend
         do jbf=shell(jshell)%istart,shell(jshell)%iend
           do ibf=shell(ishell)%istart,shell(ishell)%iend
             force_exx(:,iatom) = force_exx(:,iatom)  &
                                   + 2.0_dp / spin_fact * SUM( p_matrix(ibf,kbf,:) * p_matrix(jbf,lbf,:) ) &
                                      * shell_gradA(ibf-shell(ishell)%istart+1,       &
                                                    jbf-shell(jshell)%istart+1,       &
                                                    kbf-shell(kshell)%istart+1,       &
                                                    lbf-shell(lshell)%istart+1,:)
           enddo
         enddo
       enddo
     enddo

     if( kshell /= lshell ) then
       iatom = shell(ishell)%iatom
       do lbf=shell(lshell)%istart,shell(lshell)%iend
         do kbf=shell(kshell)%istart,shell(kshell)%iend
           do jbf=shell(jshell)%istart,shell(jshell)%iend
             do ibf=shell(ishell)%istart,shell(ishell)%iend
               force_exx(:,iatom) = force_exx(:,iatom)  & 
                                     + 2.0_dp / spin_fact * SUM( p_matrix(ibf,lbf,:) * p_matrix(jbf,kbf,:) ) &
                                        * shell_gradA(ibf-shell(ishell)%istart+1,       &
                                                      jbf-shell(jshell)%istart+1,       &
                                                      kbf-shell(kshell)%istart+1,       &
                                                      lbf-shell(lshell)%istart+1,:)
             enddo
           enddo
         enddo
       enddo
     endif

     if( ishell /= jshell ) then
       iatom = shell(jshell)%iatom
       do lbf=shell(lshell)%istart,shell(lshell)%iend
         do kbf=shell(kshell)%istart,shell(kshell)%iend
           do jbf=shell(jshell)%istart,shell(jshell)%iend
             do ibf=shell(ishell)%istart,shell(ishell)%iend
               force_exx(:,iatom) = force_exx(:,iatom)  & 
                                     + 2.0_dp / spin_fact * SUM( p_matrix(jbf,kbf,:) * p_matrix(ibf,lbf,:) ) &
                                        * shell_gradB(ibf-shell(ishell)%istart+1,       &
                                                      jbf-shell(jshell)%istart+1,       &
                                                      kbf-shell(kshell)%istart+1,       &
                                                      lbf-shell(lshell)%istart+1,:)
             enddo
           enddo
         enddo
       enddo

       if( kshell /= lshell ) then
         iatom = shell(jshell)%iatom
         do lbf=shell(lshell)%istart,shell(lshell)%iend
           do kbf=shell(kshell)%istart,shell(kshell)%iend
             do jbf=shell(jshell)%istart,shell(jshell)%iend
               do ibf=shell(ishell)%istart,shell(ishell)%iend
                 force_exx(:,iatom) = force_exx(:,iatom)  & 
                                       + 2.0_dp / spin_fact * SUM( p_matrix(jbf,lbf,:) * p_matrix(ibf,kbf,:) ) &
                                          * shell_gradB(ibf-shell(ishell)%istart+1,       &
                                                        jbf-shell(jshell)%istart+1,       &
                                                        kbf-shell(kshell)%istart+1,       &
                                                        lbf-shell(lshell)%istart+1,:)
               enddo
             enddo
           enddo
         enddo
       endif
     endif

     !
     ! When the opposite is not calculated by LIBINT:
     if( shell(ishell)%am + shell(jshell)%am /= shell(kshell)%am + shell(lshell)%am ) then
       iatom = shell(kshell)%iatom
       do lbf=shell(lshell)%istart,shell(lshell)%iend
         do kbf=shell(kshell)%istart,shell(kshell)%iend
           do jbf=shell(jshell)%istart,shell(jshell)%iend
             do ibf=shell(ishell)%istart,shell(ishell)%iend
               force_exx(:,iatom) = force_exx(:,iatom)  & 
                                     + 2.0_dp / spin_fact * SUM( p_matrix(jbf,lbf,:) * p_matrix(ibf,kbf,:) ) &
                                        * shell_gradC(ibf-shell(ishell)%istart+1,       &
                                                      jbf-shell(jshell)%istart+1,       &
                                                      kbf-shell(kshell)%istart+1,       &
                                                      lbf-shell(lshell)%istart+1,:)
             enddo
           enddo
         enddo
       enddo

       if( ishell /= jshell ) then
         iatom = shell(kshell)%iatom
         do lbf=shell(lshell)%istart,shell(lshell)%iend
           do kbf=shell(kshell)%istart,shell(kshell)%iend
             do jbf=shell(jshell)%istart,shell(jshell)%iend
               do ibf=shell(ishell)%istart,shell(ishell)%iend
                 force_exx(:,iatom) = force_exx(:,iatom)  & 
                                       + 2.0_dp / spin_fact * SUM( p_matrix(ibf,lbf,:) * p_matrix(jbf,kbf,:) ) &
                                          * shell_gradC(ibf-shell(ishell)%istart+1,       &
                                                        jbf-shell(jshell)%istart+1,       &
                                                        kbf-shell(kshell)%istart+1,       &
                                                        lbf-shell(lshell)%istart+1,:)
               enddo
             enddo
           enddo
         enddo
       endif

       if( kshell /= lshell ) then
         iatom = shell(lshell)%iatom
         do lbf=shell(lshell)%istart,shell(lshell)%iend
           do kbf=shell(kshell)%istart,shell(kshell)%iend
             do jbf=shell(jshell)%istart,shell(jshell)%iend
               do ibf=shell(ishell)%istart,shell(ishell)%iend
                 force_exx(:,iatom) = force_exx(:,iatom)  & 
                                       + 2.0_dp / spin_fact * SUM( p_matrix(jbf,kbf,:) * p_matrix(ibf,lbf,:) ) &
                                          * shell_gradD(ibf-shell(ishell)%istart+1,       &
                                                        jbf-shell(jshell)%istart+1,       &
                                                        kbf-shell(kshell)%istart+1,       &
                                                        lbf-shell(lshell)%istart+1,:)
               enddo
             enddo
           enddo
         enddo
         if( ishell /= jshell ) then
           iatom = shell(lshell)%iatom
           do lbf=shell(lshell)%istart,shell(lshell)%iend
             do kbf=shell(kshell)%istart,shell(kshell)%iend
               do jbf=shell(jshell)%istart,shell(jshell)%iend
                 do ibf=shell(ishell)%istart,shell(ishell)%iend
                   force_exx(:,iatom) = force_exx(:,iatom)  & 
                                         + 2.0_dp / spin_fact * SUM( p_matrix(ibf,kbf,:) * p_matrix(jbf,lbf,:) ) &
                                            * shell_gradD(ibf-shell(ishell)%istart+1,       &
                                                          jbf-shell(jshell)%istart+1,       &
                                                          kbf-shell(kshell)%istart+1,       &
                                                          lbf-shell(lshell)%istart+1,:)
                 enddo
               enddo
             enddo
           enddo

         endif

       endif

     endif


     deallocate(shell_gradA,shell_gradB,shell_gradC,shell_gradD)
   enddo
 enddo

 ! is_core is an inefficient way to get the Kinetic+Nucleus hamiltonian
 if( calc_type%is_core ) force_har(:,:) = 0.0_dp

 write(stdout,'(/,1x,a)') ' ====== Hartree Forces ====== '
 do iatom=1,natom
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force_har(:,iatom)
 enddo
 write(stdout,'(1x,a,/)') ' ==================== '

 write(stdout,'(/,1x,a)') ' ====== Exchange Forces ====== '
 do iatom=1,natom
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force_exx(:,iatom) * alpha_hybrid
 enddo
 write(stdout,'(1x,a,/)') ' ==================== '


 write(stdout,'(/,1x,a)') ' ====== Nucleus repulsion Forces ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'
 call nucleus_nucleus_force()
 do iatom=1,natom
   write(*,'(1x,a,i4,a,2x,3(2x,f16.8))') 'atom ',iatom,':',force_nuc_nuc(:,iatom)
 enddo
 write(stdout,'(1x,a,/)') ' ==================== '


 !
 ! Total forces
 !
 force(:,:) = force_nuc_nuc(:,:) & 
              + force_kin(:,:) + force_nuc(:,:) + force_har(:,:) + force_exx(:,:) * alpha_hybrid
#ifndef PSP
 force(:,:) = force(:,:) + force_ovp(:,:)
#endif

 write(stdout,'(/,1x,a)') ' ====== Total Forces ====== '
 do iatom=1,natom
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force(:,iatom)
 enddo
 write(stdout,'(/,1x,a)') ' ====== Hellman Feynman Forces ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'
 do iatom=1,natom
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force_hl(:,iatom)
 enddo
 write(stdout,'(/,1x,a)') ' ====== Pulay Forces ====== '
 write(*,'(1x,a)') 'Atoms                  Fx               Fy                 Fz'
 do iatom=1,natom
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force(:,iatom) - force_hl(:,iatom)
 enddo
 write(stdout,'(1x,a,/)') ' ==================== '

 write(stdout,'(/,1x,a)') ' ====== nuc_nuc  kin  nuc har exx Forces ====== '
 do iatom=1,natom
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force_nuc_nuc(:,iatom)
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force_kin(:,iatom)
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force_nuc(:,iatom)
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force_har(:,iatom)
   write(*,'(1x,a,i4,a,2x,3(2x,e16.8))') 'atom ',iatom,':',force_exx(:,iatom) * alpha_hybrid
 enddo
 write(stdout,'(1x,a,/)') ' ==================== '

 deallocate(p_matrix)
 deallocate(r_matrix)


contains

subroutine force_twobody_hartree_add(deriv,shell1,shell2,shell3,shell4,shell_grad,force_inout)
 implicit none
 integer,intent(in)     :: deriv
 integer,intent(in)     :: shell1,shell2,shell3,shell4
 real(dp),intent(in)    :: shell_grad(:,:,:,:,:)
 real(dp),intent(inout) :: force_inout(:,:)
!=====
 integer :: iatom
 integer :: ibf,jbf,kbf,lbf
!=====

 select case(deriv)
 case(1)
   iatom = shell(shell1)%iatom
 case(2)
   iatom = shell(shell2)%iatom
 case(3)
   iatom = shell(shell3)%iatom
 case(4)
   iatom = shell(shell4)%iatom
 end select

 do lbf=shell(shell4)%istart,shell(shell4)%iend
   do kbf=shell(shell3)%istart,shell(shell3)%iend
     do jbf=shell(shell2)%istart,shell(shell2)%iend
       do ibf=shell(shell1)%istart,shell(shell1)%iend
         force_inout(:,iatom) = force_inout(:,iatom)  &
                               - 2.0_dp * SUM( p_matrix(ibf,jbf,:) * p_matrix(kbf,lbf,:) ) &
                                  * shell_grad(ibf-shell(shell1)%istart+1,       &
                                               jbf-shell(shell2)%istart+1,       &
                                               kbf-shell(shell3)%istart+1,       &
                                               lbf-shell(shell4)%istart+1,:)
       enddo
     enddo
   enddo
 enddo

end subroutine force_twobody_hartree_add


subroutine force_twobody_exchange_add(deriv,shell1,shell2,shell3,shell4,shell_grad,force_inout)
 implicit none
 integer,intent(in)     :: deriv
 integer,intent(in)     :: shell1,shell2,shell3,shell4
 real(dp),intent(in)    :: shell_grad(:,:,:,:,:)
 real(dp),intent(inout) :: force_inout(:,:)
!=====
 integer :: iatom
 integer :: ibf,jbf,kbf,lbf
!=====

 select case(deriv)
 case(1)
   iatom = shell(shell1)%iatom
 case(2)
   iatom = shell(shell2)%iatom
 case(3)
   iatom = shell(shell3)%iatom
 case(4)
   iatom = shell(shell4)%iatom
 end select


 do lbf=shell(shell4)%istart,shell(shell4)%iend
   do kbf=shell(shell3)%istart,shell(shell3)%iend
     do jbf=shell(shell2)%istart,shell(shell2)%iend
       do ibf=shell(shell1)%istart,shell(shell1)%iend
         force_inout(:,iatom) = force_inout(:,iatom)  &
                               + 2.0_dp / spin_fact * SUM( p_matrix(ibf,kbf,:) * p_matrix(jbf,lbf,:) ) &
                                  * shell_grad(ibf-shell(shell1)%istart+1,       &
                                               jbf-shell(shell2)%istart+1,       &
                                               kbf-shell(shell3)%istart+1,       &
                                               lbf-shell(shell4)%istart+1,:)
       enddo
     enddo 
   enddo 
 enddo

end subroutine force_twobody_exchange_add


end subroutine calculate_force


!=========================================================================
