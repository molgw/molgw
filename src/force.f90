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
 use m_hamiltonian_onebody
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: energy(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
!=====
 real(dp),parameter      :: TOL_DENSITY_MATRIX=1.0e-2
 integer                 :: ijshellpair,klshellpair
 integer                 :: iatom
 integer                 :: ibf,jbf,kbf,lbf
 integer                 :: ishell,jshell,kshell,lshell
 integer                 :: ni,nj,nk,nl
 real(dp)                :: fact
 real(dp),allocatable    :: r_matrix(:,:)
 real(dp),allocatable    :: p_matrix(:,:,:)
 real(dp),allocatable    :: grad_onebody(:,:,:)
 real(dp),allocatable    :: grad_nucleus(:,:,:,:)
 real(dp),allocatable    :: shell_gradA(:,:,:,:,:)
 real(dp),allocatable    :: shell_gradB(:,:,:,:,:)
 real(dp),allocatable    :: shell_gradC(:,:,:,:,:)
 real(dp),allocatable    :: shell_gradD(:,:,:,:,:)
 logical,allocatable     :: skip_shellpair(:,:)
!=====

#ifndef HAVE_LIBINT_ONEBODY
 call issue_warning('calculate_force: impossible to calculate gradient if LIBINT does support the gradients')
 return
#endif

 call start_clock(timing_force)

 write(stdout,'(/,1x,a)') 'Calculate the forces'

 allocate(p_matrix(basis%nbf,basis%nbf,nspin))
 allocate(r_matrix(basis%nbf,basis%nbf))
 call setup_density_matrix(c_matrix,occupation,p_matrix)
 call setup_energy_density_matrix(c_matrix,occupation,energy,r_matrix)

 !
 ! Filter out the low density matrix shells
 allocate(skip_shellpair(basis%nshell,basis%nshell))
 skip_shellpair(:,:) = .TRUE.
 do ijshellpair=1,nshellpair
   ishell = index_shellpair(1,ijshellpair)
   jshell = index_shellpair(2,ijshellpair)
   do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
     do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
       if( ANY( ABS(p_matrix(ibf,jbf,:) / spin_fact) > TOL_DENSITY_MATRIX ) ) then
         skip_shellpair(ishell,jshell) = .FALSE.
         skip_shellpair(jshell,ishell) = .FALSE.
       endif
     enddo
   enddo
 enddo
 write(stdout,'(1x,a,i6,a,i6)') 'Shell pair skipped due to low density matrix screening:',COUNT( skip_shellpair(:,:) ),' / ',basis%nshell**2


 !
 ! Nucleus-nucleus repulsion forces
 ! 
 call nucleus_nucleus_force()


 !
 ! Energy density matrix forces
 ! 
 allocate(grad_onebody(basis%nbf,basis%nbf,3))
 call setup_overlap_grad(basis,grad_onebody)

 force_ovp(:,:) = 0.0_dp
 do ibf=1,basis%nbf
   iatom = basis%bff(ibf)%iatom
   force_ovp(:,iatom) = force_ovp(:,iatom) + 2.0_dp * MATMUL( r_matrix(ibf,:) , grad_onebody(:,ibf,:) )
 enddo
 deallocate(grad_onebody)


 !
 ! Kinetic energy forces
 ! 
 allocate(grad_onebody(basis%nbf,basis%nbf,3))
 call setup_kinetic_grad(basis,grad_onebody)

 force_kin(:,:) = 0.0_dp
 do ibf=1,basis%nbf
   iatom = basis%bff(ibf)%iatom
   force_kin(:,iatom) = force_kin(:,iatom) + 2.0_dp * MATMUL( SUM( p_matrix(ibf,:,:),DIM=2) , grad_onebody(ibf,:,:) )
 enddo
 deallocate(grad_onebody)


 !
 ! Nucleus energy forces
 ! 
 allocate(grad_nucleus(basis%nbf,basis%nbf,natom+1,3))
 call setup_nucleus_grad(basis,grad_nucleus)

 force_nuc(:,:) = 0.0_dp
 do ibf=1,basis%nbf
   iatom = basis%bff(ibf)%iatom
   force_nuc(:,iatom) = force_nuc(:,iatom) + 2.0_dp * MATMUL( SUM( p_matrix(ibf,:,:),DIM=2) , grad_nucleus(ibf,:,natom+1,:) )
 enddo

 force_hl(:,:) = 0.0_dp
 do iatom=1,natom
   force_hl(1,iatom) = force_hl(1,iatom) + SUM( SUM( p_matrix(:,:,:),DIM=3 ) * grad_nucleus(:,:,iatom,1) )
   force_hl(2,iatom) = force_hl(2,iatom) + SUM( SUM( p_matrix(:,:,:),DIM=3 ) * grad_nucleus(:,:,iatom,2) )
   force_hl(3,iatom) = force_hl(3,iatom) + SUM( SUM( p_matrix(:,:,:),DIM=3 ) * grad_nucleus(:,:,iatom,3) )
 enddo
 force_nuc(:,:) = force_nuc(:,:) + force_hl(:,:)


 deallocate(grad_nucleus)


 !
 ! Hartree-Fock energy forces
 ! 
 write(stdout,'(1x,a)') 'Calculate the Hartee-Fock part with 4-center integrals gradient (LIBINT)'


 force_har(:,:) = 0.0_dp
 force_exx(:,:) = 0.0_dp
 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)
   nk = number_basis_function_am( basis%gaussian_type , basis%shell(kshell)%am )
   nl = number_basis_function_am( basis%gaussian_type , basis%shell(lshell)%am )

   do ijshellpair=1,nshellpair
     ishell = index_shellpair(1,ijshellpair)
     jshell = index_shellpair(2,ijshellpair)
     ni = number_basis_function_am( basis%gaussian_type , basis%shell(ishell)%am )
     nj = number_basis_function_am( basis%gaussian_type , basis%shell(jshell)%am )

     if( skip_shellpair(ishell,jshell) .AND. skip_shellpair(kshell,lshell) &
             .AND.  skip_shellpair(ishell,kshell) .AND. skip_shellpair(jshell,lshell)  &
             .AND.  skip_shellpair(ishell,lshell) .AND. skip_shellpair(jshell,kshell) ) then
       cycle
     endif

     ! Libint ordering is strict!
     if( basis%shell(ishell)%am + basis%shell(jshell)%am > basis%shell(kshell)%am + basis%shell(lshell)%am ) cycle
     allocate(shell_gradA(ni,nj,nk,nl,3))
     allocate(shell_gradB(ni,nj,nk,nl,3))
     allocate(shell_gradC(ni,nj,nk,nl,3))
     allocate(shell_gradD(ni,nj,nk,nl,3))
     call calculate_eri_4center_shell_grad(basis,0.0_dp,ijshellpair,klshellpair,&
                                           shell_gradA,shell_gradB,shell_gradC,shell_gradD)

     !
     ! Hartree
     !
     if( .NOT. ( skip_shellpair(ishell,jshell) .AND. skip_shellpair(kshell,lshell) ) ) then
       if( kshell /= lshell ) then
         fact = -4.0_dp
       else
         fact = -2.0_dp
       endif

       iatom = basis%shell(ishell)%iatom
       do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
         do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
           do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
             do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
               force_har(:,iatom) = force_har(:,iatom)  &
                                     + fact * SUM( p_matrix(ibf,jbf,:) * p_matrix(kbf,lbf,:) ) &
                                        * shell_gradA(ibf-basis%shell(ishell)%istart+1,       &
                                                      jbf-basis%shell(jshell)%istart+1,       &
                                                      kbf-basis%shell(kshell)%istart+1,       &
                                                      lbf-basis%shell(lshell)%istart+1,:)
             enddo
           enddo
         enddo
       enddo

       if( ishell /= jshell ) then
         iatom = basis%shell(jshell)%iatom
         do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
           do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
             do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
               do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
                 force_har(:,iatom) = force_har(:,iatom)  &
                                       + fact * SUM( p_matrix(ibf,jbf,:) * p_matrix(kbf,lbf,:) ) &
                                          * shell_gradB(ibf-basis%shell(ishell)%istart+1,       &
                                                        jbf-basis%shell(jshell)%istart+1,       &
                                                        kbf-basis%shell(kshell)%istart+1,       &
                                                        lbf-basis%shell(lshell)%istart+1,:)
               enddo
             enddo
           enddo
         enddo
       endif

       !
       ! When the opposite is not calculated by LIBINT:
       if( basis%shell(ishell)%am + basis%shell(jshell)%am < basis%shell(kshell)%am + basis%shell(lshell)%am ) then
         if( ishell /= jshell ) then
           fact = -4.0_dp
         else
           fact = -2.0_dp
         endif

         iatom = basis%shell(kshell)%iatom
         do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
           do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
             do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
               do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
                 force_har(:,iatom) = force_har(:,iatom)  &
                                       + fact * SUM( p_matrix(kbf,lbf,:) * p_matrix(ibf,jbf,:) ) &
                                          * shell_gradC(ibf-basis%shell(ishell)%istart+1,       &
                                                        jbf-basis%shell(jshell)%istart+1,       &
                                                        kbf-basis%shell(kshell)%istart+1,       &
                                                        lbf-basis%shell(lshell)%istart+1,:)
               enddo
             enddo
           enddo
         enddo

         if( kshell /= lshell ) then
           iatom = basis%shell(lshell)%iatom
           do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
             do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
               do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
                 do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
                   force_har(:,iatom) = force_har(:,iatom)  &
                                         + fact * SUM( p_matrix(lbf,kbf,:) * p_matrix(ibf,jbf,:) ) &
                                            * shell_gradD(ibf-basis%shell(ishell)%istart+1,       &
                                                          jbf-basis%shell(jshell)%istart+1,       &
                                                          kbf-basis%shell(kshell)%istart+1,       &
                                                          lbf-basis%shell(lshell)%istart+1,:)
                 enddo
               enddo
             enddo
           enddo
         endif


       endif
     endif

     !
     ! Exchange
     !
     if( .NOT. ( skip_shellpair(ishell,kshell) .AND. skip_shellpair(jshell,lshell) ) &
         .OR. .NOT. ( skip_shellpair(ishell,lshell) .AND. skip_shellpair(jshell,kshell) ) ) then
       iatom = basis%shell(ishell)%iatom
       do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
         do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
           do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
             do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
               force_exx(:,iatom) = force_exx(:,iatom)  &
                                     + 2.0_dp / spin_fact * SUM( p_matrix(ibf,kbf,:) * p_matrix(jbf,lbf,:) ) &
                                        * shell_gradA(ibf-basis%shell(ishell)%istart+1,       &
                                                      jbf-basis%shell(jshell)%istart+1,       &
                                                      kbf-basis%shell(kshell)%istart+1,       &
                                                      lbf-basis%shell(lshell)%istart+1,:)
             enddo
           enddo
         enddo
       enddo

       if( kshell /= lshell ) then
         iatom = basis%shell(ishell)%iatom
         do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
           do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
             do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
               do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
                 force_exx(:,iatom) = force_exx(:,iatom)  & 
                                       + 2.0_dp / spin_fact * SUM( p_matrix(ibf,lbf,:) * p_matrix(jbf,kbf,:) ) &
                                          * shell_gradA(ibf-basis%shell(ishell)%istart+1,       &
                                                        jbf-basis%shell(jshell)%istart+1,       &
                                                        kbf-basis%shell(kshell)%istart+1,       &
                                                        lbf-basis%shell(lshell)%istart+1,:)
               enddo
             enddo
           enddo
         enddo
       endif

       if( ishell /= jshell ) then
         iatom = basis%shell(jshell)%iatom
         do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
           do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
             do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
               do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
                 force_exx(:,iatom) = force_exx(:,iatom)  & 
                                       + 2.0_dp / spin_fact * SUM( p_matrix(jbf,kbf,:) * p_matrix(ibf,lbf,:) ) &
                                          * shell_gradB(ibf-basis%shell(ishell)%istart+1,       &
                                                        jbf-basis%shell(jshell)%istart+1,       &
                                                        kbf-basis%shell(kshell)%istart+1,       &
                                                        lbf-basis%shell(lshell)%istart+1,:)
               enddo
             enddo
           enddo
         enddo

         if( kshell /= lshell ) then
           iatom = basis%shell(jshell)%iatom
           do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
             do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
               do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
                 do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
                   force_exx(:,iatom) = force_exx(:,iatom)  & 
                                         + 2.0_dp / spin_fact * SUM( p_matrix(jbf,lbf,:) * p_matrix(ibf,kbf,:) ) &
                                            * shell_gradB(ibf-basis%shell(ishell)%istart+1,       &
                                                          jbf-basis%shell(jshell)%istart+1,       &
                                                          kbf-basis%shell(kshell)%istart+1,       &
                                                          lbf-basis%shell(lshell)%istart+1,:)
                 enddo
               enddo
             enddo
           enddo
         endif
       endif

       !
       ! When the opposite is not calculated by LIBINT:
       if( basis%shell(ishell)%am + basis%shell(jshell)%am /= basis%shell(kshell)%am + basis%shell(lshell)%am ) then
         iatom = basis%shell(kshell)%iatom
         do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
           do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
             do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
               do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
                 force_exx(:,iatom) = force_exx(:,iatom)  & 
                                       + 2.0_dp / spin_fact * SUM( p_matrix(jbf,lbf,:) * p_matrix(ibf,kbf,:) ) &
                                          * shell_gradC(ibf-basis%shell(ishell)%istart+1,       &
                                                        jbf-basis%shell(jshell)%istart+1,       &
                                                        kbf-basis%shell(kshell)%istart+1,       &
                                                        lbf-basis%shell(lshell)%istart+1,:)
               enddo
             enddo
           enddo
         enddo

         if( ishell /= jshell ) then
           iatom = basis%shell(kshell)%iatom
           do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
             do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
               do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
                 do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
                   force_exx(:,iatom) = force_exx(:,iatom)  & 
                                         + 2.0_dp / spin_fact * SUM( p_matrix(ibf,lbf,:) * p_matrix(jbf,kbf,:) ) &
                                            * shell_gradC(ibf-basis%shell(ishell)%istart+1,       &
                                                          jbf-basis%shell(jshell)%istart+1,       &
                                                          kbf-basis%shell(kshell)%istart+1,       &
                                                          lbf-basis%shell(lshell)%istart+1,:)
                 enddo
               enddo
             enddo
           enddo
         endif

         if( kshell /= lshell ) then
           iatom = basis%shell(lshell)%iatom
           do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
             do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
               do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
                 do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
                   force_exx(:,iatom) = force_exx(:,iatom)  & 
                                         + 2.0_dp / spin_fact * SUM( p_matrix(jbf,kbf,:) * p_matrix(ibf,lbf,:) ) &
                                            * shell_gradD(ibf-basis%shell(ishell)%istart+1,       &
                                                          jbf-basis%shell(jshell)%istart+1,       &
                                                          kbf-basis%shell(kshell)%istart+1,       &
                                                          lbf-basis%shell(lshell)%istart+1,:)
                 enddo
               enddo
             enddo
           enddo
           if( ishell /= jshell ) then
             iatom = basis%shell(lshell)%iatom
             do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
               do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
                 do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
                   do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
                     force_exx(:,iatom) = force_exx(:,iatom)  & 
                                           + 2.0_dp / spin_fact * SUM( p_matrix(ibf,kbf,:) * p_matrix(jbf,lbf,:) ) &
                                              * shell_gradD(ibf-basis%shell(ishell)%istart+1,       &
                                                            jbf-basis%shell(jshell)%istart+1,       &
                                                            kbf-basis%shell(kshell)%istart+1,       &
                                                            lbf-basis%shell(lshell)%istart+1,:)
                   enddo
                 enddo
               enddo
             enddo

           endif

         endif

       endif
     endif


     deallocate(shell_gradA,shell_gradB,shell_gradC,shell_gradD)
   enddo
 enddo

 deallocate(skip_shellpair)

 ! is_core is an inefficient way to get the Kinetic+Nucleus hamiltonian
 if( calc_type%is_core ) force_har(:,:) = 0.0_dp



 !
 ! Total forces
 !
 force(:,:) = force_nuc_nuc(:,:) + force_ovp(:,:) &
              + force_kin(:,:) + force_nuc(:,:) + force_har(:,:) + force_exx(:,:) * alpha_hybrid


 write(stdout,'(/,1x,a)') ' ====== Hellman-Feynman Forces ====== '
 write(*,'(1x,a,22x,a,19x,a,19x,a)') 'Atoms','Fx','Fy','Fz'
 do iatom=1,natom
   write(*,'(3x,a,i4,a,2x,3(2x,f19.10))') 'H-F force atom   ',iatom,':',force_hl(:,iatom)
 enddo
 write(stdout,'(/,1x,a)') ' ====== Pulay Forces ====== '
 write(*,'(1x,a,22x,a,19x,a,19x,a)') 'Atoms','Fx','Fy','Fz'
 do iatom=1,natom
   write(*,'(3x,a,i4,a,2x,3(2x,f19.10))') 'Pulay force atom ',iatom,':',force(:,iatom) - force_hl(:,iatom)
 enddo
 write(stdout,'(/,1x,a)') ' ====== Total Forces ====== '
 write(*,'(1x,a,22x,a,19x,a,19x,a)') 'Atoms','Fx','Fy','Fz'
 do iatom=1,natom
   write(*,'(3x,a,i4,a,2x,3(2x,f19.10))') 'Total force atom ',iatom,':',force(:,iatom)
 enddo
 write(stdout,'(1x,a,/)') ' ==================== '


 deallocate(p_matrix)
 deallocate(r_matrix)

 call stop_clock(timing_force)


end subroutine calculate_force


!=========================================================================
