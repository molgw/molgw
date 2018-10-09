!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to evaluate the Kohn-Sham Hamiltonian
! with no distribution of the memory
!
!=========================================================================
module m_hamiltonian
 use m_definitions
 use m_timing
 use m_mpi
 use m_scalapack
 use m_warning
 use m_memory
 use m_cart_to_pure
 use m_inputparam,only: nspin,spin_fact,scalapack_block_min
 use m_basis_set
 use m_eri_calculate


contains


!=========================================================================
subroutine setup_hartree(p_matrix,hartree_ij,ehartree)
 implicit none
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: hartree_ij(:,:)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: nbf
 integer              :: ibf,jbf,kbf,lbf
 character(len=100)   :: title
!=====

 call start_clock(timing_hartree)

 write(stdout,*) 'Calculate Hartree term'

 nbf = SIZE(hartree_ij(:,:),DIM=1)

 hartree_ij(:,:)=0.0_dp

 do jbf=1,nbf
   do ibf=1,nbf
     if( negligible_basispair(ibf,jbf) ) cycle
     do lbf=1,nbf
       !
       ! symmetry k <-> l
       do kbf=1,lbf-1 ! nbf
         if( negligible_basispair(kbf,lbf) ) cycle
         !
         ! symmetry (ij|kl) = (kl|ij) has been used to loop in the fast order
         hartree_ij(ibf,jbf) = hartree_ij(ibf,jbf) &
                    + eri(kbf,lbf,ibf,jbf) * SUM( p_matrix(kbf,lbf,:) ) * 2.0_dp
       enddo
       hartree_ij(ibf,jbf) = hartree_ij(ibf,jbf) &
                  + eri(lbf,lbf,ibf,jbf) * SUM( p_matrix(lbf,lbf,:) )
     enddo
   enddo
 enddo


 title='=== Hartree contribution ==='
 call dump_out_matrix(.FALSE.,title,nbf,1,hartree_ij)

 ehartree = 0.5_dp*SUM(hartree_ij(:,:)*p_matrix(:,:,1))
 if( nspin == 2 ) then
   ehartree = ehartree + 0.5_dp*SUM(hartree_ij(:,:)*p_matrix(:,:,2))
 endif

 call stop_clock(timing_hartree)

end subroutine setup_hartree


!=========================================================================
subroutine setup_hartree_oneshell(basis,p_matrix,hartree_ij,ehartree)
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: p_matrix(:,:,:)
 real(dp),intent(out)       :: hartree_ij(:,:)
 real(dp),intent(out)       :: ehartree
!=====
 real(dp),parameter      :: TOL_DENSITY_MATRIX=1.0e-7_dp
 integer                 :: ijshellpair,klshellpair
 integer                 :: ibf,jbf,kbf,lbf
 integer                 :: ishell,jshell,kshell,lshell
 integer                 :: ni,nj,nk,nl
 real(dp)                :: fact
 real(dp),allocatable    :: shell_ijkl(:,:,:,:)
 logical,allocatable     :: skip_shellpair(:)
 real(dp)                :: cost(nshellpair)
 real(dp)                :: cost2(nshellpair)
 real(dp)                :: load(nproc_world)
 integer                 :: shellpair_cpu(nshellpair)
 logical                 :: mask(nshellpair)
!=====


 call start_clock(timing_hartree)

 write(stdout,*) 'Calculate Hartree term with out-of-core integrals'

 !
 ! Phenomenological formula to evaluate the CPU time for each shell pair
 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)
   nk = number_basis_function_am('CART', basis%shell(kshell)%am )
   nl = number_basis_function_am('CART', basis%shell(lshell)%am )
   cost(klshellpair) =   5.0e-6_dp * ( nk * nl )**1.6   &
                       + 9.0e-4_dp * ( basis%shell(kshell)%ng * basis%shell(lshell)%ng )
 enddo
 do klshellpair=1,nshellpair
   cost2(klshellpair) = SUM(cost(1:klshellpair)) * cost(klshellpair)
 enddo

 !
 ! Distribute workload among CPUs
 mask(:) = .TRUE.
 load(:) = 0.0_dp
 do klshellpair=1,nshellpair
   ijshellpair = MAXLOC(cost2(:), MASK=mask(:), DIM=1)
   mask(ijshellpair) = .FALSE.
   shellpair_cpu(ijshellpair) = MINLOC(load(:), DIM=1)
   load(shellpair_cpu(ijshellpair)) = load(shellpair_cpu(ijshellpair)) +  cost2(ijshellpair)
 enddo

 !
 ! Filter out the low density matrix shells
 !
 allocate(skip_shellpair(nshellpair))
 skip_shellpair(:) = .TRUE.
 do ijshellpair=1,nshellpair
   ishell = index_shellpair(1,ijshellpair)
   jshell = index_shellpair(2,ijshellpair)
   do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
     do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend
       if( ANY( ABS(p_matrix(ibf,jbf,:) / spin_fact) > TOL_DENSITY_MATRIX ) ) then
         skip_shellpair(ijshellpair) = .FALSE.
       endif
     enddo
   enddo
 enddo
 write(stdout,'(1x,a,i6,a,i6)') 'Shell pair skipped due to low density matrix screening:',COUNT( skip_shellpair(:) ),' / ',nshellpair


 hartree_ij(:,:) = 0.0_dp
 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)
   nk = number_basis_function_am( basis%gaussian_type , basis%shell(kshell)%am )
   nl = number_basis_function_am( basis%gaussian_type , basis%shell(lshell)%am )

   !if( MODULO(klshellpair,nproc_world) /= rank_world ) cycle
   if( shellpair_cpu(klshellpair) - 1 /= rank_world ) cycle

   if( skip_shellpair(klshellpair) ) cycle

   do ijshellpair=1,klshellpair
     ishell = index_shellpair(1,ijshellpair)
     jshell = index_shellpair(2,ijshellpair)
     ni = number_basis_function_am( basis%gaussian_type , basis%shell(ishell)%am )
     nj = number_basis_function_am( basis%gaussian_type , basis%shell(jshell)%am )

     ! I don't know if it is a good idea
     if( skip_shellpair(ijshellpair) ) cycle


     ! Libint ordering is enforced by the shell pair ordering
     call calculate_eri_4center_shell(basis,0.0_dp,ijshellpair,klshellpair,shell_ijkl)

     fact = 0.50_dp
     if( ishell /= jshell ) fact = fact * 2.0_dp
     if( kshell /= lshell ) fact = fact * 2.0_dp


     do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
       do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
         do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
           do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend

             hartree_ij(ibf,jbf) = hartree_ij(ibf,jbf)  &
                                   + SUM(p_matrix(kbf,lbf,:)) * fact                       &
                                      * shell_ijkl(ibf-basis%shell(ishell)%istart+1,       &
                                                   jbf-basis%shell(jshell)%istart+1,       &
                                                   kbf-basis%shell(kshell)%istart+1,       &
                                                   lbf-basis%shell(lshell)%istart+1)
           enddo
         enddo
       enddo
     enddo

     if( ijshellpair /= klshellpair ) then
       do lbf=basis%shell(lshell)%istart,basis%shell(lshell)%iend
         do kbf=basis%shell(kshell)%istart,basis%shell(kshell)%iend
           do jbf=basis%shell(jshell)%istart,basis%shell(jshell)%iend
             do ibf=basis%shell(ishell)%istart,basis%shell(ishell)%iend

               hartree_ij(kbf,lbf) = hartree_ij(kbf,lbf)  &
                                     + SUM(p_matrix(ibf,jbf,:)) * fact                       &
                                        * shell_ijkl(ibf-basis%shell(ishell)%istart+1,       &
                                                     jbf-basis%shell(jshell)%istart+1,       &
                                                     kbf-basis%shell(kshell)%istart+1,       &
                                                     lbf-basis%shell(lshell)%istart+1)
             enddo
           enddo
         enddo
       enddo
     endif

     deallocate(shell_ijkl)

   enddo
 enddo


 deallocate(skip_shellpair)


 hartree_ij(:,:) = hartree_ij(:,:) + TRANSPOSE( hartree_ij(:,:) )

 call xsum_world(hartree_ij)

 call dump_out_matrix(.FALSE.,'=== Hartree contribution ===',basis%nbf,1,hartree_ij)

 ehartree = 0.5_dp*SUM(hartree_ij(:,:)*p_matrix(:,:,1))
 if( nspin == 2 ) then
   ehartree = ehartree + 0.5_dp*SUM(hartree_ij(:,:)*p_matrix(:,:,2))
 endif

 call stop_clock(timing_hartree)

end subroutine setup_hartree_oneshell


!=========================================================================
subroutine setup_hartree_versatile_ri(p_matrix,hartree_ij,ehartree)
 implicit none
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: hartree_ij(:,:)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: nauxil_local,npair_local
 integer              :: nbf
 integer              :: ibf,jbf,kbf,lbf
 integer              :: ipair,ipair_local
 real(dp),allocatable :: partial_sum(:)
 real(dp)             :: rtmp
 character(len=100)   :: title
!=====

 call start_clock(timing_hartree)

 write(stdout,*) 'Calculate Hartree term with Resolution-of-Identity: versatile version'

 nbf = SIZE(hartree_ij(:,:),DIM=1)

 nauxil_local = NUMROC(nauxil_2center,block_row,iprow_3center,first_row,nprow_3center)
 npair_local  = NUMROC(npair         ,block_col,ipcol_3center,first_col,npcol_3center)


 allocate(partial_sum(nauxil_local))

 partial_sum(:) = 0.0_dp
 !$OMP PARALLEL PRIVATE(kbf,lbf,ipair)
 !$OMP DO REDUCTION(+:partial_sum)
 do ipair_local=1,npair_local
   ipair = INDXL2G(ipair_local,block_col,ipcol_3center,first_col,npcol_3center)
   kbf = index_basis(1,ipair)
   lbf = index_basis(2,ipair)
   ! Factor 2 comes from the symmetry of p_matrix
   partial_sum(:) = partial_sum(:) + eri_3center(:,ipair_local) * SUM( p_matrix(kbf,lbf,:) ) * 2.0_dp
   if( kbf == lbf ) partial_sum(:) = partial_sum(:) - eri_3center(:,ipair_local) * SUM( p_matrix(kbf,lbf,:) )
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 !FIXME ortho parallelization is not taken into account here!
 call DGSUM2D(cntxt_3center,'R',' ',nauxil_local,1,partial_sum,nauxil_local,-1,-1)

 ! Hartree potential is not sensitive to spin
 hartree_ij(:,:) = 0.0_dp

 !$OMP PARALLEL PRIVATE(ibf,jbf,ipair,rtmp)
 !$OMP DO
 do ipair_local=1,npair_local
   ipair = INDXL2G(ipair_local,block_col,ipcol_3center,first_col,npcol_3center)

   rtmp = DOT_PRODUCT( eri_3center(:,ipair_local) , partial_sum(:) )

   ibf = index_basis(1,ipair)
   jbf = index_basis(2,ipair)
   hartree_ij(ibf,jbf) = rtmp
   hartree_ij(jbf,ibf) = rtmp
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 deallocate(partial_sum)

 !
 ! Sum up the different contribution from different procs only if needed
 !FIXME ortho parallelization is not taken into account here!
 call xsum_world(hartree_ij)


 title='=== Hartree contribution ==='
 call dump_out_matrix(.FALSE.,title,nbf,1,hartree_ij)

 ehartree = 0.5_dp * SUM(hartree_ij(:,:)*p_matrix(:,:,1))
 if( nspin == 2 ) then
   ehartree = ehartree + 0.5_dp * SUM(hartree_ij(:,:)*p_matrix(:,:,2))
 endif

 call stop_clock(timing_hartree)

end subroutine setup_hartree_versatile_ri


!=========================================================================
subroutine setup_exchange(p_matrix,exchange_ij,eexchange)
 implicit none
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: exchange_ij(:,:,:)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: nbf
 integer              :: ibf,jbf,kbf,lbf,ispin
!=====

 call start_clock(timing_exchange)

 write(stdout,*) 'Calculate Exchange term'

 nbf = SIZE(exchange_ij,DIM=1)

 exchange_ij(:,:,:)=0.0_dp

 do ispin=1,nspin
   do jbf=1,nbf
     do lbf=1,nbf
       if( negligible_basispair(lbf,jbf) ) cycle
       do kbf=1,nbf
!         if( ABS(p_matrix(kbf,lbf,ispin)) <  1.0e-12_dp ) cycle
         do ibf=1,nbf
           if( negligible_basispair(ibf,kbf) ) cycle
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           exchange_ij(ibf,jbf,ispin) = exchange_ij(ibf,jbf,ispin) &
                      - eri(ibf,kbf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact
         enddo
       enddo
     enddo
   enddo
 enddo


 eexchange = 0.5_dp*SUM(exchange_ij(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange


!=========================================================================
subroutine setup_exchange_versatile_ri(occupation,c_matrix,p_matrix,exchange_ij,eexchange)
 implicit none
 real(dp),intent(in)  :: occupation(:,:)
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: exchange_ij(:,:,:)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: nauxil_local,npair_local
 integer              :: nbf,nstate
 integer              :: ibf,jbf,ispin,istate
 integer              :: nocc
 real(dp),allocatable :: tmp(:,:)
 integer              :: ipair,ipair_local
 integer              :: ibf_auxil_first,nbf_auxil_core
!=====

 call start_clock(timing_exchange)

 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity: versatile version'

 exchange_ij(:,:,:) = 0.0_dp

 nbf    = SIZE(exchange_ij,DIM=1)
 nstate = SIZE(occupation(:,:),DIM=1)

 nauxil_local = NUMROC(nauxil_2center,block_row,iprow_3center,first_row,nprow_3center)
 npair_local  = NUMROC(npair         ,block_col,ipcol_3center,first_col,npcol_3center)

 ! Prepare a sub-division of tasks after the reduction DGSUM2D
 nbf_auxil_core  = CEILING( REAL(nauxil_local,dp) / REAL(npcol_3center,dp) )
 ibf_auxil_first = ipcol_3center * nbf_auxil_core + 1
 ! Be careful when a core has fewer data or no data at all
 if( ibf_auxil_first + nbf_auxil_core - 1 > nauxil_local ) nbf_auxil_core = nauxil_local - ibf_auxil_first + 1
 if( ibf_auxil_first > nauxil_local ) nbf_auxil_core = 0


 allocate(tmp(nauxil_local,nbf))

 ! Find highest occupied state
 nocc = get_number_occupied_states(occupation)

 do ispin=1,nspin

   do istate=1,nocc
     if( MODULO( istate - 1 , nproc_ortho ) /= rank_ortho ) cycle

     if( ABS(occupation(istate,ispin)) < completely_empty ) cycle

     tmp(:,:) = 0.0_dp
     !$OMP PARALLEL PRIVATE(ibf,jbf,ipair)
     !$OMP DO REDUCTION(+:tmp)
     do ipair_local=1,npair_local
       ipair = INDXL2G(ipair_local,block_col,ipcol_3center,first_col,npcol_3center)
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix(jbf,istate,ispin) * eri_3center(:,ipair_local)
       if( ibf /= jbf) &
         tmp(:,jbf) = tmp(:,jbf) + c_matrix(ibf,istate,ispin) * eri_3center(:,ipair_local)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     call DGSUM2D(cntxt_3center,'R',' ',nauxil_local,nbf,tmp,nauxil_local,-1,-1)

     ! exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
     !                    - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
     ! C = A^T * A + C
     if( nbf_auxil_core > 0 ) &
       call DSYRK('L','T',nbf,nbf_auxil_core,-occupation(istate,ispin)/spin_fact,tmp(ibf_auxil_first,1),nauxil_local,1.0_dp,exchange_ij(:,:,ispin),nbf)

   enddo

   !
   ! Need to symmetrize exchange_ij
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ij(ibf,jbf,ispin) = exchange_ij(jbf,ibf,ispin)
     enddo
   enddo

 enddo
 deallocate(tmp)

 call xsum_world(exchange_ij)

 eexchange = 0.5_dp * SUM( exchange_ij(:,:,:) * p_matrix(:,:,:) )

 call stop_clock(timing_exchange)

end subroutine setup_exchange_versatile_ri


!=========================================================================
subroutine setup_exchange_versatile_longrange_ri(occupation,c_matrix,p_matrix,exchange_ij,eexchange)
 implicit none
 real(dp),intent(in)  :: occupation(:,:)
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: exchange_ij(:,:,:)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: nauxil_local,npair_local
 integer              :: nbf,nstate
 integer              :: ibf,jbf,ispin,istate
 integer              :: nocc
 real(dp),allocatable :: tmp(:,:)
 integer              :: ipair,ipair_local
 integer              :: ibf_auxil_first,nbf_auxil_core
!=====

 call start_clock(timing_exchange)

 write(stdout,*) 'Calculate LR Exchange term with Resolution-of-Identity: versatile version'

 exchange_ij(:,:,:) = 0.0_dp

 nbf    = SIZE(exchange_ij,DIM=1)
 nstate = SIZE(occupation(:,:),DIM=1)

 nauxil_local = NUMROC(nauxil_2center_lr,block_row,iprow_3center,first_row,nprow_3center)
 npair_local  = NUMROC(npair            ,block_col,ipcol_3center,first_col,npcol_3center)

 ! Prepare a sub-division of tasks after the reduction DGSUM2D
 nbf_auxil_core  = CEILING( REAL(nauxil_local,dp) / REAL(npcol_3center,dp) )
 ibf_auxil_first = ipcol_3center * nbf_auxil_core + 1
 ! Be careful when a core has fewer data or no data at all
 if( ibf_auxil_first + nbf_auxil_core - 1 > nauxil_local ) nbf_auxil_core = nauxil_local - ibf_auxil_first + 1
 if( ibf_auxil_first > nauxil_local ) nbf_auxil_core = 0


 allocate(tmp(nauxil_local,nbf))

 ! Find highest occupied state
 nocc = get_number_occupied_states(occupation)

 do ispin=1,nspin

   do istate=1,nocc
     if( MODULO( istate - 1 , nproc_ortho ) /= rank_ortho ) cycle

     if( ABS(occupation(istate,ispin)) < completely_empty ) cycle

     tmp(:,:) = 0.0_dp
     !$OMP PARALLEL PRIVATE(ibf,jbf,ipair)
     !$OMP DO REDUCTION(+:tmp)
     do ipair_local=1,npair_local
       ipair = INDXL2G(ipair_local,block_col,ipcol_3center,first_col,npcol_3center)
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix(jbf,istate,ispin) * eri_3center_lr(:,ipair_local)
       if( ibf /= jbf) &
         tmp(:,jbf) = tmp(:,jbf) + c_matrix(ibf,istate,ispin) * eri_3center_lr(:,ipair_local)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     call DGSUM2D(cntxt_3center,'R',' ',nauxil_local,nbf,tmp,nauxil_local,-1,-1)

     ! exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
     !                    - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
     ! C = A^T * A + C
     if( nbf_auxil_core > 0 ) &
       call DSYRK('L','T',nbf,nbf_auxil_core,-occupation(istate,ispin)/spin_fact,tmp(ibf_auxil_first,1),nauxil_local,1.0_dp,exchange_ij(:,:,ispin),nbf)

   enddo

   !
   ! Need to symmetrize exchange_ij
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ij(ibf,jbf,ispin) = exchange_ij(jbf,ibf,ispin)
     enddo
   enddo

 enddo
 deallocate(tmp)

 call xsum_world(exchange_ij)

 eexchange = 0.5_dp * SUM( exchange_ij(:,:,:) * p_matrix(:,:,:) )

 call stop_clock(timing_exchange)

end subroutine setup_exchange_versatile_longrange_ri


!=========================================================================
subroutine setup_exchange_longrange(p_matrix,exchange_ij,eexchange)
 implicit none
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: exchange_ij(:,:,:)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: nbf
 integer              :: ibf,jbf,kbf,lbf,ispin
!=====

 call start_clock(timing_exchange)

 write(stdout,*) 'Calculate Long-Range Exchange term'

 nbf = SIZE(exchange_ij,DIM=1)

 exchange_ij(:,:,:)=0.0_dp

 do ispin=1,nspin
   do jbf=1,nbf
     do ibf=1,nbf
       do lbf=1,nbf
         do kbf=1,nbf
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           exchange_ij(ibf,jbf,ispin) = exchange_ij(ibf,jbf,ispin) &
                      - eri_lr(kbf,ibf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact
         enddo
       enddo
     enddo
   enddo
 enddo


 eexchange = 0.5_dp * SUM(exchange_ij(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange


!=========================================================================
subroutine setup_density_matrix(c_matrix,occupation,p_matrix)
 implicit none
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: occupation(:,:)
 real(dp),intent(out) :: p_matrix(:,:,:)
!=====
 integer :: nbf,nstate,nocc
 integer :: ispin,ibf,jbf
 integer :: istate
 real(dp),allocatable :: c_matrix_sqrtocc(:,:)
!=====

 call start_clock(timing_density_matrix)

 write(stdout,'(1x,a)') 'Build density matrix'

 nbf    = SIZE(c_matrix(:,:,:),DIM=1)
 nstate = SIZE(c_matrix(:,:,:),DIM=2)

 if( ANY( occupation(:,:) < -1.0e-6_dp ) ) call die('setup_density_matrix: negative occupation number should not happen here.')
 ! Find the number of occupatied states
 nocc = get_number_occupied_states(occupation)

 allocate(c_matrix_sqrtocc(nbf,nocc))

 p_matrix(:,:,:) = 0.0_dp
 do ispin=1,nspin

   do istate=1,nocc
     c_matrix_sqrtocc(:,istate) = c_matrix(:,istate,ispin) * SQRT(occupation(istate,ispin))
   enddo

   call DSYRK('L','N',nbf,nocc,1.0d0,c_matrix_sqrtocc,nbf,0.0d0,p_matrix(1,1,ispin),nbf)

   ! Symmetrize
   do jbf=1,nbf
     do ibf=jbf+1,nbf
       p_matrix(jbf,ibf,ispin) = p_matrix(ibf,jbf,ispin)
     enddo
   enddo
 enddo

 deallocate(c_matrix_sqrtocc)

 call stop_clock(timing_density_matrix)


end subroutine setup_density_matrix


!=========================================================================
subroutine setup_energy_density_matrix(c_matrix,occupation,energy,q_matrix)
 implicit none
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: occupation(:,:)
 real(dp),intent(in)  :: energy(:,:)
 real(dp),intent(out) :: q_matrix(:,:)
!=====
 integer :: nbf,nstate
 integer :: ispin,ibf,jbf
 integer :: istate
!=====

 call start_clock(timing_density_matrix)

 write(stdout,'(1x,a)') 'Build energy-density matrix'

 nbf    = SIZE(c_matrix(:,:,:),DIM=1)
 nstate = SIZE(c_matrix(:,:,:),DIM=2)

 q_matrix(:,:) = 0.0_dp
 do ispin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     call DSYR('L',nbf,occupation(istate,ispin)*energy(istate,ispin),c_matrix(:,istate,ispin),1,q_matrix(:,:),nbf)
   enddo
 enddo


 ! Symmetrize
 do jbf=1,nbf
   do ibf=jbf+1,nbf
     q_matrix(jbf,ibf) = q_matrix(ibf,jbf)
   enddo
 enddo

 call stop_clock(timing_density_matrix)


end subroutine setup_energy_density_matrix


!=========================================================================
subroutine test_density_matrix(p_matrix,s_matrix)
 implicit none
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(in)  :: s_matrix(:,:)
!=====
 integer              :: nbf
 integer              :: ispin
 real(dp),allocatable :: matrix(:,:)
 character(len=100)   :: title
!=====

 nbf = SIZE(p_matrix(:,:,:),DIM=1)
 allocate(matrix(nbf,nbf))

 write(stdout,*) 'Check equality PSP = P'
 write(stdout,*) ' valid only for integer occupation numbers'

 do ispin=1,nspin

   !
   ! Calculate PSP
   matrix(:,:) = MATMUL( p_matrix(:,:,ispin), MATMUL( s_matrix(:,:) , p_matrix(:,:,ispin) ) )


   title='=== PSP ==='
   call dump_out_matrix(1,title,nbf,1,matrix)
   title='===  P  ==='
   call dump_out_matrix(1,title,nbf,1,p_matrix(:,:,ispin))

 enddo

 deallocate(matrix)

end subroutine test_density_matrix


!=========================================================================
subroutine set_occupation(nstate,temperature,electrons,magnetization,energy,occupation)
 implicit none
 integer,intent(in)   :: nstate
 real(dp),intent(in)  :: electrons,magnetization,temperature
 real(dp),intent(in)  :: energy(nstate,nspin)
 real(dp),intent(out) :: occupation(nstate,nspin)
!=====
 real(dp)             :: remaining_electrons(nspin)
 real(dp)             :: electrons_mu,mu,delta_mu,grad_electrons
 real(dp)             :: mu_change
 integer              :: istate,nlines,ilines
 logical              :: file_exists
 integer              :: occfile
 integer              :: iter
!=====

 if( temperature < 1.0e-8_dp ) then

   occupation(:,:)=0.0_dp

   inquire(file='manual_occupations',exist=file_exists)

   if(.NOT. file_exists) then
     remaining_electrons(1) = (electrons+magnetization) / REAL(nspin,dp)
     if(nspin==2) remaining_electrons(2) = (electrons-magnetization) / REAL(nspin,dp)

     do istate=1,nstate
       occupation(istate,:) = MIN(remaining_electrons(:), spin_fact)
       remaining_electrons(:)  = remaining_electrons(:) - occupation(istate,:)
     enddo
   else
     write(stdout,*)
     write(stdout,*) 'occupations are read from file: manual_occupations'
     msg='reading manual occupations from file'
     call issue_warning(msg)
     open(newunit=occfile,file='manual_occupations',status='old')
     !
     ! read nlines, all other occupations are set to zero
     read(occfile,*) nlines
     do ilines=1,nlines
       read(occfile,*) occupation(ilines,:)
     enddo
     close(occfile)
     write(stdout,*) 'occupations set, closing file'
   endif

 else

   !
   ! Finite temperature case
   !
   write(stdout,'(1x,a,f12.6,3x,f15.3)') 'Find new the occupations and Fermi level for temperature (Ha) (K): ', &
                                         temperature,temperature * Ha_K

   ! First, set mu half way between the HOMO and the LUMO
   mu = 0.50_dp * ( energy(NINT(electrons/2.0_dp)+1,1) + energy(NINT(electrons/2.0_dp),1) )

   delta_mu = 1.0e-5_dp
   electrons_mu = -1.0_dp
   iter = 0
   mu_change = 0.0_dp

   do while( ABS( electrons_mu - electrons ) > 1.0e-8_dp .AND. iter <= 100 )

     iter = iter + 1
     mu = mu + mu_change

     occupation(:,:) = fermi_dirac(energy,mu)
     electrons_mu    = SUM( occupation(:,:) )

     grad_electrons = ( SUM( fermi_dirac(energy,mu+delta_mu) ) - SUM( fermi_dirac(energy,mu-delta_mu) ) ) / ( 2.0_dp* delta_mu )

     ! Maximum change is made bounded within +/- 0.10 Hartree
     mu_change = -( electrons_mu - electrons ) / grad_electrons
     mu_change = MAX( MIN( mu_change , 0.10_dp / REAL(iter) ), -0.10_dp / REAL(iter) )

!     write(*,*) iter,mu,mu_change,0.10_dp / REAL(iter),electrons_mu

   enddo

   write(stdout,'(1x,a,f12.6)') 'Fermi level (eV): ', mu * Ha_eV

 endif

 !
 ! final check
 if( ABS( SUM(occupation(:,:)) - electrons ) > 1.0e-4_dp ) then
   write(stdout,*) 'occupation set up failed to give the right number of electrons'
   write(stdout,*) 'sum of occupations',SUM(occupation(:,:))
   write(stdout,*) 'electrons',electrons
   do istate=1,nstate
     write(stdout,*) istate,occupation(istate,:)
   enddo
   call die('FAILURE in set_occupation')
 endif

 call dump_out_occupation('=== Occupations ===',nstate,nspin,occupation)

contains

function fermi_dirac(energy,mu)
 implicit none
 real(dp),intent(in) :: energy(nstate,nspin)
 real(dp),intent(in) :: mu
 real(dp)            :: fermi_dirac(nstate,nspin)
!=====

 fermi_dirac(:,:) = spin_fact / ( 1.0_dp + EXP( ( energy(:,:) - mu ) / temperature ) )

end function fermi_dirac

end subroutine set_occupation


!=========================================================================
pure function get_number_occupied_states(occupation) result(nocc)
 implicit none

 real(dp),intent(in) :: occupation(:,:)
 integer             :: nocc
!=====
 integer :: nstate,istate,ispin,nspin_local
!=====

 nstate      = SIZE(occupation(:,:),DIM=1)
 nspin_local = SIZE(occupation(:,:),DIM=2)

 ! Find highest occupied state
 ! Take care of negative occupations, this can happen if C comes from P^{1/2}
 nocc = 0
 do ispin=1,nspin_local
   do istate=1,nstate
     if( ABS(occupation(istate,ispin)) < completely_empty )  cycle
     nocc = MAX(nocc,istate)
   enddo
 enddo


end function get_number_occupied_states


!=========================================================================
subroutine matrix_basis_to_eigen(c_matrix,matrix_inout)
 implicit none
 real(dp),intent(in)     :: c_matrix(:,:,:)
 real(dp),intent(inout)  :: matrix_inout(:,:,:)
!=====
 integer                 :: nbf,nstate
 integer                 :: ispin
 real(dp),allocatable    :: matrix_tmp(:,:)
!=====

 nbf    = SIZE(c_matrix(:,:,:),DIM=1)
 nstate = SIZE(c_matrix(:,:,:),DIM=2)


 allocate(matrix_tmp(nbf,nstate))

 do ispin=1,nspin
   !matrix_inout(1:nstate,1:nstate,ispin) = MATMUL( TRANSPOSE( c_matrix(:,:,ispin) ) , MATMUL( matrix_inout(:,:,ispin) , c_matrix(:,:,ispin) ) )
   ! H * C
   call DGEMM('N','N',nbf,nstate,nbf,1.0d0,matrix_inout(:,:,ispin),nbf, &
                                           c_matrix(:,:,ispin),nbf,     &
                                     0.0d0,matrix_tmp,nbf)
   ! C**T * (H * C)
   call DGEMM('T','N',nstate,nstate,nbf,1.0d0,c_matrix(:,:,ispin),nbf,  &
                                              matrix_tmp,nbf,           &
                                        0.0d0,matrix_inout,nbf)

 enddo
 deallocate(matrix_tmp)

end subroutine matrix_basis_to_eigen


!=========================================================================
subroutine evaluate_s2_operator(occupation,c_matrix,s_matrix)
 implicit none
 real(dp),intent(in)     :: occupation(:,:)
 real(dp),intent(in)     :: c_matrix(:,:,:)
 real(dp),intent(in)     :: s_matrix(:,:)
!=====
 integer                 :: nstate
 integer                 :: istate,jstate
 real(dp)                :: s2,s2_exact
 real(dp)                :: n1,n2,nmax,nmin
!=====

 if(nspin /= 2) return

 nstate = SIZE(occupation,DIM=1)

 n1 = SUM( occupation(:,1) )  ! Number of spin up   electrons
 n2 = SUM( occupation(:,2) )  ! Number of spin down electrons
 nmax = MAX(n1,n2)
 nmin = MIN(n1,n2)

 s2_exact = (nmax-nmin)/2.0_dp * ( (nmax-nmin)/2.0_dp + 1.0_dp )
 s2       = s2_exact + nmin
 do istate=1,nstate
   if( occupation(istate,1) < completely_empty ) cycle
   do jstate=1,nstate
     if( occupation(jstate,2) < completely_empty ) cycle

     s2 = s2 - ABS( DOT_PRODUCT( c_matrix(:,istate,1) , MATMUL( s_matrix(:,:) , c_matrix(:,jstate,2) ) )  &
                      * occupation(istate,1) * occupation(jstate,2) )**2

   enddo
 enddo


 write(stdout,'(/,a,f8.4)') ' Total Spin S**2: ',s2
 write(stdout,'(a,f8.4)')   ' Instead of:      ',s2_exact


end subroutine evaluate_s2_operator


!=========================================================================
subroutine level_shifting_up(s_matrix,c_matrix,occupation,level_shifting_energy,hamiltonian)
 implicit none
 real(dp),intent(in)    :: s_matrix(:,:)
 real(dp),intent(in)    :: c_matrix(:,:,:)
 real(dp),intent(in)    :: occupation(:,:)
 real(dp),intent(in)    :: level_shifting_energy
 real(dp),intent(inout) :: hamiltonian(:,:,:)
!=====
 integer              :: nstate
 integer              :: ispin,istate
 real(dp),allocatable :: sqrt_level_shifting(:)
 real(dp),allocatable :: matrix_tmp(:,:)
!=====

 write(stdout,'(/,a)')     ' Level shifting switched on'
 write(stdout,'(a,f12.6)') '   energy shift (eV):',level_shifting_energy * Ha_eV

 if( level_shifting_energy < 0.0_dp ) then
   call die('level_shifting_energy has to be positive!')
 endif

 nstate = SIZE(occupation,DIM=1)
 allocate(matrix_tmp,MOLD=s_matrix)
 allocate(sqrt_level_shifting(nstate))

 do ispin=1,nspin
   !
   ! Shift up empty states only
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) then
       sqrt_level_shifting(istate) = SQRT( level_shifting_energy )
     else
       sqrt_level_shifting(istate) = 0.0_dp
     endif
   enddo
   forall(istate=1:nstate)
     matrix_tmp(:,istate) =  c_matrix(:,istate,ispin) * sqrt_level_shifting(istate)
   end forall

   !
   ! M = C * E * tC
   matrix_tmp(:,:) = MATMUL( matrix_tmp(:,1:nstate) , TRANSPOSE(matrix_tmp(:,1:nstate)) )
   ! M = S * M * S
   matrix_tmp(:,:) = MATMUL( s_matrix , MATMUL( matrix_tmp , s_matrix ) )

   ! Finally update the total hamiltonian
   hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) + matrix_tmp(:,:)

 enddo

 deallocate(matrix_tmp)
 deallocate(sqrt_level_shifting)


end subroutine level_shifting_up


!=========================================================================
subroutine level_shifting_down(s_matrix,c_matrix,occupation,level_shifting_energy,energy,hamiltonian)
 implicit none
 real(dp),intent(in)    :: s_matrix(:,:)
 real(dp),intent(in)    :: c_matrix(:,:,:)
 real(dp),intent(in)    :: occupation(:,:)
 real(dp),intent(inout) :: energy(:,:)
 real(dp),intent(in)    :: level_shifting_energy
 real(dp),intent(inout) :: hamiltonian(:,:,:)
!=====
 integer              :: nstate
 integer              :: ispin,istate
 real(dp),allocatable :: sqrt_level_shifting(:)
 real(dp),allocatable :: matrix_tmp(:,:)
!=====

 if( level_shifting_energy < 0.0_dp ) then
   call die('level_shifting_energy has to be positive!')
 endif

 nstate = SIZE(occupation,DIM=1)
 allocate(matrix_tmp,MOLD=s_matrix)
 allocate(sqrt_level_shifting(nstate))

 !
 ! Shift down the energies of the virtual orbitals
 do ispin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) then
       energy(istate,ispin) = energy(istate,ispin) - level_shifting_energy
     endif
   enddo
 enddo


 do ispin=1,nspin
   !
   ! Shift down empty states only
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) then
       sqrt_level_shifting(istate) = SQRT( level_shifting_energy )
     else
       sqrt_level_shifting(istate) = 0.0_dp
     endif
   enddo
   forall(istate=1:nstate)
     matrix_tmp(:,istate) =  c_matrix(:,istate,ispin) * sqrt_level_shifting(istate)
   end forall

   !
   ! M = C * E * tC
   matrix_tmp(:,:) = MATMUL( matrix_tmp(:,1:nstate) , TRANSPOSE(matrix_tmp(:,1:nstate)) )
   ! M = S * M * S
   matrix_tmp(:,:) = MATMUL( s_matrix , MATMUL( matrix_tmp , s_matrix ) )

   ! Finally update the total hamiltonian
   hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) - matrix_tmp(:,:)

 enddo

 deallocate(matrix_tmp)
 deallocate(sqrt_level_shifting)

end subroutine level_shifting_down


!=========================================================================
subroutine setup_sqrt_overlap(TOL_OVERLAP,s_matrix,nstate,s_matrix_sqrt_inv)
 use m_tools
 implicit none

 real(dp),intent(in)                :: TOL_OVERLAP
 real(dp),intent(in)                :: s_matrix(:,:)
 integer,intent(out)                :: nstate
 real(dp),allocatable,intent(inout) :: s_matrix_sqrt_inv(:,:)
!=====
 integer  :: nbf
 integer  :: istate,jbf
 real(dp),allocatable :: s_eigval(:)
 real(dp),allocatable :: matrix_tmp(:,:)
!=====

 write(stdout,'(/,a)') ' Calculate overlap matrix square-root S^{1/2}'

 nbf = SIZE(s_matrix,DIM=1)

 allocate(matrix_tmp(nbf,nbf))
 allocate(s_eigval(nbf))

 matrix_tmp(:,:) = s_matrix(:,:)
 ! Diagonalization with or without SCALAPACK
 call diagonalize_scalapack(scalapack_block_min,nbf,matrix_tmp,s_eigval)

 nstate = COUNT( s_eigval(:) > TOL_OVERLAP )

 call clean_allocate('Overlap sqrt S^{-1/2}',s_matrix_sqrt_inv,nbf,nstate)

 write(stdout,'(/,a)')       ' Filtering basis functions that induce overcompleteness'
 write(stdout,'(a,es9.2)')   '   Lowest S eigenvalue is           ',MINVAL( s_eigval(:) )
 write(stdout,'(a,es9.2)')   '   Tolerance on overlap eigenvalues ',TOL_OVERLAP
 write(stdout,'(a,i5,a,i5)') '   Retaining ',nstate,' among ',nbf

 istate = 0
 do jbf=1,nbf
   if( s_eigval(jbf) > TOL_OVERLAP ) then
     istate = istate + 1
     s_matrix_sqrt_inv(:,istate) = matrix_tmp(:,jbf) / SQRT( s_eigval(jbf) )
   endif
 enddo

 deallocate(matrix_tmp,s_eigval)

end subroutine setup_sqrt_overlap


!=========================================================================
subroutine setup_sqrt_density_matrix(p_matrix,p_matrix_sqrt,p_matrix_occ)
 use m_tools
 implicit none

 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: p_matrix_sqrt(:,:,:)
 real(dp),intent(out) :: p_matrix_occ(:,:)
!=====
 integer              :: nbf
 integer              :: ispin,ibf
!=====

 nbf = SIZE( p_matrix(:,:,:), DIM=1 )

 write(stdout,*) 'Calculate the square root of the density matrix'
 call start_clock(timing_sqrt_density_matrix)

 do ispin=1,nspin
   p_matrix_sqrt(:,:,ispin) = p_matrix(:,:,ispin)
   ! Diagonalization with or without SCALAPACK
   call diagonalize_scalapack(scalapack_block_min,nbf,p_matrix_sqrt(:,:,ispin),p_matrix_occ(:,ispin))
   do ibf=1,nbf
     ! this is to avoid instabilities
     if( p_matrix_occ(ibf,ispin) < 1.0e-8_dp ) then
       p_matrix_occ(ibf,ispin)    = 0.0_dp
       p_matrix_sqrt(:,ibf,ispin) = 0.0_dp
     else
       p_matrix_sqrt(:,ibf,ispin) = p_matrix_sqrt(:,ibf,ispin) * SQRT( p_matrix_occ(ibf,ispin) )
     endif
   enddo
 enddo

 call stop_clock(timing_sqrt_density_matrix)

end subroutine setup_sqrt_density_matrix


!=========================================================================
subroutine get_c_matrix_from_p_matrix(p_matrix,c_matrix,occupation)
 use m_tools
 implicit none

 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: c_matrix(:,:,:)
 real(dp),intent(out) :: occupation(:,:)
!=====
 real(dp),allocatable :: p_matrix_sqrt(:,:)
 integer              :: nbf,nstate
 integer              :: ispin
!=====

 nbf    = SIZE( p_matrix(:,:,:), DIM=1 )
 nstate = SIZE( c_matrix(:,:,:), DIM=2 )
 allocate(p_matrix_sqrt(nbf,nbf))

 write(stdout,*) 'Calculate the square root of the density matrix to obtain the C matrix'
 call start_clock(timing_sqrt_density_matrix)

 do ispin=1,nspin

   ! Minus the p_matrix so that the eigenvalues are ordered from the largest to the lowest
   p_matrix_sqrt(:,:) = -p_matrix(:,:,ispin)

   ! Diagonalization with or without SCALAPACK
   call diagonalize_scalapack(scalapack_block_min,nbf,p_matrix_sqrt,occupation(:,ispin))

   c_matrix(:,:,ispin) = p_matrix_sqrt(:,1:nstate)
   occupation(:,ispin) = -occupation(:,ispin)

 enddo

 deallocate(p_matrix_sqrt)

 call stop_clock(timing_sqrt_density_matrix)

end subroutine get_c_matrix_from_p_matrix


!=========================================================================
subroutine dft_exc_vxc_batch(batch_size,basis,occupation,c_matrix,vxc_ij,exc_xc)
 use m_inputparam
 use m_dft_grid
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
 implicit none

 integer,intent(in)         :: batch_size
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(:,:)
 real(dp),intent(in)        :: c_matrix(:,:,:)
 real(dp),intent(out)       :: vxc_ij(:,:,:)
 real(dp),intent(out)       :: exc_xc
!=====
 real(dp),parameter   :: TOL_RHO=1.0e-9_dp
 integer              :: nstate
 integer              :: ibf,jbf,ispin
 integer              :: idft_xc
 integer              :: igrid_start,igrid_end,ir,nr
 real(dp)             :: normalization(nspin)
 real(dp),allocatable :: weight_batch(:)
 real(dp),allocatable :: tmp_batch(:,:)
 real(dp),allocatable :: basis_function_r_batch(:,:)
 real(dp),allocatable :: basis_function_gradr_batch(:,:,:)
 real(dp),allocatable :: exc_batch(:)
 real(dp),allocatable :: rhor_batch(:,:)
 real(dp),allocatable :: vrho_batch(:,:)
 real(dp),allocatable :: dedd_r_batch(:,:)
 real(dp),allocatable :: grad_rhor_batch(:,:,:)
 real(dp),allocatable :: sigma_batch(:,:)
 real(dp),allocatable :: dedgd_r_batch(:,:,:)
 real(dp),allocatable :: vsigma_batch(:,:)
!=====

 exc_xc = 0.0_dp
 vxc_ij(:,:,:) = 0.0_dp
 if( ndft_xc == 0 ) return

 call start_clock(timing_dft)

 nstate = SIZE(occupation,DIM=1)

#ifdef HAVE_LIBXC

 write(stdout,*) 'Calculate DFT XC potential'
 if( batch_size /= 1 ) write(stdout,*) 'Using batches of size',batch_size


 normalization(:) = 0.0_dp

 !
 ! Loop over batches of grid points
 !
 do igrid_start=1,ngrid,batch_size
   igrid_end = MIN(ngrid,igrid_start+batch_size-1)
   nr = igrid_end - igrid_start + 1

   allocate(weight_batch(nr))
   allocate(rhor_batch(nspin,nr))
   allocate(basis_function_r_batch(basis%nbf,nr))
   allocate(exc_batch(nr))
   allocate(vrho_batch(nspin,nr))
   allocate(dedd_r_batch(nspin,nr))

   if( dft_xc_needs_gradient ) then
     allocate(basis_function_gradr_batch(basis%nbf,nr,3))
     allocate(grad_rhor_batch(nspin,nr,3))
     allocate(dedgd_r_batch(3,nr,nspin))
     allocate(sigma_batch(2*nspin-1,nr))
     allocate(vsigma_batch(2*nspin-1,nr))
   endif

   weight_batch(:) = w_grid(igrid_start:igrid_end)

   call get_basis_functions_r_batch(basis,igrid_start,nr,basis_function_r_batch)
   !
   ! Get the gradient at points r
   if( dft_xc_needs_gradient ) call get_basis_functions_gradr_batch(basis,igrid_start,nr,basis_function_gradr_batch)

   !
   ! Calculate the density at points r for spin up and spin down
   ! Calculate grad rho at points r for spin up and spin down
   if( .NOT. dft_xc_needs_gradient ) then
     call calc_density_r_batch(nspin,basis%nbf,nstate,nr,occupation,c_matrix,basis_function_r_batch,rhor_batch)

   else
     call calc_density_gradr_batch(nspin,basis%nbf,nstate,nr,occupation,c_matrix, &
                                   basis_function_r_batch,basis_function_gradr_batch,rhor_batch,grad_rhor_batch)
     do ir=1,nr
       sigma_batch(1,ir) = DOT_PRODUCT( grad_rhor_batch(1,ir,:) , grad_rhor_batch(1,ir,:) )
       if( nspin == 2 ) then
         sigma_batch(2,ir) = DOT_PRODUCT( grad_rhor_batch(1,ir,:) , grad_rhor_batch(2,ir,:) )
         sigma_batch(3,ir) = DOT_PRODUCT( grad_rhor_batch(2,ir,:) , grad_rhor_batch(2,ir,:) )
       endif
     enddo

   endif

   ! Normalization
   normalization(:) = normalization(:) + MATMUL( rhor_batch(:,:) , weight_batch(:) )

   !
   ! LIBXC calls
   !

   dedd_r_batch(:,:) = 0.0_dp
   if( dft_xc_needs_gradient ) dedgd_r_batch(:,:,:) = 0.0_dp

   do idft_xc=1,ndft_xc
     if( ABS(dft_xc_coef(idft_xc)) < 1.0e-6_dp ) cycle

     select case(xc_f90_info_family(calc_type%xc_info(idft_xc)))
     case(XC_FAMILY_LDA)
       call xc_f90_lda_exc_vxc(calc_type%xc_func(idft_xc),nr,rhor_batch(1,1),exc_batch(1),vrho_batch(1,1))

     case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
       call xc_f90_gga_exc_vxc(calc_type%xc_func(idft_xc),nr,rhor_batch(1,1),sigma_batch(1,1),exc_batch(1),vrho_batch(1,1),vsigma_batch(1,1))

       ! Remove too small densities to stabilize the computation
       ! especially useful for Becke88
       do ir=1,nr
         if( ALL( rhor_batch(:,ir) < TOL_RHO ) ) then
           exc_batch(ir)      = 0.0_dp
           vrho_batch(:,ir)   = 0.0_dp
           vsigma_batch(:,ir) = 0.0_dp
         endif
       enddo

     case default
       call die('functional is not LDA nor GGA nor hybrid')
     end select

     ! XC energy
     exc_xc = exc_xc + SUM( weight_batch(:) * exc_batch(:) * SUM(rhor_batch(:,:),DIM=1) ) * dft_xc_coef(idft_xc)

     dedd_r_batch(:,:) = dedd_r_batch(:,:) + vrho_batch(:,:) * dft_xc_coef(idft_xc)

     !
     ! Set up divergence term if needed (GGA case)
     !
     if( dft_xc_needs_gradient ) then
       if(nspin==1) then

         do ir=1,nr
           dedgd_r_batch(:,ir,1) = dedgd_r_batch(:,ir,1)  &
                      + 2.0_dp * vsigma_batch(1,ir) * grad_rhor_batch(1,ir,:) * dft_xc_coef(idft_xc)
         enddo

       else

         do ir=1,nr
           dedgd_r_batch(:,ir,1) = dedgd_r_batch(:,ir,1) &
                     + ( 2.0_dp * vsigma_batch(1,ir) * grad_rhor_batch(1,ir,:) &
                                 + vsigma_batch(2,ir) * grad_rhor_batch(2,ir,:) ) * dft_xc_coef(idft_xc)

           dedgd_r_batch(:,ir,2) = dedgd_r_batch(:,ir,2) &
                     + ( 2.0_dp * vsigma_batch(3,ir) * grad_rhor_batch(2,ir,:) &
                                 + vsigma_batch(2,ir) * grad_rhor_batch(1,ir,:) ) * dft_xc_coef(idft_xc)
         enddo

       endif

     endif

   enddo ! loop on the XC functional



   !
   ! Eventually set up the vxc term
   !

   allocate(tmp_batch(basis%nbf,nr))
   do ispin=1,nspin
     !
     ! LDA and GGA
     do ir=1,nr
       tmp_batch(:,ir) = weight_batch(ir) * dedd_r_batch(ispin,ir) * basis_function_r_batch(:,ir) * 0.50_dp
     enddo
     !
     ! GGA-only
     if( dft_xc_needs_gradient ) then
       do ir=1,nr
         tmp_batch(:,ir) = tmp_batch(:,ir) &
                          +  MATMUL( basis_function_gradr_batch(:,ir,:) , dedgd_r_batch(:,ir,ispin) * weight_batch(ir) )
       enddo
     endif

     call DSYR2K('L','N',basis%nbf,nr,1.0d0,basis_function_r_batch,basis%nbf,tmp_batch,basis%nbf,1.0d0,vxc_ij(:,:,ispin),basis%nbf)
   enddo
   deallocate(tmp_batch)



   deallocate(weight_batch)
   deallocate(basis_function_r_batch)
   deallocate(exc_batch)
   deallocate(rhor_batch)
   deallocate(vrho_batch)
   deallocate(dedd_r_batch)
   if( dft_xc_needs_gradient ) then
     deallocate(basis_function_gradr_batch)
     deallocate(grad_rhor_batch)
     deallocate(sigma_batch)
     deallocate(dedgd_r_batch)
     deallocate(vsigma_batch)
   endif

 enddo ! loop on the batches




 ! Symmetrize now
 do ispin=1,nspin
   do jbf=1,basis%nbf
     do ibf=jbf+1,basis%nbf
       vxc_ij(jbf,ibf,ispin) = vxc_ij(ibf,jbf,ispin)
     enddo
   enddo
 enddo

 !
 ! Sum up the contributions from all procs only if needed
 call xsum_grid(normalization)
 call xsum_grid(vxc_ij)
 call xsum_grid(exc_xc)


#else
 write(stdout,*) 'XC energy and potential set to zero'
 write(stdout,*) 'LIBXC is not present'
#endif

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization(:)
 write(stdout,'(a,2x,f12.6,/)')    '  DFT xc energy (Ha):',exc_xc

 call stop_clock(timing_dft)

end subroutine dft_exc_vxc_batch


!=========================================================================
subroutine dft_approximate_vhxc(basis,vhxc_ij)
 use m_dft_grid
 use m_eri_calculate
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: vhxc_ij(basis%nbf,basis%nbf)
!=====
 integer,parameter    :: BATCH_SIZE=64
 real(dp),allocatable :: weight_batch(:)
 real(dp),allocatable :: basis_function_r_batch(:,:)
 real(dp),allocatable :: exc_batch(:)
 real(dp),allocatable :: rhor_batch(:)
 real(dp),allocatable :: vrho_batch(:)
 integer              :: igrid,ibf,jbf
 integer              :: igrid_start,igrid_end
 integer              :: ir,nr
 real(dp)             :: normalization
 real(dp)             :: exc
 real(dp)             :: vhgau(basis%nbf,basis%nbf)
 integer              :: iatom,ngau
 real(dp),allocatable :: alpha(:),coeff(:)
!=====

 call start_clock(timing_approx_ham)

 vhxc_ij(:,:) = 0.0_dp

 write(stdout,'(/,a)') ' Calculate approximate HXC potential with a superposition of atomic densities'

 do iatom=1,natom
   if( rank_grid /= MODULO(iatom,nproc_grid) ) cycle

   ngau = 4
   allocate(alpha(ngau),coeff(ngau))
   call element_atomicdensity(zvalence(iatom),zatom(iatom),coeff,alpha)


   call calculate_eri_approximate_hartree(basis,basis%nbf,basis%nbf,xatom(:,iatom),ngau,coeff,alpha,vhgau)
   vhxc_ij(:,:) = vhxc_ij(:,:) + vhgau(:,:)

   deallocate(alpha,coeff)
 enddo

 write(stdout,*) 'Simple LDA functional on a coarse grid'
 !
 ! Create a temporary grid with low quality
 ! This grid is to be destroyed at the end of the present subroutine
 call init_dft_grid(basis,low,.FALSE.,.FALSE.,BATCH_SIZE)


 normalization = 0.0_dp
 exc           = 0.0_dp
 do igrid_start=1,ngrid,BATCH_SIZE
   igrid_end = MIN(ngrid,igrid_start+batch_size-1)
   nr = igrid_end - igrid_start + 1

   allocate(weight_batch(nr))
   allocate(basis_function_r_batch(basis%nbf,nr))
   allocate(exc_batch(nr))
   allocate(rhor_batch(nr))
   allocate(vrho_batch(nr))

   weight_batch(:) = w_grid(igrid_start:igrid_end)

   !
   ! Get all the functions and gradients at point rr
   call get_basis_functions_r_batch(basis,igrid_start,nr,basis_function_r_batch)

   !
   ! Calculate the density at points r
   do ir=1,nr
     igrid = igrid_start + ir - 1
     call setup_atomic_density(rr_grid(:,igrid),rhor_batch(ir))
   enddo

   !
   ! Normalization
   normalization = normalization + SUM( rhor_batch(:) * weight_batch(:) )

   call teter_lda_vxc_exc(nr,rhor_batch,vrho_batch,exc_batch)

   !
   ! XC energy
   exc = exc + SUM( exc_batch(:) * weight_batch(:) * rhor_batch(:) )

   !
   ! HXC
   ! Hopefully, vrho is always negative
   do ir=1,nr
     basis_function_r_batch(:,ir) = SQRT( -weight_batch(ir) * vrho_batch(ir) ) * basis_function_r_batch(:,ir)
   enddo
   call DSYRK('L','N',basis%nbf,nr,-1.0d0,basis_function_r_batch,basis%nbf,1.0d0,vhxc_ij,basis%nbf)


   deallocate(weight_batch)
   deallocate(basis_function_r_batch)
   deallocate(exc_batch)
   deallocate(rhor_batch)
   deallocate(vrho_batch)

 enddo ! loop on the grid point
 !
 ! Sum up the contributions from all procs only if needed
 call xsum_grid(normalization)
 call xsum_grid(exc)
 call xsum_grid(vhxc_ij)

 ! Symmetrize vhxc
 do ibf=1,basis%nbf
   do jbf=ibf+1,basis%nbf
     vhxc_ij(ibf,jbf) = vhxc_ij(jbf,ibf)
   enddo
 enddo

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization
 write(stdout,  '(a,2(2x,f12.6))') '      XC energy (Ha):',exc

 !
 ! Temporary grid destroyed
 call destroy_dft_grid()

 call stop_clock(timing_approx_ham)

end subroutine dft_approximate_vhxc


end module m_hamiltonian
!=========================================================================
