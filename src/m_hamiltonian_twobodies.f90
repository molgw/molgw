!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to evaluate the Kohn-Sham Hamiltonian
!
!=========================================================================
module m_hamiltonian_twobodies
 use m_definitions
 use m_timing
 use m_mpi
 use m_scalapack
 use m_warning
 use m_memory
 use m_cart_to_pure
 use m_inputparam
 use m_basis_set
 use m_eri_calculate
 use m_density_tools
 use m_dft_grid
 use m_libxc_tools
 use m_io


contains


!=========================================================================
subroutine setup_hartree(p_matrix,hartree_ao,ehartree)
 implicit none
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: hartree_ao(:,:)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: nbf
 integer              :: ibf,jbf,kbf,lbf
 character(len=100)   :: title
 integer(kind=int8)   :: iint
 integer              :: index_ij,index_kl,stride
 real(dp)             :: fact_ij,fact_kl
#if defined(_OPENMP)
 integer,external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
!=====

 call start_clock(timing_hartree)
 nbf = SIZE(hartree_ao,DIM=1)

 write(stdout,*) 'Calculate Hartree term using the 8 permutation symmetries'
 !$OMP PARALLEL PRIVATE(index_ij,index_kl,ibf,jbf,kbf,lbf,stride,fact_ij,fact_kl)
 hartree_ao(:,:) = 0.0_dp

#if defined(_OPENMP)
 stride = OMP_GET_NUM_THREADS()
 index_ij = OMP_GET_THREAD_NUM() - stride + 1
#else
 stride = 1
 index_ij = 0
#endif
 index_kl = 1

 ! SCHEDULE(static,1) should not be modified
 ! else a race competition will occur when performing index_ij = index_ij + stride
 !$OMP DO REDUCTION(+:hartree_ao) SCHEDULE(static,1)
 do iint=1,nint_4center
   index_ij = index_ij + stride
   do while( index_ij > npair )
     index_kl = index_kl + 1
     index_ij = index_kl + index_ij - npair - 1
   enddo

   ibf = index_basis(1,index_ij)
   jbf = index_basis(2,index_ij)
   kbf = index_basis(1,index_kl)
   lbf = index_basis(2,index_kl)

   if( kbf /= lbf ) then
     fact_kl = 2.0_dp
   else
     fact_kl = 1.0_dp
   endif
   if( ibf /= jbf ) then
     fact_ij = 2.0_dp
   else
     fact_ij = 1.0_dp
   endif

   hartree_ao(ibf,jbf) = hartree_ao(ibf,jbf) &
                + eri_4center(iint) * SUM( p_matrix(kbf,lbf,:) ) * fact_kl
   if( ibf /= jbf ) then
     hartree_ao(jbf,ibf) = hartree_ao(jbf,ibf) &
                  + eri_4center(iint) * SUM( p_matrix(kbf,lbf,:) ) * fact_kl
   endif
   if( index_ij /= index_kl ) then
     hartree_ao(kbf,lbf) = hartree_ao(kbf,lbf) &
                 + eri_4center(iint) * SUM( p_matrix(ibf,jbf,:) ) * fact_ij
     if( kbf /= lbf ) then
       hartree_ao(lbf,kbf) = hartree_ao(lbf,kbf) &
                    + eri_4center(iint) * SUM( p_matrix(ibf,jbf,:) ) * fact_ij
     endif
   endif

 enddo
 !$OMP END DO

 !$OMP WORKSHARE
 ehartree = 0.5_dp * SUM( hartree_ao(:,:) * SUM( p_matrix(:,:,:),DIM=3) )
 !$OMP END WORKSHARE
 !$OMP END PARALLEL


 call stop_clock(timing_hartree)

end subroutine setup_hartree


!=========================================================================
subroutine setup_hartree_oneshell(basis,p_matrix,hartree_ao,ehartree)
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: p_matrix(:,:,:)
 real(dp),intent(out)       :: hartree_ao(:,:)
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
 write(stdout,'(1x,a,i6,a,i6)') 'Shell pair skipped due to low density matrix screening:', &
                                COUNT( skip_shellpair(:) ),' / ',nshellpair


 hartree_ao(:,:) = 0.0_dp
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

             hartree_ao(ibf,jbf) = hartree_ao(ibf,jbf)  &
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

               hartree_ao(kbf,lbf) = hartree_ao(kbf,lbf)  &
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


 hartree_ao(:,:) = hartree_ao(:,:) + TRANSPOSE( hartree_ao(:,:) )

 call xsum_world(hartree_ao)

 call dump_out_matrix(.FALSE.,'=== Hartree contribution ===',hartree_ao)

 ehartree = 0.5_dp*SUM(hartree_ao(:,:)*p_matrix(:,:,1))
 if( nspin == 2 ) then
   ehartree = ehartree + 0.5_dp*SUM(hartree_ao(:,:)*p_matrix(:,:,2))
 endif

 call stop_clock(timing_hartree)

end subroutine setup_hartree_oneshell


!=========================================================================
subroutine setup_hartree_ri(p_matrix,hartree_ao,ehartree)
 implicit none
 class(*),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: hartree_ao(:,:)
 real(dp),intent(out) :: ehartree
!=====
 integer              :: nbf
 integer              :: ibf,jbf,kbf,lbf
 integer              :: ipair
 real(dp),allocatable :: x_vector(:)
 real(dp)             :: rtmp,factor
 integer              :: timing_xxdft_hartree
 real(dp),allocatable :: pmat(:)
!=====

 nbf = SIZE(hartree_ao(:,:),DIM=1)

 select type(p_matrix)
 type is(real(dp))
   timing_xxdft_hartree   = timing_hartree
 type is(complex(dp))
   timing_xxdft_hartree   = timing_tddft_hartree
 class default
   call die("setup_hartree_ri: p_matrix is neither real nor complex")
 end select

 call start_clock(timing_xxdft_hartree)

 write(stdout,*) 'Calculate Hartree term with Resolution-of-Identity'


 hartree_ao(:,:) = 0.0_dp

 allocate(pmat(npair))
 allocate(x_vector(nauxil_3center))

 select type(p_matrix)
 type is(real(dp))
   !$OMP PARALLEL PRIVATE(kbf,lbf)
   !$OMP DO
   do ipair=1,npair
     kbf = index_basis(1,ipair)
     lbf = index_basis(2,ipair)
     pmat(ipair) = SUM(p_matrix(kbf,lbf,:)) * 2.0_dp
   enddo
   !$OMP END DO
   !$OMP END PARALLEL
 type is(complex(dp))
   !$OMP PARALLEL PRIVATE(kbf,lbf)
   !$OMP DO
   do ipair=1,npair
     kbf = index_basis(1,ipair)
     lbf = index_basis(2,ipair)
     ! As all pairs contribute twice for (k,l) and (l,k) and as P is hermitian,
     ! only the real part survives
     pmat(ipair) = SUM(p_matrix(kbf,lbf,:)%re) * 2.0_dp
   enddo
   !$OMP END DO
   !$OMP END PARALLEL
 end select

 ! X_P = \sum_{\alpha \beta} P_{\alpha \beta} * ( \alpha \beta | P )
 call DGEMV('T',npair,nauxil_3center,1.0d0,eri_3center,npair,pmat,1,0.0d0,x_vector,1)
 ! v_H_{alpha beta} = \sum_P ( alpha beta | P ) * X_P
 call DGEMV('N',npair,nauxil_3center,1.0d0,eri_3center,npair,x_vector,1,0.0d0,pmat,1)
 !$OMP PARALLEL PRIVATE(kbf,lbf)
 !$OMP DO
 do ipair=1,npair
   kbf = index_basis(1,ipair)
   lbf = index_basis(2,ipair)
   hartree_ao(kbf,lbf) = pmat(ipair)
   hartree_ao(lbf,kbf) = pmat(ipair)
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 ! Do not forget that the eri_3center(ibf,ibf | P ) included a factor 0.50
 do ibf=1,nbf
   hartree_ao(ibf,ibf) = hartree_ao(ibf,ibf) * 2.0_dp
 enddo
 deallocate(x_vector,pmat)

 !
 ! Sum up the different contribution from different procs
 call xsum_world(hartree_ao)
 hartree_ao(:,:) = hartree_ao(:,:) / REAL(nproc_ortho,dp)

 call dump_out_matrix(.FALSE.,'=== Hartree contribution ===',hartree_ao)

 select type(p_matrix)
 type is(real(dp))
   ehartree = 0.5_dp * SUM( hartree_ao(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
 type is(complex(dp))
   ehartree = 0.5_dp * SUM( hartree_ao(:,:) * SUM(REAL(p_matrix(:,:,:),dp),DIM=3) )
 end select


 call stop_clock(timing_xxdft_hartree)

end subroutine setup_hartree_ri


!=========================================================================
! calculate_density_auxilbasis
!
! Find the coefficients of the density expansion in the auxiliary basis
! rho(r) = \sum_I  \phi_I(r) R_I
!
! From the resolution-of-the-identity on the Coulomb metric
! R_I = \sum_J \sum_{\alpha\beta} ( I | 1/r12 | J )^{-1} (\alpha\beta | 1/r12 | J ) P_{\alpha\beta}
!
!=========================================================================
subroutine calculate_density_auxilbasis(p_matrix,rho_coeff)
 implicit none

 class(*),intent(in)              :: p_matrix(:,:,:)
 real(dp),allocatable,intent(out) :: rho_coeff(:,:)
 !=====
 integer              :: nbf
 integer              :: ibf,jbf,kbf,lbf
 integer              :: ipair,ispin
 real(dp),allocatable :: x_vector(:)
 real(dp)             :: factor
 integer              :: timing_xxdft_rhoauxil
 real(dp),allocatable :: pmat(:)
 !=====

 select type(p_matrix)
 type is(real(dp))
   timing_xxdft_rhoauxil = timing_rhoauxil
 type is(complex(dp))
   timing_xxdft_rhoauxil = timing_tddft_rhoauxil
 class default
   call die("calculate_density_auxilbasis: p_matrix is neither real nor complex")
 end select

 call start_clock(timing_xxdft_rhoauxil)

 write(stdout,*) 'Calculate Hartree term with Resolution-of-Identity'

 allocate(rho_coeff(nauxil_2center,nspin))

 allocate(pmat(npair))
 allocate(x_vector(nauxil_3center))

 do ispin=1,nspin

   select type(p_matrix)
   type is(real(dp))
     !$OMP PARALLEL PRIVATE(kbf,lbf)
     !$OMP DO
     do ipair=1,npair
       kbf = index_basis(1,ipair)
       lbf = index_basis(2,ipair)
       pmat(ipair) = p_matrix(kbf,lbf,ispin) * 2.0_dp
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
   type is(complex(dp))
     !$OMP PARALLEL PRIVATE(kbf,lbf)
     !$OMP DO
     do ipair=1,npair
       kbf = index_basis(1,ipair)
       lbf = index_basis(2,ipair)
       ! As all pairs contribute twice for (k,l) and (l,k) and as P is hermitian,
       ! only the real part survives
       pmat(ipair) = p_matrix(kbf,lbf,ispin)%re * 2.0_dp
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
   end select

   ! X_J = \sum_{\alpha \beta} P_{\alpha \beta} * ( \alpha \beta | J )
   call DGEMV('T',npair,nauxil_3center,1.0d0,eri_3center,npair,pmat,1,0.0d0,x_vector,1)

   ! R_I = \sum_I ( I | 1/r12 | J )^{-1} * X_J
   call DGEMV('N',nauxil_2center,nauxil_3center,1.0d0,eri_2center_inv,nauxil_2center,x_vector,1,0.0d0,rho_coeff(:,ispin),1)

 enddo

 call xsum_world(rho_coeff)

 deallocate(x_vector,pmat)

 call stop_clock(timing_xxdft_rhoauxil)


end subroutine calculate_density_auxilbasis


!=========================================================================
subroutine setup_hartree_genuine_ri(p_matrix,rho_coeff,hartree_ao,ehartree)
 implicit none
 class(*),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(in)  :: rho_coeff(:,:)
 real(dp),intent(out) :: hartree_ao(:,:)
 real(dp),intent(out) :: ehartree
 !=====
 integer              :: nbf
 integer              :: ibf,jbf,kbf,lbf
 integer              :: ipair,iauxil_local,iauxil_global
 real(dp),allocatable :: x_vector(:)
 real(dp)             :: rtmp,factor
 integer              :: timing_xxdft_hartree
 real(dp),allocatable :: vh(:)
 real(dp),allocatable :: rho_coeff_local_nospin(:)
 !=====

 if( nproc_ortho > 1 ) call die('setup_hartree_genuine_ri: ortho-parallelization not coded')

 nbf = SIZE(hartree_ao(:,:),DIM=1)

 select type(p_matrix)
 type is(real(dp))
   timing_xxdft_hartree   = timing_hartree
 type is(complex(dp))
   timing_xxdft_hartree   = timing_tddft_hartree
 class default
   call die("setup_hartree_ri: p_matrix is neither real nor complex")
 end select

 call start_clock(timing_xxdft_hartree)

 write(stdout,*) 'Calculate Hartree term with Resolution-of-Identity-expanded density'

 hartree_ao(:,:) = 0.0_dp

 allocate(vh(npair))
 allocate(rho_coeff_local_nospin(nauxil_3center))

 do iauxil_local=1,nauxil_3center
   iauxil_global = ibf_auxil_g(iauxil_local)
   rho_coeff_local_nospin(iauxil_local) = SUM(rho_coeff(iauxil_global,:))
 enddo

 ! vH_\alpha\beta = \sum_I ( \alpha \beta | 1/r12 | I ) * R_I
 call DGEMV('N',npair,nauxil_3center,1.0d0,eri_3center,npair,rho_coeff_local_nospin,1,0.0d0,vh,1)

 !$OMP PARALLEL PRIVATE(kbf,lbf)
 !$OMP DO
 do ipair=1,npair
   kbf = index_basis(1,ipair)
   lbf = index_basis(2,ipair)
   hartree_ao(kbf,lbf) = vh(ipair)
   hartree_ao(lbf,kbf) = vh(ipair)
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 ! Do not forget that the eri_3center(ibf,ibf | P ) included a factor 0.50
 do ibf=1,nbf
   hartree_ao(ibf,ibf) = hartree_ao(ibf,ibf) * 2.0_dp
 enddo
 deallocate(rho_coeff_local_nospin,vh)

 !
 ! Sum up the different contribution from different procs
 call xsum_world(hartree_ao)

 call dump_out_matrix(.FALSE.,'=== Hartree contribution ===',hartree_ao)

 select type(p_matrix)
 type is(real(dp))
   ehartree = 0.5_dp * SUM( hartree_ao(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
 type is(complex(dp))
   ehartree = 0.5_dp * SUM( hartree_ao(:,:) * SUM(REAL(p_matrix(:,:,:),dp),DIM=3) )
 end select


 call stop_clock(timing_xxdft_hartree)


end subroutine setup_hartree_genuine_ri


!=========================================================================
subroutine setup_exchange(p_matrix,exchange_ao,eexchange)
 implicit none
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: exchange_ao(:,:,:)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: nbf
 integer              :: ibf,jbf,kbf,lbf,ispin
 integer(kind=int8)   :: iint
 integer              :: index_ik,index_lj,stride
#if defined(_OPENMP)
 integer,external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
!=====

 call start_clock(timing_exchange)

 write(stdout,*) 'Calculate Exchange term using the 8 permutation symmetries'
 !$OMP PARALLEL PRIVATE(index_ik,index_lj,ibf,jbf,kbf,lbf,stride)
 exchange_ao(:,:,:) = 0.0_dp

#if defined(_OPENMP)
 stride = OMP_GET_NUM_THREADS()
 index_ik = OMP_GET_THREAD_NUM() - stride + 1
#else
 stride = 1
 index_ik = 0
#endif
 index_lj = 1

 ! SCHEDULE(static,1) should not be modified
 ! else a race competition will occur when performing index_ik = index_ik + stride
 !$OMP DO REDUCTION(+:exchange_ao) SCHEDULE(static,1)
 do iint=1,nint_4center
   index_ik = index_ik + stride
   do while( index_ik > npair )
     index_lj = index_lj + 1
     index_ik = index_lj + index_ik - npair - 1
   enddo

   ibf = index_basis(1,index_ik)
   kbf = index_basis(2,index_ik)
   lbf = index_basis(1,index_lj)
   jbf = index_basis(2,index_lj)


   exchange_ao(ibf,jbf,:) = exchange_ao(ibf,jbf,:) &
                  - eri_4center(iint) * p_matrix(kbf,lbf,:) / spin_fact
   if( ibf /= kbf ) then
     exchange_ao(kbf,jbf,:) = exchange_ao(kbf,jbf,:) &
                    - eri_4center(iint) * p_matrix(ibf,lbf,:) / spin_fact
     if( lbf /= jbf ) then
       exchange_ao(kbf,lbf,:) = exchange_ao(kbf,lbf,:) &
                      - eri_4center(iint) * p_matrix(ibf,jbf,:) / spin_fact
     endif
   endif
   if( lbf /= jbf ) then
     exchange_ao(ibf,lbf,:) = exchange_ao(ibf,lbf,:) &
                    - eri_4center(iint) * p_matrix(kbf,jbf,:) / spin_fact
   endif

   if( index_ik /= index_lj ) then
     exchange_ao(lbf,kbf,:) = exchange_ao(lbf,kbf,:) &
                    - eri_4center(iint) * p_matrix(jbf,ibf,:) / spin_fact
     if( ibf /= kbf ) then
       exchange_ao(lbf,ibf,:) = exchange_ao(lbf,ibf,:) &
                      - eri_4center(iint) * p_matrix(jbf,kbf,:) / spin_fact
       if( lbf /= jbf ) then
         exchange_ao(jbf,ibf,:) = exchange_ao(jbf,ibf,:) &
                        - eri_4center(iint) * p_matrix(lbf,kbf,:) / spin_fact
       endif
     endif
     if( lbf /= jbf ) then
       exchange_ao(jbf,kbf,:) = exchange_ao(jbf,kbf,:) &
                      - eri_4center(iint) * p_matrix(lbf,ibf,:) / spin_fact
     endif
   endif

 enddo
 !$OMP END DO

 !$OMP WORKSHARE
 eexchange = 0.5_dp * SUM( exchange_ao(:,:,:) * p_matrix(:,:,:) )
 !$OMP END WORKSHARE
 !$OMP END PARALLEL

 call stop_clock(timing_exchange)

end subroutine setup_exchange


!=========================================================================
subroutine setup_exchange_longrange(p_matrix,exchange_ao,eexchange)
 implicit none
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: exchange_ao(:,:,:)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: nbf
 integer              :: ibf,jbf,kbf,lbf,ispin
 integer(kind=int8)   :: iint
 integer              :: index_ik,index_lj,stride
#if defined(_OPENMP)
 integer,external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
!=====

 call start_clock(timing_exchange)

 write(stdout,*) 'Calculate LR-Exchange term using the 8 permutation symmetries'
 !$OMP PARALLEL PRIVATE(index_ik,index_lj,ibf,jbf,kbf,lbf,stride)
 exchange_ao(:,:,:) = 0.0_dp

#if defined(_OPENMP)
 stride = OMP_GET_NUM_THREADS()
 index_ik = OMP_GET_THREAD_NUM() - stride + 1
#else
 stride = 1
 index_ik = 0
#endif
 index_lj = 1

 ! SCHEDULE(static,1) should not be modified
 ! else a race competition will occur when performing index_ik = index_ik + stride
 !$OMP DO REDUCTION(+:exchange_ao) SCHEDULE(static,1)
 do iint=1,nint_4center
   index_ik = index_ik + stride
   do while( index_ik > npair )
     index_lj = index_lj + 1
     index_ik = index_lj + index_ik - npair - 1
   enddo

   ibf = index_basis(1,index_ik)
   kbf = index_basis(2,index_ik)
   lbf = index_basis(1,index_lj)
   jbf = index_basis(2,index_lj)


   exchange_ao(ibf,jbf,:) = exchange_ao(ibf,jbf,:) &
                  - eri_4center_lr(iint) * p_matrix(kbf,lbf,:) / spin_fact
   if( ibf /= kbf ) then
     exchange_ao(kbf,jbf,:) = exchange_ao(kbf,jbf,:) &
                    - eri_4center_lr(iint) * p_matrix(ibf,lbf,:) / spin_fact
     if( lbf /= jbf ) then
       exchange_ao(kbf,lbf,:) = exchange_ao(kbf,lbf,:) &
                      - eri_4center_lr(iint) * p_matrix(ibf,jbf,:) / spin_fact
     endif
   endif
   if( lbf /= jbf ) then
     exchange_ao(ibf,lbf,:) = exchange_ao(ibf,lbf,:) &
                    - eri_4center_lr(iint) * p_matrix(kbf,jbf,:) / spin_fact
   endif

   if( index_ik /= index_lj ) then
     exchange_ao(lbf,kbf,:) = exchange_ao(lbf,kbf,:) &
                    - eri_4center_lr(iint) * p_matrix(jbf,ibf,:) / spin_fact
     if( ibf /= kbf ) then
       exchange_ao(lbf,ibf,:) = exchange_ao(lbf,ibf,:) &
                      - eri_4center_lr(iint) * p_matrix(jbf,kbf,:) / spin_fact
       if( lbf /= jbf ) then
         exchange_ao(jbf,ibf,:) = exchange_ao(jbf,ibf,:) &
                        - eri_4center_lr(iint) * p_matrix(lbf,kbf,:) / spin_fact
       endif
     endif
     if( lbf /= jbf ) then
       exchange_ao(jbf,kbf,:) = exchange_ao(jbf,kbf,:) &
                      - eri_4center_lr(iint) * p_matrix(lbf,ibf,:) / spin_fact
     endif
   endif

 enddo
 !$OMP END DO

 !$OMP WORKSHARE
 eexchange = 0.5_dp * SUM( exchange_ao(:,:,:) * p_matrix(:,:,:) )
 !$OMP END WORKSHARE
 !$OMP END PARALLEL

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange


!=========================================================================
subroutine setup_exchange_ri(occupation,c_matrix,p_matrix,exchange_ao,eexchange)
 implicit none
 real(dp),intent(in)  :: occupation(:,:)
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: exchange_ao(:,:,:)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: nauxil_local
 integer              :: nbf,nstate
 integer              :: ibf,jbf,ispin,istate
 integer              :: nocc
 real(dp),allocatable :: tmp(:,:),c_t(:,:)
 integer              :: ipair,iauxil
 integer              :: ibf_auxil_first,nbf_auxil_core
!=====

 call start_clock(timing_exchange)

 write(stdout,*) 'Calculate Exchange term with Resolution-of-Identity'

 exchange_ao(:,:,:) = 0.0_dp

 ! Find highest occupied state
 nocc = get_number_occupied_states(occupation)

 nbf    = SIZE(exchange_ao,DIM=1)
 nstate = SIZE(occupation(:,:),DIM=1)

 nauxil_local = nauxil_3center

 allocate(tmp(nocc,nbf))
 allocate(c_t(nocc,nbf))

 do ispin=1,nspin

   !$OMP PARALLEL DO
   do ibf=1,nbf
     c_t(:,ibf) = c_matrix(ibf,1:nocc,ispin) * SQRT( occupation(1:nocc,ispin) / spin_fact )
   enddo
   !$OMP END PARALLEL DO

   do iauxil=1,nauxil_local
     if( MODULO( iauxil - 1 , nproc_ortho ) /= rank_ortho ) cycle
     tmp(:,:) = 0.0_dp
     !$OMP PARALLEL PRIVATE(ibf,jbf)
     !$OMP DO REDUCTION(+:tmp)
     do ipair=1,npair
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_t(:,jbf) * eri_3center(ipair,iauxil)
       tmp(:,jbf) = tmp(:,jbf) + c_t(:,ibf) * eri_3center(ipair,iauxil)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     ! exchange_ao(:,:,ispin) = exchange_ao(:,:,ispin) &
     !                    - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
     ! C = A^T * A + C
     call DSYRK('L','T',nbf,nocc,-1.0_dp,tmp(1,1),nocc,1.0_dp,exchange_ao(1,1,ispin),nbf)

   enddo
 enddo

 deallocate(c_t)
 deallocate(tmp)


 !
 ! Need to symmetrize exchange_ao
 do ispin=1,nspin
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ao(ibf,jbf,ispin) = exchange_ao(jbf,ibf,ispin)
     enddo
   enddo
 enddo
 call xsum_world(exchange_ao)

 eexchange = 0.5_dp * SUM( exchange_ao(:,:,:) * p_matrix(:,:,:) )

 call stop_clock(timing_exchange)

end subroutine setup_exchange_ri


!=========================================================================
subroutine setup_exchange_longrange_ri(occupation,c_matrix,p_matrix,exchange_ao,eexchange)
 implicit none
 real(dp),intent(in)  :: occupation(:,:)
 real(dp),intent(in)  :: c_matrix(:,:,:)
 real(dp),intent(in)  :: p_matrix(:,:,:)
 real(dp),intent(out) :: exchange_ao(:,:,:)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: nauxil_local
 integer              :: nbf,nstate
 integer              :: ibf,jbf,ispin,istate
 integer              :: nocc
 real(dp),allocatable :: tmp(:,:),c_t(:,:)
 integer              :: ipair,iauxil
 integer              :: ibf_auxil_first,nbf_auxil_core
!=====

 call start_clock(timing_exchange)

 write(stdout,*) 'Calculate LR Exchange term with Resolution-of-Identity'

 exchange_ao(:,:,:) = 0.0_dp

 ! Find highest occupied state
 nocc = get_number_occupied_states(occupation)

 nbf    = SIZE(exchange_ao,DIM=1)
 nstate = SIZE(occupation(:,:),DIM=1)

 nauxil_local = nauxil_3center_lr


 allocate(tmp(nocc,nbf))
 allocate(c_t(nocc,nbf))

 do ispin=1,nspin

   !$OMP PARALLEL DO
   do ibf=1,nbf
     c_t(:,ibf) = c_matrix(ibf,1:nocc,ispin) * SQRT( occupation(1:nocc,ispin) / spin_fact )
   enddo
   !$OMP END PARALLEL DO


   do iauxil=1,nauxil_local
     if( MODULO( iauxil - 1 , nproc_ortho ) /= rank_ortho ) cycle
     tmp(:,:) = 0.0_dp
     !$OMP PARALLEL PRIVATE(ibf,jbf)
     !$OMP DO REDUCTION(+:tmp)
     do ipair=1,npair
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_t(:,jbf) * eri_3center_lr(ipair,iauxil)
       tmp(:,jbf) = tmp(:,jbf) + c_t(:,ibf) * eri_3center_lr(ipair,iauxil)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     ! exchange_ao(:,:,ispin) = exchange_ao(:,:,ispin) &
     !                    - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
     ! C = A^T * A + C
     call DSYRK('L','T',nbf,nocc,-1.0_dp,tmp(1,1),nocc,1.0_dp,exchange_ao(1,1,ispin),nbf)

   enddo
 enddo

 deallocate(c_t)
 deallocate(tmp)


 !
 ! Need to symmetrize exchange_ao
 do ispin=1,nspin
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ao(ibf,jbf,ispin) = exchange_ao(jbf,ibf,ispin)
     enddo
   enddo
 enddo
 call xsum_world(exchange_ao)

 eexchange = 0.5_dp * SUM( exchange_ao(:,:,:) * p_matrix(:,:,:) )

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange_ri


!=========================================================================
subroutine setup_exchange_ri_cmplx(occupation,c_matrix,p_matrix,exchange_ao,eexchange)
 implicit none
 real(dp),intent(in)     :: occupation(:,:)
 complex(dp),intent(in)  :: c_matrix(:,:,:)
 complex(dp),intent(in)  :: p_matrix(:,:,:)
 complex(dp),intent(out) :: exchange_ao(:,:,:)
 real(dp),intent(out)    :: eexchange
!=====
 integer                 :: nbf,nstate
 integer                 :: nocc
 integer                 :: ibf,jbf,ispin,istate
 complex(dp),allocatable :: tmp_cmplx(:,:),c_t_cmplx(:,:)
 integer                 :: ipair,iauxil
 integer                 :: ibf_auxil_first,nbf_auxil_core
!=====

 call start_clock(timing_exchange)

 write(stdout,*) 'Calculate Complex Exchange term with Resolution-of-Identity'

 exchange_ao(:,:,:) = (0.0_dp, 0.0_dp)

 ! Find highest occupied state
 nocc = get_number_occupied_states(occupation)

 nbf    = SIZE(exchange_ao,DIM=1)
 nstate = SIZE(occupation(:,:),DIM=1)

 allocate(tmp_cmplx(nocc,nbf))
 allocate(c_t_cmplx(nocc,nbf))

 do ispin=1,nspin

   !$OMP PARALLEL DO
   do ibf=1,nbf
     c_t_cmplx(:,ibf) = CONJG( c_matrix(ibf,1:nocc,ispin) ) * SQRT( occupation(1:nocc,ispin) / spin_fact )
   enddo
   !$OMP END PARALLEL DO

   do iauxil=1,nauxil_3center
     if( MODULO( iauxil - 1 , nproc_ortho ) /= rank_ortho ) cycle
     tmp_cmplx(:,:) = (0.0_dp, 0.0_dp)
     !$OMP PARALLEL PRIVATE(ibf,jbf)
     !$OMP DO REDUCTION(+:tmp_cmplx)
     do ipair=1,npair
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)
       tmp_cmplx(:,ibf) = tmp_cmplx(:,ibf) + c_t_cmplx(:,jbf) * eri_3center(ipair,iauxil)
       tmp_cmplx(:,jbf) = tmp_cmplx(:,jbf) + c_t_cmplx(:,ibf) * eri_3center(ipair,iauxil)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
     ! exchange_ao(:,:,ispin) = exchange_ao(:,:,ispin) &
     !                    - MATMUL( CONJG(TRANSPOSE(tmp(:,:))) , tmp(:,:) ) / spin_fact
     ! C = A^H * A + C
     call ZHERK('L','C',nbf,nocc,-1.0_dp,tmp_cmplx,nocc,1.0_dp,exchange_ao(1,1,ispin),nbf)

   enddo
 enddo

 deallocate(c_t_cmplx)
 deallocate(tmp_cmplx)


 !
 ! Need to hermitianize exchange_ao
 do ispin=1,nspin
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ao(ibf,jbf,ispin) = CONJG( exchange_ao(jbf,ibf,ispin) )
     enddo
   enddo
 enddo
 call xsum_world(exchange_ao)

 eexchange = 0.5_dp * REAL( SUM( exchange_ao(:,:,:) * CONJG( p_matrix(:,:,:) ) ) , dp)

 call stop_clock(timing_exchange)

end subroutine setup_exchange_ri_cmplx


!=========================================================================
! Calculate the exchange-correlation potential and energy
! * subroutine works for both real and complex wavefunctions c_matrix
!   using "class" syntax of Fortran2003
subroutine dft_exc_vxc_batch(batch_size,basis,occupation,c_matrix,vxc_ao,exc_xc)
 implicit none

 integer,intent(in)         :: batch_size
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(:,:)
 class(*),intent(in)        :: c_matrix(:,:,:)
 real(dp),intent(out)       :: vxc_ao(:,:,:)
 real(dp),intent(out)       :: exc_xc
!=====
 real(dp),parameter   :: TOL_RHO=1.0e-9_dp
 integer              :: nstate
 integer              :: ibf,jbf,ispin
 integer              :: ixc
 integer              :: igrid_start,igrid_end,ir
 integer              :: timing_xxdft_xc,timing_xxdft_densities,timing_xxdft_libxc,timing_xxdft_vxc
 real(dp)             :: normalization(nspin)
 real(dp),allocatable :: weight_batch(:)
 real(dp),allocatable :: tmp_batch(:,:)
 real(dp),allocatable :: basis_function_r_batch(:,:)
 real(dp),allocatable :: bf_gradx_batch(:,:)
 real(dp),allocatable :: bf_grady_batch(:,:)
 real(dp),allocatable :: bf_gradz_batch(:,:)
 real(dp),allocatable :: dedd_r_batch(:,:)
 real(dp),allocatable :: grad_rhor_batch(:,:,:)
 real(dp),allocatable :: dedgd_r_batch(:,:,:)
 integer(C_INT)             :: nr
 real(C_DOUBLE),allocatable :: rhor_batch(:,:)
 real(C_DOUBLE),allocatable :: sigma_batch(:,:)
 real(C_DOUBLE),allocatable :: vrho_batch(:,:)
 real(C_DOUBLE),allocatable :: exc_batch(:)
 real(C_DOUBLE),allocatable :: vsigma_batch(:,:)
!=====

 exc_xc = 0.0_dp
 vxc_ao(:,:,:) = 0.0_dp

 if( dft_xc(1)%nxc == 0 ) return

 !
 ! Switch the timers: complex wavefunctions imply RT-TDDFT run
 select type(c_matrix)
 type is(real(dp))
   timing_xxdft_xc        = timing_dft_xc
   timing_xxdft_densities = timing_dft_densities
   timing_xxdft_libxc     = timing_dft_libxc
   timing_xxdft_vxc       = timing_dft_vxc
 type is(complex(dp))
   timing_xxdft_xc        = timing_tddft_xc
   timing_xxdft_densities = timing_tddft_densities
   timing_xxdft_libxc     = timing_tddft_libxc
   timing_xxdft_vxc       = timing_tddft_vxc
 class default
   call die("dft_exc_vxc_batch: c_matrix is neither real nor complex")
 end select

 call start_clock(timing_xxdft_xc)

 nstate = SIZE(occupation,DIM=1)

#if defined(HAVE_LIBXC)

 write(stdout,*) 'Calculate DFT XC potential'
 if( batch_size /= 1 ) write(stdout,*) 'Using batches of size',batch_size


 normalization(:) = 0.0_dp

 !
 ! Loop over batches of grid points
 !
 do igrid_start=1,ngrid,batch_size
   igrid_end = MIN(ngrid,igrid_start+batch_size-1)
   nr = igrid_end - igrid_start + 1

   call start_clock(timing_xxdft_densities)

   allocate(weight_batch(nr))
   allocate(rhor_batch(nspin,nr))
   allocate(basis_function_r_batch(basis%nbf,nr))
   allocate(exc_batch(nr))
   allocate(vrho_batch(nspin,nr))
   allocate(dedd_r_batch(nspin,nr))

   if( dft_xc(1)%needs_gradient ) then
     allocate(bf_gradx_batch(basis%nbf,nr))
     allocate(bf_grady_batch(basis%nbf,nr))
     allocate(bf_gradz_batch(basis%nbf,nr))
     allocate(grad_rhor_batch(nspin,nr,3))
     allocate(dedgd_r_batch(3,nr,nspin))
     allocate(sigma_batch(2*nspin-1,nr))
     allocate(vsigma_batch(2*nspin-1,nr))
   endif

   weight_batch(:) = w_grid(igrid_start:igrid_end)

   call get_basis_functions_r_batch(basis,igrid_start,basis_function_r_batch)
   !
   ! Get the gradient at points r
   if( dft_xc(1)%needs_gradient ) &
      call get_basis_functions_gradr_batch(basis,igrid_start,bf_gradx_batch,bf_grady_batch,bf_gradz_batch)

   !
   ! Calculate the density at points r for spin up and spin down
   ! Calculate grad rho at points r for spin up and spin down
   if( .NOT. dft_xc(1)%needs_gradient ) then
     call calc_density_r_batch(occupation,c_matrix,basis_function_r_batch,rhor_batch)
   else
     call calc_density_gradr_batch(occupation,c_matrix,basis_function_r_batch, &
                                   bf_gradx_batch,bf_grady_batch,bf_gradz_batch,rhor_batch,grad_rhor_batch)

     !$OMP PARALLEL DO
     do ir=1,nr
       sigma_batch(1,ir) = DOT_PRODUCT( grad_rhor_batch(1,ir,:) , grad_rhor_batch(1,ir,:) )
       if( nspin == 2 ) then
         sigma_batch(2,ir) = DOT_PRODUCT( grad_rhor_batch(1,ir,:) , grad_rhor_batch(2,ir,:) )
         sigma_batch(3,ir) = DOT_PRODUCT( grad_rhor_batch(2,ir,:) , grad_rhor_batch(2,ir,:) )
       endif
     enddo
     !$OMP END PARALLEL DO

   endif

   ! Normalization
   normalization(:) = normalization(:) + MATMUL( rhor_batch(:,:) , weight_batch(:) )

   call stop_clock(timing_xxdft_densities)


   !
   ! LIBXC calls
   !
   call start_clock(timing_xxdft_libxc)

   dedd_r_batch(:,:) = 0.0_dp
   if( dft_xc(1)%needs_gradient ) dedgd_r_batch(:,:,:) = 0.0_dp

   do ixc=1,dft_xc(1)%nxc
     if( ABS(dft_xc(ixc)%coeff) < 1.0e-6_dp ) cycle

     select case(dft_xc(ixc)%family)
     case(XC_FAMILY_LDA)
       call xc_lda_exc_vxc(dft_xc(ixc)%func,nr,rhor_batch(1,1),exc_batch(1),vrho_batch(1,1))

     case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
       call xc_gga_exc_vxc(dft_xc(ixc)%func,nr,rhor_batch(1,1),sigma_batch(1,1),exc_batch(1),vrho_batch(1,1),vsigma_batch(1,1))

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
       write(stdout,*) 'Family id:',ixc,dft_xc(ixc)%family
       call die('dft_exc_vxc_batch: functional is not LDA nor GGA nor hybrid')
     end select

     ! XC energy
     exc_xc = exc_xc + SUM( weight_batch(:) * exc_batch(:) * SUM(rhor_batch(:,:),DIM=1) ) * dft_xc(ixc)%coeff

     dedd_r_batch(:,:) = dedd_r_batch(:,:) + vrho_batch(:,:) * dft_xc(ixc)%coeff

     !
     ! Set up divergence term if needed (GGA case)
     !
     if( dft_xc(1)%needs_gradient ) then
       if(nspin==1) then

         do ir=1,nr
           dedgd_r_batch(:,ir,1) = dedgd_r_batch(:,ir,1)  &
                      + 2.0_dp * vsigma_batch(1,ir) * grad_rhor_batch(1,ir,:) * dft_xc(ixc)%coeff
         enddo

       else

         do ir=1,nr
           dedgd_r_batch(:,ir,1) = dedgd_r_batch(:,ir,1) &
                     + ( 2.0_dp * vsigma_batch(1,ir) * grad_rhor_batch(1,ir,:) &
                                 + vsigma_batch(2,ir) * grad_rhor_batch(2,ir,:) ) * dft_xc(ixc)%coeff

           dedgd_r_batch(:,ir,2) = dedgd_r_batch(:,ir,2) &
                     + ( 2.0_dp * vsigma_batch(3,ir) * grad_rhor_batch(2,ir,:) &
                                 + vsigma_batch(2,ir) * grad_rhor_batch(1,ir,:) ) * dft_xc(ixc)%coeff
         enddo

       endif

     endif

   enddo ! loop on the XC functional

   call stop_clock(timing_xxdft_libxc)



   !
   ! Eventually set up the vxc term
   !
   call start_clock(timing_xxdft_vxc)

   allocate(tmp_batch(basis%nbf,nr))
   do ispin=1,nspin
     if( dft_xc(1)%needs_gradient ) then
       !
       ! GGA
       !
       !$OMP PARALLEL DO
       do ir=1,nr
         tmp_batch(:,ir) = weight_batch(ir) * dedd_r_batch(ispin,ir) * basis_function_r_batch(:,ir) * 0.50_dp
         tmp_batch(:,ir) = tmp_batch(:,ir) &
                          +  bf_gradx_batch(:,ir) * dedgd_r_batch(1,ir,ispin) * weight_batch(ir) &
                          +  bf_grady_batch(:,ir) * dedgd_r_batch(2,ir,ispin) * weight_batch(ir) &
                          +  bf_gradz_batch(:,ir) * dedgd_r_batch(3,ir,ispin) * weight_batch(ir)
       enddo
       !$OMP END PARALLEL DO
     else
       !
       ! LDA
       !
       !$OMP PARALLEL DO
       do ir=1,nr
         tmp_batch(:,ir) = weight_batch(ir) * dedd_r_batch(ispin,ir) * basis_function_r_batch(:,ir) * 0.50_dp
       enddo
       !$OMP END PARALLEL DO
     endif
     call DSYR2K('L','N',basis%nbf,nr,1.0d0,basis_function_r_batch,basis%nbf,tmp_batch,basis%nbf,1.0d0,vxc_ao(:,:,ispin),basis%nbf)

   enddo
   deallocate(tmp_batch)

   call stop_clock(timing_xxdft_vxc)


   deallocate(weight_batch)
   deallocate(basis_function_r_batch)
   deallocate(exc_batch)
   deallocate(rhor_batch)
   deallocate(vrho_batch)
   deallocate(dedd_r_batch)
   if( dft_xc(1)%needs_gradient ) then
     deallocate(bf_gradx_batch)
     deallocate(bf_grady_batch)
     deallocate(bf_gradz_batch)
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
       vxc_ao(jbf,ibf,ispin) = vxc_ao(ibf,jbf,ispin)
     enddo
   enddo
 enddo

 !
 ! Sum up the contributions from all procs only if needed
 call xsum_grid(normalization)
 call xsum_grid(vxc_ao)
 call xsum_grid(exc_xc)


#else
 write(stdout,*) 'XC energy and potential set to zero'
 write(stdout,*) 'LIBXC is not present'
#endif

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization(:)
 write(stdout,'(a,2x,f12.6,/)')    '  DFT xc energy (Ha):',exc_xc

 call stop_clock(timing_xxdft_xc)

end subroutine dft_exc_vxc_batch


!=========================================================================
subroutine dft_approximate_vhxc(basis,vhxc_ao)
 use m_dft_grid
 use m_eri_calculate
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: vhxc_ao(basis%nbf,basis%nbf)
!=====
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

 vhxc_ao(:,:) = 0.0_dp

 write(stdout,'(/,a)') ' Calculate approximate HXC potential with a superposition of atomic densities'

 do iatom=1,natom
   if( rank_grid /= MODULO(iatom,nproc_grid) ) cycle

   ngau = 4
   allocate(alpha(ngau),coeff(ngau))
   call element_atomicdensity(zvalence(iatom),zatom(iatom),coeff,alpha)


   call calculate_eri_approximate_hartree(basis,xatom(:,iatom),coeff,alpha,vhgau)
   vhxc_ao(:,:) = vhxc_ao(:,:) + vhgau(:,:)

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
   call get_basis_functions_r_batch(basis,igrid_start,basis_function_r_batch)

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
   call DSYRK('L','N',basis%nbf,nr,-1.0d0,basis_function_r_batch(1,1),basis%nbf,1.0d0,vhxc_ao(1,1),basis%nbf)


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
 call xsum_grid(vhxc_ao)

 ! Symmetrize vhxc
 do ibf=1,basis%nbf
   do jbf=ibf+1,basis%nbf
     vhxc_ao(ibf,jbf) = vhxc_ao(jbf,ibf)
   enddo
 enddo

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization
 write(stdout,  '(a,2(2x,f12.6))') '      XC energy (Ha):',exc

 !
 ! Temporary grid destroyed
 call destroy_dft_grid()

 call stop_clock(timing_approx_ham)

end subroutine dft_approximate_vhxc


end module m_hamiltonian_twobodies
!=========================================================================
