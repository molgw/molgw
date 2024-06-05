!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to perform the Atomic Orbital to Molecular Orbital transform
!
!=========================================================================
#include "molgw.h"
module m_eri_ao_mo
  use m_definitions
  use m_mpi
  use m_scalapack
  use m_memory
  use m_warning
  use m_basis_set
  use m_timing
  use m_eri
  use m_inputparam, only: nspin, has_auxil_basis
  use m_io, only: yaml_search_keyword


  real(dp),protected,allocatable :: eri_3center_eigen(:,:,:,:)

  real(dp),protected,allocatable :: eri_3center_eigen_lr(:,:,:,:)

  real(dp),protected,allocatable :: eri_4center_eigen_uks(:,:,:,:)

  complex(dp),protected,allocatable :: eri_3center_eigen_cmplx(:,:,:,:)


contains


!=========================================================================
function eri_eigen(istate,jstate,ijspin,kstate,lstate,klspin)
  implicit none
  integer,intent(in) :: ijspin,klspin
  integer,intent(in) :: istate,jstate,kstate,lstate
  real(dp)           :: eri_eigen
  !=====

  if(has_auxil_basis) then
    eri_eigen = DOT_PRODUCT( eri_3center_eigen(:,istate,jstate,ijspin) , eri_3center_eigen(:,kstate,lstate,klspin) )
    call auxil%sum(eri_eigen)
  else
    eri_eigen = eri_4center_eigen_uks(istate,jstate,kstate,lstate)
  endif


end function eri_eigen


!=========================================================================
function eri_eigen_ri(istate,jstate,ijspin,kstate,lstate,klspin)
  implicit none
  integer,intent(in) :: ijspin,klspin
  integer,intent(in) :: istate,jstate,kstate,lstate
  real(dp)           :: eri_eigen_ri
  !=====

  eri_eigen_ri = DOT_PRODUCT( eri_3center_eigen(:,istate,jstate,ijspin) , eri_3center_eigen(:,kstate,lstate,klspin) )

  call auxil%sum(eri_eigen_ri)

end function eri_eigen_ri


!=========================================================================
function eri_eigen_ri_lr(istate,jstate,ijspin,kstate,lstate,klspin)
  implicit none
  integer,intent(in) :: ijspin,klspin
  integer,intent(in) :: istate,jstate,kstate,lstate
  real(dp)           :: eri_eigen_ri_lr
  !=====

  eri_eigen_ri_lr = DOT_PRODUCT( eri_3center_eigen_lr(:,istate,jstate,ijspin) , eri_3center_eigen_lr(:,kstate,lstate,klspin) )

  call auxil%sum(eri_eigen_ri_lr)

end function eri_eigen_ri_lr


!=========================================================================
function eri_eigen_ri_cmplx(istate,jstate,ijspin,kstate,lstate,klspin)
  implicit none
  integer,intent(in) :: ijspin,klspin
  integer,intent(in) :: istate,jstate,kstate,lstate
  complex(dp)        :: eri_eigen_ri_cmplx
  !=====

  eri_eigen_ri_cmplx = SUM( eri_3center_eigen_cmplx(:,istate,jstate,ijspin) * eri_3center_eigen_cmplx(:,kstate,lstate,klspin) )

  call auxil%sum(eri_eigen_ri_cmplx)

end function eri_eigen_ri_cmplx


!=========================================================================
pure function eri_eigen_ri_paral(istate,jstate,ijspin,kstate,lstate,klspin)
  implicit none
  integer,intent(in) :: ijspin,klspin
  integer,intent(in) :: istate,jstate,kstate,lstate
  real(dp)           :: eri_eigen_ri_paral
  !=====

  eri_eigen_ri_paral = DOT_PRODUCT( eri_3center_eigen(:,istate,jstate,ijspin) , eri_3center_eigen(:,kstate,lstate,klspin) )

end function eri_eigen_ri_paral


!=================================================================
subroutine calculate_eri_4center_eigen(c_matrix,istate,ijspin,eri_eigenstate_i)
  implicit none

  integer,intent(in)     :: istate,ijspin
  real(dp),intent(in)    :: c_matrix(:,:,:)
  real(dp),intent(inout) :: eri_eigenstate_i(:,:,:,:)
  !=====
  logical,save         :: disclaimer = .TRUE.
  integer,save         :: istate_previous=0
  integer,save         :: ijspin_previous=0
  integer              :: nbf,nstate
  integer              :: klspin
  integer              :: ibf,jbf,kbf,lbf
  integer              :: jstate,kstate,lstate
  real(dp),allocatable :: eri_tmp3(:,:,:),eri_tmp2(:,:,:),eri_tmp1(:,:)
  integer(kind=int8)   :: iint
  integer              :: index_ij,index_kl,stride
  !=====

  nbf    = SIZE(c_matrix,DIM=1)
  nstate = SIZE(c_matrix,DIM=2)

  ! Check if the calculation can be skipped
  if( istate_previous == istate .AND. ijspin_previous == ijspin .AND. ANY(ABS(eri_eigenstate_i(:,:,:,:))>1.0e-6_dp) ) then
    return
  else
    istate_previous = istate
    ijspin_previous = ijspin
  endif


  call start_clock(timing_eri_4center_eigen)

  allocate(eri_tmp1(nbf,nstate))
  call clean_allocate('TMP array',eri_tmp2,nstate,nbf,nbf)
  call clean_allocate('TMP array',eri_tmp3,nbf,nbf,nbf)

  eri_tmp3(:,:,:) = 0.0_dp
  eri_eigenstate_i(:,:,:,:) = 0.0_dp


  ! New implementation looping over the unique integrals
  ! write this disclaimer only at the first call
  if( disclaimer ) then
    write(stdout,'(1x,a)')        'AO to MO transform using the 8 permutation symmetries'
    write(stdout,'(1x,a,f9.3,a)') 'Make sure OMP_STACKSIZE is larger than ',REAL(nbf,dp)**3 * REAL(dp,dp) / 1024.0_dp**2,' (Mb)'
    disclaimer = .FALSE.
  endif
  !$OMP PARALLEL PRIVATE(index_ij,index_kl,ibf,jbf,kbf,lbf,stride)
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
  !$OMP DO REDUCTION(+:eri_tmp3) SCHEDULE(static,1)
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

    eri_tmp3(jbf,kbf,lbf) = eri_tmp3(jbf,kbf,lbf) + eri_4center(iint) * c_matrix(ibf,istate,ijspin)
    if( ibf /= jbf ) then
      eri_tmp3(ibf,kbf,lbf) = eri_tmp3(ibf,kbf,lbf) + eri_4center(iint) * c_matrix(jbf,istate,ijspin)
      if( kbf /= lbf ) then
        eri_tmp3(ibf,lbf,kbf) = eri_tmp3(ibf,lbf,kbf) + eri_4center(iint) * c_matrix(jbf,istate,ijspin)
      endif
    endif
    if( kbf /= lbf ) then
      eri_tmp3(jbf,lbf,kbf) = eri_tmp3(jbf,lbf,kbf) + eri_4center(iint) * c_matrix(ibf,istate,ijspin)
    endif
    if( index_kl /= index_ij ) then
      eri_tmp3(lbf,ibf,jbf) = eri_tmp3(lbf,ibf,jbf) + eri_4center(iint) * c_matrix(kbf,istate,ijspin)
      if( ibf /= jbf ) then
        eri_tmp3(lbf,jbf,ibf) = eri_tmp3(lbf,jbf,ibf) + eri_4center(iint) * c_matrix(kbf,istate,ijspin)
        if( kbf /= lbf ) then
          eri_tmp3(kbf,jbf,ibf) = eri_tmp3(kbf,jbf,ibf) + eri_4center(iint) * c_matrix(lbf,istate,ijspin)
        endif
      endif
      if( kbf /= lbf ) then
        eri_tmp3(kbf,ibf,jbf) = eri_tmp3(kbf,ibf,jbf) + eri_4center(iint) * c_matrix(lbf,istate,ijspin)
      endif
    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


  call DGEMM('T','N',nstate,nbf*nbf,nbf,1.0d0,c_matrix(1,1,ijspin),nbf, &
                                              eri_tmp3(1,1,1),nbf,    &
                                        0.0d0,eri_tmp2(1,1,1),nstate)

  do klspin=1,nspin
    do lbf=1,nbf
      call DGEMM('N','N',nstate,nstate,nbf,1.0d0,eri_tmp2(1,1,lbf),nstate,    &
                                                 c_matrix(1,1,klspin),nbf, &
                                           0.0d0,eri_tmp3(1,1,lbf),nbf)
    enddo

    do lstate=1,nstate
      eri_tmp1(:,:) = TRANSPOSE( eri_tmp3(1:nstate,lstate,1:nbf) )

      call DGEMM('T','N',nstate,nstate,nbf,1.0d0,eri_tmp1,nbf,              &
                                                 c_matrix(1,1,klspin),nbf,  &
                                           0.0d0,eri_eigenstate_i(1,1,lstate,klspin),nstate)

    enddo
  enddo !klspin

  deallocate(eri_tmp1)
  call clean_deallocate('TMP array',eri_tmp2)
  call clean_deallocate('TMP array',eri_tmp3)

  call stop_clock(timing_eri_4center_eigen)

end subroutine calculate_eri_4center_eigen


!=================================================================
subroutine calculate_eri_4center_eigen_uks(c_matrix,nstate_min,nstate_max)
  implicit none

  real(dp),intent(in)    :: c_matrix(:,:,:)
  integer,intent(in)     :: nstate_min,nstate_max
  !=====
  integer              :: nbf,nstate,nstate_maxmin
  integer              :: ijspin,klspin
  integer              :: ibf,jbf,kbf,lbf
  integer              :: istate,jstate,kstate,lstate
  real(dp),allocatable :: eri_tmp3(:,:,:),eri_tmp2(:,:,:),eri_tmp1(:,:),eri_tmp1b(:,:)
  integer(kind=int8)   :: iint
  integer              :: index_ij,index_kl,stride
  !=====

  nbf    = SIZE(c_matrix,DIM=1)
  nstate = SIZE(c_matrix,DIM=2)

  write(stdout,'(/,1x,a)') 'Calculate all the 4-center MO integrals at once (using the 8 permutation symmetries)'
  write(stdout,'(1x,a,f9.3,a)') 'Make sure OMP_STACKSIZE is larger than ',REAL(nbf,dp)**3 * 8.0_dp / 1024.0_dp**2,' (Mb)'

  if( nspin /= 1 ) call die('calculate_eri_4center_eigen_uks: requires spin-restricted calculation')
  ! spin-unrestricted will be coded later
  ijspin = 1

  call start_clock(timing_eri_4center_eigen)

  call clean_allocate('4-center MO integrals',eri_4center_eigen_uks, &
                      nstate_min,nstate_max,nstate_min,nstate_max,nstate_min,nstate_max,nstate_min,nstate_max)
  eri_4center_eigen_uks(:,:,:,:) = 0.0_dp
  nstate_maxmin = nstate_max - nstate_min + 1

  allocate(eri_tmp1(nstate_maxmin,nbf))
  allocate(eri_tmp1b(nbf,nbf))
  call clean_allocate('TMP array',eri_tmp2,nstate_maxmin,nbf,nbf)
  call clean_allocate('TMP array',eri_tmp3,nbf,nbf,nbf)


  do istate=nstate_min,nstate_max
    if( MODULO( istate - nstate_min , ortho%nproc ) /= ortho%rank ) cycle

    eri_tmp3(:,:,:) = 0.0_dp

    !$OMP PARALLEL PRIVATE(index_ij,index_kl,ibf,jbf,kbf,lbf,stride)
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
    !$OMP DO REDUCTION(+:eri_tmp3) SCHEDULE(static,1)
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

      eri_tmp3(jbf,kbf,lbf) = eri_tmp3(jbf,kbf,lbf) + eri_4center(iint) * c_matrix(ibf,istate,ijspin)
      if( ibf /= jbf ) then
        eri_tmp3(ibf,kbf,lbf) = eri_tmp3(ibf,kbf,lbf) + eri_4center(iint) * c_matrix(jbf,istate,ijspin)
        if( kbf /= lbf ) then
          eri_tmp3(ibf,lbf,kbf) = eri_tmp3(ibf,lbf,kbf) + eri_4center(iint) * c_matrix(jbf,istate,ijspin)
        endif
      endif
      if( kbf /= lbf ) then
        eri_tmp3(jbf,lbf,kbf) = eri_tmp3(jbf,lbf,kbf) + eri_4center(iint) * c_matrix(ibf,istate,ijspin)
      endif
      if( index_kl /= index_ij ) then
        eri_tmp3(lbf,ibf,jbf) = eri_tmp3(lbf,ibf,jbf) + eri_4center(iint) * c_matrix(kbf,istate,ijspin)
        if( ibf /= jbf ) then
          eri_tmp3(lbf,jbf,ibf) = eri_tmp3(lbf,jbf,ibf) + eri_4center(iint) * c_matrix(kbf,istate,ijspin)
          if( kbf /= lbf ) then
            eri_tmp3(kbf,jbf,ibf) = eri_tmp3(kbf,jbf,ibf) + eri_4center(iint) * c_matrix(lbf,istate,ijspin)
          endif
        endif
        if( kbf /= lbf ) then
          eri_tmp3(kbf,ibf,jbf) = eri_tmp3(kbf,ibf,jbf) + eri_4center(iint) * c_matrix(lbf,istate,ijspin)
        endif
      endif
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call DGEMM('T','N',nstate_maxmin,nbf*nbf,nbf,1.0d0,c_matrix(1,nstate_min,ijspin),nbf, &
                                                   eri_tmp3(1,1,1),nbf,             &
                                             0.0d0,eri_tmp2(1,1,1),nstate_maxmin)

    do klspin=1,nspin
      do jstate=nstate_min,nstate_max
        eri_tmp1b(:,:) = eri_tmp2(jstate-nstate_min+1,:,1:nbf)
        call DGEMM('T','N',nstate_maxmin,nbf,nbf,1.0d0,c_matrix(1,nstate_min,klspin),nbf,   &
                                                                 eri_tmp1b(1,1),nbf,                  &
                                                           0.0d0,eri_tmp1(nstate_min,1),nstate_maxmin)

        call DGEMM('N','N',nstate_maxmin,nstate_maxmin,nbf, &
                                              1.0d0,eri_tmp1(1,1),nstate_maxmin,              &
                                                    c_matrix(1,nstate_min,klspin),nbf,  &
                                              0.0d0,eri_4center_eigen_uks(nstate_min,nstate_min,jstate,istate),nstate_maxmin)
      enddo
    enddo !klspin

  enddo !istate

  deallocate(eri_tmp1)
  deallocate(eri_tmp1b)
  call clean_deallocate('TMP array',eri_tmp2)
  call clean_deallocate('TMP array',eri_tmp3)

  call ortho%sum(eri_4center_eigen_uks)

  call stop_clock(timing_eri_4center_eigen)

end subroutine calculate_eri_4center_eigen_uks


!=================================================================
subroutine destroy_eri_4center_eigen_uks()
  implicit none
  !=====
  !=====

  call clean_deallocate('4-center MO integrals',eri_4center_eigen_uks)

end subroutine destroy_eri_4center_eigen_uks


!=================================================================
subroutine calculate_eri_3center_eigen(c_matrix,mstate_min,mstate_max,nstate_min,nstate_max,timing,verbose,long_range)
  implicit none
  real(dp),intent(in)         :: c_matrix(:,:,:)
  integer,optional,intent(in) :: mstate_min,mstate_max,nstate_min,nstate_max
  integer,optional,intent(in) :: timing
  logical,optional,intent(in) :: verbose,long_range
  !=====
  logical              :: verbose_,lr_exists=.true.
  integer              :: nbf,nstate
  integer              :: mstate_min_,mstate_max_,nstate_min_,nstate_max_
  integer              :: mstate_count_,nstate_count_
  integer              :: kbf,lbf,iauxil
  integer              :: lstate
  integer              :: klspin
  real(dp),allocatable :: tmp1(:,:),tmp2(:,:),c_t(:,:)
  integer              :: ipair
  !=====

  if( PRESENT(long_range) ) then
    if( long_range .AND. (.NOT. ALLOCATED(eri_3center_eigen_lr)) ) lr_exists=.false.
  endif

  ! eri_3center_eigen(_lr) is/are already allocated, then assume that you know what you are doing
  if( ALLOCATED(eri_3center_eigen) .AND. lr_exists ) then
    return
  endif

  if( PRESENT(timing) ) then
    call start_clock(timing)
  else
    call start_clock(timing_eri_3center_eigen)
  endif

  if(PRESENT(verbose)) then
    verbose_ = verbose
  else
    verbose_ = .TRUE.
  endif
  if(verbose_) then
    write(stdout,'(/,a)') ' Calculate 3-center integrals on eigenstates: AO -> MO transform'
  endif

  nbf    = SIZE(c_matrix,DIM=1)
  nstate = SIZE(c_matrix,DIM=2)

  if( PRESENT(mstate_min) ) then
    mstate_min_ = mstate_min
  else
    mstate_min_ = 1
  endif

  if( PRESENT(mstate_max) ) then
    mstate_max_ = mstate_max
  else
    mstate_max_ = nstate
  endif

  if( PRESENT(nstate_min) ) then
    nstate_min_ = nstate_min
  else
    nstate_min_ = 1
  endif

  if( PRESENT(nstate_max) ) then
    nstate_max_ = nstate_max
  else
    nstate_max_ = nstate
  endif

  mstate_count_ = mstate_max_ - mstate_min_ + 1
  nstate_count_ = nstate_max_ - nstate_min_ + 1

  !TODO merge the 2 last indexes to save a factor 2! (i<->j symmetry)
  call clean_allocate('3-center MO integrals',eri_3center_eigen,1,nauxil_local, &
                     mstate_min_,mstate_max_,nstate_min_,nstate_max_,1,nspin,verbose_)
  if(PRESENT(long_range)) then
    call clean_allocate('3-center_lr MO integrals',eri_3center_eigen_lr,1,nauxil_local_lr, &
                     mstate_min_,mstate_max_,nstate_min_,nstate_max_,1,nspin,verbose_)
  endif
  eri_3center_eigen(:,:,:,:) = 0.0_dp

  call clean_allocate('TMP 3-center ints',tmp1,nstate_min_,nstate_max_,1,nbf,verbose_)
  call clean_allocate('TMP 3-center ints',c_t ,nstate_min_,nstate_max_,1,nbf,verbose_)
  call clean_allocate('TMP 3-center ints',tmp2,mstate_min_,mstate_max_,nstate_min_,nstate_max_,verbose_)

  do klspin=1,nspin

    c_t(:,:)  = TRANSPOSE( c_matrix(:,nstate_min_:nstate_max_,klspin) )

    do iauxil=1,nauxil_local
      if( MODULO( iauxil - 1 , ortho%nproc ) /= ortho%rank ) cycle

      tmp1(:,:) = 0.0_dp
      !$OMP PARALLEL PRIVATE(kbf,lbf)
      !$OMP DO REDUCTION(+:tmp1)
      do ipair=1,npair
        kbf = index_basis(1,ipair)
        lbf = index_basis(2,ipair)
        tmp1(:,kbf) = tmp1(:,kbf) +  c_t(:,lbf) * eri_3center(ipair,iauxil)
        tmp1(:,lbf) = tmp1(:,lbf) +  c_t(:,kbf) * eri_3center(ipair,iauxil)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      ! Transformation of the second index
      call DGEMM('T','T',mstate_count_,nstate_count_,nbf, &
                 1.0d0,c_matrix(1,mstate_min_,klspin),nbf, &
                       tmp1(nstate_min_,1),nstate_count_,   &
                 0.0d0,tmp2(mstate_min_,nstate_min_),mstate_count_)

      ! Transposition (ij, P) -> (P, ij) happens here!
      eri_3center_eigen(iauxil,mstate_min_:mstate_max_,nstate_min_:nstate_max_,klspin) &
                                  = tmp2(mstate_min_:mstate_max_,nstate_min_:nstate_max_)

    enddo !iauxil

    if(PRESENT(long_range)) then
      do iauxil=1,nauxil_local_lr
        if( MODULO( iauxil - 1 , ortho%nproc ) /= ortho%rank ) cycle

        tmp1(:,:) = 0.0_dp
        !$OMP PARALLEL PRIVATE(kbf,lbf)
        !$OMP DO REDUCTION(+:tmp1)
        do ipair=1,npair
          kbf = index_basis(1,ipair)
          lbf = index_basis(2,ipair)
          tmp1(:,kbf) = tmp1(:,kbf) +  c_t(:,lbf) * eri_3center_lr(ipair,iauxil)
          tmp1(:,lbf) = tmp1(:,lbf) +  c_t(:,kbf) * eri_3center_lr(ipair,iauxil)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        ! Transformation of the second index
        call DGEMM('T','T',mstate_count_,nstate_count_,nbf, &
                   1.0d0,c_matrix(1,mstate_min_,klspin),nbf, &
                         tmp1(nstate_min_,1),nstate_count_,   &
                   0.0d0,tmp2(mstate_min_,nstate_min_),mstate_count_)

        ! Transposition (ij, P) -> (P, ij) happens here!
        eri_3center_eigen_lr(iauxil,mstate_min_:mstate_max_,nstate_min_:nstate_max_,klspin) &
                                    = tmp2(mstate_min_:mstate_max_,nstate_min_:nstate_max_)

      enddo !iauxil_lr
    endif

  enddo !klspin

  call clean_deallocate('TMP 3-center ints',tmp2,verbose_)
  call clean_deallocate('TMP 3-center ints',c_t,verbose_)
  call clean_deallocate('TMP 3-center ints',tmp1,verbose_)

  call ortho%sum(eri_3center_eigen)
  if( PRESENT(long_range) ) call ortho%sum(eri_3center_eigen_lr)

  if( PRESENT(timing) ) then
    call stop_clock(timing)
  else
    call stop_clock(timing_eri_3center_eigen)
  endif

end subroutine calculate_eri_3center_eigen


!=================================================================
! Calculate LR Coulomb integrals and place them in eri_3center_eigen
subroutine calculate_eri_3center_eigen_lr(c_matrix,mstate_min,mstate_max,nstate_min,nstate_max,timing,verbose)
  implicit none
  real(dp),intent(in)         :: c_matrix(:,:,:)
  integer,optional,intent(in) :: mstate_min,mstate_max,nstate_min,nstate_max
  integer,optional,intent(in) :: timing
  logical,optional,intent(in) :: verbose
  !=====
  logical              :: verbose_
  integer              :: nbf,nstate
  integer              :: mstate_min_,mstate_max_,nstate_min_,nstate_max_
  integer              :: mstate_count_,nstate_count_
  integer              :: kbf,lbf,iauxil,nauxil_local_
  integer              :: lstate
  integer              :: klspin
  real(dp),allocatable :: tmp1(:,:),tmp2(:,:),c_t(:,:)
  integer              :: ipair
  !=====

  ! eri_3center_eigen is already allocated, then assume that you know what you are doing
  if( ALLOCATED(eri_3center_eigen) ) then
    return
  endif

  if( PRESENT(timing) ) then
    call start_clock(timing)
  else
    call start_clock(timing_eri_3center_eigen)
  endif

  if(PRESENT(verbose)) then
    verbose_ = verbose
  else
    verbose_ = .TRUE.
  endif
  if(verbose_) then
    write(stdout,'(/,a)') ' Calculate 3-center integrals on eigenstates: AO -> MO transform'
  endif

  nbf    = SIZE(c_matrix,DIM=1)
  nstate = SIZE(c_matrix,DIM=2)

  if( PRESENT(mstate_min) ) then
    mstate_min_ = mstate_min
  else
    mstate_min_ = 1
  endif

  if( PRESENT(mstate_max) ) then
    mstate_max_ = mstate_max
  else
    mstate_max_ = nstate
  endif

  if( PRESENT(nstate_min) ) then
    nstate_min_ = nstate_min
  else
    nstate_min_ = 1
  endif

  if( PRESENT(nstate_max) ) then
    nstate_max_ = nstate_max
  else
    nstate_max_ = nstate
  endif

  mstate_count_ = mstate_max_ - mstate_min_ + 1
  nstate_count_ = nstate_max_ - nstate_min_ + 1

  nauxil_local_ = nauxil_local_lr
  !TODO merge the 2 last indexes to save a factor 2! (i<->j symmetry)
  call clean_allocate('3-center MO integrals',eri_3center_eigen,1,nauxil_local_, &
                     mstate_min_,mstate_max_,nstate_min_,nstate_max_,1,nspin,verbose_)
  eri_3center_eigen(:,:,:,:) = 0.0_dp

  call clean_allocate('TMP 3-center ints',tmp1,mstate_min_,mstate_max_,1,nbf,verbose_)
  call clean_allocate('TMP 3-center ints',c_t ,mstate_min_,mstate_max_,1,nbf,verbose_)
  call clean_allocate('TMP 3-center ints',tmp2,mstate_min_,mstate_max_,nstate_min_,nstate_max_,verbose_)

  do klspin=1,nspin

    c_t(:,:)  = TRANSPOSE( c_matrix(:,mstate_min_:mstate_max_,klspin) )

    do iauxil=1,nauxil_local_
      if( MODULO( iauxil - 1 , ortho%nproc ) /= ortho%rank ) cycle

      tmp1(:,:) = 0.0_dp
      !$OMP PARALLEL PRIVATE(kbf,lbf)
      !$OMP DO REDUCTION(+:tmp1)
      do ipair=1,npair
        kbf = index_basis(1,ipair)
        lbf = index_basis(2,ipair)
        tmp1(:,kbf) = tmp1(:,kbf) +  c_t(:,lbf) * eri_3center_lr(ipair,iauxil)
        tmp1(:,lbf) = tmp1(:,lbf) +  c_t(:,kbf) * eri_3center_lr(ipair,iauxil)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      ! Transformation of the second index
      call DGEMM('N','N',mstate_count_,nstate_count_,nbf, &
                 1.0d0,tmp1(mstate_min_,1),mstate_count_,   &
                       c_matrix(1,nstate_min_,klspin),nbf, &
                 0.0d0,tmp2(mstate_min_,nstate_min_),mstate_count_)

      ! Transposition happens here!
      eri_3center_eigen(iauxil,mstate_min_:mstate_max_,nstate_min_:nstate_max_,klspin) &
                                  = tmp2(mstate_min_:mstate_max_,nstate_min_:nstate_max_)

    enddo !iauxil

  enddo !klspin

  call clean_deallocate('TMP 3-center ints',tmp2,verbose_)
  call clean_deallocate('TMP 3-center ints',c_t,verbose_)
  call clean_deallocate('TMP 3-center ints',tmp1,verbose_)

  call ortho%sum(eri_3center_eigen)

  if( PRESENT(timing) ) then
    call stop_clock(timing)
  else
    call stop_clock(timing_eri_3center_eigen)
  endif

end subroutine calculate_eri_3center_eigen_lr


!=================================================================
subroutine calculate_eri_3center_eigen_cmplx(c_matrix_cmplx,mstate_min,mstate_max,nstate_min,nstate_max,timing,verbose)
  implicit none
  complex(dp),intent(in)         :: c_matrix_cmplx(:,:,:)
  integer,optional,intent(in) :: mstate_min,mstate_max,nstate_min,nstate_max
  integer,optional,intent(in) :: timing
  logical,optional,intent(in) :: verbose
  !=====
  logical              :: verbose_
  integer              :: nbf,nstate
  integer              :: mstate_min_,mstate_max_,nstate_min_,nstate_max_
  integer              :: mstate_count_,nstate_count_
  integer              :: kbf,lbf,iauxil
  integer              :: lstate
  integer              :: klspin
  complex(dp),allocatable :: tmp1_cmplx(:,:),tmp2_cmplx(:,:),c_t_cmplx(:,:)
  integer              :: ipair
  !=====

  if( PRESENT(timing) ) then
    call start_clock(timing)
  else
    call start_clock(timing_eri_3center_eigen)
  endif

  if(PRESENT(verbose)) then
    verbose_ = verbose
  else
    verbose_ = .TRUE.
  endif
  if(verbose_) then
    write(stdout,'(/,a)') ' Calculate 3-center integrals on eigenstates: AO -> MO transform'
  endif

  nbf    = SIZE(c_matrix_cmplx,DIM=1)
  nstate = SIZE(c_matrix_cmplx,DIM=2)

  if( PRESENT(mstate_min) ) then
    mstate_min_ = mstate_min
  else
    mstate_min_ = 1
  endif

  if( PRESENT(mstate_max) ) then
    mstate_max_ = mstate_max
  else
    mstate_max_ = nstate
  endif

  if( PRESENT(nstate_min) ) then
    nstate_min_ = nstate_min
  else
    nstate_min_ = 1
  endif

  if( PRESENT(nstate_max) ) then
    nstate_max_ = nstate_max
  else
    nstate_max_ = nstate
  endif

  mstate_count_ = mstate_max_ - mstate_min_ + 1
  nstate_count_ = nstate_max_ - nstate_min_ + 1

  !TODO merge the 2 last indexes to save a factor 2! (i<->j symmetry)
  call clean_allocate('3-center MO integrals',eri_3center_eigen_cmplx,1,nauxil_local, &
                      mstate_min_,mstate_max_,nstate_min_,nstate_max_,1,nspin,verbose_)
  eri_3center_eigen_cmplx(:,:,:,:) = (0.0_dp,0.0_dp)

  call clean_allocate('TMP 3-center ints',tmp1_cmplx,mstate_min_,mstate_max_,1,nbf,verbose_)
  call clean_allocate('TMP 3-center ints',c_t_cmplx ,mstate_min_,mstate_max_,1,nbf,verbose_)
  call clean_allocate('TMP 3-center ints',tmp2_cmplx,mstate_min_,mstate_max_,nstate_min_,nstate_max_,verbose_)

  do klspin=1,nspin

    c_t_cmplx(:,:)  = CONJG( TRANSPOSE( c_matrix_cmplx(:,mstate_min_:mstate_max_,klspin) ) )

    do iauxil=1,nauxil_local
      if( MODULO( iauxil - 1 , ortho%nproc ) /= ortho%rank ) cycle

      tmp1_cmplx(:,:) = complex_zero
      !$OMP PARALLEL PRIVATE(kbf,lbf)
      !$OMP DO REDUCTION(+:tmp1_cmplx)
      do ipair=1,npair
        kbf = index_basis(1,ipair)
        lbf = index_basis(2,ipair)
        tmp1_cmplx(:,kbf) = tmp1_cmplx(:,kbf) +  c_t_cmplx(:,lbf) * eri_3center(ipair,iauxil)
        tmp1_cmplx(:,lbf) = tmp1_cmplx(:,lbf) +  c_t_cmplx(:,kbf) * eri_3center(ipair,iauxil)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      ! Transformation of the second index
      call ZGEMM('N','N',mstate_count_,nstate_count_,nbf, &
                 complex_one,tmp1_cmplx(mstate_min_,1),mstate_count_,   &
                       c_matrix_cmplx(1,nstate_min_,klspin),nbf, &
                 complex_zero,tmp2_cmplx(mstate_min_,nstate_min_),mstate_count_)

      ! Transposition happens here!
      eri_3center_eigen_cmplx(iauxil,mstate_min_:mstate_max_,nstate_min_:nstate_max_,klspin) &
                                  = tmp2_cmplx(mstate_min_:mstate_max_,nstate_min_:nstate_max_)

    enddo !iauxil
  enddo !klspin

  call clean_deallocate('TMP 3-center ints',tmp2_cmplx,verbose_)
  call clean_deallocate('TMP 3-center ints',c_t_cmplx,verbose_)
  call clean_deallocate('TMP 3-center ints',tmp1_cmplx,verbose_)

  call ortho%sum(eri_3center_eigen_cmplx)

  if( PRESENT(timing) ) then
    call stop_clock(timing)
  else
    call stop_clock(timing_eri_3center_eigen)
  endif

end subroutine calculate_eri_3center_eigen_cmplx


!=================================================================
subroutine destroy_eri_3center_eigen(verbose,long_range)
  implicit none

  logical,optional,intent(in) :: verbose,long_range
  !=====
  logical :: verbose_
  !=====

  if(PRESENT(verbose)) then
    verbose_ = verbose
  else
    verbose_ = .TRUE.
  endif

  if(verbose_) write(stdout,'(/,a)') ' Destroy 3-center integrals on eigenstates'
  call clean_deallocate('3-center MO integrals',eri_3center_eigen,verbose_)
  if(PRESENT(long_range)) then
    if(long_range) then
      if(verbose_) write(stdout,'(a,/)') ' Destroy 3-center_lr integrals on eigenstates'
      call clean_deallocate('3-center_lr MO integrals',eri_3center_eigen_lr,verbose_)
    endif
  endif

end subroutine destroy_eri_3center_eigen


!=================================================================
subroutine destroy_eri_3center_eigen_cmplx(verbose)
  implicit none

  logical,optional,intent(in) :: verbose
  !=====
  logical :: verbose_
  !=====

  if(PRESENT(verbose)) then
    verbose_ = verbose
  else
    verbose_ = .TRUE.
  endif

  if(verbose_) write(stdout,'(/,a)') ' Destroy 3-center integrals on eigenstates'
  call clean_deallocate('3-center MO integrals',eri_3center_eigen_cmplx,verbose_)

end subroutine destroy_eri_3center_eigen_cmplx


!=================================================================
subroutine form_erimol(nbf,nstate_tot,nstate_jkl,c_matrix,c_matrix_cmplx,ERImol,ERImol_cmplx)
  implicit none

  integer,intent(in)            :: nbf,nstate_tot,nstate_jkl
  real(dp),optional,intent(in)  :: c_matrix(nbf,nstate_tot,nspin)
  real(dp),optional,intent(out) :: ERImol(nstate_tot,nstate_jkl,nstate_jkl,nstate_jkl)
  complex(dp),optional,intent(in)  :: c_matrix_cmplx(nbf,nstate_tot,nspin)
  complex(dp),optional,intent(out) :: ERImol_cmplx(nstate_tot,nstate_jkl,nstate_jkl,nstate_jkl)
  !=====
  integer                 :: pstate,istate,jstate,kstate
  integer                 :: ibf,jbf,abf,bbf
  real(dp),allocatable    :: tmp_pxjx(:,:)
  real(dp),allocatable    :: tmp_pijx(:)
  real(dp),allocatable    :: tmp_pxxx(:,:,:)
  complex(dp),allocatable :: tmp_pxjx_cmplx(:,:)
  complex(dp),allocatable :: tmp_pijx_cmplx(:)
  complex(dp),allocatable :: tmp_pxxx_cmplx(:,:,:)
  !=====
  ! Build <ij|kl> integrals
  !=====
  if(present(c_matrix_cmplx).and.present(ERImol_cmplx)) then

    allocate(tmp_pxxx_cmplx(nbf,nbf,nbf),tmp_pxjx_cmplx(nbf,nbf),tmp_pijx_cmplx(nbf))

    do pstate=1,nstate_tot
      tmp_pxxx_cmplx(:,:,:) = COMPLEX_ZERO
      do bbf=1,nbf
        do jbf=1,nbf
          if( negligible_basispair(jbf,bbf) ) cycle
          do abf=1,nbf
            do ibf=1,nbf
              if( negligible_basispair(ibf,abf) ) cycle
              tmp_pxxx_cmplx(abf,jbf,bbf)=tmp_pxxx_cmplx(abf,jbf,bbf)+conjg(c_matrix_cmplx(ibf,pstate,1))*eri(ibf,abf,jbf,bbf)
            enddo
          enddo
        enddo
      enddo
      do jstate=1,nstate_jkl
        tmp_pxjx_cmplx(:,:) = COMPLEX_ZERO
        do bbf=1,nbf
          do abf=1,nbf
            tmp_pxjx_cmplx(abf,bbf) = SUM( conjg(c_matrix_cmplx(:,jstate,1)) * tmp_pxxx_cmplx(abf,:,bbf) )
          enddo
        enddo
        do istate=1,nstate_jkl
          tmp_pijx_cmplx(:) = COMPLEX_ZERO
          do bbf=1,nbf
            tmp_pijx_cmplx(bbf) = SUM( c_matrix_cmplx(:,istate,1) * tmp_pxjx_cmplx(:,bbf) )
          enddo
          do kstate=1,nstate_jkl
            ERImol_cmplx(pstate,jstate,istate,kstate) = SUM( tmp_pijx_cmplx(:) * c_matrix_cmplx(:,kstate,1) ) ! < pi|jk>
          enddo
        enddo
      enddo
    enddo ! pstate

    deallocate(tmp_pxxx_cmplx,tmp_pxjx_cmplx,tmp_pijx_cmplx)

    elseif(present(c_matrix).and.present(ERImol)) then

    allocate(tmp_pxxx(nbf,nbf,nbf),tmp_pxjx(nbf,nbf),tmp_pijx(nbf))

    do pstate=1,nstate_tot
      tmp_pxxx(:,:,:) = zero
      do bbf=1,nbf
        do jbf=1,nbf
          if( negligible_basispair(jbf,bbf) ) cycle
          do abf=1,nbf
            do ibf=1,nbf
              if( negligible_basispair(ibf,abf) ) cycle
              tmp_pxxx(abf,jbf,bbf)=tmp_pxxx(abf,jbf,bbf)+c_matrix(ibf,pstate,1)*eri(ibf,abf,jbf,bbf)
            enddo
          enddo
        enddo
      enddo
      do jstate=1,nstate_jkl
        tmp_pxjx(:,:) = zero
        do bbf=1,nbf
          do abf=1,nbf
            tmp_pxjx(abf,bbf) = SUM( c_matrix(:,jstate,1) * tmp_pxxx(abf,:,bbf) )
          enddo
        enddo
        do istate=1,nstate_jkl
          tmp_pijx(:) = zero
          do bbf=1,nbf
            tmp_pijx(bbf) = SUM( c_matrix(:,istate,1) * tmp_pxjx(:,bbf) )
          enddo
          do kstate=1,nstate_jkl
            ERImol(pstate,jstate,istate,kstate) = SUM( tmp_pijx(:) * c_matrix(:,kstate,1) ) ! < pi|jk>
          enddo
        enddo
      enddo
    enddo ! pstate

    deallocate(tmp_pxxx,tmp_pxjx,tmp_pijx)

  else
    ! Nth to be done
  endif

end subroutine form_erimol


!=========================================================================
subroutine read_coulombvertex()
  implicit none

  !=====
  integer :: nstate, istate, jstate
  integer :: unitcv, ierr
  complex(dp),allocatable :: coulomb_vertex(:,:,:)
  integer, allocatable :: dims(:)
#if defined(HAVE_MPI)
  integer(kind=MPI_OFFSET_KIND) :: disp
#endif
  !=====


  call yaml_search_keyword('CoulombVertex.yaml', 'length', dims)
  nauxil_global = dims(1)
  nstate        = dims(2)
  write(stdout,*) 'ALLOCATE ',dims(:)
  call distribute_auxil_basis(nauxil_global)

  write(stdout,*) "nauxil_local",nauxil_local
  write(stdout,*) "nstate",nstate
  
  allocate(coulomb_vertex(nauxil_local,nstate,nstate))

#if !defined(HAVE_MPI)
  write(stdout,*) 'Reading file CoulombVertex.elements'
  open(newunit=unitcv, file='CoulombVertex.elements', form='unformatted', access='stream', status='old', action='read')
  read(unitcv) coulomb_vertex(:,:,:)
  close(unitcv)

#else
  write(stdout,*) 'Reading file CoulombVertex.elements with MPI-IO'
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'CoulombVertex.elements', &
                    MPI_MODE_RDONLY, &
                    MPI_INFO_NULL,unitcv,ierr)
 
  disp = 0
  do jstate=1,nstate
    do istate=1,nstate
      iauxil_local = 0
      do ibf_auxil=1,nauxil_global
        disp = disp + 1
        if( auxil%rank == iproc_ibf_auxil(ibf_auxil) ) then
          iauxil_local = iauxil_local + 1
          call MPI_FILE_READ_AT(unitcv, disp, coulomb_vertex(iauxil_local, istate, jstate), &
                                1, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE,ierr)
        endif
        disp = disp + STORAGE_SIZE(coulomb_vertex(1,1,1))
      enddo
    enddo
  enddo
                  
  call MPI_FILE_CLOSE(unitcv, ierr)
#endif

  call clean_allocate('3-center MO integrals',eri_3center_eigen,1,nauxil_local,1,nstate,1,nstate,1,1)

  write(stdout,*) coulomb_vertex(1,1,1)
  write(stdout,*) coulomb_vertex(nauxil_local,nstate,nstate)

end subroutine read_coulombvertex


!=========================================================================
end module m_eri_ao_mo
