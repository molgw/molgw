!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains the calculation of the GW self-energy
! within different flavors: G0W0, GnW0, GnWn, COHSEX, QSGW
!
!=========================================================================
#include "molgw.h"
module m_gw_selfenergy_analytic
  use m_definitions
  use m_mpi
  use m_timing
  use m_inputparam
  use m_warning, only: issue_warning
  use m_spectral_function
  use m_eri_ao_mo
  use m_selfenergy_tools
  use m_scalapack
  use m_linear_algebra, only: diagonalize


contains


!=========================================================================
subroutine gw_selfenergy(selfenergy_approx, occupation, energy, c_matrix, wpol, se)
  implicit none

  integer, intent(in)                  :: selfenergy_approx
  real(dp), intent(in)                 :: occupation(:, :), energy(:, :)
  real(dp), intent(in)                 :: c_matrix(:, :, :)
  type(spectral_function), intent(in)  :: wpol
  type(selfenergy_grid), intent(inout) :: se
  !=====
  integer               :: nstate
  integer               :: iomega
  integer               :: ipstate
  integer               :: pstate
  integer               :: istate, ispin, spole
  real(dp), allocatable  :: w_s(:, :)
  real(dp)              :: fact_full_i, fact_empty_i
  real(dp)              :: energy_gw
  !=====

  call start_clock(timing_gw_self)

  nstate = SIZE(occupation, DIM=1)

  write(stdout, *)
  select case(selfenergy_approx)
  case(GW)
    write(stdout, *) 'Perform a one-shot G0W0 calculation'
  case(ONE_RING)
    write(stdout, *) 'Perform a one-shot one-ring calculation'
  case(COHSEX)
    write(stdout, *) 'Perform a COHSEX calculation'
  case(GnW0)
    write(stdout, *) 'Perform an eigenvalue self-consistent GnW0 calculation'
  case(GnWn)
    write(stdout, *) 'Perform an eigenvalue self-consistent GnWn calculation'
  case default
    write(stdout, *) 'type:', selfenergy_approx
    call die('gw_selfenergy: calculation type unknown')
  end select


  if(has_auxil_basis) then
    call calculate_eri_3center_eigen(c_matrix, nsemin, nsemax, ncore_G+1, nvirtual_G-1, timing=timing_aomo_gw)
  endif


  call clean_allocate('Temporary array', w_s, 1, wpol%npole_reso, nsemin, nsemax)


  energy_gw = 0.0_dp
  se%sigma(:, :, :) = 0.0_dp

  do ispin=1, nspin
    do istate=ncore_G+1, nvirtual_G-1 !INNER LOOP of G

      if( MODULO( istate - (ncore_G+1) , poorman%nproc) /= poorman%rank ) cycle

      !
      ! Prepare w_s with the knowledge of index istate and pstate
      if( .NOT. has_auxil_basis) then
        !$OMP PARALLEL
        !$OMP DO PRIVATE(ipstate)
        ! Here just grab the precalculated value
        do pstate=nsemin, nsemax
          ipstate = index_prodstate(istate, pstate) + (ispin-1) * index_prodstate(nvirtual_W-1, nvirtual_W-1)
          w_s(:, pstate) = wpol%residue_left(ipstate, :)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
      else
        ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
        w_s(:, nsemin:nsemax)     = MATMUL( TRANSPOSE(wpol%residue_left(:, :)) , &
                                            eri_3center_eigen(:, nsemin:nsemax, istate, ispin) )
        call auxil%sum(w_s)
      endif



      ! The application of residue theorem only retains the pole in given
      ! quadrants.
      ! The positive poles of W go with the poles of occupied states in G
      ! The negative poles of W go with the poles of empty states in G
      fact_full_i   = occupation(istate, ispin) / spin_fact
      fact_empty_i = (spin_fact - occupation(istate, ispin)) / spin_fact

      do spole=1, wpol%npole_reso

        if( MODULO( spole - 1 , auxil%nproc ) /= auxil%rank ) cycle


        select case(selfenergy_approx)

        case(GW, GnW0, GnWn, ONE_RING)
          !
          ! calculate only the diagonal !
          do pstate=nsemin, nsemax
            !$OMP PARALLEL
            !$OMP DO
            do iomega=-se%nomega, se%nomega
              se%sigma(iomega, pstate, ispin) = se%sigma(iomega, pstate, ispin) &
                      + w_s(spole, pstate) * w_s(spole, pstate)                                          &
                        * ( fact_full_i  / ( se%energy0(pstate, ispin) + se%omega(iomega) &
                                             - energy(istate, ispin) + wpol%pole(spole) - ieta )  &
                          + fact_empty_i / ( se%energy0(pstate, ispin) + se%omega(iomega) &
                                             - energy(istate, ispin) - wpol%pole(spole) + ieta ) )
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            if( (spin_fact - occupation(pstate, ispin)) / spin_fact < completely_empty) then
              energy_gw = energy_gw + fact_empty_i * occupation(pstate, ispin) &
                               * w_s(spole, pstate)**2 / ( energy(pstate, ispin) - energy(istate, ispin) - wpol%pole(spole) )
            endif
          enddo
        case(COHSEX)
          !$OMP PARALLEL
          !$OMP DO PRIVATE(pstate)
          do pstate=nsemin, nsemax
            !
            ! SEX
            !
            se%sigma(0, pstate, ispin) = se%sigma(0, pstate, ispin) &
                      + w_s(spole, pstate) * w_s(spole, pstate) &
                            * fact_full_i / wpol%pole(spole) * 2.0_dp

            !
            ! COH
            !
            se%sigma(0, pstate, ispin) = se%sigma(0, pstate, ispin) &
                      - w_s(spole, pstate) * w_s(spole, pstate) &
                            / wpol%pole(spole)

          enddo
          !$OMP END DO
          !$OMP END PARALLEL
        case default
          call die('BUG')
        end select

      enddo !spole

    enddo !istate
  enddo !ispin

  ! Sum up the contribution from different poles (= different procs)
  call world%sum(se%sigma)
  call world%sum(energy_gw)


  write(stdout, '(a)') ' Sigma_c(omega) is calculated'

  if( nsemax >= nhomo_G .AND. nsemin <= ncore_G+1 ) &
      write(stdout, '(1x,a,1x,f19.10)') 'Correlation energy 1/2 Tr[ Sig_c * G ] (Ha):', energy_gw


  call clean_deallocate('Temporary array', w_s)

  if(has_auxil_basis) then
    call destroy_eri_3center_eigen()
  endif


  call stop_clock(timing_gw_self)


end subroutine gw_selfenergy


!=========================================================================
subroutine gw_selfenergy_upfolding(selfenergy_approx, occupation, energy, c_matrix, wpol, exchange_m_vxc)
  implicit none

  integer, intent(in)                  :: selfenergy_approx
  real(dp), intent(in)                 :: occupation(:, :), energy(:, :)
  real(dp), intent(in)                 :: c_matrix(:, :, :), exchange_m_vxc(:, :, :)
  type(spectral_function), intent(in)  :: wpol
  !=====
  character(len=4)     :: ctmp
  integer              :: nstate
  integer              :: ipstate
  integer              :: pstate
  integer              :: istate, ispin
  real(dp)             :: sign_i, mu
  real(dp)             :: weight
  real(dp), allocatable :: matrix_wing(:, :), matrix_head(:, :), matrix_diag(:)
  real(dp), allocatable :: super_matrix(:, :), eigval(:)
  integer              :: nmat, mwing, imat, jmat
  integer              :: mstate, jstate
  integer              :: irecord
  integer              :: fu, info
  integer              :: desc_wing(NDEL), desc_eri(NDEL), desc_wpol(NDEL)
  integer              :: nauxil_local, ilocal, mwing_local
#if defined(HAVE_SCALAPACK)
  integer              :: desc_matrix(NDEL)
  real(dp)             :: work(1), nelect, rtmp
  integer              :: nlocal, jlocal, iglobal, jglobal
#endif
  !=====

  call start_clock(timing_gw_self)

  nstate = SIZE(energy, DIM=1)

  write(stdout, *)
  select case(selfenergy_approx)
  case(GW)
    write(stdout, *) 'Perform a one-shot G0W0 calculation with full solution of the Dyson equation'
  case default
    write(stdout, *) 'type:', selfenergy_approx
    call die('gw_selfenergy_upfolding: calculation type unknown')
  end select

  if( nspin > 1 ) call die('gw_selfenergy_upfolding: not functional for nspin>1')

  if( nsemin /= ncore_G+1 .OR. nsemax /= nvirtual_G-1 ) then
    write(stdout, '(1x,a,i5,1x,i5)') 'nsemin ?= ncore_G+1    ', nsemin, ncore_G+1
    write(stdout, '(1x,a,i5,1x,i5)') 'nsemax ?= nvirtual_G-1 ', nsemax, nvirtual_G-1
    call die('gw_selfenergy_upfolding: selfenergy state range should contain all the active states')
  endif

  if(has_auxil_basis) then
    call calculate_eri_3center_eigen(c_matrix, nsemin, nsemax, ncore_G+1, nvirtual_G-1, timing=timing_aomo_gw)
  endif

  mstate = nvirtual_G - ncore_G - 1
  nmat   = mstate * ( 1 + wpol%npole_reso)
  mwing  = mstate * wpol%npole_reso

  !
  ! Temporary descriptors
  ! desc_wpol for wpol%residue_left
  nauxil_local = NUMROC(nauxil_global, MB_eri3_mo, iprow_eri3_mo, first_row, nprow_eri3_mo)
  call DESCINIT(desc_wpol, nauxil_global, wpol%npole_reso, MB_eri3_mo, NB_eri3_mo, first_row, first_col, &
                cntxt_eri3_mo, MAX(1, nauxil_local), info)
  ! desc_eri for wpol%residue_left
  call DESCINIT(desc_eri, nauxil_global, mstate, MB_eri3_mo, NB_eri3_mo, first_row, first_col, &
                cntxt_eri3_mo, MAX(1, nauxil_local), info)

  ! desc_wing for matrix_wing
  mwing_local = NUMROC(mwing, MB_eri3_mo, iprow_eri3_mo, first_row, nprow_eri3_mo)
  call DESCINIT(desc_wing, mwing, mstate, MB_eri3_mo, NB_eri3_mo, first_row, first_col, &
                cntxt_eri3_mo, MAX(1, mwing_local), info)


  call clean_allocate('Matrix head', matrix_head, mstate, mstate)
  call clean_allocate('Matrix wing', matrix_wing, mwing_local, mstate)


  allocate(matrix_diag(mwing))
  matrix_head(:, :) = 0.0_dp
  matrix_wing(:, :) = 0.0_dp
  matrix_diag(:)    = 0.0_dp

  matrix_head(:, :) = exchange_m_vxc(ncore_G+1:nvirtual_G-1, ncore_G+1:nvirtual_G-1, 1)  ! spin index set to 1

  write(stdout, '(/,1x,a,i8,a,i8)') 'Diagonalization problem of size: ', nmat, ' x ', nmat

  do ispin=1, nspin
    do istate=ncore_G+1, nvirtual_G-1 !INNER LOOP of G

      if( MODULO( istate - (ncore_G+1) , poorman%nproc) /= poorman%rank ) cycle
      !
      ! indeces
      jstate = istate - ncore_G
      sign_i = MERGE(-1.0_dp, 1.0_dp, occupation(istate, ispin) / spin_fact > completely_empty )
      irecord = ( jstate - 1 ) * wpol%npole_reso

      !
      ! Head
      matrix_head(jstate, jstate) = matrix_head(jstate, jstate) + energy(istate, ispin)

      !
      ! Diagonal
      matrix_diag(irecord+1:irecord+wpol%npole_reso)   = energy(istate, ispin) + sign_i * wpol%pole(:)

      !
      ! Wing
      if( .NOT. has_auxil_basis) then
        do pstate=nsemin, nsemax
          ipstate = index_prodstate(istate, pstate) + (ispin-1) * index_prodstate(nvirtual_W-1, nvirtual_W-1)
          matrix_wing(irecord+1:irecord+wpol%npole_reso, pstate-ncore_G) = wpol%residue_left(ipstate, :)
        enddo
      else
        ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
#if defined(HAVE_SCALAPACK)
        call PDGEMM('T', 'N', wpol%npole_reso, mstate,nauxil_global, &
                   1.0_dp, wpol%residue_left, 1, 1, desc_wpol,        &
                          eri_3center_eigen(1, nsemin, istate, ispin), 1, 1, desc_eri,        &
                   0.0_dp, matrix_wing, irecord+1, 1, desc_wing)
#else
        matrix_wing(irecord+1:irecord+wpol%npole_reso, :) = &
             MATMUL( TRANSPOSE(wpol%residue_left(:, :)) , eri_3center_eigen(:, nsemin:nsemax, istate, ispin) )
        call auxil%sum(matrix_wing(irecord+1:irecord+wpol%npole_reso, :))
#endif
      endif


    enddo !istate
  enddo !ispin

  write(stdout, '(a)') ' Matrix is setup'

  !
  ! Dump the super_matrix on files (1 file per SCALAPACK thread)
  write(stdout, *) 'Dump the large sparse super_matrix on disk'
  write(ctmp, '(i4.4)') world%rank
  open(newunit=fu, file='MATRIX_'//ctmp, form='formatted', action='write')

  ! only master writes the head and the long diagonal
  if( is_iomaster) then
    do jmat=1, mstate
      do imat=1, jmat
        write(fu, '(1x,i7,1x,i7,1x,e16.8)') imat, jmat, matrix_head(imat, jmat)*Ha_eV
      enddo
    enddo
    do imat=1, mwing
      write(fu, '(1x,i7,1x,i7,1x,e16.8)') imat+mstate, imat+mstate, matrix_diag(imat)*Ha_eV
    enddo
  endif

  do jmat=1, mstate
    do ilocal=1, mwing_local
      imat = INDXL2G(ilocal, MB_eri3_mo, iprow_eri3_mo, first_row, nprow_eri3_mo)
      write(fu, '(1x,i7,1x,i7,1x,e16.8)') mstate+imat, jmat, matrix_wing(ilocal, jmat)*Ha_eV
    enddo
  enddo
  close(fu)


  !
  ! If the super_matrix is small enough, then diagonalize it!
  if( nmat < 13001 ) then

    mu = ( MINVAL(energy(nhomo_G+1, :)) + MAXVAL(energy(nhomo_G, :)) ) / 2.0_dp
    write(stdout, '(1x,a,f12.6)') 'Center of the HOMO-LUMO gap (eV): ', mu*Ha_eV
    allocate(eigval(nmat))

#if defined(HAVE_SCALAPACK)
    write(stdout, '(1x,a,i4,a,i4)') 'Diagonalize the large sparse super_matrix as if ' // &
                                    'it were dense with SCALAPACK with distribution: ', &
                                  nprow_sd, ' x ', npcol_sd
    mwing_local = NUMROC(nmat, block_row, iprow_sd, first_row, nprow_sd)
    nlocal = NUMROC(nmat, block_col, ipcol_sd, first_col, npcol_sd)
    call DESCINIT(desc_matrix, nmat, nmat, block_row, block_col, first_row, first_col, cntxt_sd, MAX(1, mwing_local), info)
    call clean_allocate('Super matrix', super_matrix, mwing_local, nlocal)
    super_matrix(:, :) = 0.0_dp

    write(stdout, *) 'Start copy wing block'
    call PDGEMR2D(mwing, mstate, matrix_wing, 1, 1, desc_wing, super_matrix, mstate+1, 1, desc_matrix, cntxt_sd)
    write(stdout, *) 'copy done'

    write(stdout, *) 'Set head block'
    do jglobal=1, mstate
      if( INDXG2P(jglobal, block_col, ipcol_sd, first_col, npcol_sd) /= ipcol_sd ) cycle
      jlocal = INDXG2L(jglobal, block_col, ipcol_sd, first_col, npcol_sd)
      do iglobal=1, mstate
        if( INDXG2P(iglobal, block_row, iprow_sd, first_row, nprow_sd) /= iprow_sd ) cycle
        ilocal = INDXG2L(iglobal, block_row, iprow_sd, first_row, nprow_sd)
        super_matrix(ilocal, jlocal) = matrix_head(iglobal, jglobal)
      enddo
    enddo

    write(stdout, *) 'Set large diagonal'
    do iglobal=mstate+1, nmat
      jglobal=iglobal
      if( INDXG2P(iglobal, block_row, iprow_sd, first_row, nprow_sd) /= iprow_sd ) cycle
      if( INDXG2P(jglobal, block_col, ipcol_sd, first_col, npcol_sd) /= ipcol_sd ) cycle
      jlocal = INDXG2L(jglobal, block_col, ipcol_sd, first_col, npcol_sd)
      ilocal = INDXG2L(iglobal, block_row, iprow_sd, first_row, nprow_sd)
      super_matrix(ilocal, jlocal) = matrix_diag(iglobal-mstate)
    enddo

    write(stdout, *) 'Start diago'
    call diagonalize_sca('R', super_matrix, desc_matrix,eigval)
    write(stdout, *) 'Diago done'

    write(stdout, *) '============== Poles in eV , weight ==============='
    open(newunit=fu, file='GREENS_FUNCTION', action='write')
    nelect = 0.0_dp
    do jmat=1, nmat
      weight = PDLANGE('F', mstate, 1, super_matrix, 1, jmat, desc_matrix,work)**2
      ! If eigenvalue lower than the middle of the HOMO-LUMO gap,
      ! then consider the excitation is occupied
      if( eigval(jmat) < mu ) nelect = nelect + spin_fact * weight
      if( weight > 5.0e-2_dp ) then
        call PDAMAX(mstate, rtmp, jstate, super_matrix, 1, jmat, desc_matrix, 1)
        call world%max(jstate)
        write(stdout, '(1x,a,i5.5,a,f16.6,4x,f12.6)') 'Projection on state ', jstate, ': ', eigval(jmat)*Ha_eV, weight
      endif
      write(fu, '(1x,f16.6,4x,f12.6)') eigval(jmat)*Ha_eV, weight
    enddo
    close(fu)
    write(stdout, '(1x,a,f12.6)') 'Number of electrons: ', nelect
    write(stdout, *) '==================================================='

#else
    write(stdout, *) 'Diagonalize the large sparse super_matrix as if it were dense'
    call clean_allocate('Super matrix', super_matrix, nmat, nmat)
    super_matrix(:, :)                    = 0.0_dp
    super_matrix(1:mstate, 1:mstate)      = matrix_head(:, :)
    super_matrix(1:mstate, mstate+1:nmat) = TRANSPOSE(matrix_wing(:, :))
    super_matrix(mstate+1:nmat, 1:mstate) = matrix_wing(:, :)
    do imat=mstate+1, nmat
      super_matrix(imat, imat) = matrix_diag(imat-mstate)
    enddo
    call diagonalize('R', super_matrix, eigval)

    write(stdout, '(1x,a,i8)') 'Number of non-negligible poles: ', COUNT( SUM(super_matrix(1:mstate, :)**2, DIM=1) > 1.0e-3_dp )
    write(stdout, *) '============== Poles in eV , weight ==============='
    open(newunit=fu, file='GREENS_FUNCTION', action='write')
    do jmat=1, nmat
      weight = SUM(super_matrix(1:mstate, jmat)**2)
      if( weight > 5.0e-2_dp ) then
        jstate = MAXLOC(ABS(super_matrix(1:mstate, jmat)), DIM=1)
        write(stdout, '(1x,a,i5.5,a,f16.6,4x,f12.6)') 'Projection on state ', jstate, ': ', eigval(jmat)*Ha_eV, weight
      endif
      write(fu, '(1x,f16.6,4x,f12.6)') eigval(jmat)*Ha_eV, SUM(super_matrix(1:mstate, jmat)**2)
    enddo
    close(fu)
    ! If eigenvalue lower than the middle of the HOMO-LUMO gap,
    ! then consider the excitation is occupied
    write(stdout, '(1x,a,f12.6)') 'Number of electrons: ', & 
             spin_fact * SUM( SUM(super_matrix(1:mstate, :)**2, DIM=1), MASK=(eigval(:) < mu) )
    write(stdout, *) '==================================================='

#endif
    call clean_deallocate('Super matrix', super_matrix)
    deallocate(eigval)

  endif

  call clean_deallocate('Matrix wing', matrix_wing)
  call clean_deallocate('Matrix head', matrix_head)
  deallocate(matrix_diag)
  if(has_auxil_basis) then
    call destroy_eri_3center_eigen()
  endif

  call stop_clock(timing_gw_self)

end subroutine gw_selfenergy_upfolding


!=========================================================================
subroutine gw_selfenergy_scalapack(selfenergy_approx, occupation, energy, c_matrix, wpol, se)
  implicit none

  integer, intent(in)                  :: selfenergy_approx
  real(dp), intent(in)                 :: occupation(:, :), energy(:, :)
  real(dp), intent(in)                 :: c_matrix(:, :, :)
  type(spectral_function), intent(in)  :: wpol
  type(selfenergy_grid), intent(inout) :: se
  !=====
#if defined(HAVE_SCALAPACK)
  integer                 :: pstate, pspin, nstate
  integer                 :: iomega
  integer                 :: istate, spole
  real(dp)                :: fact_full_i, fact_empty_i
  integer                 :: desc_wauxil(NDEL), desc_wsd(NDEL)
  integer                 :: desc_3auxil(NDEL), desc_3sd(NDEL)
  integer                 :: desc_w_s(NDEL)
  integer                 :: mlocal, nlocal
  integer                 :: ilocal, jlocal, jglobal
  integer                 :: info
  real(dp), allocatable    :: eri_3tmp_auxil(:, :), eri_3tmp_sd(:, :)
  real(dp), allocatable    :: wresidue_sd(:, :)
  real(dp), allocatable    :: w_s(:, :)
  complex(dp), allocatable :: sigmagw(:, :, :)
#endif
  !=====

  if(.NOT. has_auxil_basis) return


#if defined(HAVE_SCALAPACK)
  nstate = SIZE(energy, DIM=1)

  call start_clock(timing_gw_self)

  write(stdout, *)
  select case(selfenergy_approx)
  case(ONE_RING)
    write(stdout, *) 'Perform a one-shot one_ring calculation: SCALAPACK'
  case(GW)
    write(stdout, *) 'Perform a one-shot G0W0 calculation: SCALAPACK'
  case(GnW0)
    write(stdout, *) 'Perform an eigenvalue self-consistent GnW0 calculation: SCALAPACK'
  case(GnWn)
    write(stdout, *) 'Perform an eigenvalue self-consistent GnWn calculation: SCALAPACK'
  case default
    write(stdout, *) 'type:', selfenergy_approx
    call die('gw_selfenergy_scalapack: calculation type unknown')
  end select


  call calculate_eri_3center_eigen(c_matrix, ncore_G+1, nvirtual_G-1, nsemin, nsemax, timing=timing_aomo_gw)



  !
  ! SCALAPACK preparation for W
  !  wpol%residue_left
  mlocal = NUMROC(nauxil_global , MB_eri3_mo, iprow_eri3_mo, first_row, nprow_eri3_mo)
  nlocal = NUMROC(wpol%npole_reso, NB_eri3_mo, ipcol_eri3_mo, first_col, npcol_eri3_mo)
  call DESCINIT(desc_wauxil, nauxil_global, wpol%npole_reso, MB_eri3_mo, NB_eri3_mo, &
               first_row, first_col, cntxt_eri3_mo, MAX(1, mlocal), info)
  !
  ! Change data distribution
  ! from cntxt_eri3_mo to cntxt_sd
  mlocal = NUMROC(nauxil_global , block_row, iprow_sd, first_row, nprow_sd)
  nlocal = NUMROC(wpol%npole_reso, block_col, ipcol_sd, first_col, npcol_sd)
  call DESCINIT(desc_wsd, nauxil_global, wpol%npole_reso, block_row, block_col, first_row, first_col, &
                cntxt_sd, MAX(1, mlocal), info)
  call clean_allocate('TMP distributed W', wresidue_sd, mlocal, nlocal)
  call PDGEMR2D(nauxil_global, wpol%npole_reso, wpol%residue_left, 1, 1, desc_wauxil, &
                                                    wresidue_sd, 1, 1, desc_wsd, cntxt_sd)

  ! Temporary array sigmagw is created because OPENMP does not want to work directly with se%sigma
  allocate(sigmagw, MOLD=se%sigma)
  sigmagw(:, :, :)  = 0.0_dp

  do pspin=1, nspin
    do pstate=nsemin, nsemax


      !
      ! SCALAPACK preparation for the 3-center integrals
      !
      mlocal = NUMROC(nauxil_global      , MB_eri3_mo, iprow_eri3_mo, first_row, nprow_eri3_mo)
      nlocal = NUMROC(nvirtual_G-ncore_G-1, NB_eri3_mo, ipcol_eri3_mo, first_col, npcol_eri3_mo)
      call DESCINIT(desc_3auxil, nauxil_global, nvirtual_G-ncore_G-1, MB_eri3_mo, NB_eri3_mo, &
                   first_row, first_col, cntxt_eri3_mo, MAX(1, mlocal), info)

      if( cntxt_eri3_mo > 0 ) then
        call clean_allocate('TMP 3center eigen', eri_3tmp_auxil, mlocal, nlocal, verbose=.FALSE.)
        do jlocal=1, nlocal
          jglobal = INDXL2G(jlocal, NB_eri3_mo, ipcol_eri3_mo, first_col, npcol_eri3_mo) + ncore_G
          do ilocal=1, mlocal
            eri_3tmp_auxil(ilocal, jlocal) = eri_3center_eigen(ilocal, jglobal, pstate, pspin)
          enddo
        enddo
      else
        call clean_allocate('TMP 3center eigen', eri_3tmp_auxil, 1, 1, verbose=.FALSE.)
      endif
      !
      ! Change data distribution
      ! from cntxt_eri3_mo to cntxt_sd
      mlocal = NUMROC(nauxil_global      , block_row, iprow_sd, first_row, nprow_sd)
      nlocal = NUMROC(nvirtual_G-ncore_G-1, block_col, ipcol_sd, first_col, npcol_sd)
      call DESCINIT(desc_3sd, nauxil_global, nvirtual_G-ncore_G-1, block_row, block_col, &
                   first_row, first_col, cntxt_sd, MAX(1, mlocal), info)
      call clean_allocate('TMP 3center eigen', eri_3tmp_sd, mlocal, nlocal, verbose=.FALSE.)
      call PDGEMR2D(nauxil_global, nvirtual_G-ncore_G-1, eri_3tmp_auxil, 1, 1, desc_3auxil, &
                                                          eri_3tmp_sd, 1, 1, desc_3sd, cntxt_sd)
      call clean_deallocate('TMP 3center eigen', eri_3tmp_auxil, verbose=.FALSE.)


      !
      ! Prepare a SCALAPACKed w_s that is to contain  wresidue**T * v**1/2
      mlocal = NUMROC(wpol%npole_reso     , block_row, iprow_sd, first_row, nprow_sd)
      nlocal = NUMROC(nvirtual_G-ncore_G-1, block_col, ipcol_sd, first_col, npcol_sd)
      call DESCINIT(desc_w_s, wpol%npole_reso, nvirtual_G-ncore_G-1, block_row, block_col, &
                   first_row, first_col, cntxt_sd, MAX(1, mlocal), info)
      call clean_allocate('Temporary array', w_s, mlocal, nlocal, verbose=.FALSE.)

      ! And calculate it
      call PDGEMM('T', 'N', wpol%npole_reso, nvirtual_G-ncore_G-1,nauxil_global, &
                             1.0_dp, wresidue_sd, 1, 1, desc_wsd,    &
                                    eri_3tmp_sd, 1, 1, desc_3sd,    &
                             0.0_dp, w_s        , 1, 1, desc_w_s)
      call clean_deallocate('TMP 3center eigen', eri_3tmp_sd, verbose=.FALSE.)



      !$OMP PARALLEL PRIVATE(istate,spole,fact_full_i,fact_empty_i)
      !$OMP DO REDUCTION(+:sigmagw)
      do jlocal=1, nlocal
        istate = INDXL2G(jlocal, block_col, ipcol_sd, first_col, npcol_sd) + ncore_G
        do ilocal=1, mlocal
          spole = INDXL2G(ilocal, block_row, iprow_sd, first_row, nprow_sd)


          ! The application of residue theorem only retains the pole in given
          ! quadrants.
          ! The positive poles of W go with the poles of occupied states in G
          ! The negative poles of W go with the poles of empty states in G
          fact_full_i   = occupation(istate, pspin) / spin_fact
          fact_empty_i = (spin_fact - occupation(istate, pspin)) / spin_fact

          sigmagw(:, pstate, pspin) = sigmagw(:, pstate, pspin) &
                  + w_s(ilocal, jlocal) * w_s(ilocal, jlocal)                    &
                    * ( fact_full_i  / ( se%energy0(pstate, pspin) + se%omega(:) &
                                        - energy(istate, pspin) + wpol%pole(spole) - ieta )   &
                      + fact_empty_i / ( se%energy0(pstate, pspin) + se%omega(:) &
                                        - energy(istate, pspin) - wpol%pole(spole) + ieta )  )
        enddo  !ilocal -> spole
      enddo !jlocal -> istate
      !$OMP END DO
      !$OMP END PARALLEL

      call clean_deallocate('Temporary array', w_s, verbose=.FALSE.)

    enddo !pstate
  enddo !pspin

  ! Sum up the contribution from different poles (= different procs)
  call world%sum(sigmagw)

  se%sigma(:, :, :) = sigmagw(:, :, :)
  deallocate(sigmagw)

  write(stdout, '(a)') ' Sigma_c(omega) is calculated'

  call clean_deallocate('TMP distributed W', wresidue_sd)
  call destroy_eri_3center_eigen()

  call stop_clock(timing_gw_self)

#endif

end subroutine gw_selfenergy_scalapack


!=========================================================================
subroutine gw_selfenergy_qs(occupation, energy, c_matrix, s_matrix, wpol, selfenergy_ao)
  implicit none

  real(dp), intent(in)                :: occupation(:, :), energy(:, :)
  real(dp), intent(in)                :: c_matrix(:, :, :)
  real(dp), intent(in)                :: s_matrix(:, :)
  type(spectral_function), intent(in) :: wpol
  real(dp), intent(out)               :: selfenergy_ao(:, :, :)
  !=====
  integer               :: nstate
  integer               :: ipstate, pstate, qstate, istate
  integer               :: ispin, spole
  real(dp), allocatable  :: w_s(:, :), selfenergy_mo(:, :, :)
  real(dp)              :: fact_full_i, fact_empty_i
  !=====

  call start_clock(timing_gw_self)

  nstate = SIZE(energy, DIM=1)

  write(stdout, *)
  select case(calc_type%selfenergy_approx)
  case(ONE_RING)
    write(stdout, *) 'Perform a QP self-consistent one-ring calculation (QS1-ring)'
  case(GW)
    write(stdout, *) 'Perform a QP self-consistent GW calculation (QSGW)'
  case(COHSEX)
    write(stdout, *) 'Perform a self-consistent COHSEX calculation'
  case default
    call die('gw_selfenergy_qs: calculation type unknown')
  end select


  if(has_auxil_basis) then
    call calculate_eri_3center_eigen(c_matrix, nsemin, nsemax, ncore_G+1, nvirtual_G-1, timing=timing_aomo_gw)
  endif

  call clean_allocate('Temporary array', w_s, 1, wpol%npole_reso, nsemin, nsemax)

  allocate(selfenergy_mo(nstate, nstate, nspin))

  selfenergy_mo(:, :, :)  = 0.0_dp
  do ispin=1, nspin
    do istate=ncore_G+1, nvirtual_G-1 !INNER LOOP of G

      if( MODULO( istate - (ncore_G+1) , poorman%nproc) /= poorman%rank ) cycle

      !
      ! Prepare w_s with the knowledge of index istate and pstate
      if( .NOT. has_auxil_basis) then
        !$OMP PARALLEL
        !$OMP DO PRIVATE(ipstate)
        ! Here just grab the precalculated value
        do pstate=nsemin, nsemax
          ipstate = index_prodstate(istate, pstate) + (ispin-1) * index_prodstate(nvirtual_W-1, nvirtual_W-1)
          w_s(:, pstate) = wpol%residue_left(ipstate, :)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
      else
        ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
        w_s(:, nsemin:nsemax)     = MATMUL( TRANSPOSE(wpol%residue_left(:, :)), &
                                            eri_3center_eigen(:, nsemin:nsemax, istate, ispin) )
        call auxil%sum(w_s)
      endif


      ! The application of residue theorem only retains the pole in given
      ! quadrants.
      ! The positive poles of W go with the poles of occupied states in G
      ! The negative poles of W go with the poles of empty states in G
      fact_full_i   = occupation(istate, ispin) / spin_fact
      fact_empty_i = (spin_fact - occupation(istate, ispin)) / spin_fact


      do spole=1, wpol%npole_reso

        if( MODULO( spole - 1 , auxil%nproc ) /= auxil%rank ) cycle

        select case(calc_type%selfenergy_approx)

        case(GW)
          !$OMP PARALLEL
          !$OMP DO COLLAPSE(2)
          do qstate=nsemin, nsemax
            do pstate=nsemin, nsemax

              selfenergy_mo(pstate, qstate, ispin) = selfenergy_mo(pstate, qstate, ispin) &
                        + w_s(spole, pstate) * w_s(spole, qstate)                            &
                          * ( REAL(  fact_full_i  / ( energy(qstate, ispin) - ieta  &
                                                     - energy(istate, ispin) + wpol%pole(spole) ) , dp ) &
                            + REAL(  fact_empty_i / ( energy(qstate, ispin) + ieta  &
                                                     - energy(istate, ispin) - wpol%pole(spole) ) , dp ) )

            enddo
          enddo
          !$OMP END DO
          !$OMP END PARALLEL
        case(COHSEX)
          !$OMP PARALLEL
          !$OMP DO COLLAPSE(2)
          do qstate=nsemin, nsemax
            do pstate=nsemin, nsemax
              !
              ! SEX
              !
              selfenergy_mo(pstate, qstate, ispin) = selfenergy_mo(pstate, qstate, ispin) &
                        + w_s(spole, pstate) * w_s(spole, qstate)                                &
                              * fact_full_i / wpol%pole(spole) * 2.0_dp

              !
              ! COH
              !
              selfenergy_mo(pstate, qstate, ispin) = selfenergy_mo(pstate, qstate, ispin) &
                        - w_s(spole, pstate) * w_s(spole, qstate) &
                              / wpol%pole(spole)
            enddo
          enddo
          !$OMP END DO
          !$OMP END PARALLEL
        case default
          call die('BUG')
        end select

      enddo !spole

    enddo !istate
  enddo !ispin

  ! Sum up the contribution from different poles (= different procs)
  call world%sum(selfenergy_mo)

  call clean_deallocate('Temporary array', w_s)


  ! Kotani's hermitianization trick
  call apply_qs_approximation(s_matrix, c_matrix, selfenergy_mo, selfenergy_ao)


  if(has_auxil_basis) call destroy_eri_3center_eigen()


  call stop_clock(timing_gw_self)


end subroutine gw_selfenergy_qs


!=========================================================================
subroutine dump_gw_ingredients(energy, c_matrix, wpol)
  implicit none

  real(dp), intent(in)                 :: energy(:, :)
  real(dp), intent(in)                 :: c_matrix(:, :, :)
  type(spectral_function), intent(in)  :: wpol
  !=====
  integer               :: qstate, qspin, spole
  real(dp), allocatable :: wcoeff(:, :)
  integer               :: file_v, file_w, file_e, file_omega
  !=====

  if(.NOT. has_auxil_basis) return
  if( world%nproc > 1 ) call die('dump_gw_ingredients: not implemented with MPI')

  write(stdout, '(/,1x,a)') 'Dump on file the GW ingredients'

  if(has_auxil_basis) then
    call calculate_eri_3center_eigen(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1, timing=timing_aomo_gw)
  endif

  !
  ! energies.dat file that contains gKS energies
  !
  open(newunit=file_e, file='energies.dat', form='formatted', status='replace')
  write(file_e, *) nvirtual_G-ncore_G-1
  do qstate=ncore_G+1, nvirtual_G-1
    write(file_e, *) energy(qstate, 1)
  enddo
  close(file_e)
  write(stdout, '(1x,a)') 'energies.dat written'

  !
  ! omegas.dat file that contains RPA neutral excitation energies
  !
  open(newunit=file_omega, file='omegas.dat', form='formatted', status='replace')
  write(file_omega, *) wpol%npole_reso
  do spole=1, wpol%npole_reso
    write(file_omega, *) wpol%pole(spole)
  enddo
  close(file_omega)
  write(stdout, '(1x,a)') 'omegas.dat written'

  !
  ! w.bin file that contains the w coefficients:
  !   w_{pq} = \sum_{P ia} ( p q | P ) ( P | i a) [ X_{ia} + Y_{ia} ]
  !
  open(newunit=file_w, file='w.bin', form='unformatted', access='stream', status='replace')
  call clean_allocate('w coeff for dumping', wcoeff, 1, wpol%npole_reso, ncore_G+1, nvirtual_G-1)
  do qspin=1, nspin
    do qstate=ncore_G+1, nvirtual_G-1
      ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
      wcoeff(:, ncore_G+1:nvirtual_G-1) = MATMUL( TRANSPOSE(wpol%residue_left(:, :)) , &
                                                 eri_3center_eigen(:, ncore_G+1:nvirtual_G-1, qstate, qspin) )
      call auxil%sum(wcoeff)
      write(file_w) wcoeff(:, :)
    enddo
  enddo
  call clean_deallocate('w coeff for dumping', wcoeff)
  close(file_w)
  write(stdout, '(1x,a)') 'w.bin written'

  !
  ! v.bin file that contains the Coulomb v on the auxiliary basis
  ! v = ( P | p q)
  !
  open(newunit=file_v, file='v.dat', form='formatted', status='replace')
  write(file_v, *) nauxil_global
  write(file_v, *) nvirtual_G-ncore_G-1
  write(file_v, *) nvirtual_G-ncore_G-1
  close(file_v)
  open(newunit=file_v, file='v.bin', form='unformatted', access='stream', status='replace')
  do qspin=1, nspin
    do qstate=ncore_G+1, nvirtual_G-1
      write(file_v) eri_3center_eigen(:, ncore_G+1:nvirtual_G-1, qstate, qspin)
    enddo
  enddo
  close(file_v)
  write(stdout, '(1x,a)') 'v.bin written'


end subroutine dump_gw_ingredients


end module m_gw_selfenergy_analytic
!=========================================================================
