!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains the calculation of the vertex function
! within different flavors: SOX, SOSEX, 2SOSEX, G3W2, static G3W02
!
!=========================================================================
#include "molgw.h"
module m_g3w2_selfenergy
  use m_definitions
  use m_mpi
  use m_timing
  use m_inputparam
  use m_warning
  use m_basis_set
  use m_spectral_function
  use m_eri_ao_mo
  use m_selfenergy_tools
  use m_tddft_fxc
  use m_linear_response


contains


!=========================================================================
subroutine sosex_selfenergy(basis, occupation, energy, c_matrix, wpol, se)
  implicit none

  type(basis_set)                    :: basis
  real(dp), intent(in)                :: occupation(:, :), energy(:, :)
  real(dp), intent(in)                :: c_matrix(:, :, :)
  type(spectral_function), intent(inout) :: wpol
  type(selfenergy_grid), intent(inout) :: se
  !=====
  integer                 :: nstate
  integer                 :: iomega
  complex(dp), allocatable :: sigma_sosex(:, :, :)
  complex(dp), allocatable :: sigma_sox(:, :, :)
  integer                 :: astate, bstate, cstate
  integer                 :: istate, jstate, kstate, ispin, spole
  integer                 :: pstate
  real(dp), allocatable    :: w_s(:, :)
  real(dp)                :: vcoul, vcoul1, vcoul2
  real(dp)                :: pole_s
  real(dp)                :: fxc
  !=====

  call start_clock(timing_gwgamma_self)
  nstate = SIZE(energy, DIM=1)

  ! Turn dynamic into static
  !call issue_warning('HACK in SOSEX')
  !wpol%pole(:) = wpol%pole(:) * 1000.0d0
  !wpol%residue_left(:,:) = wpol%residue_left(:,:) * SQRT( 1000.0d0 )

  write(stdout, *)
  select case(calc_type%selfenergy_approx)
  case(GWSOX)
    write(stdout, *) 'Perform a one-shot GW+SOX calculation'
  case(GWSOSEX)
    if( calc_type%selfenergy_technique == EVSC ) then
      write(stdout, *) 'Perform an eigenvalue-self-consistent GW+SOSEX calculation'
    else
      write(stdout, *) 'Perform a one-shot GW+SOSEX calculation'
    endif
  case(G3W2, G3W2_NUMERICAL)
    write(stdout, *) 'Perform a one-shot GW+SOSEX calculation to prepare full G3W2'
  case default
    call die('sosex_selfenergy: calculation type unknown')
  end select


  if( gwgamma_tddft_ ) then
    write(stdout, *) 'Include a TDDFT kernel contribution to the vertex'
    write(stdout, '(1x,a,f12.4)') 'Exact-exchange amount: ', alpha_hybrid
    call prepare_tddft(.FALSE., nstate, basis, c_matrix, occupation)
  endif

  if(has_auxil_basis) then
    call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)
  else
    call calculate_eri_4center_mo_uks(c_matrix, ncore_G+1, nvirtual_G-1)
  endif


  call clean_allocate('Temporary array', w_s, ncore_G+1, nvirtual_G-1, ncore_G+1, MAX(nhomo_G, nsemax))


  !
  !
  allocate(sigma_sosex(-se%nomega:se%nomega, nsemin:nsemax, nspin))
  allocate(sigma_sox(-se%nomega:se%nomega, nsemin:nsemax, nspin))

  sigma_sosex(:, :, :)  = 0.0_dp
  sigma_sox(:, :, :)  = 0.0_dp


  write(stdout, *) 'Calculate SOX'

  do ispin=1, nspin

    !==========================
    do bstate=ncore_G+1, nvirtual_G-1
      if( (spin_fact - occupation(bstate, ispin)) / spin_fact < completely_empty) cycle
      if( MODULO( bstate-(ncore_G+1) , poorman%nproc ) /= poorman%rank ) cycle

      do istate=ncore_G+1, nvirtual_G-1
        if( occupation(istate, ispin) / spin_fact < completely_empty ) cycle
        do kstate=ncore_G+1, nvirtual_G-1
          if( occupation(kstate, ispin) / spin_fact < completely_empty ) cycle

          do pstate=nsemin, nsemax

            vcoul1 = evaluate_eri_mo(pstate, istate, ispin, bstate, kstate, ispin)
            vcoul2 = evaluate_eri_mo(istate, bstate, ispin, kstate, pstate, ispin)
            if( gwgamma_tddft_ ) then
              fxc = eval_fxc_rks_singlet(istate, bstate, ispin, kstate, pstate, ispin)
              call grid%sum(fxc)
              vcoul2 = alpha_hybrid * vcoul2 - fxc
            endif
            !
            ! calculate only the diagonal !
            do iomega=-se%nomega, se%nomega
              sigma_sox(iomega, pstate, ispin) = sigma_sox(iomega, pstate, ispin) &
                  - vcoul1 * vcoul2            &
                    / ( se%energy0(pstate, ispin) + se%omega(iomega) &
                        - energy(istate, ispin) - energy(kstate, ispin) + energy(bstate, ispin) - ieta )

            enddo
          enddo

        enddo
      enddo
    enddo

    !==========================
    do cstate=ncore_G+1, nvirtual_G-1
      if( (spin_fact - occupation(cstate, ispin)) / spin_fact < completely_empty) cycle
      if( MODULO( cstate-(ncore_G+1) , poorman%nproc ) /= poorman%rank ) cycle

      do jstate=ncore_G+1, nvirtual_G-1
        if( occupation(jstate, ispin) / spin_fact < completely_empty ) cycle
        do astate=ncore_G+1, nvirtual_G-1
          if( (spin_fact - occupation(astate, ispin)) / spin_fact < completely_empty) cycle

          do pstate=nsemin, nsemax

            vcoul1 = evaluate_eri_mo(pstate, astate, ispin, jstate, cstate, ispin)
            vcoul2 = evaluate_eri_mo(astate, jstate, ispin, cstate, pstate, ispin)
            if( gwgamma_tddft_ ) then
              fxc = eval_fxc_rks_singlet(astate, jstate, ispin, cstate, pstate, ispin)
              call grid%sum(fxc)
              vcoul2 = alpha_hybrid * vcoul2 - fxc
            endif
            !
            ! calculate only the diagonal !
            do iomega=-se%nomega, se%nomega
              sigma_sox(iomega, pstate, ispin) = sigma_sox(iomega, pstate, ispin) &
                  - vcoul1 * vcoul2            &
                    / ( se%energy0(pstate, ispin) + se%omega(iomega) &
                        - energy(astate, ispin) - energy(cstate, ispin) + energy(jstate, ispin) + ieta )
            enddo
          enddo

        enddo
      enddo
    enddo


  enddo

  call poorman%sum(sigma_sox)


  if( calc_type%selfenergy_approx == GWSOSEX .OR. calc_type%selfenergy_approx == G3W2 &
      .OR. calc_type%selfenergy_approx == G3W2_NUMERICAL ) then

    write(stdout, *) 'Calculate dynamical SOSEX'


    do ispin=1, nspin

      do spole=1, wpol%npole_reso

        if( MODULO( spole - 1 , poorman%nproc ) /= poorman%rank ) cycle
        write(stdout, *) 'SOSEX W poles:', spole, ' / ', wpol%npole_reso

        pole_s = wpol%pole(spole)

        if(has_auxil_basis) then
          do pstate=ncore_G+1, MAX(nhomo_G, nsemax)
            ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
            w_s(:, pstate)     = MATMUL( wpol%residue_left(:, spole) , eri_3center_mo(:, :, pstate, ispin) )
          enddo
          call auxil%sum(w_s)
        else
          ! Here just grab the precalculated value
          forall(istate=ncore_G+1:nvirtual_G-1, pstate=ncore_G+1:MAX(nhomo_G, nsemax))
            w_s(istate, pstate) = wpol%residue_left(index_prodstate(istate, pstate) &
                                                    + (ispin-1) * index_prodstate(nvirtual_W-1, nvirtual_W-1), &
                                                   spole)
          end forall
        endif



        !==========================
        do istate=ncore_G+1, nvirtual_G-1
          if( occupation(istate, ispin) / spin_fact < completely_empty ) cycle
          do bstate=ncore_G+1, nvirtual_G-1
            if( (spin_fact - occupation(bstate, ispin)) / spin_fact < completely_empty ) cycle
            do cstate=ncore_G+1, nvirtual_G-1
              if( (spin_fact - occupation(cstate, ispin)) / spin_fact < completely_empty ) cycle

              !
              ! calculate only the diagonal !
              do pstate=nsemin, nsemax

                vcoul = evaluate_eri_mo(istate, cstate, ispin, bstate, pstate, ispin)
                if( gwgamma_tddft_ ) then
                  fxc = eval_fxc_rks_singlet(istate, cstate, ispin, bstate, pstate, ispin)
                  call grid%sum(fxc)
                  vcoul = alpha_hybrid * vcoul - fxc
                endif

                do iomega=-se%nomega, se%nomega
                  sigma_sosex(iomega, pstate, ispin) = sigma_sosex(iomega, pstate, ispin) &
                           - w_s(cstate, pstate) * w_s(bstate, istate) * vcoul                          &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(cstate, ispin) - pole_s + ieta )    &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) &
                                 - energy(cstate, ispin) + energy(istate, ispin) - energy(bstate, ispin) + ieta )


                  sigma_sosex(iomega, pstate, ispin) = sigma_sosex(iomega, pstate, ispin) &
                           + w_s(cstate, pstate) * w_s(bstate, istate) * vcoul                          &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(bstate, ispin) &
                                 - energy(cstate, ispin) + energy(istate, ispin) + ieta )  &
                             / ( energy(bstate, ispin) - energy(istate, ispin) + pole_s - ieta )

                enddo
              enddo

            enddo
          enddo
        enddo

        !==========================
        do jstate=ncore_G+1, nvirtual_G-1
          if( occupation(jstate, ispin) / spin_fact < completely_empty ) cycle
          do astate=ncore_G+1, nvirtual_G-1
            if( (spin_fact - occupation(astate, ispin)) / spin_fact < completely_empty) cycle
            do kstate=ncore_G+1, nvirtual_G-1
              if( occupation(kstate, ispin) / spin_fact < completely_empty ) cycle

              !
              ! calculate only the diagonal !
              do pstate=nsemin, nsemax

                vcoul = evaluate_eri_mo(jstate, kstate, ispin, astate, pstate, ispin)
                if( gwgamma_tddft_ ) then
                  fxc = eval_fxc_rks_singlet(jstate, kstate, ispin, astate, pstate, ispin)
                  call grid%sum(fxc)
                  vcoul = alpha_hybrid * vcoul - fxc
                endif

                do iomega=-se%nomega, se%nomega
                  sigma_sosex(iomega, pstate, ispin) = sigma_sosex(iomega, pstate, ispin) &
                           - w_s(kstate, pstate) * w_s(astate, jstate) * vcoul                          &
                              / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(kstate, ispin) + pole_s - ieta )  &
                              / ( -pole_s + energy(jstate, ispin) - energy(astate, ispin) + ieta )
                enddo
              enddo

            enddo
          enddo
        enddo

        do astate=ncore_G+1, nvirtual_G-1
          if( (spin_fact - occupation(astate, ispin)) / spin_fact < completely_empty  ) cycle
          do jstate=ncore_G+1, nvirtual_G-1
            if( occupation(jstate, ispin) / spin_fact < completely_empty ) cycle
            do kstate=ncore_G+1, nvirtual_G-1
              if( occupation(kstate, ispin) / spin_fact < completely_empty ) cycle

              !
              ! calculate only the diagonal !
              do pstate=nsemin, nsemax

                vcoul = evaluate_eri_mo(astate, kstate, ispin, jstate, pstate, ispin)
                if( gwgamma_tddft_ ) then
                  fxc = eval_fxc_rks_singlet(astate, kstate, ispin, jstate, pstate, ispin)
                  call grid%sum(fxc)
                  vcoul = alpha_hybrid * vcoul - fxc
                endif

                do iomega=-se%nomega, se%nomega
                  sigma_sosex(iomega, pstate, ispin) = sigma_sosex(iomega, pstate, ispin) &
                           - w_s(kstate, pstate) * w_s(astate, jstate) * vcoul                          &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(kstate, ispin) &
                                + energy(astate, ispin) - energy(jstate, ispin)  - ieta )  &
                             / ( energy(jstate, ispin) - energy(astate, ispin) - pole_s + ieta )

                  sigma_sosex(iomega, pstate, ispin) = sigma_sosex(iomega, pstate, ispin) &
                           + w_s(kstate, pstate) * w_s(astate, jstate) * vcoul                          &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(kstate, ispin) &
                                 + energy(astate, ispin) - energy(jstate, ispin)  - ieta )  &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(kstate, ispin) + pole_s - ieta )


                enddo
              enddo

            enddo
          enddo
        enddo

        !==========================
        do astate=ncore_G+1, nvirtual_G-1
          if( (spin_fact - occupation(astate, ispin)) / spin_fact < completely_empty  ) cycle
          do jstate=ncore_G+1, nvirtual_G-1
            if( occupation(jstate, ispin) / spin_fact < completely_empty ) cycle
            do cstate=ncore_G+1, nvirtual_G-1
              if( (spin_fact - occupation(cstate, ispin)) / spin_fact < completely_empty ) cycle

              !
              ! calculate only the diagonal !
              do pstate=nsemin, nsemax

                vcoul = evaluate_eri_mo(astate, cstate, ispin, jstate, pstate, ispin)
                if( gwgamma_tddft_ ) then
                  fxc = eval_fxc_rks_singlet(astate, cstate, ispin, jstate, pstate, ispin)
                  call grid%sum(fxc)
                  vcoul = alpha_hybrid * vcoul - fxc
                endif

                do iomega=-se%nomega, se%nomega
                  sigma_sosex(iomega, pstate, ispin) = sigma_sosex(iomega, pstate, ispin) &
                           + w_s(cstate, pstate) * w_s(astate, jstate) * vcoul                          &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(cstate, ispin) - pole_s + ieta )  &
                             / ( pole_s + energy(astate, ispin) - energy(jstate, ispin) - ieta )

                enddo
              enddo

            enddo
          enddo
        enddo



      enddo !spole
    enddo !ispin

    call poorman%sum(sigma_sosex)

  endif


  write(stdout, '(a)') ' Sigma_c(omega) is calculated'



  forall(pstate=nsemin:nsemax)
    se%sigma(:, pstate, :) = factor_sox * sigma_sox(:, pstate, :) + factor_sosex * sigma_sosex(:, pstate, :)
  end forall


  ! if( print_sigma_) then
  !   call write_selfenergy_omega('selfenergy_sox'    ,energy,exchange_m_vxc_diag,occupation,energy,sigma_sox)
  !   call write_selfenergy_omega('selfenergy_sosex'  ,energy,exchange_m_vxc_diag,occupation,energy,sigma_sosex)
  ! endif


  if( ABS( factor_sosex - 1.0_dp ) < 0.001 ) then
    write(stdout, '(/,a)') ' GW+SOSEX self-energy contributions at E0 (eV)'
  else
    write(stdout, '(/,a)') ' GW+SOSEX2 self-energy contributions at E0 (eV)'
  endif
  write(stdout, '(a)') &
     '   #          E0        SigC_SOX   SigC_G(W(w)-v)GvG SigC_TOT'

  do pstate=nsemin, nsemax
    write(stdout, '(i4,1x,20(1x,f12.6))') pstate, se%energy0(pstate, :)*Ha_eV,          &
                                         sigma_sox(0, pstate, :)%re*Ha_eV,  &
                                         sigma_sosex(0, pstate, :)%re*Ha_eV, &
                                         se%sigma(0, pstate, :)%re*Ha_eV
  enddo




  call clean_deallocate('Temporary array', w_s)

  if(has_auxil_basis) then
    call destroy_eri_3center_mo()
  else
    call destroy_eri_4center_mo_uks()
  endif

  if( gwgamma_tddft_ ) then
    call destroy_tddft()
  endif

  call stop_clock(timing_gwgamma_self)


end subroutine sosex_selfenergy


!=========================================================================
subroutine sosex_selfenergy_analyzed(basis, occupation, energy, c_matrix, wpol, se)
  implicit none

  type(basis_set)                    :: basis
  real(dp), intent(in)               :: occupation(:, :), energy(:, :)
  real(dp), intent(in)               :: c_matrix(:, :, :)
  type(spectral_function), intent(inout) :: wpol
  type(selfenergy_grid), intent(inout) :: se
  !=====
  logical, parameter       :: SCREENED = .FALSE. ! .TRUE. ! .FALSE.
  integer                 :: nstate
  integer                 :: iomega
  complex(dp), allocatable :: sigma_sosex(:, :, :, :)
  complex(dp), allocatable :: sigma_sox(:, :, :, :)
  integer                 :: astate, bstate, cstate
  integer                 :: istate, jstate, kstate, pspin, spole, tpole
  integer                 :: pstate, qstate, ustate
  integer                 :: ibf_auxil
  real(dp), allocatable    :: w_s(:, :)
  real(dp)                :: vcoul, vcoul1, vcoul2
  real(dp)                :: Omega_s, norm
  real(dp)                :: ea, eb, ec, ei, ej, ek
  complex(dp), allocatable :: omega(:)
  ! Store weights for analysis
  ! Huge memory impact
  real(dp), allocatable    :: w_voo_hpl(:, :, :)
  real(dp), allocatable    :: w_ovv_ppl(:, :, :)
  real(dp), allocatable    :: w_vov_pph_ring(:, :, :)
  real(dp), allocatable    :: w_vov_pph_sox(:, :, :)
  real(dp), allocatable    :: w_vov_ppl(:, :, :)
  real(dp), allocatable    :: w_vov_pph(:, :, :)
  real(dp), allocatable    :: w_ovo_hhp_ring(:, :, :)
  real(dp), allocatable    :: w_ovo_hhp_sox(:, :, :)
  real(dp), allocatable    :: w_ovo_hpl(:, :, :)
  real(dp), allocatable    :: w_ovo_hhp(:, :, :)
  real(dp), allocatable    :: w_o_hpl(:, :, :)
  real(dp), allocatable    :: w_v_ppl(:, :, :)
  real(dp), allocatable    :: chi_static(:, :), chi_up(:, :, :), up(:, :, :)
  integer                 :: file
  integer                 :: state_range, nstate2
  character(len=3)        :: ctmp
  !=====

  call start_clock(timing_gwgamma_self)
  nstate = SIZE(energy, DIM=1)

  write(stdout, *)
  select case(calc_type%selfenergy_approx)
  case(GWSOX)
    write(stdout, *) 'Perform a one-shot GW+SOX calculation'
  case(GWSOSEX)
    if( calc_type%selfenergy_technique == EVSC ) then
      write(stdout, *) 'Perform an eigenvalue-self-consistent GW+SOSEX calculation'
    else
      write(stdout, *) 'Perform a one-shot GW+SOSEX calculation'
    endif
  case(G3W2, G3W2_NUMERICAL)
    write(stdout, *) 'Perform a one-shot GW+SOSEX calculation to prepare full G3W2'
  case default
    call die('sosex_selfenergy_analyzed: calculation type unknown')
  end select


  if(has_auxil_basis) then
    call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)
  else
    call die('sosex_selfenergy_analyzed: auxil_basis is compulsory')
  endif

  if( SCREENED ) then
    call issue_warning('hardcoded: use a statically screened Coulomb interaction')
    call clean_allocate('chi static', chi_static, nauxil_global, nauxil_global)
    call wpol%evaluate((0.0_dp, 0.0_dp), chi_static)
    do ibf_auxil=1, nauxil_global
      chi_static(ibf_auxil, ibf_auxil) = chi_static(ibf_auxil, ibf_auxil) + 1.0_dp
    enddo
    allocate(up(nauxil_global, ncore_G+1:nvirtual_G-1, nsemin:nsemax))
    allocate(chi_up(nauxil_global, ncore_G+1:nvirtual_G-1, nsemin:nsemax))
    state_range = nvirtual_G - ncore_G-1
    nstate2 = state_range * ( nsemax - nsemin +1 )
    do pstate=nsemin, nsemax
      do ustate=ncore_G+1, nvirtual_G-1
        up(:, ustate, pstate) = eri_3center_mo(:, ustate, pstate, 1)
      enddo
    enddo
    call DGEMM('T', 'N', nauxil_global, nstate2, nauxil_global, &
                1.0_dp, chi_static, nauxil_global, &
                       up, nauxil_global, &
                0.0_dp, chi_up, nauxil_global)
    deallocate(up)
    call clean_deallocate('chi static', chi_static)
  else
    allocate(chi_up(nauxil_global, ncore_G+1:nvirtual_G-1, nsemin:nsemax))
    do pstate=nsemin, nsemax
      do ustate=ncore_G+1, nvirtual_G-1
        chi_up(:, ustate, pstate) = eri_3center_mo(:, ustate, pstate, 1)
      enddo
    enddo
  endif

  ! Only poorman parallelization is authorized here
  if( poorman%nproc /= world%nproc ) then
    call die('sosex_selfenergy_analyzed: poorman parallelization is compulsory')
  endif


  call clean_allocate('Temporary array', w_s, ncore_G+1, nvirtual_G-1, ncore_G+1, MAX(nhomo_G, nsemax))

  allocate(omega(-se%nomega:se%nomega))
  !
  !
  allocate(sigma_sosex(-se%nomega:se%nomega, nsemin:nsemax, nspin, 6))
  allocate(sigma_sox(-se%nomega:se%nomega, nsemin:nsemax, nspin, 2))
  sigma_sosex(:, :, :, :) = 0.0_dp
  sigma_sox(:, :, :, :)   = 0.0_dp

  !
  ! Allocate weights
  ! hpl stands for hole + plasmon
  ! hhp stands for hole + hole + plasmon
  allocate(w_o_hpl(ncore_G+1:nhomo_G, wpol%npole_reso, nsemin:nsemax))
  allocate(w_voo_hpl(ncore_G+1:nhomo_G, wpol%npole_reso, nsemin:nsemax))
  allocate(w_ovo_hpl(ncore_G+1:nhomo_G, wpol%npole_reso, nsemin:nsemax))
  allocate(w_ovo_hhp_ring(ncore_G+1:nhomo_G, wpol%npole_reso, nsemin:nsemax))
  allocate(w_ovo_hhp_sox(ncore_G+1:nhomo_G, wpol%npole_reso, nsemin:nsemax))
  allocate(w_ovo_hhp(ncore_G+1:nhomo_G, wpol%npole_reso, nsemin:nsemax))
  w_o_hpl(:, :, :)       = 0.0_dp
  w_voo_hpl(:, :, :)     = 0.0_dp
  w_ovo_hpl(:, :, :)     = 0.0_dp
  w_ovo_hhp_ring(:, :, :) = 0.0_dp
  w_ovo_hhp_sox(:, :, :) = 0.0_dp
  w_ovo_hhp(:, :, :)     = 0.0_dp

  ! ppl stands for particle + plasmon
  ! pph stands for particle + particle + hole
  allocate(w_v_ppl(nhomo_G+1:nvirtual_G-1, wpol%npole_reso, nsemin:nsemax))
  allocate(w_ovv_ppl(nhomo_G+1:nvirtual_G-1, wpol%npole_reso, nsemin:nsemax))
  allocate(w_vov_ppl(nhomo_G+1:nvirtual_G-1, wpol%npole_reso, nsemin:nsemax))
  allocate(w_vov_pph_ring(nhomo_G+1:nvirtual_G-1, wpol%npole_reso, nsemin:nsemax))
  allocate(w_vov_pph_sox(nhomo_G+1:nvirtual_G-1, wpol%npole_reso, nsemin:nsemax))
  allocate(w_vov_pph(nhomo_G+1:nvirtual_G-1, wpol%npole_reso, nsemin:nsemax))

  w_v_ppl(:, :, :)       = 0.0_dp
  w_ovv_ppl(:, :, :)     = 0.0_dp
  w_vov_ppl(:, :, :)     = 0.0_dp
  w_vov_pph_ring(:, :, :) = 0.0_dp
  w_vov_pph_sox(:, :, :) = 0.0_dp
  w_vov_pph(:, :, :)     = 0.0_dp


  write(stdout, *) 'Calculate SOX'

  do pspin=1, nspin

    !==========================
    ! SOX ovo

    ! Looping over resonant poles is an handy way to run over (occ, virt) pairs
    do tpole=1, wpol%npole_reso
      if( MODULO( tpole - 1, poorman%nproc ) /= poorman%rank ) cycle
      istate = wpol%transition_table(1, tpole)
      bstate = wpol%transition_table(2, tpole)
      ei = energy(istate, pspin)
      eb = energy(bstate, pspin)

      do kstate=ncore_G+1, nhomo_G
        ek = energy(kstate, pspin)
        do pstate=nsemin, nsemax
          omega(:) = se%energy0(pstate, pspin) + se%omega(:)
          !
          ! calculate only the diagonal!
          qstate = pstate

          !vcoul1 = evaluate_eri_mo(pstate,istate,pspin,bstate,kstate,pspin)
          vcoul1 = DOT_PRODUCT(chi_up(:, istate, pstate), eri_3center_mo(:, bstate, kstate, pspin))
          vcoul2 = evaluate_eri_mo(istate, bstate, pspin, kstate, qstate, pspin)

          sigma_sox(:, pstate, pspin, 1) = sigma_sox(:, pstate, pspin, 1) &
              - vcoul1 * vcoul2 / ( omega(:) - ei + eb - ek - ieta )
          w_ovo_hhp_sox(kstate, tpole, pstate) = w_ovo_hhp_sox(kstate, tpole, pstate) &
                        - vcoul1 * vcoul2

        enddo
      enddo
    enddo

    !==========================
    ! RING ovo

    ! Looping over resonant poles is an handy way to run over (occ, virt) pairs
    do tpole=1, wpol%npole_reso
      if( MODULO( tpole - 1, poorman%nproc ) /= poorman%rank ) cycle
      istate = wpol%transition_table(1, tpole)
      bstate = wpol%transition_table(2, tpole)
      ei = energy(istate, pspin)
      eb = energy(bstate, pspin)

      do kstate=ncore_G+1, nhomo_G
        ek = energy(kstate, pspin)
        do pstate=nsemin, nsemax
          omega(:) = se%energy0(pstate, pspin) + se%omega(:)
          !
          ! calculate only the diagonal!
          qstate = pstate

          vcoul1 = evaluate_eri_mo(pstate, istate, pspin, bstate, kstate, pspin)
          vcoul2 = evaluate_eri_mo(kstate, bstate, pspin, istate, qstate, pspin)

          w_ovo_hhp_ring(kstate, tpole, pstate) = w_ovo_hhp_ring(kstate, tpole, pstate) &
                        + spin_fact * vcoul1 * vcoul2

        enddo
      enddo
    enddo

    !==========================
    ! SOX vov

    ! Looping over resonant poles is an handy way to run over (occ, virt) pairs
    do tpole=1, wpol%npole_reso
      if( MODULO( tpole - 1, poorman%nproc ) /= poorman%rank ) cycle
      jstate = wpol%transition_table(1, tpole)
      astate = wpol%transition_table(2, tpole)
      ea = energy(astate, pspin)
      ej = energy(jstate, pspin)

      do cstate=nhomo_G+1, nvirtual_G-1
        ec = energy(cstate, pspin)
        do pstate=nsemin, nsemax
          omega(:) = se%energy0(pstate, pspin) + se%omega(:)
          !
          ! calculate only the diagonal!
          qstate = pstate

          !vcoul1 = evaluate_eri_mo(pstate,astate,pspin,jstate,cstate,pspin)
          vcoul1 = DOT_PRODUCT(chi_up(:, astate, pstate), eri_3center_mo(:, jstate, cstate, pspin))
          vcoul2 = evaluate_eri_mo(astate, jstate, pspin, cstate, qstate, pspin)

          sigma_sox(:, pstate, pspin, 2) = sigma_sox(:, pstate, pspin, 2) &
              - vcoul1 * vcoul2 / ( omega(:) - ea + ej - ec + ieta )
          w_vov_pph_sox(cstate, tpole, pstate) = w_vov_pph_sox(cstate, tpole, pstate) &
                        - vcoul1 * vcoul2
        enddo
      enddo
    enddo

    !==========================
    ! RING vov

    ! Looping over resonant poles is an handy way to run over (occ, virt) pairs
    do tpole=1, wpol%npole_reso
      if( MODULO( tpole - 1, poorman%nproc ) /= poorman%rank ) cycle
      jstate = wpol%transition_table(1, tpole)
      astate = wpol%transition_table(2, tpole)
      ea = energy(astate, pspin)
      ej = energy(jstate, pspin)

      do cstate=nhomo_G+1, nvirtual_G-1
        ec = energy(cstate, pspin)
        do pstate=nsemin, nsemax
          omega(:) = se%energy0(pstate, pspin) + se%omega(:)
          !
          ! calculate only the diagonal!
          qstate = pstate

          vcoul1 = evaluate_eri_mo(pstate, astate, pspin, jstate, cstate, pspin)
          vcoul2 = evaluate_eri_mo(cstate, jstate, pspin, astate, qstate, pspin)

          w_vov_pph_ring(cstate, tpole, pstate) = w_vov_pph_ring(cstate, tpole, pstate) &
                        + spin_fact * vcoul1 * vcoul2
        enddo
      enddo
    enddo

  enddo ! pspin

  call poorman%sum(sigma_sox)


  write(stdout, *) 'Calculate dynamical SOSEX'


  do pspin=1, nspin

    do spole=1, wpol%npole_reso

      if( MODULO( spole - 1 , poorman%nproc ) /= poorman%rank ) cycle
      write(stdout, *) 'SOSEXanalyzed W poles:', spole, ' / ', wpol%npole_reso

      Omega_s = wpol%pole(spole)

      if( has_auxil_basis ) then
        do pstate=ncore_G+1, MAX(nhomo_G, nsemax)
          ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
          w_s(:, pstate)     = MATMUL( wpol%residue_left(:, spole) , eri_3center_mo(:, :, pstate, pspin) )
        enddo
        call auxil%sum(w_s)
      endif

      !==========================
      ! GW occupied
      do kstate=ncore_G+1, nhomo_G
        do pstate=nsemin, nsemax
          !
          ! calculate only the diagonal !
          qstate = pstate

          w_o_hpl(kstate, spole, pstate) = w_s(kstate, pstate) * w_s(kstate, qstate)
        enddo
      enddo

      !==========================
      ! GW virtual
      do cstate=nhomo_G+1, nvirtual_G-1
        do pstate=nsemin, nsemax
          !
          ! calculate only the diagonal !
          qstate = pstate
          
          w_v_ppl(cstate, spole, pstate) = w_s(cstate, pstate) * w_s(cstate, qstate)
        enddo
      enddo


      !==========================
      ! SOSEX \Sigma^{ovo}
      !do istate=ncore_G+1,nhomo_G
      !  do bstate=nhomo_G+1,nvirtual_G-1
      ! Looping over resonant poles is an handy way to run over (occ, virt) pairs
      do tpole=1, wpol%npole_reso
          istate = wpol%transition_table(1, tpole)
          bstate = wpol%transition_table(2, tpole)
          ei = energy(istate, pspin)
          eb = energy(bstate, pspin)

          do kstate=ncore_G+1, nhomo_G
            ek = energy(kstate, pspin)

            do pstate=nsemin, nsemax
              omega(:) = se%energy0(pstate, pspin) + se%omega(:)
              !
              ! calculate only the diagonal !
              qstate = pstate

              !vcoul = evaluate_eri_mo(istate,pstate,pspin,bstate,kstate,pspin)
              vcoul = DOT_PRODUCT(chi_up(:, istate, pstate), eri_3center_mo(:, bstate, kstate, pspin))

!             sigma_sosex(:,pstate,pspin,1) = sigma_sosex(:,pstate,pspin,1) &
!                       + w_s(kstate,qstate) * w_s(bstate,istate) * vcoul &
!                         / ( se%energy0(pstate,pspin) + se%omega(:) &
!                             - energy(istate,pspin) + energy(bstate,pspin) - energy(kstate,pspin) - ieta ) &
!                      * ( &
!                          -1.0_dp / ( energy(istate,pspin) - energy(bstate,pspin) - Omega_s + ieta )    &
!                          +1.0_dp / ( se%energy0(pstate,pspin) + se%omega(:) - energy(kstate,pspin) + Omega_s - ieta ) &
!                        )

            sigma_sosex(:, pstate, pspin, 1) = sigma_sosex(:, pstate, pspin, 1) &
                       + w_s(kstate, qstate) * w_s(bstate, istate) * vcoul &
                         * ( &
                            - 1.0_dp / ( omega(:) - ek + Omega_s - ieta ) &
                                     / ( ei - eb + Omega_s + ieta ) &
                            )
            sigma_sosex(:, pstate, pspin, 5) = sigma_sosex(:, pstate, pspin, 5) &
                       + w_s(kstate, qstate) * w_s(bstate, istate) * vcoul &
                         * ( &
                              1.0_dp / ( omega(:) - ei + eb - ek - ieta ) &
                                * (-2.0_dp) * Omega_s / ( ei - eb - Omega_s + ieta ) &
                                                      / ( ei - eb + Omega_s + ieta ) &
                            )
              w_ovo_hhp(kstate, tpole, pstate) = w_ovo_hhp(kstate, tpole, pstate) &
                       + w_s(kstate, qstate) * w_s(bstate, istate) * vcoul &
                         * REAL( -2.0_dp * Omega_s / ( ei - eb - Omega_s + ieta ) &
                                                   / ( ei - eb + Omega_s + ieta ) )

              w_ovo_hpl(kstate, spole, pstate) = w_ovo_hpl(kstate, spole, pstate) &
                       + w_s(kstate, qstate) * w_s(bstate, istate) * vcoul &
                         * REAL( -1.0_dp / ( ei - eb + Omega_s + ieta ) )


            enddo
          enddo
        !enddo
      enddo

      !==========================
      ! SOSEX \Sigma^{vov}

      ! Looping over resonant poles is an handy way to run over (occ, virt) pairs
      do tpole=1, wpol%npole_reso
          jstate = wpol%transition_table(1, tpole)
          astate = wpol%transition_table(2, tpole)
          ea = energy(astate, pspin)
          ej = energy(jstate, pspin)

          do cstate=nhomo_G+1, nvirtual_G-1
            ec = energy(cstate, pspin)

            do pstate=nsemin, nsemax
              omega(:) = se%energy0(pstate, pspin) + se%omega(:)
              !
              ! calculate only the diagonal !
              qstate = pstate

              !vcoul = evaluate_eri_mo(astate,pstate,pspin,jstate,cstate,pspin)
              vcoul = DOT_PRODUCT(chi_up(:, astate, pstate), eri_3center_mo(:, jstate, cstate, pspin))

              !sigma_sosex(:,pstate,pspin,2) = sigma_sosex(:,pstate,pspin,2) &
              !         + w_s(cstate,qstate) * w_s(astate,jstate) * vcoul &
              !           / ( se%energy0(pstate,pspin) + se%omega(:) &
              !               - energy(astate,pspin) + energy(jstate,pspin) - energy(cstate,pspin) + ieta ) &
              !        * (  1.0_dp / ( energy(astate,pspin) - energy(jstate,pspin) + Omega_s - ieta )    &
              !            -1.0_dp / ( se%energy0(pstate,pspin) + se%omega(:) - energy(cstate,pspin) - Omega_s + ieta ) )
              sigma_sosex(:, pstate, pspin, 2) = sigma_sosex(:, pstate, pspin, 2) &
                         + w_s(cstate, qstate) * w_s(astate, jstate) * vcoul &
                           * ( &
                                1.0_dp / ( omega(:) - ec - Omega_s + ieta ) &
                                       / ( ea - ej - Omega_s - ieta ) &
                              )
              sigma_sosex(:, pstate, pspin, 6) = sigma_sosex(:, pstate, pspin, 6) &
                         + w_s(cstate, qstate) * w_s(astate, jstate) * vcoul &
                           * ( &
                                1.0_dp / ( omega(:) - ea + ej - ec + ieta ) &
                               * (-2.0_dp) * Omega_s / ( ea - ej + Omega_s + ieta ) &
                                                     / ( ea - ej - Omega_s + ieta ) &
                              )
              w_vov_pph(cstate, tpole, pstate) = w_vov_pph(cstate, tpole, pstate) &
                         + w_s(cstate, qstate) * w_s(astate, jstate) * vcoul &
                 * REAL( -2.0_dp * Omega_s / ( ea - ej + Omega_s + ieta ) &
                                           / ( ea - ej - Omega_s + ieta )  )
              w_vov_ppl(cstate, spole, pstate) = w_vov_ppl(cstate, spole, pstate) &
                         + w_s(cstate, qstate) * w_s(astate, jstate) * vcoul &
                 * REAL( 1.0_dp / ( ea - ej - Omega_s - ieta )  )

            enddo

          enddo
        !enddo
      enddo

      !==========================
      ! SOSEX \Sigma^{ovv}
      do istate=ncore_G+1, nhomo_G
        ei = energy(istate, pspin)
        do bstate=nhomo_G+1, nvirtual_G-1
          eb = energy(bstate, pspin)
          do cstate=nhomo_G+1, nvirtual_G-1
            ec = energy(cstate, pspin)

            do pstate=nsemin, nsemax
              omega(:) = se%energy0(pstate, pspin) + se%omega(:)
              !
              ! calculate only the diagonal !
              qstate = pstate

              !vcoul = evaluate_eri_mo(istate,pstate,pspin,bstate,cstate,pspin)
              vcoul = DOT_PRODUCT(chi_up(:, istate, pstate), eri_3center_mo(:, bstate, cstate, pspin))

              sigma_sosex(:, pstate, pspin, 3) = sigma_sosex(:, pstate, pspin, 3) &
                       + w_s(cstate, qstate) * w_s(bstate, istate) * vcoul &
                          / ( omega(:) - ec - Omega_s + ieta ) &
                          * ( -1.0_dp / ( ei - eb - Omega_s + ieta ) )

              w_ovv_ppl(cstate, spole, pstate) = w_ovv_ppl(cstate, spole, pstate) &
                      + w_s(cstate, qstate) * w_s(bstate, istate) * vcoul &
                        * REAL( -1.0_dp / ( ei - eb - Omega_s + ieta ) )

            enddo

          enddo
        enddo
      enddo
      !==========================
      ! SOSEX \Sigma^{voo}
      do astate=nhomo_G+1, nvirtual_G-1
        ea = energy(astate, pspin)
        do jstate=ncore_G+1, nhomo_G
          ej = energy(jstate, pspin)
          do kstate=ncore_G+1, nhomo_G
            ek = energy(kstate, pspin)

            do pstate=nsemin, nsemax
              omega(:) = se%energy0(pstate, pspin) + se%omega(:)
              !
              ! calculate only the diagonal !
              qstate = pstate

              !vcoul = evaluate_eri_mo(astate,pstate,pspin,jstate,kstate,pspin)
              vcoul = DOT_PRODUCT(chi_up(:, astate, pstate), eri_3center_mo(:, jstate, kstate, pspin))

              sigma_sosex(:, pstate, pspin, 4) = sigma_sosex(:, pstate, pspin, 4) &
                       + w_s(kstate, qstate) * w_s(astate, jstate) * vcoul &
                          / ( omega(:) - ek + Omega_s - ieta ) &
                          / ( ea - ej + Omega_s - ieta )

              w_voo_hpl(kstate, spole, pstate) = w_voo_hpl(kstate, spole, pstate) &
                      + w_s(kstate, qstate) * w_s(astate, jstate) * vcoul &
                        * REAL( 1.0_dp / ( ea - ej + Omega_s - ieta ) )

            enddo

          enddo
        enddo
      enddo


    enddo !spole
  enddo !pspin

  call poorman%sum(sigma_sosex)


  write(stdout, '(a)') ' Sigma_c(omega) is calculated'

  do pstate=nsemin, nsemax
    write(ctmp, '(i3.3)') pstate

    open(newunit=file, file='weights_GWSOSEX_occ_state' // ctmp // '.dat', action='write')
    write(file, '(a,*(a14))') '#state excit', 'GW_hpl   ', 'SOSEXvoo_hpl', 'SOSEXovo_hpl', &
                             'SOSEXovo_hhp', 'RING_hhp', 'SOX_hhp', 'PT2_hhp',&
                             'hpl_pole', 'hhp_pole'
    
    do kstate=ncore_G+1, nhomo_G
      ek = energy(kstate, 1)
      do spole=1, wpol%npole_reso
        istate = wpol%transition_table(1, spole)
        astate = wpol%transition_table(2, spole)
        ei = energy(istate, 1)
        ea = energy(astate, 1)
        norm = ABS(w_o_hpl(kstate, spole, pstate)) &
              + ABS(w_voo_hpl(kstate, spole, pstate)) &
              + ABS(w_ovo_hpl(kstate, spole, pstate)) &
              + ABS(w_ovo_hhp(kstate, spole, pstate)) &
              + ABS(w_ovo_hhp_ring(kstate, spole, pstate)) &
              + ABS(w_ovo_hhp_sox(kstate, spole, pstate))
        if( norm < 1.0e-5_dp ) cycle
        write(file, '(i4,1x,i4,1x,*(f14.6))') kstate, spole, &
                                  w_o_hpl(kstate, spole, pstate), &
                                  w_voo_hpl(kstate, spole, pstate), &
                                  w_ovo_hpl(kstate, spole, pstate), &
                                  w_ovo_hhp(kstate, spole, pstate), &
                                  w_ovo_hhp_ring(kstate, spole, pstate), &
                                  w_ovo_hhp_sox(kstate, spole, pstate), &
                                  w_ovo_hhp_ring(kstate, spole, pstate) &
                                    + w_ovo_hhp_sox(kstate, spole, pstate), &
                                  ek - wpol%pole(spole), &
                                  ek - (ea - ei)
      enddo
    enddo
    close(file)

    open(newunit=file, file='weights_GWSOSEX_virt_state' // ctmp // '.dat', action='write')

    write(file, '(a,*(a14))') '#state excit', 'GW_ppl   ', 'SOSEXovv_ppl', 'SOSEXvov_ppl', &
                             'SOSEXvov_pph', 'RING_pph', 'SOX_pph', 'PT2_pph', &
                             'ppl_pole', 'pph_pole'
    do cstate=nhomo_G+1, nvirtual_G-1
      ec = energy(cstate, 1)
      do spole=1, wpol%npole_reso
        istate = wpol%transition_table(1, spole)
        astate = wpol%transition_table(2, spole)
        ei = energy(istate, 1)
        ea = energy(astate, 1)
        norm = ABS(w_v_ppl(cstate, spole, pstate)) &
              + ABS(w_ovv_ppl(cstate, spole, pstate)) &
              + ABS(w_vov_ppl(cstate, spole, pstate)) &
              + ABS(w_vov_pph(cstate, spole, pstate)) &
              + ABS(w_vov_pph_ring(cstate, spole, pstate)) &
              + ABS(w_vov_pph_sox(cstate, spole, pstate))
        if( norm < 1.0e-5_dp ) cycle
        write(file, '(i4,1x,i4,1x,*(f14.6))') cstate, spole, &
                                  w_v_ppl(cstate, spole, pstate), &
                                  w_ovv_ppl(cstate, spole, pstate), &
                                  w_vov_ppl(cstate, spole, pstate), &
                                  w_vov_pph(cstate, spole, pstate), &
                                  w_vov_pph_ring(cstate, spole, pstate), &
                                  w_vov_pph_sox(cstate, spole, pstate), &
                                  w_vov_pph_ring(cstate, spole, pstate) &
                                    + w_vov_pph_sox(cstate, spole, pstate), &
                                  ec + wpol%pole(spole), &
                                  ec + (ea - ei)
      enddo
    enddo
    close(file)

    open(newunit=file, file='selfenergy_GWSOSEX_parts_state' // ctmp // '.dat', action='write')
    write(file, '(*(a13))') '# omega', 'SOXovo  ', 'SOXvov  ', &
                           'SOSEXovo_pl', 'SOSEXvov_pl', &
                           'SOSEXovv_pl', 'SOSEXvoo_pl', &
                           'SOSEXovo_hp', 'SOSEXvov_hp'
    do iomega=-se%nomega, se%nomega
      write(file, '(*(1x,f12.6))') REAL(se%energy0(pstate, 1) + se%omega(iomega) ) * Ha_eV, &
                                  sigma_sox(iomega, pstate, 1, :)%re * Ha_eV, &
                                  sigma_sosex(iomega, pstate, 1, :)%re * Ha_eV
    enddo
    close(file)
  enddo

  se%sigma(:, :, :) = factor_sox * SUM(sigma_sox(:, :, :, :), DIM=4 ) &
                     + factor_sosex * SUM(sigma_sosex(:, :, :, :), DIM=4)


  ! if( print_sigma_) then
  !   call write_selfenergy_omega('selfenergy_sox'    ,energy,exchange_m_vxc_diag,occupation,energy,sigma_sox)
  !   call write_selfenergy_omega('selfenergy_sosex'  ,energy,exchange_m_vxc_diag,occupation,energy,sigma_sosex)
  ! endif


  if( ABS( factor_sosex - 1.0_dp ) < 0.001 ) then
    write(stdout, '(/,a)') ' GW+SOSEX self-energy contributions at E0 (eV)'
  else
    write(stdout, '(/,a)') ' GW+2SOSEX self-energy contributions at E0 (eV)'
  endif
  write(stdout, '(a)') &
     '   #          E0        SigC_SOX   SigC_G(W(w)-v)GvG SigC_TOT'

  do pstate=nsemin, nsemax
    write(stdout, '(i4,1x,*(1x,f12.6))') pstate, se%energy0(pstate, :)*Ha_eV,          &
                                         sigma_sox(0, pstate, :, :)%re*Ha_eV,  &
                                         sigma_sosex(0, pstate, :, :)%re*Ha_eV, &
                                         se%sigma(0, pstate, :)%re*Ha_eV
  enddo




  call clean_deallocate('Temporary array', w_s)

  if(has_auxil_basis) then
    call destroy_eri_3center_mo()
  else
    call destroy_eri_4center_mo_uks()
  endif

  if( gwgamma_tddft_ ) then
    call destroy_tddft()
  endif

  call stop_clock(timing_gwgamma_self)


end subroutine sosex_selfenergy_analyzed


!=========================================================================
subroutine gwgw0g_selfenergy(occupation, energy, c_matrix, wpol, se)
  implicit none

  real(dp), intent(in)                :: occupation(:, :), energy(:, :)
  real(dp), intent(in)                :: c_matrix(:, :, :)
  type(spectral_function), intent(in) :: wpol
  type(selfenergy_grid), intent(inout) :: se
  !=====
  integer                 :: nstate
  integer                 :: iomega, ibf_auxil
  complex(dp), allocatable :: sigma_gwgw0g(:, :, :)
  complex(dp), allocatable :: sigma_gw0gw0g(:, :, :)
  complex(dp), allocatable :: sigma_gvgw0g(:, :, :)
  integer                 :: astate, bstate, cstate
  integer                 :: istate, jstate, kstate, ispin, spole
  integer                 :: pstate
  real(dp), allocatable    :: w_s(:, :)
  real(dp)                :: v_1, w0_1, w0_2
  real(dp)                :: pole_s
  real(dp), allocatable    :: chi_static(:, :)
  real(dp), allocatable    :: ip(:), bk(:), ib(:), kp(:)
  real(dp), allocatable    :: pa(:), jc(:), aj(:), cp(:)
  real(dp), allocatable    :: ik(:), ac(:), bp(:), ak(:), jp(:), ic(:)
  integer                 :: ustate, state_range, nstate2
  real(dp), allocatable    :: chi_up(:, :, :), up(:, :, :)
  type(spectral_function) :: wpol_static_rpa
  !=====

  call start_clock(timing_gwgamma_self)
  if( .NOT. has_auxil_basis ) call die('gwgw0g_selfenergy: not implemented without an auxiliary basis')

  nstate = SIZE(occupation, DIM=1)

  write(stdout, *)
  select case(calc_type%selfenergy_approx)
  case(GW0GW0G)
    write(stdout, *) 'Perform a one-shot GW+GW0GW0G calculation'
  case(GWGW0G)
      write(stdout, *) 'Perform a one-shot GW+GWGW0G calculation'
  case(GWGW0RPAG)
      write(stdout, *) 'Perform a one-shot G * W + G * W * G * W_0^RPA * G calculation'
  case default
    call die('gwgw0g_selfenergy: calculation type unknown')
  end select


  call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)

  call clean_allocate('Temporary array', w_s, ncore_G+1, nvirtual_G-1, ncore_G+1, MAX(nhomo_G, nsemax))


  call clean_allocate('chi static', chi_static, nauxil_global, nauxil_global)
  !
  ! Calculate a new RPA static screening or use the existing one
  !
  if( calc_type%selfenergy_approx == GWGW0RPAG ) then
    call wpol_static_rpa%init(nstate, occupation, 1, grid_type=STATIC)
    call wpol_static_rpa%vsqrt_chi_vsqrt_rpa(occupation, energy, c_matrix, verbose=.FALSE.)
    chi_static(:, :) = wpol_static_rpa%chi(:, :, 1)
  else
    call wpol%evaluate((0.0_dp, 0.0_dp), chi_static)
  endif

  ! Turn dynamic into static for debug purposes
  !call issue_warning('hack to recover GW+SOX for debug purposes')
  !chi_static(:,:) = 0.0d0   ! to recover GW+SOX
  do ibf_auxil=1, nauxil_global
    chi_static(ibf_auxil, ibf_auxil) = chi_static(ibf_auxil, ibf_auxil) + 1.0_dp
  enddo
  allocate(ib(nauxil_global))
  allocate(kp(nauxil_global))
  allocate(ip(nauxil_global))
  allocate(bk(nauxil_global))
  allocate(aj(nauxil_global))
  allocate(cp(nauxil_global))
  allocate(pa(nauxil_global))
  allocate(jc(nauxil_global))
  allocate(ac(nauxil_global))
  allocate(ik(nauxil_global))
  allocate(bp(nauxil_global))
  allocate(ak(nauxil_global))
  allocate(jp(nauxil_global))
  allocate(ic(nauxil_global))

  allocate(up(nauxil_global, ncore_G+1:nvirtual_G-1, nsemin:nsemax))
  allocate(chi_up(nauxil_global, ncore_G+1:nvirtual_G-1, nsemin:nsemax))
  state_range=nvirtual_G-ncore_G-1
  nstate2 = state_range * ( nsemax - nsemin +1 )
  do pstate=nsemin, nsemax
    do ustate=ncore_G+1, nvirtual_G-1
      up(:, ustate, pstate) = eri_3center_mo(:, ustate, pstate, 1)
    enddo
  enddo
  call DGEMM('T', 'N', nauxil_global, nstate2, nauxil_global, &
              1.0_dp, chi_static, nauxil_global, &
                     up, nauxil_global, &
              0.0_dp, chi_up, nauxil_global)

  !
  !
  allocate(sigma_gwgw0g(-se%nomega:se%nomega, nsemin:nsemax, nspin))
  allocate(sigma_gvgw0g(-se%nomega:se%nomega, nsemin:nsemax, nspin))
  allocate(sigma_gw0gw0g(-se%nomega:se%nomega, nsemin:nsemax, nspin))

  sigma_gwgw0g(:, :, :)  = (0.0_dp, 0.0_dp)
  sigma_gvgw0g(:, :, :)  = (0.0_dp, 0.0_dp)
  sigma_gw0gw0g(:, :, :) = (0.0_dp, 0.0_dp)


  write(stdout, *) 'Calculate two static terms analog to SOX'

  do ispin=1, nspin

    !==========================
    do bstate=nhomo_G+1, nvirtual_G-1
      if( MODULO( bstate-(ncore_G+1) , poorman%nproc ) /= poorman%rank ) cycle

      do istate=ncore_G+1, nhomo_G
        do kstate=ncore_G+1, nhomo_G

          do pstate=nsemin, nsemax

            !v_1 = evaluate_eri_mo(pstate,istate,ispin,bstate,kstate,ispin)
            !v_2 = evaluate_eri_mo(istate,bstate,ispin,kstate,pstate,ispin)

            ip(:) = eri_3center_mo(:, pstate, istate, ispin)
            bk(:) = eri_3center_mo(:, bstate, kstate, ispin)
            v_1  = DOT_PRODUCT( ip(:) , bk(:) )
            !w0_1 = DOT_PRODUCT( ip(:) , MATMUL( chi_static(:,:) , bk(:) ) )
            w0_1 = DOT_PRODUCT( bk(:) , chi_up(:, istate, pstate) )

            ib(:) = eri_3center_mo(:, istate, bstate, ispin)
            kp(:) = eri_3center_mo(:, kstate, pstate, ispin)
            !w0_2 = DOT_PRODUCT( ib(:) , MATMUL( chi_static(:,:) , kp(:) ) )
            w0_2 = DOT_PRODUCT( ib(:) , chi_up(:, kstate, pstate) )

            !
            ! calculate only the diagonal !
            do iomega=-se%nomega, se%nomega
              sigma_gvgw0g(iomega, pstate, ispin) = sigma_gvgw0g(iomega, pstate, ispin) &
                  - v_1 * w0_2            &
                    / ( se%energy0(pstate, ispin) + se%omega(iomega) &
                        - energy(istate, ispin) - energy(kstate, ispin) + energy(bstate, ispin) - ieta )
              sigma_gw0gw0g(iomega, pstate, ispin) = sigma_gw0gw0g(iomega, pstate, ispin) &
                  - w0_1 * w0_2            &
                    / ( se%energy0(pstate, ispin) + se%omega(iomega) &
                        - energy(istate, ispin) - energy(kstate, ispin) + energy(bstate, ispin) - ieta )
            enddo
          enddo

        enddo
      enddo
    enddo

    !==========================
    do cstate=nhomo_G+1, nvirtual_G-1
      if( MODULO( cstate-(ncore_G+1) , poorman%nproc ) /= poorman%rank ) cycle

      do jstate=ncore_G+1, nhomo_G
        do astate=nhomo_G+1, nvirtual_G-1

          do pstate=nsemin, nsemax

            !v_1 = evaluate_eri_mo(pstate,astate,ispin,jstate,cstate,ispin)
            !v_2 = evaluate_eri_mo(astate,jstate,ispin,cstate,pstate,ispin)

            pa(:) = eri_3center_mo(:, pstate, astate, ispin)
            jc(:) = eri_3center_mo(:, jstate, cstate, ispin)
            v_1  = DOT_PRODUCT( pa(:) , jc(:) )
            !w0_1 = DOT_PRODUCT( pa(:) , MATMUL( chi_static(:,:) , jc(:) ) )
            w0_1 = DOT_PRODUCT( jc(:) , chi_up(:, astate, pstate) )

            aj(:) = eri_3center_mo(:, astate, jstate, ispin)
            cp(:) = eri_3center_mo(:, cstate, pstate, ispin)
            !w0_2 = DOT_PRODUCT( aj(:) , MATMUL( chi_static(:,:) , cp(:) ) )
            w0_2 = DOT_PRODUCT( aj(:) , chi_up(:, cstate, pstate) )

            !
            ! calculate only the diagonal !
            do iomega=-se%nomega, se%nomega
              sigma_gvgw0g(iomega, pstate, ispin) = sigma_gvgw0g(iomega, pstate, ispin) &
                  - v_1 * w0_2            &
                    / ( se%energy0(pstate, ispin) + se%omega(iomega) &
                        - energy(astate, ispin) - energy(cstate, ispin) + energy(jstate, ispin) + ieta )
              sigma_gw0gw0g(iomega, pstate, ispin) = sigma_gw0gw0g(iomega, pstate, ispin) &
                  - w0_1 * w0_2            &
                    / ( se%energy0(pstate, ispin) + se%omega(iomega) &
                        - energy(astate, ispin) - energy(cstate, ispin) + energy(jstate, ispin) + ieta )
            enddo
          enddo

        enddo
      enddo
    enddo


  enddo


  call poorman%sum(sigma_gvgw0g)


  if( calc_type%selfenergy_approx == GWGW0G &
     .OR. calc_type%selfenergy_approx == GWGW0RPAG ) then

    write(stdout, *) 'Calculate dynamical term analog to SOSEX'


    do ispin=1, nspin

      do spole=1, wpol%npole_reso

        if( MODULO( spole - 1 , poorman%nproc ) /= poorman%rank ) cycle
        write(stdout, *) 'GWGW0G W poles:', spole, ' / ', wpol%npole_reso

        pole_s = wpol%pole(spole)

        if(has_auxil_basis) then
          do pstate=ncore_G+1, MAX(nhomo_G, nsemax)
            ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
            w_s(:, pstate)     = MATMUL( wpol%residue_left(:, spole) , eri_3center_mo(:, :, pstate, ispin) )
          enddo
          call auxil%sum(w_s)
        else
          ! Here just grab the precalculated value
          forall(istate=ncore_G+1:nvirtual_G-1, pstate=ncore_G+1:MAX(nhomo_G, nsemax))
            w_s(istate, pstate) = wpol%residue_left(index_prodstate(istate, pstate) &
                                     + (ispin-1) * index_prodstate(nvirtual_W-1, nvirtual_W-1), spole)
          end forall
        endif


        !==========================
        do istate=ncore_G+1, nvirtual_G-1
          if( occupation(istate, ispin) / spin_fact < completely_empty ) cycle
          do bstate=ncore_G+1, nvirtual_G-1
            if( (spin_fact - occupation(bstate, ispin)) / spin_fact < completely_empty) cycle
            do kstate=ncore_G+1, nvirtual_G-1
              if( occupation(kstate, ispin) / spin_fact < completely_empty ) cycle

              !
              ! calculate only the diagonal !
              do pstate=nsemin, nsemax

                !v_2 = evaluate_eri_mo(istate,kstate,ispin,bstate,pstate,ispin)
                ik(:) = eri_3center_mo(:, istate, kstate, ispin)
                !bp(:) = eri_3center_mo(:,bstate,pstate,ispin)
                !w0_2 = DOT_PRODUCT( ik(:) , MATMUL( chi_static(:,:) , bp(:) ) )
                w0_2 = DOT_PRODUCT( ik(:) , chi_up(:, bstate, pstate) )

                do iomega=-se%nomega, se%nomega
                  sigma_gwgw0g(iomega, pstate, ispin) = sigma_gwgw0g(iomega, pstate, ispin) &
                           - w_s(kstate, pstate) * w_s(bstate, istate) * w0_2                          &
                              / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(kstate, ispin) + pole_s - ieta )  &
                              / ( -pole_s + energy(istate, ispin) - energy(bstate, ispin) + ieta )
                enddo
              enddo

            enddo
          enddo
        enddo

        !==========================
        do istate=ncore_G+1, nvirtual_G-1
          if( occupation(istate, ispin) / spin_fact < completely_empty ) cycle
          do bstate=ncore_G+1, nvirtual_G-1
            if( (spin_fact - occupation(bstate, ispin)) / spin_fact < completely_empty ) cycle
            do cstate=ncore_G+1, nvirtual_G-1
              if( (spin_fact - occupation(cstate, ispin)) / spin_fact < completely_empty ) cycle

              !
              ! calculate only the diagonal !
              do pstate=nsemin, nsemax

                !v_2 = evaluate_eri_mo(istate,cstate,ispin,bstate,pstate,ispin)
                ic(:) = eri_3center_mo(:, istate, cstate, ispin)
                !bp(:) = eri_3center_mo(:,bstate,pstate,ispin)
                !w0_2 = DOT_PRODUCT( ic(:) , MATMUL( chi_static(:,:) , bp(:) ) )
                w0_2 = DOT_PRODUCT( ic(:) , chi_up(:, bstate, pstate) )

                do iomega=-se%nomega, se%nomega
                  sigma_gwgw0g(iomega, pstate, ispin) = sigma_gwgw0g(iomega, pstate, ispin) &
                           - w_s(cstate, pstate) * w_s(bstate, istate) * w0_2                          &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(cstate, ispin) - pole_s + ieta )    &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) &
                                 - energy(cstate, ispin) + energy(istate, ispin) - energy(bstate, ispin) + ieta )


                  sigma_gwgw0g(iomega, pstate, ispin) = sigma_gwgw0g(iomega, pstate, ispin) &
                           + w_s(cstate, pstate) * w_s(bstate, istate) * w0_2                          &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(bstate, ispin) &
                                 - energy(cstate, ispin) + energy(istate, ispin) + ieta )  &
                             / ( energy(bstate, ispin) - energy(istate, ispin) + pole_s - ieta )

                enddo
              enddo

            enddo
          enddo
        enddo

        !==========================
        do astate=ncore_G+1, nvirtual_G-1
          if( (spin_fact - occupation(astate, ispin)) / spin_fact < completely_empty  ) cycle
          do jstate=ncore_G+1, nvirtual_G-1
            if( occupation(jstate, ispin) / spin_fact < completely_empty ) cycle
            do kstate=ncore_G+1, nvirtual_G-1
              if( occupation(kstate, ispin) / spin_fact < completely_empty ) cycle

              !
              ! calculate only the diagonal !
              do pstate=nsemin, nsemax

                !v_2 = evaluate_eri_mo(astate,kstate,ispin,jstate,pstate,ispin)
                ak(:) = eri_3center_mo(:, astate, kstate, ispin)
                !jp(:) = eri_3center_mo(:,jstate,pstate,ispin)
                !w0_2 = DOT_PRODUCT( ak(:) , MATMUL( chi_static(:,:) , jp(:) ) )
                w0_2 = DOT_PRODUCT( ak(:) , chi_up(:, jstate, pstate) )

                do iomega=-se%nomega, se%nomega
                  sigma_gwgw0g(iomega, pstate, ispin) = sigma_gwgw0g(iomega, pstate, ispin) &
                           - w_s(kstate, pstate) * w_s(astate, jstate) * w0_2                          &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(kstate, ispin) &
                                + energy(astate, ispin) - energy(jstate, ispin)  - ieta )  &
                             / ( energy(jstate, ispin) - energy(astate, ispin) - pole_s + ieta )

                  sigma_gwgw0g(iomega, pstate, ispin) = sigma_gwgw0g(iomega, pstate, ispin) &
                           + w_s(kstate, pstate) * w_s(astate, jstate) * w0_2                          &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(kstate, ispin) &
                                 + energy(astate, ispin) - energy(jstate, ispin)  - ieta )  &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(kstate, ispin) + pole_s - ieta )


                enddo
              enddo

            enddo
          enddo
        enddo

        !==========================
        do astate=ncore_G+1, nvirtual_G-1
          if( (spin_fact - occupation(astate, ispin)) / spin_fact < completely_empty  ) cycle
          do jstate=ncore_G+1, nvirtual_G-1
            if( occupation(jstate, ispin) / spin_fact < completely_empty ) cycle
            do cstate=ncore_G+1, nvirtual_G-1
              if( (spin_fact - occupation(cstate, ispin)) / spin_fact < completely_empty ) cycle

              !
              ! calculate only the diagonal !
              do pstate=nsemin, nsemax

                !v_2 = evaluate_eri_mo(astate,cstate,ispin,jstate,pstate,ispin)
                ac(:) = eri_3center_mo(:, astate, cstate, ispin)
                !jp(:) = eri_3center_mo(:,jstate,pstate,ispin)
                !w0_2 = DOT_PRODUCT( ac(:) , MATMUL( chi_static(:,:) , jp(:) ) )
                w0_2 = DOT_PRODUCT( ac(:) , chi_up(:, jstate, pstate) )

                do iomega=-se%nomega, se%nomega
                  sigma_gwgw0g(iomega, pstate, ispin) = sigma_gwgw0g(iomega, pstate, ispin) &
                           + w_s(cstate, pstate) * w_s(astate, jstate) * w0_2             &
                             / ( se%energy0(pstate, ispin) + se%omega(iomega) - energy(cstate, ispin) - pole_s + ieta )  &
                             / ( pole_s + energy(astate, ispin) - energy(jstate, ispin) - ieta )

                enddo
              enddo

            enddo
          enddo
        enddo



      enddo !spole
    enddo !ispin

    call poorman%sum(sigma_gwgw0g)

  endif


  write(stdout, '(a)') ' Sigma_c(omega) is calculated'



  select case(calc_type%selfenergy_approx)
  case(GW0GW0G)
    forall(pstate=nsemin:nsemax)
      se%sigma(:, pstate, :) = sigma_gw0gw0g(:, pstate, :)
    end forall
    write(stdout, '(/,a)') ' G * W(w=0) * G * W(w=0) * G self-energy contributions at E0 (eV)'
    write(stdout, '(a)') &
      '   #       E0           GW0GW0G'
    do pstate=nsemin, nsemax
      write(stdout, '(i4,1x,*(1x,f12.6))') pstate, se%energy0(pstate, :)*Ha_eV, &
                                          se%sigma(0, pstate, :)%re*Ha_eV
    enddo

  case(GWGW0G)
    forall(pstate=nsemin:nsemax)
      se%sigma(:, pstate, :) = 2.0_dp * sigma_gvgw0g(:, pstate, :) &
                           - sigma_gw0gw0g(:, pstate, :) + 2.0_dp * sigma_gwgw0g(:, pstate, :)
    end forall
    write(stdout, '(/,a)') ' G * W(w) * G * W(w=0) * G self-energy contributions at E0 (eV)'
    write(stdout, '(a)') &
      '   #       E0          GvGW0G        GW0GW0G     G(W-v)GW0G    G(W-W0)GW0G    Total'

    do pstate=nsemin, nsemax
      write(stdout, '(i4,1x,*(1x,f12.6))') pstate, se%energy0(pstate, :)*Ha_eV, &
                                          sigma_gvgw0g(0, pstate, :)%re*Ha_eV, &
                                          sigma_gw0gw0g(0, pstate, :)%re*Ha_eV, &
                                          sigma_gwgw0g(0, pstate, :)%re*Ha_eV, &
              (sigma_gwgw0g(0, pstate, 1)%re+sigma_gvgw0g(0, pstate, :)%re-sigma_gw0gw0g(0, pstate, 1)%re)*Ha_eV, &
                                          se%sigma(0, pstate, :)%re*Ha_eV
    enddo
  case(GWGW0RPAG)
    forall(pstate=nsemin:nsemax)
      se%sigma(:, pstate, :) = sigma_gvgw0g(:, pstate, :) + sigma_gwgw0g(:, pstate, :)
    end forall
    write(stdout, '(/,a)') ' G * W(w) * G * W(w=0)^RPA * G self-energy contributions at E0 (eV)'
    write(stdout, '(a)') &
      '   #       E0          GvGW0G          G(W-v)GW0G      Total'
    do pstate=nsemin, nsemax
      write(stdout, '(i4,1x,*(1x,f12.6))') pstate, se%energy0(pstate, :)*Ha_eV, &
                                          sigma_gvgw0g(0, pstate, :)%re*Ha_eV, &
                                          sigma_gwgw0g(0, pstate, :)%re*Ha_eV, &
                                          se%sigma(0, pstate, :)%re*Ha_eV
    enddo
  case default
    call die('gwgw0g_selfenergy: calculation type unknown')
  end select


  call clean_deallocate('Temporary array', w_s)

  call destroy_eri_3center_mo()

  call clean_deallocate('chi static', chi_static)

  call stop_clock(timing_gwgamma_self)


end subroutine gwgw0g_selfenergy


!=========================================================================
subroutine g3w2_selfenergy(occupation, energy, c_matrix, wpol, se)
  implicit none

  real(dp), intent(in)                :: occupation(:, :), energy(:, :)
  real(dp), intent(in)                :: c_matrix(:, :, :)
  type(spectral_function), intent(in) :: wpol
  type(selfenergy_grid), intent(inout) :: se
  !=====
  integer                 :: nstate
  integer                 :: iomega
  complex(dp), allocatable :: sigma_rest(:, :, :)
  complex(dp), allocatable :: sigma_g3w2(:, :, :, :)
  integer                 :: astate, bstate, cstate
  integer                 :: istate, jstate, kstate, pqspin, spole, tpole
  integer                 :: pstate, qstate
  real(dp), allocatable    :: w_t(:, :), w_s(:, :)
  real(dp)                :: Omega_s, Omega_t
  complex(dp)             :: denom1, denom2, denom3, denom4, num3
  real(dp)                :: omega, num1, num2, ei, ej, ek, ea, eb, ec
  !=====

  call start_clock(timing_gwgamma_self)


  nstate = SIZE(c_matrix, DIM=2)

  write(stdout, *)
  write(stdout, *) 'Perform a one-shot G3W2 calculation'

  if( g3w2_skip_vvv_ ) then
    call issue_warning('g3w2_selfenergy: g3w2_skip_vvv has been triggered')    
  endif
  if( g3w2_skip_vv_ ) then
    call issue_warning('g3w2_selfenergy: g3w2_skip_vv has been triggered')    
  endif

  if(has_auxil_basis) then
    call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)
  else
    call die('not implemented')
  endif


  call clean_allocate('Temporary array', w_s, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)
  call clean_allocate('Temporary array', w_t, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)



  !
  !
  allocate(sigma_g3w2(-se%nomega:se%nomega, nsemin:nsemax, nspin, 6))
  allocate(sigma_rest(-se%nomega:se%nomega, nsemin:nsemax, nspin))

  sigma_g3w2(:, :, :, :)  = 0.0_dp


  do pqspin=1, nspin

    do spole=1, wpol%npole_reso
      if( MODULO( spole - 1 , poorman%nproc ) /= poorman%rank ) cycle
      write(stdout, *) 'W poles for G3W2:', spole, ' / ', wpol%npole_reso

      Omega_s = wpol%pole(spole)

      if(has_auxil_basis) then
        do pstate=ncore_G+1, nvirtual_G-1
          ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v) in MO
          w_s(:, pstate)     = MATMUL( wpol%residue_left(:, spole) , eri_3center_mo(:, :, pstate, pqspin) )
        enddo
      endif

      do tpole=1, wpol%npole_reso

        Omega_t = wpol%pole(tpole)
        if(has_auxil_basis) then
          do pstate=ncore_G+1, nvirtual_G-1
            ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v) in MO
            w_t(:, pstate)     = MATMUL( wpol%residue_left(:, tpole) , eri_3center_mo(:, :, pstate, pqspin) )
          enddo
        endif



        do pstate=nsemin, nsemax
          qstate=pstate
          do iomega=-se%nomega, se%nomega
            omega = se%energy0(pstate, pqspin) + se%omega(iomega)%re

            !
            ! 000
            !

            ! occ R occ R occ
            !  i  t  j  s  k
            do istate=ncore_G+1, nhomo_G
              ei = energy(istate, pqspin)
              do jstate=ncore_G+1, nhomo_G
                ej = energy(jstate, pqspin)
                do kstate=ncore_G+1, nhomo_G
                  ek = energy(kstate, pqspin)
                  num1 = w_t(pstate, istate) * w_t(jstate, kstate)
                  num2 = w_s(qstate, kstate) * w_s(istate, jstate)

                  denom1 = omega + Omega_t - ei -2.0_dp*ieta
                  denom2 = omega + Omega_t + Omega_s - ej -3.0_dp*ieta
                  denom3 = omega + Omega_s - ek -2.0_dp*ieta

                  sigma_g3w2(iomega, pstate, pqspin, 1) = sigma_g3w2(iomega, pstate, pqspin, 1) &
                            + num1 * num2 / denom1 / denom2 / denom3
                enddo
              enddo
            enddo

            !
            ! 001
            !

            ! occ R occ R emp + occ R occ AR emp
            ! + emp AR occ R occ + emp AR occ R occ
            !  i  t  j  s  c
            do astate=nhomo_G+1, nvirtual_G-1
              ea = energy(astate, pqspin)
              do jstate=ncore_G+1, nhomo_G
                ej = energy(jstate, pqspin)
                do kstate=ncore_G+1, nhomo_G
                  ek = energy(kstate, pqspin)
                  num1 = w_t(pstate, astate) * w_t(jstate, kstate)
                  num2 = w_s(qstate, kstate) * w_s(astate, jstate)

                  ! emp R occ R occ
                  denom1 = omega - ek + Omega_s - 2.0_dp*ieta
                  denom2 = omega - ej + Omega_s + Omega_t - 3.0_dp*ieta
                  denom3 = Omega_s + ea  - ej - 3.0_dp*ieta

                  sigma_g3w2(iomega, pstate, pqspin, 2) = sigma_g3w2(iomega, pstate, pqspin, 2) &
                            - 2.0_dp * num1 * num2 / denom1 / denom2 / denom3

                  ! emp AR occ AR emp
                  !denom1 = omega - ek + Omega_s - 2.0_dp*ieta
                  denom2 = omega - ea - Omega_t + 2.0_dp*ieta
                  !denom3 = Omega_s + ea  - ej - 3.0_dp*ieta

                  sigma_g3w2(iomega, pstate, pqspin, 2) = sigma_g3w2(iomega, pstate, pqspin, 2) &
                            + 2.0_dp * num1 * num2 / denom1 / denom2 / denom3

                enddo
              enddo
            enddo


            !
            ! 010
            !

            ! occ R emp R occ + occ AR emp R occ + occ AR emp AR occ + occ R emp AR occ
            !  i  t  b  s  k
            do istate=ncore_G+1, nhomo_G
              ei = energy(istate, pqspin)
              do bstate=nhomo_G+1, nvirtual_G-1
                eb = energy(bstate, pqspin)
                do kstate=ncore_G+1, nhomo_G
                  ek = energy(kstate, pqspin)
                  num1 = w_t(pstate, istate) * w_t(bstate, kstate)
                  num2 = w_s(qstate, kstate) * w_s(istate, bstate)

                  ! occ R emp R occ
                  denom1 = omega - ei + Omega_t - 2.0_dp*ieta
                  denom2 = omega - ei - ek + eb - 3.0_dp*ieta
                  denom3 = omega + Omega_s - ek - 2.0_dp*ieta

                  sigma_g3w2(iomega, pstate, pqspin, 3) = sigma_g3w2(iomega, pstate, pqspin, 3) &
                            - num1 * num2 / denom1 / denom2 / denom3

                  ! occ AR emp R occ
                  ! occ R emp AR occ
                  denom1 = Omega_t + eb - ek - 3.0_dp*ieta
                  !denom2 = omega - ei - ek + eb - 3.0_dp*ieta
                  !denom3 = omega + Omega_s - ek - 2.0_dp*ieta

                  sigma_g3w2(iomega, pstate, pqspin, 3) = sigma_g3w2(iomega, pstate, pqspin, 3) &
                            - 2.0_dp * num1 * num2 / denom1 / denom2 / denom3

                  ! occ AR emp AR occ
                  num3   = (2.0_dp * eb - ei - ek + Omega_s + Omega_t - 6.0_dp * ieta )
                  denom1 = eb - ei + Omega_s - 3.0_dp*ieta
                  !denom2 = omega - ei - ek + eb - 3.0_dp*ieta
                  denom3 = eb - ek + Omega_t - 3.0_dp*ieta
                  denom4 = omega - eb - Omega_s - Omega_t + 3.0_dp*ieta

                  sigma_g3w2(iomega, pstate, pqspin, 3) = sigma_g3w2(iomega, pstate, pqspin, 3) &
                            + num1 * num2 * num3 / denom1 / denom2 / denom3 / denom4

                enddo
              enddo
            enddo

            !
            ! 011 + 110
            !

            ! occ AR emp AR emp + occ R  emp AR emp
            ! emp AR emp R  occ + emp AR emp AR occ
            !  i  t  b  s  c
            if( .NOT. g3w2_skip_vv_ ) then
              do istate=ncore_G+1, nhomo_G
                ei = energy(istate, pqspin)
                do bstate=nhomo_G+1, nvirtual_G-1
                  eb = energy(bstate, pqspin)
                  do cstate=nhomo_G+1, nvirtual_G-1
                    ec = energy(cstate, pqspin)
                    num1 = w_t(pstate, istate) * w_t(bstate, cstate)
                    num2 = w_s(qstate, cstate) * w_s(istate, bstate)

                    ! occ AR emp AR emp
                    denom1 = omega - eb - Omega_s - Omega_t + 3.0_dp*ieta
                    denom2 = Omega_s - ei + eb - 3.0_dp*ieta
                    denom3 = omega - Omega_s - ec + 2.0_dp*ieta

                    sigma_g3w2(iomega, pstate, pqspin, 4) = sigma_g3w2(iomega, pstate, pqspin, 4) &
                              + 2.0_dp * num1 * num2 / denom1 / denom2 / denom3

                    ! occ R emp AR emp
                    denom1 = omega - ei + Omega_t - 2.0_dp*ieta
                    !denom2 = Omega_s - ei + eb - 3.0_dp*ieta
                    !denom3 = omega - Omega_s - ec + 2.0_dp*ieta

                    sigma_g3w2(iomega, pstate, pqspin, 4) = sigma_g3w2(iomega, pstate, pqspin, 4) &
                              - 2.0_dp * num1 * num2 / denom1 / denom2 / denom3
                  enddo
                enddo
              enddo
            endif

            !
            ! 101
            !

            ! emp R occ R emp
            ! emp AR occ R emp
            ! emp AR occ AR emp
            ! emp R occ AR emp
            !  a  t  j  s  c
            if( .NOT. g3w2_skip_vv_ ) then
              do astate=nhomo_G+1, nvirtual_G-1
                ea = energy(astate, pqspin)
                do jstate=ncore_G+1, nhomo_G
                  ej = energy(jstate, pqspin)
                  do cstate=nhomo_G+1, nvirtual_G-1
                    ec = energy(cstate, pqspin)
                    num1 = w_t(pstate, astate) * w_t(jstate, cstate)
                    num2 = w_s(qstate, cstate) * w_s(astate, jstate)

                    ! emp R occ R emp
                    num3   = 2.0_dp * ej  - ea - ec - Omega_s - Omega_t + 6.0_dp * ieta 
                    denom1 = Omega_s - ej + ea - 3.0_dp*ieta
                    denom2 = Omega_t - ej + ec - 3.0_dp*ieta
                    denom3 = omega - ea - ec + ej + 3.0_dp*ieta
                    denom4 = omega - ej + Omega_s + Omega_t - 3.0_dp*ieta

                    sigma_g3w2(iomega, pstate, pqspin, 5) = sigma_g3w2(iomega, pstate, pqspin, 5) &
                              + num1 * num2 * num3 / denom1 / denom2 / denom3 / denom4

                    ! emp AR occ R  emp
                    ! emp R  occ AR emp
                    denom1 = omega - ea - Omega_t + 2.0_dp*ieta
                    denom2 = Omega_s + ea - ej - 3.0_dp*ieta
                    !denom3 = omega - ea  - ec + ej + 3.0_dp*ieta

                    sigma_g3w2(iomega, pstate, pqspin, 5) = sigma_g3w2(iomega, pstate, pqspin, 5) &
                              + 2.0_dp * num1 * num2 / denom1 / denom2 / denom3
                    ! emp AR occ AR emp
                    !denom1 = omega - ea - Omega_t + 2.0_dp*ieta
                    denom2 = omega - ec - Omega_s + 2.0_dp*ieta
                    !denom3 = omega - ea - ec + ej + 3.0_dp*ieta

                    sigma_g3w2(iomega, pstate, pqspin, 5) = sigma_g3w2(iomega, pstate, pqspin, 5) &
                              - num1 * num2 / denom1 / denom2 / denom3

                  enddo
                enddo
              enddo
            endif

            !
            ! 111
            !

            ! emp AR emp AR emp
            !  a  t  b  s  c
            if( .NOT. g3w2_skip_vvv_ ) then
              do astate=nhomo_G+1, nvirtual_G-1
                ea = energy(astate, pqspin)
                do bstate=nhomo_G+1, nvirtual_G-1
                  eb = energy(bstate, pqspin)
                  do cstate=nhomo_G+1, nvirtual_G-1
                    ec = energy(cstate, pqspin)
                    num1 = w_t(pstate, astate) * w_t(bstate, cstate)
                    num2 = w_s(qstate, cstate) * w_s(astate, bstate)

                    denom1 = omega - Omega_t - ea + 2.0_dp*ieta
                    denom2 = omega - eb - Omega_s - Omega_t + 3.0_dp*ieta
                    denom3 = omega - Omega_s - ec + 2.0_dp*ieta

                    sigma_g3w2(iomega, pstate, pqspin, 6) = sigma_g3w2(iomega, pstate, pqspin, 6) &
                              + num1 * num2 / denom1 / denom2 / denom3
                  enddo
                enddo
              enddo
            endif

          enddo ! iomega
        enddo ! pstate


      enddo !tpole
    enddo !spole
  enddo !pqspin

  call poorman%sum(sigma_g3w2)


  write(stdout, '(a)') ' Sigma_c(omega) is calculated'


  !
  ! The input sigma contains the 2SOSEX selfenergy
  sigma_rest(:, :, :) = se%sigma(:, :, :)


  forall(pstate=nsemin:nsemax)
    se%sigma(:, pstate, :) = sigma_rest(:, pstate, :) + SUM(sigma_g3w2(:, pstate, :, :), DIM=3)
  end forall


  ! if( print_sigma_) then
  !   call write_selfenergy_omega('selfenergy_g3w2',energy,exchange_m_vxc_diag,occupation,energy,sigma_g3w2)
  ! endif


  write(stdout, '(/,a)') ' G3W2 self-energy contributions at E0 (eV)'
  write(stdout, '(a)') &
     '   #  SigC_G(W(w)-v)G(W(w)-v)G'
  write(stdout, '(a3,7(a13))') ' ', 'ooo', 'oov+voo', 'ovo','ovv+vvo','vov','vvv','tot'

  do pstate=nsemin, nsemax
    write(stdout, '(i4,1x,*(1x,f12.6))') pstate, &
                                        sigma_g3w2(0, pstate, :, :)%re*Ha_eV, &
                                        SUM(sigma_g3w2(0, pstate, :, :)%re, DIM=2)*Ha_eV
  enddo
  write(stdout, '(a)') &
     '   #          E0       SigC_2SOSEX SigC_G(W-v)G(W-v)G SigC_TOT'
  do pstate=nsemin, nsemax
    write(stdout, '(i4,1x,*(1x,f12.6))') pstate, se%energy0(pstate, :)*Ha_eV,          &
                                        sigma_rest(0, pstate, :)%re*Ha_eV,   &
                                        SUM(sigma_g3w2(0, pstate, :, :)%re, DIM=2)*Ha_eV, &
                                        se%sigma(0, pstate, :)%re*Ha_eV
  enddo

  deallocate(sigma_rest)
  deallocate(sigma_g3w2)


  call clean_deallocate('Temporary array', w_s)
  call clean_deallocate('Temporary array', w_t)

  call destroy_eri_3center_mo()


  call stop_clock(timing_gwgamma_self)


end subroutine g3w2_selfenergy


!=========================================================================
subroutine g3w2_selfenergy_real_grid(basis, occupation, energy, c_matrix, se)
  implicit none

  type(basis_set), intent(in)          :: basis
  real(dp), intent(in)                 :: energy(:, :), occupation(:, :)
  real(dp), intent(in)                 :: c_matrix(:, :, :)
  type(selfenergy_grid), intent(inout) :: se
  !=====
  real(dp)             :: domega
  real(dp)             :: mu
  integer              :: nstate
  integer              :: iomega_sigma, iomegap, iomegapp
  integer              :: pstate, qstate, pqspin
  integer              :: ustate, vstate, wstate
  integer              :: first_omega, last_omega
  complex(dp), allocatable :: sigmag3w2(:, :, :)
  complex(dp)          :: num1, num2
  type(spectral_function) :: wpol_analytic
  real(dp) :: vw(nauxil_global), up(nauxil_global), uv(nauxil_global), qw(nauxil_global)
  real(dp) :: erpa, egw, omegap, omegapp
  complex(dp) :: ieta_u, ieta_v, ieta_w, g_u, g_v, g_w
  complex(dp), allocatable :: chi_omegap(:, :), chi_omegapp(:, :)
  !=====
  domega = eta * 0.5_dp

  call start_clock(timing_gwgamma_self)

  if( nspin > 1 ) then
    call die('g3w2_selfenergy_real_grid only for spin restricted')
  endif
  if( .NOT. has_auxil_basis ) then
    call die('g3w2_selfenergy_real_grid requires an auxiliary basis')
  endif
  first_omega = -se%nomega
  last_omega  = se%nomega
  write(stdout, *) first_omega, last_omega


  write(stdout, '(/,1x,a)') 'G3W2 self-energy on a grid of real frequencies with double integration'
  write(stdout, '(/,1x,a)') '========= Sigma evaluated at frequencies (eV): ========='
  do iomega_sigma=first_omega, last_omega
    write(stdout, '(1x,i4,1x,f14.4,1x,f14.4)') iomega_sigma, se%omega(iomega_sigma)*Ha_eV
  enddo
  write(stdout, '(1x,a)') '========================================================'

  nstate = SIZE(energy, DIM=1)

  mu = 0.5_dp * ( MAXVAL(energy(nhomo_G, :)) + MINVAL(energy(nhomo_G+1, :)) )
  write(stdout, '(1x,a,f12.3)') 'Fermi energy mu (eV): ', mu*Ha_eV


  call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1, timing=timing_aomo_gw)

  if( analytic_chi_ ) then
    call wpol_analytic%init(nstate, occupation, 0, grid_type=NO_GRID)
    call polarizability(.TRUE., .TRUE., basis, occupation, energy, c_matrix, erpa, egw, wpol_analytic)
  else
    call die('g3w2_selfenergy_real_grid: analytic_chi is compulsory')
  endif

  allocate(chi_omegap(nauxil_global, nauxil_global))
  allocate(chi_omegapp(nauxil_global, nauxil_global))

  allocate(sigmag3w2(first_omega:last_omega, nsemin:nsemax, nspin))
  sigmag3w2(:, :, :) = 0.0_dp
  pqspin=1

  do iomegap=-nomega_chi_real, nomega_chi_real
    if( MODULO( iomegap - 1 , poorman%nproc) /= poorman%rank ) cycle
    omegap = domega * iomegap
    write(stdout, '(1x,a,i4,es12.4)') 'External omega loop (eV): ', iomegap, omegap*Ha_eV
    call wpol_analytic%evaluate(omegap, chi_omegap)

    do iomegapp=-nomega_chi_real, nomega_chi_real
      omegapp = domega * (iomegapp-0.5_dp)
      call wpol_analytic%evaluate(omegapp, chi_omegapp)

      do pstate=nsemin, nsemax
        qstate=pstate ! diagonal only
        do ustate=ncore_G+1, nvirtual_G-1
        !do ustate=ncore_G+1,nhomo_G
        !do ustate=nhomo_G+1,nvirtual_G-1
          ieta_u = (0.0, 1.0_dp) * SIGN(eta, energy(ustate, pqspin) - mu)
          do vstate=ncore_G+1, nvirtual_G-1
          !do vstate=ncore_G+1,nhomo_G
          !do vstate=nhomo_G+1,nvirtual_G-1
            ieta_v = (0.0, 1.0_dp) * SIGN(eta, energy(vstate, pqspin) - mu)
            do wstate=ncore_G+1, nvirtual_G-1
            !do wstate=ncore_G+1,nhomo_G
            !do wstate=nhomo_G+1,nvirtual_G-1
              ieta_w = (0.0, 1.0_dp) * SIGN(eta, energy(wstate, pqspin) - mu)

              qw(:) = eri_3center_mo(:, qstate, wstate, pqspin)
              uv(:) = eri_3center_mo(:, ustate, vstate, pqspin)
              up(:) = eri_3center_mo(:, ustate, pstate, pqspin)
              vw(:) = eri_3center_mo(:, vstate, wstate, pqspin)
              num1 = DOT_PRODUCT( vw(:) , MATMUL( chi_omegap(:, :), up(:) ) )
              num2 = DOT_PRODUCT( uv(:) , MATMUL( chi_omegapp(:, :) , qw(:) ) )

              do iomega_sigma=first_omega, last_omega
                g_u = 1.0_dp / ( se%energy0(pstate, pqspin) + se%omega(iomega_sigma) + omegap &
                                 -energy(ustate, pqspin) + ieta_u )
                g_v = 1.0_dp / ( se%energy0(pstate, pqspin) + se%omega(iomega_sigma) + omegap &
                                 + omegapp -energy(vstate, pqspin) + ieta_v )
                g_w = 1.0_dp / ( se%energy0(pstate, pqspin) + se%omega(iomega_sigma) + omegapp &
                                 -energy(wstate, pqspin) + ieta_w )

                sigmag3w2(iomega_sigma, pstate, pqspin) = sigmag3w2(iomega_sigma, pstate, pqspin) &
                                                    -1.0_dp / (2.0_dp * pi)**2 * g_u * num1 * g_v * num2 * g_w &
                                                       * domega**2
              enddo

            enddo
          enddo
        enddo
      enddo


    enddo
  enddo
  call poorman%sum(sigmag3w2)
  write(stdout, *) 'Self-energy correction (eV): '
  do pstate=nsemin, nsemax
    do iomega_sigma=first_omega, last_omega
      write(stdout, '(4x,i4,*(4x,f14.6))') pstate, (se%energy0(pstate, 1) + se%omega(iomega_sigma)%re)*Ha_eV, &
                                          sigmag3w2(iomega_sigma, pstate, 1)%re*Ha_eV
    enddo
  enddo

  se%sigma(:, :, :) = se%sigma(:, :, :) + sigmag3w2(:, :, :)
  deallocate(sigmag3w2)


  call destroy_eri_3center_mo()

  call stop_clock(timing_gwgamma_self)



end subroutine g3w2_selfenergy_real_grid


!=========================================================================
subroutine g3w2_selfenergy_imag_grid(basis, occupation, energy, c_matrix, se)
  implicit none

  type(basis_set), intent(in)          :: basis
  real(dp), intent(in)                 :: energy(:, :), occupation(:, :)
  real(dp), intent(in)                 :: c_matrix(:, :, :)
  type(selfenergy_grid), intent(inout) :: se
  !=====
  integer              :: nstate
  integer              :: iomega_sigma, iomegap, iomegapp, iomega_pair
  real(dp)             :: braket1, braket2
  integer              :: pstate, qstate, pqspin
  integer              :: ustate, vstate, wstate
  integer              :: first_omega, last_omega
  complex(dp), allocatable :: sigmag3w2(:, :, :)
  complex(dp)          :: denom1, denom2, denom3, denoms
  type(spectral_function) :: wpol_imag, wpol_analytic
  !type(spectral_function) :: wpol_one
  real(dp) :: vw(nauxil_global), uv(nauxil_global)
  real(dp) :: erpa, egw
  real(dp), allocatable :: chi_omegap_up(:, :), chi_omegapp_wq(:, :)
  !=====


  if( .NOT. has_auxil_basis ) then
    call die('g3w2_selfenergy_imag_grid: it requires an auxiliary basis')
  endif
  first_omega = LBOUND(se%omega_calc(:), DIM=1)
  last_omega  = UBOUND(se%omega_calc(:), DIM=1)

  call start_clock(timing_gwgamma_self)

  write(stdout, '(/,1x,a)') 'G3W2 self-energy on a grid of imaginary frequencies centered on the HOMO-LUMO gap'
  write(stdout, '(/,1x,a)') '========= Sigma evaluated at frequencies (eV): ========='
  do iomega_sigma=first_omega, last_omega
    write(stdout, '(1x,i4,1x,f14.4,1x,f14.4)') iomega_sigma, se%omega_calc(iomega_sigma)*Ha_eV
  enddo
  write(stdout, '(1x,a)') '========================================================'

  nstate = SIZE(energy, DIM=1)


  call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1, timing=timing_aomo_gw)

  call wpol_imag%init(nstate, occupation, nomega_chi_imag, grid_type=IMAGINARY_QUAD)
  if( analytic_chi_ ) then
    call wpol_analytic%init(nstate, occupation, 0, grid_type=NO_GRID)
    call polarizability(.TRUE., .TRUE., basis, occupation, energy, c_matrix, erpa, egw, wpol_analytic)
    call clean_allocate('Chi', wpol_imag%chi, nauxil_global, nauxil_global, wpol_imag%nomega, verbose=.FALSE.)
    call wpol_analytic%evaluate(wpol_imag%omega, wpol_imag%chi)
  else
    call wpol_imag%vsqrt_chi_vsqrt_rpa(occupation, energy, c_matrix, low_rank=.FALSE.)
  endif

  allocate(chi_omegap_up(nauxil_global, ncore_G+1:nvirtual_G-1))
  allocate(chi_omegapp_wq(nauxil_global, ncore_G+1:nvirtual_G-1))


  allocate(sigmag3w2(first_omega:last_omega, nsemin:nsemax, nspin))
  sigmag3w2(:, :, :) = 0.0_dp

  do pqspin=1, nspin
    do pstate=nsemin, nsemax
      qstate = pstate ! only the diagonal

      iomega_pair = 0

      !
      ! First imaginary axis integral: i omega'
      !
      do iomegap=1, wpol_imag%nomega
        write(stdout, '(1x,a,i4,a,i4)') 'Quadrature on omega prime: ', iomegap, ' / ', wpol_imag%nomega

        !
        ! Second imaginary axis integral: i omega''
        !
        do iomegapp=1, wpol_imag%nomega
          iomega_pair = iomega_pair + 1
          if( MODULO( iomega_pair - 1 , poorman%nproc) /= poorman%rank ) cycle

          call DGEMM('N', 'N', nauxil_global, nvirtual_G-ncore_G-1,nauxil_global, &
                         1.0_dp, wpol_imag%chi(:, :, iomegap), nauxil_global, &
                                eri_3center_mo(:, ncore_G+1:nvirtual_G-1, pstate, pqspin), nauxil_global, &
                         0.0_dp, chi_omegap_up, nauxil_global)

          call DGEMM('N', 'N', nauxil_global, nvirtual_G-ncore_G-1,nauxil_global, &
                         1.0_dp, wpol_imag%chi(:, :, iomegapp), nauxil_global, &
                                eri_3center_mo(:, ncore_G+1:nvirtual_G-1, qstate, pqspin), nauxil_global, &
                         0.0_dp, chi_omegapp_wq, nauxil_global)

          do ustate=ncore_G+1, nvirtual_G-1
            do vstate=ncore_G+1, nvirtual_G-1
              do wstate=ncore_G+1, nvirtual_G-1


                ! v * chi( +/- iw') * v
                vw(:) = eri_3center_mo(:, vstate, wstate, pqspin)
                braket1 = DOT_PRODUCT( vw(:), chi_omegap_up(:, ustate) )


                ! v * chi( +/- iw'') * v
                uv(:) = eri_3center_mo(:, ustate, vstate, pqspin)
                braket2 = DOT_PRODUCT( uv(:) , chi_omegapp_wq(:, wstate) )


                do iomega_sigma=first_omega, last_omega
                  ! +w' +w''
                  denom1 = se%omega_calc(iomega_sigma) + wpol_imag%omega(iomegap)  - energy(ustate, pqspin)
                  denom2 = se%omega_calc(iomega_sigma) + wpol_imag%omega(iomegap)  &
                                                       + wpol_imag%omega(iomegapp) - energy(vstate, pqspin)
                  denom3 = se%omega_calc(iomega_sigma) + wpol_imag%omega(iomegapp) - energy(wstate, pqspin)

                  denoms = 1.0_dp / ( denom1 * denom2 * denom3 )

                  ! -w' +w''
                  denom1 = se%omega_calc(iomega_sigma) - wpol_imag%omega(iomegap)  - energy(ustate, pqspin)
                  denom2 = se%omega_calc(iomega_sigma) - wpol_imag%omega(iomegap)  &
                                                       + wpol_imag%omega(iomegapp) - energy(vstate, pqspin)
                  denom3 = se%omega_calc(iomega_sigma) + wpol_imag%omega(iomegapp) - energy(wstate, pqspin)

                  denoms = denoms + 1.0_dp / ( denom1 * denom2 * denom3 )

                  ! +w' -w''
                  denom1 = se%omega_calc(iomega_sigma) + wpol_imag%omega(iomegap)  - energy(ustate, pqspin)
                  denom2 = se%omega_calc(iomega_sigma) + wpol_imag%omega(iomegap)  &
                                                       - wpol_imag%omega(iomegapp) - energy(vstate, pqspin)
                  denom3 = se%omega_calc(iomega_sigma) - wpol_imag%omega(iomegapp) - energy(wstate, pqspin)

                  denoms = denoms + 1.0_dp / ( denom1 * denom2 * denom3 )

                  ! -w' -w''
                  denom1 = se%omega_calc(iomega_sigma) - wpol_imag%omega(iomegap)  - energy(ustate, pqspin)
                  denom2 = se%omega_calc(iomega_sigma) - wpol_imag%omega(iomegap)  &
                                                       - wpol_imag%omega(iomegapp) - energy(vstate, pqspin)
                  denom3 = se%omega_calc(iomega_sigma) - wpol_imag%omega(iomegapp) - energy(wstate, pqspin)

                  denoms = denoms + 1.0_dp / ( denom1 * denom2 * denom3 )



                  sigmag3w2(iomega_sigma, pstate, pqspin) = sigmag3w2(iomega_sigma, pstate, pqspin) &
                                + wpol_imag%weight_quad(iomegap) * wpol_imag%weight_quad(iomegapp)  &
                                   * ( braket1 * braket2 ) * denoms * ( -1.0_dp / (2.0_dp * pi) ) **2

                enddo !iomega_sigma

              enddo
            enddo
          enddo


        enddo !iomegapp
      enddo !iomegap

    enddo !pstate
  enddo !pqspin
  call poorman%sum(sigmag3w2)

  deallocate(chi_omegap_up)
  deallocate(chi_omegapp_wq)

  write(stdout, *) 'G (W-v) G (W-v) G'
  write(stdout, *) ' state index  imag frequency index  rest self-energy  G3W2 self-energy'
  do pstate=nsemin, nsemax
    do iomega_sigma=first_omega, last_omega
      write(stdout, '(2(4x,i4),4x,f16.6,1x,f16.6)')  pstate, iomega_sigma, &
                                                    se%sigma_calc(iomega_sigma, pstate, 1)%re*Ha_eV, &
                                                    sigmag3w2(iomega_sigma, pstate, 1)%re*Ha_eV
    enddo
  enddo
  se%sigma_calc(:, :, :) = se%sigma_calc(:, :, :) + sigmag3w2(:, :, :)

  deallocate(sigmag3w2)
  call wpol_imag%destroy()

  call destroy_eri_3center_mo()

  call stop_clock(timing_gwgamma_self)

end subroutine g3w2_selfenergy_imag_grid


!=========================================================================
subroutine sosex_selfenergy_imag_grid(basis, occupation, energy, c_matrix, se)
  implicit none

  type(basis_set), intent(in)          :: basis
  real(dp), intent(in)                 :: energy(:, :), occupation(:, :)
  real(dp), intent(in)                 :: c_matrix(:, :, :)
  type(selfenergy_grid), intent(inout) :: se
  !=====
  logical, parameter    :: SOSEX_INCLUDING_SOX=.FALSE.
  integer              :: nstate
  integer              :: iomega_sigma, iomegap
  real(dp)             :: braket1, braket2
  integer              :: mstate, rstate, mpspin, iauxil
  integer              :: prange, isignp
  integer              :: first_omega, last_omega
  complex(dp), allocatable :: sigma_sosex(:, :, :)
  complex(dp)          :: denom1, denom2, omega_cmplx
  type(spectral_function) :: wpol_imag, wpol_one, wpol_analytic
  real(dp) :: erpa, egw
  ! DGEMM
  integer :: astate, istate
  real(dp), allocatable :: eri3_i_m(:, :), eri3_i_a(:, :), eri3_r_m(:, :), eri3_r_a(:, :)
  real(dp), allocatable :: tmp(:, :), braket1_ri(:, :), braket2_ri(:, :)
  real(dp), allocatable :: eri3_a_m(:, :), eri3_a_i(:, :), eri3_r_i(:, :), braket1_ra(:, :), braket2_ra(:, :)
  real(dp) :: chi_wp(nauxil_global, nauxil_global), chi_wwp(nauxil_global, nauxil_global)
  !=====


  if( .NOT. has_auxil_basis ) then
    call die('sosex_selfenergy_imag_grid: requires an auxiliary basis')
  endif
  first_omega = LBOUND(se%omega_calc(:), DIM=1)
  last_omega  = UBOUND(se%omega_calc(:), DIM=1)

  call start_clock(timing_gwgamma_self)

  write(stdout, '(/,1x,a)') 'SOSEX self-energy on a grid of imaginary frequencies centered on the HOMO-LUMO gap'
  write(stdout, '(/,1x,a)') '========= Sigma evaluated at frequencies (eV): ========='
  do iomega_sigma=first_omega, last_omega
    write(stdout, '(1x,i4,1x,f14.4,1x,f14.4)') iomega_sigma, se%omega_calc(iomega_sigma)*Ha_eV
  enddo
  write(stdout, '(1x,a)') '========================================================'

  if( SOSEX_INCLUDING_SOX ) then
    write(stdout, '(1x,a)') 'Calculate SOSEX including SOX: G*W*G*v*G'
  else
    write(stdout, '(1x,a)') 'Calculate SOSEX excluding SOX: G*(W-v)*G*v*G (SOX is calculated elsewhere)'
  endif

  nstate = SIZE(energy, DIM=1)


  call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1, timing=timing_aomo_gw)

  !
  ! Initialize wpol_imag any way to obtain the quadrature grid points and weights
  call wpol_imag%init(nstate, occupation, nomega_chi_imag, grid_type=IMAGINARY_QUAD)
  if( analytic_chi_ ) then
    call wpol_analytic%init(nstate, occupation, 0, grid_type=NO_GRID)
    call polarizability(.TRUE., .TRUE., basis, occupation, energy, c_matrix, erpa, egw, wpol_analytic)
  else
    call wpol_imag%vsqrt_chi_vsqrt_rpa(occupation, energy, c_matrix, low_rank=.FALSE.)
  endif


  prange = nvirtual_G - ncore_G - 1


  allocate(sigma_sosex(first_omega:last_omega, nsemin:nsemax, nspin))
  sigma_sosex(:, :, :) = 0.0_dp

  do mpspin=1, nspin
    do mstate=nsemin, nsemax


      !
      ! Imaginary axis integral
      !
      do iomegap=1, wpol_imag%nomega
        if( analytic_chi_ ) then
          call wpol_analytic%evaluate(wpol_imag%omega(iomegap), chi_wp)
        else
          chi_wp(:, :) = wpol_imag%chi(:, :, iomegap)
        endif
        if( SOSEX_INCLUDING_SOX ) then
          do iauxil=1, nauxil_global
            chi_wp(iauxil, iauxil) = chi_wp(iauxil, iauxil) + 1.0_dp
          enddo
        endif
        write(stdout, '(1x,a,i4,a,i4)') 'Quadrature on omega prime: ', iomegap, ' / ', wpol_imag%nomega

        ! TODO eliminate the isignp loop and move iomega_sigma to become the most inner one
        !
        ! positive and negative omega'
        !
        do isignp=1, -1, -2

          !
          ! loop on Sigma(omega)
          !
          do iomega_sigma=first_omega, last_omega

            omega_cmplx = ABS(se%omega_calc(iomega_sigma)%im+isignp*wpol_imag%omega(iomegap)%im)*im
            chi_wwp(:, :) = 0.0_dp
            do iauxil=1, nauxil_global
              chi_wwp(iauxil, iauxil) = chi_wwp(iauxil, iauxil) + 1.0_dp
            enddo

            allocate(eri3_i_m(nauxil_global, ncore_G+1:nhomo_G))
            allocate(eri3_r_m(nauxil_global, ncore_G+1:nvirtual_G-1))
            allocate(eri3_i_a(nauxil_global, ncore_G+1:nhomo_G))
            allocate(eri3_r_a(nauxil_global, ncore_G+1:nvirtual_G-1))
            allocate(tmp(nauxil_global, ncore_G+1:nhomo_G))
            allocate(braket1_ri(ncore_G+1:nvirtual_G-1, ncore_G+1:nhomo_G))
            allocate(braket2_ri(ncore_G+1:nvirtual_G-1, ncore_G+1:nhomo_G))
            eri3_i_m(:, :) = eri_3center_mo(:, ncore_G+1:nhomo_G, mstate, mpspin)
            eri3_r_m(:, :) = eri_3center_mo(:, ncore_G+1:nvirtual_G-1, mstate, mpspin)

            do astate=nhomo_G+1, nvirtual_G-1
              if( MODULO( astate - (nhomo_G+1) , poorman%nproc) /= poorman%rank ) cycle

              eri3_i_a(:, :) = eri_3center_mo(:, ncore_G+1:nhomo_G, astate, mpspin)
              eri3_r_a(:, :) = eri_3center_mo(:, ncore_G+1:nvirtual_G-1, astate, mpspin)

              !
              ! Fix mstate and astate and then use BLAS level 3 for rstate, istate
              !

              !
              ! Chemist's notation:
              ! ( phi_r phi_m | W( +/- iw') | phi_i phi_a )
              call DGEMM('N', 'N', nauxil_global, nhomo_G-ncore_G,nauxil_global, &
                         1.0_dp, chi_wp(:, :), nauxil_global, &
                                eri3_i_a, nauxil_global, &
                         0.0_dp, tmp, nauxil_global)
              call DGEMM('T', 'N', nvirtual_G-ncore_G-1, nhomo_G-ncore_G,nauxil_global, &
                         1.0_dp, eri3_r_m(:, :), nauxil_global, &
                                tmp(:, :), nauxil_global, &
                         0.0_dp, braket1_ri(:, :), nvirtual_G-ncore_G-1)

              !
              ! Chemist's notation:
              ! ( phi_r phi_a | W( iw +/- iw') | phi_i phi_m )
              call DGEMM('N', 'N', nauxil_global, nhomo_G-ncore_G,nauxil_global, &
                         1.0_dp, chi_wwp(:, :), nauxil_global, &
                                eri3_i_m, nauxil_global, &
                         0.0_dp, tmp, nauxil_global)
              call DGEMM('T', 'N', nvirtual_G-ncore_G-1, nhomo_G-ncore_G,nauxil_global, &
                         1.0_dp, eri3_r_a(:, :), nauxil_global, &
                                tmp(:, :), nauxil_global, &
                         0.0_dp, braket2_ri(:, :), nvirtual_G-ncore_G-1)


              do rstate=ncore_G+1, nvirtual_G-1
                do istate=ncore_G+1, nhomo_G

                  braket1 = braket1_ri(rstate, istate)
                  braket2 = braket2_ri(rstate, istate)

                  denom1 = se%omega_calc(iomega_sigma) + isignp * wpol_imag%omega(iomegap) - energy(rstate, mpspin)
                  denom2 = isignp * wpol_imag%omega(iomegap) + energy(istate, mpspin) - energy(astate, mpspin)

                  sigma_sosex(iomega_sigma, mstate, mpspin) = sigma_sosex(iomega_sigma, mstate, mpspin) &
                                + wpol_imag%weight_quad(iomegap) &
                                    * braket1 / denom1 * braket2 / denom2 / (2.0_dp * pi)

                enddo
              enddo

            enddo

            deallocate(eri3_i_m)
            deallocate(eri3_i_a)
            deallocate(eri3_r_a)
            deallocate(tmp)
            deallocate(braket1_ri)
            deallocate(braket2_ri)


            allocate(eri3_a_m(nauxil_global, nhomo_G+1:nvirtual_G-1))
            allocate(eri3_a_i(nauxil_global, nhomo_G+1:nvirtual_G-1))
            allocate(eri3_r_i(nauxil_global, ncore_G+1:nvirtual_G-1))
            allocate(braket1_ra(ncore_G+1:nvirtual_G-1, nhomo_G+1:nvirtual_G-1))
            allocate(braket2_ra(ncore_G+1:nvirtual_G-1, nhomo_G+1:nvirtual_G-1))
            allocate(tmp(nauxil_global, nhomo_G+1:nvirtual_G-1))
            eri3_a_m(:, :) = eri_3center_mo(:, nhomo_G+1:nvirtual_G-1, mstate, mpspin)

            do istate=ncore_G+1, nhomo_G
              if( MODULO( istate - (ncore_G+1) , poorman%nproc) /= poorman%rank ) cycle

              eri3_a_i(:, :) = eri_3center_mo(:, nhomo_G+1:nvirtual_G-1, istate, mpspin)
              eri3_r_i(:, :) = eri_3center_mo(:, ncore_G+1:nvirtual_G-1, istate, mpspin)

              !
              ! Fix mstate and istate and then use BLAS level 3 for rstate, astate
              !

              !
              ! Chemist's notation:
              ! ( phi_r phi_m | W( +/-iw') | phi_i phi_a )
              call DGEMM('N', 'N', nauxil_global, nvirtual_G-nhomo_G-1,nauxil_global, &
                         1.0_dp, chi_wp(:, :), nauxil_global, &
                                eri3_a_i(:, :), nauxil_global, &
                         0.0_dp, tmp(:, :), nauxil_global)
              call DGEMM('T', 'N', nvirtual_G-ncore_G-1, nvirtual_G-nhomo_G-1,nauxil_global, &
                         1.0_dp, eri3_r_m(:, :), nauxil_global, &
                                tmp(:, :), nauxil_global, &
                         0.0_dp, braket1_ra(:, :), nvirtual_G-ncore_G-1)
              !braket1_ra(:,:) = MATMUL( TRANSPOSE(eri3_r_m(:,ncore_G+1:nvirtual_G-1)), tmp(:,nhomo_G+1:nvirtual_G-1) )

              ! v + v * chi(iw +/- iw') * v
              !
              ! Chemist's notation:
              ! ( phi_r phi_i | W( iw +/- iw') | phi_a phi_m )
              call DGEMM('N', 'N', nauxil_global, nvirtual_G-nhomo_G-1,nauxil_global, &
                         1.0_dp, chi_wwp(:, :), nauxil_global, &
                                eri3_a_m, nauxil_global, &
                         0.0_dp, tmp, nauxil_global)
              call DGEMM('T', 'N', nvirtual_G-ncore_G-1, nvirtual_G-nhomo_G-1,nauxil_global, &
                         1.0_dp, eri3_r_i(:, :), nauxil_global, &
                                tmp(:, :), nauxil_global, &
                         0.0_dp, braket2_ra(:, :), nvirtual_G-ncore_G-1)
              !braket2_ra(:,:) = MATMUL( TRANSPOSE(eri3_r_i(:,ncore_G+1:nvirtual_G-1)), tmp(:,nhomo_G+1:nvirtual_G-1) )


              do astate=nhomo_G+1, nvirtual_G-1
                do rstate=ncore_G+1, nvirtual_G-1


                  ! v + v * chi( +/- iw') * v
                  !braket1 = DOT_PRODUCT( mr(:), MATMUL( wpol_imag%chi(:,:,iomegap), ai(:) ) ) + DOT_PRODUCT( mr, ai )
                  braket1 = braket1_ra(rstate, astate)
                  braket2 = braket2_ra(rstate, astate)
                  !braket1 = DOT_PRODUCT( mr, ai )  !SOX

                  ! v + v * chi(iw +/- iw') * v
                  !braket2 = DOT_PRODUCT( ir(:) , MATMUL( wpol_one%chi(:,:,1) , ma(:) ) ) + DOT_PRODUCT( ir, ma)
                  !braket2 = DOT_PRODUCT( ir, ma)  !SOX or SOSEX

                  denom1 = se%omega_calc(iomega_sigma) + isignp * wpol_imag%omega(iomegap) - energy(rstate, mpspin)
                  denom2 = isignp * wpol_imag%omega(iomegap) + energy(astate, mpspin) - energy(istate, mpspin)

                  sigma_sosex(iomega_sigma, mstate, mpspin) = sigma_sosex(iomega_sigma, mstate, mpspin) &
                                - wpol_imag%weight_quad(iomegap) &
                                    * braket1 / denom1 * braket2 / denom2 / (2.0_dp * pi)

                enddo
              enddo

            enddo

            deallocate(eri3_a_m)
            deallocate(eri3_a_i)
            deallocate(eri3_r_i)
            deallocate(braket1_ra)
            deallocate(braket2_ra)
            deallocate(tmp)

            deallocate(eri3_r_m)

            if( .NOT. analytic_chi_ ) call wpol_one%destroy(verbose=.FALSE.)
          enddo !iomega_sigma
        enddo ! isignp
      enddo !iomegap

    enddo
  enddo
  call poorman%sum(sigma_sosex)

  se%sigma_calc(:, :, :) = se%sigma_calc(:, :, :) + factor_sosex * sigma_sosex(:, :, :)

  deallocate(sigma_sosex)
  call wpol_imag%destroy()
  if( analytic_chi_ ) call wpol_analytic%destroy()

  call destroy_eri_3center_mo()

  call stop_clock(timing_gwgamma_self)

end subroutine sosex_selfenergy_imag_grid


!=========================================================================
subroutine sox_selfenergy_imag_grid(occupation, energy, c_matrix, se)
  implicit none

  real(dp), intent(in)                :: occupation(:, :), energy(:, :)
  real(dp), intent(in)                :: c_matrix(:, :, :)
  type(selfenergy_grid), intent(inout) :: se
  !=====
  integer                 :: iomega_sigma
  integer                 :: first_omega, last_omega
  complex(dp), allocatable :: sigma_gw(:, :, :)
  complex(dp), allocatable :: sigma_sox(:, :, :)
  integer                 :: astate, bstate, cstate
  integer                 :: istate, jstate, kstate, ispin
  integer                 :: pstate
  real(dp)                :: vcoul1, vcoul2
  !=====

  call start_clock(timing_gwgamma_self)

  if( .NOT. has_auxil_basis ) then
    call die('sex_selfenergy_imag_grid: requires an auxiliary basis')
  endif
  first_omega = LBOUND(se%omega_calc(:), DIM=1)
  last_omega  = UBOUND(se%omega_calc(:), DIM=1)

  call start_clock(timing_gwgamma_self)

  write(stdout, '(/,1x,a)') 'SOX self-energy on a grid of imaginary frequencies centered on the HOMO-LUMO gap'
  write(stdout, '(/,1x,a)') '========= Sigma evaluated at frequencies (eV): ========='
  do iomega_sigma=first_omega, last_omega
    write(stdout, '(1x,i4,1x,f14.4,1x,f14.4)') iomega_sigma, se%omega_calc(iomega_sigma)*Ha_eV
  enddo
  write(stdout, '(1x,a)') '========================================================'


  call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)



  !
  !
  allocate(sigma_sox(first_omega:last_omega, nsemin:nsemax, nspin))
  allocate(sigma_gw(first_omega:last_omega, nsemin:nsemax, nspin))

  sigma_sox(:, :, :)  = 0.0_dp


  write(stdout, *) 'Calculate SOX'

  do ispin=1, nspin

    !==========================
    do bstate=nhomo_G+1, nvirtual_G-1
      if( MODULO( bstate-(nhomo_G+1) , poorman%nproc ) /= poorman%rank ) cycle

      do istate=ncore_G+1, nhomo_G
        do kstate=ncore_G+1, nhomo_G

          do pstate=nsemin, nsemax

            vcoul1 = evaluate_eri_mo(pstate, istate, ispin, bstate, kstate, ispin)
            vcoul2 = evaluate_eri_mo(istate, bstate, ispin, kstate, pstate, ispin)

            do iomega_sigma=first_omega, last_omega
              sigma_sox(iomega_sigma, pstate, ispin) = sigma_sox(iomega_sigma, pstate, ispin) &
                  - vcoul1 * vcoul2            &
                    / ( se%omega_calc(iomega_sigma) - energy(istate, ispin) - energy(kstate, ispin) + energy(bstate, ispin) )
            enddo
          enddo

        enddo
      enddo
    enddo

    !==========================
    do cstate=nhomo_G+1, nvirtual_G-1
      if( MODULO( cstate-(nhomo_G+1) , poorman%nproc ) /= poorman%rank ) cycle

      do jstate=ncore_G+1, nhomo_G
        do astate=nhomo_G+1, nvirtual_G-1

          do pstate=nsemin, nsemax

            vcoul1 = evaluate_eri_mo(pstate, astate, ispin, jstate, cstate, ispin)
            vcoul2 = evaluate_eri_mo(astate, jstate, ispin, cstate, pstate, ispin)

            do iomega_sigma=first_omega, last_omega
              sigma_sox(iomega_sigma, pstate, ispin) = sigma_sox(iomega_sigma, pstate, ispin) &
                  - vcoul1 * vcoul2            &
                    / ( se%omega_calc(iomega_sigma) - energy(astate, ispin) - energy(cstate, ispin) + energy(jstate, ispin) )
            enddo
          enddo

        enddo
      enddo
    enddo

  enddo

  call poorman%sum(sigma_sox)
  !
  ! The input sigma contains the GW selfenergy
  sigma_gw(:, :, :) = se%sigma_calc(:, :, :)
  se%sigma_calc(:, :, :) = sigma_gw(:, :, :) + sigma_sox(:, :, :)


  deallocate(sigma_gw, sigma_sox)

  call destroy_eri_3center_mo()

  call stop_clock(timing_gwgamma_self)

end subroutine sox_selfenergy_imag_grid


!=========================================================================
subroutine psd_gw2sosex_selfenergy(energy, c_matrix, wpol, se)
  implicit none

  real(dp), intent(in)                 :: energy(:, :)
  real(dp), intent(in)                 :: c_matrix(:, :, :)
  type(spectral_function), intent(in)  :: wpol
  type(selfenergy_grid), intent(inout) :: se
  !=====
  integer               :: qspin, spole, pstate, istate, astate, iastate
  integer               :: kstate, cstate
  integer               :: nstate2
  complex(dp), allocatable :: sigma(:, :, :)
  real(dp), allocatable :: w_s(:, :, :), w_s_tilde(:, :, :)
  real(dp)              :: v_paik, v_piak, v_paic, v_piac, ei, ea, Omega_s
  ! DEBUG flag
  character(len=32), parameter :: selfenergy_switch = 'PSD' ! 'GW+2SOX+2SOSEX' ! 'GW+SOSEX' ! 'PSD'  ! 'GW'
  !=====

  if(.NOT. has_auxil_basis) return
  if( world%nproc /= poorman%nproc ) call die('psd_2sosex: only implemented with poorman MPI parallelization')
  if( nspin > 1 ) call die('psd_2sosex: not implemented with spin')

  call start_clock(timing_gwgamma_self)

  write(stdout, *) 'Perform a one-shot GW+2SOSEX PDF calculation'

  if(has_auxil_basis) then
    call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1, timing=timing_aomo_gw)
  endif

  call clean_allocate('GW Lehman amplitudes w_s', w_s, 1, wpol%npole_reso, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)

  do qspin=1, nspin

    !do qstate=ncore_G+1, nvirtual_G-1
    !  ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
    !  w_s(:, ncore_G+1:nvirtual_G-1, qstate) = MATMUL( TRANSPOSE(wpol%residue_left(:, :)) , &
    !                                             eri_3center_mo(:, ncore_G+1:nvirtual_G-1, qstate, qspin) )
    !  call auxil%sum(w_s)
    !enddo

    ! w_s^{mn} = \sum_P w_s^P * ( P | m n )
    ! Collapse the 2nd and 3rd indices
    nstate2 = ( nvirtual_G - ncore_G - 1 )**2
    call DGEMM('T', 'N', wpol%npole_reso, nstate2, nauxil_local, &
                1.0_dp, wpol%residue_left(1, 1), nauxil_local, &
                        eri_3center_mo(1, ncore_G+1, ncore_G+1, 1), nauxil_local, &
                0.0_dp, w_s(1, ncore_G+1, ncore_G+1), wpol%npole_reso)
    call auxil%sum(w_s)

  enddo
  qspin = 1


  call clean_allocate('GW+2SOSEX_PSD Lehman amplitudes ~w_s', w_s_tilde, 1, wpol%npole_reso, ncore_G+1, nvirtual_G-1, nsemin, nsemax)


  w_s_tilde(:, :, :) = 0.0_dp
  do pstate=nsemin, nsemax
    !$OMP PARALLEL DO PRIVATE(ei, ea, v_paik, v_piak, Omega_s)
    do kstate=ncore_G+1, nhomo_G

      iastate = 0
      do istate=ncore_G+1, nhomo_G
        ei = energy(istate, qspin)
        do astate=nhomo_G+1, nvirtual_G-1
          ea = energy(astate, qspin)
          iastate = iastate + 1
          if( MODULO( iastate - 1 , poorman%nproc ) /= poorman%rank ) cycle

          v_paik = DOT_PRODUCT( eri_3center_mo(:, pstate, astate, qspin), eri_3center_mo(:, istate, kstate, qspin) )
          v_piak = DOT_PRODUCT( eri_3center_mo(:, pstate, istate, qspin), eri_3center_mo(:, astate, kstate, qspin) )

          do spole=1, wpol%npole_reso
            Omega_s = wpol%pole(spole)
            w_s_tilde(spole, kstate, pstate) = w_s_tilde(spole, kstate, pstate) &
                + w_s(spole, astate, istate) &
                      * (  v_paik * REAL( 1.0_dp / ( ea - ei + Omega_s - ieta) ) &
                         + v_piak * REAL( 1.0_dp / ( ea - ei - Omega_s - ieta) ) )
          enddo

        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(ei, ea, v_paic, v_piac, Omega_s)
    do cstate=nhomo_G+1, nvirtual_G-1

      iastate = 0
      do istate=ncore_G+1, nhomo_G
        ei = energy(istate, qspin)
        do astate=nhomo_G+1, nvirtual_G-1
          ea = energy(astate, qspin)
          iastate = iastate + 1
          if( MODULO( iastate - 1 , poorman%nproc ) /= poorman%rank ) cycle

          v_paic = DOT_PRODUCT( eri_3center_mo(:, pstate, astate, qspin), eri_3center_mo(:, istate, cstate, qspin) )
          v_piac = DOT_PRODUCT( eri_3center_mo(:, pstate, istate, qspin), eri_3center_mo(:, astate, cstate, qspin) )

          do spole=1, wpol%npole_reso
            Omega_s = wpol%pole(spole)
            w_s_tilde(spole, cstate, pstate) = w_s_tilde(spole, cstate, pstate) &
                + w_s(spole, astate, istate) &
                      * (  v_paic * REAL( 1.0_dp / ( ea - ei - Omega_s - ieta) ) &
                         + v_piac * REAL( 1.0_dp / ( ea - ei + Omega_s - ieta) ) )
          enddo

        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  enddo

  call poorman%sum(w_s_tilde)



  allocate(sigma, MOLD=se%sigma)
  sigma(:, :, :) = 0.0_dp

  select case(TRIM(selfenergy_switch))
  case('GW')
    !$OMP PARALLEL PRIVATE(Omega_s)
    !$OMP DO REDUCTION(+:sigma)
    do spole=1, wpol%npole_reso
      Omega_s = wpol%pole(spole)
      do pstate=nsemin, nsemax
        do kstate=ncore_G+1, nhomo_G
          sigma(:, pstate, qspin) = sigma(:, pstate, qspin) &
                 + w_s(spole, kstate, pstate)**2  &
                           / ( se%energy0(pstate, qspin) + se%omega(:) - energy(kstate, qspin) + Omega_s - ieta )
        enddo
      enddo
      do pstate=nsemin, nsemax
        do cstate=nhomo_G+1, nvirtual_G-1
          sigma(:, pstate, qspin) = sigma(:, pstate, qspin) &
                 + w_s(spole, cstate, pstate)**2 &
                           / ( se%energy0(pstate, qspin) + se%omega(:) - energy(cstate, qspin) - Omega_s + ieta )
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  case('GW+2SOX+2SOSEX')
    ! GW + 2SOX + 2SOSEX
    !$OMP PARALLEL PRIVATE(Omega_s)
    !$OMP DO REDUCTION(+:sigma)
    do spole=1, wpol%npole_reso
      Omega_s = wpol%pole(spole)
      do pstate=nsemin, nsemax
        do kstate=ncore_G+1, nhomo_G
          sigma(:, pstate, qspin) = sigma(:, pstate, qspin) &
                 + ( w_s(spole, kstate, pstate) + 2.0_dp * w_s_tilde(spole, kstate, pstate) ) &
                        * w_s(spole, kstate, pstate)  &
                           / ( se%energy0(pstate, qspin) + se%omega(:) - energy(kstate, qspin) + Omega_s - ieta )
        enddo
      enddo
      do pstate=nsemin, nsemax
        do cstate=nhomo_G+1, nvirtual_G-1
          sigma(:, pstate, qspin) = sigma(:, pstate, qspin) &
                 + ( w_s(spole, cstate, pstate) + 2.0_dp * w_s_tilde(spole, cstate, pstate) ) &
                        * w_s(spole, cstate, pstate)  &
                           / ( se%energy0(pstate, qspin) + se%omega(:) - energy(cstate, qspin) - Omega_s + ieta )
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  case('GW+SOSEX')
    ! GW + SOSEX
    !$OMP PARALLEL PRIVATE(Omega_s)
    !$OMP DO REDUCTION(+:sigma)
    do spole=1, wpol%npole_reso
      Omega_s = wpol%pole(spole)
      do pstate=nsemin, nsemax
        do kstate=ncore_G+1, nhomo_G
          sigma(:, pstate, qspin) = sigma(:, pstate, qspin) &
                 + ( w_s(spole, kstate, pstate) + w_s_tilde(spole, kstate, pstate) ) &
                        * w_s(spole, kstate, pstate)  &
                           / ( se%energy0(pstate, qspin) + se%omega(:) - energy(kstate, qspin) + Omega_s - ieta )
        enddo
      enddo
      do pstate=nsemin, nsemax
        do cstate=nhomo_G+1, nvirtual_G-1
          sigma(:, pstate, qspin) = sigma(:, pstate, qspin) &
                 + ( w_s(spole, cstate, pstate) + w_s_tilde(spole, cstate, pstate) ) &
                        * w_s(spole, cstate, pstate)  &
                           / ( se%energy0(pstate, qspin) + se%omega(:) - energy(cstate, qspin) - Omega_s + ieta )
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  case default
    !
    ! This is the full PSD case. The other one are DEBUG
    !
    !$OMP PARALLEL PRIVATE(Omega_s)
    !$OMP DO REDUCTION(+:sigma)
    do spole=1, wpol%npole_reso
      Omega_s = wpol%pole(spole)
      do pstate=nsemin, nsemax
        do kstate=ncore_G+1, nhomo_G
          sigma(:, pstate, qspin) = sigma(:, pstate, qspin) &
                 + ( w_s(spole, kstate, pstate) + w_s_tilde(spole, kstate, pstate) )**2  &
                           / ( se%energy0(pstate, qspin) + se%omega(:) - energy(kstate, qspin) + Omega_s - ieta )
        enddo
      enddo
      do pstate=nsemin, nsemax
        do cstate=nhomo_G+1, nvirtual_G-1
          sigma(:, pstate, qspin) = sigma(:, pstate, qspin) &
                 + ( w_s(spole, cstate, pstate) + w_s_tilde(spole, cstate, pstate) )**2 &
                           / ( se%energy0(pstate, qspin) + se%omega(:) - energy(cstate, qspin) - Omega_s + ieta )
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL


  end select

  call clean_deallocate('GW Lehman amplitudes w_s', w_s)
  call clean_deallocate('GW+2SOSEX_PSD Lehman amplitudes ~w_s', w_s_tilde)

  se%sigma(:, :, :) = sigma(:, :, :)
  deallocate(sigma)


  if(has_auxil_basis) then
    call destroy_eri_3center_mo()
  endif

  call stop_clock(timing_gwgamma_self)


end subroutine psd_gw2sosex_selfenergy


!=========================================================================
subroutine psd_gw2sosex_selfenergy_upfolding(occupation, energy, c_matrix, wpol, exchange_m_vxc, p_matrix)
  implicit none

  real(dp), intent(in)                 :: occupation(:, :), energy(:, :)
  real(dp), intent(in)                 :: c_matrix(:, :, :), exchange_m_vxc(:, :, :)
  type(spectral_function), intent(in)  :: wpol
  real(dp), intent(out), optional      :: p_matrix(:, :, :)
  !=====
  character(len=4)     :: ctmp
  integer              :: nstate, nstate2, nbf
  integer              :: pstate
  integer              :: qstate, qspin
  real(dp)             :: sign_i, mu
  real(dp)             :: weight
  real(dp), allocatable :: matrix_wing(:, :), matrix_head(:, :), matrix_diag(:)
  real(dp), allocatable :: super_matrix(:, :), eigval(:)
  integer              :: nmat, mwing, imat, jmat
  integer              :: mstate, astate, iastate, cstate, istate, kstate, spole
  integer              :: irecord
  integer              :: fu
  real(dp)             :: ea, ei, v_paik, v_piak, Omega_s, v_paic, v_piac, ecorr
  real(dp), allocatable :: w_s(:, :, :), w_s_tilde(:, :, :)
  ! DEBUG flag
  character(len=32), parameter :: selfenergy_switch = 'PSD' ! 'PSD' ! 'GW'
  !=====

  call start_clock(timing_gwgamma_self)

  nstate = SIZE(energy, DIM=1)
  nbf = SIZE(c_matrix, DIM=1)

  write(stdout, *)
  select case(TRIM(selfenergy_switch))
  case('GW')
    write(stdout, *) 'Perform a one-shot GW calculation with super matrix formulation'
  case('PSD')
    write(stdout, *) 'Perform a one-shot PSD-ized GW+2SOSEX calculation with super matrix formulation'
  case default
    call die('psd_gw2sosex_selfenergy_upfolding: invalid choice')
  end select

  if( nspin > 1 ) call die('psd_gw2sosex_selfenergy_upfolding: not functional for nspin>1')

  if( nsemin /= ncore_G+1 .OR. nsemax /= nvirtual_G-1 ) then
    write(stdout, '(1x,a,i5,1x,i5)') 'nsemin ?= ncore_G+1    ', nsemin, ncore_G+1
    write(stdout, '(1x,a,i5,1x,i5)') 'nsemax ?= nvirtual_G-1 ', nsemax, nvirtual_G-1
    call die('psd_gw2sosex_selfenergy_upfolding: selfenergy state range should contain all the active states')
  endif

  if(has_auxil_basis) then
    call calculate_eri_3center_mo(c_matrix, nsemin, nsemax, ncore_G+1, nvirtual_G-1, timing=timing_aomo_gw)
  else
    call die('psd_gw2sosex_selfenergy_upfolding: only works with an auxiliary basis')
  endif


  !
  ! Firstly, prepare w_s^{mn}
  !
  call clean_allocate('GW Lehman amplitudes w_s', w_s, 1, wpol%npole_reso, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)

  do qspin=1, nspin
    ! w_s^{mn} = \sum_P w_s^P * ( P | m n)
    ! Collapse the 2nd and 3rd indices
    nstate2 = ( nvirtual_G - ncore_G - 1 )**2
    call DGEMM('T', 'N', wpol%npole_reso, nstate2, nauxil_local, &
                1.0_dp, wpol%residue_left(1, 1), nauxil_local, &
                        eri_3center_mo(1, ncore_G+1, ncore_G+1, 1), nauxil_local, &
                0.0_dp, w_s(1, ncore_G+1, ncore_G+1), wpol%npole_reso)
    call auxil%sum(w_s)
  enddo

  !
  ! Secondly, prepare tilde w_s^{mn}
  !
  call clean_allocate('GW+2SOSEX_PSD Lehman amplitudes ~w_s', w_s_tilde, 1, wpol%npole_reso, ncore_G+1, nvirtual_G-1, nsemin, nsemax)

  w_s_tilde(:, :, :) = 0.0_dp
  do qspin=1, nspin
    do pstate=nsemin, nsemax
      !$OMP PARALLEL DO PRIVATE(ei, ea, v_paik, v_piak, Omega_s)
      do kstate=ncore_G+1, nhomo_G

        iastate = 0
        do istate=ncore_G+1, nhomo_G
          ei = energy(istate, qspin)
          do astate=nhomo_G+1, nvirtual_G-1
            ea = energy(astate, qspin)
            iastate = iastate + 1
            if( MODULO( iastate - 1 , poorman%nproc ) /= poorman%rank ) cycle

            v_paik = DOT_PRODUCT( eri_3center_mo(:, pstate, astate, qspin), eri_3center_mo(:, istate, kstate, qspin) )
            v_piak = DOT_PRODUCT( eri_3center_mo(:, pstate, istate, qspin), eri_3center_mo(:, astate, kstate, qspin) )

            do spole=1, wpol%npole_reso
              Omega_s = wpol%pole(spole)
              w_s_tilde(spole, kstate, pstate) = w_s_tilde(spole, kstate, pstate) &
                  + w_s(spole, astate, istate) &
                        * (  v_paik * REAL( 1.0_dp / ( ea - ei + Omega_s - ieta) ) &
                           + v_piak * REAL( 1.0_dp / ( ea - ei - Omega_s - ieta) ) )
            enddo

          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ei, ea, v_paic, v_piac, Omega_s)
      do cstate=nhomo_G+1, nvirtual_G-1

        iastate = 0
        do istate=ncore_G+1, nhomo_G
          ei = energy(istate, qspin)
          do astate=nhomo_G+1, nvirtual_G-1
            ea = energy(astate, qspin)
            iastate = iastate + 1
            if( MODULO( iastate - 1 , poorman%nproc ) /= poorman%rank ) cycle

            v_paic = DOT_PRODUCT( eri_3center_mo(:, pstate, astate, qspin), eri_3center_mo(:, istate, cstate, qspin) )
            v_piac = DOT_PRODUCT( eri_3center_mo(:, pstate, istate, qspin), eri_3center_mo(:, astate, cstate, qspin) )

            do spole=1, wpol%npole_reso
              Omega_s = wpol%pole(spole)
              w_s_tilde(spole, cstate, pstate) = w_s_tilde(spole, cstate, pstate) &
                  + w_s(spole, astate, istate) &
                        * (  v_paic * REAL( 1.0_dp / ( ea - ei - Omega_s - ieta) ) &
                           + v_piac * REAL( 1.0_dp / ( ea - ei + Omega_s - ieta) ) )
            enddo

          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

    enddo
  enddo

  call poorman%sum(w_s_tilde)


  mstate = nvirtual_G - ncore_G - 1
  nmat   = mstate * ( 1 + wpol%npole_reso)
  mwing  = mstate * wpol%npole_reso

  !    | M_head | M_wing**T |
  !    | M_wing | M_diag    |

  call clean_allocate('Matrix head', matrix_head, mstate, mstate)
  call clean_allocate('Matrix wing', matrix_wing, mwing, mstate)


  allocate(matrix_diag(mwing))
  matrix_head(:, :) = 0.0_dp
  matrix_wing(:, :) = 0.0_dp
  matrix_diag(:)    = 0.0_dp

  matrix_head(:, :) = exchange_m_vxc(ncore_G+1:nvirtual_G-1, ncore_G+1:nvirtual_G-1, 1)  ! spin index set to 1

  write(stdout, '(/,1x,a,i8,a,i8)') 'Diagonalization problem of size: ', nmat, ' x ', nmat

  do qspin=1, nspin
    do qstate=ncore_G+1, nvirtual_G-1 ! INNER LOOP of G

      if( MODULO( qstate - (ncore_G+1) , poorman%nproc) /= poorman%rank ) cycle
      !
      ! indeces
      pstate = qstate - ncore_G
      sign_i = MERGE(-1.0_dp, 1.0_dp, occupation(qstate, qspin) / spin_fact > completely_empty )
      irecord = ( pstate - 1 ) * wpol%npole_reso

      !
      ! Head
      matrix_head(pstate, pstate) = matrix_head(pstate, pstate) + energy(qstate, qspin)

      !
      ! Diagonal
      matrix_diag(irecord+1:irecord+wpol%npole_reso) = energy(qstate, qspin) + sign_i * wpol%pole(:)

      !
      ! Wing
      ! Here transform (sqrt(v) * chi * sqrt(v)) into  (v * chi * v)
      select case(TRIM(selfenergy_switch))
      case('GW')
        matrix_wing(irecord+1:irecord+wpol%npole_reso, :) = w_s(1:wpol%npole_reso, qstate, ncore_G+1:nvirtual_G-1)
      case('PSD')
        matrix_wing(irecord+1:irecord+wpol%npole_reso, :) = w_s(1:wpol%npole_reso, qstate, ncore_G+1:nvirtual_G-1) &
                                                           + w_s_tilde(1:wpol%npole_reso, qstate, ncore_G+1:nvirtual_G-1)
      end select


    enddo ! qstate
  enddo ! qspin

  ecorr = 0.0_dp
  do astate=nhomo_G+1, nvirtual_G-1
    do istate=ncore_G+1, nhomo_G
      ecorr = ecorr - spin_fact * SUM( w_s(:, istate, astate)**2 / ( energy(astate, 1) - energy(istate, 1) + wpol%pole(:) ) )
    enddo
  enddo
  write(stdout,*) 'FBFB E correlation GW: ', ecorr

  ecorr = 0.0_dp
  do astate=nhomo_G+1, nvirtual_G-1
    do istate=ncore_G+1, nhomo_G
      ecorr = ecorr - spin_fact * SUM( ( w_s(:, istate, astate) + w_s_tilde(:, istate, astate) )**2 &
                     / ( energy(astate, 1) - energy(istate, 1) + wpol%pole(:) ) )
    enddo
  enddo
  write(stdout,*) 'FBFB E correlation PSD: ', ecorr


  write(stdout, '(a)') ' Matrix is setup'

  if( .FALSE. ) then
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
        write(fu, '(1x,i7,1x,i7,1x,e16.8)') imat + mstate, imat+mstate, matrix_diag(imat)*Ha_eV
      enddo
    endif

    do jmat=1, mstate
      do imat=1, mwing
        write(fu, '(1x,i7,1x,i7,1x,e16.8)') mstate + imat, jmat, matrix_wing(imat, jmat)*Ha_eV
      enddo
    enddo
    close(fu)
  endif

  !
  ! If the super_matrix is small enough, then diagonalize it!
  !
  if( nmat < 30001 ) then

    mu = ( MINVAL(energy(nhomo_G+1, :)) + MAXVAL(energy(nhomo_G, :)) ) / 2.0_dp
    write(stdout, '(1x,a,f12.6)') 'Center of the HOMO-LUMO gap (eV): ', mu * Ha_eV
    allocate(eigval(nmat))

    write(stdout, *) 'Diagonalize the large sparse super_matrix as if it were dense'
    call clean_allocate('Super matrix', super_matrix, nmat, nmat)
    super_matrix(:, :)                    = 0.0_dp
    super_matrix(1:mstate, 1:mstate)      = matrix_head(:, :)
    super_matrix(1:mstate, mstate+1:nmat) = TRANSPOSE(matrix_wing(:, :))
    super_matrix(mstate+1:nmat, 1:mstate) = matrix_wing(:, :)
    do imat=mstate+1, nmat
      super_matrix(imat, imat) = matrix_diag(imat-mstate)
    enddo
    call diagonalize(postscf_diago_flavor, super_matrix, eigval)

    write(stdout, '(1x,a,i8)') 'Number of non-negligible poles: ', COUNT( SUM(super_matrix(1:mstate, :)**2, DIM=1) > 1.0e-3_dp )

    write(stdout, '(/,a)') '============== Poles in eV , weight ==============='
    open(newunit=fu, file='GREENS_FUNCTION', action='write')
    do jmat=1, nmat
      weight = SUM(super_matrix(1:mstate, jmat)**2)
      if( weight > 1.0e-2_dp ) then
        pstate = MAXLOC(ABS(super_matrix(1:mstate, jmat)), DIM=1)
        write(stdout, '(1x,a,i5.5,a,f16.6,4x,f12.6)') 'Projection on state ', pstate, ': ', eigval(jmat)*Ha_eV, weight
      endif
      write(fu, '(1x,f16.6,4x,f12.6)') eigval(jmat)*Ha_eV, SUM(super_matrix(1:mstate, jmat)**2)
    enddo
    close(fu)
    ! If eigenvalue lower than the middle of the HOMO-LUMO gap,
    ! then consider the excitation is occupied
    write(stdout, '(1x,a,f12.6)') 'Number of electrons: ', & 
             spin_fact * SUM( SUM(super_matrix(1:mstate, :)**2, DIM=1), MASK=(eigval(:) < mu) )
    write(stdout, *) '==================================================='


    if( PRESENT(p_matrix) ) then
      write(stdout, '(/,1x,a)') 'Evaluate the density-matrix P'
      call greensfunction_supermatrix_to_density_matrix(occupation, energy, c_matrix, super_matrix, eigval, p_matrix)
    endif

    call clean_deallocate('Super matrix', super_matrix)
    deallocate(eigval)

  endif

  call clean_deallocate('GW+2SOSEX_PSD Lehman amplitudes ~w_s', w_s_tilde)
  call clean_deallocate('GW Lehman amplitudes w_s', w_s)

  call clean_deallocate('Matrix wing', matrix_wing)
  call clean_deallocate('Matrix head', matrix_head)
  deallocate(matrix_diag)
  if(has_auxil_basis) then
    call destroy_eri_3center_mo()
  endif

  call stop_clock(timing_gwgamma_self)


end subroutine psd_gw2sosex_selfenergy_upfolding


end module m_g3w2_selfenergy
!=========================================================================
