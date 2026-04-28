!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the many-body perturbation theory to obtain the (perturbative) density matrix
!
!=========================================================================
#include "molgw.h"
module m_mbpt_density_matrix

  use m_definitions
  use m_mpi
  use m_warning
  use m_timing
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam
  use m_hamiltonian_onebody
  use m_selfenergy_tools
  use m_hamiltonian_tools
  use m_spectral_function
  use m_io

  implicit none



contains


!=========================================================================
subroutine pt2_density_matrix(occupation, energy, c_matrix, p_matrix_pt2)

  real(dp), intent(in)  :: occupation(:, :), energy(:, :)
  real(dp), intent(in)  :: c_matrix(:, :, :)
  real(dp), intent(out) :: p_matrix_pt2(:, :, :)
  !=====
  integer              :: nstate
  integer              :: istate, jstate, kstate
  integer              :: astate, bstate, cstate
  integer              :: pqspin
  real(dp)             :: denom1, denom2
  real(dp)             :: num1, num2
  !=====

  nstate = SIZE(occupation, DIM=1)

  call start_clock(timing_mbpt_dm)

  write(stdout, '(/,a)') ' Calculate the PT2 density matrix'

  if( nspin /= 1 ) call die('pt2_density_matrix: only implemented for spin restricted calculations')

  if(has_auxil_basis) then
    call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)
  else
    call calculate_eri_4center_mo_uks(c_matrix, ncore_G+1, nvirtual_G-1)
  endif


  ! Full calculation of the PT2 density matrix
  p_matrix_pt2(:, :, :) = 0.0_dp
  ! so far, only spin-restricted calculation are possible
  pqspin = 1

  ! A1 P_ij sum over k and a,b
  do istate=ncore_G+1, nhomo_G
    do jstate=ncore_G+1, nhomo_G
      do astate=nhomo_G+1, nvirtual_G-1
        if( MODULO( astate - (nhomo_G+1), poorman%nproc ) /= poorman%rank ) cycle
        do bstate=nhomo_G+1, nvirtual_G-1
          do kstate=ncore_G+1, nhomo_G

            denom1 = energy(astate, pqspin) + energy(bstate, pqspin) - energy(istate, pqspin) - energy(kstate, pqspin)
            denom2 = energy(astate, pqspin) + energy(bstate, pqspin) - energy(jstate, pqspin) - energy(kstate, pqspin)

            num1 = 2.0_dp * evaluate_eri_mo(istate, astate, pqspin, kstate, bstate, pqspin) &
                       - evaluate_eri_mo(istate, bstate, pqspin, kstate, astate, pqspin)
            num2 = 2.0_dp * evaluate_eri_mo(jstate, astate, pqspin, kstate, bstate, pqspin)

            p_matrix_pt2(istate, jstate, pqspin) = p_matrix_pt2(istate, jstate, pqspin)  &
                  - num1 * num2 / ( denom1 * denom2 )

          enddo
        enddo
      enddo
    enddo
  enddo

  ! A2 P_ab sum over  i,j  and c
  do astate=nhomo_G+1, nvirtual_G-1
    do bstate=nhomo_G+1, nvirtual_G-1
      do cstate=nhomo_G+1, nvirtual_G-1
        if( MODULO( cstate - (nhomo_G+1), poorman%nproc ) /= poorman%rank ) cycle
        do istate=ncore_G+1, nhomo_G
          do jstate=ncore_G+1, nhomo_G

            denom1 = energy(istate, pqspin) + energy(jstate, pqspin) - energy(astate, pqspin) - energy(cstate, pqspin)
            denom2 = energy(istate, pqspin) + energy(jstate, pqspin) - energy(bstate, pqspin) - energy(cstate, pqspin)

            num1 = 2.0_dp * evaluate_eri_mo(astate, istate, pqspin, cstate, jstate, pqspin) &
                       - evaluate_eri_mo(astate, jstate, pqspin, cstate, istate, pqspin)
            num2 = 2.0_dp * evaluate_eri_mo(bstate, istate, pqspin, cstate, jstate, pqspin)

            p_matrix_pt2(astate, bstate, pqspin) = p_matrix_pt2(astate, bstate, pqspin)  &
                  + num1 * num2 / ( denom1 * denom2 )

          enddo
        enddo
      enddo
    enddo
  enddo

  ! A3    P_cj  sum over i, and a,b
  ! A4    P_jc  sum over i, and a,b
  do cstate=nhomo_G+1, nvirtual_G-1
    do jstate=ncore_G+1, nhomo_G
      do astate=nhomo_G+1, nvirtual_G-1
        if( MODULO( astate - (nhomo_G+1), poorman%nproc ) /= poorman%rank ) cycle
        do bstate=nhomo_G+1, nvirtual_G-1
          do istate=ncore_G+1, nhomo_G
            denom1 = energy(jstate, pqspin) + energy(istate, pqspin) - energy(astate, pqspin) - energy(bstate, pqspin)
            denom2 = energy(jstate, pqspin) - energy(cstate, pqspin)
            num1 = 2.0_dp * evaluate_eri_mo(jstate, astate, pqspin, istate, bstate, pqspin) &
                - evaluate_eri_mo(jstate, bstate, pqspin, istate, astate, pqspin)
            num2 = 2.0_dp * evaluate_eri_mo(astate, cstate, pqspin, bstate, istate, pqspin)

            p_matrix_pt2(cstate, jstate, pqspin) = p_matrix_pt2(cstate, jstate, pqspin)  &
                                            + num1 * num2 / ( denom1 * denom2 )
            p_matrix_pt2(jstate, cstate, pqspin) = p_matrix_pt2(jstate, cstate, pqspin)  &
                                            + num1 * num2 / ( denom1 * denom2 )
          enddo
        enddo
      enddo
    enddo
  enddo

  ! A5   P_bk  sum over i,j,a
  ! A6   P_kb  sum over i,j,a
  do bstate=nhomo_G+1, nvirtual_G-1
    do kstate=ncore_G+1, nhomo_G
      do astate=nhomo_G+1, nvirtual_G-1
        if( MODULO( astate - (nhomo_G+1), poorman%nproc ) /= poorman%rank ) cycle
        do istate=ncore_G+1, nhomo_G
          do jstate=ncore_G+1, nhomo_G
            denom1 = energy(jstate, pqspin) + energy(istate, pqspin) - energy(astate, pqspin) - energy(bstate, pqspin)
            denom2 = energy(kstate, pqspin) - energy(bstate, pqspin)
            num1 = 2.0_dp * evaluate_eri_mo(jstate, astate, pqspin, istate, bstate, pqspin) &
                - evaluate_eri_mo(jstate, bstate, pqspin, istate, astate, pqspin)
            num2 = 2.0_dp * evaluate_eri_mo(istate, kstate, pqspin, jstate, astate, pqspin)

            p_matrix_pt2(bstate, kstate, pqspin) = p_matrix_pt2(bstate, kstate, pqspin)  &
                                            - num1 * num2 / ( denom1 * denom2 )
            p_matrix_pt2(kstate, bstate, pqspin) = p_matrix_pt2(kstate, bstate, pqspin)  &
                                            - num1 * num2 / ( denom1 * denom2 )
          enddo
        enddo
      enddo
    enddo
  enddo
  call poorman%sum(p_matrix_pt2)


  if(has_auxil_basis) then
    call destroy_eri_3center_mo()
  else
    call destroy_eri_4center_mo_uks()
  endif

  call stop_clock(timing_mbpt_dm)

end subroutine pt2_density_matrix


!=========================================================================
subroutine onering_density_matrix(occupation, energy, c_matrix, p_matrix_pt2)

  real(dp), intent(in)   :: occupation(:, :), energy(:, :)
  real(dp), intent(in)   :: c_matrix(:, :, :)
  real(dp), intent(out)  :: p_matrix_pt2(:, :, :)
  !=====
  integer              :: nstate
  integer              :: istate, jstate, kstate
  integer              :: astate, bstate, cstate
  integer              :: pqspin
  real(dp)             :: denom1, denom2
  real(dp)             :: num1, num2
  !=====

  call start_clock(timing_mbpt_dm)

  nstate = SIZE(occupation, DIM=1)

  write(stdout, '(/,a)') ' Calculate the 1-ring density matrix'

  if( nspin /= 1 ) call die('pt2_density_matrix: only implemented for spin restricted calculations')

  if(has_auxil_basis) then
    call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)
  else
    call calculate_eri_4center_mo_uks(c_matrix, ncore_G+1, nvirtual_G-1)
  endif


  ! Full calculation of the 1-ring density matrix

  p_matrix_pt2(:, :, :) = 0.0_dp
  pqspin = 1

  ! A1
  do istate=ncore_G+1, nhomo_G
    do jstate=ncore_G+1, nhomo_G
      do kstate=ncore_G+1, nhomo_G
        do astate=nhomo_G+1, nvirtual_G-1
          do bstate=nhomo_G+1, nvirtual_G-1

            denom1 = energy(astate, pqspin) + energy(bstate, pqspin) - energy(istate, pqspin) - energy(kstate, pqspin)
            denom2 = energy(astate, pqspin) + energy(bstate, pqspin) - energy(jstate, pqspin) - energy(kstate, pqspin)

            num1 = 2.0_dp * evaluate_eri_mo(istate, astate, pqspin, kstate, bstate, pqspin)
            num2 = 2.0_dp * evaluate_eri_mo(jstate, astate, pqspin, kstate, bstate, pqspin)

            p_matrix_pt2(istate, jstate, pqspin) = p_matrix_pt2(istate, jstate, pqspin)  &
                  - num1 * num2 / ( denom1 * denom2 )

          enddo
        enddo
      enddo
    enddo
  enddo

  ! A2
  do astate=nhomo_G+1, nvirtual_G-1
    do bstate=nhomo_G+1, nvirtual_G-1
      do cstate=nhomo_G+1, nvirtual_G-1
        do istate=ncore_G+1, nhomo_G
          do jstate=ncore_G+1, nhomo_G

            denom1 = energy(istate, pqspin) + energy(jstate, pqspin) - energy(astate, pqspin) - energy(cstate, pqspin)
            denom2 = energy(istate, pqspin) + energy(jstate, pqspin) - energy(bstate, pqspin) - energy(cstate, pqspin)

            num1 = 2.0_dp * evaluate_eri_mo(astate, istate, pqspin, cstate, jstate, pqspin)
            num2 = 2.0_dp * evaluate_eri_mo(bstate, istate, pqspin, cstate, jstate, pqspin)

            p_matrix_pt2(astate, bstate, pqspin) = p_matrix_pt2(astate, bstate, pqspin)  &
                  + num1 * num2 / ( denom1 * denom2 )

          enddo
        enddo
      enddo
    enddo
  enddo

  ! A3    P_cj  sum over i,a,b
  ! A4    P_jc  sum over i,a,b
  do cstate=nhomo_G+1, nvirtual_G-1
    do jstate=ncore_G+1, nhomo_G
      do astate=nhomo_G+1, nvirtual_G-1
        do bstate=nhomo_G+1, nvirtual_G-1
          do istate=ncore_G+1, nhomo_G
            denom1 = energy(jstate, pqspin) + energy(istate, pqspin) - energy(astate, pqspin) - energy(bstate, pqspin)
            denom2 = energy(jstate, pqspin) - energy(cstate, pqspin)
            num1 = 2.0_dp * evaluate_eri_mo(jstate, astate, pqspin, istate, bstate, pqspin)
            num2 = 2.0_dp * evaluate_eri_mo(astate, cstate, pqspin, bstate, istate, pqspin)

            p_matrix_pt2(cstate, jstate, pqspin) = p_matrix_pt2(cstate, jstate, pqspin)  &
                                            + num1 * num2 / ( denom1 * denom2 )
            p_matrix_pt2(jstate, cstate, pqspin) = p_matrix_pt2(jstate, cstate, pqspin)  &
                                            + num1 * num2 / ( denom1 * denom2 )
          enddo
        enddo
      enddo
    enddo
  enddo

  ! A5   P_bk  sum over i,j,a
  ! A6   P_kb  sum over i,j,a
  do bstate=nhomo_G+1, nvirtual_G-1
    do kstate=ncore_G+1, nhomo_G
      do astate=nhomo_G+1, nvirtual_G-1
        do istate=ncore_G+1, nhomo_G
          do jstate=ncore_G+1, nhomo_G
            denom1 = energy(jstate, pqspin) + energy(istate, pqspin) - energy(astate, pqspin) - energy(bstate, pqspin)
            denom2 = energy(kstate, pqspin) - energy(bstate, pqspin)
            num1 = 2.0_dp * evaluate_eri_mo(jstate, astate, pqspin, istate, bstate, pqspin)
            num2 = 2.0_dp * evaluate_eri_mo(istate, kstate, pqspin, jstate, astate, pqspin)

            p_matrix_pt2(bstate, kstate, pqspin) = p_matrix_pt2(bstate, kstate, pqspin)  &
                                            - num1 * num2 / ( denom1 * denom2 )
            p_matrix_pt2(kstate, bstate, pqspin) = p_matrix_pt2(kstate, bstate, pqspin)  &
                                            - num1 * num2 / ( denom1 * denom2 )
          enddo
        enddo
      enddo
    enddo
  enddo


  if(has_auxil_basis) then
    call destroy_eri_3center_mo()
  else
    call destroy_eri_4center_mo_uks()
  endif

  call stop_clock(timing_mbpt_dm)

end subroutine onering_density_matrix


!=========================================================================
subroutine gw_density_matrix(occupation, energy, c_matrix, wpol, p_matrix_gw)

  real(dp), intent(in)                :: occupation(:, :), energy(:, :)
  real(dp), intent(in)                :: c_matrix(:, :, :)
  type(spectral_function), intent(in) :: wpol
  real(dp), intent(out)               :: p_matrix_gw(:, :, :)
  !=====
  integer  :: nstate
  integer  :: pstate, qstate
  integer  :: istate, jstate
  integer  :: astate, bstate
  integer  :: pqspin
  integer  :: spole
  integer  :: npole_local, spole_local
  integer  :: nstate_occ, nstate_virt
  real(dp), allocatable :: w_s_occ(:, :), w_s_virt(:, :)
  real(dp), allocatable :: w_s_occ_local(:, :), w_s_virt_local(:, :)
  !=====

  call start_clock(timing_mbpt_dm)

  nstate = SIZE(occupation, DIM=1)

  write(stdout, '(/,a)') ' Calculate the GW density matrix'

  if( nspin /= 1 ) call die('gw_density_matrix: only implemented for spin restricted calculations')
  if( .NOT. has_auxil_basis)  call die('gw_density_matrix: only implemented without RI')

  if(has_auxil_basis) then
    call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)
  else
    call calculate_eri_4center_mo_uks(c_matrix, ncore_G+1, nvirtual_G-1)
  endif


  ! First order calculation of the GW density matrix
  p_matrix_gw(:, :, :) = 0.0_dp
  pqspin = 1

  nstate_occ  = nhomo_G - ncore_G
  nstate_virt = nvirtual_G - nhomo_G - 1
  allocate(w_s_occ(wpol%npole_reso, ncore_G+1:nhomo_G))
  allocate(w_s_virt(wpol%npole_reso, nhomo_G+1:nvirtual_G-1))
  npole_local = NUMROC(wpol%npole_reso, 1, auxil%rank, 0, auxil%nproc)
  allocate(w_s_occ_local(npole_local, ncore_G+1:nhomo_G))
  allocate(w_s_virt_local(npole_local, nhomo_G+1:nvirtual_G-1))

  do astate=nhomo_G+1, nvirtual_G-1
    if( MODULO( astate - (nhomo_G+1) , poorman%nproc ) /= poorman%rank ) cycle

    ! A1
    !w_s_occ(:,ncore_G+1:nhomo_G) = MATMUL( TRANSPOSE(wpol%w_s(:,:)) , eri_3center_mo(:,ncore_G+1:nhomo_G,astate,pqspin) )
    call DGEMM('T', 'N', wpol%npole_reso, nstate_occ,nauxil_local, &
                          1.0d0, wpol%w_s, nauxil_local, &
                                eri_3center_mo(1, ncore_G+1, astate, pqspin), nauxil_local, &
                          0.0_dp, w_s_occ(1, ncore_G+1), wpol%npole_reso)
    call auxil%sum(w_s_occ)

    spole_local = 0
    do spole=1, wpol%npole_reso
      if( MODULO( spole-1 , auxil%nproc ) /= auxil%rank ) cycle
      spole_local = spole_local + 1
      do jstate=ncore_G+1, nhomo_G
        w_s_occ_local(spole_local, jstate) = w_s_occ(spole, jstate) &
                                             / ( energy(jstate, pqspin) - energy(astate, pqspin) - wpol%pole(spole) )
      enddo
    enddo

    call DSYRK('U', 'T', nstate_occ, npole_local, -2.0d0, w_s_occ_local,npole_local,&
               1.0d0, p_matrix_gw(ncore_G+1, ncore_G+1, pqspin), nstate)


    ! A3    P_cj  sum over i,a,b
    ! A4    P_jc  sum over i,a,b     ! not actually calculated, but included through the symmetrization step
    !w_s_virt(:,nhomo_G+1:nvirtual_G-1) = MATMUL( TRANSPOSE(wpol%w_s(:,:)) , eri_3center_mo(:,nhomo_G+1:nvirtual_G-1,astate,pqspin) )
    call DGEMM('T', 'N', wpol%npole_reso, nstate_virt,nauxil_local, &
                          1.0d0, wpol%w_s, nauxil_local,  &
                                eri_3center_mo(1, nhomo_G+1, astate, pqspin), nauxil_local, &
                          0.0_dp, w_s_virt(1, nhomo_G+1), wpol%npole_reso)
    call auxil%sum(w_s_virt)

    spole_local = 0
    do spole=1, wpol%npole_reso
      if( MODULO( spole-1 , auxil%nproc ) /= auxil%rank ) cycle
      spole_local = spole_local + 1
      do bstate=nhomo_G+1, nvirtual_G-1
        w_s_virt_local(spole_local, bstate) = w_s_virt(spole, bstate)
      enddo
    enddo

    call DGEMM('T', 'N', nstate_occ, nstate_virt, npole_local, 2.0d0, w_s_occ_local,npole_local, &
                                                                w_s_virt_local, npole_local, &
                                                          1.0d0, p_matrix_gw(ncore_G+1, nhomo_G+1, pqspin), nstate)
  enddo

  do istate=ncore_G+1, nhomo_G
    if( MODULO( istate - (ncore_G+1) , poorman%nproc ) /= poorman%rank ) cycle

    ! A2
    !w_s_virt(:,nhomo_G+1:nvirtual_G-1) = MATMUL( TRANSPOSE(wpol%w_s(:,:)) , eri_3center_mo(:,nhomo_G+1:nvirtual_G-1,istate,pqspin) )
    call DGEMM('T', 'N', wpol%npole_reso, nstate_virt,nauxil_local, &
                          1.0d0, wpol%w_s, nauxil_local,  &
                                eri_3center_mo(1, nhomo_G+1, istate, pqspin), nauxil_local, &
                          0.0_dp, w_s_virt(1, nhomo_G+1), wpol%npole_reso)
    call auxil%sum(w_s_virt)

    spole_local = 0
    do spole=1, wpol%npole_reso
      if( MODULO( spole-1 , auxil%nproc ) /= auxil%rank ) cycle
      spole_local = spole_local + 1
      do bstate=nhomo_G+1, nvirtual_G-1
        w_s_virt_local(spole_local, bstate) = w_s_virt(spole, bstate) &
                                             / ( energy(istate, pqspin) - energy(bstate, pqspin) - wpol%pole(spole) )
      enddo
    enddo

    call DSYRK('U', 'T', nstate_virt, npole_local, 2.0d0, w_s_virt_local,npole_local, &
               1.0d0, p_matrix_gw(nhomo_G+1, nhomo_G+1, pqspin), nstate)

    ! A5   P_bk  sum over i,j,a
    ! A6   P_kb  sum over i,j,a   ! not actually calculated, but included through the symmetrization step
    !w_s_occ(:,ncore_G+1:nhomo_G)       = MATMUL( TRANSPOSE(wpol%w_s(:,:)) , eri_3center_mo(:,ncore_G+1:nhomo_G,istate,pqspin) )
    call DGEMM('T', 'N', wpol%npole_reso, nstate_occ,nauxil_local, &
                          1.0d0, wpol%w_s, nauxil_local, &
                                eri_3center_mo(1, ncore_G+1, istate, pqspin), nauxil_local, &
                          0.0_dp, w_s_occ(1, ncore_G+1), wpol%npole_reso)
    call auxil%sum(w_s_occ)

    spole_local = 0
    do spole=1, wpol%npole_reso
      if( MODULO( spole-1 , auxil%nproc ) /= auxil%rank ) cycle
      spole_local = spole_local + 1
      do jstate=ncore_G+1, nhomo_G
        w_s_occ_local(spole_local, jstate) = w_s_occ(spole, jstate)
      enddo
    enddo

    call DGEMM('T', 'N', nstate_occ, nstate_virt, npole_local, -2.0d0, w_s_occ_local,npole_local, &
                                                                 w_s_virt_local, npole_local, &
                                                           1.0d0, p_matrix_gw(ncore_G+1, nhomo_G+1, pqspin), nstate)
  enddo

  ! A common factor 1/(e_j-e_c) is to be added for the occupied-virtual block (terms A3,A4,A5,A6)
  do bstate=nhomo_G+1, nvirtual_G-1
    do jstate=ncore_G+1, nhomo_G
      p_matrix_gw(jstate, bstate, pqspin) = p_matrix_gw(jstate, bstate, pqspin) &
                                            / ( energy(jstate, pqspin) - energy(bstate, pqspin) )
    enddo
  enddo

  call world%sum(p_matrix_gw)

  deallocate(w_s_occ)
  deallocate(w_s_virt)
  deallocate(w_s_occ_local)
  deallocate(w_s_virt_local)

  ! Symmetrization of the p_matrix_gw here
  ! Only the upper triangle was set up before
  do pstate=1, nstate
    do qstate=pstate+1, nstate
      p_matrix_gw(qstate, pstate, pqspin) = p_matrix_gw(pstate, qstate, pqspin)
    enddo
  enddo


  if(has_auxil_basis) then
    call destroy_eri_3center_mo()
  else
    call destroy_eri_4center_mo_uks()
  endif

  call stop_clock(timing_mbpt_dm)


end subroutine gw_density_matrix


!=========================================================================
subroutine gw_density_matrix_imag(occupation, energy, c_matrix, wpol, p_matrix_gw)

  real(dp), intent(in)                :: occupation(:, :), energy(:, :)
  real(dp), intent(in)                :: c_matrix(:, :, :)
  type(spectral_function), intent(in) :: wpol
  real(dp), intent(out)               :: p_matrix_gw(:, :, :)
  !=====
  real(dp), parameter  :: alpha_quad = 1.0_dp
  real(dp), parameter  :: beta_quad = 1.0_dp
  integer              :: nstate
  integer              :: iomegas, iomega
  integer              :: info
  real(dp), allocatable :: eri3_sca_p(:, :), eri3_sca_q(:, :)
  real(dp), allocatable :: chi_eri3_sca_q(:, :)
  real(dp), allocatable :: omega_sigma(:), weight_sigma(:)
  real(dp)             :: v_chi_v_pq, mu
  integer              :: desc_eri3_t(NDEL)
  integer              :: iprow, ipcol, nprow, npcol
  integer              :: desc_eri3_final(NDEL)
  integer              :: meri3, neri3
  integer              :: mstate, pstate, qstate, pqspin
  integer              :: mrange, mlocal
  !=====


  if( .NOT. has_auxil_basis ) then
    call die('gw_density_matrix_imag: requires an auxiliary basis')
  endif

  call start_clock(timing_mbpt_dm)

  nstate = SIZE(occupation, DIM=1)

  write(stdout, '(/,1x,a)') 'GW density matrix from a grid of imaginary frequencies'

  nprow = 1
  npcol = 1
#if defined(HAVE_SCALAPACK)
  ! Get the processor grid included in the input wpol%desc_chi
  call BLACS_GRIDINFO(wpol%desc_chi(CTXT_), nprow, npcol, iprow, ipcol)
  write(stdout, '(1x,a,i4,a,i4)') 'SCALAPACK grid', nprow, ' x ', npcol
#endif


  if( has_auxil_basis ) call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)


  mrange = nvirtual_G - ncore_G - 1

  meri3 = NUMROC(nauxil_global, wpol%desc_chi(MB_), iprow, wpol%desc_chi(RSRC_), nprow)
  neri3 = NUMROC(mrange        , wpol%desc_chi(NB_), ipcol, wpol%desc_chi(CSRC_), npcol)
  call DESCINIT(desc_eri3_final, nauxil_global, mrange, wpol%desc_chi(MB_), wpol%desc_chi(NB_), &
                wpol%desc_chi(RSRC_), wpol%desc_chi(CSRC_), wpol%desc_chi(CTXT_), MAX(1, meri3), info)

  call clean_allocate('TMP 3-center MO integrals', eri3_sca_p, meri3, neri3)
  call clean_allocate('TMP 3-center MO integrals', eri3_sca_q, meri3, neri3)
  call clean_allocate('TMP 3-center MO integrals', chi_eri3_sca_q, meri3, neri3)

  call DESCINIT(desc_eri3_t, nauxil_global, mrange, MB_eri3_mo, NB_eri3_mo, &
                first_row, first_col, cntxt_eri3_mo, MAX(1, nauxil_local), info)


  !
  ! Set up the imaginary frequency grid
  !
  allocate(omega_sigma(nomega_sigma), weight_sigma(nomega_sigma))
  call coeffs_gausslegint(0.0_dp, 1.0_dp, omega_sigma, weight_sigma, nomega_sigma)

  write(stdout, '(/,1x,a)') 'Numerical integration on a grid along the imaginary axis'
  ! Variable change [0,1] -> [0,+\inf[
  write(stdout, '(a)') '    #    Frequencies (eV)    Quadrature weights'
  do iomegas=1, nomega_sigma
    weight_sigma(iomegas) = weight_sigma(iomegas) / ( 2.0_dp**alpha_quad - 1.0_dp ) &
                           * alpha_quad * (1.0_dp -  omega_sigma(iomegas))**(-alpha_quad-1.0_dp) * beta_quad
    omega_sigma(iomegas)  = 1.0_dp / ( 2.0_dp**alpha_quad - 1.0_dp ) &
                             * ( 1.0_dp / (1.0_dp - omega_sigma(iomegas))**alpha_quad - 1.0_dp ) * beta_quad
    write(stdout, '(i5,2(2x,f14.6))') iomegas, omega_sigma(iomegas)*Ha_eV, weight_sigma(iomegas)
  enddo

  !
  ! Find the HOMO-LUMO gap
  mu = ( MINVAL(energy(nhomo_G+1, :)) + MAXVAL(energy(nhomo_G, :)) ) / 2.0_dp
  write(stdout, '(1x,a,f12.6)') 'Center of the HOMO-LUMO gap (eV): ', mu*Ha_eV

  p_matrix_gw(:, :, :) = 0.0_dp

  do pqspin=1, nspin
    do qstate=nsemin, nsemax

      eri3_sca_q(:, 1:mrange) = eri_3center_mo(:, ncore_G+1:nvirtual_G-1, qstate, pqspin)


      do iomega=1, wpol%nomega

        call DGEMM('N', 'N', nauxil_global, mrange,nauxil_global,  &
                   1.0_dp, wpol%chi(:, :, iomega), nauxil_global,    &
                          eri3_sca_q          , nauxil_global,    &
                   0.0_dp, chi_eri3_sca_q      , nauxil_global)

        do pstate=qstate, nsemax
          eri3_sca_p(:, 1:mrange) = eri_3center_mo(:, ncore_G+1:nvirtual_G-1, pstate, pqspin)

          do mlocal=1, neri3
            mstate = INDXL2G(mlocal, wpol%desc_chi(NB_), ipcol, wpol%desc_chi(CSRC_), npcol) + ncore_G

            v_chi_v_pq = DOT_PRODUCT( eri3_sca_p(:, mlocal) , chi_eri3_sca_q(:, mlocal) )

            ! Factor 2 and the real part come from the integration of Sigma
            ! for negative AND positive imaginary frequencies
            ! Sigma_pq(-iw) = Sigma_qp(iw)^*
            p_matrix_gw(pstate, qstate, pqspin) = p_matrix_gw(pstate, qstate, pqspin) &
                 - spin_fact * 2.0_dp  / (2.0_dp *  pi)**2  &
                       * SUM( wpol%weight_quad(iomega) * weight_sigma(:) * v_chi_v_pq  &
                    * REAL( (  1.0_dp / ( im * omega_sigma(:) - (energy(mstate, pqspin)-mu) + wpol%omega(iomega) )    &
                             + 1.0_dp / ( im * omega_sigma(:) - (energy(mstate, pqspin)-mu) - wpol%omega(iomega) )  ) &
                            / ( im * omega_sigma(:) - (energy(pstate, pqspin)-mu) )   &
                            / ( im * omega_sigma(:) - (energy(qstate, pqspin)-mu) ) , dp ) )
          enddo

        enddo

      enddo

    enddo
  enddo
  call world%sum(p_matrix_gw)

  !debug
  call dump_out_matrix(.FALSE., 'P matrix', p_matrix_gw)
  deallocate(omega_sigma, weight_sigma)


  ! Symmetrization of the p_matrix_gw here
  ! Only the upper triangle was set up before
  do pqspin=1, nspin
    do qstate=1, nstate
      do pstate=qstate+1, nstate
        p_matrix_gw(qstate, pstate, pqspin) = p_matrix_gw(pstate, qstate, pqspin)
      enddo
    enddo
  enddo


  call clean_deallocate('TMP 3-center MO integrals', eri3_sca_p)
  call clean_deallocate('TMP 3-center MO integrals', eri3_sca_q)
  call clean_deallocate('TMP 3-center MO integrals', chi_eri3_sca_q)

  call destroy_eri_3center_mo()

  call stop_clock(timing_mbpt_dm)

end subroutine gw_density_matrix_imag


!=========================================================================
subroutine gw_density_matrix_dyson_imag(occupation, energy, c_matrix, wpol, p_matrix_gw)

  real(dp), intent(in)                :: occupation(:, :), energy(:, :)
  real(dp), intent(in)                :: c_matrix(:, :, :)
  type(spectral_function), intent(in) :: wpol
  real(dp), intent(out)               :: p_matrix_gw(:, :, :)
  !=====
  real(dp), parameter  :: alpha_quad = 1.0_dp
  real(dp), parameter  :: beta_quad = 1.0_dp
  integer              :: nstate
  integer              :: iomegas, iomega
  integer              :: info
  real(dp), allocatable :: eri3_sca_p(:, :), eri3_sca_q(:, :)
  real(dp), allocatable :: chi_eri3_sca_q(:, :)
  real(dp), allocatable :: omega_sigma(:), weight_sigma(:)
  complex(dp), allocatable :: m_matrix(:, :, :, :), sigma_c_g0(:, :, :)
  real(dp)             :: v_chi_v_pq, mu, ec_gm
  integer              :: desc_eri3_t(NDEL)
  integer              :: iprow, ipcol, nprow, npcol
  integer              :: desc_eri3_final(NDEL)
  integer              :: meri3, neri3
  integer              :: mstate, pstate, qstate, pqspin
  integer              :: mrange, mlocal
  !=====


  if( .NOT. has_auxil_basis ) then
    call die('gw_density_matrix_imag: requires an auxiliary basis')
  endif

  call start_clock(timing_mbpt_dm)

  nstate = SIZE(occupation, DIM=1)

  write(stdout, '(/,1x,a)') 'GW density matrix from a grid of imaginary frequencies'

  nprow = 1
  npcol = 1
#if defined(HAVE_SCALAPACK)
  ! Get the processor grid included in the input wpol%desc_chi
  call BLACS_GRIDINFO(wpol%desc_chi(CTXT_), nprow, npcol, iprow, ipcol)
  write(stdout, '(1x,a,i4,a,i4)') 'SCALAPACK grid', nprow, ' x ', npcol
#endif


  if( has_auxil_basis ) call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)


  mrange = nvirtual_G - ncore_G - 1

  meri3 = NUMROC(nauxil_global, wpol%desc_chi(MB_), iprow, wpol%desc_chi(RSRC_), nprow)
  neri3 = NUMROC(mrange        , wpol%desc_chi(NB_), ipcol, wpol%desc_chi(CSRC_), npcol)
  call DESCINIT(desc_eri3_final, nauxil_global, mrange, wpol%desc_chi(MB_), wpol%desc_chi(NB_), &
                wpol%desc_chi(RSRC_), wpol%desc_chi(CSRC_), wpol%desc_chi(CTXT_), MAX(1, meri3), info)

  call clean_allocate('TMP 3-center MO integrals', eri3_sca_p, meri3, neri3)
  call clean_allocate('TMP 3-center MO integrals', eri3_sca_q, meri3, neri3)
  call clean_allocate('TMP 3-center MO integrals', chi_eri3_sca_q, meri3, neri3)

  call DESCINIT(desc_eri3_t, nauxil_global, mrange, MB_eri3_mo, NB_eri3_mo, &
                first_row, first_col, cntxt_eri3_mo, MAX(1, nauxil_local), info)


  !
  ! Set up the imaginary frequency grid
  !
  allocate(omega_sigma(nomega_sigma), weight_sigma(nomega_sigma))
  call coeffs_gausslegint(0.0_dp, 1.0_dp, omega_sigma, weight_sigma, nomega_sigma)

  write(stdout, '(/,1x,a)') 'Numerical integration on a grid along the imaginary axis'
  ! Variable change [0,1] -> [0,+\inf[
  write(stdout, '(a)') '    #    Frequencies (eV)    Quadrature weights'
  do iomegas=1, nomega_sigma
    weight_sigma(iomegas) = weight_sigma(iomegas) / ( 2.0_dp**alpha_quad - 1.0_dp ) &
                           * alpha_quad * (1.0_dp -  omega_sigma(iomegas))**(-alpha_quad-1.0_dp) * beta_quad
    omega_sigma(iomegas)  = 1.0_dp / ( 2.0_dp**alpha_quad - 1.0_dp ) &
                           * ( 1.0_dp / (1.0_dp - omega_sigma(iomegas))**alpha_quad - 1.0_dp ) * beta_quad
    write(stdout, '(i5,2(2x,f14.6))') iomegas, omega_sigma(iomegas)*Ha_eV, weight_sigma(iomegas)
  enddo

  !
  ! Find the HOMO-LUMO gap
  mu = ( MINVAL(energy(nhomo_G+1, :)) + MAXVAL(energy(nhomo_G, :)) ) / 2.0_dp
  write(stdout, '(1x,a,f12.6)') 'Center of the HOMO-LUMO gap (eV): ', mu*Ha_eV

  allocate(m_matrix(nsemin:nsemax, nsemin:nsemax, nomega_sigma, nspin))
  allocate(sigma_c_g0(nomega_sigma, nsemin:nsemax, nspin))
  m_matrix(:, :, :, :) = (0.0_dp, 0.0_dp)
  sigma_c_g0(:, :, :) = (0.0_dp, 0.0_dp)

  do pqspin=1, nspin
    do qstate=nsemin, nsemax

      eri3_sca_q(:, 1:mrange) = eri_3center_mo(:, ncore_G+1:nvirtual_G-1, qstate, pqspin)


      do iomega=1, wpol%nomega

        call DGEMM('N', 'N', nauxil_global, mrange,nauxil_global,  &
                   1.0_dp, wpol%chi(:, :, iomega), nauxil_global,    &
                          eri3_sca_q          , nauxil_global,    &
                   0.0_dp, chi_eri3_sca_q      , nauxil_global)

        do pstate=nsemin, nsemax
          eri3_sca_p(:, 1:mrange) = eri_3center_mo(:, ncore_G+1:nvirtual_G-1, pstate, pqspin)

          do mlocal=1, neri3
            mstate = INDXL2G(mlocal, wpol%desc_chi(NB_), ipcol, wpol%desc_chi(CSRC_), npcol) + ncore_G

            v_chi_v_pq = DOT_PRODUCT( eri3_sca_p(:, mlocal) , chi_eri3_sca_q(:, mlocal) )

            ! M will contain (1 - G_0 * Sigma)
            ! Sigma has a -1/(2*pi) prefactor
            ! then -G_0 * Sigma should be multiplied by 1/(2*pi)
            m_matrix(pstate, qstate, :, pqspin) = m_matrix(pstate, qstate, :, pqspin) &
                 + wpol%weight_quad(iomega) * v_chi_v_pq /  ( 2 * pi )  &
                    * (  1.0_dp / ( im * omega_sigma(:) - (energy(mstate, pqspin)-mu) + wpol%omega(iomega) )    &
                       + 1.0_dp / ( im * omega_sigma(:) - (energy(mstate, pqspin)-mu) - wpol%omega(iomega) )  ) &
                    / ( im * omega_sigma(:) - (energy(pstate, pqspin)-mu) )
          enddo

        enddo

      enddo

      sigma_c_g0(:, qstate, pqspin) = -m_matrix(qstate, qstate, :, pqspin)
      m_matrix(qstate, qstate, :, pqspin) = m_matrix(qstate, qstate, :, pqspin) + (1.0_dp, 0.0_dp)
    enddo
  enddo

  do pqspin=1, nspin
    do iomegas=1, nomega_sigma
      call invert(m_matrix(:, :, iomegas, pqspin))
    enddo
    !
    ! Remove G_0 since p_matrix_gw is to contain the *correction* to p_matrix
    do qstate=nsemin, nsemax
      m_matrix(qstate, qstate, :, pqspin) = m_matrix(qstate, qstate, :, pqspin) - (1.0_dp, 0.0_dp)
    enddo
  enddo

  !
  ! Perform the final omega integrals:
  ! delta P_MO = spin_fact /(2 pi i) int_{-inf}^{+inf} d(i omega)  (1-G_0*Sigma)^{-1} * G_0(\mu + i omega)
  !            = spin_fact /(2 pi i) int_{  0 }^{+inf} d(i omega) 2 * Re{ (1-G_0*Sigma)^{-1} * G_0(\mu + i omega) }
  ! Ec_GM      = spin_fact /(4 pi i) int_{-inf}^{+inf} d(i omega) Sigma(i omega) * G0(i omega)
  !            = spin_fact /(4 pi i) int_{  0 }^{+inf} d(i omega) 2 * Re{ Sigma(i omega) * G0(i omega) }
  ec_gm = 0.0_dp
  p_matrix_gw(:, :, :) = 0.0_dp

  do pqspin=1, nspin
    do qstate=nsemin, nsemax
      do pstate=nsemin, nsemax
        p_matrix_gw(pstate, qstate, pqspin) = &
           2.0_dp * spin_fact / ( 2.0_dp * pi )  &
            * REAL( SUM( m_matrix(pstate, qstate, :, pqspin) &
                    / ( im * omega_sigma(:) - (energy(qstate, pqspin)-mu) ) * weight_sigma(:) ) , dp)
      enddo
      ec_gm = ec_gm + &
         2.0_dp * spin_fact / ( 2.0_dp * pi )  &
         * REAL( SUM( sigma_c_g0(:, qstate, pqspin) * weight_sigma(:) ), dp)
    enddo
  enddo

  ec_gm = 0.5_dp * ec_gm
  write(stdout, '(/,a)')       ' Galitskii-Migdal formula on the Imag. axis 1/(4*pi) int [G0(iw) * Sigma_c(iw)] dw:'
  write(stdout, '(a,f19.10,/)') ' GM correlation energy (Ha): ', ec_gm

  deallocate(omega_sigma, weight_sigma, sigma_c_g0)


  ! Symmetrization of the p_matrix here
  ! Only the upper triangle was set up before
  do pqspin=1, nspin
    do qstate=1, nstate
      do pstate=qstate+1, nstate
        p_matrix_gw(qstate, pstate, pqspin) = p_matrix_gw(pstate, qstate, pqspin)
      enddo
    enddo
  enddo


  call clean_deallocate('TMP 3-center MO integrals', eri3_sca_p)
  call clean_deallocate('TMP 3-center MO integrals', eri3_sca_q)
  call clean_deallocate('TMP 3-center MO integrals', chi_eri3_sca_q)

  call destroy_eri_3center_mo()

  call stop_clock(timing_mbpt_dm)

end subroutine gw_density_matrix_dyson_imag


!=========================================================================
subroutine update_density_matrix(occupation, c_matrix, p_matrix_mo, p_matrix)

  real(dp), intent(in)    :: occupation(:, :)
  real(dp), intent(in)    :: c_matrix(:, :, :)
  real(dp), intent(in)    :: p_matrix_mo(:, :, :)
  real(dp), intent(inout) :: p_matrix(:, :, :)
  !=====
  integer              :: pstate, pqspin
  integer              :: nbf, nstate
  integer              :: file_density_matrix
  real(dp), allocatable :: p_matrix_tmp(:, :, :)
  !=====

  nbf    = SIZE(c_matrix, DIM=1)
  nstate = SIZE(c_matrix, DIM=2)

  ! Input density matrix (p_matrix_mo) is the change in the density matrix on the MO
  ! Input density matrix (p_matrix) is the Fock density matrix on the AO
  ! Output density matrix (p_matrix) is the full density matrix on the AO

  !
  ! If input Fock density (p_matrix) is zero, then assume an Hartree-Fock SCF calculation
  if( ALL( ABS(p_matrix(:, :, :)) < 1.0e-6_dp ) ) then
    ! Add the SCF density matrix to get to the total density matrix
    allocate(p_matrix_tmp(nstate, nstate, nspin))
    p_matrix_tmp(:, :, :) = p_matrix_mo(:, :, :)
    do pstate=1, nstate
      p_matrix_tmp(pstate, pstate, :) = p_matrix_tmp(pstate, pstate, :) + occupation(pstate, :)
    enddo
    ! Transform from MO to AO
    call p_mo_to_ao(c_matrix, p_matrix_tmp, p_matrix)
    deallocate(p_matrix_tmp)

  else

    write(stdout, *) 'An input Fock density matrix was provided. Use it now!'
    allocate(p_matrix_tmp(nbf, nbf, nspin))
    ! Transform from MO to AO
    call p_mo_to_ao(c_matrix, p_matrix_mo, p_matrix_tmp)
    p_matrix(:, :, :) = p_matrix(:, :, :) + p_matrix_tmp(:, :, :)
    deallocate(p_matrix_tmp)

  endif

  !
  ! Write down the DENSITY_MATRIX file
  ! that contains the complete p_matrix on the AO basis
  !
  if( print_density_matrix_ .AND. is_iomaster ) then
    write(stdout, '(1x,a)') 'Write DENSITY_MATRIX file'
    open(newunit=file_density_matrix, file='DENSITY_MATRIX', form='unformatted', action='write')
    do pqspin=1, nspin
      write(file_density_matrix) p_matrix(:, :, pqspin)
    enddo
    close(file_density_matrix)
  endif


end subroutine update_density_matrix


!=========================================================================
! in input:
!   p_matrix: Fock density-matrix in AO
!   p_matrix_mo: correlated density-matrix correction in MO
! 
! in output:
!   p_matrix: 0.0
!   p_matrix_mo: Complete relaxed density-matrix in MO
!
subroutine prepare_cederbaum(c_matrix, s_matrix, occupation, p_matrix, p_matrix_mo)
  real(dp), intent(in) :: c_matrix(:, :, :), s_matrix(:, :), occupation(:, :)
  real(dp), intent(inout) :: p_matrix(:, :, :), p_matrix_mo(:, :, :)
  !=====
  integer :: nstate, istate
  real(dp), allocatable :: p_matrix_mo_tmp(:, :, :)
  !=====

  nstate = SIZE(c_matrix, DIM=2)

  ! First convert the Fock density matrix from AO -> MO
  allocate(p_matrix_mo_tmp(nstate, nstate, nspin))
  call p_ao_to_mo(c_matrix, s_matrix, p_matrix, p_matrix_mo_tmp)

  ! Remove the mean-field part since we just want the Fock correction
  do istate=1, nhomo_G
    p_matrix_mo_tmp(istate, istate, :) = p_matrix_mo_tmp(istate, istate, :) - occupation(istate, :)
  enddo
  ! Add p_matrix_mo_tmp to p_matrix_mo
  p_matrix_mo(:, :, :) = p_matrix_mo(:, :, :) + p_matrix_mo_tmp(:, :, :)

  deallocate(p_matrix_mo_tmp)

  ! Set p_matrix to zero here so that it is not added again later
  ! in "update_density_matrix" 
  p_matrix(:, :, :) = 0.0_dp


end subroutine prepare_cederbaum


!=========================================================================
subroutine cederbaum_blas(energy, c_matrix, hfock_ao, p_matrix_mo)

  real(dp), intent(in)    :: energy(:, :)
  real(dp), intent(in)    :: c_matrix(:, :, :)
  real(dp), intent(in)    :: hfock_ao(:, :, :)
  real(dp), intent(inout) :: p_matrix_mo(:, :, :)
  !=====
  integer :: kstate, cstate, nstate
  integer :: it, jt, nocc, nvirt, nmo, nt, imo, nneg
  real(dp), allocatable :: a11(:, :), b1(:), sigma_inf1(:)
  integer :: istate, jstate, astate, bstate, pstate, qstate
  integer :: pqspin
  real(dp) :: deltae
  real(dp) :: eri_iajb, eri_ibja, eri_ijab
  integer :: info
  integer, allocatable :: ipiv(:)
  real(dp), allocatable :: hfock_mo(:, :, :)
  real(dp), allocatable :: hartree_tmp(:), exchange_tmp(:, :, :)
  real(dp), allocatable :: eigval(:), eigvec(:, :)
  !=====

  write(stdout, '(/,1x,a)') &
      'Renormalization of the occupied-virtual coupling block following Cederbaum Appendix B'
  pqspin=1

  nstate = SIZE(energy, DIM=1)
  nvirt = nvirtual_G-1 - nhomo_G
  nocc  = nhomo_G - ncore_G
  nmo   = nocc + nvirt
  nt = nvirt * nocc
  allocate(b1(nt))
  allocate(sigma_inf1(nt))

  ! First create the Fock hamiltonian in MO basis
  call clean_allocate('Fock matrix F_MO', hfock_mo, nstate, nstate, nspin)
  call h_ao_to_mo(c_matrix, hfock_ao, hfock_mo)

  ! And fill b1 with the occ-virt Fock hamiltonian in MO basis
  ! Note that
  ! < i | Σₓ - vₓ𞁞 | a > = < i | h + Σₓ - (h+vₓ𞁞) | a >
  !                      = < i | F | a > - δᵢₐ εₐ
  !                      = < i | F | a >
  it = 0
  do istate=ncore_G+1, nhomo_G
    do astate=nhomo_G+1, nvirtual_G-1
      it = it + 1
      b1(it) = hfock_mo(istate, astate, pqspin)
    enddo
  enddo

  call clean_deallocate('Fock matrix F_MO', hfock_mo)

  if(has_auxil_basis) then
    call calculate_eri_3center_mo(c_matrix, ncore_G+1, nvirtual_G-1, ncore_G+1, nvirtual_G-1)
  else
    call die("cederbaum_blas: not coded without RI")
  endif

  !
  ! Step 1: Build B= (B₁
  !                   B₂)
  !
  ! Block 1: occ-virt block
  ! Equation B2.b
  !
  call start_clock(timing_tmp2)
  !
  ! B₁ᵢₐ = B₁ᴴᵢₐ + B₁ˣᵢₐ
  !

  !
  ! Hartree term
  !
  ! B₁ᴴᵢₐ = Σ_pq 2 (i a| p q ) γpq
  write(stdout, *) 'Direct term in B1'

  allocate(hartree_tmp(nauxil_local))
  hartree_tmp(:) = 0.0_dp
  do pstate=ncore_G+1, nvirtual_G-1
    do qstate=ncore_G+1, nvirtual_G-1
      hartree_tmp(:) = hartree_tmp(:) + eri_3center_mo(:, pstate, qstate, pqspin) &
                                        * p_matrix_mo(pstate, qstate, pqspin) / spin_fact
    enddo
  enddo
  it = 0
  do istate=ncore_G+1, nhomo_G
    do astate=nhomo_G+1, nvirtual_G-1
      it = it + 1
      b1(it) = b1(it) + spin_fact * DOT_PRODUCT(eri_3center_mo(:, istate, astate, 1), hartree_tmp(:))
    enddo
  enddo
  deallocate(hartree_tmp)

  !
  ! Exchange term
  !
  ! B₁ˣᵢₐ = -Σ_pq (i p| a q ) γpq
  write(stdout, *) 'Exchange term in B1'

  allocate(eigval(nmo), eigvec(nmo, nmo))
  eigvec(:, :) = p_matrix_mo(ncore_G+1:nvirtual_G-1, ncore_G+1:nvirtual_G-1, pqspin) / spin_fact
  ! in-place diago
  call diagonalize(' ', eigvec, eigval)
  ! negative and positive eigenvalues are treated separately because the square-root
  nneg = COUNT(eigval(:) < 0.0d0)

  do imo=1, nmo
    eigvec(:, imo) = eigvec(:, imo) * SQRT(ABS(eigval(imo)))
  enddo

  allocate(exchange_tmp(nauxil_local, nmo, ncore_G+1:nvirtual_G-1))
  do pstate=ncore_G+1, nvirtual_G-1
    !exchange_tmp(:, :, pstate) = MATMUL( eri_3center_mo(:, ncore_G+1:nvirtual_G-1, pstate, 1), eigvec(:, :) )
    call DGEMM('N', 'N', nauxil_local, nmo, nmo, &
               1.0d0, eri_3center_mo(1, ncore_G+1, pstate, 1), nauxil_local, &
               eigvec, nmo, &
               0.0d0, exchange_tmp(1, 1, pstate), nauxil_local)
  enddo
  deallocate(eigvec)

  it = 0
  do istate=ncore_G+1, nhomo_G
    do astate=nhomo_G+1, nvirtual_G-1
      it = it + 1
      b1(it) = b1(it) + SUM(exchange_tmp(:, :nneg  , istate) * exchange_tmp(:, :nneg, astate) )
      b1(it) = b1(it) - SUM(exchange_tmp(:, nneg+1:, istate) * exchange_tmp(:, nneg+1:, astate) )
    enddo
  enddo
  deallocate(exchange_tmp)
  call auxil%sum(b1)

  write(stdout, *) 'B1 is set up'

  ! Naive implementation
  !it = 0
  !do istate=ncore_G+1, nhomo_G
  !  do astate=nhomo_G+1, nvirtual_G-1
  !    it = it + 1
  !    do pstate=ncore_G+1, nvirtual_G-1
  !      do qstate=ncore_G+1, nvirtual_G-1
  !        eri_iqpa = DOT_PRODUCT(eri_3center_mo(:, istate, qstate, 1), eri_3center_mo(:, pstate, astate, 1) )
    !        b1(it) = b1(it) - eri_iqpa * p_matrix_mo(pstate, qstate, 1) / spin_fact
  !      enddo
  !    enddo
  !  enddo
  !enddo

  call stop_clock(timing_tmp2)

  !
  ! Step 2: Build A = ( A₁₁  0 )
  !                   ( A₂₁  0 )
  !
  ! to obtain
  ! A⁻¹ = ( (1-A₁₁)⁻¹       0 )
  !       ( A₂₁·(1-A₁₁)⁻¹   1 )

  ! Step 2.1: do (1 - A₁₁)
  ! Block 11: indices  (ia, jb)
  call start_clock(timing_tmp3)
  call clean_allocate('A11 matrix', a11, nt, nt)

  it = 0
  do istate=ncore_G+1, nhomo_G
    do astate=nhomo_G+1, nvirtual_G-1
      it = it + 1
      jt = 0
      do jstate=ncore_G+1, nhomo_G
        do bstate=nhomo_G+1, nvirtual_G-1
          jt = jt + 1
          deltae = energy(jstate, pqspin) - energy(bstate, pqspin)
          eri_iajb = DOT_PRODUCT(eri_3center_mo(:, istate, astate, pqspin), eri_3center_mo(:, jstate, bstate, pqspin) )
          eri_ibja = DOT_PRODUCT(eri_3center_mo(:, istate, bstate, pqspin), eri_3center_mo(:, jstate, astate, pqspin) )
          eri_ijab = DOT_PRODUCT(eri_3center_mo(:, istate, jstate, pqspin), eri_3center_mo(:, astate, bstate, pqspin) )
          a11(it, jt) = ( 4.0_dp * eri_iajb - eri_ibja - eri_ijab ) / deltae
        enddo
      enddo
    enddo
  enddo
  call stop_clock(timing_tmp3)
  call auxil%sum(a11)
  a11(:, :) = -a11(:, :)
  do it=1, nt
    a11(it, it) = a11(it, it) + 1.0d0
  enddo

  !
  ! a11 now contains (1-A₁₁)
  !
  ! Solve linear system:
  ! (1-A₁₁) · Σ₁(∞) = (1-A₁₁)⁻¹ · B₁
  write(stdout, '(1x,a)') 'Solve (1-A₁₁) · Σ₁(∞) = B₁'
  !sigma_inf1(:) = b1(:)
  call move_alloc(b1, sigma_inf1)
  allocate(ipiv(nt))
  call DGESV(nt, 1, a11, nt, ipiv, sigma_inf1, nt, info)
  deallocate(ipiv)

  if(info /= 0) then
    call die("cederbaum_blas: error in the linear system solution")
  endif
  call clean_deallocate('A11 matrix', a11)

  !
  ! Step 3: Add the density-matrix correction to the GW density-matrix
  ! Note: only the occ-virt block is involved
  it = 0
  do istate=ncore_G+1, nhomo_G
    do astate=nhomo_G+1, nvirtual_G-1
      it = it + 1
      deltae = energy(istate, pqspin) - energy(astate, pqspin)
      p_matrix_mo(istate, astate, pqspin) = p_matrix_mo(istate, astate, pqspin) + sigma_inf1(it) * spin_fact / deltae
      p_matrix_mo(astate, istate, pqspin) = p_matrix_mo(istate, astate, pqspin)
    enddo
  enddo
  deallocate(sigma_inf1)

  call destroy_eri_3center_mo()

end subroutine cederbaum_blas


!=========================================================================
end module m_mbpt_density_matrix
!=========================================================================
