!=========================================================================
! This file is part of MOLGW.
! Authors: Fabien Bruneval, Mauricio Rodriguez-Mayorga
!
! This module contains
! the routines to calculate the ACFD total energy
!
!=========================================================================
#include "molgw.h"
module m_acfd
  use m_definitions
  use m_timing
  use m_warning
  use m_memory
  use m_inputparam
  use m_mpi
  use m_scalapack
  use m_cart_to_pure
  use m_basis_set
  use m_dft_grid
  use m_hamiltonian_twobodies
  use m_spectral_function
  use m_spectra
  use m_scf
  use m_selfenergy_tools
  use m_gw_selfenergy_grid
  use m_linear_response
  use m_numerical_tools, only: coeffs_gausslegint

#if !defined(NO_LIBXC)
#include <xc_funcs.h>
#endif

contains


!=========================================================================
subroutine acfd_total_energy(basis, nstate, occupation, energy, c_matrix, en_mbpt)
  type(basis_set), intent(in)  :: basis
  integer, intent(in)          :: nstate
  real(dp), intent(in)         :: occupation(:, :), energy(:, :)
  real(dp), intent(in)         :: c_matrix(:, :, :)
  real(dp), allocatable        :: matrix_tmp(:, :, :)
  type(energy_contributions), intent(inout) :: en_mbpt
  !=====
  type(dft_xc_info), allocatable  :: dft_xc_tmp(:)
  type(spectral_function)        :: wpol
  real(dp)                       :: ebse_singlet, erpa_singlet, erpa_triplet, egw_tmp, erpa_sie_KP, erpa_state
  real(dp)                       :: wlambda(acfd_nlambda), lambda(acfd_nlambda)
  real(dp), allocatable           :: x_matrix(:, :), y_matrix(:, :)
  real(dp), allocatable           :: a_matrix(:, :), b_matrix(:, :)
  integer                        :: nmat, desc_x(NDEL)
  !=====

  !
  ! First the "+" part of RPA+
  !
  erpa_sie_KP = 0.0_dp

  if( TRIM(postscf) == 'RPAP' .OR. TRIM(postscf) == 'RPAP_IM' .OR. TRIM(postscf) == 'RPA+' .OR. TRIM(postscf) == 'RPA+_IM') then
    allocate(dft_xc_tmp(2))
    dft_xc_tmp(:)%id = 0
    dft_xc_tmp(:)%nspin = nspin
    dft_xc_tmp(1)%id = XC_LDA_C_PW      ! HEG
    dft_xc_tmp(2)%id = XC_LDA_C_PW_RPA  ! RPA@HEG
    dft_xc_tmp(1)%coeff = one
    dft_xc_tmp(2)%coeff = -one
    call init_libxc_info(dft_xc_tmp)
    call init_dft_grid(basis, grid_level, dft_xc_tmp(1)%needs_gradient, .TRUE., BATCH_SIZE)
    call clean_allocate('XC operator RPA+', matrix_tmp, basis%nbf, basis%nbf, nspin)
    call dft_exc_vxc_batch(BATCH_SIZE, basis, occupation, c_matrix, matrix_tmp, erpa_sie_KP, dft_xc_in=dft_xc_tmp)
    call destroy_dft_grid()
    call destroy_libxc_info(dft_xc_tmp)
    call clean_deallocate('XC operator RPA+', matrix_tmp)
    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout, '(/,a,f16.10)') ' RPA+ Energy scaled by :', kappa_hybrid
      erpa_sie_KP = kappa_hybrid * erpa_sie_KP
    endif
    write(stdout, '(/,a,f19.10)') ' E(LHEG) - E(LRPA) correlation energy (Ha):', erpa_sie_KP
  endif

  ! Should we add the missing short-range part?
  !if( TRIM(postscf) == 'RPALR' ) then
  !  allocate(dft_xc_tmp(2))
  !  dft_xc_tmp(:)%id = 0
  !  dft_xc_tmp(:)%nspin = nspin
  !  dft_xc_tmp(1)%id = XC_LDA_C_PW      ! HEG
  !  dft_xc_tmp(2)%id = XC_LDA_C_PMGB06  ! RPALR@HEG
  !  dft_xc_tmp(1)%coeff = one
  !  dft_xc_tmp(2)%coeff = -one
  !  dft_xc_tmp(2)%gamma = 1.0_dp / rcut_mbpt
  !  call init_libxc_info(dft_xc_tmp)
  !  call init_dft_grid(basis,grid_level,dft_xc_tmp(1)%needs_gradient,.TRUE.,BATCH_SIZE)
  !  call clean_allocate('XC operator RPA+',matrix_tmp,basis%nbf,basis%nbf,nspin)
  !  call dft_exc_vxc_batch(BATCH_SIZE,basis,occupation,c_matrix,matrix_tmp,erpa_sie_KP,dft_xc_in=dft_xc_tmp)
  !  call destroy_dft_grid()
  !  call destroy_libxc_info(dft_xc_tmp)
  !  call clean_deallocate('XC operator RPA+',matrix_tmp)
  !  if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
  !    write(stdout,'(/,a,f16.10)') ' RPA+ Energy scaled by :',kappa_hybrid
  !    erpa_sie_KP = kappa_hybrid * erpa_sie_KP
  !  endif
  !  write(stdout,'(/,a,f19.10)') ' Ec(HEG) - Ec(LR-RPA) correlation energy (Ha):',erpa_sie_KP
  !endif


  select case(TRIM(postscf))
  case('RPA', 'RPAP', 'RPA+','RPALR')
    call wpol%init(nstate, occupation, 0)
    call polarizability(.FALSE., .FALSE., basis, occupation, energy, c_matrix, &
                        erpa_singlet, egw_tmp, wpol, enforce_spin_multiplicity=1)
    call wpol%destroy()
    en_mbpt%rpa = erpa_singlet
    write(stdout, '(a,2x,f19.10)') ' RPA Energy      (Ha):', en_mbpt%rpa

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout, '(/,a,f16.10)') ' RPA Energy scaled by :', kappa_hybrid
      en_mbpt%rpa = kappa_hybrid * en_mbpt%rpa
    endif

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      en_mbpt%total = en_mbpt%total + en_mbpt%rpa + erpa_sie_KP
    else
      en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus &
                    + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa + erpa_sie_KP
    endif

    if( TRIM(postscf) == 'RPAP' .OR. TRIM(postscf) == 'RPA+' ) then
      write(stdout, '(/,a,f19.10)') ' RPA+ correlation energy (Ha):', en_mbpt%rpa + erpa_sie_KP
      write(stdout, '(a,2x,f19.10)') ' RPA+ Total Energy (Ha):', en_mbpt%total
    else
      write(stdout, '(a,2x,f19.10)') ' RPA Total Energy (Ha):', en_mbpt%total
    endif

  case('RPA_IM', 'RPAP_IM', 'RPA+_IM')
    call wpol%init(nstate, occupation, nomega_chi_imag, grid_type=IMAGINARY_QUAD)
    call polarizability_grid_scalapack(occupation, energy, c_matrix, erpa_state, egw_tmp, wpol)
    call wpol%destroy()
    en_mbpt%rpa = erpa_state

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout, '(/,a,f16.10)') ' RPA Energy scaled by :', kappa_hybrid
      en_mbpt%rpa = kappa_hybrid * en_mbpt%rpa
    endif

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      en_mbpt%total = en_mbpt%total + en_mbpt%rpa + erpa_sie_KP
    else
      en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus &
                    + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa + erpa_sie_KP
    endif

    if( TRIM(postscf) == 'RPAP_IM' .OR. TRIM(postscf) == 'RPA+_IM' ) then
      write(stdout, '(/,a,f19.10)') ' RPA+ correlation energy (Ha):', en_mbpt%rpa + erpa_sie_KP
      write(stdout, '(a,2x,f19.10)') ' RPA+ Total Energy (Ha):', en_mbpt%total
    else
      write(stdout, '(a,2x,f19.10)') ' RPA Total Energy (Ha):', en_mbpt%total
    endif

  case('RPAX', 'RPAX-II')
    call wpol%init(nstate, occupation, 0)
    call polarizability(.FALSE., .FALSE., basis, occupation, energy, c_matrix, &
                        erpa_singlet, egw_tmp, wpol, enforce_spin_multiplicity=1)
    call wpol%destroy()
    en_mbpt%rpa = 0.50_dp * erpa_singlet

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout, '(/,a,f16.10)') ' Singlet RPAx Energy scaled by :', kappa_hybrid
      en_mbpt%rpa=kappa_hybrid*en_mbpt%rpa
    endif
    write(stdout, '(a,2x,f19.10)') ' Singlet RPAx Energy contribution      (Ha):', en_mbpt%rpa

    call wpol%init(nstate, occupation, 0)
    call polarizability(.FALSE., .FALSE., basis, occupation, energy, c_matrix, erpa_triplet, &
                        egw_tmp, wpol, enforce_spin_multiplicity=3)
    call wpol%destroy()
    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout, '(/,a,f16.10)') ' Triplet RPAx Energy scaled by :', kappa_hybrid
      erpa_triplet=kappa_hybrid*erpa_triplet
    endif
    write(stdout, '(a,2x,f19.10)') ' Triplet RPAx Energy contribution      (Ha):', 1.50_dp * erpa_triplet
    en_mbpt%rpa = en_mbpt%rpa + 1.50_dp * erpa_triplet
    write(stdout, '(a,2x,f19.10)') ' RPAx Energy      (Ha):', en_mbpt%rpa

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      en_mbpt%total = en_mbpt%total + en_mbpt%rpa
    else
      en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa
    endif

    write(stdout, *)
    write(stdout, '(a,2x,f19.10)') ' RPAx Total Energy (Ha):', en_mbpt%total
    write(stdout, *)

  case('RPA-I')
    write(stdout, '(/,1x,a,i4)') 'RPA with integration over lambda with acfd_nlambda: ', acfd_nlambda
    call coeffs_gausslegint(0.0_dp, 1.0_dp, lambda, wlambda, acfd_nlambda)

    call wpol%init(nstate, occupation, 0)
    nmat = wpol%npole_reso
    m_x = NUMROC(nmat, block_row, iprow_sd, first_row, nprow_sd)
    n_x = NUMROC(nmat, block_col, ipcol_sd, first_col, npcol_sd)
    call DESCINIT(desc_x, nmat, nmat, block_row, block_col, first_row, first_col, cntxt_sd, MAX(1, m_x), info)

    call clean_allocate('X matrix', x_matrix, m_x, n_x)
    call clean_allocate('Y matrix', y_matrix, m_x, n_x)

    call clean_allocate('A matrix', a_matrix, m_x, n_x)
    call clean_allocate('B matrix', b_matrix, m_x, n_x)
    ! Get A and B
    call polarizability(.FALSE., .FALSE., basis, occupation, energy, c_matrix, erpa_singlet, egw_tmp, wpol, &
                        enforce_spin_multiplicity=1, lambda=1.0_dp, a_matrix=a_matrix, b_matrix=b_matrix)
    call remove_a_energy_diag(energy, wpol, a_matrix)
    call wpol%destroy()

    en_mbpt%rpa = 0.0_dp
    do ilambda=1, acfd_nlambda
      write(stdout, '(/,1x,a,i4,a,i4)') '=== Lambda', ilambda, ' / ', acfd_nlambda
      write(stdout, '(1x,a,f12.8)')   'lambda: ', lambda(ilambda)
      call wpol%init(nstate, occupation, 0)
      call polarizability(.FALSE., .FALSE., basis, occupation, energy, c_matrix, erpa_singlet, egw_tmp, wpol, &
                          enforce_spin_multiplicity=1, lambda=lambda(ilambda), x_matrix=x_matrix, y_matrix=y_matrix)
      call wpol%destroy()
   
      call calculate_ec_acft(desc_x, a_matrix, b_matrix, x_matrix, y_matrix, erpa_singlet)
   
      write(stdout, '(1x,a,f12.8)')   'Energy(lambda): ', erpa_singlet
      en_mbpt%rpa = en_mbpt%rpa + erpa_singlet * wlambda(ilambda)
    enddo
    if(has_auxil_basis) call destroy_eri_3center_mo()

    call clean_deallocate('X matrix', x_matrix)
    call clean_deallocate('Y matrix', y_matrix)
    call clean_deallocate('A matrix', a_matrix)
    call clean_deallocate('B matrix', b_matrix)

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout, '(/,a,f16.10)') ' RPA Energy scaled by :', kappa_hybrid
      en_mbpt%rpa=kappa_hybrid*en_mbpt%rpa
    endif
    write(stdout, '(a,2x,f19.10)') ' RPA Energy      (Ha):', en_mbpt%rpa

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      en_mbpt%total = en_mbpt%total + en_mbpt%rpa
    else
      en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa
    endif
    write(stdout, *)
    write(stdout, '(a,2x,f19.10)') ' RPA Total Energy (Ha):', en_mbpt%total
    write(stdout, *)

  case('RPAX-I')

    write(stdout, '(/,1x,a,i4)') 'RPAx-I with integration over lambda with acfd_nlambda: ', acfd_nlambda
    call coeffs_gausslegint(0.0_dp, 1.0_dp, lambda, wlambda, acfd_nlambda)

    call wpol%init(nstate, occupation, 0)
    nmat = wpol%npole_reso
    m_x = NUMROC(nmat, block_row, iprow_sd, first_row, nprow_sd)
    n_x = NUMROC(nmat, block_col, ipcol_sd, first_col, npcol_sd)
    call DESCINIT(desc_x, nmat, nmat, block_row, block_col, first_row, first_col, cntxt_sd, MAX(1, m_x), info)

    call clean_allocate('X matrix', x_matrix, m_x, n_x)
    call clean_allocate('Y matrix', y_matrix, m_x, n_x)

    call clean_allocate('A matrix', a_matrix, m_x, n_x)
    call clean_allocate('B matrix', b_matrix, m_x, n_x)
    !
    ! Only singlet needed for RPAx-I
    !
    ! Get A and B
    call polarizability(.TRUE., .FALSE., basis, occupation, energy, c_matrix, erpa_singlet, egw_tmp, wpol, &
                        enforce_spin_multiplicity=1, lambda=1.0_dp, a_matrix=a_matrix, b_matrix=b_matrix)
    call remove_a_energy_diag(energy, wpol, a_matrix)
    call wpol%destroy()

    en_mbpt%rpa = 0.0_dp
    do ilambda=1, acfd_nlambda
      write(stdout, '(/,1x,a,i4,a,i4)') '=== Lambda', ilambda, ' / ', acfd_nlambda
      write(stdout, '(1x,a,f12.8)')   'lambda: ', lambda(ilambda)
      call wpol%init(nstate, occupation, 0)
      call polarizability(.FALSE., .FALSE., basis, occupation, energy, c_matrix, erpa_singlet, egw_tmp, wpol, &
                          enforce_spin_multiplicity=1, lambda=lambda(ilambda), x_matrix=x_matrix, y_matrix=y_matrix)
      call wpol%destroy()

      call calculate_ec_acft(desc_x, a_matrix, b_matrix, x_matrix, y_matrix, erpa_singlet)

      write(stdout, '(1x,a,f12.8)')   'Energy(lambda): ', erpa_singlet
      en_mbpt%rpa = en_mbpt%rpa + erpa_singlet * wlambda(ilambda)
    enddo
    call wpol%destroy()

    if(has_auxil_basis) call destroy_eri_3center_mo()

    call clean_deallocate('X matrix', x_matrix)
    call clean_deallocate('Y matrix', y_matrix)
    call clean_deallocate('A matrix', a_matrix)
    call clean_deallocate('B matrix', b_matrix)

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      write(stdout, '(/,a,f16.10)') ' RPAx-I Energy scaled by :', kappa_hybrid
      en_mbpt%rpa=kappa_hybrid*en_mbpt%rpa
    endif
    write(stdout, '(a,2x,f19.10)') ' RPAx-I Energy      (Ha):', en_mbpt%rpa

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using RPA (and RPA versions)
      en_mbpt%total = en_mbpt%total + en_mbpt%rpa
    else
      en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%rpa
    endif
    write(stdout, *)
    write(stdout, '(a,2x,f19.10)') ' RPAx-I Total Energy (Ha):', en_mbpt%total
    write(stdout, *)
  
  case('BSE-I')
    !  P.-F. Loos et al. J. Chem. Phys. Lett. 2020, 11, 9, 3536
    !  Pros and Cons of the Bethe-Salpeter Formalism for Ground-State Energies
    write(stdout, '(/,1x,a,i4)') 'BSE-I with integration over lambda with acfd_nlambda: ', acfd_nlambda
    call coeffs_gausslegint(0.0_dp, 1.0_dp, lambda, wlambda, acfd_nlambda)

    call wpol%init(nstate, occupation, 0)
    nmat = wpol%npole_reso
    m_x = NUMROC(nmat, block_row, iprow_sd, first_row, nprow_sd)
    n_x = NUMROC(nmat, block_col, ipcol_sd, first_col, npcol_sd)
    call DESCINIT(desc_x, nmat, nmat, block_row, block_col, first_row, first_col, cntxt_sd, MAX(1, m_x), info)

    call clean_allocate('X matrix', x_matrix, m_x, n_x)
    call clean_allocate('Y matrix', y_matrix, m_x, n_x)

    call clean_allocate('A matrix', a_matrix, m_x, n_x)
    call clean_allocate('B matrix', b_matrix, m_x, n_x)
    ! Get A and B
    call polarizability(.TRUE., .FALSE., basis, occupation, energy, c_matrix, ebse_singlet, egw_tmp, wpol, &
                        enforce_spin_multiplicity=1, lambda=1.0_dp, a_matrix=a_matrix, b_matrix=b_matrix)
    call remove_a_energy_diag(energy, wpol, a_matrix)
    call wpol%destroy()

    en_mbpt%bse = 0.0_dp
    do ilambda=1, acfd_nlambda
      write(stdout, '(/,1x,a,i4,a,i4)') '=== Lambda', ilambda, ' / ', acfd_nlambda
      write(stdout, '(1x,a,f12.8)')   'lambda: ', lambda(ilambda)
      call wpol%init(nstate, occupation, 0, grid_type=0)
      call polarizability(.FALSE., .FALSE., basis, occupation, energy, c_matrix, ebse_singlet, egw_tmp, wpol, &
                          enforce_spin_multiplicity=1, lambda=lambda(ilambda), x_matrix=x_matrix, y_matrix=y_matrix)
      call wpol%destroy()

      call calculate_ec_acft(desc_x, a_matrix, b_matrix, x_matrix, y_matrix, ebse_singlet)

      write(stdout, '(1x,a,f12.8)')   'Energy(lambda): ', ebse_singlet
      en_mbpt%bse = en_mbpt%bse + ebse_singlet * wlambda(ilambda)
    enddo
    if(has_auxil_basis) call destroy_eri_3center_mo()

    call clean_deallocate('X matrix', x_matrix)
    call clean_deallocate('Y matrix', y_matrix)
    call clean_deallocate('A matrix', a_matrix)
    call clean_deallocate('B matrix', b_matrix)

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using BSE energy
      write(stdout, '(/,a,f16.10)') ' BSE-I Energy scaled by :', kappa_hybrid
      en_mbpt%bse=kappa_hybrid*en_mbpt%bse
    endif
    write(stdout, '(a,2x,f19.10)') ' BSE-I Energy      (Ha):', en_mbpt%bse

    if( abs(kappa_hybrid) > 1.0e-10_dp ) then ! Double-hybrids using BSE energy
      en_mbpt%total = en_mbpt%total + en_mbpt%bse
    else
      en_mbpt%total = en_mbpt%nuc_nuc + en_mbpt%kinetic + en_mbpt%nucleus + en_mbpt%hartree + en_mbpt%exx + en_mbpt%bse
    endif
    write(stdout, *)
    write(stdout, '(a,2x,f19.10)') ' BSE-I Total Energy (Ha):', en_mbpt%total
    write(stdout, *)

  case default
    call die('acfd_total_energy: postscf option not recognized')
  end select

  call print_energy_yaml('mbpt energy', en_mbpt)

end subroutine acfd_total_energy


!=========================================================================
subroutine calculate_ec_acft(desc_x, a_matrix, b_matrix, x_matrix, y_matrix, erpa)
  implicit none

  integer, intent(in)   :: desc_x(NDEL)
  real(dp), intent(in)  :: x_matrix(:, :), y_matrix(:, :), a_matrix(:, :), b_matrix(:, :)
  real(dp), intent(out) :: erpa
  !=====
  integer :: nmat
  real(dp), allocatable :: m_matrix(:, :), tmp_matrix(:, :)
#if defined(HAVE_SCALAPACK)
  integer :: nprow, npcol, myprow, mypcol
#endif
  !=====

  ! Beware that a_matrix and b_matrix are symmetric and in the SCALAPACK case, only the lower part is filled.

  allocate(m_matrix, MOLD=a_matrix)
  allocate(tmp_matrix, MOLD=a_matrix)

#if defined(HAVE_SCALAPACK)
  call BLACS_GRIDINFO(desc_x(CTXT_), nprow, npcol, myprow, mypcol)
  ! if only one proc, then use default coding
  if( nprow * npcol > 1 ) then

    ! all global matrices are square and share the same descriptor desc_x
    nmat = desc_x(M_)

    erpa = -0.5_dp * PDLATRA(nmat, a_matrix, 1, 1, desc_x)

    !
    ! TMP = A * X
    call PDSYMM('L', 'L', nmat, nmat, 1.0d0, a_matrix,1,1,desc_x, &
                                        x_matrix, 1, 1, desc_x, &
                                  0.0d0, tmp_matrix, 1, 1, desc_x)
    !
    ! TMP = (A * X) + B * Y
    call PDSYMM('L', 'L', nmat, nmat, 1.0d0, b_matrix,1,1,desc_x, &
                                        y_matrix, 1, 1, desc_x, &
                                  1.0d0, tmp_matrix, 1, 1, desc_x)
    !
    ! M = X**T * (A * X + B * Y)
    call PDGEMM('T', 'N', nmat, nmat, nmat, 1.0d0, x_matrix,1,1,desc_x, &
                                             tmp_matrix, 1, 1, desc_x, &
                                       0.0d0, m_matrix, 1, 1, desc_x)
    erpa = erpa + 0.5_dp * PDLATRA(nmat, m_matrix, 1, 1, desc_x)

    !
    ! TMP = A * Y
    call PDSYMM('L', 'L', nmat, nmat, 1.0d0, a_matrix,1,1,desc_x, &
                                        y_matrix, 1, 1, desc_x, &
                                  0.0d0, tmp_matrix, 1, 1, desc_x)
    !
    ! TMP = (A * Y) + B * X
    call PDSYMM('L', 'L', nmat, nmat, 1.0d0, b_matrix,1,1,desc_x, &
                                        x_matrix, 1, 1, desc_x, &
                                  1.0d0, tmp_matrix, 1, 1, desc_x)
    !
    ! M = Y**T * (A * Y + B * X)
    call PDGEMM('T', 'N', nmat, nmat, nmat, 1.0d0, y_matrix,1,1,desc_x, &
                                             tmp_matrix, 1, 1, desc_x, &
                                       0.0d0, m_matrix, 1, 1, desc_x)
    erpa = erpa + 0.5_dp * PDLATRA(nmat, m_matrix, 1, 1, desc_x)
  else

#endif
    nmat=SIZE(a_matrix, DIM=1)

    erpa = -0.5_dp * matrix_trace(a_matrix)
    !
    ! TMP = A * X
    call DSYMM('L', 'L', nmat, nmat, 1.0d0, a_matrix, nmat,x_matrix,nmat,  &
               0.0d0, tmp_matrix, nmat)
    !
    ! TMP = (A * X) + B * Y
    call DSYMM('L', 'L', nmat, nmat, 1.0d0, b_matrix, nmat,y_matrix,nmat,  &
               1.0d0, tmp_matrix, nmat)
    ! M = X**T * (A * X + B * Y)
    call DGEMM('T', 'N', nmat, nmat, nmat, 1.0d0, x_matrix, nmat,tmp_matrix,nmat, &
               0.0d0, m_matrix, nmat)
    erpa = erpa + 0.5_dp * matrix_trace(m_matrix)

    !
    ! TMP = A * Y
    call DSYMM('L', 'L', nmat, nmat, 1.0d0, a_matrix, nmat,y_matrix,nmat,  &
               0.0d0, tmp_matrix, nmat)
    !
    ! TMP = (A * Y) + B * X
    call DSYMM('L', 'L', nmat, nmat, 1.0d0, b_matrix, nmat,x_matrix,nmat,  &
               1.0d0, tmp_matrix, nmat)
    ! M = Y**T * (A * Y + B * X)
    call DGEMM('T', 'N', nmat, nmat, nmat, 1.0d0, y_matrix, nmat,tmp_matrix,nmat, &
               0.0d0, m_matrix, nmat)
    erpa = erpa + 0.5_dp * matrix_trace(m_matrix)
#if defined(HAVE_SCALAPACK)
  endif
#endif

  !m_matrix(:,:) = -a_matrix(:,:)
  !m_matrix(:,:) = m_matrix(:,:) + MATMUL( TRANSPOSE(x_matrix) , MATMUL(a_matrix,x_matrix) )
  !m_matrix(:,:) = m_matrix(:,:) + MATMUL( TRANSPOSE(y_matrix) , MATMUL(a_matrix,y_matrix) )
  !m_matrix(:,:) = m_matrix(:,:) + MATMUL( TRANSPOSE(x_matrix) , MATMUL(b_matrix,y_matrix) )
  !m_matrix(:,:) = m_matrix(:,:) + MATMUL( TRANSPOSE(y_matrix) , MATMUL(b_matrix,x_matrix) )
  !
  !erpa = 0.0_dp
  !do imat=1,nmat
  !  erpa = erpa + 0.5_dp * matrix_trace(m_matrix)
  !enddo
  deallocate(m_matrix)
  deallocate(tmp_matrix)


end subroutine calculate_ec_acft


end module m_acfd
!=========================================================================
