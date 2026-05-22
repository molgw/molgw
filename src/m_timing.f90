!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the timings of the code
!
!=========================================================================
#include "molgw.h"
module m_timing
  use m_definitions
  use m_warning, only: die
#if defined(HAVE_NVTX)
  use m_nvtx
#endif

  implicit none

  !
  ! Calculations are divided in 3 "stages"
  ! 1. preSCF
  ! 2. SCF
  ! 3. postSCF
  !
  integer, parameter, private :: nstage = 3

  integer, parameter :: TIMING_PRESCF  = 1
  integer, parameter :: TIMING_SCF     = 2
  integer, parameter :: TIMING_POSTSCF = 3

  integer, private :: current_stage
  integer, private :: id_last_used
  integer, private :: count_rate, count_max

  type timer
    character(len=64) :: label
    integer           :: id
    logical           :: running
    real(dp)          :: start_time
    real(dp)          :: timing(nstage)
    integer           :: calls(nstage)
    logical           :: has_been_printed(nstage)
  contains
    procedure :: init  => timer_init
    procedure :: start => timer_start
    procedure :: stop  => timer_stop
    procedure :: print => timer_print
  end type

  ! Top-level stage timers
  type(timer) :: timer_molgw
  type(timer) :: timer_prescf
  type(timer) :: timer_scf
  type(timer) :: timer_postscf

  ! Temporary timers for development use
  type(timer) :: timer_tmp0, timer_tmp1, timer_tmp2, timer_tmp3, timer_tmp4
  type(timer) :: timer_tmp5, timer_tmp6, timer_tmp7, timer_tmp8, timer_tmp9

  ! Pre-SCF: integrals
  type(timer) :: timer_auto_auxil
  type(timer) :: timer_eri_screening
  type(timer) :: timer_eri_4center
  type(timer) :: timer_eri_2center
  type(timer) :: timer_eri_2center_ints
  type(timer) :: timer_eri_2center_invert
  type(timer) :: timer_eri_2center_inverse_sqrt
  type(timer) :: timer_eri_3center
  type(timer) :: timer_eri_3center_ints
  type(timer) :: timer_eri_3center_matmul
  type(timer) :: timer_eri_3center_ao2mo
  type(timer) :: timer_eri_4center_ao2mo

  ! Pre-SCF: Hamiltonian
  type(timer) :: timer_overlap
  type(timer) :: timer_overlap_grad
  type(timer) :: timer_approx_ham
  type(timer) :: timer_hamiltonian_kin
  type(timer) :: timer_hamiltonian_nuc
  type(timer) :: timer_hamiltonian_ecp
  type(timer) :: timer_ecp
  type(timer) :: timer_sqrt_density_matrix
  type(timer) :: timer_relativistic

  ! SCF
  type(timer) :: timer_grid_init
  type(timer) :: timer_grid_generation
  type(timer) :: timer_grid_wfn
  type(timer) :: timer_density_matrix
  type(timer) :: timer_rhoauxil
  type(timer) :: timer_hartree
  type(timer) :: timer_exchange
  type(timer) :: timer_dft_xc
  type(timer) :: timer_dft_density
  type(timer) :: timer_dft_libxc
  type(timer) :: timer_dft_vxc
  type(timer) :: timer_diago_hamiltonian
  type(timer) :: timer_diis
  type(timer) :: timer_restart_file
  type(timer) :: timer_fno
  type(timer) :: timer_force

  ! PBC
  type(timer) :: timer_pbc_eval_bf
  type(timer) :: timer_pbc_density
  type(timer) :: timer_pbc_potential_to_hao
  type(timer) :: timer_pbc_nuclei_density

  ! Post-SCF: linear response
  type(timer) :: timer_x_m_vxc
  type(timer) :: timer_pola
  type(timer) :: timer_rpa_dynamic
  type(timer) :: timer_rpa_static
  type(timer) :: timer_aomo_pola
  type(timer) :: timer_build_h2p
  type(timer) :: timer_build_common
  type(timer) :: timer_build_tddft
  type(timer) :: timer_build_bse
  type(timer) :: timer_diago_h2p
  type(timer) :: timer_vchiv
  type(timer) :: timer_spectrum
  type(timer) :: timer_stopping
  type(timer) :: timer_read_coulombvertex
  type(timer) :: timer_write_coulombvertex

  ! Post-SCF: self-energies
  type(timer) :: timer_gw_self
  type(timer) :: timer_aomo_gw
  type(timer) :: timer_pt_self
  type(timer) :: timer_vertex_selfenergy
  type(timer) :: timer_mp2_energy
  type(timer) :: timer_mbpt_dm

  ! CI
  type(timer) :: timer_full_ci
  type(timer) :: timer_ci_config
  type(timer) :: timer_zeroes_ci
  type(timer) :: timer_ham_ci
  type(timer) :: timer_ci_diago
  type(timer) :: timer_ci_write
  type(timer) :: timer_ci_selfenergy

  ! NOFT
  type(timer) :: timer_noft_energy

  ! RT-TDDFT
  type(timer) :: timer_tddft_loop
  type(timer) :: timer_tddft_one_iter
  type(timer) :: timer_tddft_fourier
  type(timer) :: timer_tddft_propagation
  type(timer) :: timer_propagate_diago
  type(timer) :: timer_propagate_matmul
  type(timer) :: timer_propagate_inverse
  type(timer) :: timer_update_basis_eri
  type(timer) :: timer_tddft_eri_2center
  type(timer) :: timer_tddft_eri_2center_ints
  type(timer) :: timer_tddft_eri_2center_invert
  type(timer) :: timer_tddft_eri_3center
  type(timer) :: timer_tddft_eri_3center_ints
  type(timer) :: timer_update_overlaps
  type(timer) :: timer_tddft_grid_init
  type(timer) :: timer_tddft_grid_generation
  type(timer) :: timer_tddft_grid_wfn
  type(timer) :: timer_tddft_frozen_core
  type(timer) :: timer_tddft_q_matrix
  type(timer) :: timer_tddft_hamiltonian
  type(timer) :: timer_tddft_kin
  type(timer) :: timer_density_matrix_cmplx
  type(timer) :: timer_density_matrix_MO
  type(timer) :: timer_tddft_hamiltonian_nuc
  type(timer) :: timer_tddft_rhoauxil
  type(timer) :: timer_tddft_hartree
  type(timer) :: timer_tddft_exchange
  type(timer) :: timer_tddft_xc
  type(timer) :: timer_grad_kin
  type(timer) :: timer_grad_nuc
  type(timer) :: timer_tddft_densities
  type(timer) :: timer_tddft_libxc
  type(timer) :: timer_tddft_vxc
  type(timer) :: timer_tddft_ham_orthobasis
  type(timer) :: timer_restart_tddft_file
  type(timer) :: timer_print_cube_rho_tddft
  type(timer) :: timer_print_line_rho_tddft
  type(timer) :: timer_calc_dens_disc


contains


!=========================================================================
subroutine init_timers()

  call system_clock(COUNT_RATE=count_rate, COUNT_MAX=count_max)
  current_stage = TIMING_PRESCF
  id_last_used  = 0

  ! Stage timers
  call timer_molgw%init('Total')
  call timer_prescf%init('Total pre SCF')
  call timer_scf%init('Total SCF')
  call timer_postscf%init('Total post SCF')

  ! Temporary timers
  call timer_tmp0%init('temporary timer 0')
  call timer_tmp1%init('temporary timer 1')
  call timer_tmp2%init('temporary timer 2')
  call timer_tmp3%init('temporary timer 3')
  call timer_tmp4%init('temporary timer 4')
  call timer_tmp5%init('temporary timer 5')
  call timer_tmp6%init('temporary timer 6')
  call timer_tmp7%init('temporary timer 7')
  call timer_tmp8%init('temporary timer 8')
  call timer_tmp9%init('temporary timer 9')

  call timer_relativistic%init('Total Relativistic')

  ! Pre-SCF: integrals
  call timer_auto_auxil%init('Automatic auxiliary basis')
  call timer_eri_screening%init('Integral pre-screening')
  call timer_eri_4center%init('4-center integrals')
  call timer_eri_2center%init('2-center integrals')
  call timer_eri_2center_ints%init('Integrals evaluation')
  call timer_eri_2center_invert%init('Matrix inversion')
  call timer_eri_2center_inverse_sqrt%init('Matrix inverse sqrt')
  call timer_eri_3center%init('3-center integrals')
  call timer_eri_3center_ints%init('Integrals evaluation')
  call timer_eri_3center_matmul%init('Matrix multiplication')
  call timer_eri_3center_ao2mo%init('3-center AO to MO transform')
  call timer_eri_4center_ao2mo%init('4-center AO to MO transform')

  ! Pre-SCF: Hamiltonian
  call timer_overlap%init('Overlap matrix S')
  call timer_overlap_grad%init('Overlap gradient')
  call timer_approx_ham%init('Approximate guess Hamiltonian')
  call timer_hamiltonian_kin%init('Kinetic Hamiltonian')
  call timer_hamiltonian_nuc%init('Electron-nucleus Hamiltonian')
  call timer_hamiltonian_ecp%init('ECP Hamiltonian')
  call timer_ecp%init('Effective core potential Hamiltonian')
  call timer_sqrt_density_matrix%init('Square root density matrix')

  ! SCF
  call timer_grid_init%init('DFT grid initialization')
  call timer_grid_generation%init('Grid generation')
  call timer_grid_wfn%init('Wavefunction evaluation')
  call timer_density_matrix%init('Density matrix')
  call timer_rhoauxil%init('Auxiliary basis density')
  call timer_hartree%init('Hartree potential')
  call timer_exchange%init('Exchange operator')
  call timer_dft_xc%init('DFT xc potential')
  call timer_dft_density%init('Densities on a grid')
  call timer_dft_libxc%init('LIBXC calls')
  call timer_dft_vxc%init('Setting up Vxc')
  call timer_diago_hamiltonian%init('Hamiltonian diagonalization')
  call timer_diis%init('Pulay DIIS mixing')
  call timer_restart_file%init('RESTART file writing')
  call timer_fno%init('Virtual FNO generation')
  call timer_force%init('Forces')

  ! PBC
  call timer_pbc_eval_bf%init('PBC: basis functions on grid')
  call timer_pbc_density%init('PBC: density on FFT grid')
  call timer_pbc_potential_to_hao%init('PBC: from v(r) to H_AO')
  call timer_pbc_nuclei_density%init('PBC: nuclei density on FFT grid')

  ! Post-SCF: linear response
  call timer_x_m_vxc%init('Sigma_x - Vxc')
  call timer_pola%init('Response function chi')
  call timer_rpa_dynamic%init('Response function chi on grid')
  call timer_rpa_static%init('Static polarization for BSE')
  call timer_aomo_pola%init('3-center AO to MO transform in chi')
  call timer_build_h2p%init('Build 2-particle Hamiltonian')
  call timer_build_common%init('RPA part')
  call timer_build_tddft%init('TDDFT part')
  call timer_build_bse%init('BSE part')
  call timer_diago_h2p%init('Diago 2 particle H')
  call timer_vchiv%init('Build W')
  call timer_spectrum%init('Optical spectrum')
  call timer_stopping%init('Stopping power')
  call timer_read_coulombvertex%init('Reading Coulomb vertex file')
  call timer_write_coulombvertex%init('Writing Coulomb vertex file')

  ! Post-SCF: self-energies
  call timer_gw_self%init('GW self-energy')
  call timer_aomo_gw%init('3-center AO to MO transform in GW')
  call timer_pt_self%init('PT self-energy')
  call timer_vertex_selfenergy%init('Vertex corrected self-energy')
  call timer_mp2_energy%init('MP2 energy')
  call timer_mbpt_dm%init('MBPT density matrix')

  ! CI
  call timer_full_ci%init('Full CI for few electrons')
  call timer_ci_config%init('Setup CI configurations')
  call timer_zeroes_ci%init('Screen CI Hamiltonian zeroes')
  call timer_ham_ci%init('Build CI Hamiltonian')
  call timer_ci_diago%init('CI diagonalization')
  call timer_ci_write%init('CI eigenvector file writing')
  call timer_ci_selfenergy%init('CI self-energy')

  ! NOFT
  call timer_noft_energy%init('NOFT calculation')

  ! RT-TDDFT
  call timer_tddft_loop%init('Real-time TDDFT')
  call timer_tddft_one_iter%init('TDDFT one iteration')
  call timer_tddft_fourier%init('TDDFT Fourier')
  call timer_tddft_propagation%init('TDDFT propagator')
  call timer_propagate_diago%init('TDDFT propagator diago')
  call timer_propagate_matmul%init('TDDFT propagator matmul')
  call timer_propagate_inverse%init('TDDFT propagator invert')
  call timer_update_basis_eri%init('Update basis, auxiliary and eri')
  call timer_tddft_eri_2center%init('2-center ERI (TDDFT)')
  call timer_tddft_eri_2center_ints%init('Update 2-center ERI')
  call timer_tddft_eri_2center_invert%init('Matrix inversion')
  call timer_tddft_eri_3center%init('Update 3-center ERI')
  call timer_tddft_eri_3center_ints%init('Integrals evaluation')
  call timer_update_overlaps%init('Update S and D matrices')
  call timer_tddft_grid_init%init('Grid initialization')
  call timer_tddft_grid_generation%init('Grid generation')
  call timer_tddft_grid_wfn%init('Wavefunction evaluation')
  call timer_tddft_frozen_core%init('TDDFT frozen core')
  call timer_tddft_q_matrix%init('TDDFT q_matrix')
  call timer_tddft_hamiltonian%init('Hamiltonian calculation')
  call timer_tddft_kin%init('Kinetic energy')
  call timer_density_matrix_cmplx%init('Complex density matrix')
  call timer_density_matrix_MO%init('P in MO basis')
  call timer_tddft_hamiltonian_nuc%init('Electron-Nucleus potential')
  call timer_tddft_rhoauxil%init('Auxiliary basis density')
  call timer_tddft_hartree%init('Hartree potential')
  call timer_tddft_exchange%init('Exchange operator')
  call timer_tddft_xc%init('XC potential')
  call timer_grad_kin%init('Kinetic energy gradient')
  call timer_grad_nuc%init('Nucleus energy gradient')
  call timer_tddft_densities%init('Densities on a grid')
  call timer_tddft_libxc%init('LIBXC calls')
  call timer_tddft_vxc%init('Setting up Vxc')
  call timer_tddft_ham_orthobasis%init('Orthogonal basis')
  call timer_restart_tddft_file%init('RESTART_TDDFT file writing')
  call timer_print_cube_rho_tddft%init('Cube density file writing')
  call timer_print_line_rho_tddft%init('Line density file writing')
  call timer_calc_dens_disc%init('Electronic density in discs')

end subroutine init_timers


!=========================================================================
subroutine timers_setstage(stage)
  integer, intent(in) :: stage
  current_stage = stage
end subroutine timers_setstage




!=========================================================================
subroutine timer_init(t, label, id)
  class(timer), intent(inout)  :: t
  character(len=*), intent(in) :: label
  integer, optional            :: id
  !=====

  if( PRESENT(id) ) then
    t%id = id
  else
    id_last_used = id_last_used + 1
    t%id = id_last_used
  endif
  t%running              = .FALSE.
  t%label                = label
  t%timing(:)            = 0.0_dp
  t%calls(:)             = 0
  t%has_been_printed(:)  = .FALSE.

end subroutine timer_init


!=========================================================================
subroutine timer_start(t)
  class(timer), intent(inout) :: t
  !=====
  integer :: count_tmp
  !=====

  if( t%running ) then
    write(stdout, *) 'Timer is already started:', t%label
    call die('timer_start: timer already started: ' // TRIM(t%label))
  endif

  t%running = .TRUE.

#if defined(HAVE_NVTX)
  call nvtxStartRange(t%label, t%id)
#endif

  call system_clock(COUNT=count_tmp)
  t%start_time = count_tmp
  t%calls(current_stage) = t%calls(current_stage) + 1

end subroutine timer_start


!=========================================================================
subroutine timer_stop(t)
  class(timer), intent(inout) :: t
  !=====
  integer :: count_tmp
  !=====

  if( .NOT. t%running ) then
    write(stdout, *) 'Timer had never been started:', t%label
    call die('timer_stop: timer had never been started: ' // TRIM(t%label))
  endif

  t%running = .FALSE.

#if defined(HAVE_NVTX)
  call nvtxEndRange()
#endif

  call system_clock(COUNT=count_tmp)
  t%timing(current_stage) = t%timing(current_stage) &
       + MODULO( count_tmp - NINT(t%start_time), count_max ) / REAL(count_rate, KIND=dp)

end subroutine timer_stop


!=========================================================================
function timer_get(t, stage) RESULT(time)
  class(timer), intent(inout) :: t
  integer, optional           :: stage
  real(dp)                    :: time
  !=====

  if( PRESENT(stage) ) then
    time = t%timing(stage)
  else
    time = t%timing(current_stage)
  endif

end function timer_get


!=========================================================================
subroutine output_timers()

  !=====
  !=====

  write(stdout, '(/,a,/)') '                 --- Timings in (s) and # of calls ---'

  call timer_molgw%print(0)

  write(stdout, '(/,a,/)') '                 -------------------------------------'

  call timer_prescf%print(0)
  call timer_scf%print(0)
  call timer_postscf%print(0)

  write(stdout, '(/,a,/)') '                 -------------------------------------'
  write(stdout, '(a,/)')   '                             Pre SCF'

  call timer_auto_auxil%print(1)
  call timer_eri_screening%print(1)
  call timer_eri_4center%print(1)
  call timer_eri_2center%print(1)
  call timer_eri_2center_ints%print(2)
  call timer_eri_2center_invert%print(2)
  call timer_eri_2center_inverse_sqrt%print(2)
  call timer_eri_3center%print(1)
  call timer_eri_3center_ints%print(2)
  call timer_eri_3center_matmul%print(2)
  call timer_overlap%print(1)
  call timer_approx_ham%print(1)
  call timer_hamiltonian_kin%print(1)
  call timer_hamiltonian_nuc%print(1)
  call timer_ecp%print(1)
  call timer_hamiltonian_ecp%print(1)
  call timer_pbc_eval_bf%print(1)
  call timer_pbc_nuclei_density%print(1)

  write(stdout, '(/,a,/)') '                 -------------------------------------'
  write(stdout, '(a,/)')   '                                 SCF'

  call timer_grid_init%print(1)
  call timer_grid_generation%print(2)
  call timer_grid_wfn%print(2)
  call timer_density_matrix%print(1)
  call timer_rhoauxil%print(1)
  call timer_hartree%print(1)
  call timer_pbc_density%print(2)
  call timer_exchange%print(1)
  call timer_dft_xc%print(1)
  call timer_dft_density%print(2)
  call timer_dft_libxc%print(2)
  call timer_dft_vxc%print(2)
  call timer_pbc_potential_to_hao%print(1)
  call timer_diago_hamiltonian%print(1)
  call timer_diis%print(1)
  call timer_restart_file%print(1)
  call timer_fno%print(1)
  call timer_force%print(1)

  write(stdout, '(/,a,/)') '                 -------------------------------------'
  write(stdout, '(a,/)')   '                            Post SCF'

  call timer_read_coulombvertex%print(1)
  call timer_write_coulombvertex%print(1)
  call timer_x_m_vxc%print(1)

  ! Linear response: RPA, TDDFT, BSE
  call timer_eri_3center_ao2mo%print(1)
  call timer_rpa_dynamic%print(1)
  call timer_pola%print(1)
  call timer_aomo_pola%print(2)
  call timer_eri_4center_ao2mo%print(2)
  call timer_rpa_static%print(2)
  call timer_build_h2p%print(2)
  call timer_build_common%print(3)
  call timer_build_tddft%print(3)
  call timer_build_bse%print(3)
  call timer_diago_h2p%print(2)
  call timer_vchiv%print(2)
  call timer_spectrum%print(2)
  call timer_stopping%print(2)

  ! Self-energies
  call timer_mbpt_dm%print(1)
  call timer_gw_self%print(1)
  call timer_aomo_gw%print(2)
  call timer_pt_self%print(1)
  call timer_vertex_selfenergy%print(1)
  call timer_mp2_energy%print(1)

  ! CI
  call timer_full_ci%print(1)
  call timer_ci_config%print(2)
  call timer_zeroes_ci%print(2)
  call timer_ham_ci%print(2)
  call timer_ci_diago%print(2)
  call timer_ci_write%print(2)
  call timer_ci_selfenergy%print(2)

  ! NOFT
  call timer_noft_energy%print(1)

  ! RT-TDDFT
  call timer_tddft_loop%print(1)
  call timer_tddft_propagation%print(2)
  call timer_propagate_diago%print(3)
  call timer_propagate_matmul%print(3)
  call timer_propagate_inverse%print(3)
  call timer_update_basis_eri%print(2)
  call timer_tddft_eri_2center_ints%print(3)
  call timer_tddft_eri_2center_invert%print(4)
  call timer_tddft_eri_3center%print(3)
  call timer_tddft_eri_3center_ints%print(4)
  call timer_update_overlaps%print(2)
  call timer_overlap_grad%print(3)
  call timer_tddft_grid_init%print(2)
  call timer_tddft_grid_generation%print(3)
  call timer_tddft_grid_wfn%print(3)
  call timer_tddft_frozen_core%print(2)
  call timer_tddft_q_matrix%print(2)
  call timer_tddft_hamiltonian%print(2)
  call timer_tddft_kin%print(3)
  call timer_density_matrix_cmplx%print(3)
  call timer_density_matrix_MO%print(3)
  call timer_tddft_hamiltonian_nuc%print(3)
  call timer_tddft_rhoauxil%print(3)
  call timer_tddft_hartree%print(3)
  call timer_tddft_exchange%print(3)
  call timer_tddft_xc%print(3)
  call timer_grad_kin%print(3)
  call timer_grad_nuc%print(3)
  call timer_tddft_densities%print(4)
  call timer_tddft_libxc%print(4)
  call timer_tddft_vxc%print(4)
  call timer_tddft_ham_orthobasis%print(3)
  call timer_restart_tddft_file%print(2)
  call timer_print_cube_rho_tddft%print(2)
  call timer_print_line_rho_tddft%print(2)
  call timer_calc_dens_disc%print(2)

  write(stdout, '(/,a,/)') '                 -------------------------------------'

  ! Developer temporary timers
  call timer_tmp0%print(1)
  call timer_tmp1%print(1)
  call timer_tmp2%print(1)
  call timer_tmp3%print(1)
  call timer_tmp4%print(1)
  call timer_tmp5%print(1)
  call timer_tmp6%print(1)
  call timer_tmp7%print(1)
  call timer_tmp8%print(1)
  call timer_tmp9%print(1)

  write(stdout, '(/,a,/)') '                 -------------------------------------'

end subroutine output_timers


!=========================================================================
subroutine timer_print(t, level, stage)
  class(timer), intent(inout)   :: t
  integer, intent(in)           :: level
  integer, intent(in), optional :: stage
  !=====
  integer, parameter           :: max_length = 46
  character(len=max_length+10) :: prepend
  integer                      :: lt, lp
  character(len=3)             :: key
  real(dp)                     :: timing
  integer                      :: calls
  !=====

  lt = LEN_TRIM(t%label)

  if( lt > max_length ) &
      call die('timer_print: label too long. Shorten it or increase max_length.')

  select case(level)
  case(0)
    prepend = ''
  case(1)
    prepend = '|'
  case(2)
    prepend = '      |'
  case(3)
    prepend = '            |'
  case(4)
    prepend = '                   |'
  end select

  lp = LEN_TRIM(prepend)
  prepend = TRIM(prepend) // REPEAT('-', max_length - lt - lp)

  write(key, '(i3.3)') max_length + 1

  if( PRESENT(stage) ) then
    timing = t%timing(stage)
    calls  = t%calls(stage)
    t%has_been_printed(stage) = .TRUE.
  else
    timing = SUM(t%timing(:))
    calls  = SUM(t%calls(:))
    t%has_been_printed(:) = .TRUE.
  endif

  ! No output if this timer was never used
  if( calls < 1 ) return

  if( level == 0 ) then
    write(stdout, '(1x,a'//key//',4x,f12.2)') TRIM(t%label), timing
  else
    write(stdout, '(1x,a,1x,a,4x,f12.2,2x,i8)') TRIM(prepend), TRIM(t%label), timing, calls
  endif

end subroutine timer_print


!=========================================================================
end module m_timing
!=========================================================================
