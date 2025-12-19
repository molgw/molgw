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

  integer, parameter, private :: nstage = 3
  
  integer, private :: current_stage
  integer, private :: id_last_used

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
    procedure :: print  => timer_print
  end type

  type(timer) :: timer_molgw
  type(timer) :: timer_prescf
  type(timer) :: timer_scf
  type(timer) :: timer_postscf

  type(timer) :: timer_tmp0
  type(timer) :: timer_tmp1
  type(timer) :: timer_tmp2
  type(timer) :: timer_tmp3
  type(timer) :: timer_tmp4
  type(timer) :: timer_tmp5
  type(timer) :: timer_tmp6
  type(timer) :: timer_tmp7
  type(timer) :: timer_tmp8
  type(timer) :: timer_tmp9






  integer, parameter :: NTIMING=150

  integer, parameter :: timing_total               =  1

  integer, parameter :: timing_relativistic        = 80
  integer, parameter :: timing_prescf              = 81
  integer, parameter :: timing_scf                 = 82
  integer, parameter :: timing_postscf             = 83

  integer, parameter :: timing_dft_xc              =  2
  integer, parameter :: timing_pola                =  3
  integer, parameter :: timing_gw_self             =  4
  integer, parameter :: timing_overlap             =  5
  integer, parameter :: timing_eri_4center         =  6
  integer, parameter :: timing_exchange            =  7
  integer, parameter :: timing_hartree             =  8
  integer, parameter :: timing_sqrt_density_matrix =  9
  integer, parameter :: timing_diago_h2p           = 10
  integer, parameter :: timing_rpa_static          = 11
  integer, parameter :: timing_mp2_energy          = 12
  integer, parameter :: timing_pt_self             = 13
  integer, parameter :: timing_eri_4center_ao2mo   = 14
  integer, parameter :: timing_overlap_grad        = 15
  integer, parameter :: timing_eri_2center         = 16
  integer, parameter :: timing_eri_3center         = 17
  integer, parameter :: timing_eri_3center_ao2mo   = 18
  integer, parameter :: timing_vchiv               = 19
  integer, parameter :: timing_build_h2p           = 20
  integer, parameter :: timing_restart_file        = 21
  integer, parameter :: timing_diago_hamiltonian   = 22
  integer, parameter :: timing_hamiltonian_nuc     = 23
  integer, parameter :: timing_hamiltonian_kin     = 24
  integer, parameter :: timing_build_common        = 25
  integer, parameter :: timing_build_tddft         = 26
  integer, parameter :: timing_build_bse           = 27
  integer, parameter :: timing_spectrum            = 28
  integer, parameter :: timing_eri_screening       = 29
  integer, parameter :: timing_hamiltonian_ecp     = 30
  integer, parameter :: timing_grid_init           = 32
  integer, parameter :: timing_diis                = 33
  integer, parameter :: timing_approx_ham          = 34
  integer, parameter :: timing_fno                 = 36
  integer, parameter :: timing_full_ci             = 37
  integer, parameter :: timing_vertex_selfenergy   = 38
  integer, parameter :: timing_ecp                 = 39
  integer, parameter :: timing_density_matrix      = 40
  integer, parameter :: timing_force               = 41
  integer, parameter :: timing_rpa_dynamic         = 42
  integer, parameter :: timing_ci_selfenergy       = 43
  integer, parameter :: timing_ham_ci              = 44
  integer, parameter :: timing_ci_diago            = 45
  integer, parameter :: timing_ci_write            = 46
  integer, parameter :: timing_ci_config           = 47
  integer, parameter :: timing_zeroes_ci           = 48
  integer, parameter :: timing_density_matrix_cmplx= 49
  integer, parameter :: timing_aomo_pola           = 50
  integer, parameter :: timing_aomo_gw             = 51
  integer, parameter :: timing_mbpt_dm             = 52
  integer, parameter :: timing_eri_3center_ints    = 53
  integer, parameter :: timing_eri_3center_matmul  = 54
  integer, parameter :: timing_dft_densities       = 55
  integer, parameter :: timing_dft_libxc           = 56
  integer, parameter :: timing_dft_vxc             = 57
  integer, parameter :: timing_x_m_vxc             = 58
  integer, parameter :: timing_auto_auxil          = 59
  integer, parameter :: timing_stopping            = 60
  integer, parameter :: timing_noft_energy         = 61
  integer, parameter :: timing_rhoauxil            = 62
  integer, parameter :: timing_eri_2center_ints    = 63
  integer, parameter :: timing_eri_2center_invert  = 64
  integer, parameter :: timing_eri_2center_inverse_sqrt = 65
  integer, parameter :: timing_grid_generation     = 66
  integer, parameter :: timing_grid_wfn            = 67
  integer, parameter :: timing_read_coulombvertex  = 68
  integer, parameter :: timing_write_coulombvertex = 69
  integer, parameter :: timing_density_matrix_MO   = 70

  !
  ! PBC related timers
  !
  integer, parameter :: timing_pbc_eval_bf          = 85
  integer, parameter :: timing_pbc_density          = 86
  integer, parameter :: timing_pbc_potential_to_hao = 87
  integer, parameter :: timing_pbc_nuclei_density   = 88

  !
  ! temporary timers to be used when coding and then freed again
  !
  integer, parameter :: timing_tmp0                = 90
  integer, parameter :: timing_tmp1                = 91
  integer, parameter :: timing_tmp2                = 92
  integer, parameter :: timing_tmp3                = 93
  integer, parameter :: timing_tmp4                = 94
  integer, parameter :: timing_tmp5                = 95
  integer, parameter :: timing_tmp6                = 96
  integer, parameter :: timing_tmp7                = 97
  integer, parameter :: timing_tmp8                = 98
  integer, parameter :: timing_tmp9                = 99

  !
  ! RT-TDDFT timers
  !
  integer, parameter :: timing_tddft_loop             = 110
  integer, parameter :: timing_tddft_fourier          = 111
  integer, parameter :: timing_tddft_one_iter         = 112
  integer, parameter :: timing_tddft_propagation      = 113
  integer, parameter :: timing_tddft_hamiltonian      = 114
  integer, parameter :: timing_tddft_xc               = 115
  integer, parameter :: timing_tddft_exchange         = 116
  integer, parameter :: timing_tddft_hartree          = 117
  integer, parameter :: timing_tddft_hamiltonian_nuc  = 118
  integer, parameter :: timing_tddft_ham_orthobasis   = 119
  integer, parameter :: timing_tddft_eri_2center      = 120
  integer, parameter :: timing_tddft_eri_2center_ints = 121
  integer, parameter :: timing_tddft_eri_2center_invert = 122
  integer, parameter :: timing_update_basis_eri       = 123
  integer, parameter :: timing_update_overlaps        = 124
  integer, parameter :: timing_print_cube_rho_tddft   = 125
  integer, parameter :: timing_restart_tddft_file     = 126
  integer, parameter :: timing_propagate_diago        = 127
  integer, parameter :: timing_propagate_matmul       = 128
  integer, parameter :: timing_print_line_rho_tddft   = 129
  integer, parameter :: timing_calc_dens_disc         = 130
  integer, parameter :: timing_tddft_densities        = 131
  integer, parameter :: timing_tddft_libxc            = 132
  integer, parameter :: timing_tddft_vxc              = 133
  integer, parameter :: timing_tddft_frozen_core      = 134
  integer, parameter :: timing_tddft_q_matrix         = 135
  integer, parameter :: timing_tddft_rhoauxil         = 136
  integer, parameter :: timing_propagate_inverse      = 137
  integer, parameter :: timing_tddft_grid_init        = 139
  integer, parameter :: timing_tddft_grid_generation  = 140
  integer, parameter :: timing_tddft_grid_wfn         = 141
  integer, parameter :: timing_tddft_eri_3center      = 142
  integer, parameter :: timing_tddft_eri_3center_ints = 143
  integer, parameter :: timing_tddft_kin              = 144
  integer, parameter :: timing_grad_kin               = 145
  integer, parameter :: timing_grad_nuc               = 146


  integer, private     :: count_rate, count_max


contains


!=========================================================================
subroutine init_timers()

  call system_clock(COUNT_RATE=count_rate, COUNT_MAX=count_max)

  ! General timers
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


end subroutine init_timers


!=========================================================================
subroutine timers_setstage(stage)
  integer, intent(in) :: stage
  !=====

  current_stage = stage

end subroutine timers_setstage


!=========================================================================
subroutine timer_init(t, label, id)
  class(timer), intent(inout) :: t
  character(len=*), intent(in) :: label
  integer, optional :: id
  !=====

  if( PRESENT(id) ) then
    t%id = id
  else
    id_last_used = id_last_used + 1
    t%id = id_last_used
  endif
  t%running   = .FALSE.
  t%label     = label
  t%timing(:) = 0.0_dp
  t%calls(:)  = 0
  t%has_been_printed(:) = .FALSE.

end subroutine timer_init


!=========================================================================
subroutine timer_start(t)
  class(timer), intent(inout) :: t
  !=====
  integer :: count_tmp
  !=====

  if( t%running ) then
    write(stdout, *) 'Timer is already started:', t%label
    call die('timer_start: timer already started:' // t%label)
  endif

  t%running =.TRUE.

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

  if( .NOT. t%running) then
    write(stdout, *) 'Timer had never been started:', t%label
    call die('timer_stop: timer had never been started:' // t%label)
  endif


  t%running = .FALSE.

#if defined(HAVE_NVTX)
  call nvtxEndRange()
#endif

  call system_clock(COUNT=count_tmp)
  t%timing(current_stage) = t%timing(current_stage) &
       + MODULO( count_tmp - NINT(t%start_time) , count_max) / REAL(count_rate, KIND=dp)

end subroutine timer_stop


!=========================================================================
function timer_get(t, stage) RESULT(time)
  class(timer), intent(inout) :: t
  integer, optional :: stage
  real(dp)          :: time
  !=====
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

  !call output_timing_line('Total Relativistic' , timing_relativistic , 0)
  call timer_prescf%print(0)
  call timer_scf%print(0)
  call timer_postscf%print(0)

  write(stdout, '(/,a,/)') '                 -------------------------------------'
  write(stdout, '(a,/)')   '                             Pre SCF'

  call output_timing_line('Automatic auxiliary basis', timing_auto_auxil, 1)
  call output_timing_line('Integral pre-screening', timing_eri_screening, 1)
  call output_timing_line('4-center integrals', timing_eri_4center, 1)
  call output_timing_line('2-center integrals', timing_eri_2center, 1)
  call output_timing_line('Integrals evaluation', timing_eri_2center_ints, 2)
  call output_timing_line('Matrix inversion', timing_eri_2center_invert, 2)
  call output_timing_line('Matrix inverse sqrt', timing_eri_2center_inverse_sqrt, 2)
  call output_timing_line('3-center integrals', timing_eri_3center, 1)
  call output_timing_line('Integrals evaluation', timing_eri_3center_ints, 2)
  call output_timing_line('Matrix multiplication', timing_eri_3center_matmul, 2)
  call output_timing_line('Overlap matrix S', timing_overlap, 1)
  call output_timing_line('Approximate guess Hamiltonian', timing_approx_ham, 1)
  call output_timing_line('Kinetic Hamiltonian', timing_hamiltonian_kin, 1)
  call output_timing_line('Electron-nucleus Hamiltonian', timing_hamiltonian_nuc, 1)
  call output_timing_line('Effective core potential Hamiltonian', timing_ecp, 1)
  call output_timing_line('ECP Hamiltonian', timing_hamiltonian_ecp, 1)
  call output_timing_line('PBC: basis functions on grid', timing_pbc_eval_bf, 1)
  call output_timing_line('PBC: nuclei density on FFT grid', timing_pbc_density, 1)

  write(stdout, '(/,a,/)') '                 -------------------------------------'
  write(stdout, '(a,/)')   '                                 SCF'

  call output_timing_line('DFT grid initialization', timing_grid_init, 1)
  call output_timing_line('Grid generation', timing_grid_generation, 2)
  call output_timing_line('Wavefunction evaluation', timing_grid_wfn, 2)
  call output_timing_line('Density matrix', timing_density_matrix, 1)
  call output_timing_line('Auxiliary basis density', timing_rhoauxil, 1)
  call output_timing_line('Hartree potential', timing_hartree, 1)
  call output_timing_line('PBC: density on FFT grid', timing_pbc_density, 2)
  call output_timing_line('Exchange operator', timing_exchange, 1)
  call output_timing_line('DFT xc potential', timing_dft_xc, 1)
  call output_timing_line('Densities on a grid', timing_dft_densities, 2)
  call output_timing_line('LIBXC calls', timing_dft_libxc, 2)
  call output_timing_line('Setting up Vxc ', timing_dft_vxc, 2)
  call output_timing_line('PBC: from v(r) to H_AO', timing_pbc_potential_to_hao, 1)
  call output_timing_line('Hamiltonian diagonalization', timing_diago_hamiltonian, 1)
  call output_timing_line('Pulay DIIS mixing', timing_diis, 1)
  call output_timing_line('RESTART file writing', timing_restart_file, 1)
  call output_timing_line('Virtual FNO generation', timing_fno, 1)
  call output_timing_line('Forces', timing_force,1)


  write(stdout, '(/,a,/)') '                 -------------------------------------'
  write(stdout, '(a,/)')   '                            Post SCF'

  ! Prepare post scf
  call output_timing_line('Reading Coulomb vertex file', timing_read_coulombvertex, 1)
  call output_timing_line('Writing Coulomb vertex file', timing_write_coulombvertex, 1)
  call output_timing_line('Sigma_x - Vxc', timing_x_m_vxc, 1)

  ! Linear response polarization RPA or TDDFT or BSE
  call output_timing_line('3-center AO to MO transform', timing_eri_3center_ao2mo, 1)
  call output_timing_line('Response function chi on grid', timing_rpa_dynamic, 1)
  call output_timing_line('Response function chi', timing_pola, 1)
  call output_timing_line('3-center AO to MO transform in chi', timing_aomo_pola, 2)
  call output_timing_line('4-center AO to MO transform', timing_eri_4center_ao2mo, 2)
  call output_timing_line('Static polarization for BSE', timing_rpa_static, 2)
  call output_timing_line('Build 2-particle Hamiltonian', timing_build_h2p, 2)
  call output_timing_line('RPA part', timing_build_common, 3)
  call output_timing_line('TDDFT part', timing_build_tddft, 3)
  call output_timing_line('BSE part', timing_build_bse, 3)
  call output_timing_line('Diago 2 particle H', timing_diago_h2p, 2)
  call output_timing_line('Build W', timing_vchiv,2)
  call output_timing_line('Optical spectrum', timing_spectrum, 2)
  call output_timing_line('Stopping power', timing_stopping, 2)

  ! Self-energies
  call output_timing_line('MBPT density matrix', timing_mbpt_dm, 1)
  call output_timing_line('GW self-energy', timing_gw_self, 1)
  call output_timing_line('3-center AO to MO transform in GW', timing_aomo_gw, 2)
  call output_timing_line('PT self-energy', timing_pt_self, 1)
  call output_timing_line('Vertex corrected self-energy', timing_vertex_selfenergy, 1)
  call output_timing_line('MP2 energy', timing_mp2_energy, 1)

  ! CI
  call output_timing_line('Full CI for few electrons', timing_full_ci, 1)
  call output_timing_line('Setup CI configurations', timing_ci_config, 2)
  call output_timing_line('Screen CI Hamiltonian zeroes', timing_zeroes_ci, 2)
  call output_timing_line('Build CI Hamiltonian', timing_ham_ci, 2)
  call output_timing_line('CI diagonalization', timing_ci_diago, 2)
  call output_timing_line('CI eigenvector file writing', timing_ci_write, 2)
  call output_timing_line('CI self-energy', timing_ci_selfenergy, 2)

  ! NOFT
  call output_timing_line('NOFT calculation', timing_noft_energy, 1)

  ! RT-TDDFT
  call output_timing_line('Real-time TDDFT', timing_tddft_loop, 1)
  call output_timing_line('TDDFT propagator', timing_tddft_propagation, 2)
  call output_timing_line('TDDFT propagator diago', timing_propagate_diago, 3)
  call output_timing_line('TDDFT propagator matmul', timing_propagate_matmul, 3)
  call output_timing_line('TDDFT propagator invert', timing_propagate_inverse, 3)

  call output_timing_line('Update basis, auxilary and eri', timing_update_basis_eri, 2)
  call output_timing_line('Update 2-center ERI', timing_tddft_eri_2center_ints, 3)
  call output_timing_line('Integrals evaluation', timing_tddft_eri_2center_ints, 4)
  call output_timing_line('Matrix inversion', timing_tddft_eri_2center_invert, 4)
  call output_timing_line('Update 3-center ERI', timing_tddft_eri_3center, 3)
  call output_timing_line('Integrals evaluation', timing_tddft_eri_3center_ints, 4)
  call output_timing_line('Update S and D matrices', timing_update_overlaps, 2)
  call output_timing_line('Overlap gradient', timing_overlap_grad, 3)
  call output_timing_line('Grid initialization', timing_tddft_grid_init, 2)
  call output_timing_line('Grid generation', timing_tddft_grid_generation, 3)
  call output_timing_line('Wavefunction evaluation', timing_tddft_grid_wfn, 3)

  call output_timing_line('TDDFT frozen core', timing_tddft_frozen_core, 2)
  call output_timing_line('TDDFT q_matrix', timing_tddft_q_matrix, 2)

  call output_timing_line('Hamiltonian calculation', timing_tddft_hamiltonian, 2)
  call output_timing_line('Kinetic energy', timing_tddft_kin, 3)
  call output_timing_line('Complex density matrix', timing_density_matrix_cmplx, 3)
  call output_timing_line('P in MO basis', timing_density_matrix_MO, 3)
  call output_timing_line('Electron-Nucleus potential', timing_tddft_hamiltonian_nuc, 3)
  call output_timing_line('Auxiliary basis density', timing_tddft_rhoauxil, 3)
  call output_timing_line('Hartree potential', timing_tddft_hartree, 3)
  call output_timing_line('Exchange operator', timing_tddft_exchange, 3)
  call output_timing_line('XC potential', timing_tddft_xc, 3)
  call output_timing_line('Kinetic energy gradient', timing_grad_kin, 3)
  call output_timing_line('Nucleus energy gradient', timing_grad_nuc, 3)
  call output_timing_line('Densities on a grid', timing_tddft_densities, 4)
  call output_timing_line('LIBXC calls', timing_tddft_libxc, 4)
  call output_timing_line('Setting up Vxc ', timing_tddft_vxc, 4)
  call output_timing_line('Orthogonal basis', timing_tddft_ham_orthobasis, 3)

  call output_timing_line('RESTART_TDDFT file writing', timing_restart_tddft_file, 2)
  call output_timing_line('Cube density file writing', timing_print_cube_rho_tddft, 2)
  call output_timing_line('Line density file writing', timing_print_line_rho_tddft, 2)
  call output_timing_line('Electronic density in discs', timing_calc_dens_disc, 2)

  write(stdout, '(/,a,/)') '                 -------------------------------------'



  !
  ! Developer's timings for temporary use only!
  !
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
      call die('timer_print: label string too long. Shorten it or increase the string length. Ask developers')

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

  prepend = TRIM(prepend) // REPEAT('-', max_length-lt-lp)

  write(key, '(i3.3)') max_length+1

  if( PRESENT(stage) ) then
    timing = t%timing(stage)
    calls  = t%calls(stage)
    t%has_been_printed(stage) = .TRUE.
  else
    timing = SUM(t%timing(:))
    calls  = SUM(t%calls(:))
    t%has_been_printed(:) = .TRUE.
  endif

  ! No writing if the timing counter has never been used
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
