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
  use m_warning,only: die

  integer,parameter :: NTIMING=150

  integer,parameter :: timing_total               =  1

  integer,parameter :: timing_prescf              = 81
  integer,parameter :: timing_scf                 = 82
  integer,parameter :: timing_postscf             = 83

  integer,parameter :: timing_dft_xc              =  2
  integer,parameter :: timing_pola                =  3
  integer,parameter :: timing_gw_self             =  4
  integer,parameter :: timing_overlap             =  5
  integer,parameter :: timing_eri_4center         =  6
  integer,parameter :: timing_exchange            =  7
  integer,parameter :: timing_hartree             =  8
  integer,parameter :: timing_sqrt_density_matrix =  9
  integer,parameter :: timing_diago_h2p           = 10
  integer,parameter :: timing_rpa_static          = 11
  integer,parameter :: timing_mp2_energy          = 12
  integer,parameter :: timing_pt_self             = 13
  integer,parameter :: timing_eri_4center_eigen   = 14
  integer,parameter :: timing_overlap_grad        = 15
  integer,parameter :: timing_eri_2center         = 16
  integer,parameter :: timing_eri_3center         = 17
  integer,parameter :: timing_eri_3center_eigen   = 18
  integer,parameter :: timing_vchiv               = 19
  integer,parameter :: timing_build_h2p           = 20
  integer,parameter :: timing_restart_file        = 21
  integer,parameter :: timing_diago_hamiltonian   = 22
  integer,parameter :: timing_hamiltonian_nuc     = 23
  integer,parameter :: timing_hamiltonian_kin     = 24
  integer,parameter :: timing_build_common        = 25
  integer,parameter :: timing_build_tddft         = 26
  integer,parameter :: timing_build_bse           = 27
  integer,parameter :: timing_spectrum            = 28
  integer,parameter :: timing_eri_screening       = 29
  integer,parameter :: timing_hamiltonian_ecp     = 30
  integer,parameter :: timing_sca_distr1          = 31
  integer,parameter :: timing_grid_init           = 32
  integer,parameter :: timing_diis                = 33
  integer,parameter :: timing_approx_ham          = 34
  integer,parameter :: timing_sca_distr2          = 35
  integer,parameter :: timing_fno                 = 36
  integer,parameter :: timing_full_ci             = 37
  integer,parameter :: timing_gwgamma_self        = 38
  integer,parameter :: timing_ecp                 = 39
  integer,parameter :: timing_density_matrix      = 40
  integer,parameter :: timing_force               = 41
  integer,parameter :: timing_rpa_dynamic         = 42
  integer,parameter :: timing_ci_selfenergy       = 43
  integer,parameter :: timing_ham_ci              = 44
  integer,parameter :: timing_ci_diago            = 45
  integer,parameter :: timing_ci_write            = 46
  integer,parameter :: timing_ci_config           = 47
  integer,parameter :: timing_zeroes_ci           = 48
  integer,parameter :: timing_density_matrix_cmplx= 49
  integer,parameter :: timing_aomo_pola           = 50
  integer,parameter :: timing_aomo_gw             = 51
  integer,parameter :: timing_mbpt_dm             = 52
  integer,parameter :: timing_eri_3center_ints    = 53
  integer,parameter :: timing_eri_3center_matmul  = 54
  integer,parameter :: timing_dft_densities       = 55
  integer,parameter :: timing_dft_libxc           = 56
  integer,parameter :: timing_dft_vxc             = 57
  integer,parameter :: timing_x_m_vxc             = 58
  integer,parameter :: timing_auto_auxil          = 59
  integer,parameter :: timing_stopping            = 60
  integer,parameter :: timing_noft_energy         = 61
  integer,parameter :: timing_rhoauxil            = 62
  integer,parameter :: timing_eri_2center_ints    = 63
  integer,parameter :: timing_eri_2center_invert  = 64
  integer,parameter :: timing_eri_2center_inverse_sqrt = 65
  integer,parameter :: timing_grid_generation     = 66
  integer,parameter :: timing_grid_wfn            = 67

  integer,parameter :: timing_density_matrix_MO_cmplx = 70

  integer,parameter :: timing_tmp0                = 90
  integer,parameter :: timing_tmp1                = 91
  integer,parameter :: timing_tmp2                = 92
  integer,parameter :: timing_tmp3                = 93
  integer,parameter :: timing_tmp4                = 94
  integer,parameter :: timing_tmp5                = 95
  integer,parameter :: timing_tmp6                = 96
  integer,parameter :: timing_tmp7                = 97
  integer,parameter :: timing_tmp8                = 98
  integer,parameter :: timing_tmp9                = 99

  integer,parameter :: timing_tddft_loop             = 110
  integer,parameter :: timing_tddft_fourier          = 111
  integer,parameter :: timing_tddft_one_iter         = 112
  integer,parameter :: timing_tddft_propagation      = 113
  integer,parameter :: timing_tddft_hamiltonian      = 114
  integer,parameter :: timing_tddft_xc               = 115
  integer,parameter :: timing_tddft_exchange         = 116
  integer,parameter :: timing_tddft_hartree          = 117
  integer,parameter :: timing_tddft_hamiltonian_nuc  = 118
  integer,parameter :: timing_tddft_ham_orthobasis   = 119
  integer,parameter :: timing_tddft_eri_2center      = 120
  integer,parameter :: timing_tddft_eri_2center_ints = 121
  integer,parameter :: timing_tddft_eri_2center_invert = 122
  integer,parameter :: timing_update_basis_eri       = 123
  integer,parameter :: timing_update_overlaps        = 124
  integer,parameter :: timing_print_cube_rho_tddft   = 125
  integer,parameter :: timing_restart_tddft_file     = 126
  integer,parameter :: timing_propagate_diago        = 127
  integer,parameter :: timing_propagate_matmul       = 128
  integer,parameter :: timing_print_line_rho_tddft   = 129
  integer,parameter :: timing_calc_dens_disc         = 130
  integer,parameter :: timing_tddft_densities        = 131
  integer,parameter :: timing_tddft_libxc            = 132
  integer,parameter :: timing_tddft_vxc              = 133
  integer,parameter :: timing_tddft_frozen_core      = 134
  integer,parameter :: timing_tddft_q_matrix         = 135
  integer,parameter :: timing_tddft_rhoauxil         = 136
  integer,parameter :: timing_propagate_inverse      = 137
  integer,parameter :: timing_tddft_grid_init        = 139
  integer,parameter :: timing_tddft_grid_generation  = 140
  integer,parameter :: timing_tddft_grid_wfn         = 141

  integer,parameter :: timing_tddft_eri_3center      = 142
  integer,parameter :: timing_tddft_eri_3center_ints = 143
  integer,parameter :: timing_tddft_kin              = 144


  integer,private     :: count_rate,count_max
  logical,private     :: time_running(NTIMING)
  real(dp),private    :: time_start(NTIMING)
  real(dp),private    :: timing(NTIMING)
  integer(dp),private :: calls(NTIMING)

  logical,protected   :: in_rt_tddft = .FALSE.


contains


!=========================================================================
subroutine init_timing()
  implicit none
  !=====
  !=====

  time_running(:) = .FALSE.
  timing(:)       = 0.0_dp
  calls(:)        = 0

  call system_clock(COUNT_RATE=count_rate,COUNT_MAX=count_max)

end subroutine init_timing


!=========================================================================
function get_timing(itiming) RESULT(time)
  implicit none

  integer,intent(in) :: itiming
  real(dp)           :: time
  !=====
  !=====

  time = timing(itiming)

end function get_timing


!=========================================================================
subroutine start_clock(itiming)
  implicit none
  integer,intent(in) :: itiming
  !=====
  integer            :: count_tmp
  character(len=5)   :: msg
  !=====

  ! 0 means no timing
  if(itiming == 0) return

  if(time_running(itiming)) then
    write(stdout,*) 'clock # is already started:',itiming
    write(msg,'(i05)') itiming
    call die('start_clock: clock already started:' // msg)
  endif

  time_running(itiming)=.TRUE.

  call system_clock(COUNT=count_tmp)
  time_start(itiming) = count_tmp
  calls(itiming) = calls(itiming) + 1

end subroutine start_clock


!=========================================================================
subroutine stop_clock(itiming)
  implicit none

  integer,intent(in) :: itiming
  !=====
  integer            :: count_tmp
  character(len=5)   :: msg
  !=====

  ! 0 means no timing
  if(itiming == 0) return

  if(.NOT.time_running(itiming)) then
    write(stdout,*) 'clock # has not been started:',itiming
    write(msg,'(i05)') itiming
    call die('stop_clock: clock was never started' // msg)
  endif


  time_running(itiming)=.FALSE.

  call system_clock(COUNT=count_tmp)
  timing(itiming) = timing(itiming) + MODULO( count_tmp - NINT(time_start(itiming)) , count_max) / REAL(count_rate,dp)

end subroutine stop_clock


!=========================================================================
subroutine output_timing()

  implicit none
  !=====
  !=====

  write(stdout,'(/,a,/)') '                 --- Timings in (s) and # of calls ---'

  call output_timing_line('Total time' ,timing_total ,0)

  write(stdout,'(/,a,/)') '                 -------------------------------------'

  call output_timing_line('Total pre SCF' ,timing_prescf ,0)
  call output_timing_line('Total SCF'     ,timing_scf    ,0)
  call output_timing_line('Total post SCF',timing_postscf,0)

  write(stdout,'(/,a,/)') '                 -------------------------------------'
  write(stdout,'(a,/)')   '                             Pre SCF'

  call output_timing_line('Automatic auxiliary basis',timing_auto_auxil,1)
  call output_timing_line('Integral pre-screening',timing_eri_screening,1)
  call output_timing_line('4-center integrals',timing_eri_4center,1)
  call output_timing_line('2-center integrals',timing_eri_2center,1)
  call output_timing_line('Integrals evaluation',timing_eri_2center_ints,2)
  call output_timing_line('Matrix inversion',timing_eri_2center_invert,2)
  call output_timing_line('Matrix inverse sqrt',timing_eri_2center_inverse_sqrt,2)
  call output_timing_line('3-center integrals',timing_eri_3center,1)
  call output_timing_line('Integrals evaluation',timing_eri_3center_ints,2)
  call output_timing_line('Matrix multiplication',timing_eri_3center_matmul,2)
  call output_timing_line('Overlap matrix S',timing_overlap,1)
  call output_timing_line('Approximate guess Hamiltonian',timing_approx_ham,1)
  call output_timing_line('Kinetic Hamiltonian',timing_hamiltonian_kin,1)
  call output_timing_line('Electron-nucleus Hamiltonian',timing_hamiltonian_nuc,1)
  call output_timing_line('ECP Hamiltonian',timing_hamiltonian_ecp,1)

  write(stdout,'(/,a,/)') '                 -------------------------------------'
  write(stdout,'(a,/)')   '                                 SCF'

  call output_timing_line('DFT grid initialization',timing_grid_init,1)
  call output_timing_line('Grid generation',timing_grid_generation,2)
  call output_timing_line('Wavefunction evaluation',timing_grid_wfn,2)
  call output_timing_line('Density matrix',timing_density_matrix,1)
  call output_timing_line('Auxiliary basis density',timing_rhoauxil,1)
  call output_timing_line('Hartree potential',timing_hartree,1)
  call output_timing_line('Exchange operator',timing_exchange,1)
  call output_timing_line('DFT xc potential',timing_dft_xc,1)
  call output_timing_line('Densities on a grid',timing_dft_densities,2)
  call output_timing_line('LIBXC calls',timing_dft_libxc,2)
  call output_timing_line('Setting up Vxc ',timing_dft_vxc,2)
  call output_timing_line('Hamiltonian diagonalization',timing_diago_hamiltonian,1)
  call output_timing_line('Pulay DIIS mixing',timing_diis,1)
  call output_timing_line('RESTART file writing',timing_restart_file,1)
  call output_timing_line('Virtual FNO generation',timing_fno,1)
  call output_timing_line('Forces',timing_force,1)


  write(stdout,'(/,a,/)') '                 -------------------------------------'
  write(stdout,'(a,/)')   '                            Post SCF'

  ! Prepare post scf
  call output_timing_line('Sigma_x - Vxc',timing_x_m_vxc,1)

  ! Linear response polarization RPA or TDDFT or BSE
  call output_timing_line('3-center AO to MO transform',timing_eri_3center_eigen,1)
  call output_timing_line('Response function chi on grid',timing_rpa_dynamic,1)
  call output_timing_line('Response function chi',timing_pola,1)
  call output_timing_line('3-center AO to MO transform in chi',timing_aomo_pola,2)
  call output_timing_line('4-center AO to MO transform',timing_eri_4center_eigen,2)
  call output_timing_line('Static polarization for BSE',timing_rpa_static,2)
  call output_timing_line('Build 2-particle Hamiltonian',timing_build_h2p,2)
  call output_timing_line('RPA part',timing_build_common,3)
  call output_timing_line('TDDFT part',timing_build_tddft,3)
  call output_timing_line('BSE part',timing_build_bse,3)
  call output_timing_line('Diago 2 particle H',timing_diago_h2p,2)
  call output_timing_line('Build W',timing_vchiv,2)
  call output_timing_line('Optical spectrum',timing_spectrum,2)
  call output_timing_line('Stopping power',timing_stopping,2)

  ! Self-energies
  call output_timing_line('MBPT density matrix',timing_mbpt_dm,1)
  call output_timing_line('GW self-energy',timing_gw_self,1)
  call output_timing_line('3-center AO to MO transform in GW',timing_aomo_gw,2)
  call output_timing_line('PT self-energy',timing_pt_self,1)
  call output_timing_line('GWGamma self-energy',timing_gwgamma_self,1)
  call output_timing_line('MP2 energy',timing_mp2_energy,1)

  ! CI
  call output_timing_line('Full CI for few electrons',timing_full_ci,1)
  call output_timing_line('Setup CI configurations',timing_ci_config,2)
  call output_timing_line('Screen CI Hamiltonian zeroes',timing_zeroes_ci,2)
  call output_timing_line('Build CI Hamiltonian',timing_ham_ci,2)
  call output_timing_line('CI diagonalization',timing_ci_diago,2)
  call output_timing_line('CI eigenvector file writing',timing_ci_write,2)
  call output_timing_line('CI self-energy',timing_ci_selfenergy,2)

  ! NOFT
  call output_timing_line('NOFT calculation',timing_noft_energy,1)

  ! RT-TDDFT
  call output_timing_line('TDDFT loop',timing_tddft_loop,1)
  call output_timing_line('TDDFT propagator',timing_tddft_propagation,2)
  call output_timing_line('TDDFT propagator diago',timing_propagate_diago,3)
  call output_timing_line('TDDFT propagator matmul',timing_propagate_matmul,3)
  call output_timing_line('TDDFT propagator invert',timing_propagate_inverse,3)

  call output_timing_line('Update basis, auxilary and eri',timing_update_basis_eri,2)
  call output_timing_line('Update 2-center ERI',timing_tddft_eri_2center_ints,3)
  call output_timing_line('Integrals evaluation',timing_tddft_eri_2center_ints,4)
  call output_timing_line('Matrix inversion',timing_tddft_eri_2center_invert,4)
  call output_timing_line('Update 3-center ERI',timing_tddft_eri_3center,3)
  call output_timing_line('Integrals evaluation',timing_tddft_eri_3center_ints,4)
  call output_timing_line('Update S and D matrices',timing_update_overlaps,2)
  call output_timing_line('Gradient of overlap',timing_overlap_grad,3)
  call output_timing_line('Grid initialization',timing_tddft_grid_init,2)
  call output_timing_line('Grid generation',timing_tddft_grid_generation,3)
  call output_timing_line('Wavefunction evaluation',timing_tddft_grid_wfn,3)

  call output_timing_line('TDDFT frozen core',timing_tddft_frozen_core,2)
  call output_timing_line('TDDFT q_matrix',timing_tddft_q_matrix,2)

  call output_timing_line('Hamiltonian calculation',timing_tddft_hamiltonian,2)
  call output_timing_line('Kinetic energy',timing_tddft_kin,3)
  call output_timing_line('Complex density matrix',timing_density_matrix_cmplx,3)
  call output_timing_line('Complex P in MO basis',timing_density_matrix_MO_cmplx,3)
  call output_timing_line('Electron-Nucleus potential',timing_tddft_hamiltonian_nuc,3)
  call output_timing_line('Auxiliary basis density',timing_tddft_rhoauxil,3)
  call output_timing_line('Hartree potential',timing_tddft_hartree,3)
  call output_timing_line('Exchange operator',timing_tddft_exchange,3)
  call output_timing_line('XC potential',timing_tddft_xc,3)
  call output_timing_line('Densities on a grid',timing_tddft_densities,4)
  call output_timing_line('LIBXC calls',timing_tddft_libxc,4)
  call output_timing_line('Setting up Vxc ',timing_tddft_vxc,4)
  call output_timing_line('Orthogonal basis',timing_tddft_ham_orthobasis,3)

  call output_timing_line('RESTART_TDDFT file writing',timing_restart_tddft_file,2)
  call output_timing_line('Cube density file writing',timing_print_cube_rho_tddft,2)
  call output_timing_line('Line density file writing',timing_print_line_rho_tddft,2)
  call output_timing_line('Electronic density in discs',timing_calc_dens_disc,2)

  write(stdout,'(/,a,/)') '                 -------------------------------------'



  !
  ! Developer's timings for temporary use only!
  !
  if( ANY( calls(timing_tmp0:timing_tmp9) > 0 ) .OR. calls(timing_sca_distr1) > 0 .OR. calls(timing_sca_distr2) > 0 ) then
    write(stdout,'(a,/)')   '                            Testing'
    call output_timing_line('timing SCALAPACK tmp1',timing_sca_distr1,1)
    call output_timing_line('timing SCALAPACK tmp2',timing_sca_distr2,1)
    call output_timing_line('Tmp timing 0',timing_tmp0,1)
    call output_timing_line('Tmp timing 1',timing_tmp1,1)
    call output_timing_line('Tmp timing 2',timing_tmp2,1)
    call output_timing_line('Tmp timing 3',timing_tmp3,1)
    call output_timing_line('Tmp timing 4',timing_tmp4,1)
    call output_timing_line('Tmp timing 5',timing_tmp5,1)
    call output_timing_line('Tmp timing 6',timing_tmp6,1)
    call output_timing_line('Tmp timing 7',timing_tmp7,1)
    call output_timing_line('Tmp timing 8',timing_tmp8,1)
    call output_timing_line('Tmp timing 9',timing_tmp9,1)
    write(stdout,'(/,a,/)') '                 -------------------------------------'
  endif


end subroutine output_timing


!=========================================================================
subroutine output_timing_line(title,itiming,level)
  implicit none

  character(len=*),intent(in) :: title
  integer,intent(in)          :: itiming
  integer,intent(in)          :: level
  !=====
  integer,parameter            :: max_length = 46
  character(len=max_length+10) :: prepend
  integer                      :: lt,lp
  character(len=3)             :: key
  !=====

  ! No writing if the timing counter has never been used
  if( calls(itiming) < 1 ) return

  lt = LEN_TRIM(title)

  if( lt > max_length ) &
      call die('output_timing_line: title string too long. Shorten it or increase the string length. Ask developers')

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

  prepend = TRIM(prepend) // REPEAT('-',max_length-lt-lp)

  write(key,'(i3.3)') max_length+1

  if( level == 0 ) then
    write(stdout,'(1x,a'//key//',4x,f12.2)') TRIM(title),timing(itiming)
  else
    write(stdout,'(1x,a,1x,a,4x,f12.2,2x,i8)') TRIM(prepend),TRIM(title),timing(itiming),calls(itiming)
  endif


end subroutine output_timing_line

!=========================================================================
subroutine switch_on_rt_tddft_timers()
  implicit none

  in_rt_tddft = .TRUE.

end subroutine switch_on_rt_tddft_timers


!=========================================================================
subroutine switch_off_rt_tddft_timers()
  implicit none

  in_rt_tddft = .FALSE.

end subroutine switch_off_rt_tddft_timers


!=========================================================================
end module m_timing
!=========================================================================
