!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to set up and store the input parameters from the input file
!
!=========================================================================
#include "molgw.h"
#if !defined(NO_LIBINT)
#include<libint2/libint2_params.h>
#endif
module m_inputparam
  use m_definitions
  use m_mpi
  use m_warning
  use m_atoms
  use m_elements
  use m_ecp
  use m_string_tools,only: capitalize,yesno_to_logical,yesno_to_TrueFalse
  use m_libxc_tools

#if !defined(NO_LIBXC)
#include <xc_funcs.h>
#endif

  !
  ! Self-energy evaluation technique
  integer,parameter :: one_shot                = 101
  integer,parameter :: QS                      = 102
  integer,parameter :: EVSC                    = 103
  integer,parameter :: imaginary_axis_pade     = 104
  integer,parameter :: static_selfenergy       = 105
  integer,parameter :: imaginary_axis_integral = 106
  integer,parameter :: exact_dyson             = 107
  integer,parameter :: imaginary_axis_homolumo = 108
  integer,parameter :: contour_deformation     = 109

  !
  ! Self-energy approximation
  integer,parameter :: CI              =-201
  integer,parameter :: COHSEX          = 204
  integer,parameter :: GnW0            = 205
  integer,parameter :: GnWn            = 206
  integer,parameter :: GW              = 207
  integer,parameter :: GWSOSEX         = 217
  integer,parameter :: GWSOX           = 219
  integer,parameter :: PT2             = 220
  integer,parameter :: ONE_RING        = 221
  integer,parameter :: SOX             = 222
  integer,parameter :: PT3             = 223
  integer,parameter :: TWO_RINGS       = 224
  integer,parameter :: GWPT3           = 226
  integer,parameter :: GWGWG           = 229
  integer,parameter :: GW0GW0G         = 231
  integer,parameter :: GWGW0G          = 232
  integer,parameter :: GWGWG_NUMERICAL = 233
  integer,parameter :: GWTILDE         = 234
  integer,parameter :: GWGW0RPAG       = 235

  !
  ! TDDFT variables
  integer,parameter :: EXCIT_NO                 = 501
  integer,parameter :: EXCIT_LIGHT              = 502
  integer,parameter :: EXCIT_PROJECTILE         = 503
  integer,parameter :: EXCIT_PROJECTILE_W_BASIS = 504

  integer,protected           :: unit_yaml
  character(len=10),parameter :: filename_yaml = 'molgw.yaml'


  type calculation_type
    character(len=100) :: scf_name
    character(len=100) :: postscf_name
    logical            :: is_core
    logical            :: is_dft
    logical            :: is_real_time
    logical            :: need_exchange
    logical            :: need_exchange_lr
    logical            :: need_rpa
    logical            :: is_lr_mbpt
    logical            :: is_gw
    logical            :: is_mp2
    logical            :: is_mp3
    logical            :: is_noft
    logical            :: is_selfenergy
    logical            :: is_ci
    logical            :: is_bse,no_bse_kernel,include_tddft_kernel,include_tdhf_kernel
    integer            :: selfenergy_technique      ! perturbative or quasiparticle self-consistent or eigenvalue-sc
    integer            :: selfenergy_approx         ! GW, COHSEX, PT2
  end type calculation_type

  type excitation_type
    character(len=256)   :: name
    integer              :: form
    real(dp)             :: kappa, omega, time0
    real(dp)             :: dir(3)
  end type


  ! There are the input variables of MOLGW
  ! They should not be modified anywhere else in the code.
  ! Declare them as protected and work on copies if absolutely necessary.
  type(calculation_type),protected :: calc_type
  type(dft_xc_info),allocatable    :: dft_xc(:)
  type(excitation_type),protected  :: excit_type
  logical,protected                :: is_frozencore
  logical,protected                :: is_tddft_frozencore
  logical,protected                :: is_tda,is_triplet
  logical,protected                :: is_virtual_fno
  real(dp),protected               :: spin_fact
  real(dp),protected               :: electrons

  character(len=100),allocatable,protected :: basis_name(:)
  character(len=100),allocatable,protected :: auxil_basis_name(:)
  character(len=100),allocatable,protected :: small_basis_name(:)
  character(len=100),allocatable,protected :: ecp_basis_name(:)
  character(len=100),allocatable,protected :: ecp_auxil_basis_name(:)
  character(len=100),allocatable,protected :: ecp_small_basis_name(:)

  integer,protected                :: tddft_grid_level
  integer,protected                :: grid_level
  integer,protected                :: ecp_level
  integer,protected                :: integral_level
  logical,protected                :: has_auxil_basis
  logical,protected                :: has_small_basis
  logical,protected                :: incore_
  !
  ! the boring small complex number eta: (0.0_dp,0.001_dp) is typically over converged
  ! Having a larger ieta value smoothen the oscillation far from the HOMO-LUMO gap
  complex(dp),protected            :: ieta

  logical,protected                :: use_correlated_density_matrix_
  logical,protected                :: gwgamma_tddft_
  logical,protected                :: memory_evaluation_
  logical,protected                :: read_restart_
  logical,protected                :: ignore_bigrestart_
  logical,protected                :: force_energy_qp_
  logical,protected                :: print_eri_
  logical,protected                :: print_wfn_
  logical,protected                :: print_w_
  logical,protected                :: print_sigma_
  logical,protected                :: print_restart_
  logical,protected                :: print_bigrestart_
  logical,protected                :: print_pdos_
  logical,protected                :: print_spatial_extension_
  logical,protected                :: print_cube_
  logical,protected                :: print_wfn_files_
  logical,protected                :: print_all_MO_wfn_file_
  logical,protected                :: print_multipole_
  logical,protected                :: print_hartree_
  logical,protected                :: print_density_matrix_
  logical,protected                :: print_rho_grid_
  logical,protected                :: print_tddft_matrices_
  logical,protected                :: print_cube_rho_tddft_
  logical,protected                :: print_cube_diff_tddft_
  logical,protected                :: print_line_rho_tddft_
  logical,protected                :: print_line_rho_diff_tddft_
  logical,protected                :: print_dens_traj_tddft_
  logical,protected                :: print_dens_traj_
  logical,protected                :: print_dens_traj_points_set_
  logical,protected                :: print_charge_tddft_
  logical,protected                :: print_transition_density_
  logical,protected                :: cphf_cpks_0_
  logical,protected                :: calc_q_matrix_
  logical,protected                :: calc_dens_disc_
  logical,protected                :: calc_spectrum_
  logical,protected                :: read_tddft_restart_
  logical,protected                :: print_tddft_restart_
  logical,protected                :: print_yaml_
  logical,protected                :: assume_scf_converged_
  logical,protected                :: analytic_chi_
  logical,protected                :: eri3_genuine_
  logical,protected                :: auto_occupation_
  logical,protected                :: gwgwg_skip_vvv_
  logical,protected                :: gwgwg_skip_vv_
  logical,protected                :: gwgwg_static_approximation_

  real(dp),protected               :: rcut         = 0.0_dp
  real(dp),protected               :: factor_sosex = 1.0_dp

  ! Here we call the fortran code that was generated by the python script
  ! Any new variable should be added through the python script
#include "input_variable_declaration.f90"

contains


!=========================================================================
subroutine init_calculation_type(scf,postscf)
  implicit none

  character(len=*),intent(in) :: scf,postscf
  !=====
  !=====

  !
  ! default values
  calc_type%is_dft               = .FALSE.
  calc_type%need_rpa             = .FALSE.
  calc_type%is_lr_mbpt           = .FALSE.
  calc_type%is_gw                = .FALSE.
  calc_type%is_mp2               = .FALSE.
  calc_type%is_mp3               = .FALSE.
  calc_type%is_noft              = .FALSE.
  calc_type%is_ci                = .FALSE.
  calc_type%is_bse               = .FALSE.
  calc_type%no_bse_kernel        = .FALSE.
  calc_type%include_tddft_kernel = .FALSE.
  calc_type%include_tdhf_kernel  = .FALSE.
  calc_type%is_real_time         = .FALSE.
  calc_type%selfenergy_technique = one_shot
  calc_type%selfenergy_approx    = 0
  calc_type%postscf_name         = 'None'
  calc_type%is_selfenergy        = .FALSE.
  calc_type%is_core              = .FALSE.

  !
  ! If it exists, first read the last part of the calculation specifier
  if( LEN(TRIM(postscf)) > 0 ) then

    calc_type%postscf_name =  TRIM(postscf)

    select case(TRIM(postscf))
    case('GNW0')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GnW0
      calc_type%selfenergy_technique = EVSC
    case('GNWN','EVGW')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GnWn
      calc_type%selfenergy_technique = EVSC
    case('GW','G0W0')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GW
    case('G0W0_DYSON')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GW
      calc_type%selfenergy_technique = exact_dyson
    case('COHSEX')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = COHSEX
      calc_type%selfenergy_technique = static_selfenergy
    case('G0W0_AC','GW_AC','G0W0_PADE','GW_PADE')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx    = GW
      calc_type%selfenergy_technique = imaginary_axis_pade
    case('G0W0_CONTOUR','GW_CONTOUR')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx    = GW
      calc_type%selfenergy_technique = contour_deformation
    case('G0W0_HOMOLUMO','GW_HOMOLUMO')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx    = GW
      calc_type%selfenergy_technique = imaginary_axis_homolumo
    case('GWSOX','GW+SOX')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWSOX
    case('GWSOX_PADE','GW+SOX_PADE')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWSOX
      calc_type%selfenergy_technique = imaginary_axis_pade
    case('GWPT3')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWPT3
    case('GWSOSEX','GW+SOSEX','GW+GWGVG')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWSOSEX
    case('GWSOSEX_PADE','GW+SOSEX_PADE')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWSOSEX
      calc_type%selfenergy_technique = imaginary_axis_pade
    case('GWSOSEX2','GW+SOSEX2','GW+2SOSEX')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWSOSEX
      factor_sosex = 2.0_dp
    case('GWSOSEX2_PADE','GW+SOSEX2_PADE','GW+2SOSEX_PADE')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWSOSEX
      factor_sosex = 2.0_dp
      calc_type%selfenergy_technique = imaginary_axis_pade
    case('GWGWGWG','GW+GWGWG','GW+G3W2')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWGWG
      factor_sosex = 2.0_dp
    case('GWGWGWG_PADE','GW+GWGWG_PADE','GW+G3W2_PADE')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWGWG
      factor_sosex = 2.0_dp
      calc_type%selfenergy_technique = imaginary_axis_pade
    case('GWGWGWG_NUMERICAL','GW+GWGWG_NUMERICAL','GW+G3W2_NUMERICAL')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWGWG_NUMERICAL
      factor_sosex = 2.0_dp
    case('GW+GW0GW0G')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GW0GW0G
    case('GW+GWGW0G')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWGW0G
    case('GW+GWGW0RPAG')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWGW0RPAG
    case('GWTILDE')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWTILDE
    case('EVGWGAMMA','GNW0GAMMAN','GNW0SOSEX','EVGWSOSEX')
      calc_type%is_gw    =.TRUE.
      calc_type%selfenergy_approx = GWSOSEX
      calc_type%selfenergy_technique = EVSC
    case('LRGW')
      calc_type%is_gw      =.TRUE.
      calc_type%selfenergy_approx = GW
      calc_type%is_lr_mbpt = .TRUE.
    case('MP2')
      calc_type%is_mp2   =.TRUE.
    case('MP3')
      calc_type%is_mp2   =.TRUE.
      calc_type%is_mp3   =.TRUE.
    case('MP2_SELFENERGY','PT2','GF2')
      calc_type%selfenergy_approx = PT2
    case('MP3_SELFENERGY','PT3','GF3')
      calc_type%selfenergy_approx = PT3
    case('EVMP2_SELFENERGY','EVPT2','EVGF2')
      calc_type%selfenergy_approx = PT2
      calc_type%selfenergy_technique = EVSC
    case('EVMP3_SELFENERGY','EVPT3','EVGF3')
      calc_type%selfenergy_approx = PT3
      calc_type%selfenergy_technique = EVSC
    case('NOFT')
      calc_type%is_noft   =.TRUE.
    case('TWO_RINGS','TWO-RINGS','TWORINGS','2RINGS')
      calc_type%selfenergy_approx = TWO_RINGS
    case('ONE_RING','ONE-RING','ONERING','1RING')
      calc_type%selfenergy_approx = ONE_RING
    case('SOX')
      calc_type%selfenergy_approx = SOX
    case('CI')
      calc_type%is_ci = .TRUE.
    case('CI_SELFENERGY')
      calc_type%is_ci = .TRUE.
      calc_type%selfenergy_approx = CI
    case('BSE')
      calc_type%is_bse = .TRUE.
    case('BSE_RPA','BSE-RPA') ! debug only
      calc_type%is_bse        = .TRUE.
      calc_type%no_bse_kernel = .TRUE.
    case('TD')
      calc_type%include_tddft_kernel = .TRUE.
    case('CPHF','CPKS')
      calc_type%include_tddft_kernel = .TRUE.
    case('REAL_TIME')
      calc_type%is_real_time = .TRUE.
    case('RPA','RPAP','RPA_IM','RPAP_IM','RPA+','RPA+_IM','RPA-I')
      ! nothing to declare
    case('BSE-I')
      calc_type%is_bse        = .TRUE.
      ! nothing to declare
    case('RPALR')
      calc_type%is_lr_mbpt = .TRUE.
    case('RPAX','RPAX-I','RPAX-II')
      calc_type%include_tdhf_kernel =.TRUE.
    case default
      write(stdout,*) 'postscf: ',TRIM(postscf)
      call die('init_calculation_type: Error reading keyword: postscf')
    end select
  endif

  calc_type%scf_name =  TRIM(scf)

  !
  ! Then read the first part of the calculation specifier
  select case(TRIM(scf))
  case('CI')
    calc_type%is_ci         = .TRUE.
    alpha_hybrid            = 1.00_dp
  case('CORE')
    alpha_hybrid            = 0.00_dp
    calc_type%is_core       = .TRUE.
  case('H','HARTREE')
    alpha_hybrid            = 0.0_dp
  case('HF')
    alpha_hybrid            = 1.00_dp
  case('MP2','PT2')
    calc_type%selfenergy_approx = PT2
    calc_type%selfenergy_technique = QS
    alpha_hybrid            = 1.00_dp
  case('GW','QSGW')
    calc_type%is_gw                = .TRUE.
    calc_type%selfenergy_approx    = GW
    calc_type%selfenergy_technique = QS
    alpha_hybrid            = 1.00_dp
  case('COHSEX')
    calc_type%is_gw                = .TRUE.           !ABCD
    calc_type%selfenergy_approx    = COHSEX
    calc_type%selfenergy_technique = QS
    alpha_hybrid            = 1.00_dp
  case default
    !
    ! If the calculation type is none of the above, let's assume it is DFT-type
    calc_type%is_dft=.TRUE.
    call init_dft_type(scf)
  end select

  !
  ! Do we need Coulomb integrals?
  ! Do we need LR Coulomb integrals?
  !
  calc_type%need_exchange    = ( alpha_hybrid > 1.0e-6 )
  calc_type%need_exchange_lr = ( rcut > 1.0e-6 ) .AND. ( ABS(beta_hybrid) > 1.0e-6 )

  calc_type%is_selfenergy = ( calc_type%selfenergy_approx /= 0 )

end subroutine init_calculation_type


!=========================================================================
subroutine init_excitation_type()
  implicit none
  !=====

  excit_type%name  = excit_name
  excit_type%kappa = excit_kappa
  excit_type%omega = excit_omega
  excit_type%time0 = excit_time0
  excit_type%dir   = excit_dir

  if( LEN(TRIM(excit_name)) /= 0 ) then
    select case (excit_type%name)
    case("ION","ANTIION")
      excit_type%form = EXCIT_PROJECTILE_W_BASIS
    case("NUCLEUS","ANTINUCLEUS")
      excit_type%form = EXCIT_PROJECTILE
    case("NO")
      excit_type%form = EXCIT_NO
    case("GAU","HSW","STEP","DEL")
      excit_type%form = EXCIT_LIGHT
    case default
      write(stdout,*) 'error reading excitation name (excit_name variable)'
      write(stdout,*) TRIM(excit_name)
      call die('excit_name is unknown')
    end select
  end if

end subroutine init_excitation_type


!=========================================================================
subroutine init_dft_type(key)
  implicit none

  character(len=*),intent(in)          :: key
  !=====
  integer              :: ixc,off1,off2
  real(C_DOUBLE)       :: globalx_libxc,srx_libxc,omega_libxc
  !=====


  if( ALLOCATED(dft_xc) ) then
    call destroy_libxc_info(dft_xc)
  endif

  !
  ! Prepare the object dft_xc
  allocate(dft_xc(3))
  dft_xc(:)%nspin = nspin
  ! default is one, otherwise it is modified later
  dft_xc(:)%coeff = 1.0_dp
  ! id = 0 means not xc
  dft_xc(:)%id    = 0



  select case(TRIM(key))
#if !defined(NO_LIBXC)
  !
  ! LDA functionals
  case('LDAX')
    dft_xc(1)%id = XC_LDA_X
    alpha_hybrid   = 0.00_dp
  case('SPL')
    dft_xc(1)%id = XC_LDA_X
    dft_xc(2)%id = XC_LDA_C_PZ
    alpha_hybrid   = 0.00_dp
  case('LDA')
    dft_xc(1)%id = XC_LDA_X
    dft_xc(2)%id = XC_LDA_C_PW
    alpha_hybrid   = 0.00_dp
  case('VWN')
    dft_xc(1)%id = XC_LDA_X
    dft_xc(2)%id = XC_LDA_C_VWN
    alpha_hybrid   = 0.00_dp
  case('VWN_RPA')
    dft_xc(1)%id = XC_LDA_X
    dft_xc(2)%id = XC_LDA_C_VWN_RPA
    alpha_hybrid   = 0.00_dp
  case('VWNC')
    dft_xc(1)%id = XC_LDA_C_VWN
    alpha_hybrid   = 0.00_dp
  !
  ! GGA functionals
  case('PBEX')
    dft_xc(1)%id    = XC_GGA_X_PBE
    alpha_hybrid    = 0.00_dp
  case('PBE')
    dft_xc(1)%id    = XC_GGA_X_PBE
    dft_xc(2)%id    = XC_GGA_C_PBE
    alpha_hybrid    = 0.00_dp
  case('PBE_SOL','PBESOL')
    dft_xc(1)%id    = XC_GGA_X_PBE_SOL
    dft_xc(2)%id    = XC_GGA_C_PBE_SOL
    alpha_hybrid    = 0.00_dp
  case('WPBEHX')
    dft_xc(1)%id    = XC_GGA_X_WPBEH
    alpha_hybrid    = 0.00_dp
    dft_xc(1)%gamma = gamma_hybrid
  case('WPBEH')
    dft_xc(1)%id    = XC_GGA_X_WPBEH
    dft_xc(2)%id    = XC_GGA_C_PBE
    alpha_hybrid    = 0.00_dp
    dft_xc(1)%gamma = gamma_hybrid
  case('BX')
    dft_xc(1)%id = XC_GGA_X_B88
    alpha_hybrid   = 0.00_dp
  case('BLYP')
    dft_xc(1)%id = XC_GGA_X_B88
    dft_xc(2)%id = XC_GGA_C_LYP
    alpha_hybrid   = 0.00_dp
  case('PW91X')
    dft_xc(1)%id = XC_GGA_X_PW91
    alpha_hybrid   = 0.00_dp
  case('PW91')
    dft_xc(1)%id = XC_GGA_X_PW91
    dft_xc(2)%id = XC_GGA_C_PW91
    alpha_hybrid   = 0.00_dp
  case('HCTH')
    dft_xc(1)%id = XC_GGA_XC_HCTH_407
    alpha_hybrid   = 0.00_dp
  case('TH')
    dft_xc(1)%id = XC_GGA_XC_TH1
    alpha_hybrid   = 0.00_dp
  case('HJSX')
    dft_xc(1)%id = XC_GGA_X_HJS_PBE
    alpha_hybrid   = 0.00_dp
    dft_xc(1)%gamma = gamma_hybrid
  case('LYP')
    dft_xc(1)%id = XC_GGA_C_LYP
    alpha_hybrid   = 0.00_dp
  case('PBEC')
    dft_xc(1)%id = XC_GGA_C_PBE
    alpha_hybrid   = 0.00_dp
  !
  ! Meta-GGA functionals
  case('RPPX')
    dft_xc(1)%id = XC_MGGA_X_RPP09
    alpha_hybrid = 0.00_dp
  !
  ! Hybrid functionals
  case('HFPBE')
    dft_xc(1)%id = XC_GGA_C_PBE
    alpha_hybrid = 1.00_dp
  case('BHANDH')
    dft_xc(1)%id = XC_HYB_GGA_XC_BHANDH
    alpha_hybrid = 0.50_dp
  case('BHANDHLYP','BHLYP')
    dft_xc(1)%id = XC_HYB_GGA_XC_BHANDHLYP
    alpha_hybrid = 0.50_dp
  case('B3LYP')
    dft_xc(1)%id = XC_HYB_GGA_XC_B3LYP
    alpha_hybrid = 0.20_dp
  case('B3LYP5')
    dft_xc(1)%id = XC_HYB_GGA_XC_B3LYP5
    alpha_hybrid = 0.20_dp
  case('PBE0')
    dft_xc(1)%id = XC_HYB_GGA_XC_PBEH
    alpha_hybrid   = 0.25_dp
  case('PBE50')
    dft_xc(1)%id = XC_HYB_GGA_XC_PBE50
    alpha_hybrid   = 0.50_dp
  case('WB97')
    dft_xc(1)%id = XC_HYB_GGA_XC_WB97
    alpha_hybrid  = 0.0_dp
    beta_hybrid   = 1.0_dp
    gamma_hybrid  = 0.40_dp
  case('WB97X')
    dft_xc(1)%id = XC_HYB_GGA_XC_WB97X
    alpha_hybrid  = 0.157706_dp
    beta_hybrid   = 1.0_dp - alpha_hybrid
    gamma_hybrid  = 0.30_dp
  case('HSE03')
    dft_xc(1)%id  = XC_HYB_GGA_XC_HSE03
    alpha_hybrid  = 0.25_dp
    beta_hybrid   = -alpha_hybrid
    gamma_hybrid  = 0.15_dp / SQRT(2.0_dp)
  case('HSE06')
    dft_xc(1)%id  = XC_HYB_GGA_XC_HSE06
    alpha_hybrid  = 0.25_dp
    beta_hybrid   = -alpha_hybrid
    gamma_hybrid  = 0.11_dp
  case('HSE08')
    dft_xc(1)%id  = XC_HYB_GGA_XC_HJS_PBE
    alpha_hybrid  = 0.25_dp
    beta_hybrid   = -alpha_hybrid
    gamma_hybrid  = 0.11_dp
  case('LC-BLYP') ! Notice that it is built to be also used as a double-hybrid functional
    dft_xc(1)%id  =  XC_GGA_X_ITYH
    dft_xc(2)%id  =  XC_GGA_C_LYP
    if( abs(kappa_hybrid) < tol8 ) kappa_hybrid=0.00_dp
    if( abs( gamma_hybrid - 1000000.0_dp ) < tol8 ) gamma_hybrid=0.47_dp
    alpha_hybrid  = 0.00_dp
    beta_hybrid   = 1.00_dp
    dft_xc(1)%coeff = 1.00_dp
    dft_xc(2)%coeff = 1.00_dp - kappa_hybrid
    dft_xc(1)%gamma = gamma_hybrid
  case('CAM-B3LYP')
    dft_xc(1)%id  = XC_HYB_GGA_XC_CAM_B3LYP
    alpha_hybrid  = 0.19_dp
    beta_hybrid   = 0.46_dp
    gamma_hybrid  = 0.33_dp
  case('TUNED-CAM-B3LYP')
    dft_xc(1)%id  = XC_HYB_GGA_XC_TUNED_CAM_B3LYP
    alpha_hybrid  = 0.0799_dp
    beta_hybrid   = 0.9201_dp
    gamma_hybrid  = 0.150_dp
    dft_xc(2)%gamma = gamma_hybrid
  case('RSHX')
    dft_xc(1)%id = XC_GGA_X_PBE
    dft_xc(2)%id = XC_GGA_X_HJS_PBE
    dft_xc(1)%coeff = 1.00_dp - (alpha_hybrid + beta_hybrid)
    dft_xc(2)%coeff = beta_hybrid
    dft_xc(2)%gamma = gamma_hybrid
  case('LDA0')
    alpha_hybrid = 0.25_dp
    dft_xc(1)%id = XC_LDA_X
    dft_xc(2)%id = XC_LDA_C_PW
    dft_xc(1)%coeff = 1.00_dp - alpha_hybrid
    dft_xc(2)%coeff = 1.00_dp
  case('PBEH')
    dft_xc(1)%id = XC_GGA_X_PBE
    dft_xc(2)%id = XC_GGA_C_PBE
    dft_xc(1)%coeff = 1.00_dp - alpha_hybrid
    dft_xc(2)%coeff = 1.00_dp
  case('RSH')   ! This one can also be used for double hybrid functionals (e.g. PBEQIDH, PBE0-DH, and their RPA+ versions).
    dft_xc(1)%id = XC_GGA_X_PBE
    dft_xc(2)%id = XC_GGA_X_HJS_PBE
    dft_xc(3)%id = XC_GGA_C_PBE
    dft_xc(1)%coeff = 1.00_dp - (alpha_hybrid + beta_hybrid)
    dft_xc(2)%coeff = beta_hybrid
    dft_xc(3)%coeff = 1.00_dp - kappa_hybrid
    dft_xc(2)%gamma = gamma_hybrid

  !
  ! Double Hybrid functionals (used with postscf='MP2' or 'RPA')
  case('PBE0-DH')
    if( abs(alpha_hybrid) < tol8 ) alpha_hybrid=0.5_dp
    if( abs(kappa_hybrid) < tol8 ) kappa_hybrid=0.125_dp
    dft_xc(1)%id = XC_GGA_X_PBE
    dft_xc(2)%id = XC_GGA_C_PBE
    dft_xc(1)%coeff = 1.00_dp - alpha_hybrid
    dft_xc(2)%coeff = 1.00_dp - kappa_hybrid
  case('PBE-QIDH')
    if( abs(alpha_hybrid) < tol8 ) alpha_hybrid=0.693361274_dp
    if( abs(kappa_hybrid) < tol8 ) kappa_hybrid=0.333333333_dp
    dft_xc(1)%id = XC_GGA_X_PBE
    dft_xc(2)%id = XC_GGA_C_PBE
    dft_xc(1)%coeff = 1.00_dp - alpha_hybrid
    dft_xc(2)%coeff = 1.00_dp - kappa_hybrid
  case('B2PLYP')
    if( abs(alpha_hybrid) < tol8 ) alpha_hybrid=0.53_dp
    if( abs(kappa_hybrid) < tol8 ) kappa_hybrid=0.27_dp
    dft_xc(1)%id = XC_GGA_X_B88
    dft_xc(2)%id = XC_GGA_C_LYP
    dft_xc(1)%coeff = 1.00_dp - alpha_hybrid
    dft_xc(2)%coeff = 1.00_dp - kappa_hybrid
  case('RSX-QIDH') ! Using beta = 1.0 - alpha in RSX-QIDH.
    if( abs(alpha_hybrid) < tol8 ) alpha_hybrid=0.693361274_dp
    if( abs(kappa_hybrid) < tol8 ) kappa_hybrid=0.333333333_dp
    if( abs( gamma_hybrid - 1000000.0_dp ) < tol8 ) gamma_hybrid=0.27_dp
    beta_hybrid  = 1.00_dp - alpha_hybrid
    dft_xc(1)%id = XC_GGA_X_ITYH_PBE
    dft_xc(2)%id = XC_GGA_C_PBE
    dft_xc(1)%coeff = beta_hybrid
    dft_xc(2)%coeff = 1.00_dp - kappa_hybrid
    dft_xc(1)%gamma = gamma_hybrid
  case('RSH-QIDH') ! Forcing beta = 1.0 - alpha in RSH (and XC_GGA_X_PBE cancels).
    if( abs(alpha_hybrid) < tol8 ) alpha_hybrid=0.693361274_dp
    if( abs(kappa_hybrid) < tol8 ) kappa_hybrid=0.333333333_dp
    if( abs( gamma_hybrid - 1000000.0_dp ) < tol8 ) gamma_hybrid=0.27_dp
    beta_hybrid = 1.00_dp - alpha_hybrid
    dft_xc(1)%id = XC_GGA_X_HJS_PBE
    dft_xc(2)%id = XC_GGA_C_PBE
    dft_xc(1)%coeff = beta_hybrid
    dft_xc(2)%coeff = 1.00_dp - kappa_hybrid
    dft_xc(1)%gamma = gamma_hybrid
#endif
  case default

    !
    ! Possibility offered to enter directly the LIBXC index of the functional
    ! using  LIBXC:101+130 for PBE
    !
    ! Check whether the key starts with "LIBXC:")
    off1 = INDEX(key,"LIBXC:")
    if( off1 > 0 ) then
      write(stdout,'(1x,a,a)') 'Reading a manual input of the LIBXC indeces: ',TRIM(key)
      off1 = off1 + 6

      ixc = 0
      off2 = off1
      do while( off2 == off1 )
        ! Look for "+" sign to assemble LIBXC functionals
        off2 = INDEX(key(off1:),"+") - 2 + off1
        ixc = ixc + 1
        if( ixc > SIZE(dft_xc) ) &
          call die('init_dft_type: too many functionals in the input string. Maximum is 3')
        if( off2 < off1 ) then
          read(key(off1:),'(i4)') dft_xc(ixc)%id
        else
          read(key(off1:off2),'(i4)') dft_xc(ixc)%id
          off2 = off2+2
          off1 = off2
        endif
      enddo

    else
      call die('init_dft_type: Error reading keyword scf. DFT type not understood')
    endif
  end select

  call init_libxc_info(dft_xc)

#if !defined(NO_LIBXC)
  !
  ! If "LIBXC:" syntax, obtain exact-exchange parameters from LIBXC
  if( INDEX(key,"LIBXC:") > 0 ) then
    ! Assumes hybrid functionals consists of only one XC id
    call xc_hyb_cam_coef(dft_xc(1)%func,omega_libxc,globalx_libxc,srx_libxc)
    gamma_hybrid = omega_libxc
    alpha_hybrid = globalx_libxc + srx_libxc
    beta_hybrid  = -srx_libxc
  endif
#endif

  if( gamma_hybrid > 1.0e-12 ) then
    rcut = 1.0_dp / gamma_hybrid
  endif

end subroutine init_dft_type


!=========================================================================
subroutine summary_input()
  implicit none

  !=====

  !
  ! Summarize some important input parameters
  write(stdout,'(/,a,/)')    ' Summary of important input parameters '
  write(stdout,'(a30,2x,a)') '         SCF type: ',calc_type%scf_name
  write(stdout,'(a30,2x,a)') '    Post SCF type: ',calc_type%postscf_name
  write(stdout,'(a30,i4)')   ' number of atoms: ',ncenter_nuclei
  write(stdout,'(a30,i4)')   ' number of basis centers: ',ncenter_basis
  write(stdout,'(a30,f10.4)') ' electrons: ',electrons
  write(stdout,'(a30,f10.4)') ' charge: ',charge
  write(stdout,'(a30,i4)')   ' spin polarization: ',nspin
  write(stdout,'(a30,f10.4)') ' magnetization: ',magnetization
  write(stdout,'(a30,2x,a)') ' basis file path:',TRIM(basis_path)
  write(stdout,'(a30,2x,a)') ' basis set: ',basis_name(1)
  write(stdout,'(a30,2x,a)') ' auxiliary basis set: ',auxil_basis_name(1)
  write(stdout,'(a30,2x,a)') ' gaussian type: ',TRIM(gaussian_type)
  write(stdout,'(a30,f10.4)') ' global exchange: ',alpha_hybrid
  write(stdout,'(a30,f10.4)') ' long-range-only exchange: ',beta_hybrid
  if( rcut > 1.0e-6_dp ) then
    write(stdout,'(a30,f8.4)') ' range-separation (bohr-1): ',gamma_hybrid
  else
    write(stdout,'(a30,a)') ' range-separation (bohr-1): ',' none'
  endif
  write(stdout,*)

  call output_positions()

  write(stdout,'(a,i5)') ' Number of bonds ',nbond
  if(inversion) then
    write(stdout,*) 'Molecule has inversion symmetry'
  else
    write(stdout,*) 'Molecule does not have inversion symmetry'
  endif
  if(linear) then
    write(stdout,*) 'Molecule is linear'
  else
    write(stdout,*) 'Molecule is not linear'
    if(planar) then
      write(stdout,*) 'Molecule is planar'
    else
      write(stdout,*) 'Molecule is not planar'
    endif
  endif
  write(stdout,*)

end subroutine summary_input


!=========================================================================
subroutine read_inputfile_namelist()
  implicit none

  !=====
  integer              :: idot,ichart,ninput_argument
  character(len=140)   :: input_file_name
  integer              :: inputfile
  logical              :: file_exists

  character(len=140)   :: read_restart
  character(len=140)   :: basis
  character(len=140)   :: auxil_basis
  character(len=140)   :: small_basis
  character(len=140)   :: ecp_basis
  character(len=140)   :: ecp_auxil_basis
  character(len=140)   :: ecp_small_basis
  character(len=256)   :: default_basis_path
  character(len=256)   :: basis_path_env
  integer              :: istatus,icenter

  ! Here we call the fortran code that was generated by the python script
  ! Any new variable should be added through the python script
#include "input_variables.f90"

  ! Here we call the fortran code that was generated by the python script
#include "basis_path.f90"
  !=====


  ! If no value is given to basis_path in the input file,
  ! then try to obtain it from the environment variable MOLGW_BASIS_PATH
  ! and if not defined, then get the default value in the source
  ! basis_path is set in this priority order:
  ! input file > environment variable > default from the sources
  if( LEN(TRIM(basis_path)) == 0 ) then
    call GET_ENVIRONMENT_VARIABLE('MOLGW_BASIS_PATH',basis_path_env,status=istatus)
    if( istatus == 0 ) then
      write(stdout,'(1x,a,a)') 'Found environment variable: MOLGW_BASIS_PATH=',TRIM(basis_path_env)
      basis_path = TRIM(basis_path_env)
    else
      basis_path = default_basis_path
    endif
  endif


  ! Get the number of inline arguments with the new Fortran 2003 statement
  ninput_argument = COMMAND_ARGUMENT_COUNT()

  select case(ninput_argument)
  case(1)
    call GET_COMMAND_ARGUMENT(1,VALUE=input_file_name)
    write(stdout,'(a,a)') ' Opening input file: ',TRIM(input_file_name)
    inquire(file=TRIM(input_file_name),exist=file_exists)
    if( .NOT. file_exists) then
      write(stdout,*) 'Tried to open file:',TRIM(input_file_name)
      call die('Input file not found')
    endif
    output_name=TRIM(input_file_name)
    idot=1
    ichart=1
    do
      if(output_name(ichart:ichart+2)==".in") then
        idot=ichart
        exit
      endif
      ichart=ichart+1
      if(ichart+2>140) exit
    enddo
    output_name=TRIM(output_name(1:idot))
    open(newunit=inputfile,file=TRIM(input_file_name),status='old')
  case(0)
    inputfile = 5
    call issue_warning('Deprecated reading from stdin. Please use instead the newer syntax ./molgw inputfile > outfile')
  case default
    call die('input file name not understood')
  end select

  !
  ! Read all the input file in one statement!
  !
  read(inputfile,molgw)


  ! Here we call the fortran code that was generated by the python script
  write(stdout,'(/,1x,a,/)') ' === Echo the entire list of input variables  ==='
#include "echo_input_variables.f90"
  write(stdout,'(/,1x,a,/)') ' ================================================'


  ieta = (0.0_dp,1.0_dp) * eta

  scf                = capitalize(scf)
  postscf            = capitalize(postscf)
  gaussian_type      = capitalize(gaussian_type)
  tddft_grid_quality = capitalize(tddft_grid_quality)
  grid_quality       = capitalize(grid_quality)
  ecp_quality        = capitalize(ecp_quality)
  integral_quality   = capitalize(integral_quality)
  mixing_scheme      = capitalize(mixing_scheme)
  length_unit        = capitalize(length_unit)
  init_hamiltonian   = capitalize(init_hamiltonian)
  prop_type          = capitalize(prop_type)
  excit_name         = capitalize(excit_name)
  pred_corr          = capitalize(pred_corr)
  ci_greens_function = capitalize(ci_greens_function)
  ci_type            = capitalize(ci_type)
  read_fchk          = capitalize(read_fchk)
  pt_density_matrix  = capitalize(pt_density_matrix)
  pt3_a_diagrams     = capitalize(pt3_a_diagrams)
  partition_scheme   = capitalize(partition_scheme)
  w_screening        = capitalize(w_screening)

  memory_evaluation_        = yesno_to_logical(memory_evaluation)
  read_restart_             = yesno_to_logical(read_restart)
  ignore_bigrestart_        = yesno_to_logical(ignore_bigrestart)
  force_energy_qp_          = yesno_to_logical(force_energy_qp)
  is_tda                    = yesno_to_logical(tda)
  is_triplet                = yesno_to_logical(triplet)
  is_frozencore             = yesno_to_logical(frozencore)
  is_tddft_frozencore       = yesno_to_logical(tddft_frozencore)
  is_virtual_fno            = yesno_to_logical(virtual_fno)
  incore_                   = yesno_to_logical(incore)

  gwgwg_skip_vvv_           = yesno_to_logical(gwgwg_skip_vvv)
  gwgwg_skip_vv_            = yesno_to_logical(gwgwg_skip_vv)
  gwgwg_static_approximation_ = yesno_to_logical(gwgwg_static_approximation)
  print_eri_                = yesno_to_logical(print_eri)
  print_wfn_                = yesno_to_logical(print_wfn)
  print_w_                  = yesno_to_logical(print_w)
  print_sigma_              = yesno_to_logical(print_sigma)
  print_restart_            = yesno_to_logical(print_restart)
  print_bigrestart_         = yesno_to_logical(print_bigrestart)
  print_pdos_               = yesno_to_logical(print_pdos)
  print_charge_tddft_       = yesno_to_logical(print_charge_tddft)
  print_spatial_extension_  = yesno_to_logical(print_spatial_extension)
  print_multipole_          = yesno_to_logical(print_multipole)
  print_cube_               = yesno_to_logical(print_cube)
  print_wfn_files_          = yesno_to_logical(print_wfn_files)
  print_all_MO_wfn_file_    = yesno_to_logical(print_all_MO_wfn_file)
  print_hartree_            = yesno_to_logical(print_hartree)
  print_density_matrix_     = yesno_to_logical(print_density_matrix)
  print_rho_grid_           = yesno_to_logical(print_rho_grid)
  gwgamma_tddft_            = yesno_to_logical(gwgamma_tddft)
  use_correlated_density_matrix_ = yesno_to_logical(use_correlated_density_matrix)
  print_tddft_matrices_       = yesno_to_logical(print_tddft_matrices)
  print_cube_rho_tddft_       = yesno_to_logical(print_cube_rho_tddft)
  print_cube_diff_tddft_      = yesno_to_logical(print_cube_diff_tddft)
  print_line_rho_tddft_       = yesno_to_logical(print_line_rho_tddft)
  print_line_rho_diff_tddft_  = yesno_to_logical(print_line_rho_diff_tddft)
  print_dens_traj_tddft_      = yesno_to_logical(print_dens_traj_tddft)
  print_dens_traj_            = yesno_to_logical(print_dens_traj)
  print_dens_traj_points_set_ = yesno_to_logical(print_dens_traj_points_set)
  print_transition_density_   = yesno_to_logical(print_transition_density)
  calc_q_matrix_              = yesno_to_logical(calc_q_matrix)
  calc_dens_disc_             = yesno_to_logical(calc_dens_disc)
  calc_spectrum_              = yesno_to_logical(calc_spectrum)
  read_tddft_restart_         = yesno_to_logical(read_tddft_restart)
  print_tddft_restart_        = yesno_to_logical(print_tddft_restart)
  print_yaml_                 = yesno_to_logical(print_yaml)
  assume_scf_converged_       = yesno_to_logical(assume_scf_converged)
  cphf_cpks_0_                = yesno_to_logical(cphf_cpks_0)
  analytic_chi_               = yesno_to_logical(analytic_chi)
  eri3_genuine_               = yesno_to_logical(eri3_genuine)
  auto_occupation_            = yesno_to_logical(auto_occupation)

  tddft_grid_level   = interpret_quality(tddft_grid_quality)
  grid_level         = interpret_quality(grid_quality)
  ecp_level          = interpret_quality(ecp_quality)
  integral_level     = interpret_quality(integral_quality)


  select case(TRIM(mixing_scheme))
  case('SIMPLE','PULAY','DIIS','ADIIS','EDIIS')
  case default
    write(stdout,*) TRIM(mixing_scheme)
    call die('mixing scheme not recognized')
  end select

  select case(TRIM(pt3_a_diagrams))
  case('YES','NO','ONLY')
  case default
    call die('pt3_a_diagram input variable can only be yes, no, or only')
  end select

  !
  ! A few consistency checks
  !
  if(alpha_mixing<0.0 .OR. alpha_mixing > 1.0 ) call die('alpha_mixing should be inside [0,1]')
  if(ncoreg<0) call die('negative ncoreg is meaningless')
  if(ncorew<0) call die('negative ncorew is meaningless')
  if(ncore_tddft<0) call die('negative ncore_tddft is meaningless')
  if(nvirtualg<0) call die('negative nvirtualg is meaningless')
  if(nvirtualw<0) call die('negative nvirtualw is meaningless')
  if(nvirtualg<ncoreg) call die('too small nvirtualg is meaningless')
  if(nvirtualw<ncorew) call die('too small nvirtualw is meaningless')
  if(nspin/=1 .AND. nspin/=2) call die('nspin in incorrect')
  if(magnetization<-1.e-5)    call die('magnetization is negative')
  if(magnetization>1.e-5 .AND. nspin==1) call die('magnetization is non-zero and nspin is 1')
  if(nomega_sigma<0)    call die('nomega_sigma < 0')
  if(step_sigma<0.0_dp) call die('step_sigma < 0.0')
  if(auto_auxil_fsam<1.00001_dp) call die('auto_auxil_fsam should be strictly greater to 1. Increase it a bit please')


#if !defined(LIBINT2_DERIV_ONEBODY_ORDER) || (LIBINT2_DERIV_ONEBODY_ORDER == 0) || !defined(LIBINT2_DERIV_ERI_ORDER) || (LIBINT2_DERIV_ERI_ORDER == 0)
  if( move_nuclei /= 'no' ) then
    call die('LIBINT does not contain the gradients of the integrals that are needed when move_nuclei is different from no')
  endif
#endif

  if( mpi_nproc_ortho > world%nproc ) then
    mpi_nproc_ortho = world%nproc
    call issue_warning('mpi_nproc_ortho has been resized to the max number of processors')
    write(stdout,'(1x,a,i4)') 'Now mpi_nproc_ortho = ',mpi_nproc_ortho
  endif
  if( MODULO( world%nproc , mpi_nproc_ortho) /= 0 ) then
    write(stdout,'(1x,a,i6,a,i6)') 'mpi_nproc_ortho must be a divisor of nproc ',mpi_nproc_ortho,' / ',world%nproc
    mpi_nproc_ortho = 1
    call issue_warning('mpi_nproc_ortho value is invalid. Override it and set mpi_nproc_ortho=1')
  endif

  call init_excitation_type()
  nprojectile = MERGE(1,0,excit_type%form==EXCIT_PROJECTILE .OR. excit_type%form == EXCIT_PROJECTILE_W_BASIS)

  !
  ! If no nuclei motion is requested, then override nstep and set it to 1
  if( move_nuclei == 'no' ) then
    nstep = 1
  endif

  call setup_nuclei(inputfile,basis,auxil_basis,small_basis,ecp_basis,ecp_auxil_basis,ecp_small_basis)

  has_auxil_basis = TRIM(auxil_basis_name(1)) /= '' .OR. TRIM(ecp_auxil_basis_name(1)) /= ''
  has_small_basis = TRIM(small_basis_name(1)) /= '' .OR. TRIM(ecp_small_basis_name(1)) /= ''

  ! To avoid SCREENED_COULOMB bug:
  ! SCREENED_COULOMB file unfortunately depends on the diagonalization of the 2-center Coulomb repulsion of the the auxiliary basis.
  ! It is therefore runtime-dependent and should not be used!
  if( print_w_ .AND. has_auxil_basis ) then
    call die('input check: print_w is not numerically stable when using an auxiliary basis.' // &
             ' Do not use this keyword and everything is gonna be alright')
  endif

  !
  ! Interpret the scf and postscf input parameters
  call init_calculation_type(scf,postscf)

  !
  ! Some additional checks
  !
  if(calc_type%selfenergy_technique == imaginary_axis_pade .AND. nomega_chi_imag<1) &
    call die('when asking for a numerical evaluation of the self-energy, one needs nomega_chi_imag > 0')
  if(calc_type%selfenergy_technique == imaginary_axis_pade .AND. nomega_sigma_calc==1) &
    call issue_warning('when asking for a numerical evaluation of the self-energy,' &
                    // ' consider more frequencies than just one for sigma.' &
                    // ' nomega_sigma_calc > 1 advised')
  if( nexcitation /=0 .AND. calc_type%is_gw ) then
    call die('Davidson diago is not compatible with GW. Set nexcitation to 0')
  endif
  if( nstep_gw > 1 .AND. calc_type%selfenergy_technique /= EVSC ) then
    call die('nstep_gw > 1 is only valid when performing ev-GW. Change either postscf or nstep_gw')
  endif
  if( eri3_genuine_ .AND. ( calc_type%need_exchange .OR. calc_type%need_exchange_lr ) ) then
    call die('eri3_genuine does not work with exact-exchange')
  endif
  if( excit_type%form == EXCIT_PROJECTILE_W_BASIS .AND. .NOT.(eri3_genuine_) ) then
    call die('eri3_genuine is required for moving basis (=excit_name=ion)')
  endif
  if( excit_type%form == EXCIT_PROJECTILE_W_BASIS .AND. .NOT.(pred_corr(1:2)=='MB') ) then
    call die('Predictor-correction scheme is not valid for moving basis. Use instead MB_PC2B for instance')
  endif

  spin_fact = REAL(-nspin+3,dp)
  electrons = SUM(zvalence(:)) - charge

  ! Echo the interpreted input variables
  call summary_input()

  !
  ! Here we open a YAML file and keep it open until "this_is_the_end"
  ! Here we call the fortran code that was generated by the python script
  !
  if( print_yaml_ .AND. is_iomaster ) then
    open(newunit=unit_yaml,file=filename_yaml,action='write')
    write(unit_yaml,'(a)') '---'
    write(unit_yaml,'(a)') 'input parameters:'
#include "echo_input_variables_yaml.f90"

    write(unit_yaml,'(/,a)') 'physical system:'
    write(unit_yaml,'(4x,a,1x,es18.8)') 'electrons:',electrons
    write(unit_yaml,'(4x,a)') 'length unit: bohr'
    write(unit_yaml,'(4x,a)') 'atom list:'
    do icenter=1,ncenter_nuclei
      write(unit_yaml,'(8x,a,"[ ",a2,", ",es18.8,", ",es18.8,", ",es18.8,"]")') '- ', &
                  element_name(REAL(zatom(icenter),dp)),xatom(:,icenter)
    enddo

    write(unit_yaml,'(4x,a)') 'basis list:'
    do icenter=1,ncenter_basis
      write(unit_yaml,'(8x,a,"[ ",a2,", ",a,", ",a,"]")') '- ', &
              element_name(REAL(zbasis(icenter),dp)),TRIM(basis_name(icenter)),TRIM(auxil_basis_name(icenter))
    end do

  endif


end subroutine read_inputfile_namelist


!=========================================================================
function interpret_quality(quality) result(quality_level)
  implicit none

  character(len=*),intent(in) :: quality
  integer                     :: quality_level
  !=====
  !=====

  select case(TRIM(quality))
  case('LOW','L')
    quality_level = low
  case('MEDIUM','MED','M')
    quality_level = medium
  case('HIGH','HI','H')
    quality_level = high
  case('VERY HIGH','VERYHIGH','VH')
    quality_level = very_high
  case('INSANE','I')
    quality_level = insane
  end select

end function interpret_quality


!=========================================================================
function standardize_basis_name(basis_name_in) result(basis_name_out)
  implicit none

  character(len=*),intent(in)       :: basis_name_in
  character(len=LEN(basis_name_in)) :: basis_name_out
  !=====
  integer :: istring
  !=====

  do istring=1,LEN(basis_name_in)

    basis_name_out(istring:istring) = basis_name_in(istring:istring)

    if( basis_name_in(istring:istring) == '*' ) then
      basis_name_out(istring:istring) = 's'
    endif

    if( basis_name_in(istring:istring) == '+' ) then
      basis_name_out(istring:istring) = 'p'
    endif

  enddo

end function standardize_basis_name


!=========================================================================
subroutine setup_nuclei(inputfile,basis,auxil_basis,small_basis,ecp_basis,ecp_auxil_basis,ecp_small_basis)
  implicit none

  integer,intent(in)          :: inputfile
  character(len=*),intent(in) :: basis,auxil_basis,small_basis
  character(len=*),intent(in) :: ecp_basis,ecp_small_basis,ecp_auxil_basis
  !=====
  integer              :: natom_read
  integer              :: element,iatom,ielement_ecp,icenter
  integer              :: info,info1,info2,xyzfile
  character(len=12)    :: element_symbol
  real(dp),allocatable :: zatom_read(:),x_read(:,:)
  character(len=140)   :: ctmp1,ctmp2
  character(len=256)   :: line_char
  logical              :: file_exists
  real(dp)             :: length_factor
  integer              :: ncenter_basis_max
  logical,allocatable  :: nucleus_wo_basis(:)
  !=====

  select case(TRIM(length_unit))
  case('A','ANGSTROM')
    length_factor = 1.0_dp/bohr_A
  case('BOHR','AU','A.U','A.U.')
    length_factor = 1.0_dp
  case default
    call die('units for lengths in input file not understood')
  end select
  if( LEN(TRIM(xyz_file)) > 0 ) then
    if( ABS(length_factor - 1.0_dp ) < 1.0e-6_dp  ) then
      write(stdout,*) 'xyz files are always in Angstrom. However, length_unit was set to bohr'
      call die('Please set length_unit to Angstrom')
    endif
  endif

  !
  ! Read the atom positions if no xyz file is specified
  if( LEN(TRIM(xyz_file)) == 0 ) then
    !
    ! In this case, natom must be set to a positive value
    !if(natom<1) call die('natom<1')

    if(excit_type%form == EXCIT_PROJECTILE_W_BASIS) then
      ncenter_basis_max = natom + nghost + nprojectile
    else
      ncenter_basis_max = natom + nghost
    endif
    ncenter_nuclei = natom + nprojectile

    natom_read     = natom + nghost + nprojectile
    allocate(nucleus_wo_basis(natom_read))
    nucleus_wo_basis(:) = .FALSE.

    !
    ! Need to know the number of atoms to allocate the basis arrays
    allocate(auxil_basis_name(ncenter_basis_max))
    allocate(small_basis_name(ncenter_basis_max))
    allocate(basis_name(ncenter_basis_max))
    allocate(ecp_basis_name(ncenter_basis_max))
    allocate(ecp_auxil_basis_name(ncenter_basis_max))
    allocate(ecp_small_basis_name(ncenter_basis_max))
    basis_name(:)           = standardize_basis_name(basis)
    auxil_basis_name(:)     = standardize_basis_name(auxil_basis)
    small_basis_name(:)     = standardize_basis_name(small_basis)
    ecp_basis_name(:)       = standardize_basis_name(ecp_basis)
    ecp_auxil_basis_name(:) = standardize_basis_name(ecp_auxil_basis)
    ecp_small_basis_name(:) = standardize_basis_name(ecp_small_basis)

    allocate(x_read(3,natom_read),zatom_read(natom_read))
    ncenter_basis = 0

    do iatom=1,natom_read
      ! First, read the full line
      read(inputfile,'(a)') line_char

      ! Then, try to interpret it
      read(line_char,*,iostat=info2) element_symbol,x_read(:,iatom),ctmp1,ctmp2
      if( info2 == 0 ) then
        if( TRIM(capitalize(ctmp1)) == 'NONE' ) then
          nucleus_wo_basis(iatom) = .TRUE.
        else
          ncenter_basis = ncenter_basis + 1
          basis_name(ncenter_basis)           = standardize_basis_name(TRIM(ctmp1))
          ecp_basis_name(ncenter_basis)       = standardize_basis_name(TRIM(ctmp1))
          auxil_basis_name(ncenter_basis)     = standardize_basis_name(TRIM(ctmp2))
          ecp_auxil_basis_name(ncenter_basis) = standardize_basis_name(TRIM(ctmp2))
        endif

      else
        read(line_char,*,iostat=info1) element_symbol,x_read(:,iatom),ctmp1
        if( info1 == 0 ) then
          if( TRIM(capitalize(ctmp1)) == 'NONE' ) then
            nucleus_wo_basis(iatom) = .TRUE.
          else
            ncenter_basis = ncenter_basis + 1
            basis_name(ncenter_basis) = standardize_basis_name(TRIM(ctmp1))
            ecp_basis_name(iatom)     = standardize_basis_name(TRIM(ctmp1))
          endif
        else
          ! when excitation is a bare projectile, assume the last atom has no basis functions
          if( excit_type%form == EXCIT_PROJECTILE .AND. iatom == natom_read) then
            nucleus_wo_basis(iatom) = .TRUE.
          else
            ncenter_basis = ncenter_basis + 1
          endif
          read(line_char,*) element_symbol,x_read(:,iatom)

        endif
      endif

      !
      ! First, try to interpret element_symbol as an integer
      read(element_symbol,*,iostat=info) zatom_read(iatom)
      ! If it fails, then assumes it is a character
      if( info /= 0 ) then
        zatom_read(iatom) = element_number(element_symbol)
      endif
    enddo

  else
    !
    ! Try to open the xyz file
    write(stdout,'(a,a)') ' Opening xyz file: ',TRIM(xyz_file)
    inquire(file=TRIM(xyz_file),exist=file_exists)
    if( .NOT. file_exists) then
      write(stdout,*) 'Tried to open the requested xyz file:',TRIM(xyz_file)
      call die('xyz file not found')
    endif
    open(newunit=xyzfile,file=TRIM(xyz_file),status='old')
    read(xyzfile,*) natom_read
    if( natom /= 0 .AND. natom+nghost+nprojectile /= natom_read ) then
      call die('the number of atoms in the input file does not correspond to the number of atoms in the xyz file')
    endif
    read(xyzfile,*)

    !natom_read read from xyz file but not from input file
    ncenter_nuclei = natom_read - nghost
    ncenter_basis  = natom_read - nprojectile
    allocate(nucleus_wo_basis(natom_read))
    nucleus_wo_basis(:) = .FALSE.

    !
    ! Need to know the number of atoms to allocate the basis arrays
    allocate(auxil_basis_name(ncenter_basis))
    allocate(small_basis_name(ncenter_basis))
    allocate(basis_name(ncenter_basis))
    allocate(ecp_basis_name(ncenter_basis))
    allocate(ecp_auxil_basis_name(ncenter_basis))
    allocate(ecp_small_basis_name(ncenter_basis))
    basis_name(:)           = standardize_basis_name(basis)
    auxil_basis_name(:)     = standardize_basis_name(auxil_basis)
    small_basis_name(:)     = standardize_basis_name(small_basis)
    ecp_basis_name(:)       = standardize_basis_name(ecp_basis)
    ecp_auxil_basis_name(:) = standardize_basis_name(ecp_auxil_basis)
    ecp_small_basis_name(:) = standardize_basis_name(ecp_small_basis)

    allocate(x_read(3,natom_read),zatom_read(natom_read))
    do iatom=1,natom_read
      read(xyzfile,*) element_symbol,x_read(:,iatom)
      !
      ! First, try to interpret element_symbol as an integer
      read(element_symbol,'(i2)',iostat=info) element
      ! If it fails, then assumes it is a character
      if( info /=0 ) then
        element = element_number(element_symbol)
      endif
      zatom_read(iatom) = element
    enddo

    close(xyzfile)



  endif

  ! Conversion to bohr if necessary
  x_read(:,:) = x_read(:,:) * length_factor
  ! vel_projectile(:) = vel_projectile(:) * length_factor
  ! For the moment, velocity in input file must be in bohrs per a.u.[time] for any length_factor

  call init_atoms(natom,nghost,nucleus_wo_basis,zatom_read,x_read,vel_projectile, &
                  (move_nuclei/='no'),excit_type%name,projectile_charge_scaling)
  deallocate(x_read,zatom_read)

  call init_ecp(ecp_elements,basis_path,ecp_type,ecp_level)
  ! If ECP are used, tweak the nuclei charges here
  zvalence(:) = zatom(:)
  do icenter=1,ncenter_nuclei
    do ielement_ecp=1,nelement_ecp
      if( ABS( element_ecp(ielement_ecp) - zatom(icenter) ) < 1.0e-5_dp ) then
        zvalence(icenter) = zatom(icenter) - REAL( ecp(ielement_ecp)%ncore , dp )
        exit
      endif
    enddo
  enddo

end subroutine setup_nuclei


!=========================================================================
end module m_inputparam
!=========================================================================
