!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to set up and store the input parameters from the input file
!
!=========================================================================
module m_inputparam
 use m_definitions
 use m_mpi
 use m_warning
 use m_atoms
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif

 !
 ! Method definitions
 integer,parameter :: perturbative = 101
 integer,parameter :: QS           = 102
 integer,parameter :: COHSEX       = 103
 integer,parameter :: QSCOHSEX     = 104
 integer,parameter :: GnW0         = 105
 integer,parameter :: GnWn         = 106
 integer,parameter :: G0W0         = 107
 integer,parameter :: GV           = 108   ! perturbative HF
 integer,parameter :: GSIGMA       = 109   ! Total energy calc
 integer,parameter :: LW           = 110   ! Luttinger-Ward log term
 integer,parameter :: GSIGMA3      = 112   ! Total energy calc
 integer,parameter :: LW2          = 113   ! Luttinger-Ward log term
 integer,parameter :: COHSEX_DEVEL = 114
 integer,parameter :: TUNED_COHSEX = 115
 integer,parameter :: G0W0_IOMEGA  = 116
 integer,parameter :: G0W0GAMMA0   = 117
 integer,parameter :: GWTILDE      = 118
 integer,parameter :: G0W0SOX0     = 119

 type calculation_type
 character(len=100) :: calc_name
 character(len=100) :: scf_name
 character(len=100) :: postscf_name
 logical            :: is_dft
 logical            :: need_exchange
 logical            :: need_exchange_lr
 logical            :: need_rpa
 logical            :: is_lr_mbpt
 logical            :: is_gw
 logical            :: is_mp2
 logical            :: is_mp2_selfenergy
 logical            :: is_ci
 logical            :: read_potential
 logical            :: is_bse,is_td
 logical            :: read_energy_qp
 integer            :: gwmethod                    ! perturbative or quasiparticle self-consistent
#ifdef HAVE_LIBXC
 type(xc_f90_pointer_t),allocatable :: xc_func(:)
 type(xc_f90_pointer_t),allocatable :: xc_info(:)
#endif
 end type calculation_type

 type(calculation_type),public    :: calc_type                   ! Sometimes we may need to tune this
 integer,protected                :: selfenergy_state_min
 integer,protected                :: selfenergy_state_max
 integer,protected                :: selfenergy_state_range
 integer,protected                :: ncoreg 
 integer,protected                :: ncorew 
 integer,protected                :: nvirtualg 
 integer,protected                :: nvirtualw 
 integer,protected                :: nvirtualspa
 logical,protected                :: is_frozencore
 logical,protected                :: is_tda,is_triplet
 logical,protected                :: is_virtual_fno
 integer,protected                :: nexcitation
 integer,protected                :: nspin
 real(dp),protected               :: spin_fact
 integer,protected                :: nscf
 integer,protected                :: mixing_first_nscf
 real(dp),protected               :: alpha_mixing
 character(len=100),protected     :: xyz_file
 character(len=100),protected     :: basis_path
 character(len=100),protected     :: basis_name
 character(len=100),protected     :: auxil_basis_name
 character(len=4),protected       :: gaussian_type
 character(len=12),protected      :: mixing_scheme
 character(len=12),protected      :: partition_scheme
 real(dp),protected               :: tolscf
 real(dp),protected               :: toldav
 real(dp),protected               :: min_overlap
 real(dp),protected               :: electrons,charge
 real(dp),protected               :: temperature
 real(dp),protected               :: magnetization
 integer,protected                :: tddft_grid_level
 integer,protected                :: grid_level
 integer,protected                :: integral_level
 logical,protected                :: has_auxil_basis
 logical,protected                :: is_full_auxil
 !
 ! the boring small complex number eta: (0.0_dp,0.001_dp) is typically over converged
 ! Having a larger ieta value smoothen the oscillation far from the HOMO-LUMO gap
 complex(dp),protected            :: ieta

 integer,protected                :: nomega_sigma
 real(dp),protected               :: step_sigma
 real(dp),protected               :: level_shifting_energy
 real(dp),protected               :: scissor
 integer,protected                :: dft_core
 integer,protected                :: npulay_hist
 integer,protected                :: scalapack_block_min
 integer,protected                :: scalapack_nprow
 integer,protected                :: scalapack_npcol
 integer,protected                :: mpi_nproc_ortho
 real(dp),protected               :: alpha_cohsex,beta_cohsex,gamma_cohsex,delta_cohsex,epsilon_cohsex

 logical,protected                :: gwgamma_tddft_
 logical,protected                :: ignore_restart_
 logical,protected                :: ignore_bigrestart_
 logical,protected                :: print_matrix_
 logical,protected                :: print_eri_
 logical,protected                :: print_wfn_
 logical,protected                :: print_w_
 logical,protected                :: print_sigma_
 logical,protected                :: print_restart_
 logical,protected                :: print_bigrestart_
 logical,protected                :: print_pdos_
 logical,protected                :: print_cube_
 real(dp),protected               :: rcut_mbpt

 real(dp),protected               :: alpha_hybrid    = 0.0_dp
 real(dp),protected               :: alpha_hybrid_lr = 0.0_dp
 real(dp),protected               :: rcut            = 0.0_dp
 real(dp),protected               :: gamma_hybrid  

!This ones are not protected, as we sometimes need to hack them
 integer,protected                 :: ndft_xc = 0
 integer,protected,allocatable     :: dft_xc_type(:)
 real(dp),protected,allocatable    :: dft_xc_coef(:)


contains


!=========================================================================
subroutine init_calculation_type(calc_type,input_key)
 implicit none
!=====
 type(calculation_type),intent(out)   :: calc_type
 character(len=100),intent(in)        :: input_key
!=====
 integer                              :: ipos
 character(len=100)                   :: key1,key2
!=====

! msg='calculation name: '//TRIM(input_key)
! call issue_warning(msg)
 !
 ! default values
 calc_type%calc_name           =  TRIM(input_key)
 calc_type%is_dft              = .FALSE.
 calc_type%need_rpa            = .FALSE.
 calc_type%is_lr_mbpt          = .FALSE.
 calc_type%is_gw               = .FALSE.
 calc_type%is_mp2              = .FALSE.
 calc_type%is_mp2_selfenergy   = .FALSE.
 calc_type%is_ci               = .FALSE.
 calc_type%is_bse              = .FALSE.
 calc_type%is_td               = .FALSE.
 calc_type%gwmethod            = 0
 calc_type%read_potential      = .FALSE.
 calc_type%postscf_name        = 'None'
 calc_type%read_energy_qp      = .FALSE.
 

 ipos=index(input_key,'+',.TRUE.)

 key1=''
 key2=''

 !
 ! If it exists, first read the last part of the calculation specifier
 if(ipos/=0) then
   key1(:ipos-1) = input_key(:ipos-1)
   key2(1:) = input_key(ipos+1:)

   calc_type%postscf_name =  TRIM(key2)

   select case(TRIM(key2))
   case('LW')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = LW
   case('LW2')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = LW2
   case('GSIGMA3')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = GSIGMA3
     calc_type%read_energy_qp = .TRUE.
   case('GSIGMA')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = GSIGMA
   case('GV')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = GV
   case('GNW0')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = GnW0
     calc_type%read_energy_qp = .TRUE.
   case('GNWN')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = GnWn
     calc_type%read_energy_qp = .TRUE.
   case('GW','G0W0')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = G0W0
   case('COHSEX')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = COHSEX
   case('COHSEX_DEVEL')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = COHSEX_DEVEL
   case('G0W0_IOMEGA')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = G0W0_IOMEGA
   case('GWTILDE')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = GWTILDE
   case('G0W0SOX0')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = G0W0SOX0
   case('G0W0GAMMA0','GWGAMMA')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = G0W0GAMMA0
   case('EVGWGAMMA','GNW0GAMMAN')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = G0W0GAMMA0
     calc_type%read_energy_qp    =.TRUE.
   case('TUNED_COHSEX')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = TUNED_COHSEX
     calc_type%is_lr_mbpt = .TRUE.
   case('LRGW')
     calc_type%is_gw      =.TRUE.
     calc_type%gwmethod   = G0W0
     calc_type%is_lr_mbpt = .TRUE.
   case('MP2')
     calc_type%is_mp2   =.TRUE.
     calc_type%gwmethod = perturbative
   case('MP2_SELFENERGY')
     calc_type%is_mp2_selfenergy =.TRUE.
     calc_type%gwmethod = perturbative
   case('EVMP2')
     calc_type%is_mp2_selfenergy =.TRUE.
     calc_type%read_energy_qp    =.TRUE.
     calc_type%gwmethod = perturbative
   case('CI')
     calc_type%is_ci =.TRUE.
   case('BSE')
     calc_type%is_bse     =.TRUE.
   case('TD')
     calc_type%is_td      =.TRUE.
   case default
     call die('Error reading keyword: postscf')
   end select
 else
   key1 = input_key
 endif

 calc_type%scf_name =  TRIM(key1)

 !
 ! Then read the first part of the calculation specifier
 select case(TRIM(key1))
 case('CI')
   calc_type%is_ci         = .TRUE.
   alpha_hybrid            = 1.00_dp
 case('H','HARTREE')
   alpha_hybrid            = 0.0_dp
 case('HF')
   alpha_hybrid            = 1.00_dp
 case('MP2')
   calc_type%is_mp2_selfenergy = .TRUE.
   calc_type%gwmethod      = QS
   alpha_hybrid            = 1.00_dp
 case('GW')
   calc_type%is_gw         = .TRUE.
   calc_type%gwmethod      = QS
   alpha_hybrid            = 1.00_dp
 case('COHSEX')
   calc_type%is_gw         = .TRUE.
   calc_type%gwmethod      = QSCOHSEX
   alpha_hybrid            = 1.00_dp
 case('VIN')
   calc_type%read_potential= .TRUE.  
 case default
   !
   ! If the calculation type is none of the above, let's assume it is DFT-type
   calc_type%is_dft=.TRUE.
   call init_dft_type(key1,calc_type)
 end select

 !
 ! Do we need Coulomb integrals?
 ! Do we need LR Coulomb integrals?
 !
 calc_type%need_exchange    = ( alpha_hybrid > 1.0e-6 )
 calc_type%need_exchange_lr = ( rcut > 1.0e-6 )


end subroutine init_calculation_type


!=========================================================================
subroutine init_dft_type(key,calc_type)
 implicit none
!=====
 character(len=*),intent(in)          :: key
 type(calculation_type),intent(inout) :: calc_type
!=====
 integer              :: idft_xc
 character(len=256)   :: string
!=====


 select case(TRIM(key))
 case('LDAX','HFPBE','PBEX','PBEHX','BX','PW91X','RPPX',&
      'BHANDH','BHANDHLYP','BHLYP','B3LYP','PBE0','HSE03','HSE06','HSE08','HCTH','CAM-B3LYP','TUNED-CAM-B3LYP','HJSX')
   ndft_xc=1
 case('LDA','SPL','VWN','VWN_RPA','PBE','PBEH','BLYP','PW91')
   ndft_xc=2
 case('RSH')
   ndft_xc=3
 case('TESTPBE0','TESTLDA0')
   ndft_xc=2
 case('TESTHSE')
   ndft_xc=3
 case default
   write(stdout,*) 'error reading calculation type'
   write(stdout,*) TRIM(key)
   call die('DFT xc is unknown')
 end select

 if( ALLOCATED(dft_xc_type) ) deallocate(dft_xc_type)
 if( ALLOCATED(dft_xc_coef) ) deallocate(dft_xc_coef)

 allocate(dft_xc_type(ndft_xc))
 allocate(dft_xc_coef(ndft_xc))
 !
 ! default is one, otherwise it is modified later
 dft_xc_coef(:) = 1.0_dp

 select case(TRIM(key))
#ifdef HAVE_LIBXC
 !
 ! LDA functionals
 case('LDAX')
   dft_xc_type(1) = XC_LDA_X
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('SPL')
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_PZ
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('LDA')
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_PW
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('VWN')
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_VWN
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('VWN_RPA')
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_VWN_RPA
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 !
 ! GGA functionals
 case('PBEX')
   dft_xc_type(1) = XC_GGA_X_PBE
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('PBE')
   dft_xc_type(1) = XC_GGA_X_PBE
   dft_xc_type(2) = XC_GGA_C_PBE
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('PBEHX')
   dft_xc_type(1) = XC_GGA_X_WPBEH
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('PBEH')
   dft_xc_type(1) = XC_GGA_X_WPBEH
   dft_xc_type(2) = XC_GGA_C_PBE
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('BX')
   dft_xc_type(1) = XC_GGA_X_B88
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('BLYP')
   dft_xc_type(1) = XC_GGA_X_B88
   dft_xc_type(2) = XC_GGA_C_LYP
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('PW91X')
   dft_xc_type(1) = XC_GGA_X_PW91
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('PW91')
   dft_xc_type(1) = XC_GGA_X_PW91
   dft_xc_type(2) = XC_GGA_C_PW91
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('HCTH')
   dft_xc_type(1) = XC_GGA_XC_HCTH_407
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('TH')
   dft_xc_type(1) = XC_GGA_XC_TH1
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('HJSX')
   dft_xc_type(1) = XC_GGA_X_HJS_PBE
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 !
 ! Meta-GGA functionals
 case('RPPX')
   dft_xc_type(1) = XC_MGGA_X_RPP09
   alpha_hybrid   = 0.00_dp
   alpha_hybrid_lr= 0.00_dp
 !
 ! Hybrid functionals
 case('HFPBE')
   dft_xc_type(1) = XC_GGA_C_PBE
   alpha_hybrid   = 1.00_dp
   alpha_hybrid_lr= 0.00_dp
 case('BHANDH')
   dft_xc_type(1) = XC_HYB_GGA_XC_BHANDH
   alpha_hybrid   = 0.50_dp
   alpha_hybrid_lr= 0.00_dp
 case('BHANDHLYP','BHLYP')
   dft_xc_type(1) = XC_HYB_GGA_XC_BHANDHLYP
   alpha_hybrid   = 0.50_dp
   alpha_hybrid_lr= 0.00_dp
 case('B3LYP')
   dft_xc_type(1) = XC_HYB_GGA_XC_B3LYP
   alpha_hybrid   = 0.20_dp
   alpha_hybrid_lr= 0.00_dp
 case('PBE0')
   dft_xc_type(1) = XC_HYB_GGA_XC_PBEH
   alpha_hybrid   = 0.25_dp
   alpha_hybrid_lr= 0.00_dp
 case('HSE03')
   dft_xc_type(1)  = XC_HYB_GGA_XC_HSE03
   alpha_hybrid    = 0.25_dp
   alpha_hybrid_lr = -alpha_hybrid
   rcut            = 1.0_dp / ( 0.15_dp / SQRT(2.0_dp) )
 case('HSE06')
   dft_xc_type(1)  = XC_HYB_GGA_XC_HSE06
   alpha_hybrid    = 0.25_dp
   alpha_hybrid_lr = -alpha_hybrid
   gamma_hybrid    = 0.11_dp
   rcut            = 1.0_dp / 0.11_dp
 case('HSE08')
   dft_xc_type(1)  = XC_HYB_GGA_XC_HJS_PBE
   alpha_hybrid    = 0.25_dp
   alpha_hybrid_lr = -alpha_hybrid
   gamma_hybrid    = 0.11_dp
   rcut            = 1.0_dp / 0.11_dp
 case('CAM-B3LYP')
   dft_xc_type(1)  = XC_HYB_GGA_XC_CAM_B3LYP
   alpha_hybrid    =  0.19_dp 
   alpha_hybrid_lr =  0.46_dp 
   rcut            =  1.0_dp / 0.33_dp  
 case('TUNED-CAM-B3LYP')
   dft_xc_type(1)  = XC_HYB_GGA_XC_TUNED_CAM_B3LYP
   alpha_hybrid    =  0.0799_dp 
   alpha_hybrid_lr =  0.9201_dp
   rcut            =  1.0_dp / 0.150_dp  
 case('RSH')
   dft_xc_type(1) = XC_GGA_X_PBE
   dft_xc_type(2) = XC_GGA_X_HJS_PBE  ! HJS is not correct in Libxc <= 2.2.2
   dft_xc_type(3) = XC_GGA_C_PBE
   dft_xc_coef(1) = 1.00_dp - (alpha_hybrid + alpha_hybrid_lr)
   dft_xc_coef(2) = alpha_hybrid_lr
   dft_xc_coef(3) = 1.00_dp
   rcut           = 1.0_dp / gamma_hybrid
 ! Testing
 case('TESTHSE')
   dft_xc_type(1) = XC_GGA_X_PBE
   dft_xc_type(2) = XC_GGA_X_HJS_PBE ! XC_GGA_X_WPBEH ! 2001  
   dft_xc_type(3) = XC_GGA_C_PBE
   alpha_hybrid   =  0.25_dp
   dft_xc_coef(1) =  1.00_dp 
   dft_xc_coef(2) = -0.25_dp
   dft_xc_coef(3) =  1.00_dp
   alpha_hybrid_lr = -alpha_hybrid
   gamma_hybrid    = 0.11_dp
   rcut           = 1.0_dp / gamma_hybrid
 case('TESTLDA0')
   alpha_hybrid   = 0.25_dp
   alpha_hybrid_lr= 0.00_dp
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_PW
   dft_xc_coef(1) =  1.00_dp - alpha_hybrid
   dft_xc_coef(2) =  1.00_dp
 case('TESTPBE0')
   alpha_hybrid   = 0.25_dp
   alpha_hybrid_lr= 0.00_dp
   dft_xc_type(1) = XC_GGA_X_PBE
   dft_xc_type(2) = XC_GGA_C_PBE
   dft_xc_coef(1) =  1.00_dp - alpha_hybrid
   dft_xc_coef(2) =  1.00_dp
#endif
 case default
   call die('Error reading keyword scf')
 end select


#ifdef HAVE_LIBXC
 !
 ! Initialize the DFT objects for LIBXC
 !
 if( ALLOCATED(calc_type%xc_func) ) then
   do idft_xc=1,SIZE(calc_type%xc_func(:))
     call xc_f90_func_end(calc_type%xc_func(idft_xc))
   enddo
   deallocate(calc_type%xc_func)
   deallocate(calc_type%xc_info)
 endif

 allocate(calc_type%xc_func(ndft_xc))
 allocate(calc_type%xc_info(ndft_xc))

 do idft_xc=1,ndft_xc

   if( nspin == 1 ) then
     call xc_f90_func_init(calc_type%xc_func(idft_xc),calc_type%xc_info(idft_xc),dft_xc_type(idft_xc),XC_UNPOLARIZED)
   else
     call xc_f90_func_init(calc_type%xc_func(idft_xc),calc_type%xc_info(idft_xc),dft_xc_type(idft_xc),XC_POLARIZED)
   endif


   if( dft_xc_type(idft_xc) < 1000 ) then
     call xc_f90_info_name(calc_type%xc_info(idft_xc),string)
     write(stdout,'(a,i4,a,a)') '   XC functional ',idft_xc,' :  ',TRIM(string)
   else
     write(stdout,'(a,i4,a,a)') '   XC functional ',idft_xc,' :  ','INTERNAL EXCHANGE-CORRELATION'
   endif

   !
   ! Tune the range for range separated hybrids
   if( dft_xc_type(idft_xc) == XC_GGA_X_HJS_PBE ) then
     call xc_f90_gga_x_hjs_set_par(calc_type%xc_func(idft_xc), gamma_hybrid )
   endif
   if( dft_xc_type(idft_xc) == XC_GGA_X_WPBEH ) then
     call xc_f90_gga_x_wpbeh_set_par(calc_type%xc_func(idft_xc),gamma_hybrid )
   endif
 enddo
#endif

end subroutine init_dft_type


!=========================================================================
subroutine summary_input(grid_quality,integral_quality)
 implicit none

 character(len=12),intent(in) :: grid_quality
 character(len=12),intent(in) :: integral_quality
!=====
 integer :: iatom,ighost
!=====

 !
 ! Summarize input parameters
 write(stdout,'(/,a,/)')    ' Summary of the input parameters '
 write(stdout,'(a25,2x,a)') ' Calculation type: ',calc_type%calc_name
 write(stdout,'(a25,2x,a)') '         SCF type: ',calc_type%scf_name
 write(stdout,'(a25,2x,a)') '    Post SCF type: ',calc_type%postscf_name
 write(stdout,'(a25,i3)')   ' Natom: ',natom
 write(stdout,'(a25,i3)')   ' Nghost:',nghost
 write(stdout,'(a25,f8.4)') ' Electrons: ',electrons
 write(stdout,'(a25,f8.4)') ' Charge: ',charge
 write(stdout,'(a25,f8.4)') ' Magnetization: ',magnetization
 write(stdout,'(a25,2x,a)') ' Basis set: ',basis_name
 write(stdout,'(a25,2x,a)') ' Auxiliary basis set: ',auxil_basis_name
 write(stdout,'(a25,2x,a)') ' Gaussian type: ',gaussian_type
 write(stdout,'(a25,2x,a)') ' Basis file path:',basis_path
 write(stdout,'(a25,i3)')   ' Spin polarization: ',nspin
 write(stdout,'(a25,i3)')   ' SCF steps: ',nscf
 write(stdout,'(a25,f8.4)') ' Mixing: ',alpha_mixing
 write(stdout,'(a25,2x,a)') ' Grid quality: ',grid_quality
 write(stdout,'(a25,2x,a)') ' Integral quality: ',integral_quality
 write(stdout,*)
 write(stdout,'(a19)')      ' IO options:'
 write(stdout,'(a30,l3)')   ' - matrices details:   ',print_matrix_       
 write(stdout,'(a30,l3)')   ' - ERI file:           ',print_eri_          
 write(stdout,'(a30,l3)')   ' - ignore big RESTART: ',ignore_bigrestart_
 write(stdout,'(a30,l3)')   ' - plot some wfns:     ',print_wfn_          
 write(stdout,'(a30,l3)')   ' - dump spectral functs',print_w_
 write(stdout,'(a30,l3)')   ' - dump self-energy    ',print_sigma_
 write(stdout,'(a30,l3)')   ' - RESTART files       ',print_restart_
 write(stdout,'(a30,l3)')   ' - big RESTART file    ',print_bigrestart_


 write(stdout,*)
 write(stdout,*) '================================'
 write(stdout,*) '      atom list'
 write(stdout,*) '                       bohr                                        angstrom'
 do iatom=1,natom
   write(stdout,'(2x,a2,3(1x,f12.6),6x,3(1x,f12.6))') element_name(zatom(iatom)),x(:,iatom),x(:,iatom)*bohr_A
 enddo
 if( nghost>0) write(stdout,'(a)') ' == ghost list'
 do ighost=1,nghost
   write(stdout,'(2x,a2,3(1x,f12.6),6x,3(1x,f12.6))') element_name(REAL(basis_element(natom+ighost),dp)),x(:,natom+ighost),x(:,natom+ighost)*bohr_A
 enddo

 write(stdout,*) '================================'
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
 use m_elements
 use m_tools,only: capitalize
 implicit none

!=====
 integer              :: ninput_argument
 character(len=100)   :: input_file_name
 integer              :: inputfile,xyzfile
 logical              :: file_exists

 character(len=100)   :: input_key
 character(len=24)    :: scf
 character(len=24)    :: postscf
 character(len=100)   :: basis
 character(len=100)   :: auxil_basis
 character(len=100)   :: default_basis_path
 character(len=12)    :: length_unit
 character(len=3)     :: ignore_restart,ignore_bigrestart,no_4center
 character(len=3)     :: print_matrix,print_eri,print_wfn,print_w,print_sigma
 character(len=3)     :: print_restart,print_bigrestart
 character(len=3)     :: print_pdos,print_cube
 character(len=3)     :: tda,triplet,frozencore,virtual_fno
 character(len=3)     :: gwgamma_tddft
 real(dp)             :: length_factor,eta
 integer              :: natom_read
 integer              :: atom_number,info,iatom
 character(len=2)     :: atom_symbol
 real(dp),allocatable :: zatom_read(:),x_read(:,:)
 real(dp)             :: beta_hybrid
 character(len=12)    :: tddft_grid_quality
 character(len=12)    :: grid_quality
 character(len=12)    :: integral_quality
 character(len=100)   :: ctmp
!=====

! Here we call the fortran code that was generated by the python script
! Any new variable should be added through the python script
#include "input_variables.f90"

! Here we call the fortran code that was generated by the python script
#include "basis_path.f90"

!=====

 ! If no value is given to basis_path in the input file,
 ! then get the default value in the source
 if( LEN(TRIM(basis_path)) == 0 ) then
   basis_path = default_basis_path
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
   open(newunit=inputfile,file=TRIM(input_file_name),status='old')
 case(0)
   inputfile = 5
   call issue_warning('Deprecated reading from stdin. Please use instead the newer syntax ./molgw inputfile > outfile')
 case default
   call die('input file name not understood')
 end select


 ! Read all the input file in one statement!
 read(inputfile,molgw)

! Here we call the fortran code that was generated by the python script
 write(stdout,'(/,1x,a,/)') ' === Echo the entire list of input variables  ==='
#include "echo_input_variables.f90"
 write(stdout,'(/,1x,a,/)') ' ================================================'


 basis_name = basis
 auxil_basis_name = auxil_basis
 has_auxil_basis = TRIM(auxil_basis) /= ''
 ieta = (0.0_dp,1.0_dp) * eta 
 alpha_hybrid_lr = beta_hybrid
 

 scf                = capitalize(scf)
 postscf            = capitalize(postscf)
 gaussian_type      = capitalize(gaussian_type)
 tddft_grid_quality = capitalize(tddft_grid_quality)
 grid_quality       = capitalize(grid_quality)
 integral_quality   = capitalize(integral_quality)
 mixing_scheme      = capitalize(mixing_scheme)
 length_unit        = capitalize(length_unit)

 ignore_restart_    = yesno(ignore_restart)
 ignore_bigrestart_ = yesno(ignore_bigrestart)
 is_full_auxil      = yesno(no_4center)
 is_tda             = yesno(tda)
 is_triplet         = yesno(triplet)
 is_frozencore      = yesno(frozencore)
 is_virtual_fno     = yesno(virtual_fno)

 print_matrix_      = yesno(print_matrix)
 print_eri_         = yesno(print_eri)
 print_wfn_         = yesno(print_wfn)
 print_w_           = yesno(print_w)
 print_sigma_       = yesno(print_sigma)
 print_restart_     = yesno(print_restart)
 print_bigrestart_  = yesno(print_bigrestart)
 print_pdos_        = yesno(print_pdos)
 print_cube_        = yesno(print_cube)
 gwgamma_tddft_     = yesno(gwgamma_tddft)

 tddft_grid_level   = interpret_quality(tddft_grid_quality)
 grid_level         = interpret_quality(grid_quality)
 integral_level     = interpret_quality(integral_quality)

 select case(TRIM(mixing_scheme))
 case('SIMPLE','PULAY')
 case default
   write(stdout,*) TRIM(mixing_scheme)
   call die('mixing scheme not recognized')
 end select

 select case(TRIM(length_unit))
 case('A','ANGSTROM')
   length_factor=1.0_dp/bohr_A
 case('BOHR','AU','A.U','A.U.')
   length_factor=1.0_dp
 case default
   call die('units for lengths in input file not understood')
 end select


 ! A few consistency checks
 if(alpha_mixing<0.0 .OR. alpha_mixing > 1.0 ) call die('alpha_mixing should be inside [0,1]')
 if(ncoreg<0) call die('negative ncoreg is meaningless')
 if(ncorew<0) call die('negative ncorew is meaningless')
 if(nvirtualg<0) call die('negative nvirtualg is meaningless')
 if(nvirtualw<0) call die('negative nvirtualw is meaningless')
 if(nvirtualspa<0) call die('negative nvirtualspa is meaningless')
 if(nvirtualg<ncoreg) call die('too small nvirtualg is meaningless')
 if(nvirtualw<ncorew) call die('too small nvirtualw is meaningless')
 if(nvirtualspa<ncorew) call die('too small nvirtualspa is meaningless')
 if(nspin/=1 .AND. nspin/=2) call die('nspin in incorrect')
 if(magnetization<-1.e-5)    call die('magnetization is negative')
 if(magnetization>1.e-5 .AND. nspin==1) call die('magnetization is non-zero and nspin is 1')
 if(nomega_sigma<0)    call die('nomega_sigma < 0')
 if(step_sigma<0.0_dp) call die('step_sigma < 0.0')

 if( is_full_auxil .AND. .NOT. has_auxil_basis) then
   write(stdout,*) 'A calculation is no 4 center integrals has been requested'
   write(stdout,*) 'However no auxiliary basis has been provided in the input file'
   call die('Please provide MOLGW with an auxiliary basis set')
 endif

 if( .NOT. has_auxil_basis .AND. nproc_world > 1 ) then
   write(stdout,*) 'Parallelization is not available without an auxiliary basis'
   call die('Please run with one CPU only or provide MOLGW with an auxiliary basis')
 endif
 if( scalapack_nprow * scalapack_npcol > nproc_world ) then
   write(stdout,'(1x,a,i4,a,i4)') 'The requested number of processors in the SCALAPACK grid: ',scalapack_nprow,' x ',scalapack_npcol
   write(stdout,'(1x,a,i5)') 'is larger than the number of total processors: ',nproc_world
   scalapack_nprow = FLOOR( SQRT( REAL(nproc_world,dp) ) )
   scalapack_npcol = nproc_world / scalapack_nprow
   write(stdout,'(1x,a,i4,a,i4)') 'Continue with a reduced SCALAPACK grid: ',scalapack_nprow,' x ',scalapack_npcol
   write(ctmp,'(a,i4,a,i4)') 'scalapack_nprow or scalapack_npcol was decreased automatically to ',scalapack_nprow,' x ',scalapack_npcol
   call issue_warning(ctmp)
 endif
 if( MODULO( nproc_world , mpi_nproc_ortho) /= 0 ) then
   write(stdout,'(1x,a,i6,a,i6)') 'mpi_nproc_ortho must be a divisor of nproc ',mpi_nproc_ortho,' / ',nproc_world
   mpi_nproc_ortho = 1
   call issue_warning('mpi_nproc_ortho value is invalid. Override it and set mpi_nproc_ortho=1')
 endif


 !
 ! Read the atom positions if no xyz file is specified
 if( LEN(TRIM(xyz_file)) == 0 ) then
    !
    ! In this case, natom must be set to a positive value
    if(natom<1) call die('natom<1')

   allocate(x_read(3,natom+nghost),zatom_read(natom+nghost))
   do iatom=1,natom+nghost
     read(inputfile,*) atom_symbol,x_read(:,iatom)
     !
     ! First, try to interpret atom_symbol as an integer
     read(atom_symbol,'(i2)',iostat=info) atom_number
     ! If it fails, then assumes it is a character
     if( info /=0 ) then
       atom_number = element_number(atom_symbol)
     endif
     zatom_read(iatom) = atom_number
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
   if( natom /= 0 .AND. natom+nghost /= natom_read ) then
     call die('the number of atoms in the input file does not correspond to the number of atoms in the xyz file')
   endif
   read(xyzfile,*) 

   natom = natom_read-nghost

   allocate(x_read(3,natom+nghost),zatom_read(natom+nghost))
   do iatom=1,natom+nghost
     read(xyzfile,*) atom_symbol,x_read(:,iatom)
     !
     ! First, try to interpret atom_symbol as an integer
     read(atom_symbol,'(i2)',iostat=info) atom_number
     ! If it fails, then assumes it is a character
     if( info /=0 ) then
       atom_number = element_number(atom_symbol)
     endif
     zatom_read(iatom) = atom_number
   enddo

   close(xyzfile)
   

   if( ABS(length_factor - 1.0_dp ) < 1.0e-6_dp  ) then
     write(stdout,*) 'xyz files are always in Angstrom. However, length_unit was set to bohr'
     call die('Please set length_unit to Angstrom')
   endif

 endif

 x_read(:,:) = x_read(:,:) * length_factor
 call init_atoms(natom,nghost,zatom_read,x_read)
 deallocate(x_read,zatom_read)

 !
 ! Interpret the scf and postscf input parameters
 if( TRIM(postscf) =='' ) then
   input_key=TRIM(scf)
 else
   input_key=TRIM(scf)//'+'//TRIM(postscf)
 endif
 call init_calculation_type(calc_type,input_key)

 ! Some additional checks
 if( nexcitation /=0 .AND. calc_type%is_gw ) then
   call die('Davidson diago is not compatible with GW. Set nexcitation to 0')
 endif

 spin_fact = REAL(-nspin+3,dp)
 electrons = SUM(zatom(:)) - charge


 ! Echo the interpreted input variables
 call summary_input(grid_quality,integral_quality)


contains


function interpret_quality(quality) result(quality_level)
 implicit none
 character(len=12),intent(inout) :: quality
 integer                         :: quality_level
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


function yesno(char3)
 implicit none
 character(len=3),intent(inout) :: char3
 logical                        :: yesno
!=====
 
 char3 = capitalize(char3)

 select case(TRIM(char3))
 case('YES','Y')
   yesno=.TRUE.
 case('NO','N')
   yesno=.FALSE.
 case default
  call die('Yes or No, I cannot interpret this input')
 end select
 
end function yesno

end subroutine read_inputfile_namelist


!=========================================================================
end module m_inputparam
