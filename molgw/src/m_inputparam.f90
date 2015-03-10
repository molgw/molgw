!=========================================================================
#include "macros.h"
!=========================================================================
module m_inputparam
 use m_definitions
 use m_mpi
 use m_atoms
 use m_basis_set
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
 use iso_c_binding,only: C_INT,C_DOUBLE

 !
 ! Method definitions
 integer,parameter :: perturbative = 101
 integer,parameter :: QS           = 102
 integer,parameter :: COHSEX       = 103
 integer,parameter :: QSCOHSEX     = 104
 integer,parameter :: GnW0         = 105
 integer,parameter :: GnWn         = 106
 integer,parameter :: G0W0         = 107

 type calculation_type
   character(len=100) :: calc_name
   character(len=100) :: scf_name
   character(len=100) :: postscf_name
   logical            :: is_dft
   logical            :: need_exchange
   logical            :: need_rpa
   logical            :: is_lr_mbpt
   logical            :: is_screened_hybrid
   logical            :: is_gw
   logical            :: is_mp2
   logical            :: is_ci
   logical            :: read_potential
   logical            :: is_bse,is_td
   integer            :: gwmethod                    ! perturbative or quasiparticle self-consistent
 end type calculation_type

 integer,protected                :: ncoreg 
 integer,protected                :: ncorew 
 integer,protected                :: nvirtualg 
 integer,protected                :: nvirtualw 
 logical,protected                :: is_frozencore
 logical,protected                :: is_tda,is_triplet
 integer,protected                :: nspin
 real(dp),protected               :: spin_fact
 integer,protected                :: nscf
 real(dp),protected               :: alpha_mixing
 character(len=100),protected     :: basis_path
 character(len=100),protected     :: basis_name
 character(len=100),protected     :: auxil_basis_name
 character(len=4),protected       :: gaussian_type
 character(len=12),protected      :: mixing_scheme
 real(dp),protected               :: tolscf
 real(dp),protected               :: electrons,charge
 real(dp),protected               :: magnetization
 type(calculation_type),protected :: calc_type
 integer,protected                :: grid_level
 integer,protected                :: integral_level
 logical,protected                :: has_auxil_basis
 logical,protected                :: is_full_auxil
 real(dp),protected               :: pole_eta
 integer,protected                :: nomega_sigma
 real(dp),protected               :: step_sigma

 logical,protected                :: no_restart   
 logical,protected                :: ignore_big_restart
 logical,protected                :: print_matrix_
 logical,protected                :: print_eri_
 logical,protected                :: print_wfn_
 logical,protected                :: print_w_
 logical,protected                :: print_sigma_

 real(dp),protected               :: alpha_hybrid    = 0.0_dp
 real(dp),protected               :: alpha_hybrid_lr = 0.0_dp
 real(dp),protected               :: rcut            = 0.0_dp
 real(dp),protected               :: gamma_hybrid  

 integer,protected                    :: ndft_xc      = 0
 integer(C_INT),protected,allocatable :: dft_xc_type(:)
 real(C_DOUBLE),protected,allocatable :: dft_xc_coef(:)


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

 msg='calculation name: '//TRIM(input_key)
 call issue_warning(msg)
 !
 ! default values
 calc_type%calc_name           =  TRIM(input_key)
 calc_type%is_dft              = .FALSE.
 calc_type%need_exchange       = .FALSE.
 calc_type%need_rpa            = .FALSE.
 calc_type%is_lr_mbpt          = .FALSE.
 calc_type%is_screened_hybrid  = .FALSE.
 calc_type%is_gw               = .FALSE.
 calc_type%is_mp2              = .FALSE.
 calc_type%is_ci               = .FALSE.
 calc_type%is_bse              = .FALSE.
 calc_type%is_td               = .FALSE.
 calc_type%gwmethod            = 0
 calc_type%read_potential      = .FALSE.
 calc_type%postscf_name        = 'None'
 

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
   case('GnW0')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = GnW0
   case('GnWn')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = GnWn
   case('GW','G0W0')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = G0W0
   case('COHSEX')
     calc_type%is_gw    =.TRUE.
     calc_type%gwmethod = COHSEX
   case('LRGW')
     calc_type%is_gw      =.TRUE.
     calc_type%gwmethod   = G0W0
     calc_type%is_lr_mbpt = .TRUE.
   case('MP2')
     calc_type%is_mp2   =.TRUE.
     calc_type%gwmethod = perturbative
   case('CI')
     calc_type%is_ci =.TRUE.

   case('BSE')
     calc_type%is_bse     =.TRUE.

   case('TD')
     calc_type%is_td      =.TRUE.

   case default
     stop'error reading calculation type part 2'
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
   calc_type%need_exchange = .TRUE.
   alpha_hybrid            = 1.00_dp
 case('HF')
   calc_type%need_exchange = .TRUE.  
   alpha_hybrid            = 1.00_dp
 case('MP2')
   calc_type%need_exchange = .TRUE.  
   calc_type%is_mp2        = .TRUE.
   calc_type%gwmethod      = QS
   alpha_hybrid            = 1.00_dp
 case('GW')
   calc_type%need_exchange = .TRUE.  
   calc_type%is_gw         = .TRUE.
   calc_type%gwmethod      = QS
   alpha_hybrid            = 1.00_dp
 case('COHSEX')
   calc_type%need_exchange = .TRUE.  
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

end subroutine init_calculation_type


!=========================================================================
subroutine init_dft_type(key,calc_type)
 implicit none
!=====
 character(len=100),intent(in)          :: key
 type(calculation_type),intent(inout)   :: calc_type
!=====


 select case(TRIM(key))
 case('LDAx','PBEx','PBEhx','Bx','PW91x','BJx','RPPx',&
      'BHANDH','BHANDHLYP','B3LYP','PBE0','HSE03','HSE06','HSE08','HCTH','CAM-B3LYP','TUNED-CAM-B3LYP')
   ndft_xc=1
 case('LDA','SPL','VWN','VWN_RPA','PBE','PBEh','BLYP','PW91')
   ndft_xc=2
 case('RSH')
   ndft_xc=3
 case('TESTPBE0','TESTLDA0')
   ndft_xc=2
 case('TESTHSE')
   ndft_xc=3
 case default
   WRITE_MASTER(*,*) 'error reading calculation type'
   WRITE_MASTER(*,*) TRIM(key)
   stop' is unknown'
 end select

 allocate(dft_xc_type(ndft_xc))
 allocate(dft_xc_coef(ndft_xc))
 !
 ! default is one, otherwise it is modified later
 dft_xc_coef(:) = 1.0_dp

 select case(TRIM(key))
#ifdef HAVE_LIBXC
 !
 ! LDA functionals
 case('LDAx')
   dft_xc_type(1) = XC_LDA_X
 case('SPL')
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_PZ
 case('LDA')
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_PW
 case('VWN')
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_VWN
 case('VWN_RPA')
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_VWN_RPA
 !
 ! GGA functionals
 case('PBEx')
   dft_xc_type(1) = XC_GGA_X_PBE
 case('PBE')
   dft_xc_type(1) = XC_GGA_X_PBE
   dft_xc_type(2) = XC_GGA_C_PBE
 case('PBEhx')
   dft_xc_type(1) = XC_GGA_X_WPBEH
 case('PBEh')
   dft_xc_type(1) = XC_GGA_X_WPBEH
   dft_xc_type(2) = XC_GGA_C_PBE
 case('Bx')
   dft_xc_type(1) = XC_GGA_X_B88
 case('BLYP')
   dft_xc_type(1) = XC_GGA_X_B88
   dft_xc_type(2) = XC_GGA_C_LYP
 case('PW91x')
   dft_xc_type(1) = XC_GGA_X_PW91
 case('PW91')
   dft_xc_type(1) = XC_GGA_X_PW91
   dft_xc_type(2) = XC_GGA_C_PW91
 case('HCTH')
   dft_xc_type(1) = XC_GGA_XC_HCTH_407
 case('TH')
   dft_xc_type(1) = XC_GGA_XC_TH1
 !
 ! Meta-GGA functionals
 case('BJx')
   dft_xc_type(1) = XC_MGGA_X_BJ06
 case('RPPx')
   dft_xc_type(1) = XC_MGGA_X_RPP09
 !
 ! Hybrid functionals
 case('BHANDH')
   calc_type%need_exchange = .TRUE.  
   dft_xc_type(1) = XC_HYB_GGA_XC_BHANDH
   alpha_hybrid = 0.50_dp
 case('BHANDHLYP')
   calc_type%need_exchange = .TRUE.  
   dft_xc_type(1) = XC_HYB_GGA_XC_BHANDHLYP
   alpha_hybrid = 0.50_dp
 case('B3LYP')
   calc_type%need_exchange = .TRUE.  
   dft_xc_type(1) = XC_HYB_GGA_XC_B3LYP
   alpha_hybrid = 0.20_dp
 case('PBE0')
   calc_type%need_exchange = .TRUE.  
   dft_xc_type(1) = XC_HYB_GGA_XC_PBEH
   alpha_hybrid = 0.25_dp
 case('HSE03')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.  
   dft_xc_type(1) = XC_HYB_GGA_XC_HSE03
   alpha_hybrid    = 0.25_dp
   alpha_hybrid_lr = -alpha_hybrid
   rcut         = 1.0_dp / ( 0.15_dp / SQRT(2.0_dp) )
 case('HSE06')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.  
   dft_xc_type(1) = XC_HYB_GGA_XC_HSE06
   alpha_hybrid    = 0.25_dp
   alpha_hybrid_lr = -alpha_hybrid
   rcut         = 1.0_dp / 0.11_dp
 case('HSE08')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.
   dft_xc_type(1) = XC_HYB_GGA_XC_HJS_PBE
   alpha_hybrid    = 0.25_dp
   alpha_hybrid_lr = -alpha_hybrid
   rcut           = 1.0_dp / 0.11_dp
 case('CAM-B3LYP')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.
   dft_xc_type(1)  = XC_HYB_GGA_XC_CAM_B3LYP
   alpha_hybrid    =  0.19_dp 
   alpha_hybrid_lr =  0.46_dp 
   rcut            =  1.0_dp / 0.33_dp  
 case('TUNED-CAM-B3LYP')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.
   dft_xc_type(1)  = XC_HYB_GGA_XC_TUNED_CAM_B3LYP
   alpha_hybrid    =  0.0799_dp 
   alpha_hybrid_lr =  0.9201_dp
   rcut            =  1.0_dp / 0.150_dp  
 case('RSH')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.  
   dft_xc_type(1) = XC_GGA_X_PBE
   dft_xc_type(2) = XC_GGA_X_HJS_PBE 
   dft_xc_type(3) = XC_GGA_C_PBE
   dft_xc_coef(1) =  1.00_dp - (alpha_hybrid + alpha_hybrid_lr)
   dft_xc_coef(2) =  -alpha_hybrid
   dft_xc_coef(3) =  1.00_dp
   rcut           = 1.0_dp / gamma_hybrid
 ! Testing
 case('TESTHSE')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.  
   dft_xc_type(1) = XC_GGA_X_PBE
   dft_xc_type(2) = XC_GGA_X_HJS_PBE 
   dft_xc_type(3) = XC_GGA_C_PBE
   alpha_hybrid   =  0.25_dp
   dft_xc_coef(1) =  1.00_dp 
   dft_xc_coef(2) = -0.25_dp
   dft_xc_coef(3) =  1.00_dp
   alpha_hybrid_lr = -alpha_hybrid
   gamma_hybrid    = 0.11_dp
   rcut           = 1.0_dp / gamma_hybrid
 case('TESTLDA0')
   calc_type%need_exchange       = .TRUE.
   alpha_hybrid   = 0.25_dp
   dft_xc_type(1) = XC_LDA_X
   dft_xc_type(2) = XC_LDA_C_PW
   dft_xc_coef(1) =  1.00_dp - alpha_hybrid
   dft_xc_coef(2) =  1.00_dp
 case('TESTPBE0')
   calc_type%need_exchange       = .TRUE.
   alpha_hybrid   = 0.25_dp
   dft_xc_type(1) = XC_GGA_X_PBE
   dft_xc_type(2) = XC_GGA_C_PBE
   dft_xc_coef(1) =  1.00_dp - alpha_hybrid
   dft_xc_coef(2) =  1.00_dp
#endif
 case default
   stop'error reading calculation type part 1'
 end select


end subroutine init_dft_type


!=========================================================================
subroutine output_calculation_type(calc_type)
 implicit none
 type(calculation_type),intent(in)   :: calc_type
!=====

  WRITE_MASTER(*,*) 'need_exchange               ',calc_type%need_exchange
  WRITE_MASTER(*,*) 'need_rpa                    ',calc_type%need_rpa
  WRITE_MASTER(*,*) 'is_gw                       ',calc_type%is_gw
  WRITE_MASTER(*,*) 'is_mp2                      ',calc_type%is_mp2
  WRITE_MASTER(*,*) 'is_ci                       ',calc_type%is_ci
  WRITE_MASTER(*,*) 'is_td                       ',calc_type%is_td
  WRITE_MASTER(*,*) 'is_bse                      ',calc_type%is_bse
  WRITE_MASTER(*,*) 'method                      ',calc_type%gwmethod  

end subroutine output_calculation_type


!=========================================================================
subroutine summary_input(grid_quality,integral_quality)
 implicit none

 character(len=12),intent(in) :: grid_quality
 character(len=12),intent(in) :: integral_quality
!=====
 integer :: iatom
!=====

 !
 ! Summarize input parameters
 WRITE_MASTER(*,'(/,a,/)')    ' Summary of the input parameters '
 WRITE_MASTER(*,'(a25,2x,a)') ' Calculation type: ',calc_type%calc_name
 WRITE_MASTER(*,'(a25,2x,a)') '         SCF type: ',calc_type%scf_name
 WRITE_MASTER(*,'(a25,2x,a)') '    Post SCF type: ',calc_type%postscf_name
 WRITE_MASTER(*,'(a25,i3)')   ' Natom: ',natom
 WRITE_MASTER(*,'(a25,f8.4)') ' Electrons: ',electrons
 WRITE_MASTER(*,'(a25,f8.4)') ' Charge: ',charge
 WRITE_MASTER(*,'(a25,f8.4)') ' Magnetization: ',magnetization
 WRITE_MASTER(*,'(a25,2x,a)') ' Basis set: ',basis_name
 WRITE_MASTER(*,'(a25,2x,a)') ' Auxiliary basis set: ',auxil_basis_name
 WRITE_MASTER(*,'(a25,2x,a)') ' Gaussian type: ',gaussian_type
 WRITE_MASTER(*,'(a25,i3)')   ' Spin polarization: ',nspin
 WRITE_MASTER(*,'(a25,i3)')   ' SCF steps: ',nscf
 WRITE_MASTER(*,'(a25,f8.4)') ' Mixing: ',alpha_mixing
 WRITE_MASTER(*,'(a25,2x,a)') ' Grid quality: ',grid_quality
 WRITE_MASTER(*,'(a25,2x,a)') ' Integral quality: ',integral_quality
 WRITE_MASTER(*,*)
 WRITE_MASTER(*,'(a19)')      ' IO options:'
 WRITE_MASTER(*,'(a30,l3)')   ' - matrices details:   ',print_matrix_       
 WRITE_MASTER(*,'(a30,l3)')   ' - ERI file:           ',print_eri_          
 WRITE_MASTER(*,'(a30,l3)')   ' - ignore big RESTART: ',ignore_big_restart
 WRITE_MASTER(*,'(a30,l3)')   ' - plot some wfns:     ',print_wfn_          
 WRITE_MASTER(*,'(a30,l3)')   ' - dump spectral functs',print_w_
 WRITE_MASTER(*,'(a30,l3)')   ' - dump self-energy    ',print_sigma_


 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) '================================'
 WRITE_MASTER(*,*) '      atom list'
 WRITE_MASTER(*,*) '                       bohr                                        angstrom'
 do iatom=1,natom
   WRITE_MASTER(*,'(2x,a2,3(x,f12.6),6x,3(x,f12.6))') element_name(zatom(iatom)),x(:,iatom),x(:,iatom)*bohr_A
 enddo

 WRITE_MASTER(*,*) '================================'
 WRITE_MASTER(*,'(a,i5)') ' Number of bonds ',nbond
 if(inversion) then
   WRITE_MASTER(*,*) 'Molecule has inversion symmetry'
 else
   WRITE_MASTER(*,*) 'Molecule does not have inversion symmetry'
 endif
 if(linear) then
   WRITE_MASTER(*,*) 'Molecule is linear'
 else
   WRITE_MASTER(*,*) 'Molecule is not linear'
   if(planar) then
     WRITE_MASTER(*,*) 'Molecule is planar'
   else
     WRITE_MASTER(*,*) 'Molecule is not planar'
   endif
 endif
 WRITE_MASTER(*,*)


end subroutine summary_input


!=========================================================================
subroutine read_inputfile_namelist()
 use m_elements
 use m_tools,only: capitalize
 implicit none

!=====
 character(len=100)   :: input_key
 character(len=12)    :: scf
 character(len=12)    :: postscf
 character(len=100)   :: basis
 character(len=100)   :: auxil_basis
 character(len=12)    :: length_unit
 character(len=3)     :: ignore_restart,ignore_bigrestart,no_4center
 character(len=3)     :: print_matrix,print_eri,print_wfn,print_w,print_sigma
 character(len=3)     :: tda,triplet,frozencore
 real(dp)             :: length_factor,eta
 integer              :: atom_number,info,iatom
 character(len=2)     :: atom_symbol
 real(dp),allocatable :: zatom_read(:),x_read(:,:)
 real(dp)             :: beta_hybrid
 character(len=12)    :: grid_quality
 character(len=12)    :: integral_quality

 namelist /molgw/ scf,postscf,                                                                          &
                  alpha_hybrid,beta_hybrid,gamma_hybrid,                                                &
                  basis,auxil_basis,basis_path,gaussian_type,no_4center,                                &
                  nspin,charge,magnetization,                                                           &
                  grid_quality,integral_quality,                                                        &
                  nscf,alpha_mixing,mixing_scheme,tolscf,                                               &
                  tda,triplet,eta,frozencore,ncoreg,ncorew,nvirtualg,nvirtualw,nomega_sigma,step_sigma, &
                  ignore_restart,print_matrix,print_eri,print_wfn,print_w,print_sigma,                  &
                  length_unit,natom
!=====


 scf              = ''
 postscf          = ''

 alpha_hybrid     = 0.0_dp
 beta_hybrid      = 0.0_dp
 gamma_hybrid     = HUGE(1.0_dp)

 basis            = ''
 auxil_basis      = ''
 basis_path       = '.'
 gaussian_type    = 'PURE'
 no_4center       = 'NO'

 charge           = 0.0_dp
 nspin            = 1

 grid_quality     = 'HIGH'
 integral_quality = 'HIGH'

 nscf             = 30
 alpha_mixing     = 0.70_dp
 mixing_scheme    = 'PULAY'
 tolscf           = 1.0e-7_dp

 tda               = 'NO'
 eta               = 1.0e-3_dp
 triplet           = 'NO'
 frozencore        = 'NO'
 ncoreg            = 0
 ncorew            = 0
 nvirtualg         = HUGE(1)
 nvirtualw         = HUGE(1)
 nomega_sigma      = 51
 step_sigma        = 0.01_dp

 ignore_restart    = 'NO'
 ignore_bigrestart = 'NO'
 print_matrix      = 'NO'
 print_eri         = 'NO'
 print_wfn         = 'NO'
 print_w           = 'NO'
 print_sigma       = 'NO'

 natom             = 0
 length_unit       = 'ANGSTROM'


 ! Read all the input file in one statement!
 read(*,molgw)

 basis_name = basis
 auxil_basis_name = auxil_basis
 has_auxil_basis = TRIM(auxil_basis) /= ''
 pole_eta = eta
 alpha_hybrid_lr = beta_hybrid
 

 scf                = capitalize(scf)
 postscf            = capitalize(postscf)
 gaussian_type      = capitalize(gaussian_type)
 grid_quality       = capitalize(grid_quality)
 integral_quality   = capitalize(integral_quality)
 mixing_scheme      = capitalize(mixing_scheme)
 length_unit        = capitalize(length_unit)

 no_restart         = yesno(ignore_restart)
 ignore_big_restart = yesno(ignore_bigrestart)
 is_full_auxil      = yesno(no_4center)
 is_tda             = yesno(tda)
 is_triplet         = yesno(triplet)
 is_frozencore      = yesno(frozencore)

 print_matrix_      = yesno(print_matrix)
 print_eri_         = yesno(print_eri)
 print_wfn_         = yesno(print_wfn)
 print_w_           = yesno(print_w)
 print_sigma_       = yesno(print_sigma)

 grid_level     = interpret_quality(grid_quality)
 integral_level = interpret_quality(integral_quality)

 select case(TRIM(mixing_scheme))
 case('SIMPLE','PULAY')
 case default
   WRITE_MASTER(*,*) TRIM(mixing_scheme)
   stop'mixing scheme not recognized'
 end select

 select case(TRIM(length_unit))
 case('A','ANGSTROM')
   length_factor=1.0_dp/bohr_A
 case('BOHR','AU','A.U','A.U.')
   length_factor=1.0_dp
 case default
   stop'units for lengths in input file not understood'
 end select


 ! A few consistency checks
 if(natom<1) stop'natom<1'
 if(alpha_mixing<0.0 .OR. alpha_mixing > 1.0 ) stop'alpha_mixing should be inside [0,1]'
 if(ncoreg<0) stop'negative ncoreg is meaningless'
 if(ncorew<0) stop'negative ncorew is meaningless'
 if(nvirtualg<0) stop'negative nvirtualg is meaningless'
 if(nvirtualw<0) stop'negative nvirtualw is meaningless'
 if(nvirtualg<ncoreg) stop'too small nvirtualg is meaningless'
 if(nvirtualw<ncorew) stop'too small nvirtualw is meaningless'
 if(nspin/=1 .AND. nspin/=2) stop'nspin in incorrect'
 if(magnetization<-1.e-5)    stop'magnetization is negative'
 if(magnetization>1.e-5 .AND. nspin==1) stop'magnetization is non-zero and nspin is 1'
 if(nomega_sigma<0)    stop'nomega_sigma < 0'
 if(step_sigma<0.0_dp) stop'step_sigma < 0.0'

 if( is_full_auxil .AND. .NOT. has_auxil_basis) then
   WRITE_MASTER(*,*) 'A calculation is no 4 center integrals has been requested'
   WRITE_MASTER(*,*) 'However no auxiliary basis has been provided in the input file'
   stop'STOP HERE'
 endif


 !
 ! Read the atom positions
 allocate(x_read(3,natom),zatom_read(natom))
 do iatom=1,natom
   read(*,*) atom_symbol,x_read(:,iatom)
   !
   ! First, try to interpret atom_symbol as an integer
   read(atom_symbol,'(i2)',iostat=info) atom_number
   ! If it fails, then assumes it is a character
   if( info /=0 ) then
     atom_number = element_number(atom_symbol)
   endif
   zatom_read(iatom) = atom_number
 enddo
 x_read(:,:) = x_read(:,:) * length_factor
 call init_atoms(natom,zatom_read,x_read)
 deallocate(x_read,zatom_read)


 !
 ! Interpret the scf and postscf input parameters
 if( TRIM(postscf) =='' ) then
   input_key=TRIM(scf)
 else
   input_key=TRIM(scf)//'+'//TRIM(postscf)
 endif
 call init_calculation_type(calc_type,input_key)


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
   quality_level = 10
 case('MEDIUM','MED','M')
   quality_level = 20
 case('HIGH','HI','H')
   quality_level = 30
 case('VERY HIGH','VERYHIGH','VH')
   quality_level = 40
 case('INSANE','I')
   quality_level = 50
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
  stop'Yes or No, I cannot interpret this input'
 end select
 
end function yesno

end subroutine read_inputfile_namelist


!=========================================================================
end module m_inputparam
