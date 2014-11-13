!=========================================================================
#include "macros.h"
!=========================================================================
module m_calculation_type
 use m_definitions
 use m_mpi
 use m_warning
#ifdef HAVE_LIBXC
 use libxc_funcs_m
#endif

 !
 ! Method definitions
 integer,parameter :: perturbative = 101
 integer,parameter :: QS           = 102
 integer,parameter :: COHSEX       = 103
 integer,parameter :: QSCOHSEX     = 104


 real(dp)          :: alpha_hybrid    = 1.0_dp
 real(dp)          :: alpha_hybrid_lr = 0.0_dp
 real(dp)          :: rcut            = 0.0_dp

 integer                   :: ndft_xc      = 0
 integer,pointer           :: dft_xc_type(:)
 real(dp),pointer          :: dft_xc_coef(:)

 type calculation_type
   integer :: type
   logical :: need_exchange
   logical :: need_final_exchange
   logical :: need_rpa
   logical :: is_lr_mbpt
   logical :: is_screened_hybrid
   logical :: is_gw
   logical :: is_mp2
   logical :: is_ci
   logical :: read_potential
   logical :: is_bse,is_td
   integer :: method                    ! perturbative or quasiparticle self-consistent
 end type calculation_type

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
 calc_type%need_exchange       = .FALSE.
 calc_type%need_final_exchange = .FALSE.
 calc_type%need_rpa            = .FALSE.
 calc_type%is_lr_mbpt          = .FALSE.
 calc_type%is_screened_hybrid  = .FALSE.
 calc_type%is_gw               = .FALSE.
 calc_type%is_mp2              = .FALSE.
 calc_type%is_ci               = .FALSE.
 calc_type%is_bse              = .FALSE.
 calc_type%is_td               = .FALSE.
 calc_type%method              = 0
 calc_type%read_potential      = .FALSE.

 ipos=index(input_key,'+',.TRUE.)

 key1=''
 key2=''

 !
 ! If it exists, first read the last part of the calculation specifier
 if(ipos/=0) then
   key1(:ipos-1) = input_key(:ipos-1)
   key2(1:) = input_key(ipos+1:)
   select case(TRIM(key2))
   case('GW')
     calc_type%is_gw  =.TRUE.
     calc_type%method = perturbative
   case('COHSEX')
     calc_type%is_gw  =.TRUE.
     calc_type%method = COHSEX
   case('LRGW')
     calc_type%is_gw  =.TRUE.
     calc_type%method = perturbative
     calc_type%is_lr_mbpt        = .TRUE.
   case('MP2')
     calc_type%is_mp2 =.TRUE.
     calc_type%method = perturbative
   case('CI')
     calc_type%is_ci =.TRUE.
   case('BSE')
     calc_type%is_bse =.TRUE.
   case('TD')
     calc_type%is_td =.TRUE.
   case default
     stop'error reading calculation type part 2'
   end select
 else
   key1 = input_key
 endif

 !
 ! Then read the first part of the calculation specifier
 select case(TRIM(key1))
 case('CI')
   calc_type%is_ci =.TRUE.
   calc_type%need_exchange = .TRUE.
 case('HF')
   calc_type%need_exchange = .TRUE.  
 case('MP2')
   calc_type%need_exchange = .TRUE.  
   calc_type%is_mp2        = .TRUE.
   calc_type%method        = QS
 case('GW')
   calc_type%need_exchange = .TRUE.  
   calc_type%is_gw         = .TRUE.
   calc_type%method        = QS
 case('VIN')
   calc_type%read_potential= .TRUE.  
 case default
   !
   ! If the calculation type is none of the above, let's assume it is DFT-type
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

 alpha_hybrid = 0.0_dp
 if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.

 select case(TRIM(key))
 case('LDAx','PBEx','PBEhx','Bx','PW91x','BJx','RPPx',&
      'BHANDH','BHANDHLYP','B3LYP','PBE0','HSE03','HSE06','HSE08','HCTH','CAM-B3LYP','TD-CAM-B3LYP')
   ndft_xc=1
 case('LDA','VWN','VWN_RPA','PBE','PBEh','BLYP','PW91')
   ndft_xc=2
 case('TESTHSE')
   ndft_xc=1
 case('TESTPBE0')
   ndft_xc=2
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
 case('TD-CAM-B3LYP')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.
   dft_xc_type(1)  = XC_HYB_GGA_XC_TUNED_CAM_B3LYP
   alpha_hybrid    =  0.0799_dp 
   alpha_hybrid_lr =  0.9201_dp
   rcut            =  1.0_dp / 0.150_dp  
 ! Testing
 case('TESTHSE')
   calc_type%is_screened_hybrid  = .TRUE.
   calc_type%need_exchange       = .TRUE.  
   dft_xc_type(1) = XC_GGA_X_wPBEh
   dft_xc_coef(2) = -0.25_dp
   alpha_hybrid   = 0.25_dp
   rcut           = 1.0_dp / 0.11_dp
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

  WRITE_MASTER(*,*) 'type                        ',calc_type%type
  WRITE_MASTER(*,*) 'need_exchange               ',calc_type%need_exchange
  WRITE_MASTER(*,*) 'need_final_exchange         ',calc_type%need_final_exchange
  WRITE_MASTER(*,*) 'need_rpa                    ',calc_type%need_rpa
  WRITE_MASTER(*,*) 'is_gw                       ',calc_type%is_gw
  WRITE_MASTER(*,*) 'is_mp2                      ',calc_type%is_mp2
  WRITE_MASTER(*,*) 'is_ci                       ',calc_type%is_ci
  WRITE_MASTER(*,*) 'method                      ',calc_type%method  

end subroutine output_calculation_type

end module m_calculation_type
