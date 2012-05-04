module m_calculation_type
 use m_warning
#ifdef HAVE_LIBXC
 use libxc_funcs_m
#endif

 integer,parameter :: HF=1
 integer,parameter :: GW=2
 integer,parameter :: LDA=3
 integer,parameter :: CI=4
 integer,parameter :: MP2=5
 integer,parameter :: HARTREE=6
 integer,parameter :: RPA=7
 integer,parameter :: G0W0=8
 integer,parameter :: COHSEX=9
 integer,parameter :: HYBRID=10
 integer,parameter :: LRGW=12
! Method definitions
 integer,parameter :: perturbative=101
 integer,parameter :: QS          =102

 real(dp)          :: alpha_hybrid = 1.0_dp

 type calculation_type
   integer :: type
   logical :: need_exchange
   logical :: need_final_exchange
   logical :: need_dft_xc
   logical :: need_rpa
   logical :: need_lr_integrals
   logical :: is_gw
   logical :: is_mp2
   logical :: is_ci
   integer :: method                    ! perturbative or quasiparticle self-consistent
   integer :: dft_x
   integer :: dft_c
 end type calculation_type

contains

subroutine init_calculation_type(calc_type,input_key)
 implicit none
!=====
 type(calculation_type),intent(out)   :: calc_type
 character(len=100),intent(in)             :: input_key
!=====
 integer                       :: ipos
 character(len=100)            :: key1,key2
!=====

 msg='calculation name: '//TRIM(input_key)
 call issue_warning(msg)
 !
 ! default values
 calc_type%need_exchange       = .FALSE.
 calc_type%need_final_exchange = .FALSE.
 calc_type%need_dft_xc         = .FALSE.
 calc_type%need_rpa            = .FALSE.
 calc_type%need_lr_integrals   = .FALSE.
 calc_type%is_gw               = .FALSE.
 calc_type%is_mp2              = .FALSE.
 calc_type%is_ci               = .FALSE.
 calc_type%method              = 0
 calc_type%dft_x               = 0
 calc_type%dft_c               = 0

 ipos=index(input_key,'+',.TRUE.)

 key1=''
 key2=''
 if(ipos/=0) then
   key1(:ipos-1) = input_key(:ipos-1)
   key2(1:) = input_key(ipos+1:)
   select case(TRIM(key2))
   case('GW')
     calc_type%is_gw  =.TRUE.
     calc_type%method = perturbative
   case('LRGW')
     calc_type%is_gw  =.TRUE.
     calc_type%method = perturbative
     calc_type%need_lr_integrals = .TRUE.
   case('MP2')
     calc_type%is_mp2 =.TRUE.
     calc_type%method = perturbative
   case('CI')
     calc_type%is_ci =.TRUE.
   case default
     stop'error reading calculation type part 2'
   end select
 else
   key1 = input_key
 endif


 select case(TRIM(key1))
 case('CI')
   calc_type%type          = CI
   calc_type%need_exchange = .TRUE.
 case('HF')
   calc_type%type          = HF
   calc_type%need_exchange = .TRUE.  
 case('MP2')
   calc_type%type          = MP2
   calc_type%need_exchange = .TRUE.  
   calc_type%is_mp2        = .TRUE.
   calc_type%method        = QS
 case('GW')
   calc_type%type          = GW
   calc_type%need_exchange = .TRUE.  
   calc_type%is_gw         = .TRUE.
   calc_type%method        = QS
#ifdef HAVE_LIBXC
 !
 ! LDA functionals
 case('LDA')
   calc_type%type          = LDA
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = XC_LDA_X
   calc_type%dft_c = XC_LDA_C_PW
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.0_dp
 case('VWN')
   calc_type%type          = LDA
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = XC_LDA_X
   calc_type%dft_c = XC_LDA_C_VWN
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.0_dp
 case('VWN_RPA')
   calc_type%type          = LDA
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = XC_LDA_X
   calc_type%dft_c = XC_LDA_C_VWN_RPA
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.0_dp
 !
 ! GGA functionals
 case('PBEx')
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = XC_GGA_X_PBE
   calc_type%dft_c = 0
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.0_dp
 case('PBE')
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = XC_GGA_X_PBE
   calc_type%dft_c = XC_GGA_C_PBE
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.0_dp
 case('Bx')
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = XC_GGA_X_B88
   calc_type%dft_c = 0
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.0_dp
 case('BLYP')
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = XC_GGA_X_B88
   calc_type%dft_c = XC_GGA_C_LYP
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.0_dp
 case('PW91')
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = XC_GGA_X_PW91
   calc_type%dft_c = XC_GGA_C_PW91
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.0_dp
 case('PW91x')
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = XC_GGA_X_PW91
   calc_type%dft_c = 0
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.0_dp
 !
 ! Meta-GGA functionals
 case('BJ')
   calc_type%need_dft_xc   = .TRUE.
   calc_type%dft_x = XC_MGGA_X_BJ06
   calc_type%dft_c = 0
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.0_dp
 case('RPP')
   calc_type%need_dft_xc   = .TRUE.
   calc_type%dft_x = XC_MGGA_X_RPP09
   calc_type%dft_c = 0
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.0_dp
 !
 ! Hybrid functionals
 case('B3LYP')
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%need_exchange = .TRUE.  
   calc_type%dft_x = XC_HYB_GGA_XC_B3LYP
   calc_type%dft_c = 0
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.20_dp
 case('PBE0','PBE1PBE')
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%need_exchange = .TRUE.  
   calc_type%dft_x = XC_HYB_GGA_XC_PBEH
   calc_type%dft_c = 0
   if(calc_type%is_gw .OR. calc_type%is_mp2) calc_type%need_final_exchange=.TRUE.
   alpha_hybrid = 0.25_dp
#endif
 case default
   stop'error reading calculation type part 1'
 end select


end subroutine init_calculation_type

!=========================================================================
subroutine output_calculation_type(calc_type)
 implicit none
 type(calculation_type),intent(in)   :: calc_type
!=====

  write(*,*) 'type                        ',calc_type%type
  write(*,*) 'need_exchange               ',calc_type%need_exchange
  write(*,*) 'need_final_exchange         ',calc_type%need_final_exchange
  write(*,*) 'need_dft_xc                 ',calc_type%need_dft_xc
  write(*,*) 'need_rpa                    ',calc_type%need_rpa
  write(*,*) 'is_gw                       ',calc_type%is_gw
  write(*,*) 'is_mp2                      ',calc_type%is_mp2
  write(*,*) 'is_ci                       ',calc_type%is_ci
  write(*,*) 'method                      ',calc_type%method  
  write(*,*) 'dft_x                       ',calc_type%dft_x
  write(*,*) 'dft_c                       ',calc_type%dft_c                       

end subroutine output_calculation_type

end module m_calculation_type
