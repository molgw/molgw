module m_calculation_type
 use m_warning

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
 integer,parameter :: GWmu=11
 integer,parameter :: LDA0=12
! Method definitions
 integer,parameter :: perturbative=101
 integer,parameter :: QS          =102

 type calculation_type
   integer :: type
   logical :: need_exchange
   logical :: need_final_exchange
   logical :: need_dft_xc
   logical :: need_rpa
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
 case('LDA')
   calc_type%type          = LDA
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = 1
!   calc_type%dft_c = 7  ! VWN5
   calc_type%dft_c = 12 ! PW
   if(calc_type%is_gw) calc_type%need_final_exchange=.TRUE.
 case('PBEx')
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = 101
   calc_type%dft_c = 0
   if(calc_type%is_gw) calc_type%need_final_exchange=.TRUE.
 case('PBE')
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = 101
   calc_type%dft_c = 130
   if(calc_type%is_gw) calc_type%need_final_exchange=.TRUE.
 case('PW91')
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = 109
   calc_type%dft_c = 134
   if(calc_type%is_gw) calc_type%need_final_exchange=.TRUE.
 case('GW')
   calc_type%need_exchange = .TRUE.  
   calc_type%is_gw         = .TRUE.
   calc_type%method        = QS
 case('PW0') ! LDA-based Hybrid functional
   calc_type%need_exchange = .TRUE.  
   calc_type%need_dft_xc   = .TRUE.  
   calc_type%dft_x = 1
   calc_type%dft_c = 12 ! PW
 case default
   stop'error reading calculation type part 1'
 end select


end subroutine init_calculation_type

end module m_calculation_type
