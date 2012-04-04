module m_definitions

 integer,parameter  :: dp=kind(0.d0)
 integer,parameter  :: sp=kind(0.)
 integer,parameter  :: prec_eri=dp

 !
 ! Physical constants
 real(dp),parameter :: Ha_eV =27.2113961_dp
 real(dp),parameter :: bohr_A=0.529177249_dp
 real(dp),parameter :: pi    =3.14159265358979323_dp
 real(dp),parameter :: pi2   =pi**2
 complex(dp)        :: im=(0.0_dp,1.0_dp)

 real(dp),parameter :: completely_empty=1.0d-5

 integer,parameter  :: NOMEGA=12 !TODO remove this parameter

contains
!==========================================
 function gamma_function(rin)
 implicit none
 real(dp),intent(in) :: rin
 real(dp) :: gamma_function
!=====
 integer :: nlocal
!=====

 !
 ! just hard coded for some small half-integers
 if( ABS( rin - NINT(rin) ) - 0.5  > 1.d-6 ) stop'GAMMA FUNCTION NOT CODED'

 nlocal = FLOOR(rin-0.499999)

 select case(nlocal)
 case(0)
  gamma_function = 1.0_dp
 case(1)
  gamma_function = 0.5_dp
 case(2)
  gamma_function = 0.75_dp
 case(3)
  gamma_function = 15.0_dp / 8.0_dp
 case(4)
  gamma_function = 105.0_dp / 16.0_dp
 case(5)
  gamma_function = 945.0_dp / 32.0_dp
 case(6)
  gamma_function = 10395.0_dp / 64.0_dp
 case default
  gamma_function = double_factorial(nlocal+1)/ 2**( 0.5_dp * (nlocal+1)/2)
 end select

 gamma_function = gamma_function * SQRT(pi)

 end function gamma_function

!==========================================
 function double_factorial(intin)
 implicit none
 integer,intent(in) :: intin
! integer :: double_factorial
 real(dp) :: double_factorial
!=====
 ! just hard coded for some small integers

 select case (intin)
 case(-1) 
   double_factorial = 1
 case( 0) 
   double_factorial = 1
 case( 1) 
   double_factorial = 1
 case( 2) 
   double_factorial = 2
 case( 3) 
   double_factorial = 3
 case( 4) 
   double_factorial = 8
 case( 5) 
   double_factorial = 15
 case( 6) 
   double_factorial = 48
 case( 7) 
   double_factorial = 105
 case( 8) 
   double_factorial = 384
 case( 9) 
   double_factorial = 945
 case(10) 
   double_factorial = 3840
 case(11) 
   double_factorial = 10395
 case(12) 
   double_factorial = 46080
 case(13) 
   double_factorial = 135135
 case(14) 
   double_factorial = 645120
 case(15)
   double_factorial = 2027025
 case(16) 
   double_factorial = 10321920
 case(17) 
   double_factorial = 34459425
 case(18) 
   double_factorial = 185794560
 case(19) 
   double_factorial = 654729075
 case(20) 
   double_factorial = 3715891200
 case(21) 
   double_factorial = 13749310575
 case(22)
   double_factorial = 81749606400
 case(23) 
   double_factorial = 316234143225
 case(25) 
   double_factorial = 7905853580625
 case(27) 
   double_factorial = 213458046676875
 case(29) 
   double_factorial = 6190283353629375
 case(31) 
   double_factorial = 191898783962510625
 case default
   write(*,*) 'integer =',intin
   !stop'double factorial not coded for this integer value'
   write(*,*) 'double factorial not coded for this integer value'
   double_factorial = 1
 end select

 end function double_factorial

end module m_definitions

