!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the matrix to transform CART gaussian to PURE (=spherical) gaussian
!
!=========================================================================
module m_cart_to_pure
 use m_definitions
 use m_warning

 integer,parameter              :: LMAX_TRANSFORM      = 7
 integer,parameter              :: LMAX_TRANSFORM_PURE = 5

 type transform
   real(dp),allocatable         :: matrix(:,:)
 end type

 type(transform)                :: cart_to_pure     (0:LMAX_TRANSFORM,2)
 type(transform)                :: cart_to_pure_norm(0:LMAX_TRANSFORM,2)

 integer,parameter              :: CARTG=1
 integer,parameter              :: PUREG=2


contains


!=========================================================================
function get_gaussian_type_tag(gaussian_type)
 implicit none
 character(len=4),intent(in) :: gaussian_type
 integer                     :: get_gaussian_type_tag
!=====

 select case(gaussian_type)
 case('CART')
   get_gaussian_type_tag = CARTG
 case('PURE')
   get_gaussian_type_tag = PUREG
 end select
 
end function get_gaussian_type_tag


!=========================================================================
function number_basis_function_am(gaussian_type,am)
 implicit none
 character(len=4),intent(in) :: gaussian_type
 integer,intent(in)          :: am
 integer                     :: number_basis_function_am
!=====

 ! Above LMAX_TRANSFORM_PURE only Cartesian gaussian are implemented
 if( am <= LMAX_TRANSFORM_PURE ) then
   select case(gaussian_type)
   case('CART')
     select case(am)
     case(0)
       number_basis_function_am = 1
     case(1)
       number_basis_function_am = 3
     case(2)
       number_basis_function_am = 6
     case(3)
       number_basis_function_am = 10
     case(4)
       number_basis_function_am = 15
     case(5)
       number_basis_function_am = 21
     case(6)
       number_basis_function_am = 28
     case(7)
       number_basis_function_am = 36
     case default
       write(stdout,*) 'am=',am
       call die('number_basis_function_am: not implemented')
     end select
   case('PURE')
     number_basis_function_am = 2 * am + 1
   end select
 else
   select case(am)
   case(0)
     number_basis_function_am = 1
   case(1)
     number_basis_function_am = 3
   case(2)
     number_basis_function_am = 6
   case(3)
     number_basis_function_am = 10
   case(4)
     number_basis_function_am = 15
   case(5)
     number_basis_function_am = 21
   case(6)
     number_basis_function_am = 28
   case(7)
     number_basis_function_am = 36
   case default
     write(stdout,*) 'am=',am
     call die('number_basis_function_am: not implemented')
   end select
 endif

end function number_basis_function_am


!==========================================
function double_factorial(intin)
 implicit none
 integer,intent(in) :: intin
 real(dp) :: double_factorial
!=====
 ! just hard coded for some small integers

 select case (intin)
 case(-1) 
   double_factorial = 1.0_dp
 case( 0) 
   double_factorial = 1.0_dp
 case( 1) 
   double_factorial = 1.0_dp
 case( 2) 
   double_factorial = 2.0_dp
 case( 3) 
   double_factorial = 3.0_dp
 case( 4) 
   double_factorial = 8.0_dp
 case( 5) 
   double_factorial = 15.0_dp
 case( 6) 
   double_factorial = 48.0_dp
 case( 7) 
   double_factorial = 105.0_dp
 case( 8) 
   double_factorial = 384.0_dp
 case( 9) 
   double_factorial = 945.0_dp
 case(10) 
   double_factorial = 3840.0_dp
 case(11) 
   double_factorial = 10395.0_dp
 case(12) 
   double_factorial = 46080.0_dp
 case(13) 
   double_factorial = 135135.0_dp
 case(14) 
   double_factorial = 645120.0_dp
 case(15)
   double_factorial = 2027025.0_dp
 case(16) 
   double_factorial = 10321920.0_dp
 case(17) 
   double_factorial = 34459425.0_dp
 case(18) 
   double_factorial = 185794560.0_dp
 case(19) 
   double_factorial = 654729075.0_dp
 case(20) 
   double_factorial = 3715891200.0_dp
 case(21) 
   double_factorial = 13749310575.0_dp
 case(22)
   double_factorial = 81749606400.0_dp
 case(23) 
   double_factorial = 316234143225.0_dp
 case(25) 
   double_factorial = 7905853580625.0_dp
 case(27) 
   double_factorial = 213458046676875.0_dp
 case(29) 
   double_factorial = 6190283353629375.0_dp
 case(31) 
   double_factorial = 191898783962510625.0_dp
 case default
   write(stdout,*) 'integer =',intin
   write(stdout,*) 'double factorial not coded for this integer value'
   double_factorial = 1
 end select

end function double_factorial


!=========================================================================
subroutine setup_cart_to_pure_transforms()
 implicit none

!=====
 integer  :: il,ni,ii,jj,kk
 integer  :: nx,ny,nz
!=====


 write(stdout,'(/,1x,a)') 'Setting up the cartesian to pure transforms'


 !
 ! First setup trivial transforms in the case of CARTESIAN gaussians
 !
 do il=0,LMAX_TRANSFORM
   ni = number_basis_function_am('CART',il)
   allocate(cart_to_pure     (il,CARTG)%matrix(ni,ni))
   allocate(cart_to_pure_norm(il,CARTG)%matrix(ni,ni))
   cart_to_pure(il,CARTG)%matrix(:,:) = 0.0_dp
   do ii=1,ni
     cart_to_pure(il,CARTG)%matrix(ii,ii) = 1.0_dp
   enddo
 enddo


 !
 ! Second setup the complicated transforms in the case of PURE gaussians
 !
 ! Formula were read from Ref. 
 ! H.B. Schlegel and M.J. Frisch, INTERNATIONAL JOURNAL OF QUANTUM CHEMISTRY  54, 83-87 (1995).
 ! Formulas available up to l = 5, then use cartesian gaussian anyway

 !
 ! Transform for momentum S
 allocate(cart_to_pure     (0,PUREG)%matrix(1,1))
 allocate(cart_to_pure_norm(0,PUREG)%matrix(1,1))
 cart_to_pure(0,PUREG)%matrix(1,1) = 1.0_dp

 !
 ! Transform for momentum P
 allocate(cart_to_pure     (1,PUREG)%matrix(3,3))
 allocate(cart_to_pure_norm(1,PUREG)%matrix(3,3))
 cart_to_pure(1,PUREG)%matrix(:,:) = 0.0_dp
 cart_to_pure(1,PUREG)%matrix(3,2) = 1.0_dp
 cart_to_pure(1,PUREG)%matrix(1,3) = 1.0_dp
 cart_to_pure(1,PUREG)%matrix(2,1) = 1.0_dp

 !
 ! Transform for momentum D
 allocate(cart_to_pure     (2,PUREG)%matrix(6,5))
 allocate(cart_to_pure_norm(2,PUREG)%matrix(6,5))
 cart_to_pure(2,PUREG)%matrix(:,:) =  0.0_dp
 cart_to_pure(2,PUREG)%matrix(2,1) =  1.0_dp
 cart_to_pure(2,PUREG)%matrix(5,2) =  1.0_dp
 cart_to_pure(2,PUREG)%matrix(6,3) =  1.0_dp
 cart_to_pure(2,PUREG)%matrix(1,3) = -0.5_dp
 cart_to_pure(2,PUREG)%matrix(4,3) = -0.5_dp
 cart_to_pure(2,PUREG)%matrix(3,4) =  1.0_dp
 cart_to_pure(2,PUREG)%matrix(1,5) =  SQRT(3.0_dp/4.0_dp)
 cart_to_pure(2,PUREG)%matrix(4,5) = -SQRT(3.0_dp/4.0_dp)

 !
 ! Transform for momentum F
 allocate(cart_to_pure     (3,PUREG)%matrix(10,7))
 allocate(cart_to_pure_norm(3,PUREG)%matrix(10,7))
 cart_to_pure(3,PUREG)%matrix( :,:) =  0.0_dp
 cart_to_pure(3,PUREG)%matrix(10,4) =  1.0_dp
 cart_to_pure(3,PUREG)%matrix( 3,4) = -3.0_dp /(2.0_dp *SQRT(5.0_dp))
 cart_to_pure(3,PUREG)%matrix( 8,4) = -3.0_dp /(2.0_dp *SQRT(5.0_dp))
 cart_to_pure(3,PUREG)%matrix( 6,5) =  SQRT(6.0_dp/5.0_dp)
 cart_to_pure(3,PUREG)%matrix( 1,5) = -SQRT(6.0_dp)/4.0_dp
 cart_to_pure(3,PUREG)%matrix( 4,5) = -SQRT(6.0_dp/5.0_dp)/4.0_dp
 cart_to_pure(3,PUREG)%matrix( 9,3) =  SQRT(6.0_dp/5.0_dp)
 cart_to_pure(3,PUREG)%matrix( 7,3) = -SQRT(6.0_dp)/4.0_dp
 cart_to_pure(3,PUREG)%matrix( 2,3) = -SQRT(6.0_dp/5.0_dp)/4.0_dp
 cart_to_pure(3,PUREG)%matrix( 3,6) =  SQRT(3.0_dp/4.0_dp)
 cart_to_pure(3,PUREG)%matrix( 8,6) = -SQRT(3.0_dp/4.0_dp)
 cart_to_pure(3,PUREG)%matrix( 5,2) = -1.0_dp
 cart_to_pure(3,PUREG)%matrix( 1,7) =  SQRT(10.0_dp)/4.0_dp
 cart_to_pure(3,PUREG)%matrix( 4,7) = -SQRT(2.0_dp)*3.0_dp/4.0_dp
 cart_to_pure(3,PUREG)%matrix( 7,1) = -SQRT(10.0_dp)/4.0_dp
 cart_to_pure(3,PUREG)%matrix( 2,1) =  SQRT(2.0_dp)*3.0_dp/4.0_dp

 !
 ! Transform for momentum G
 allocate(cart_to_pure     (4,PUREG)%matrix(15,9))
 allocate(cart_to_pure_norm(4,PUREG)%matrix(15,9))
 cart_to_pure(4,PUREG)%matrix( :,:) =  0.0_dp
 cart_to_pure(4,PUREG)%matrix(15,5) =  1.0_dp
 cart_to_pure(4,PUREG)%matrix( 1,5) =  3.0/8.0
 cart_to_pure(4,PUREG)%matrix(11,5) =  3.0/8.0
 cart_to_pure(4,PUREG)%matrix( 6,5) = -3.0*SQRT(3.0/35.0)
 cart_to_pure(4,PUREG)%matrix(13,5) = -3.0*SQRT(3.0/35.0)
 cart_to_pure(4,PUREG)%matrix( 4,5) =  3.0*SQRT(3.0/35.0)/4.0
 cart_to_pure(4,PUREG)%matrix(10,6) =  SQRT(10.0/7.0)
 cart_to_pure(4,PUREG)%matrix( 3,6) = -0.75*SQRT(10.0/7.0)
 cart_to_pure(4,PUREG)%matrix( 8,6) = -0.75*SQRT( 2.0/7.0)
 cart_to_pure(4,PUREG)%matrix(14,4) =  SQRT(10.0/7.0)
 cart_to_pure(4,PUREG)%matrix(12,4) = -0.75*SQRT(10.0/7.0)
 cart_to_pure(4,PUREG)%matrix( 5,4) = -0.75*SQRT( 2.0/7.0)
 cart_to_pure(4,PUREG)%matrix( 6,7) =  1.5*SQRT( 3.0/7.0)
 cart_to_pure(4,PUREG)%matrix(13,7) = -1.5*SQRT( 3.0/7.0)
 cart_to_pure(4,PUREG)%matrix( 1,7) = -0.25*SQRT(5.0)
 cart_to_pure(4,PUREG)%matrix(11,7) =  0.25*SQRT(5.0)
 cart_to_pure(4,PUREG)%matrix( 9,3) =  3.0/SQRT(7.0)
 cart_to_pure(4,PUREG)%matrix( 2,3) = -0.5*SQRT(5.0/7.0)
 cart_to_pure(4,PUREG)%matrix( 7,3) = -0.5*SQRT(5.0/7.0)
 cart_to_pure(4,PUREG)%matrix( 3,8) =  0.25*SQRT(10.0)
 cart_to_pure(4,PUREG)%matrix( 8,8) = -0.75*SQRT(2.0)
 cart_to_pure(4,PUREG)%matrix(12,2) = -0.25*SQRT(10.0)
 cart_to_pure(4,PUREG)%matrix( 5,2) =  0.75*SQRT(2.0)
 cart_to_pure(4,PUREG)%matrix( 1,9) =  0.125*SQRT(35.0)
 cart_to_pure(4,PUREG)%matrix(11,9) =  0.125*SQRT(35.0)
 cart_to_pure(4,PUREG)%matrix( 4,9) = -0.75*SQRT(3.0)
 cart_to_pure(4,PUREG)%matrix( 2,1) =  SQRT(5.0/4.0)
 cart_to_pure(4,PUREG)%matrix( 7,1) = -SQRT(5.0/4.0)

 !
 ! Transform for momentum H
 allocate(cart_to_pure     (5,PUREG)%matrix(21,11))
 allocate(cart_to_pure_norm(5,PUREG)%matrix(21,11))
 cart_to_pure(5,PUREG)%matrix( :, :) =  0.0_dp
 cart_to_pure(5,PUREG)%matrix(21, 6) =  1.0_dp
 cart_to_pure(5,PUREG)%matrix(10, 6) = -5.0*SQRT(2.0/21.0)
 cart_to_pure(5,PUREG)%matrix(19, 6) = -5.0*SQRT(2.0/21.0)
 cart_to_pure(5,PUREG)%matrix( 3, 6) =  5.0/8.0*SQRT(2.0)
 cart_to_pure(5,PUREG)%matrix(17, 6) = -5.0/8.0*SQRT(2.0)
 cart_to_pure(5,PUREG)%matrix( 8, 6) =  0.25*SQRT(30.0/7.0)
 cart_to_pure(5,PUREG)%matrix(15, 7) =  SQRT(5.0/3.0)
 cart_to_pure(5,PUREG)%matrix( 6, 7) = -1.5*SQRT(5.0/7.0)
 cart_to_pure(5,PUREG)%matrix(13, 7) = -1.5/SQRT(7.0)
 cart_to_pure(5,PUREG)%matrix( 1, 7) = -0.125*SQRT(15.0)
 cart_to_pure(5,PUREG)%matrix(11, 7) =  0.125*SQRT(5.0/3.0)
 cart_to_pure(5,PUREG)%matrix( 4, 7) =  0.25*SQRT(5.0/7.0)
 cart_to_pure(5,PUREG)%matrix(20, 5) =  SQRT(5.0/3.0)
 cart_to_pure(5,PUREG)%matrix(18, 5) = -1.5*SQRT(5.0/7.0)
 cart_to_pure(5,PUREG)%matrix( 9, 5) = -1.5/SQRT(7.0)
 cart_to_pure(5,PUREG)%matrix(16, 5) = -0.125*SQRT(15.0)
 cart_to_pure(5,PUREG)%matrix( 2, 5) =  0.125*SQRT(5.0/3.0)
 cart_to_pure(5,PUREG)%matrix( 7, 5) =  0.25*SQRT(5.0/7.0)
 cart_to_pure(5,PUREG)%matrix(10, 8) =  SQRT(5.0/4.0)
 cart_to_pure(5,PUREG)%matrix(19, 8) = -SQRT(5.0/4.0)
 cart_to_pure(5,PUREG)%matrix( 3, 8) = -0.25*SQRT(35.0/3.0)
 cart_to_pure(5,PUREG)%matrix(17, 8) =  0.25*SQRT(35.0/3.0)
 cart_to_pure(5,PUREG)%matrix(14, 4) =  SQRT(5.0/3.0)
 cart_to_pure(5,PUREG)%matrix( 5, 4) =  -0.5*SQRT(5.0/3.0)
 cart_to_pure(5,PUREG)%matrix(12, 4) =  -0.5*SQRT(5.0/3.0)
 cart_to_pure(5,PUREG)%matrix( 6, 9) =  0.5*SQRT(10.0/3.0)
 cart_to_pure(5,PUREG)%matrix(13, 9) = -0.5*SQRT(6.0)
 cart_to_pure(5,PUREG)%matrix( 1, 9) = -SQRT(70.0)/16.0
 cart_to_pure(5,PUREG)%matrix(11, 9) =  SQRT(70.0)/16.0
 cart_to_pure(5,PUREG)%matrix( 4, 9) =  0.125*SQRT(10.0/3.0)
 cart_to_pure(5,PUREG)%matrix(18, 3) = -0.5*SQRT(10.0/3.0)
 cart_to_pure(5,PUREG)%matrix( 9, 3) =  0.5*SQRT(6.0)
 cart_to_pure(5,PUREG)%matrix(16, 3) =  SQRT(70.0)/16.0
 cart_to_pure(5,PUREG)%matrix( 2, 3) = -SQRT(70.0)/16.0
 cart_to_pure(5,PUREG)%matrix( 7, 3) = -0.125*SQRT(10.0/3.0)
 cart_to_pure(5,PUREG)%matrix( 3,10) =  0.125*SQRT(35.0)
 cart_to_pure(5,PUREG)%matrix(17,10) =  0.125*SQRT(35.0)
 cart_to_pure(5,PUREG)%matrix( 8,10) = -0.75*SQRT(3.0)
 cart_to_pure(5,PUREG)%matrix( 5, 2) =  0.5*SQRT(5.0)
 cart_to_pure(5,PUREG)%matrix(12, 2) = -0.5*SQRT(5.0)
 cart_to_pure(5,PUREG)%matrix( 1,11) =  3.0*SQRT(14.0)/16.0
 cart_to_pure(5,PUREG)%matrix(11,11) = -5.0*SQRT(14.0)/16.0
 cart_to_pure(5,PUREG)%matrix( 4,11) = -5.0*SQRT(6.0)/8.0
 cart_to_pure(5,PUREG)%matrix(16, 1) =  3.0*SQRT(14.0)/16.0
 cart_to_pure(5,PUREG)%matrix( 2, 1) = -5.0*SQRT(14.0)/16.0
 cart_to_pure(5,PUREG)%matrix( 7, 1) =  5.0*SQRT(6.0)/8.0


 !
 ! Complement with diagonal if necessary (formulas not available in the paper)
 do il=LMAX_TRANSFORM_PURE+1,LMAX_TRANSFORM
   ni = number_basis_function_am('CART',il)
   allocate(cart_to_pure     (il,PUREG)%matrix(ni,ni))
   allocate(cart_to_pure_norm(il,PUREG)%matrix(ni,ni))
   cart_to_pure(il,PUREG)%matrix(:,:) = cart_to_pure(il,CARTG)%matrix(:,:)
 enddo



 !
 ! Introduce the normalization coefficient part that depends on (nx,ny,nz)
 do il=0,LMAX_TRANSFORM
   kk=0
   do ii=0,il
     nx = il - ii
     do jj=0,ii
       kk = kk + 1
       ny = ii - jj
       nz = jj
       cart_to_pure_norm(il,CARTG)%matrix(kk,:) = cart_to_pure(il,CARTG)%matrix(kk,:) &
                 / SQRT( REAL( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) , dp ) )

       cart_to_pure_norm(il,PUREG)%matrix(kk,:) = cart_to_pure(il,PUREG)%matrix(kk,:) &
                 / SQRT( REAL( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) , dp ) )
         
     enddo
   enddo
 enddo


 write(stdout,*) 'Transformations set up completed for both CARTESIAN and PURE Gaussians'
 write(stdout,*) 

end subroutine setup_cart_to_pure_transforms


!=========================================================================
end module m_cart_to_pure
