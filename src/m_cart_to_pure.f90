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

 type transform
   real(dp),allocatable    :: matrix(:,:)
 end type

 type(transform),allocatable,protected :: cart_to_pure     (:,:)
 type(transform),allocatable,protected :: cart_to_pure_norm(:,:)

 integer,parameter         :: CARTG=1
 integer,parameter         :: PUREG=2


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
pure function number_basis_function_am(gaussian_type,am)
 implicit none
 character(len=4),intent(in) :: gaussian_type
 integer,intent(in)          :: am
 integer                     :: number_basis_function_am
!=====

 select case(gaussian_type)
 case('CART')
   number_basis_function_am = ( ( am + 1 ) * ( am + 2 ) ) / 2
 case('PURE')
   number_basis_function_am = 2 * am + 1
 end select

end function number_basis_function_am


!==========================================
pure function double_factorial(intin)
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
 case(24)
   double_factorial = 1961990553600.0_dp
 case(25)
   double_factorial = 7905853580625.0_dp
 case(26)
   double_factorial = 51011754393600.0_dp
 case(27)
   double_factorial = 213458046676875.0_dp
 case(28)
   double_factorial = 1428329123020800.0_dp
 case(29)
   double_factorial = 6190283353629375.0_dp
 case(31)
   double_factorial = 191898783962510625.0_dp
 end select

end function double_factorial


!=========================================================================
subroutine setup_cart_to_pure_transforms()
 implicit none

!=====
 integer  :: ni,nic
 integer  :: ii,jj,kk
 integer  :: nx,ny,nz
 integer  :: il,im
 integer  :: it,iu,is
 real(dp) :: rtmp
!=====

 write(stdout,'(/,1x,a,i2)') 'Setting up the cartesian to pure transforms up to l= ',MOLGW_LMAX

 allocate(cart_to_pure     (0:MOLGW_LMAX,2))
 allocate(cart_to_pure_norm(0:MOLGW_LMAX,2))

 !
 ! First setup trivial transforms in the case of CARTESIAN gaussians
 !
 do il=0,MOLGW_LMAX
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
 do il=0,MOLGW_LMAX
   nic = number_basis_function_am('CART',il)
   ni  = number_basis_function_am('PURE',il)
   allocate(cart_to_pure(il,PUREG)%matrix(nic,ni))
   allocate(cart_to_pure_norm(il,PUREG)%matrix(nic,ni))
   cart_to_pure_norm(il,PUREG)%matrix(:,:) = 0.0_dp

   kk=0
   do ii=0,il
     nx = il - ii
     do jj=0,ii
       kk = kk + 1
       ny = ii - jj
       nz = jj

       rtmp = 0.0_dp
       do it=0,il/2
         do iu=0,it
           if( 2*it-2*iu == nx .AND. 2*iu == ny .AND. il - 2*it == nz ) then
             rtmp = rtmp + (-1)**it * 0.50_dp**(2*it) * cnk(il-it,it) * cnk(il,it) * cnk(it,iu)
           endif
         enddo
       enddo
       cart_to_pure_norm(il,PUREG)%matrix(kk,il+1) = rtmp / SQRT( double_factorial(2*il-1) )

       do im=1,il
         rtmp = 0.0_dp
         do it=0,(il-im)/2
           do iu=0,it
             do is=0,(im-1)/2
               if( im + 2*it - 2*iu - 2*is -1 == nx &
                 .AND. 2*iu + 2*is + 1 == ny &
                 .AND. il - im - 2*it == nz ) then
                 rtmp = rtmp + (-1)**(it+is) * 0.50_dp**(im+2*it) * SQRT( 2.0_dp * ank(il+im,il) / ank(il,il-im) ) &
                           * cnk(it,iu) * cnk(im,2*is+1) * cnk(il-it,im+it) * cnk(il,it)
               endif

             enddo
           enddo
         enddo
         cart_to_pure_norm(il,PUREG)%matrix(kk,il+1-im) = rtmp / SQRT( double_factorial(2*il-1) )
       enddo

       do im=1,il
         rtmp = 0.0_dp
         do it=0,(il-im)/2
           do iu=0,it
             do is=0,im/2
               if( im + 2*it - 2*iu - 2*is == nx &
                 .AND. 2*iu + 2*is == ny &
                 .AND. il - im - 2*it == nz ) then
                 rtmp = rtmp + (-1)**(it+is) * 0.50_dp**(im+2*it) * SQRT( 2.0_dp * ank(il+im,il) / ank(il,il-im) ) &
                           * cnk(it,iu) * cnk(im,2*is) * cnk(il-it,im+it) * cnk(il,it)
               endif

             enddo
           enddo
         enddo
         cart_to_pure_norm(il,PUREG)%matrix(kk,il+1+im) = rtmp / SQRT( double_factorial(2*il-1) )
       enddo


     enddo
   enddo
 enddo


 !
 ! Introduce the normalization coefficient part that depends on (nx,ny,nz)
 do il=0,MOLGW_LMAX
   kk=0
   do ii=0,il
     nx = il - ii
     do jj=0,ii
       kk = kk + 1
       ny = ii - jj
       nz = jj
       cart_to_pure_norm(il,CARTG)%matrix(kk,:) = cart_to_pure(il,CARTG)%matrix(kk,:) &
                 / SQRT( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) )

       cart_to_pure(il,PUREG)%matrix(kk,:) = cart_to_pure_norm(il,PUREG)%matrix(kk,:) &
                     * SQRT( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) )

     enddo
   enddo
 enddo


 write(stdout,*) 'Transformations set up completed for both CARTESIAN and PURE Gaussians'
 write(stdout,*)

end subroutine setup_cart_to_pure_transforms


!=========================================================================
subroutine destroy_cart_to_pure_transforms()
 implicit none

!=====
 integer :: il
!=====

 do il=0,MOLGW_LMAX
   deallocate(cart_to_pure_norm(il,CARTG)%matrix)
   deallocate(cart_to_pure_norm(il,PUREG)%matrix)
   deallocate(cart_to_pure(il,CARTG)%matrix)
   deallocate(cart_to_pure(il,PUREG)%matrix)
 enddo
 deallocate(cart_to_pure_norm)
 deallocate(cart_to_pure)

end subroutine destroy_cart_to_pure_transforms


!=========================================================================
function cnk(n,k)
 implicit none

 integer,intent(in) :: n,k
 real(dp)           :: cnk
!=====
 integer  :: i
 real(dp) :: num,denom
!=====

 num   = 1.0_dp
 denom = 1.0_dp
 do i=0,k-1
   num   = num   * REAL(n-i,dp)
   denom = denom * ( REAL(i,dp) + 1.0_dp)
 enddo
 cnk = num / denom

end function cnk


!=========================================================================
function ank(n,k)
 implicit none

 integer,intent(in) :: n,k
 real(dp)           :: ank
!=====
 integer  :: i
!=====

 ank   = 1.0_dp
 do i=n,k+1,-1
   ank   = ank   * REAL(i,dp)
 enddo

end function ank


!=========================================================================
end module m_cart_to_pure
