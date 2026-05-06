!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the matrix to transform CART gaussian to PURE (=spherical) gaussian
!
!=========================================================================
#include "molgw.h"
module m_cart_to_pure
  use m_definitions
  use m_warning

  implicit none

  type transform
    real(dp), allocatable :: matrix(:, :)
  end type

  type reorder
    integer, allocatable :: reindex(:)
  end type

  type(transform), allocatable, protected :: cart_to_pure     (:, :)
  type(transform), allocatable, protected :: cart_to_pure_norm(:, :)

  ! reordering from Gaussian order to molgw/libint/libcint
  type(reorder), allocatable, protected :: g2m_cart(:)
  ! reordering from molgw/libint/libcint order to Gaussian
  type(reorder), allocatable, protected :: m2g_cart(:)
  ! reordering from molgw/libint/libcint order to MOLDEN
  type(reorder), allocatable, protected :: m2molden_cart(:)
  type(reorder), allocatable, protected :: molden2m_cart(:)

  integer, parameter :: CARTG = 1
  integer, parameter :: PUREG = 2


contains


!=========================================================================
function get_gaussian_type_tag(gaussian_type)
  character(len=4), intent(in) :: gaussian_type
  integer                      :: get_gaussian_type_tag
  !=====

  select case(gaussian_type)
  case('CART')
    get_gaussian_type_tag = CARTG
  case('PURE')
    get_gaussian_type_tag = PUREG
  end select

end function get_gaussian_type_tag


!=========================================================================
pure function number_basis_function_am(gaussian_type, am)
  character(len=4), intent(in) :: gaussian_type
  integer, intent(in)          :: am
  integer                      :: number_basis_function_am
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
  integer, intent(in) :: intin
  real(dp)            :: double_factorial
  !=====
  ! just hard-coded for some small integers

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
subroutine setup_cart_to_pure_transforms(pypzpx_order_in)

  logical, intent(in) :: pypzpx_order_in
  !=====
  integer  :: ni, nic
  integer  :: ii, jj, kk
  integer  :: nx, ny, nz
  integer  :: il, im
  integer  :: it, iu, is
  real(dp) :: rtmp
  !=====

  write(stdout, '(/,1x,a,i2)') 'Setting up the cartesian to pure transforms up to l= ', MOLGW_LMAX

  allocate(cart_to_pure     (0:MOLGW_LMAX, 2))
  allocate(cart_to_pure_norm(0:MOLGW_LMAX, 2))

  !
  ! First setup trivial transforms in the case of CARTESIAN gaussians
  !
  do il=0, MOLGW_LMAX
    nic = number_basis_function_am('CART', il)
    allocate(cart_to_pure     (il, CARTG)%matrix(nic, nic))
    allocate(cart_to_pure_norm(il, CARTG)%matrix(nic, nic))
    cart_to_pure(il, CARTG)%matrix(:, :) = 0.0_dp
    do ii=1, nic
      cart_to_pure(il, CARTG)%matrix(ii, ii) = 1.0_dp
    enddo
  enddo


  !
  ! Second setup the complicated transforms in the case of PURE gaussians
  !
  do il=0, MOLGW_LMAX
    nic = number_basis_function_am('CART', il)
    ni  = number_basis_function_am('PURE', il)
    allocate(cart_to_pure(il, PUREG)%matrix(nic, ni))
    allocate(cart_to_pure_norm(il, PUREG)%matrix(nic, ni))
    cart_to_pure_norm(il, PUREG)%matrix(:, :) = 0.0_dp

    kk=0
    do ii=0, il
      nx = il - ii
      do jj=0, ii
        kk = kk + 1
        ny = ii - jj
        nz = jj

        rtmp = 0.0_dp
        do it=0, il/2
          do iu=0, it
            if( 2*it-2*iu == nx .AND. 2*iu == ny .AND. il - 2*it == nz ) then
              rtmp = rtmp + (-1)**it * 0.50_dp**(2*it) * cnk(il-it, it) * cnk(il, it) * cnk(it, iu)
            endif
          enddo
        enddo
        cart_to_pure_norm(il, PUREG)%matrix(kk, il+1) = rtmp / SQRT( double_factorial(2*il-1) )

        do im=1, il
          rtmp = 0.0_dp
          do it=0, (il-im)/2
            do iu=0, it
              do is=0, (im-1)/2
                if( im + 2*it - 2*iu - 2*is -1 == nx &
                 .AND. 2*iu + 2*is + 1 == ny &
                 .AND. il - im - 2*it == nz ) then
                  rtmp = rtmp + (-1)**(it+is) * 0.50_dp**(im+2*it) * SQRT( 2.0_dp * ank(il+im, il) / ank(il, il-im) ) &
                           * cnk(it, iu) * cnk(im, 2*is+1) * cnk(il-it, im+it) * cnk(il, it)
                endif

              enddo
            enddo
          enddo
          cart_to_pure_norm(il, PUREG)%matrix(kk, il+1-im) = rtmp / SQRT( double_factorial(2*il-1) )
        enddo

        do im=1, il
          rtmp = 0.0_dp
          do it=0, (il-im)/2
            do iu=0, it
              do is=0, im/2
                if( im + 2*it - 2*iu - 2*is == nx &
                 .AND. 2*iu + 2*is == ny &
                 .AND. il - im - 2*it == nz ) then
                  rtmp = rtmp + (-1)**(it+is) * 0.50_dp**(im+2*it) * SQRT( 2.0_dp * ank(il+im, il) / ank(il, il-im) ) &
                           * cnk(it, iu) * cnk(im, 2*is) * cnk(il-it, im+it) * cnk(il, it)
                endif

              enddo
            enddo
          enddo
          cart_to_pure_norm(il, PUREG)%matrix(kk, il+1+im) = rtmp / SQRT( double_factorial(2*il-1) )
        enddo


      enddo
    enddo
  enddo

  ! Fix the p-orbital ordering if necessary
  if( .NOT. pypzpx_order_in ) then
    cart_to_pure_norm(1, PUREG)%matrix(1, :) = 0.0_dp
    cart_to_pure_norm(1, PUREG)%matrix(1, 3) = 1.0_dp / SQRT( double_factorial(1) )
    cart_to_pure_norm(1, PUREG)%matrix(2, 1) = 1.0_dp / SQRT( double_factorial(1) )
    cart_to_pure_norm(1, PUREG)%matrix(3, 2) = 1.0_dp / SQRT( double_factorial(1) )
  endif


  !
  ! Introduce the normalization coefficient part that depends on (nx,ny,nz)
  do il=0, MOLGW_LMAX
    kk=0
    do ii=0, il
      nx = il - ii
      do jj=0, ii
        kk = kk + 1
        ny = ii - jj
        nz = jj
        cart_to_pure_norm(il, CARTG)%matrix(kk, :) = cart_to_pure(il, CARTG)%matrix(kk, :) &
                 / SQRT( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) )

        cart_to_pure(il, PUREG)%matrix(kk, :) = cart_to_pure_norm(il, PUREG)%matrix(kk, :) &
                     * SQRT( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) )

      enddo
    enddo
  enddo

  !
  ! Store the transforms for CARTESIAN from GAUSSIAN/MOLDEN to MOLGW/LIBCINT
  !
  allocate(g2m_cart(0:MOLGW_LMAX))
  allocate(m2g_cart(0:MOLGW_LMAX))
  allocate(m2molden_cart(0:MOLGW_LMAX))
  allocate(molden2m_cart(0:MOLGW_LMAX))
  ! s
  allocate(g2m_cart(0)%reindex(1))
  allocate(m2g_cart(0)%reindex(1))
  allocate(m2molden_cart(0)%reindex(1))
  allocate(molden2m_cart(0)%reindex(1))
  g2m_cart(0)%reindex(1) = 1
  m2g_cart(0)%reindex(1) = 1
  m2molden_cart(0)%reindex(:) = m2g_cart(0)%reindex(:)
  molden2m_cart(0)%reindex(:) = g2m_cart(0)%reindex(:)
  ! p
  allocate(g2m_cart(1)%reindex(3))
  allocate(m2g_cart(1)%reindex(3))
  allocate(m2molden_cart(1)%reindex(3))
  allocate(molden2m_cart(1)%reindex(3))
  g2m_cart(1)%reindex(1) = 1
  g2m_cart(1)%reindex(2) = 2
  g2m_cart(1)%reindex(3) = 3
  m2g_cart(1)%reindex(1) = 1
  m2g_cart(1)%reindex(2) = 2
  m2g_cart(1)%reindex(3) = 3
  m2molden_cart(1)%reindex(:) = m2g_cart(1)%reindex(:)
  molden2m_cart(1)%reindex(:) = g2m_cart(1)%reindex(:)
  ! d
  ! gaussian d orbital order is xx, yy, zz, xy, xz, yz
  ! libint   d orbital order is xx, xy, xz, yy, yz, zz
  allocate(g2m_cart(2)%reindex(6))
  allocate(m2g_cart(2)%reindex(6))
  allocate(m2molden_cart(2)%reindex(6))
  allocate(molden2m_cart(2)%reindex(6))
  g2m_cart(2)%reindex(1) = 1
  g2m_cart(2)%reindex(2) = 4
  g2m_cart(2)%reindex(3) = 5
  g2m_cart(2)%reindex(4) = 2
  g2m_cart(2)%reindex(5) = 6
  g2m_cart(2)%reindex(6) = 3
  m2g_cart(2)%reindex(1) = 1
  m2g_cart(2)%reindex(2) = 4
  m2g_cart(2)%reindex(3) = 6
  m2g_cart(2)%reindex(4) = 2
  m2g_cart(2)%reindex(5) = 3
  m2g_cart(2)%reindex(6) = 5
  m2molden_cart(2)%reindex(:) = m2g_cart(2)%reindex(:)
  molden2m_cart(2)%reindex(:) = g2m_cart(2)%reindex(:)

  ! f
  ! gaussian f orbital order is XXX , YYY , ZZZ , XYY , XXY , XXZ , XZZ , YZZ , YYZ , XYZ
  ! libint   f orbital order is xxx , xxy , xxz , xyy , xyz , xzz , yyy , yyz , yzz , zzz
  allocate(g2m_cart(3)%reindex(10))
  allocate(m2g_cart(3)%reindex(10))
  allocate(m2molden_cart(3)%reindex(10))
  allocate(molden2m_cart(3)%reindex(10))
  g2m_cart(3)%reindex( 1) = 1
  g2m_cart(3)%reindex( 2) = 5
  g2m_cart(3)%reindex( 3) = 6
  g2m_cart(3)%reindex( 4) = 4
  g2m_cart(3)%reindex( 5) =10
  g2m_cart(3)%reindex( 6) = 7
  g2m_cart(3)%reindex( 7) = 2
  g2m_cart(3)%reindex( 8) = 9
  g2m_cart(3)%reindex( 9) = 8
  g2m_cart(3)%reindex(10) = 3
  m2g_cart(3)%reindex( 1) = 1
  m2g_cart(3)%reindex( 2) = 7
  m2g_cart(3)%reindex( 3) =10
  m2g_cart(3)%reindex( 4) = 4
  m2g_cart(3)%reindex( 5) = 2
  m2g_cart(3)%reindex( 6) = 3
  m2g_cart(3)%reindex( 7) = 6
  m2g_cart(3)%reindex( 8) = 9
  m2g_cart(3)%reindex( 9) = 8
  m2g_cart(3)%reindex(10) = 5
  m2molden_cart(3)%reindex(:) = m2g_cart(3)%reindex(:)
  molden2m_cart(3)%reindex(:) = g2m_cart(3)%reindex(:)

  ! g and following are like this:
  ! gaussian g orbital order is ZZZZ YZZZ YYZZ YYYZ YYYY XZZZ XYZZ XYYZ XYYY XXZZ XXYZ XXYY XXXZ XXXY XXXX
  ! libint   g orbital order is xxxx xxxy xxxz xxyy xxyz xxzz xyyy xyyz xyzz xzzz yyyy yyyz yyzz yzzz zzzz 
  do il=4, MOLGW_LMAX
    ni = number_basis_function_am('CART', il)
    allocate(g2m_cart(il)%reindex(ni))
    allocate(m2g_cart(il)%reindex(ni))
    allocate(m2molden_cart(il)%reindex(ni))
    do ii=1, ni
      g2m_cart(il)%reindex(ii) = ni + 1 - ii
      m2g_cart(il)%reindex(ii) = ni + 1 - ii
    enddo
  enddo

  molden2m_cart(4)%reindex( 1) = 1
  molden2m_cart(4)%reindex( 2) = 4
  molden2m_cart(4)%reindex( 3) = 5
  molden2m_cart(4)%reindex( 4) =10
  molden2m_cart(4)%reindex( 5) =13
  molden2m_cart(4)%reindex( 6) =11
  molden2m_cart(4)%reindex( 7) = 6
  molden2m_cart(4)%reindex( 8) =14 
  molden2m_cart(4)%reindex( 9) =15
  molden2m_cart(4)%reindex(10) = 8
  molden2m_cart(4)%reindex(11) = 2
  molden2m_cart(4)%reindex(12) = 7
  molden2m_cart(4)%reindex(13) =12
  molden2m_cart(4)%reindex(14) = 9
  molden2m_cart(4)%reindex(15) = 3

  m2molden_cart(4)%reindex( 1)= 1
  m2molden_cart(4)%reindex( 2)= 11
  m2molden_cart(4)%reindex( 3)= 15
  m2molden_cart(4)%reindex( 4)= 2
  m2molden_cart(4)%reindex( 5)= 3
  m2molden_cart(4)%reindex( 6)= 7
  m2molden_cart(4)%reindex( 7)= 12
  m2molden_cart(4)%reindex( 8)= 10
  m2molden_cart(4)%reindex( 9)= 14
  m2molden_cart(4)%reindex(10)= 4
  m2molden_cart(4)%reindex(11)= 6
  m2molden_cart(4)%reindex(12)= 13
  m2molden_cart(4)%reindex(13)= 5
  m2molden_cart(4)%reindex(14)= 8
  m2molden_cart(4)%reindex(15)= 9

  write(stdout, *) 'Transformations set up completed for both CARTESIAN and PURE Gaussians'
  write(stdout, *)

end subroutine setup_cart_to_pure_transforms


!=========================================================================
subroutine destroy_cart_to_pure_transforms()

  !=====
  integer :: il
  !=====

  do il=0, MOLGW_LMAX
    deallocate(cart_to_pure_norm(il, CARTG)%matrix)
    deallocate(cart_to_pure_norm(il, PUREG)%matrix)
    deallocate(cart_to_pure(il, CARTG)%matrix)
    deallocate(cart_to_pure(il, PUREG)%matrix)
    deallocate(g2m_cart(il)%reindex)
    deallocate(m2g_cart(il)%reindex)
    deallocate(m2molden_cart(il)%reindex)
    deallocate(molden2m_cart(il)%reindex)
  enddo
  deallocate(cart_to_pure_norm)
  deallocate(cart_to_pure)
  deallocate(g2m_cart)
  deallocate(m2g_cart)
  deallocate(m2molden_cart)
  deallocate(molden2m_cart)

end subroutine destroy_cart_to_pure_transforms


!=========================================================================
function cnk(n, k)

  integer, intent(in) :: n, k
  real(dp)            :: cnk
  !=====
  integer  :: i
  real(dp) :: num, denom
  !=====

  num   = 1.0_dp
  denom = 1.0_dp
  do i=0, k-1
    num   = num   * REAL(n-i, KIND=dp)
    denom = denom * ( REAL(i, KIND=dp) + 1.0_dp)
  enddo
  cnk = num / denom

end function cnk


!=========================================================================
function ank(n, k)

  integer, intent(in) :: n, k
  real(dp)            :: ank
  !=====
  integer  :: i
  !=====

  ank   = 1.0_dp
  do i=n, k+1, -1
    ank   = ank   * REAL(i, KIND=dp)
  enddo

end function ank


!=========================================================================
end module m_cart_to_pure
