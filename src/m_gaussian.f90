!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to perform operations on an individual Cartesian Gaussian functions
!
!=========================================================================
module m_gaussian
  use m_definitions
  use m_mpi
  use m_cart_to_pure
  use m_string_tools

  ! type containing all the information for one unnormalized cartesian gaussian
  ! (x-x0)**nx * (y-y0)**ny * (z-z0)**nz * exp( - alpha * ( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 ) )
  type gaussian
    integer          :: am
    character(len=1) :: amc
    integer          :: nx,ny,nz
    real(dp)         :: alpha
    real(dp)         :: x0(3)                ! center of the gaussian
    real(dp)         :: norm_factor          ! normalization factor for the gaussian squared
    real(dp)         :: common_norm_factor   ! normalization factor for the gaussian squared
                                             ! without the nx,ny,nz dependence
  end type gaussian

  interface
    subroutine boys_function_c(fnt,n,t) BIND(C,NAME='boys_function_c')
      import :: C_INT,C_DOUBLE
      integer(C_INT),value       :: n
      real(C_DOUBLE),value       :: t
      real(C_DOUBLE),intent(out) :: fnt(0:n)
    end subroutine boys_function_c
  end interface

contains


!=========================================================================
subroutine init_gaussian_cart(nx,ny,nz,alpha,x0,ga)
  implicit none

  integer,intent(in) :: nx,ny,nz
  real(dp),intent(in) :: alpha,x0(3)
  type(gaussian),intent(out) :: ga
  !=====
  !=====

  ga%nx = nx
  ga%ny = ny
  ga%nz = nz
  ga%am = nx + ny + nz
  ga%amc = orbital_momentum_name(ga%am)

  ga%alpha = alpha

  ga%common_norm_factor = ( 2.0_dp / pi )**0.75_dp &
                  * 2.0_dp**ga%am * ga%alpha**( 0.25_dp * ( 2.0_dp*ga%am + 3.0_dp ) )

  ga%norm_factor = ga%common_norm_factor &
                    / SQRT( REAL( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) , dp ) )

  ga%x0(:) = x0(:)

end subroutine init_gaussian_cart


!=========================================================================
function compare_gaussian(g1,g2) result(same_gaussian)
  implicit none

  logical                   :: same_gaussian
  type(gaussian),intent(in) :: g1,g2
  !=====
  !=====

  same_gaussian = .TRUE.

  if( g1%am /= g2%am                          ) same_gaussian = .FALSE.
  if( g1%nx /= g2%nx                          ) same_gaussian = .FALSE.
  if( g1%ny /= g2%ny                          ) same_gaussian = .FALSE.
  if( g1%nz /= g2%nz                          ) same_gaussian = .FALSE.
  if( ABS(g1%alpha - g2%alpha) > 1.0e-5       ) same_gaussian = .FALSE.
  if( ANY(ABS(g1%x0(:) - g2%x0(:)) > 1.0e-5 ) ) same_gaussian = .FALSE.

end function compare_gaussian


!=========================================================================
function eval_gaussian(ga,x)
  implicit none

  type(gaussian),intent(in) :: ga
  real(dp),intent(in) :: x(3)
  real(dp) :: eval_gaussian
  !=====
  real(dp) :: dx(3),dx2
  !=====

  dx(:) = x(:) - ga%x0(:)
  dx2 = SUM( dx(:)**2 )

  eval_gaussian =  dx(1)**ga%nx &
                 * dx(2)**ga%ny &
                 * dx(3)**ga%nz &
                 * EXP( - ga%alpha * dx2 )

  !
  ! normalize the cartesian gaussian
  eval_gaussian = eval_gaussian * ga%norm_factor

end function eval_gaussian


!=========================================================================
function eval_gaussian_fourier(ga,k) RESULT(gaussian_fourier)
  implicit none
  type(gaussian),intent(in) :: ga
  real(dp),intent(in) :: k(3)
  complex(dp) :: gaussian_fourier
  !=====
  real(dp) :: k2,kx2,ky2,kz2
  complex(dp) :: gx,gy,gz
  integer :: il
  !=====

  ! Fourier transform phi(k) = 1/(2pi)**3 int dr e^{-ikr} phi(r)

  kx2 = k(1)**2
  ky2 = k(2)**2
  kz2 = k(3)**2
  k2 = SUM(k(:)**2)

  ! X direction
  select case(ga%nx)
  case(0)
    gx = 1.0_dp
  case(2)
    gx = -kx2 + 2.0_dp * ga%alpha
  case(4)
    gx = kx2**2 - 12.0_dp * ga%alpha * kx2 + 12.0_dp * ga%alpha**2
  case(6)
    gx = -kx2**3 + 30.0_dp * ga%alpha * kx2**2 - 180.0_dp * ga%alpha**2 * kx2 &
                         + 120.0_dp * ga%alpha**3
  case(8)
    gx = kx2**4 - 56.0_dp * ga%alpha * kx2**3 + 840.0_dp * ga%alpha**2 * kx2**2 &
                                - 3360.0_dp * ga%alpha**3 * kx2 + 1680.0_dp * ga%alpha**4
  case(1)
    gx = -im * k(1)
  case(3)
    gx =  im * k(1)**3 -  6.0_dp * ga%alpha    * im * k(1)
  case(5)
    gx = -im * k(1)**5 + 20.0_dp * ga%alpha    * im * k(1)**3 &
                       - 60.0_dp * ga%alpha**2 * im * k(1)
  case(7)
    gx =  im * k(1)**7 -  42.0_dp * ga%alpha    * im * k(1)**5 &
                       + 420.0_dp * ga%alpha**2 * im * k(1)**3 &
                       - 840.0_dp * ga%alpha**3 * im * k(1)
  case default
    call die('eval_gaussian_fourier: too high angular momentum')
  end select

  ! Y direction
  select case(ga%ny)
  case(0)
    gy = 1.0_dp
  case(2)
    gy = -ky2 + 2.0_dp * ga%alpha
  case(4)
    gy = ky2**2 - 12.0_dp * ga%alpha * ky2 + 12.0_dp * ga%alpha**2
  case(6)
    gy = -ky2**3 + 30.0_dp * ga%alpha * ky2**2 - 180.0_dp * ga%alpha**2 * ky2 &
                         + 120.0_dp * ga%alpha**3
  case(8)
    gy = ky2**4 - 56.0_dp * ga%alpha * ky2**3 + 840.0_dp * ga%alpha**2 * ky2**2 &
                                - 3360.0_dp * ga%alpha**3 * ky2 + 1680.0_dp * ga%alpha**4
  case(1)
    gy = -im * k(2)
  case(3)
    gy =  im * k(2)**3 -  6.0_dp * ga%alpha    * im * k(2)
  case(5)
    gy = -im * k(2)**5 + 20.0_dp * ga%alpha    * im * k(2)**3 &
                       - 60.0_dp * ga%alpha**2 * im * k(2)
  case(7)
    gy = im * k(2)**7 - 42.0_dp * ga%alpha     * im * k(2)**5 &
                      + 420.0_dp * ga%alpha**2 * im * k(2)**3 &
                      - 840.0_dp * ga%alpha**3 * im * k(2)
  case default
    call die('eval_gaussian_fourier: too high angular momentum')
  end select

  ! Z direction
  select case(ga%nz)
  case(0)
    gz = 1.0_dp
  case(2)
    gz = -kz2 + 2.0_dp * ga%alpha
  case(4)
    gz = kz2**2 - 12.0_dp * ga%alpha * kz2 + 12.0_dp * ga%alpha**2
  case(6)
    gz = -kz2**3 + 30.0_dp * ga%alpha * kz2**2 - 180.0_dp * ga%alpha**2 * kz2 &
                         + 120.0_dp * ga%alpha**3
  case(8)
    gz = kz2**4 - 56.0_dp * ga%alpha * kz2**3 + 840.0_dp * ga%alpha**2 * kz2**2 &
                                - 3360.0_dp * ga%alpha**3 * kz2 + 1680.0_dp * ga%alpha**4
  case(1)
    gz = -im * k(3)
  case(3)
    gz =  im * k(3)**3 - 6.0_dp * ga%alpha     * im * k(3)
  case(5)
    gz = -im * k(3)**5 + 20.0_dp * ga%alpha    * im * k(3)**3 &
                       - 60.0_dp * ga%alpha**2 * im * k(3)
  case(7)
    gz =  im * k(3)**7 - 42.0_dp * ga%alpha    * im * k(3)**5 &
                      + 420.0_dp * ga%alpha**2 * im * k(3)**3 &
                      - 840.0_dp * ga%alpha**3 * im * k(3)
  case default
    call die('eval_gaussian_fourier: too high angular momentum')
  end select

  gaussian_fourier = gx * gy * gz * (pi/ga%alpha)**(1.5_dp) / (2*ga%alpha)**(ga%nx+ga%ny+ga%nz)   &
                          * EXP( -0.25_dp * k2 / ga%alpha ) &
                          * EXP( -im * DOT_PRODUCT(k,ga%x0) ) &
                          / ( 2.0_dp * pi )**3
  !
  ! normalize the cartesian gaussian
  gaussian_fourier = gaussian_fourier * ga%norm_factor

end function eval_gaussian_fourier


!=========================================================================
function eval_gaussian_grad(ga,x)
  implicit none

  type(gaussian),intent(in) :: ga
  real(dp),intent(in) :: x(3)
  real(dp) :: eval_gaussian_grad(3)
  !=====
  real(dp) :: dx(3),dx2
  real(dp) :: g_x,g_y,g_z,gp_x,gp_y,gp_z
  !=====

  dx(:) = x(:) - ga%x0(:)
  dx2 = SUM( dx(:)**2 )

  g_x = dx(1)**ga%nx
  g_y = dx(2)**ga%ny
  g_z = dx(3)**ga%nz

  gp_x = -2.0_dp * ga%alpha * dx(1)**(ga%nx+1)
  if(ga%nx>0) gp_x = gp_x + REAL(ga%nx,dp) * dx(1)**(ga%nx-1)
  gp_y = -2.0_dp * ga%alpha * dx(2)**(ga%ny+1)
  if(ga%ny>0) gp_y = gp_y + REAL(ga%ny,dp) * dx(2)**(ga%ny-1)
  gp_z = -2.0_dp * ga%alpha * dx(3)**(ga%nz+1)
  if(ga%nz>0) gp_z = gp_z + REAL(ga%nz,dp) * dx(3)**(ga%nz-1)

  eval_gaussian_grad(1) = gp_x * g_y  * g_z
  eval_gaussian_grad(2) = g_x  * gp_y * g_z
  eval_gaussian_grad(3) = g_x  * g_y  * gp_z

  !
  ! multiply by the common exponential factor
  ! and normalize the cartesian gaussian
  eval_gaussian_grad(:) = eval_gaussian_grad(:) * ga%norm_factor * EXP( - ga%alpha * dx2 )


end function eval_gaussian_grad


!=========================================================================
function eval_gaussian_lapl(ga,x)
  implicit none

  type(gaussian),intent(in) :: ga
  real(dp),intent(in) :: x(3)
  real(dp) :: eval_gaussian_lapl(3)
  !=====
  real(dp) :: dx(3),dx2
  real(dp) :: g_x,g_y,g_z,gpp_x,gpp_y,gpp_z
  !=====

  dx(:) = x(:) - ga%x0(:)
  dx2 = SUM( dx(:)**2 )

  g_x = dx(1)**ga%nx
  g_y = dx(2)**ga%ny
  g_z = dx(3)**ga%nz

  gpp_x = 4.0_dp * ga%alpha**2 * dx(1)**(ga%nx+2) - 2.0_dp * ga%alpha * REAL(ga%nx+1,dp) * dx(1)**ga%nx
  if(ga%nx>0) gpp_x = gpp_x - 2.0_dp * ga%alpha * REAL(ga%nx,dp) * dx(1)**(ga%nx)
  if(ga%nx>1) gpp_x = gpp_x + 2.0_dp * ga%alpha * REAL(ga%nx*(ga%nx-1),dp) * dx(1)**(ga%nx-1)

  gpp_y = 4.0_dp * ga%alpha**2 * dx(2)**(ga%ny+2) - 2.0_dp * ga%alpha * REAL(ga%ny+1,dp) * dx(2)**ga%ny
  if(ga%ny>0) gpp_y = gpp_y - 2.0_dp * ga%alpha * REAL(ga%ny,dp) * dx(2)**(ga%ny)
  if(ga%ny>1) gpp_y = gpp_y + 2.0_dp * ga%alpha * REAL(ga%ny*(ga%ny-1),dp) * dx(2)**(ga%ny-1)

  gpp_z = 4.0_dp * ga%alpha**2 * dx(3)**(ga%nz+2) - 2.0_dp * ga%alpha * REAL(ga%nz+1,dp) * dx(3)**ga%nz
  if(ga%nz>0) gpp_z = gpp_z - 2.0_dp * ga%alpha * REAL(ga%nz,dp) * dx(3)**(ga%nz)
  if(ga%nz>1) gpp_z = gpp_z + 2.0_dp * ga%alpha * REAL(ga%nz*(ga%nz-1),dp) * dx(3)**(ga%nz-1)

  eval_gaussian_lapl(1) = gpp_x * g_y   * g_z
  eval_gaussian_lapl(2) = g_x   * gpp_y * g_z
  eval_gaussian_lapl(3) = g_x   * g_y   * gpp_z

  !
  ! multiply by the common exponential factor
  ! and normalize the cartesian gaussian
  eval_gaussian_lapl(:) = eval_gaussian_lapl(:) * ga%norm_factor * EXP( - ga%alpha * dx2 )


end function eval_gaussian_lapl


!=========================================================================
subroutine print_gaussian(ga)
  implicit none

  type(gaussian),intent(in) :: ga
  !=====
  !=====

  write(stdout,*)
  write(stdout,*) '*** output information of a cartesian gaussian function ***'
  write(stdout,'(a30,5x,a1)')      'orbital momentum',ga%amc
  write(stdout,'(a30,2x,3(1x,i3))') 'momentum composition',ga%nx,ga%ny,ga%nz
  write(stdout,'(a30,(1x,f12.6))')  'alpha coeff',ga%alpha
  write(stdout,'(a30,3(1x,f12.6))') 'center of the function',ga%x0(:)
  write(stdout,*)

end subroutine print_gaussian


!=========================================================================
subroutine overlap_recurrence(ga,gb,s_ab)
  implicit none

  type(gaussian),intent(in) :: ga,gb
  real(dp),intent(out) :: s_ab
  !=====
  real(dp)             :: zeta_ab,ksi_ab,ab2,fact
  real(dp)             :: p(3),ap(3),bp(3)
  real(dp)             :: s_tmp_x(0:ga%nx,0:gb%nx)
  real(dp)             :: s_tmp_y(0:ga%ny,0:gb%ny)
  real(dp)             :: s_tmp_z(0:ga%nz,0:gb%nz)
  integer              :: ix,iy,iz
  integer              :: ixa,iya,iza
  integer              :: ixb,iyb,izb
  integer              :: ixap,iyap,izap
  integer              :: ixbp,iybp,izbp
  !=====

  ! Follow the notation of Obara and Saika, JCP  87 3963 (1986)
  zeta_ab = ga%alpha + gb%alpha
  ksi_ab   = ga%alpha * gb%alpha / zeta_ab
  ab2    = SUM( ( ga%x0(:)-gb%x0(:) )**2 )
  p(:)    = ( ga%alpha * ga%x0(:) + gb%alpha * gb%x0(:) ) / zeta_ab
  ap(:) =  p(:) - ga%x0(:)
  bp(:) =  p(:) - gb%x0(:)
  fact = 0.5_dp / zeta_ab


  !
  ! direction X
  !
  s_tmp_x(0,0) = (pi/zeta_ab)**1.5_dp * EXP( - ksi_ab * ab2 )

  do ix=1,ga%nx+gb%nx

    do ixa=0,MIN(ga%nx,ix)
      ixb=ix-ixa
      if(ixb>gb%nx) cycle

      if(ixa>0) then
        ixap=ixa-1
        s_tmp_x(ixap+1,ixb) = ap(1) * s_tmp_x(ixap,ixb)
        if(ixap>0)  s_tmp_x(ixap+1,ixb) = s_tmp_x(ixap+1,ixb) + fact * ixap * s_tmp_x(ixap-1,ixb)
        if(ixb>0)   s_tmp_x(ixap+1,ixb) = s_tmp_x(ixap+1,ixb) + fact * ixb  * s_tmp_x(ixap,ixb-1)
      else
        ixbp=ixb-1
        s_tmp_x(ixa,ixbp+1) = bp(1) * s_tmp_x(ixa,ixbp)
        if(ixbp>0) s_tmp_x(ixa,ixbp+1) = s_tmp_x(ixa,ixbp+1) + fact * ixbp * s_tmp_x(ixa,ixbp-1)
        if(ixa>0)  s_tmp_x(ixa,ixbp+1) = s_tmp_x(ixa,ixbp+1) + fact * ixa  * s_tmp_x(ixa-1,ixbp)
      endif

    enddo

  enddo

  !
  ! direction Y
  !
  s_tmp_y(0,0) = s_tmp_x(ga%nx,gb%nx)

  do iy=1,ga%ny+gb%ny

    do iya=0,MIN(ga%ny,iy)
      iyb=iy-iya
      if(iyb>gb%ny) cycle

      if(iya>0) then
        iyap=iya-1
        s_tmp_y(iyap+1,iyb) = ap(2) * s_tmp_y(iyap,iyb)
        if(iyap>0)  s_tmp_y(iyap+1,iyb) = s_tmp_y(iyap+1,iyb) + fact * iyap * s_tmp_y(iyap-1,iyb)
        if(iyb>0)   s_tmp_y(iyap+1,iyb) = s_tmp_y(iyap+1,iyb) + fact * iyb  * s_tmp_y(iyap,iyb-1)
      else
        iybp=iyb-1
        s_tmp_y(iya,iybp+1) = bp(2) * s_tmp_y(iya,iybp)
        if(iybp>0) s_tmp_y(iya,iybp+1) = s_tmp_y(iya,iybp+1) + fact * iybp * s_tmp_y(iya,iybp-1)
        if(iya>0)  s_tmp_y(iya,iybp+1) = s_tmp_y(iya,iybp+1) + fact * iya  * s_tmp_y(iya-1,iybp)
      endif

    enddo

  enddo

  !
  ! direction Z
  !
  s_tmp_z(0,0) = s_tmp_y(ga%ny,gb%ny)

  do iz=1,ga%nz+gb%nz

    do iza=0,MIN(ga%nz,iz)
      izb=iz-iza
      if(izb>gb%nz) cycle

      if(iza>0) then
        izap=iza-1
        s_tmp_z(izap+1,izb) = ap(3) * s_tmp_z(izap,izb)
        if(izap>0)  s_tmp_z(izap+1,izb) = s_tmp_z(izap+1,izb) + fact * izap * s_tmp_z(izap-1,izb)
        if(izb>0)   s_tmp_z(izap+1,izb) = s_tmp_z(izap+1,izb) + fact * izb  * s_tmp_z(izap,izb-1)
      else
        izbp=izb-1
        s_tmp_z(iza,izbp+1) = bp(3) * s_tmp_z(iza,izbp)
        if(izbp>0) s_tmp_z(iza,izbp+1) = s_tmp_z(iza,izbp+1) + fact * izbp * s_tmp_z(iza,izbp-1)
        if(iza>0)  s_tmp_z(iza,izbp+1) = s_tmp_z(iza,izbp+1) + fact * iza  * s_tmp_z(iza-1,izbp)
      endif

    enddo

  enddo


  s_ab = s_tmp_z(ga%nz,gb%nz) * ga%norm_factor * gb%norm_factor


end subroutine overlap_recurrence


!=========================================================================
subroutine kinetic_recurrence(ga,gb,k_ab)
  implicit none

  type(gaussian),intent(in)     :: ga,gb
  real(dp),intent(out)          :: k_ab
  !=====
  real(dp)             :: zeta_ab,ksi_ab,ab2,fact
  real(dp)             :: p(3),ap(3),bp(3)
  real(dp)             :: s_tmp_x(0:ga%nx,0:gb%nx)
  real(dp)             :: s_tmp_y(0:ga%ny,0:gb%ny)
  real(dp)             :: s_tmp_z(0:ga%nz,0:gb%nz)
  real(dp)             :: k_tmp_x(0:ga%nx,0:gb%nx)
  real(dp)             :: k_tmp_y(0:ga%ny,0:gb%ny)
  real(dp)             :: k_tmp_z(0:ga%nz,0:gb%nz)
  integer              :: ix,iy,iz
  integer              :: ixa,iya,iza
  integer              :: ixb,iyb,izb
  integer              :: ixap,iyap,izap
  integer              :: ixbp,iybp,izbp
  real(dp)             :: s_ab
  !=====

  ! Follow the notation of Obara and Saika, JCP  87 3963 (1986)
  zeta_ab = ga%alpha + gb%alpha
  ksi_ab   = ga%alpha * gb%alpha / zeta_ab
  ab2    = SUM( ( ga%x0(:)-gb%x0(:) )**2 )
  p(:)    = ( ga%alpha * ga%x0(:) + gb%alpha * gb%x0(:) ) / zeta_ab
  ap(:) =  p(:) - ga%x0(:)
  bp(:) =  p(:) - gb%x0(:)
  fact = 0.5_dp / zeta_ab


  !
  ! direction X
  !
  s_tmp_x(0,0) = (pi/zeta_ab)**1.5_dp * EXP( - ksi_ab * ab2 )
  k_tmp_x(0,0) = ksi_ab * ( 3.0_dp - 2.0_dp * ksi_ab * ab2 ) * s_tmp_x(0,0)

  do ix=1,ga%nx+gb%nx

    do ixa=0,MIN(ga%nx,ix)
      ixb=ix-ixa
      if(ixb>gb%nx) cycle

      if(ixa>0) then
        ixap=ixa-1

        s_tmp_x(ixap+1,ixb) = ap(1) * s_tmp_x(ixap,ixb)
        if(ixap>0)  s_tmp_x(ixap+1,ixb) = s_tmp_x(ixap+1,ixb) + fact * ixap * s_tmp_x(ixap-1,ixb)
        if(ixb>0)   s_tmp_x(ixap+1,ixb) = s_tmp_x(ixap+1,ixb) + fact * ixb  * s_tmp_x(ixap,ixb-1)

        k_tmp_x(ixap+1,ixb) = ap(1) * k_tmp_x(ixap,ixb) + 2.0_dp * ksi_ab * s_tmp_x(ixap+1,ixb)
        if(ixap>0)  k_tmp_x(ixap+1,ixb) = k_tmp_x(ixap+1,ixb) + fact * ixap * k_tmp_x(ixap-1,ixb) &
                     -  ksi_ab / ga%alpha * ixap * s_tmp_x(ixap-1,ixb)
        if(ixb>0)   k_tmp_x(ixap+1,ixb) = k_tmp_x(ixap+1,ixb) + fact * ixb  * k_tmp_x(ixap,ixb-1)

      else
        ixbp=ixb-1

        s_tmp_x(ixa,ixbp+1) = bp(1) * s_tmp_x(ixa,ixbp)
        if(ixbp>0) s_tmp_x(ixa,ixbp+1) = s_tmp_x(ixa,ixbp+1) + fact * ixbp * s_tmp_x(ixa,ixbp-1)
        if(ixa>0)  s_tmp_x(ixa,ixbp+1) = s_tmp_x(ixa,ixbp+1) + fact * ixa  * s_tmp_x(ixa-1,ixbp)

        k_tmp_x(ixa,ixbp+1) = bp(1) * k_tmp_x(ixa,ixbp) +  2.0_dp * ksi_ab * s_tmp_x(ixa,ixbp+1)
        if(ixbp>0) k_tmp_x(ixa,ixbp+1) = k_tmp_x(ixa,ixbp+1) + fact * ixbp * k_tmp_x(ixa,ixbp-1) &
                     -  ksi_ab / gb%alpha * ixbp * s_tmp_x(ixa,ixbp-1)
        if(ixa>0)  k_tmp_x(ixa,ixbp+1) = k_tmp_x(ixa,ixbp+1) + fact * ixa  * k_tmp_x(ixa-1,ixbp)

      endif

    enddo

  enddo

  !
  ! direction Y
  !
  s_tmp_y(0,0) = s_tmp_x(ga%nx,gb%nx)
  k_tmp_y(0,0) = k_tmp_x(ga%nx,gb%nx)

  do iy=1,ga%ny+gb%ny

    do iya=0,MIN(ga%ny,iy)
      iyb=iy-iya
      if(iyb>gb%ny) cycle

      if(iya>0) then
        iyap=iya-1

        s_tmp_y(iyap+1,iyb) = ap(2) * s_tmp_y(iyap,iyb)
        if(iyap>0)  s_tmp_y(iyap+1,iyb) = s_tmp_y(iyap+1,iyb) + fact * iyap * s_tmp_y(iyap-1,iyb)
        if(iyb>0)   s_tmp_y(iyap+1,iyb) = s_tmp_y(iyap+1,iyb) + fact * iyb  * s_tmp_y(iyap,iyb-1)

        k_tmp_y(iyap+1,iyb) = ap(2) * k_tmp_y(iyap,iyb) + 2.0_dp * ksi_ab * s_tmp_y(iyap+1,iyb)
        if(iyap>0)  k_tmp_y(iyap+1,iyb) = k_tmp_y(iyap+1,iyb) + fact * iyap * k_tmp_y(iyap-1,iyb) &
                     -  ksi_ab / ga%alpha * iyap * s_tmp_y(iyap-1,iyb)
        if(iyb>0)   k_tmp_y(iyap+1,iyb) = k_tmp_y(iyap+1,iyb) + fact * iyb  * k_tmp_y(iyap,iyb-1)

      else
        iybp=iyb-1

        s_tmp_y(iya,iybp+1) = bp(2) * s_tmp_y(iya,iybp)
        if(iybp>0) s_tmp_y(iya,iybp+1) = s_tmp_y(iya,iybp+1) + fact * iybp * s_tmp_y(iya,iybp-1)
        if(iya>0)  s_tmp_y(iya,iybp+1) = s_tmp_y(iya,iybp+1) + fact * iya  * s_tmp_y(iya-1,iybp)

        k_tmp_y(iya,iybp+1) = bp(2) * k_tmp_y(iya,iybp) +  2.0_dp * ksi_ab * s_tmp_y(iya,iybp+1)
        if(iybp>0) k_tmp_y(iya,iybp+1) = k_tmp_y(iya,iybp+1) + fact * iybp * k_tmp_y(iya,iybp-1) &
                     -  ksi_ab / gb%alpha * iybp * s_tmp_y(iya,iybp-1)
        if(iya>0)  k_tmp_y(iya,iybp+1) = k_tmp_y(iya,iybp+1) + fact * iya  * k_tmp_y(iya-1,iybp)

      endif

    enddo

  enddo

  !
  ! direction Z
  !
  s_tmp_z(0,0) = s_tmp_y(ga%ny,gb%ny)
  k_tmp_z(0,0) = k_tmp_y(ga%ny,gb%ny)

  do iz=1,ga%nz+gb%nz

    do iza=0,MIN(ga%nz,iz)
      izb=iz-iza
      if(izb>gb%nz) cycle

      if(iza>0) then
        izap=iza-1

        s_tmp_z(izap+1,izb) = ap(3) * s_tmp_z(izap,izb)
        if(izap>0)  s_tmp_z(izap+1,izb) = s_tmp_z(izap+1,izb) + fact * izap * s_tmp_z(izap-1,izb)
        if(izb>0)   s_tmp_z(izap+1,izb) = s_tmp_z(izap+1,izb) + fact * izb  * s_tmp_z(izap,izb-1)

        k_tmp_z(izap+1,izb) = ap(3) * k_tmp_z(izap,izb) + 2.0_dp * ksi_ab * s_tmp_z(izap+1,izb)
        if(izap>0)  k_tmp_z(izap+1,izb) = k_tmp_z(izap+1,izb) + fact * izap * k_tmp_z(izap-1,izb) &
                     -  ksi_ab / ga%alpha * izap * s_tmp_z(izap-1,izb)
        if(izb>0)   k_tmp_z(izap+1,izb) = k_tmp_z(izap+1,izb) + fact * izb  * k_tmp_z(izap,izb-1)

      else
        izbp=izb-1

        s_tmp_z(iza,izbp+1) = bp(3) * s_tmp_z(iza,izbp)
        if(izbp>0) s_tmp_z(iza,izbp+1) = s_tmp_z(iza,izbp+1) + fact * izbp * s_tmp_z(iza,izbp-1)
        if(iza>0)  s_tmp_z(iza,izbp+1) = s_tmp_z(iza,izbp+1) + fact * iza  * s_tmp_z(iza-1,izbp)

        k_tmp_z(iza,izbp+1) = bp(3) * k_tmp_z(iza,izbp) +  2.0_dp * ksi_ab * s_tmp_z(iza,izbp+1)
        if(izbp>0) k_tmp_z(iza,izbp+1) = k_tmp_z(iza,izbp+1) + fact * izbp * k_tmp_z(iza,izbp-1) &
                     -  ksi_ab / gb%alpha * izbp * s_tmp_z(iza,izbp-1)
        if(iza>0)  k_tmp_z(iza,izbp+1) = k_tmp_z(iza,izbp+1) + fact * iza  * k_tmp_z(iza-1,izbp)

      endif

    enddo

  enddo


  !
  ! overlap is a by-product
  !
  s_ab = s_tmp_z(ga%nz,gb%nz) * ga%norm_factor * gb%norm_factor

  !
  ! final result
  !
  k_ab = k_tmp_z(ga%nz,gb%nz) * ga%norm_factor * gb%norm_factor


end subroutine kinetic_recurrence


!=========================================================================
subroutine nucleus_recurrence(zatom,c,ga,gb,v_ab)
  implicit none

  real(dp),intent(in)       :: zatom,c(3)
  type(gaussian),intent(in) :: ga,gb
  real(dp),intent(out)      :: v_ab
  !=====
  real(dp)             :: zeta_ab,ksi_ab,ab2,fact
  real(dp)             :: p(3),ap(3),bp(3),cp(3)
  real(dp)             :: v_tmp_x_m(0:ga%nx,0:gb%nx)
  real(dp)             :: v_tmp_y_m(0:ga%ny,0:gb%ny)
  real(dp)             :: v_tmp_z_m(0:ga%nz,0:gb%nz)
  real(dp)             :: v_tmp_x_mp1(0:ga%nx,0:gb%nx)
  real(dp)             :: v_tmp_y_mp1(0:ga%ny,0:gb%ny)
  real(dp)             :: v_tmp_z_mp1(0:ga%nz,0:gb%nz)
  integer              :: ix,iy,iz
  integer              :: ixa,iya,iza
  integer              :: ixb,iyb,izb
  integer              :: ixap,iyap,izap
  integer              :: ixbp,iybp,izbp
  integer              :: ixam,iyam,izam
  integer              :: ixbm,iybm,izbm
  integer              :: mm
  real(dp),allocatable :: fmu(:)
  real(dp)             :: bigu
  !=====

  ! Follow the notation of Obara and Saika, JCP  87 3963 (1986)
  zeta_ab = ga%alpha + gb%alpha
  ksi_ab   = ga%alpha * gb%alpha / zeta_ab
  ab2    = SUM( ( ga%x0(:)-gb%x0(:) )**2 )
  p(:)    = ( ga%alpha * ga%x0(:) + gb%alpha * gb%x0(:) ) / zeta_ab
  ap(:) =  p(:) - ga%x0(:)
  bp(:) =  p(:) - gb%x0(:)
  cp(:) =  p(:) - c(:)
  fact  = 0.5_dp / zeta_ab
  bigu  = zeta_ab * SUM( cp(:)**2 )


  v_tmp_x_mp1(:,:) =  0.0_dp
  v_tmp_y_mp1(:,:) =  0.0_dp
  v_tmp_z_mp1(:,:) =  0.0_dp

  do mm=ga%am+gb%am,0,-1
    allocate(fmu(0:mm))
    call boys_function_c(fmu,mm,bigu)

    !
    ! direction X
    !
    v_tmp_x_m(0,0) = 2.0 * fmu(mm) * (pi/zeta_ab) * EXP( - ksi_ab * ab2 )

    ixam = 0
    ixbm = 0
    do ix=1,ga%am+gb%am-mm
      do ixa=0,MIN(ga%nx,ix)
        ixb = ix - ixa
        if( ixb > gb%nx ) cycle
        ixam = MAX(ixam,ixa)
        ixbm = MAX(ixbm,ixb)

        if( ixa > 0 ) then
          ixap = ixa - 1
          v_tmp_x_m(ixap+1,ixb) = ap(1) * v_tmp_x_m(ixap,ixb) - cp(1) * v_tmp_x_mp1(ixap,ixb)
          if( ixap > 0 ) v_tmp_x_m(ixap+1,ixb) = v_tmp_x_m(ixap+1,ixb) &
                                                + fact * ixap * ( v_tmp_x_m(ixap-1,ixb) -  v_tmp_x_mp1(ixap-1,ixb) )
          if( ixb > 0 )  v_tmp_x_m(ixap+1,ixb) = v_tmp_x_m(ixap+1,ixb) &
                                                + fact * ixb  * ( v_tmp_x_m(ixap,ixb-1) -  v_tmp_x_mp1(ixap,ixb-1) )
        else
          ixbp = ixb - 1
          v_tmp_x_m(ixa,ixbp+1) = bp(1) * v_tmp_x_m(ixa,ixbp) - cp(1) * v_tmp_x_mp1(ixa,ixbp)
          if( ixbp > 0 ) v_tmp_x_m(ixa,ixbp+1) = v_tmp_x_m(ixa,ixbp+1) &
                                                + fact * ixbp * ( v_tmp_x_m(ixa,ixbp-1) -  v_tmp_x_mp1(ixa,ixbp-1) )
          if( ixa > 0 )  v_tmp_x_m(ixa,ixbp+1) = v_tmp_x_m(ixa,ixbp+1) &
                                                + fact * ixa  * ( v_tmp_x_m(ixa-1,ixbp) -  v_tmp_x_mp1(ixa-1,ixbp) )
        endif

      enddo

    enddo

    !
    ! direction Y
    !
    v_tmp_y_m(0,0) = v_tmp_x_m(ixam,ixbm)

    iyam = 0
    iybm = 0
    do iy=1,ga%am+gb%am-mm
      do iya=0,MIN(ga%ny,iy)
        iyb = iy - iya
        if( iyb > gb%ny ) cycle
        iyam = MAX(iyam,iya)
        iybm = MAX(iybm,iyb)

        if( iya > 0 ) then
          iyap = iya - 1
          v_tmp_y_m(iyap+1,iyb) = ap(2) * v_tmp_y_m(iyap,iyb) - cp(2) * v_tmp_y_mp1(iyap,iyb)
          if( iyap > 0 ) v_tmp_y_m(iyap+1,iyb) = v_tmp_y_m(iyap+1,iyb) &
                                                + fact * iyap * ( v_tmp_y_m(iyap-1,iyb) -  v_tmp_y_mp1(iyap-1,iyb) )
          if( iyb > 0 )  v_tmp_y_m(iyap+1,iyb) = v_tmp_y_m(iyap+1,iyb) &
                                                + fact * iyb  * ( v_tmp_y_m(iyap,iyb-1) -  v_tmp_y_mp1(iyap,iyb-1) )
        else
          iybp = iyb - 1
          v_tmp_y_m(iya,iybp+1) = bp(2) * v_tmp_y_m(iya,iybp) - cp(2) * v_tmp_y_mp1(iya,iybp)
          if( iybp > 0 ) v_tmp_y_m(iya,iybp+1) = v_tmp_y_m(iya,iybp+1) &
                                                + fact * iybp * ( v_tmp_y_m(iya,iybp-1) -  v_tmp_y_mp1(iya,iybp-1) )
          if( iya > 0 )  v_tmp_y_m(iya,iybp+1) = v_tmp_y_m(iya,iybp+1) &
                                                + fact * iya  * ( v_tmp_y_m(iya-1,iybp) -  v_tmp_y_mp1(iya-1,iybp) )
        endif

      enddo

    enddo

    !
    ! direction Z
    !
    v_tmp_z_m(0,0) = v_tmp_y_m(iyam,iybm)

    izam = 0
    izbm = 0
    do iz=1,ga%am+gb%am-mm
      do iza=0,MIN(ga%nz,iz)
        izb = iz - iza
        if( izb > gb%nz ) cycle
        izam = MAX(izam,iza)
        izbm = MAX(izbm,izb)

        if( iza > 0 ) then
          izap = iza - 1
          v_tmp_z_m(izap+1,izb) = ap(3) * v_tmp_z_m(izap,izb) - cp(3) * v_tmp_z_mp1(izap,izb)
          if( izap > 0 ) v_tmp_z_m(izap+1,izb) = v_tmp_z_m(izap+1,izb) &
                                                + fact * izap * ( v_tmp_z_m(izap-1,izb) -  v_tmp_z_mp1(izap-1,izb) )
          if( izb > 0 )  v_tmp_z_m(izap+1,izb) = v_tmp_z_m(izap+1,izb) &
                                                + fact * izb  * ( v_tmp_z_m(izap,izb-1) -  v_tmp_z_mp1(izap,izb-1) )
        else
          izbp = izb - 1
          v_tmp_z_m(iza,izbp+1) = bp(3) * v_tmp_z_m(iza,izbp) - cp(3) * v_tmp_z_mp1(iza,izbp)
          if( izbp > 0 ) v_tmp_z_m(iza,izbp+1) = v_tmp_z_m(iza,izbp+1) &
                                                + fact * izbp * ( v_tmp_z_m(iza,izbp-1) -  v_tmp_z_mp1(iza,izbp-1) )
          if( iza > 0 )  v_tmp_z_m(iza,izbp+1) = v_tmp_z_m(iza,izbp+1) &
                                                + fact * iza  * ( v_tmp_z_m(iza-1,izbp) -  v_tmp_z_mp1(iza-1,izbp) )
        endif

      enddo

    enddo

    v_tmp_x_mp1(0:ixam,0:ixbm) = v_tmp_x_m(0:ixam,0:ixbm)
    v_tmp_y_mp1(0:iyam,0:iybm) = v_tmp_y_m(0:iyam,0:iybm)
    v_tmp_z_mp1(0:izam,0:izbm) = v_tmp_z_m(0:izam,0:izbm)

    deallocate(fmu)
  enddo


  v_ab = - zatom * v_tmp_z_m(izam,izbm) * ga%norm_factor * gb%norm_factor


end subroutine nucleus_recurrence


!=========================================================================
subroutine evaluate_gos(ga,gb,qvec,gos_ab)
  implicit none

  type(gaussian),intent(in) :: ga,gb
  real(dp),intent(in)       :: qvec(3)
  complex(dp),intent(out)   :: gos_ab
  !=====
  complex(dp) :: intx,inty,intz,qq(3),factor
  real(dp) :: pp,a2,b2
  integer :: ii,jj,kk
  !=====

  !
  ! Calculate < \phi_\alpha | e^{iq.r} | \phi_beta >
  !
  ! Formula from EH Mortensen et al.  Nuclear Instruments and Methods in Physics Research B69 24 (1992)
  !
  pp = ga%alpha + gb%alpha
  a2 = SUM( ga%x0(:)**2 )
  b2 = SUM( gb%x0(:)**2 )

  qq(:) = ga%alpha * ga%x0(:) + gb%alpha * gb%x0(:) + 0.5_dp * im * qvec(:)


  factor = ( pi / pp )**1.5_dp * EXP( SUM(qq(:)**2) / pp - ga%alpha * a2 - gb%alpha * b2 )


  intx = 0.0_dp
  do ii=0,ga%nx+gb%nx
    intx = intx + g(ii,pp,qq(1)) * f(ii,ga%nx,gb%nx,-ga%x0(1),-gb%x0(1))
  enddo
  inty = 0.0_dp
  do ii=0,ga%ny+gb%ny
    inty = inty + g(ii,pp,qq(2)) * f(ii,ga%ny,gb%ny,-ga%x0(2),-gb%x0(2))
  enddo
  intz = 0.0_dp
  do ii=0,ga%nz+gb%nz
    intz = intz + g(ii,pp,qq(3)) * f(ii,ga%nz,gb%nz,-ga%x0(3),-gb%x0(3))
  enddo

  gos_ab = factor * intx * inty * intz * ga%norm_factor * gb%norm_factor

contains

function g(ii,p,q)
  integer,intent(in) :: ii
  real(dp),intent(in) :: p
  complex(dp),intent(in) :: q
  complex(dp) :: g
  !=====
  integer :: ik
  !=====

  select case(ii)
  case(0)
    g = 1.0_dp
  case(1)
    g = q / p
  case(2)
    g =   1.0_dp / ( 2.0_dp * p ) + ( q / p )**2
  case(3)
    g = ( 3.0_dp * q ) / (2.0_dp * p**2 ) + ( q / p )**3
  case(4)
    g = ( 3.0_dp     ) / (4.0_dp * p**2 ) + ( 3.0_dp * q**2 ) / ( p**3 ) + ( q / p )**4
  case(5)
    g = ( 15.0_dp * q) / (4.0_dp * p**3 ) + ( 5.0_dp * q**3 ) / ( p**4 ) + ( q / p )**5
  case(6)
    g = ( 15.0_dp    ) / (8.0_dp * p**3 ) + ( 45.0_dp * q**2 ) / ( 4.0_dp * p**4 ) &
              + ( 15.0_dp * q**4 ) / ( 2.0_dp * p**5 ) + ( q / p )**6
  case default
    ! If not precalculated, then do it for real
    ! Remember Gamma(n+1) = n!
    g = (0.0_dp,0.0_dp)
    do ik=0,ii/2
      g = g + ( 0.25_dp * p )**ik * q**( ii -2 * ik )  / ( Gamma(REAL(ii+1-2*ik,dp)) * Gamma(REAL(ik+1,dp)) )
    enddo
    g = g * Gamma(REAL(ii+1,dp)) / p**ii
  end select

end function g

function f(j,l,m,a,b)
  integer,intent(in) :: j,l,m
  real(dp),intent(in) :: a,b
  real(dp) :: f
  !=====
  integer :: i
  !=====

  f = 0.0_dp
  do i=MAX(0,j-m),MIN(j,l)
    f = f + cnk(l,i) * cnk(m,j-i) * a**(l-i) * b**(m+i-j)
  enddo

 end function f

end subroutine evaluate_gos


!=========================================================================
end module m_gaussian
