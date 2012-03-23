!=========================================================================
module m_gaussian
 use m_definitions

 ! type containing all the information for one unnormalized cartesian gaussian
 ! x**nx * y**ny * z**nz * exp( - alpha * ( x**2 + y**2 + z**2 ) )
 type gaussian
   integer          :: am
   character(len=1) :: amc
   integer          :: nx,ny,nz
   real(dp)         :: alpha
   real(dp)         :: x0(3)         ! center of the gaussian *** NOT USED (YET) ***
   real(dp)         :: norm_factor   ! normalization factor for the gaussian squared
 end type gaussian


contains

!=========================================================================
subroutine init_gaussian(nx,ny,nz,alpha,ga)
 implicit none
 integer,intent(in) :: nx,ny,nz
 real(dp),intent(in) :: alpha
 type(gaussian),intent(inout) :: ga
!=====

 ga%nx = nx
 ga%ny = ny
 ga%nz = nz
 ga%am = nx + ny + nz
 ga%amc = orbital_momentum_name(ga%am)

 ga%alpha = alpha

 ga%norm_factor = ( 2.0_dp / pi )**0.75_dp &
                 * 2.0_dp**ga%am * ga%alpha**( 0.25_dp * ( 2.0_dp*ga%am + 3.0_dp ) ) &
                 / SQRT( REAL( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) , dp ) )

 ga%x0(:) = 0.0_dp

end subroutine init_gaussian

!=========================================================================
subroutine init_gaussian_general(nx,ny,nz,alpha,x0,ga)
 implicit none
 integer,intent(in) :: nx,ny,nz
 real(dp),intent(in) :: alpha,x0(3)
 type(gaussian),intent(inout) :: ga
!=====

 ga%nx = nx
 ga%ny = ny
 ga%nz = nz
 ga%am = nx + ny + nz
 ga%amc = orbital_momentum_name(ga%am)

 ga%alpha = alpha

 ga%norm_factor = ( 2.0_dp / pi )**0.75_dp &
                 * 2.0_dp**ga%am * ga%alpha**( 0.25_dp * ( 2.0_dp*ga%am + 3.0_dp ) ) &
                 / SQRT( REAL( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) , dp ) )

 ga%x0(:) = x0(:)

end subroutine init_gaussian_general

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

 eval_gaussian = dx(1)**ga%nx &
                *dx(2)**ga%ny &
                *dx(3)**ga%nz &
                *exp( - ga%alpha * dx2 )

 !
 ! normalize the cartesian gaussian
 eval_gaussian = eval_gaussian * ga%norm_factor

end function eval_gaussian

!=========================================================================
function eval_gaussian_derivative(ga,x)
 implicit none
 type(gaussian),intent(in) :: ga
 real(dp),intent(in) :: x(3)
 real(dp) :: eval_gaussian_derivative(3)
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

 eval_gaussian_derivative(1) = gp_x * g_y  * g_z 
 eval_gaussian_derivative(2) = g_x  * gp_y * g_z 
 eval_gaussian_derivative(3) = g_x  * g_y  * gp_z

 !
 ! multiply by the common exponential factor
 ! and normalize the cartesian gaussian
 eval_gaussian_derivative(:) = eval_gaussian_derivative(:) * ga%norm_factor * exp( - ga%alpha * dx2 )


end function eval_gaussian_derivative

!=========================================================================
subroutine print_gaussian(ga)
 implicit none
 type(gaussian),intent(in) :: ga
!=====

 write(*,*)
 write(*,*) '*** output information of a cartesian gaussian function ***'
 write(*,'(a30,5x,a1)')      'orbital momentum',ga%amc
 write(*,'(a30,2x,3(x,i3))') 'momentum composition',ga%nx,ga%ny,ga%nz
 write(*,'(a30,(x,f12.6))')  'alpha coeff',ga%alpha
 write(*,'(a30,3(x,f12.6))') 'center of the function',ga%x0(:)
 write(*,*)

end subroutine print_gaussian

!=========================================================================
subroutine product_gaussian(ga,gb,gprod)
 implicit none
 type(gaussian),intent(in) :: ga,gb
 type(gaussian),intent(out) :: gprod
!=====

 if( ANY( ABS(ga%x0(:) - gb%x0(:)) > 1.d-6 ) ) stop'different positions not implemented'

 call init_gaussian_general(ga%nx+gb%nx,ga%ny+gb%ny,ga%nz+gb%nz,ga%alpha+gb%alpha,ga%x0,gprod)

end subroutine product_gaussian

!=========================================================================
subroutine overlap_normalized(ga,gb,s_ab)
 implicit none
 type(gaussian),intent(in) :: ga,gb
 real(dp),intent(out) :: s_ab
!=====

 call overlap(ga,gb,s_ab)
 s_ab = s_ab * ga%norm_factor * gb%norm_factor

end subroutine overlap_normalized

!=========================================================================
subroutine overlap(ga,gb,s_ab)
 implicit none
 type(gaussian),intent(in) :: ga,gb
 real(dp),intent(out) :: s_ab
!=====
 type(gaussian) :: gprod
!=====

 if( ANY( ABS(ga%x0(:) - gb%x0(:)) > 1.d-6 ) ) stop'different positions not implemented'

 call product_gaussian(ga,gb,gprod)

 s_ab = integral_1d(gprod%alpha,gprod%nx) &
       *integral_1d(gprod%alpha,gprod%ny) &
       *integral_1d(gprod%alpha,gprod%nz)

contains

  function integral_1d(alpha,nx)
   real(dp),intent(in) :: alpha
   integer,intent(in) :: nx
   real(dp) :: integral_1d
  
   !
   ! formula obtained from Wolfram online integrator!
   !
   if( mod(nx,2) == 1 ) then
     integral_1d=0.0_dp
   else
     integral_1d = alpha**( -0.5_dp * ( nx + 1 ) ) * gamma_function( 0.5_dp * ( nx + 1 ) )
   endif
  
  end function integral_1d
end subroutine overlap

!=========================================================================
subroutine kinetic_gaussian(ga,gb,kin_me)
 implicit none
 type(gaussian),intent(in) :: ga,gb
 real(dp),intent(out) :: kin_me
!=====
 type(gaussian) :: gbtmp
 integer :: nx_t,ny_t,nz_t
 real(dp) :: beta,rtmp,rtmpx,rtmpy,rtmpz
!=====

 beta = gb%alpha

 !
 ! the kinetic energy matrix element is a sum of 9 terms
 call overlap(ga,gb,rtmp)
 rtmpx = -2.0_dp * beta * ( 2.0 * gb%nx + 1 ) * rtmp
 rtmpy = -2.0_dp * beta * ( 2.0 * gb%ny + 1 ) * rtmp
 rtmpz = -2.0_dp * beta * ( 2.0 * gb%nz + 1 ) * rtmp

 gbtmp=gb
 gbtmp%nx=gbtmp%nx+2
 call overlap(ga,gbtmp,rtmp)
 rtmpx = rtmpx + 4.0_dp * beta**2 * rtmp

 gbtmp=gb
 gbtmp%ny=gbtmp%ny+2
 call overlap(ga,gbtmp,rtmp)
 rtmpy = rtmpy + 4.0_dp * beta**2 * rtmp

 gbtmp=gb
 gbtmp%nz=gbtmp%nz+2
 call overlap(ga,gbtmp,rtmp)
 rtmpz = rtmpz + 4.0_dp * beta**2 * rtmp

 if(gb%nx>=2) then
   gbtmp=gb
   gbtmp%nx=gbtmp%nx-2
   call overlap(ga,gbtmp,rtmp)
   rtmpx = rtmpx + gb%nx * ( gb%nx - 1 ) * rtmp
 endif

 if(gb%ny>=2) then
   gbtmp=gb
   gbtmp%ny=gbtmp%ny-2
   call overlap(ga,gbtmp,rtmp)
   rtmpy = rtmpy + gb%ny * ( gb%ny - 1 ) * rtmp
 endif

 if(gb%nz>=2) then
   gbtmp=gb
   gbtmp%nz=gbtmp%nz-2
   call overlap(ga,gbtmp,rtmp)
   rtmpz = rtmpz + gb%nz * ( gb%nz - 1 ) * rtmp
 endif

 !
 ! sum up the nine terms
 kin_me = -0.5_dp * ( rtmpx + rtmpy + rtmpz )
 !
 ! normalize everything
 kin_me = kin_me * ga%norm_factor * gb%norm_factor

end subroutine kinetic_gaussian

!=========================================================================
subroutine nucleus_pot_gaussian(ga,gb,zatom,pot_me)
 use m_tools, only: coeffs_gausslegint
 implicit none
 type(gaussian),intent(in) :: ga,gb
 real(dp),intent(in) :: zatom
 real(dp),intent(out) :: pot_me
!====
 integer,parameter :: NSAMP=100
 type(gaussian) :: gprod,gu
 integer :: isamp
 real(dp) :: u(NSAMP),wu(NSAMP)
 real(dp) :: u2,integral,rtmp
!=====

 !
 ! poor man version of the nucleus - electron integral
 ! only works when the gaussian and the atoms are located in the same place
 !
 ! based on the identity (found in Obara and Saika):
 ! 1 / |r_1 - r_2| = 2 / sqrt(pi) \int_0^{+\infty} du exp( - u^2 (r_1-r_2)^2 )
 !

 !
 ! set up the gauss-legendre integration scheme from 0 to infinity
 !
 ! first get the coefficient from 0 to 1
 call coeffs_gausslegint(0.0_dp,1.0_dp,u,wu,NSAMP)
 ! second apply variable changei
 !  u' =  u / (1-u)
 ! wu' = wu / (1-u)^2
 wu(:) = wu(:) / ( 1 - u(:) )**2
 u(:)  =  u(:) / ( 1 - u(:) )

 call product_gaussian(ga,gb,gprod)

 integral = 0.0_dp
 do isamp=1,NSAMP

   u2 = u(isamp)**2

   call init_gaussian(0,0,0,u2,gu)
   call overlap(gprod,gu,rtmp)
   integral = integral + wu(isamp) * rtmp

 enddo
 integral = integral * 2.0_dp / SQRT(pi)

 pot_me = - zatom * integral

 pot_me = pot_me * ga%norm_factor * gb%norm_factor

end subroutine nucleus_pot_gaussian

!=========================================================================
function orbital_momentum_name(am)
 integer,intent(in) :: am
 character(len=1) :: orbital_momentum_name
!=====

 select case(am)
 case(0)
   orbital_momentum_name='s'
 case(1)
   orbital_momentum_name='p'
 case(2)
   orbital_momentum_name='d'
 case(3)
   orbital_momentum_name='f'
 case(4)
   orbital_momentum_name='g'
 case(5)
   orbital_momentum_name='h'
 case(6)
   orbital_momentum_name='i'
 case(7)
   orbital_momentum_name='k'
 case(8)
   orbital_momentum_name='l'
 case(9)
   orbital_momentum_name='m'
 case(10)
   orbital_momentum_name='n'
 case(11)
   orbital_momentum_name='o'
 case(12)
   orbital_momentum_name='q'
 case default
   orbital_momentum_name='tmp'
!   stop'angular momentum not implemented'
 end select

end function orbital_momentum_name


#ifdef MOLECULES
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
!=====
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
 use m_tools, only: coeffs_gausslegint
 implicit none
 real(dp),intent(in)       :: zatom,c(3)
 type(gaussian),intent(in) :: ga,gb
 real(dp),intent(out)      :: v_ab
!=====
 integer,parameter    :: NSAMP=20
 real(dp)             :: zeta_ab,ksi_ab,ab2,fact
 real(dp)             :: p(3),ap(3),bp(3),cp(3)
 real(dp)             :: v_tmp_x(0:ga%nx,0:gb%nx)
 real(dp)             :: v_tmp_y(0:ga%ny,0:gb%ny)
 real(dp)             :: v_tmp_z(0:ga%nz,0:gb%nz)
 integer              :: ix,iy,iz
 integer              :: ixa,iya,iza
 integer              :: ixb,iyb,izb
 integer              :: ixap,iyap,izap
 integer              :: ixbp,iybp,izbp
 integer              :: isamp
 real(dp)             :: u(NSAMP),wu(NSAMP)
 real(dp)             :: u2,fact2,bigu
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

 !
 ! based on the identity (found in Obara and Saika):
 ! 1 / |r_1 - r_2| = 2 / sqrt(pi) \int_0^{+\infty} du exp( - u^2 (r_1-r_2)^2 )
 !

 !
 ! set up the gauss-legendre integration scheme from 0 to infinity
 !
 ! first get the coefficient from 0 to 1
 call coeffs_gausslegint(0.0_dp,1.0_dp,u,wu,NSAMP)
 ! second apply variable changei
 !  u' =  u / (1-u)
 ! wu' = wu / (1-u)^2
 wu(:) = wu(:) / ( 1 - u(:) )**2
 u(:)  =  u(:) / ( 1 - u(:) )

 v_ab = 0.0_dp
 do isamp=1,NSAMP

   u2 = u(isamp)**2
   fact2 = u2 / ( zeta_ab + u2 )



   !
   ! direction X
   !
!   v_tmp_x(0,0) = 2.0 * Fm(u) * (pi/zeta_ab) * EXP( - ksi_ab * ab2 )
   stop'FIXME implement Fm(u)'
  
   do ix=1,ga%nx+gb%nx
  
     do ixa=0,MIN(ga%nx,ix)
       ixb=ix-ixa
       if(ixb>gb%nx) cycle
  
       if(ixa>0) then
         ixap=ixa-1
         v_tmp_x(ixap+1,ixb) = ap(1) * v_tmp_x(ixap,ixb)
         if(ixap>0)  v_tmp_x(ixap+1,ixb) = v_tmp_x(ixap+1,ixb) + fact * ixap * v_tmp_x(ixap-1,ixb)
         if(ixb>0)   v_tmp_x(ixap+1,ixb) = v_tmp_x(ixap+1,ixb) + fact * ixb  * v_tmp_x(ixap,ixb-1)
       else
         ixbp=ixb-1
         v_tmp_x(ixa,ixbp+1) = bp(1) * v_tmp_x(ixa,ixbp)
         if(ixbp>0) v_tmp_x(ixa,ixbp+1) = v_tmp_x(ixa,ixbp+1) + fact * ixbp * v_tmp_x(ixa,ixbp-1)
         if(ixa>0)  v_tmp_x(ixa,ixbp+1) = v_tmp_x(ixa,ixbp+1) + fact * ixa  * v_tmp_x(ixa-1,ixbp)
       endif
  
     enddo
  
   enddo
  
   !
   ! direction Y
   !
   v_tmp_y(0,0) = v_tmp_x(ga%nx,gb%nx)
  
   !
   ! direction Z
   !
   v_tmp_z(0,0) = v_tmp_y(ga%ny,gb%ny)






   v_ab = v_ab + wu(isamp) * v_tmp_z(ga%nz,gb%nz)

 enddo !isamp
 v_ab = - zatom * v_ab * 2.0_dp / SQRT(pi) * ga%norm_factor * gb%norm_factor


end subroutine nucleus_recurrence


#endif


!=========================================================================
subroutine numerical_overlap(ga,gb)
 implicit none
 type(gaussian),intent(in) :: ga,gb
!=====
 integer,parameter  :: nx=100
 real(dp),parameter :: rmax=10.
 real(dp)           :: dx,rtmp,x(3)
 integer            :: ix,iy,iz
!=====

 dx = rmax/REAL(nx,dp)

 rtmp=0.0d0
 do ix=1,nx
   x(1) = ( REAL(ix,dp)/DBLE(nx) - 0.5 ) * rmax
   do iy=1,nx
     x(2) = ( REAL(iy,dp)/DBLE(nx) - 0.5 ) * rmax
     do iz=1,nx
       x(3) = ( REAL(iz,dp)/DBLE(nx) - 0.5 ) * rmax
  
       rtmp = rtmp + eval_gaussian(ga,x) * eval_gaussian(gb,x) * dx**3
  
     enddo
   enddo
 enddo

 write(*,*) 'check S_ab',rtmp

end subroutine numerical_overlap

!=========================================================================
subroutine numerical_kinetic(ga,gb)
 implicit none
 type(gaussian),intent(in) :: ga,gb
!=====
 integer,parameter  :: nx=100
 real(dp),parameter :: rmax=10.0_dp
 real(dp),parameter :: dh=1.d-4
 real(dp)           :: dx,rtmp,x(3),dhx(3),dhy(3),dhz(3)
 integer            :: ix,iy,iz
!=====

 dx = rmax/REAL(nx,dp)
 dhx(:) = 0.0_dp
 dhx(1) = dh
 dhy(:) = 0.0_dp
 dhy(2) = dh
 dhz(:) = 0.0_dp
 dhz(3) = dh

 rtmp=0.0d0
 do ix=1,nx
   x(1) = ( REAL(ix,dp)/DBLE(nx) - 0.5 ) * rmax
   do iy=1,nx
     x(2) = ( REAL(iy,dp)/DBLE(nx) - 0.5 ) * rmax
     do iz=1,nx
       x(3) = ( REAL(iz,dp)/DBLE(nx) - 0.5 ) * rmax
  
       rtmp = rtmp - 0.5 * eval_gaussian(ga,x) * dx**3 &
                    * ( eval_gaussian(gb,x+dhx) + eval_gaussian(gb,x-dhx) &
                       +eval_gaussian(gb,x+dhy) + eval_gaussian(gb,x-dhy) &
                       +eval_gaussian(gb,x+dhz) + eval_gaussian(gb,x-dhz) &
                       - 6.0 * eval_gaussian(gb,x) ) / dh**2 
  
     enddo
   enddo
 enddo

 write(*,*) 'check K_ab',rtmp

end subroutine numerical_kinetic

!=========================================================================
subroutine numerical_nucleus(ga,gb)
 implicit none
 type(gaussian),intent(in) :: ga,gb
!=====
 integer,parameter  :: nx=100
 real(dp),parameter :: rmax=10.
 real(dp)           :: dx,rtmp,x(3)
 integer            :: ix,iy,iz
!=====

 dx = rmax/REAL(nx,dp)

 write(*,*) 'hydrogen atom in zero'

 rtmp=0.0d0
 do ix=1,nx
   x(1) = ( REAL(ix,dp)/DBLE(nx) - 0.5 ) * rmax
   do iy=1,nx
     x(2) = ( REAL(iy,dp)/DBLE(nx) - 0.5 ) * rmax
     do iz=1,nx
       x(3) = ( REAL(iz,dp)/DBLE(nx) - 0.5 ) * rmax

       if( SUM(x(:)**2) < 1.d-8 ) cycle ! skip the singularity

       rtmp = rtmp + eval_gaussian(ga,x) * eval_gaussian(gb,x) * dx**3 * -1.0 / SQRT(SUM(x(:)**2))
  
     enddo
   enddo
 enddo

 write(*,*) 'check V_ab',rtmp

end subroutine numerical_nucleus

!=========================================================================
end module m_gaussian

