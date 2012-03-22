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
subroutine overlap_recurrence(ga,gb,s_ab)
 implicit none
 type(gaussian),intent(in) :: ga,gb
 real(dp),intent(out) :: s_ab
!=====
 real(dp)             :: alpha_ab,ksi_ab,d2_ab
 real(dp)             :: xp(3)
 real(dp)             :: s_tmp(0:ga%nx,0:ga%ny,0:ga%nz,0:gb%nx,0:gb%ny,0:gb%nz)
!=====

 alpha_ab = ga%alpha + gb%alpha
 ksi_ab   = ga%alpha * gb%alpha / alpha_ab
 d2_ab    = SUM( ( ga%x0(:)-gb%x0(:) )**2 )
 xp(:)    = ( ga%alpha * ga%x0(:) + gb%alpha * gb%x0(:) ) / alpha_ab

 s_tmp(0,0,0,0,0,0) = (pi/alpha_ab)**1.5_dp * EXP( - ksi_ab * d2_ab )

 s_tmp(1,0,0,0,0,0) = ( xp(1) - ga%x0(1) ) * s_tmp(0,0,0,0,0,0)

 s_tmp(2,0,0,0,0,0) = ( xp(1) - ga%x0(1) ) * s_tmp(1,0,0,0,0,0) + 1.0 / ( 2.0 * alpha_ab) * s_tmp(0,0,0,0,0,0) 

 write(*,*) 's_tmp',s_tmp(0,0,0,0,0,0)
 write(*,*) 's_tmp',s_tmp(1,0,0,0,0,0)
 write(*,*) 's_tmp',s_tmp(2,0,0,0,0,0)

 stop'ENOUGH'
end subroutine overlap_recurrence


!=========================================================================
subroutine shift_gaussian(ga,x0,ngshifted,gshifted)
 implicit none
 type(gaussian),intent(in)            :: ga
 real(dp),intent(in)                  :: x0(3)
 integer,intent(out)                  :: ngshifted
 type(gaussian),pointer,intent(inout) :: gshifted(:)
!=====
 integer                              :: ngx,ngy,ngz
 integer                              :: igx,igy,igz
 integer                              :: igshifted
 integer,allocatable                  :: nx(:),ny(:),nz(:)
 real(dp),allocatable                 :: factx(:),facty(:),factz(:)
!=====

 !
 ! Cartesian component: X
 select case(ga%nx)
 case(0)
   ngx=1
   allocate(nx(ngx),factx(ngx))
   nx(1)    = 0
   factx(1) = 1.0_dp
 case(1)
   ngx=2
   allocate(nx(ngx),factx(ngx))
   nx(1)    = 1
   factx(1) = 1.0_dp
   nx(2)    = 0
   factx(2) = x0(1)-ga%x0(1)
 case(2)
   ngx=3
   allocate(nx(ngx),factx(ngx))
   nx(1)    = 2
   factx(1) = 1.0_dp
   nx(2)    = 1
   factx(2) = 2.0_dp * ( x0(1)-ga%x0(1) )
   nx(3)    = 0
   factx(3) = ( x0(1)-ga%x0(1) )**2
 case default
  stop'shift_gaussian: NOT IMPLEMENTED'
 end select
 !
 ! Cartesian component: Y
 select case(ga%ny)
 case(0)
   ngy=1
   allocate(ny(ngy),facty(ngy))
   ny(1)    = 0
   facty(1) = 1.0_dp
 case(1)
   ngy=2
   allocate(ny(ngy),facty(ngy))
   ny(1)    = 1
   facty(1) = 1.0_dp
   ny(2)    = 0
   facty(2) = x0(2)-ga%x0(2)
 case(2)
   ngy=3
   allocate(ny(ngy),facty(ngy))
   ny(1)    = 2
   facty(1) = 1.0_dp
   ny(2)    = 1
   facty(2) = 2.0_dp * ( x0(2)-ga%x0(2) )
   ny(3)    = 0
   facty(3) = ( x0(2)-ga%x0(2) )**2
 case default
  stop'shift_gaussian: NOT IMPLEMENTED'
 end select
 !
 ! Cartesian component: Z
 select case(ga%nz)
 case(0)
   ngz=1
   allocate(nz(ngz),factz(ngz))
   nz(1)    = 0
   factz(1) = 1.0_dp
 case(1)
   ngz=2
   allocate(nz(ngz),factz(ngz))
   nz(1)    = 1
   factz(1) = 1.0_dp
   nz(2)    = 0
   factz(2) = x0(3)-ga%x0(3)
 case(2)
   ngz=3
   allocate(nz(ngz),factz(ngz))
   nz(1)    = 2
   factz(1) = 1.0_dp
   nz(2)    = 1
   factz(2) = 2.0_dp * ( x0(3)-ga%x0(3) )
   nz(3)    = 0
   factz(3) = ( x0(3)-ga%x0(3) )**2
 case default
  stop'shift_gaussian: NOT IMPLEMENTED'
 end select

 ngshifted=ngx*ngy*ngz
 allocate(gshifted(ngshifted))
 igshifted=0
 do igx=1,ngx
   do igy=1,ngy
     do igz=1,ngz
       igshifted=igshifted+1
       call init_gaussian_general(nx(igx),ny(igy),nz(igz),ga%alpha,x0,gshifted(igshifted))
       gshifted%norm_factor = gshifted%norm_factor * factx(igx) * facty(igy) * factz(igz)
     enddo
   enddo
 enddo

 deallocate(nx,ny,nz)
 deallocate(factx,facty,factz)

end subroutine shift_gaussian

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
subroutine product_gaussian_general(ga,gb,gprod)
 implicit none
 type(gaussian),intent(in)  :: ga,gb
 type(gaussian),pointer,intent(out) :: gprod(:)
!=====
 integer                    :: ngprod,igprod
 integer                    :: iga,igb
 integer                    :: nx_c,ny_c,nz_c
 real(dp)                   :: alpha_c,xc(3),fact_c
 integer                    :: ngshifted_a,ngshifted_b
 type(gaussian),pointer     :: gshifted_a(:)
 type(gaussian),pointer     :: gshifted_b(:)
!=====

 if(ga%am>2 .OR. gb%am>2) stop'gaussian product with different centers not yet implemented for l>2'

 alpha_c = ga%alpha + gb%alpha
 xc(:)   = ( ga%alpha * ga%x0(:) + gb%alpha * gb%x0(:) ) / alpha_c
 fact_c  = EXP( alpha_c * SUM(xc(:)**2) - ga%alpha * SUM(ga%x0(:)**2) - gb%alpha * SUM(gb%x0(:)**2) )

 if( ALL( ABS(ga%x0(:) - gb%x0(:)) < 1.d-6 ) ) then
   ngprod = 1
   allocate(gprod(ngprod))
   call init_gaussian(ga%nx+gb%nx,ga%ny+gb%ny,ga%nz+gb%nz,alpha_c,gprod(1))
   gprod(igprod)%norm_factor = gprod(igprod)%norm_factor * fact_c
 else
   call shift_gaussian(ga,xc,ngshifted_a,gshifted_a)
   call shift_gaussian(gb,xc,ngshifted_b,gshifted_b)
   ngprod = ngshifted_a * ngshifted_b
   allocate(gprod(ngprod))
   igprod=0
   do iga=1,ngshifted_a
     do igb=1,ngshifted_b
       igprod=igprod+1
       call product_gaussian(gshifted_a(iga),gshifted_b(igb),gprod(igprod))
       write(*,*) iga,igb
       gprod(igprod)%norm_factor = gprod(igprod)%norm_factor * fact_c
     enddo
   enddo
   deallocate(gshifted_a)
   deallocate(gshifted_b)
 endif


end subroutine product_gaussian_general


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
subroutine overlap_normalized_general(ga,gb,s_ab)
 implicit none
 type(gaussian),intent(in) :: ga,gb
 real(dp),intent(out) :: s_ab
!=====
 type(gaussian),pointer :: gprod(:)
 integer                :: igprod
!=====

 call product_gaussian_general(ga,gb,gprod)
 write(*,*) 'SIZE',SIZE(gprod(:))
 s_ab=0.0_dp
 do igprod=1,SIZE(gprod(:))
   s_ab = s_ab &
         +integral_1d(gprod(igprod)%alpha,gprod(igprod)%nx) &
         *integral_1d(gprod(igprod)%alpha,gprod(igprod)%ny) &
         *integral_1d(gprod(igprod)%alpha,gprod(igprod)%nz) * gprod(igprod)%norm_factor
 enddo
 deallocate(gprod)

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
  
  end function
end subroutine overlap_normalized_general

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
  
  end function
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

!=========================================================================
end module m_gaussian

