!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! some basic numerical algorithms
!
!=========================================================================
module m_numerical_tools
  use m_definitions
  use m_warning,only: die


contains


!=========================================================================
!
! Compute the coefficients (supports and weights)
! for Gauss-Legendre integration.
! Inspired by a routine due to G. Rybicki.
!
! Input:
! xmin=lower bound of integration
! xmax=upper bound of integration
! n=order of integration
!
! Output:
! x(n)=array of support points
! weights(n)=array of integration weights
!
subroutine coeffs_gausslegint(xmin,xmax,x,weights,n)

  implicit none
 
  integer,intent(in) :: n
  real(dp),intent(in) :: xmin,xmax
  real(dp),intent(out) :: x(n),weights(n)
  !=====
  real(dp) :: tol,xl,xmean,z,p1,p2,p3,pp,z1
  integer  :: i,j
  !=====
 
  tol=1.d-13
 
  xl=(xmax-xmin)*0.50_dp
  xmean=(xmax+xmin)*0.50_dp
 
  do i=1,(n+1)/2
   z = COS(pi*(i-0.250_dp)/(n+0.50_dp))
 
   do
 
     p1=1.0_dp
     p2=0.0_dp
 
     do j=1,n
 
      p3=p2
      p2=p1
      p1=((2.0_dp*j - 1.0_dp)*z*p2 - (j-1.0_dp)*p3)/j
 
     enddo
 
     pp=n*(p2-z*p1)/(1.0_dp - z**2)
     z1=z
     z=z1-p1/pp
 
     if(abs(z-z1) < tol) exit
 
   enddo
 
   x(i)=xmean-xl*z
   x(n+1-i)=xmean+xl*z
   weights(i)=2.0_dp * xl/((1.0_dp-z**2)*pp**2)
   weights(n+1-i)=weights(i)
 
  enddo

end subroutine coeffs_gausslegint


!=========================================================================
function pade(z_in,f_in,z_out) RESULT(f_out)
  implicit none
 
  complex(dp),intent(in) :: z_in(:),f_in(:)
  complex(dp),intent(in) :: z_out
  complex(dp) :: f_out
  !=====
  integer     :: ip,np
  complex(dp),allocatable :: aa(:)
  complex(dp),allocatable :: Az(:), Bz(:)
  !=====
 
  np = SIZE(z_in)
  allocate(aa(np),Az(0:np),Bz(0:np))
 
  call calculate_pade_a(aa,z_in,f_in)
 
  Az(0) = (0.0_dp,0.0_dp)
  Az(1) = aa(1)
  Bz(0) = (1.0_dp,0.0_dp)
  Bz(1) = (1.0_dp,0.0_dp)
 
  do ip=1,np-1
    Az(ip+1) = Az(ip) + ( z_out - z_in(ip) ) * aa(ip+1) * Az(ip-1)
    Bz(ip+1) = Bz(ip) + ( z_out - z_in(ip) ) * aa(ip+1) * Bz(ip-1)
  enddo
 
  f_out = Az(np) / Bz(np)
 
  deallocate(aa,Az,Bz)

end function pade


!=========================================================================
subroutine calculate_pade_a(aa,z_in,f_in)
  implicit none
 
  complex(dp),intent(in)  :: z_in(:),f_in(:)
  complex(dp),intent(out) :: aa(:)
  !=====
  integer     :: np,ip,jp
  complex(dp),allocatable :: gtmp(:,:)
  !=====
 
  np = SIZE(z_in)
  allocate(gtmp(np,np))
 
  gtmp(1,1:np) = f_in(1:np)
 
  do ip=2,np
    do jp=ip,np
      gtmp(ip,jp) = (gtmp(ip-1,ip-1) - gtmp(ip-1,jp)) / ( (z_in(jp) - z_in(ip-1)) * gtmp(ip-1,jp) )
    enddo
  enddo
  do ip=1,np
    aa(ip) = gtmp(ip,ip)
  enddo
 
  deallocate(gtmp)

end subroutine calculate_pade_a


!=========================================================================
end module m_numerical_tools
!=========================================================================
