!=========================================================================
module m_dft
 use m_definitions

contains

subroutine dft_exc_vxc(nspin,basis,dft_xc,p_matrix,vxc_ij,exc_xc)
 use m_tools,only: coeffs_gausslegint
 use m_basis_set
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
#ifdef OPENMP
 use omp_lib
#endif
 implicit none

 integer,intent(in)         :: dft_xc(2)
 integer,intent(in)         :: nspin
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: p_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       :: vxc_ij(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       ::exc_xc
!=====
 integer,parameter :: nx=40
 integer,parameter :: nangular=38 ! 86
 integer,parameter :: nr=nx*nangular
 real(dp)          :: weight
 real(dp)          :: x1(nangular)
 real(dp)          :: y1(nangular)
 real(dp)          :: z1(nangular)
 real(dp)          :: w1(nangular)
 integer           :: n1,iangular

 integer :: ix,iy,iz,ibf,jbf,ispin,ir
 real(dp) :: rhor_r(nspin),rhor(nr,nspin),grad_rhor(3,nspin),sigma2(2*nspin-1),rr(3)
 real(dp) :: x(nx),wx(nx)
 real(dp) :: normalization(nspin)
! integer,parameter :: ixc=1  !Slater exchange
 integer,parameter :: ixc=2  !Teter exchange-correlation

#ifdef HAVE_LIBXC
 type(xc_f90_pointer_t) :: xc_func1,xc_func2
 type(xc_f90_pointer_t) :: xc_info1,xc_info2
#endif

 real(dp) :: vxc1(nspin),vxc2(nspin),exc1,exc2,vsigma1(2*nspin-1),vsigma2(2*nspin-1)
 real(dp) :: vxc_av(nspin)
 real(dp) :: dedd(nr,nspin)
 real(dp) :: dedgd(nr,3,nspin)
 real(dp) :: tmpx(nr,nspin),tmpy(nr,nspin),tmpz(nr,nspin)
 real(dp),parameter :: dx=1.0d-5

!=====

#ifdef HAVE_LIBXC

 write(*,*) 'Evaluate DFT integrals'
 write(*,'(a,i4,x,i4)') ' discretization grid [radial points , angular points] ',nx,nangular

 if(nspin==1) then
   call xc_f90_func_init(xc_func1, xc_info1, dft_xc(1), XC_UNPOLARIZED)
   call xc_f90_func_init(xc_func2, xc_info2, dft_xc(2), XC_UNPOLARIZED)
 else
   call xc_f90_func_init(xc_func1, xc_info1, dft_xc(1), XC_POLARIZED)
   call xc_f90_func_init(xc_func2, xc_info2, dft_xc(2), XC_POLARIZED)
 endif

 !
 ! spherical integration
 ! radial part with Gauss-Legendre
 call coeffs_gausslegint(-1.0_dp,1.0_dp,x,wx,nx)
 !
 ! Transformation from [-1;1] to [0;+\infty[
 ! taken from M. Krack JCP 1998
 wx(:) = wx(:) * ( 1.0_dp / log(2.0_dp) / ( 1.0_dp - x(:) ) )
 x(:)  = log( 2.0_dp / (1.0_dp - x(:) ) ) / log(2.0_dp)

 ! angular part with Lebedev - Laikov
 ! (x1,y1,z1) on the unit sphere
 n1=nangular
 select case(nangular)
 case(6)
   call ld0006(x1,y1,z1,w1,n1)
 case(14)
   call ld0014(x1,y1,z1,w1,n1)
 case(26)
   call ld0026(x1,y1,z1,w1,n1)
 case(38)
   call ld0038(x1,y1,z1,w1,n1)
 case(50)
   call ld0050(x1,y1,z1,w1,n1)
 case(74)
   call ld0074(x1,y1,z1,w1,n1)
 case(86)
   call ld0086(x1,y1,z1,w1,n1)
 case default
   write(*,*) 'grid points: ',nangular
   stop'Lebedev grid is not available'
 end select
 

 exc_xc=0.0_dp
 normalization(:)=0.0_dp

!$OMP PARALLEL DEFAULT(SHARED)

! write(*,*) 'get_num_thread()/max:',OMP_get_num_threads(),OMP_get_max_threads()

!$OMP DO REDUCTION(+:exc_xc,normalization) PRIVATE(rr,weight,ir,rhor_r,grad_rhor,exc1,exc2,vxc1,vxc2)
 do ix=1,nx
   do iangular=1,nangular
     rr(1) = x(ix) * x1(iangular)
     rr(2) = x(ix) * y1(iangular)
     rr(3) = x(ix) * z1(iangular)
     weight = wx(ix) * w1(iangular) * x(ix)**2 * 4.0_dp * pi
     ir=iangular+(ix-1)*nangular

     rhor_r(:)=0.0_dp
     do ispin=1,nspin
       do jbf=1,basis%nbf
         do ibf=1,basis%nbf
           rhor_r(ispin)=rhor_r(ispin)+p_matrix(ibf,jbf,ispin)&
                             * eval_basis_function(basis%bf(ibf),rr) * eval_basis_function(basis%bf(jbf),rr)
         enddo
       enddo
     enddo

     normalization(:) = normalization(:) + rhor_r(:) * weight

     if(xc_f90_info_family(xc_info1) == XC_FAMILY_GGA) then
       grad_rhor(:,:)=0.0_dp
       do ispin=1,nspin
         do jbf=1,basis%nbf
           do ibf=1,basis%nbf
             grad_rhor(:,ispin)=grad_rhor(:,ispin)+p_matrix(ibf,jbf,ispin)&
&                 *( eval_basis_function_derivative(basis%bf(ibf),rr) * eval_basis_function(basis%bf(jbf),rr) &
&                  + eval_basis_function_derivative(basis%bf(jbf),rr) * eval_basis_function(basis%bf(ibf),rr) ) 

           enddo
         enddo
       enddo

       sigma2(1) = SUM( grad_rhor(:,1)**2 )
       if(nspin==2) then
         sigma2(2) = SUM( grad_rhor(:,1)*grad_rhor(:,2) )
         sigma2(3) = SUM( grad_rhor(:,2)**2 )
       endif
     endif

!LIBXC called
     if(xc_f90_info_family(xc_info1) == XC_FAMILY_LDA) then 
       if(dft_xc(1)/=0) then
         call xc_f90_lda_exc_vxc(xc_func1,1,rhor_r(1),exc1,vxc1(1))
       else
         exc1=0.0_dp
         vxc1=0.0_dp
       endif
       if(dft_xc(2)/=0) then
         call xc_f90_lda_exc_vxc(xc_func2,1,rhor_r(1),exc2,vxc2(1))
       else
         exc2=0.0_dp
         vxc2=0.0_dp
       endif
     else if(xc_f90_info_family(xc_info1) == XC_FAMILY_GGA) then
       if(dft_xc(1)/=0) then
         call xc_f90_gga_exc_vxc(xc_func1,1,rhor_r(1),sigma2(1),exc1,vxc1(1),vsigma1(1))
       else
         exc1=0.0_dp
         vxc1=0.0_dp
         vsigma1=0.0_dp
       endif
       if(dft_xc(2)/=0) then
         call xc_f90_gga_exc_vxc(xc_func2,1,rhor_r(1),sigma2(1),exc2,vxc2(1),vsigma2(1))
       else
         exc2=0.0_dp
         vxc2=0.0_dp
         vsigma2=0.0_dp
       endif
     endif

!     rhor(ir,:) = rhor_r(:)
     exc_xc = exc_xc + weight * (exc1+exc2) * SUM( rhor_r(:) )

     dedd(ir,:) = vxc1(:) + vxc2(:)

     if(xc_f90_info_family(xc_info1) == XC_FAMILY_GGA) then
       if(nspin==1) then
         dedgd(ir,:,1) = 2.0_dp * ( vsigma1(1) + vsigma2(1) ) * grad_rhor(:,1) 
       else
         dedgd(ir,:,1) = 2.0_dp * ( vsigma1(1) + vsigma2(1) ) * grad_rhor(:,1) + ( vsigma1(2) + vsigma2(2) ) * grad_rhor(:,2)
         dedgd(ir,:,2) = 2.0_dp * ( vsigma1(3) + vsigma2(3) ) * grad_rhor(:,2) + ( vsigma1(2) + vsigma2(2) ) * grad_rhor(:,1)
       endif

       call xc_f90_gga_exc_vxc(xc_func1,1,rhor_r(1),sigma2(1),exc1,vxc1(1),vsigma1(1))
       if(dft_xc(2)/=0) then
         call xc_f90_gga_exc_vxc(xc_func2,1,rhor_r(1),sigma2(1),exc2,vxc2(1),vsigma2(1))
       else
         exc2=0.0_dp
         vxc2=0.0_dp
         vsigma2=0.0_dp
       endif

       stop'to be implemented in spherical coord'
!       tmpx(:) = ( dedgd(ir+nz*ny,1,:) - dedgd(ir-nz*ny,1,:) ) / ( x(ix+1)-x(ix-1) )
!       tmpy(:) = ( dedgd(ir+nz   ,2,:) - dedgd(ir-nz   ,2,:) ) / ( y(iy+1)-y(iy-1) )
!       tmpz(:) = ( dedgd(ir+1    ,3,:) - dedgd(ir-1    ,3,:) ) / ( z(iz+1)-z(iz-1) )

     else
       tmpx(ir,:) = 0.0_dp
       tmpy(ir,:) = 0.0_dp
       tmpz(ir,:) = 0.0_dp
     endif


   enddo
 enddo ! loop on the radial grid
!$OMP END DO

!$OMP END PARALLEL

 ir=0
 vxc_ij(:,:,:)=0.0_dp
 do ix=1,nx
   do iangular=1,nangular
     rr(1) = x(ix) * x1(iangular)
     rr(2) = x(ix) * y1(iangular)
     rr(3) = x(ix) * z1(iangular)
     weight = wx(ix) * w1(iangular) * x(ix)**2 * 4.0_dp * pi
     ir=ir+1

     do jbf=1,basis%nbf
       do ibf=1,basis%nbf
         vxc_ij(ibf,jbf,:) =  vxc_ij(ibf,jbf,:) + weight &
             * eval_basis_function(basis%bf(ibf),rr) * eval_basis_function(basis%bf(jbf),rr) &
             * ( dedd(ir,:) - tmpx(ir,:) - tmpy(ir,:) - tmpz(ir,:) )
       enddo
     enddo

   enddo
 enddo ! loop on the radial grid




#else
 exc_xc=0.0_dp
 vxc_ij(:,:,:)=0.0_dp
#endif

! write(*,*)
! write(*,'(a,2(2x,f12.6))') 'Average Vxc',vxc_av
 write(*,*)
 write(*,'(a,2(2x,f12.6))') 'number of electrons',normalization(:)
 write(*,'(a,2x,f12.6)') 'DFT xc energy [Ha]:',exc_xc
 write(*,*)

contains
 function rho_hydrogen(x,y,z)
 implicit none
 real(dp) :: rho_hydrogen
 real(dp) :: x,y,z

 rho_hydrogen =  exp( -2.0 * sqrt( (x-5.)**2 + (y-5.)**2 + (z-5.)**2) ) / (4.0_dp*pi) * 4.0_dp *0.9999 

 end function
 function vlda(ixc,rhor)
 implicit none
 real(dp) :: vlda
 real(dp),intent(in) :: rhor
 integer,intent(in) :: ixc

! real(dp),parameter :: a0= 0.00033959499_dp
! real(dp),parameter :: a1= 0.1912460_dp
! real(dp),parameter :: a2= 0.8008790_dp
! real(dp),parameter :: a3= 0.092956297_dp
! real(dp),parameter :: b1= 1.0_dp
! real(dp),parameter :: b2= 8.710930_dp
! real(dp),parameter :: b3= 3.945600_dp
! real(dp),parameter :: b4= 0.087989897_dp

 real(dp),parameter :: a0=.4581652932831429_dp,a1=2.40875407_dp,a2=.88642404_dp
 real(dp),parameter :: a3=.02600342_dp,b1=1.0_dp,b2=4.91962865_dp
 real(dp),parameter :: b3=1.34799453_dp,b4=.03120453_dp

 real(dp),parameter :: c1=4._dp*a0*b1/3.0_dp
 real(dp),parameter :: c2=5.0_dp*a0*b2/3.0_dp+a1*b1
 real(dp),parameter :: c3=2.0_dp*a0*b3+4.0_dp*a1*b2/3.0_dp+2.0_dp*a2*b1/3.0_dp
 real(dp),parameter :: c4=7.0_dp*a0*b4/3.0_dp+5.0_dp*a1*b3/3.0_dp+a2*b2+a3*b1/3.0_dp
 real(dp),parameter :: c5=2.0_dp*a1*b4+4.0_dp*a2*b3/3.0_dp+2.0_dp*a3*b2/3.0_dp
 real(dp),parameter :: c6=5.0_dp*a2*b4/3.0_dp+a3*b3,c7=4.0_dp*a3*b4/3.0_dp

 real(dp),parameter :: vfac=(1.5_dp/pi)**(2.0_dp/3.0_dp)

 real(dp) :: rs,num,den

 rs = ( 3.0_dp / ( 4.0_dp * pi * rhor ) )**(1.0_dp/3.0_dp)

 if(ixc==1) then

     vlda = -vfac / rs

 else

   num = -( c1*rs + c2*rs**2 + c3*rs**3 + c4*rs**4 + c5*rs**5 + c6*rs**6 + c7 *rs**7 )
   den =    b1*rs + b2*rs**2 + b3*rs**3 + b4*rs**4 
  
   vlda = num/den**2

 endif

 end function vlda
 function elda(ixc,rhor)
 implicit none
 real(dp) :: elda
 real(dp),intent(in) :: rhor
 integer,intent(in) :: ixc

! real(dp),parameter :: a0= 0.00033959499_dp
! real(dp),parameter :: a1= 0.1912460_dp
! real(dp),parameter :: a2= 0.8008790_dp
! real(dp),parameter :: a3= 0.092956297_dp
! real(dp),parameter :: b1= 1.0_dp
! real(dp),parameter :: b2= 8.710930_dp
! real(dp),parameter :: b3= 3.945600_dp
! real(dp),parameter :: b4= 0.087989897_dp

 real(dp),parameter :: a0=.4581652932831429_dp,a1=2.40875407_dp,a2=.88642404_dp
 real(dp),parameter :: a3=.02600342_dp,b1=1.0_dp,b2=4.91962865_dp
 real(dp),parameter :: b3=1.34799453_dp,b4=.03120453_dp

 real(dp),parameter :: c1=4._dp*a0*b1/3.0_dp
 real(dp),parameter :: c2=5.0_dp*a0*b2/3.0_dp+a1*b1
 real(dp),parameter :: c3=2.0_dp*a0*b3+4.0_dp*a1*b2/3.0_dp+2.0_dp*a2*b1/3.0_dp
 real(dp),parameter :: c4=7.0_dp*a0*b4/3.0_dp+5.0_dp*a1*b3/3.0_dp+a2*b2+a3*b1/3.0_dp
 real(dp),parameter :: c5=2.0_dp*a1*b4+4.0_dp*a2*b3/3.0_dp+2.0_dp*a3*b2/3.0_dp
 real(dp),parameter :: c6=5.0_dp*a2*b4/3.0_dp+a3*b3,c7=4.0_dp*a3*b4/3.0_dp

 real(dp),parameter :: efac=0.75_dp*(1.5_dp/pi)**(2.0_dp/3.0_dp)

 real(dp) :: rs,num,den

 rs = ( 3.0_dp / ( 4.0_dp * pi * rhor ) )**(1.0_dp/3.0_dp)

 if(ixc==1) then

   elda = -efac / rs

 else

   num = -( a0 + a1*rs + a2*rs**2 + a3*rs**3 )
   den =    b1*rs + b2*rs**2 + b3*rs**3 + b4*rs**4
  
   elda = num/den

 endif

 end function elda

end subroutine dft_exc_vxc

end module m_dft
