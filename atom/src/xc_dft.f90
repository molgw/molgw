!=========================================================================
subroutine dft_exc_vxc(nspin,basis,dft_xc,p_matrix,vxc_ij,exc_xc)
 use m_tools,only: coeffs_gausslegint
 use m_timing
 use m_atoms
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
 real(dp),intent(out)       :: exc_xc
!=====
 integer,parameter :: nx=40
 integer,parameter :: nangular=38 ! 86
! integer,parameter :: nr=nx*nangular
 real(dp)          :: weight
 real(dp)          :: x1(nangular)
 real(dp)          :: y1(nangular)
 real(dp)          :: z1(nangular)
 real(dp)          :: w1(nangular)
 integer           :: n1,iangular

 integer :: ix,iy,iz,ibf,jbf,ispin,iatom,jatom,katom
 real(dp) :: rhor_r(nspin),grad_rhor(3,nspin),sigma2(2*nspin-1),rr(3)
 real(dp) :: xa(nx),wxa(nx)
 real(dp) :: normalization(nspin)
! integer,parameter :: ixc=1  !Slater exchange
 integer,parameter :: ixc=2  !Teter exchange-correlation

#ifdef HAVE_LIBXC
 type(xc_f90_pointer_t) :: xc_func1,xc_func2
 type(xc_f90_pointer_t) :: xc_info1,xc_info2
#endif

 real(dp) :: vxc1(nspin),vxc2(nspin),exc1,exc2,vsigma1(2*nspin-1),vsigma2(2*nspin-1)
 real(dp) :: vxc_av(nspin)
! real(dp) :: dedd(nr,nspin)
! real(dp) :: dedgd(nr,3,nspin)
! real(dp) :: tmpx(nr,nspin),tmpy(nr,nspin),tmpz(nr,nspin)
 real(dp) :: dedd_r(nspin)
 real(dp) :: dedgd_r(3,nspin)
 real(dp) :: tmpx_r(nspin),tmpy_r(nspin),tmpz_r(nspin)
 real(dp),parameter :: dx=1.0d-5

 real(dp) :: mu,s_becke(natom,natom),p_becke(natom),fact_becke
!=====

#ifdef HAVE_LIBXC

 write(*,*) 'Evaluate DFT integrals'
 write(*,'(a,i4,x,i4)') '   discretization grid per atom [radial points , angular points] ',nx,nangular
 write(*,'(a,i8)') '   total number of real-space points',nx*nangular*natom

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
 call coeffs_gausslegint(-1.0_dp,1.0_dp,xa,wxa,nx)
 !
 ! Transformation from [-1;1] to [0;+\infty[
 ! taken from M. Krack JCP 1998
 wxa(:) = wxa(:) * ( 1.0_dp / log(2.0_dp) / ( 1.0_dp - xa(:) ) )
 xa(:)  = log( 2.0_dp / (1.0_dp - xa(:) ) ) / log(2.0_dp)

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
 vxc_ij(:,:,:)=0.0_dp
 normalization(:)=0.0_dp


 do iatom=1,natom

!!!! !$OMP PARALLEL DEFAULT(SHARED)
!!!! !$OMP DO REDUCTION(+:exc_xc,normalization,vxc_ij) PRIVATE(rr,weight,mu,s_becke,p_becke,fact_becke,rhor_r,grad_rhor,exc1,exc2,vxc1,vxc2,dedd_r,dedgd_r) COLLAPSE(2)
   do ix=1,nx
     do iangular=1,nangular
       rr(1) = xa(ix) * x1(iangular) + x(1,iatom) 
       rr(2) = xa(ix) * y1(iangular) + x(2,iatom)
       rr(3) = xa(ix) * z1(iangular) + x(3,iatom)
       weight = wxa(ix) * w1(iangular) * xa(ix)**2 * 4.0_dp * pi
  
       !
       ! Partitionning scheme of Axel Becke, J. Chem. Phys. 88, 2547 (1988).
       !
       s_becke(:,:) = 0.0_dp
       do katom=1,natom
         do jatom=1,natom
           if(katom==jatom) cycle
           mu = ( SQRT( SUM( (rr(:)-x(:,katom) )**2 ) ) - SQRT( SUM( (rr(:)-x(:,jatom) )**2 ) ) ) &
                     / SQRT( SUM( (x(:,katom)-x(:,jatom))**2 ) )
           s_becke(katom,jatom) = 0.5_dp * ( 1.0_dp - smooth_step(smooth_step(smooth_step(mu))) )
         enddo
       enddo
       p_becke(:) = 1.0_dp
       do katom=1,natom
         do jatom=1,natom
           if(katom==jatom) cycle
           p_becke(katom) = p_becke(katom) * s_becke(katom,jatom)
         enddo
       enddo
       fact_becke = p_becke(iatom) / SUM( p_becke(:) )
  
  
 call start_clock(timing_tmp5)
       rhor_r(:)=0.0_dp
       do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO REDUCTION(+:rhor_r) COLLAPSE(2)
         do jbf=1,basis%nbf
           do ibf=1,basis%nbf
             rhor_r(ispin)=rhor_r(ispin)+p_matrix(ibf,jbf,ispin)&
                               * eval_basis_function(basis%bf(ibf),rr) * eval_basis_function(basis%bf(jbf),rr)
           enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
       enddo
 call stop_clock(timing_tmp5)
  
       normalization(:) = normalization(:) + rhor_r(:) * weight * fact_becke
  
       if(xc_f90_info_family(xc_info1) == XC_FAMILY_GGA) then
         grad_rhor(:,:)=0.0_dp
         do ispin=1,nspin
           do jbf=1,basis%nbf
             do ibf=1,basis%nbf
               grad_rhor(:,ispin)=grad_rhor(:,ispin)+p_matrix(ibf,jbf,ispin)&
                    *( eval_basis_function_derivative(basis%bf(ibf),rr) * eval_basis_function(basis%bf(jbf),rr) &
                     + eval_basis_function_derivative(basis%bf(jbf),rr) * eval_basis_function(basis%bf(ibf),rr) ) 
  
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
  
       exc_xc = exc_xc + weight * fact_becke * (exc1+exc2) * SUM( rhor_r(:) )
  
       dedd_r(:) = vxc1(:) + vxc2(:)
  
       if(xc_f90_info_family(xc_info1) == XC_FAMILY_GGA) then
         if(nspin==1) then
           dedgd_r(:,1) = 2.0_dp * ( vsigma1(1) + vsigma2(1) ) * grad_rhor(:,1) 
         else
           dedgd_r(:,1) = 2.0_dp * ( vsigma1(1) + vsigma2(1) ) * grad_rhor(:,1) + ( vsigma1(2) + vsigma2(2) ) * grad_rhor(:,2)
           dedgd_r(:,2) = 2.0_dp * ( vsigma1(3) + vsigma2(3) ) * grad_rhor(:,2) + ( vsigma1(2) + vsigma2(2) ) * grad_rhor(:,1)
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
!       tmpx(:) = ( dedgd_r(1,:) - dedgd_r(1,:) ) / ( x(ix+1)-x(ix-1) )
!       tmpy(:) = ( dedgd_r(2,:) - dedgd_r(2,:) ) / ( y(iy+1)-y(iy-1) )
!       tmpz(:) = ( dedgd_r(3,:) - dedgd_r(3,:) ) / ( z(iz+1)-z(iz-1) )
  
       else
         tmpx_r(:) = 0.0_dp
         tmpy_r(:) = 0.0_dp
         tmpz_r(:) = 0.0_dp
       endif
  
 call start_clock(timing_tmp1)
       do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO COLLAPSE(2)
         do jbf=1,basis%nbf
           do ibf=1,basis%nbf
             vxc_ij(ibf,jbf,ispin) =  vxc_ij(ibf,jbf,ispin) + weight * fact_becke &
                 * eval_basis_function(basis%bf(ibf),rr) * eval_basis_function(basis%bf(jbf),rr) &
                 * ( dedd_r(ispin) - tmpx_r(ispin) - tmpy_r(ispin) - tmpz_r(ispin) )
           enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
       enddo
 call stop_clock(timing_tmp1)
  
     enddo
   enddo ! loop on the radial grid
!!!! !$OMP END DO
!!!! !$OMP END PARALLEL
  
 enddo ! loop on the atoms


#else
 write(*,*) 'XC energy and potential set to zero'
 write(*,*) 'libxc is not present'
 exc_xc=0.0_dp
 vxc_ij(:,:,:)=0.0_dp
#endif

! write(*,*)
! write(*,'(a,2(2x,f12.6))') 'Average Vxc',vxc_av
 write(*,*)
 write(*,'(a,2(2x,f12.6))') ' number of electrons:',normalization(:)
 write(*,'(a,2x,f12.6)')    '  DFT xc energy [Ha]:',exc_xc
 write(*,*)

contains

 function smooth_step(mu)
 real(dp) :: smooth_step
 real(dp),intent(in) :: mu
!=====

 smooth_step = 0.5_dp * mu * ( 3.0_dp - mu**2 )

 end function smooth_step

end subroutine dft_exc_vxc

