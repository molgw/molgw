!=========================================================================
subroutine dft_exc_vxc(nspin,basis,dft_xc,p_matrix,ehomo,vxc_ij,exc_xc)
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
 real(dp),intent(in)        :: ehomo(nspin)
 real(dp),intent(out)       :: vxc_ij(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       :: exc_xc
!=====
 real(dp),parameter :: shift=1.d-6 ! bohr  some shift used
                                   ! to evaluate numerically the divergence of the gradient
 integer,parameter :: nx= 40
 integer,parameter :: nangular= 38 ! 86

 real(dp) :: weight
 real(dp) :: x1(nangular)
 real(dp) :: y1(nangular)
 real(dp) :: z1(nangular)
 real(dp) :: w1(nangular)
 integer  :: n1,iangular
 integer  :: ix,iy,iz,ibf,jbf,ispin,iatom,jatom,katom
 real(dp) :: rr(3),rr_shift(3)
 real(dp) :: xa(nx),wxa(nx)
 real(dp) :: normalization(nspin)

#ifdef HAVE_LIBXC
 type(xc_f90_pointer_t) :: xc_func1,xc_func2,xc_functest
 type(xc_f90_pointer_t) :: xc_info1,xc_info2,xc_infotest
#endif

 real(dp) :: basis_function_r       (basis%nbf)
 real(dp) :: basis_function_r_shiftx(basis%nbf)
 real(dp) :: basis_function_r_shifty(basis%nbf)
 real(dp) :: basis_function_r_shiftz(basis%nbf)
 real(dp) :: basis_function_gradr       (3,basis%nbf)
 real(dp) :: basis_function_laplr       (3,basis%nbf)
 real(dp) :: basis_function_gradr_shiftx(3,basis%nbf)
 real(dp) :: basis_function_gradr_shifty(3,basis%nbf)
 real(dp) :: basis_function_gradr_shiftz(3,basis%nbf)
 real(dp) :: rhor_r       (nspin)
 real(dp) :: rhor_r_shiftx(nspin)
 real(dp) :: rhor_r_shifty(nspin)
 real(dp) :: rhor_r_shiftz(nspin)
 real(dp) :: grad_rhor       (3,nspin)
 real(dp) :: grad_rhor_shiftx(3,nspin)
 real(dp) :: grad_rhor_shifty(3,nspin)
 real(dp) :: grad_rhor_shiftz(3,nspin)
 real(dp) :: sigma2       (2*nspin-1)
 real(dp) :: sigma2_shiftx(2*nspin-1)
 real(dp) :: sigma2_shifty(2*nspin-1)
 real(dp) :: sigma2_shiftz(2*nspin-1)
 real(dp) :: tau(nspin),lapl_rhor(nspin)
 real(dp) :: vxc1(nspin),vxc2(nspin)
 real(dp) :: vxc_dummy(nspin)
 real(dp) :: exc1(1),exc2(1)
 real(dp) :: vsigma1       (2*nspin-1)
 real(dp) :: vsigma1_shiftx(2*nspin-1)
 real(dp) :: vsigma1_shifty(2*nspin-1)
 real(dp) :: vsigma1_shiftz(2*nspin-1)
 real(dp) :: vsigma2       (2*nspin-1)
 real(dp) :: vsigma2_shiftx(2*nspin-1)
 real(dp) :: vsigma2_shifty(2*nspin-1)
 real(dp) :: vsigma2_shiftz(2*nspin-1)
 real(dp) :: vlapl_rho(nspin),vtau(nspin)
 real(dp) :: vxc_av(nspin)
 real(dp) :: dedd_r(nspin)
 real(dp) :: dedgd_r       (3,nspin)
 real(dp) :: dedgd_r_shiftx(3,nspin)
 real(dp) :: dedgd_r_shifty(3,nspin)
 real(dp) :: dedgd_r_shiftz(3,nspin)
 real(dp) :: div(nspin)
 real(dp) :: mu,s_becke(natom,natom),p_becke(natom),fact_becke
 real(dp) :: rtmp
 character(len=256) :: string
!=====

 exc_xc = 0.0_dp
 vxc_ij(:,:,:) = 0.0_dp
 if( ALL(dft_xc(:)==0) ) return

#ifdef HAVE_LIBXC

#if 0
 call xc_f90_func_init(xc_functest, xc_infotest, XC_LDA_X, XC_UNPOLARIZED)
! call xc_f90_func_init(xc_functest, xc_infotest, XC_LDA_C_VWN_RPA , XC_UNPOLARIZED)
 do ix=1,200
   exc2(1) = exp(0.08*(DBLE(ix)-1.0))*0.05
   rhor_r(1)= 3.0/ (4.0*pi*exc2(1)**3)
   call xc_f90_lda_exc_vxc(xc_functest,1,rhor_r(1),exc1(1),vxc1(1))
   write(105,'(10(e16.8,2x))') exc2(1),rhor_r(1),exc1(1),vxc1(1)
 enddo 
 stop'ENOUGH'
#endif

 write(*,*) 'Evaluate DFT integrals'
 write(*,'(a,i4,x,i4)') '   discretization grid per atom [radial points , angular points] ',nx,nangular
 write(*,'(a,i8)')      '   total number of real-space points',nx*nangular*natom
 
 if( dft_xc(1) < 1000 ) then
   if(nspin==1) then
     call xc_f90_func_init(xc_func1, xc_info1, dft_xc(1), XC_UNPOLARIZED)
     call xc_f90_func_init(xc_func2, xc_info2, dft_xc(2), XC_UNPOLARIZED)
   else
     call xc_f90_func_init(xc_func1, xc_info1, dft_xc(1), XC_POLARIZED)
     call xc_f90_func_init(xc_func2, xc_info2, dft_xc(2), XC_POLARIZED)
   endif
 else
   write(*,*) 'Home-made functional'
   if(nspin==1) then
     call xc_f90_func_init(xc_func1, xc_info1, XC_LDA_X, XC_UNPOLARIZED)
     call xc_f90_func_init(xc_func2, xc_info2, 0, XC_UNPOLARIZED)
   else
     call xc_f90_func_init(xc_func1, xc_info1, XC_LDA_X, XC_POLARIZED)
     call xc_f90_func_init(xc_func2, xc_info2, 0, XC_POLARIZED)
   endif
 endif
! write(*,*) 'LIBXC functional index',dft_xc(:)
! write(*,*) xc_f90_info_kind(xc_info1)
! write(*,*) 'name   ',TRIM(string)
! call xc_f90_hyb_gga_exx_coef(xc_func1,rtmp)
! write(*,*) 'exx',rtmp
 write(*,'(/,a)') ' LIBXC info'
 if( dft_xc(1) /=0 ) then
   call xc_f90_info_name(xc_info1,string)
   write(*,'(a,i6,5x,a)') '   XC functional 1: ', xc_f90_info_number(xc_info1),TRIM(string)
 endif
 if( dft_xc(2) /=0 ) then
   call xc_f90_info_name(xc_info2,string)
   write(*,'(a,i6,5x,a)') '   XC functional 2: ', xc_f90_info_number(xc_info2),TRIM(string)
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
  
  
       !
       ! first calculate all the needed basis function evaluations at point rr
       do ibf=1,basis%nbf
         basis_function_r(ibf) = eval_basis_function(basis%bf(ibf),rr)
       enddo

       if(xc_f90_info_family(xc_info1) == XC_FAMILY_GGA .OR. xc_f90_info_family(xc_info1) == XC_FAMILY_HYB_GGA) then 
         do ibf=1,basis%nbf

           basis_function_gradr(:,ibf)        = eval_basis_function_grad(basis%bf(ibf),rr)

           rr_shift(:) = rr(:) + (/ shift , 0.0_dp , 0.0_dp /)
           basis_function_r_shiftx(ibf)       = eval_basis_function(basis%bf(ibf),rr_shift)
           basis_function_gradr_shiftx(:,ibf) = eval_basis_function_grad(basis%bf(ibf),rr_shift)

           rr_shift(:) = rr(:) + (/ 0.0_dp , shift , 0.0_dp /)
           basis_function_r_shifty(ibf)       = eval_basis_function(basis%bf(ibf),rr_shift)
           basis_function_gradr_shifty(:,ibf) = eval_basis_function_grad(basis%bf(ibf),rr_shift)

           rr_shift(:) = rr(:) + (/ 0.0_dp , 0.0_dp , shift /)
           basis_function_r_shiftz(ibf)       = eval_basis_function(basis%bf(ibf),rr_shift)
           basis_function_gradr_shiftz(:,ibf) = eval_basis_function_grad(basis%bf(ibf),rr_shift)

         enddo
       endif

       if(xc_f90_info_family(xc_info1) == XC_FAMILY_MGGA ) then
         do ibf=1,basis%nbf
           basis_function_gradr(:,ibf)        = eval_basis_function_grad(basis%bf(ibf),rr)
           basis_function_laplr(:,ibf)        = eval_basis_function_lapl(basis%bf(ibf),rr)
         enddo
       endif


       rhor_r(:)=0.0_dp
       do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO REDUCTION(+:rhor_r) COLLAPSE(2)
         do jbf=1,basis%nbf
           do ibf=1,basis%nbf
             rhor_r(ispin)=rhor_r(ispin)+p_matrix(ibf,jbf,ispin)&
                               * basis_function_r(ibf) &
                               * basis_function_r(jbf) 
           enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
       enddo
  
       normalization(:) = normalization(:) + rhor_r(:) * weight * fact_becke
  
       if(xc_f90_info_family(xc_info1) == XC_FAMILY_GGA .OR. xc_f90_info_family(xc_info1) == XC_FAMILY_HYB_GGA) then
         grad_rhor       (:,:)=0.0_dp
         grad_rhor_shiftx(:,:)=0.0_dp
         grad_rhor_shifty(:,:)=0.0_dp
         grad_rhor_shiftz(:,:)=0.0_dp
         rhor_r_shiftx(:)     =0.0_dp
         rhor_r_shifty(:)     =0.0_dp
         rhor_r_shiftz(:)     =0.0_dp
         do ispin=1,nspin
           do jbf=1,basis%nbf
             do ibf=1,basis%nbf

               rhor_r_shiftx(ispin) = rhor_r_shiftx(ispin) + p_matrix(ibf,jbf,ispin)&
                                        * basis_function_r_shiftx(ibf) &
                                        * basis_function_r_shiftx(jbf) 
               rhor_r_shifty(ispin) = rhor_r_shifty(ispin) + p_matrix(ibf,jbf,ispin)&
                                        * basis_function_r_shifty(ibf) &
                                        * basis_function_r_shifty(jbf) 
               rhor_r_shiftz(ispin) = rhor_r_shiftz(ispin) + p_matrix(ibf,jbf,ispin)&
                                        * basis_function_r_shiftz(ibf) &
                                        * basis_function_r_shiftz(jbf) 

               grad_rhor(:,ispin)        = grad_rhor(:,ispin)        + p_matrix(ibf,jbf,ispin) &
                    *( basis_function_gradr       (:,ibf) * basis_function_r(jbf) &
                     + basis_function_gradr       (:,jbf) * basis_function_r(ibf) ) 
               grad_rhor_shiftx(:,ispin) = grad_rhor_shiftx(:,ispin) + p_matrix(ibf,jbf,ispin) &
                    *( basis_function_gradr_shiftx(:,ibf) * basis_function_r_shiftx(jbf) &
                     + basis_function_gradr_shiftx(:,jbf) * basis_function_r_shiftx(ibf) ) 
               grad_rhor_shifty(:,ispin) = grad_rhor_shifty(:,ispin) + p_matrix(ibf,jbf,ispin) &
                    *( basis_function_gradr_shifty(:,ibf) * basis_function_r_shifty(jbf) &
                     + basis_function_gradr_shifty(:,jbf) * basis_function_r_shifty(ibf) ) 
               grad_rhor_shiftz(:,ispin) = grad_rhor_shiftz(:,ispin) + p_matrix(ibf,jbf,ispin) &
                    *( basis_function_gradr_shiftz(:,ibf) * basis_function_r_shiftz(jbf) &
                     + basis_function_gradr_shiftz(:,jbf) * basis_function_r_shiftz(ibf) ) 

             enddo
           enddo
         enddo

         sigma2(1)        = SUM( grad_rhor       (:,1)**2 )
         sigma2_shiftx(1) = SUM( grad_rhor_shiftx(:,1)**2 )
         sigma2_shifty(1) = SUM( grad_rhor_shifty(:,1)**2 )
         sigma2_shiftz(1) = SUM( grad_rhor_shiftz(:,1)**2 )
         if(nspin==2) then
           sigma2(2)        = SUM( grad_rhor       (:,1) * grad_rhor       (:,2) )
           sigma2_shiftx(2) = SUM( grad_rhor_shiftx(:,1) * grad_rhor_shiftx(:,2) )
           sigma2_shifty(2) = SUM( grad_rhor_shifty(:,1) * grad_rhor_shifty(:,2) )
           sigma2_shiftz(2) = SUM( grad_rhor_shiftz(:,1) * grad_rhor_shiftz(:,2) )

           sigma2(3)        = SUM( grad_rhor       (:,2)**2 )
           sigma2_shiftx(3) = SUM( grad_rhor_shiftx(:,2)**2 )
           sigma2_shifty(3) = SUM( grad_rhor_shifty(:,2)**2 )
           sigma2_shiftz(3) = SUM( grad_rhor_shiftz(:,2)**2 )
         endif

       endif
  
       if(xc_f90_info_family(xc_info1) == XC_FAMILY_MGGA) then

         grad_rhor(:,:)=0.0_dp
         tau(:)        =0.0_dp
         lapl_rhor(:)  =0.0_dp
         do ispin=1,nspin
           do jbf=1,basis%nbf
             do ibf=1,basis%nbf

               grad_rhor(:,ispin) = grad_rhor(:,ispin)        + p_matrix(ibf,jbf,ispin) &
                    *( basis_function_gradr       (:,ibf) * basis_function_r(jbf) &
                     + basis_function_gradr       (:,jbf) * basis_function_r(ibf) ) 

               tau(ispin)        = tau(ispin)        + p_matrix(ibf,jbf,ispin) &
                    * DOT_PRODUCT( basis_function_gradr(:,ibf) , basis_function_gradr(:,jbf) )

               lapl_rhor(ispin)  = lapl_rhor(ispin)  + p_matrix(ibf,jbf,ispin) &
                                  * (  SUM( basis_function_laplr(:,ibf) ) * basis_function_r(jbf)               &
                                     + basis_function_r(ibf)              * SUM( basis_function_laplr(:,jbf) )  &
                                     + 2.0_dp * DOT_PRODUCT( basis_function_gradr(:,ibf) , basis_function_gradr(:,jbf) ) )

             enddo
           enddo
         enddo
         sigma2(1)          = SUM( grad_rhor       (:,1)**2 )
         if(nspin==2) then
           sigma2(2)        = SUM( grad_rhor       (:,1) * grad_rhor       (:,2) )
           sigma2(3)        = SUM( grad_rhor       (:,2)**2 )
         endif

       endif

       !
       ! LIBXC call
       !
       select case(xc_f90_info_family(xc_info1))
       case(XC_FAMILY_LDA)
         if(dft_xc(1)/=0) then
           if( dft_xc(1) < 1000 ) then 
             call xc_f90_lda_exc_vxc(xc_func1,1,rhor_r(1),exc1(1),vxc1(1))
           else
             call my_lda_exc_vxc(nspin,dft_xc(1),rhor_r,exc1(1),vxc1)
           endif
         else
           exc1(1)=0.0_dp
           vxc1(1)=0.0_dp
         endif
         if(dft_xc(2)/=0) then
           call xc_f90_lda_exc_vxc(xc_func2,1,rhor_r(1),exc2(1),vxc2(1))
         else
           exc2(1)=0.0_dp
           vxc2(1)=0.0_dp
         endif

       case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
         if(dft_xc(1)/=0) then
           call xc_f90_gga_exc_vxc(xc_func1,1,rhor_r(1)       ,sigma2(1)       ,exc1(1),vxc1(1),vsigma1(1)       )
           call xc_f90_gga_vxc    (xc_func1,1,rhor_r_shiftx(1),sigma2_shiftx(1),vxc_dummy(1),vsigma1_shiftx(1))
           call xc_f90_gga_vxc    (xc_func1,1,rhor_r_shifty(1),sigma2_shifty(1),vxc_dummy(1),vsigma1_shifty(1))
           call xc_f90_gga_vxc    (xc_func1,1,rhor_r_shiftz(1),sigma2_shiftz(1),vxc_dummy(1),vsigma1_shiftz(1))
         else
           exc1(1)=0.0_dp
           vxc1(1)=0.0_dp
           vsigma1=0.0_dp
           vsigma1_shiftx=0.0_dp
           vsigma1_shifty=0.0_dp
           vsigma1_shiftz=0.0_dp
         endif
         if(dft_xc(2)/=0) then
           call xc_f90_gga_exc_vxc(xc_func2,1,rhor_r(1)       ,sigma2(1)       ,exc2(1),vxc2(1),vsigma2(1)       )
           call xc_f90_gga_vxc    (xc_func2,1,rhor_r_shiftx(1),sigma2_shiftx(1),vxc_dummy(1),vsigma2_shiftx(1))
           call xc_f90_gga_vxc    (xc_func2,1,rhor_r_shifty(1),sigma2_shifty(1),vxc_dummy(1),vsigma2_shifty(1))
           call xc_f90_gga_vxc    (xc_func2,1,rhor_r_shiftz(1),sigma2_shiftz(1),vxc_dummy(1),vsigma2_shiftz(1))
         else
           exc2(1)=0.0_dp
           vxc2(1)=0.0_dp
           vsigma2=0.0_dp
           vsigma2_shiftx=0.0_dp
           vsigma2_shifty=0.0_dp
           vsigma2_shiftz=0.0_dp
         endif

       case(XC_FAMILY_MGGA)
         if(dft_xc(1)/=0) then
           call xc_f90_mgga_vxc(xc_func1,1,rhor_r(1),sigma2(1),lapl_rhor(1),tau(1),vxc1(1),vsigma1(1),vlapl_rho(1),vtau(1))
           exc1(1)=0.0_dp
         else
           exc1(1)=0.0_dp
           vxc1(1)=0.0_dp
           vsigma1=0.0_dp
         endif
         if(dft_xc(2)/=0) then
           call xc_f90_mgga_vxc(xc_func2,1,rhor_r(1),sigma2(1),lapl_rhor(1),tau(1),vxc2(1),vsigma1(1),vlapl_rho(1),vtau(1))
           exc2(1)=0.0_dp
         else
           exc2(1)=0.0_dp
           vxc2(1)=0.0_dp
           vsigma2=0.0_dp
         endif

       case default
         stop'functional is not LDA nor GGA nor hybrid nor meta-GGA'
       end select
  
       exc_xc = exc_xc + weight * fact_becke * ( exc1(1) + exc2(1) ) * SUM( rhor_r(:) )
  
       dedd_r(:) = vxc1(:) + vxc2(:)

       if(xc_f90_info_family(xc_info1) == XC_FAMILY_GGA .OR. xc_f90_info_family(xc_info1) == XC_FAMILY_HYB_GGA) then
         if(nspin==1) then

           dedgd_r       (:,1) = 2.0_dp * ( vsigma1(1)        + vsigma2(1)        ) * grad_rhor       (:,1) 
           dedgd_r_shiftx(:,1) = 2.0_dp * ( vsigma1_shiftx(1) + vsigma2_shiftx(1) ) * grad_rhor_shiftx(:,1) 
           dedgd_r_shifty(:,1) = 2.0_dp * ( vsigma1_shifty(1) + vsigma2_shifty(1) ) * grad_rhor_shifty(:,1) 
           dedgd_r_shiftz(:,1) = 2.0_dp * ( vsigma1_shiftz(1) + vsigma2_shiftz(1) ) * grad_rhor_shiftz(:,1) 

         else

           dedgd_r(:,1)        = 2.0_dp * ( vsigma1(1)        + vsigma2(1)        ) * grad_rhor       (:,1) &
                                        + ( vsigma1(2)        + vsigma2(2)        ) * grad_rhor       (:,2)
           dedgd_r_shiftx(:,1) = 2.0_dp * ( vsigma1_shiftx(1) + vsigma2_shiftx(1) ) * grad_rhor_shiftx(:,1) &
                                        + ( vsigma1_shiftx(2) + vsigma2_shiftx(2) ) * grad_rhor_shiftx(:,2)
           dedgd_r_shifty(:,1) = 2.0_dp * ( vsigma1_shifty(1) + vsigma2_shifty(1) ) * grad_rhor_shifty(:,1) &
                                        + ( vsigma1_shifty(2) + vsigma2_shifty(2) ) * grad_rhor_shifty(:,2)
           dedgd_r_shiftz(:,1) = 2.0_dp * ( vsigma1_shiftz(1) + vsigma2_shiftz(1) ) * grad_rhor_shiftz(:,1) &
                                        + ( vsigma1_shiftz(2) + vsigma2_shiftz(2) ) * grad_rhor_shiftz(:,2)

           dedgd_r(:,2)        = 2.0_dp * ( vsigma1(3)        + vsigma2(3)        ) * grad_rhor       (:,2) &
                                        + ( vsigma1(2)        + vsigma2(2)        ) * grad_rhor       (:,1)
           dedgd_r_shiftx(:,2) = 2.0_dp * ( vsigma1_shiftx(3) + vsigma2_shiftx(3) ) * grad_rhor_shiftx(:,2) &
                                        + ( vsigma1_shiftx(2) + vsigma2_shiftx(2) ) * grad_rhor_shiftx(:,1)
           dedgd_r_shifty(:,2) = 2.0_dp * ( vsigma1_shifty(3) + vsigma2_shifty(3) ) * grad_rhor_shifty(:,2) &
                                        + ( vsigma1_shifty(2) + vsigma2_shifty(2) ) * grad_rhor_shifty(:,1)
           dedgd_r_shiftz(:,2) = 2.0_dp * ( vsigma1_shiftz(3) + vsigma2_shiftz(3) ) * grad_rhor_shiftz(:,2) &
                                        + ( vsigma1_shiftz(2) + vsigma2_shiftz(2) ) * grad_rhor_shiftz(:,1)

         endif
  
  
         div(:) = ( dedgd_r_shiftx(1,:) - dedgd_r(1,:) ) / shift &
                + ( dedgd_r_shifty(2,:) - dedgd_r(2,:) ) / shift &
                + ( dedgd_r_shiftz(3,:) - dedgd_r(3,:) ) / shift

  
       else
         div(:) = 0.0_dp
       endif
 

       !
       ! In the case of the BJ06 meta-GGA functional, a spin-dependent shift is applied
       ! since the potential does not vanish at infinity
       !
       if( dft_xc(1) == XC_MGGA_X_BJ06 ) then
         dedd_r(:) = dedd_r(:) - SQRT( 5.0_dp * ABS(ehomo(:)) / 6.0_dp ) / pi
       endif

  
       do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO COLLAPSE(2)
         do jbf=1,basis%nbf
           do ibf=1,basis%nbf
             vxc_ij(ibf,jbf,ispin) =  vxc_ij(ibf,jbf,ispin) + weight * fact_becke &
                 * ( dedd_r(ispin) - div(ispin) ) * basis_function_r(ibf) * basis_function_r(jbf)  
           enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
       enddo
  

     enddo ! loop on the angular grid
   enddo ! loop on the radial grid
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

 !
 ! Destroy operations
 if( dft_xc(1) /= 0 ) call xc_f90_func_end(xc_func1)
 if( dft_xc(2) /= 0 ) call xc_f90_func_end(xc_func2)

contains

 function smooth_step(mu)
 real(dp) :: smooth_step
 real(dp),intent(in) :: mu
!=====

 smooth_step = 0.5_dp * mu * ( 3.0_dp - mu**2 )

 end function smooth_step

end subroutine dft_exc_vxc

!=========================================================================
subroutine my_lda_exc_vxc(nspin,ixc,rhor,exc,vxc)
 use m_definitions
 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: nspin,ixc
 real(dp),intent(in) :: rhor(nspin)
!arrays
 real(dp),intent(out) :: exc,vxc(nspin)

 real(dp) :: a0p
 real(dp) :: a1p
 real(dp) :: a2p
 real(dp) :: a3p
 real(dp) :: b1p
 real(dp) :: b2p
 real(dp) :: b3p
 real(dp) :: b4p
 real(dp) :: da0
 real(dp) :: da1
 real(dp) :: da2
 real(dp) :: da3
 real(dp) :: db1
 real(dp) :: db2
 real(dp) :: db3
 real(dp) :: db4

 real(dp),parameter :: alpha_zeta=1.0_dp-1.0d-6
 real(dp),parameter :: ft=4._dp/3._dp,rsfac=0.6203504908994000_dp
 real(dp),parameter :: rsfacm3=rsfac**(-3)
 real(dp) :: a0,a1,a2,a3,b1,b2,b3,b4,d1,d1m1,d2d1drs2,d2d1drsdf,d2excdf2
 real(dp) :: d2excdrs2,d2excdrsdf,d2excdz2,d2fxcdz2,d2n1drs2,d2n1drsdf,dd1df
 real(dp) :: dd1drs,dexcdf,dexcdrs,dexcdz,dfxcdz,dn1df,dn1drs,dvxcdrs
 real(dp) :: dvxcpdrho,dvxcpdz,fact,fxc,n1
 real(dp) :: rhom1,rs,vxcp,zet,zetm,zetm_third
 real(dp) :: zetp,zetp_third
!no_abirules
!Set a minimum rho below which terms are 0
 real(dp),parameter :: rhotol=1.d-28
!real(dp) :: delta,rho,rho_dn,rho_dnm,rho_dnp,rho_up,rho_upm,rho_upp,zeta_mean

! *************************************************************************

 select case(ixc)
 case(1100)
   !
   ! The usual full LDA parameters of Teter
   a0p=.4581652932831429_dp
   a1p=2.217058676663745_dp
   a2p=0.7405551735357053_dp
   a3p=0.01968227878617998_dp
   b1p=1.0_dp
   b2p=4.504130959426697_dp
   b3p=1.110667363742916_dp
   b4p=0.02359291751427506_dp
   da0=.119086804055547_dp
   da1=0.6157402568883345_dp
   da2=0.1574201515892867_dp
   da3=0.003532336663397157_dp
   db1=0.0_dp
   db2=0.2673612973836267_dp
   db3=0.2052004607777787_dp
   db4=0.004200005045691381_dp
 case(1000)
   !
   ! full range RPA 
   a0p= 0.00033959499_dp
   a1p= 0.1912460_dp
   a2p= 0.8008790_dp
   a3p= 0.092956297_dp
   b1p=1.0_dp
   b2p= 8.710930_dp
   b3p= 3.945600_dp
   b4p= 0.087989897_dp
   da0=-0.00015974799_dp
   da1=-0.082753003_dp
   da2=-0.3675560_dp
   da3=-0.044320997_dp
   db1=0.0_dp
   db2=-1.113190_dp
   db3=-1.221470_dp
   db4=-0.040220797_dp
 case(1020)
   !
   ! Long-range only RPA parameters with rc=2.0
   a0p=-0.000012396600_dp
   a1p= 0.0014478000_dp
   a2p= 0.021771301_dp
   a3p= 0.00025572101_dp
   b1p=1.0_dp
   b2p= 0.3820980_dp
   b3p= 0.042663701_dp
   b4p= 0.00010346600_dp
   da0= 0.0000018310002_dp
   da1=-0.00021740992_dp
   da2=-0.0077045010_dp
   da3=-0.00020484751_dp
   db1=0.0_dp
   db2= 0.021046013_dp
   db3=-0.018320801_dp
   db4=-0.00012402580_dp
 case(1010)
   !
   ! Long-range only RPA parameters with rc=1.0
   a0p=-0.000028384500_dp
   a1p= 0.0037404201_dp
   a2p= 0.063176401_dp
   a3p= 0.0023404600_dp
   b1p=1.0_dp
   b2p= 0.8482450_dp
   b3p= 0.1845470_dp
   b4p= 0.0016536200_dp
   da0= 0.0000059325994_dp
   da1=-0.00076668011_dp
   da2=-0.024234399_dp
   da3=-0.0014384059_dp
   db1=0.0_dp
   db2= 0.025729001_dp
   db3=-0.071891010_dp
   db4=-0.0010981541_dp
 case(1005)
   !
   ! Long-range only RPA parameters with rc=0.5
   a0p=-5.8032401E-05
   a1p= 8.9607602E-03
   a2p= 0.1718570
   a3p= 1.3439300E-02
   b1p=1.0_dp
   b2p= 1.849290
   b3p= 0.7096860
   b4p= 1.1346900E-02
   da0= 1.3963599E-05
   da1= -2.1155602E-03
   da2= -7.3816001E-02
   da3= -7.0218993E-03
   db1=0.0_dp
   db2=-7.2860003E-02
   db3=-0.2463360
   db4=-5.8700801E-03
 end select

!Although fact is parameter value, some compilers are not able to evaluate
!it at compile time.
 fact=1.0_dp/(2.0_dp**(4.0_dp/3.0_dp)-2.0_dp)


 if (nspin==1) then

       rs=( 3.0_dp / (4.0_dp*pi*rhor(1)) )**(1.0_dp/3.0_dp) 
       n1=a0p+rs*(a1p+rs*(a2p+rs*a3p))
       d1=rs*(b1p+rs*(b2p+rs*(b3p+rs*b4p)))
       d1m1=1.0_dp/d1

!      Exchange-correlation energy
       exc=-n1*d1m1

!      Exchange-correlation potential
       dn1drs=a1p+rs*(2._dp*a2p+rs*(3._dp*a3p))
       dd1drs=b1p+rs*(2._dp*b2p+rs*(3._dp*b3p+rs*(4._dp*b4p)))

!      dexcdrs is d(exc)/d(rs)
       dexcdrs=-(dn1drs+exc*dd1drs)*d1m1
       vxc(1)=exc-rs*dexcdrs/3.0_dp

 else 

!    Allow for spin polarization. This part could be optimized for speed.

       rs=( 3.0_dp / (4.0_dp*pi*SUM(rhor(:))) )**(1.0_dp/3.0_dp) 
       zet= ( rhor(1) - rhor(2) ) / SUM( rhor(:) )
       zetp=1.0_dp+zet*alpha_zeta
       zetm=1.0_dp-zet*alpha_zeta
       zetp_third=zetp**(1.0_dp/3.0_dp)
       zetm_third=zetm**(1.0_dp/3.0_dp)
!      Exchange energy spin interpolation function f(zeta)
       fxc=( zetp*zetp_third + zetm*zetm_third - 2.0_dp ) *fact

       a0=a0p+fxc*da0
       a1=a1p+fxc*da1
       a2=a2p+fxc*da2
       a3=a3p+fxc*da3
       b1=b1p+fxc*db1
       b2=b2p+fxc*db2
       b3=b3p+fxc*db3
       b4=b4p+fxc*db4

       n1= a0+rs*(a1+rs*(a2+rs*a3))
       d1=rs*(b1+rs*(b2+rs*(b3+rs*b4)))
       d1m1=1.0_dp/d1

!      Exchange-correlation energy
       exc=-n1*d1m1

!      Exchange-correlation potential
       dn1drs=a1+rs*(2._dp*a2+rs*(3._dp*a3))
       dd1drs=b1+rs*(2._dp*b2+rs*(3._dp*b3+rs*(4._dp*b4)))
!      dexcdrs is d(exc)/d(rs)
       dexcdrs=-(dn1drs+exc*dd1drs)*d1m1

!      Only vxcp contributes when paramagnetic
       vxcp=exc-rs*dexcdrs/3.0_dp

!      d(fxc)/d(zeta)  (which is 0 at zeta=0)
       dfxcdz=ft*alpha_zeta*(zetp_third-zetm_third)*fact

!      dn1df=d(n1)/d(fxc) and dd1df=d(d1)/d(fxc)
       dn1df=da0+rs*(da1+rs*(da2+rs*da3))
       dd1df=rs*(db1+rs*(db2+rs*(db3+rs*db4)))

!      dexcdz is d(exc)/d(zeta)
       dexcdf=-(dn1df+exc*dd1df)*d1m1
       dexcdz=dfxcdz*dexcdf

!      Compute Vxc for both spin channels

       vxc(1)=vxcp - (zet-1.0_dp)*dexcdz
       vxc(2)=vxcp - (zet+1.0_dp)*dexcdz

 end if





end subroutine my_lda_exc_vxc
!=========================================================================
