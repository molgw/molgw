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
 real(dp),parameter :: shift=1.d-5 ! bohr  some shift used
                                   ! to evaluate numerically the divergence of the gradient
 integer,parameter :: nx=40
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
 type(xc_f90_pointer_t) :: xc_func1,xc_func2
 type(xc_f90_pointer_t) :: xc_info1,xc_info2
#endif

 real(dp) :: basis_function_r       (basis%nbf)
 real(dp) :: basis_function_r_shiftx(basis%nbf)
 real(dp) :: basis_function_r_shifty(basis%nbf)
 real(dp) :: basis_function_r_shiftz(basis%nbf)
 real(dp) :: basis_function_gradr       (3,basis%nbf)
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
 real(dp) :: vxc1(nspin),vxc2(nspin)
 real(dp) :: vxc_dummy(nspin)
 real(dp) :: exc1,exc2
 real(dp) :: vsigma1       (2*nspin-1)
 real(dp) :: vsigma1_shiftx(2*nspin-1)
 real(dp) :: vsigma1_shifty(2*nspin-1)
 real(dp) :: vsigma1_shiftz(2*nspin-1)
 real(dp) :: vsigma2       (2*nspin-1)
 real(dp) :: vsigma2_shiftx(2*nspin-1)
 real(dp) :: vsigma2_shifty(2*nspin-1)
 real(dp) :: vsigma2_shiftz(2*nspin-1)
 real(dp) :: vxc_av(nspin)
 real(dp) :: dedd_r(nspin)
 real(dp) :: dedgd_r       (3,nspin)
 real(dp) :: dedgd_r_shiftx(3,nspin)
 real(dp) :: dedgd_r_shifty(3,nspin)
 real(dp) :: dedgd_r_shiftz(3,nspin)
 real(dp) :: div(nspin)

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
       if(xc_f90_info_family(xc_info1) == XC_FAMILY_GGA) then 
         do ibf=1,basis%nbf

           basis_function_gradr(:,ibf)        = eval_basis_function_grad(basis%bf(ibf),rr)

           rr_shift(:) = rr(:) + (/ shift , 0.0_dp , 0.0_dp /)
           basis_function_r_shiftx(ibf)       = eval_basis_function(basis%bf(ibf),rr_shift)
           basis_function_gradr_shiftx(:,ibf) = eval_basis_function_grad(basis%bf(ibf),rr_shift)
!!           write(*,*) '____________________'
!!           write(*,*) rr(:), basis_function_r(ibf)
!!           write(*,*) rr_shift(:), basis_function_r_shiftx(ibf)
!!           write(*,*) '____________________'

           rr_shift(:) = rr(:) + (/ 0.0_dp , shift , 0.0_dp /)
           basis_function_r_shifty(ibf)       = eval_basis_function(basis%bf(ibf),rr_shift)
           basis_function_gradr_shifty(:,ibf) = eval_basis_function_grad(basis%bf(ibf),rr_shift)

           rr_shift(:) = rr(:) + (/ 0.0_dp , 0.0_dp , shift /)
           basis_function_r_shiftz(ibf)       = eval_basis_function(basis%bf(ibf),rr_shift)
           basis_function_gradr_shiftz(:,ibf) = eval_basis_function_grad(basis%bf(ibf),rr_shift)

!           write(*,*) '==========',ibf
!           write(*,*) (  basis_function_r_shiftx(ibf) - basis_function_r(ibf) ) / shift
!           write(*,*) (  basis_function_r_shifty(ibf) - basis_function_r(ibf) ) / shift
!           write(*,*) (  basis_function_r_shiftz(ibf) - basis_function_r(ibf) ) / shift
!           write(*,*) basis_function_gradr(:,ibf)
         enddo
       endif
!       stop'ENOUGHr'


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
  
       if(xc_f90_info_family(xc_info1) == XC_FAMILY_GGA) then
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

               grad_rhor(:,ispin)       = grad_rhor(:,ispin)         + p_matrix(ibf,jbf,ispin) &
                    *( basis_function_gradr(:,ibf) * basis_function_r(jbf) &
                     + basis_function_gradr(:,jbf) * basis_function_r(ibf) ) 
               grad_rhor_shiftx(:,ispin) = grad_rhor_shiftx(:,ispin) + p_matrix(ibf,jbf,ispin) &
                    *( basis_function_gradr_shiftx(:,ibf) * basis_function_r(jbf) &
                     + basis_function_gradr_shiftx(:,jbf) * basis_function_r(ibf) ) 
               grad_rhor_shifty(:,ispin) = grad_rhor_shifty(:,ispin) + p_matrix(ibf,jbf,ispin) &
                    *( basis_function_gradr_shifty(:,ibf) * basis_function_r(jbf) &
                     + basis_function_gradr_shifty(:,jbf) * basis_function_r(ibf) ) 
               grad_rhor_shiftz(:,ispin) = grad_rhor_shiftz(:,ispin) + p_matrix(ibf,jbf,ispin) &
                    *( basis_function_gradr_shiftz(:,ibf) * basis_function_r(jbf) &
                     + basis_function_gradr_shiftz(:,jbf) * basis_function_r(ibf) ) 
             enddo
           enddo
         enddo
  
         sigma2(1)        = SUM( grad_rhor(:,1)**2 )
         sigma2_shiftx(1) = SUM( grad_rhor_shiftx(:,1)**2 )
         sigma2_shifty(1) = SUM( grad_rhor_shifty(:,1)**2 )
         sigma2_shiftz(1) = SUM( grad_rhor_shiftz(:,1)**2 )
         if(nspin==2) then
           sigma2(2)        = SUM( grad_rhor(:,1)*grad_rhor(:,2) )
           sigma2_shiftx(2) = SUM( grad_rhor_shiftx(:,1)*grad_rhor_shiftx(:,2) )
           sigma2_shifty(2) = SUM( grad_rhor_shifty(:,1)*grad_rhor_shifty(:,2) )
           sigma2_shiftz(2) = SUM( grad_rhor_shiftz(:,1)*grad_rhor_shiftz(:,2) )
           sigma2(3)        = SUM( grad_rhor(:,2)**2 )
           sigma2_shiftx(3) = SUM( grad_rhor_shiftx(:,2)**2 )
           sigma2_shifty(3) = SUM( grad_rhor_shifty(:,2)**2 )
           sigma2_shiftz(3) = SUM( grad_rhor_shiftz(:,2)**2 )
         endif

       endif
  
       !
       ! LIBXC call
       !
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
           call xc_f90_gga_exc_vxc(xc_func1,1,rhor_r(1)       ,sigma2(1)       ,exc1,vxc1(1),vsigma1(1)       )
           call xc_f90_gga_vxc    (xc_func1,1,rhor_r_shiftx(1),sigma2_shiftx(1),vxc_dummy(1),vsigma1_shiftx(1))
           call xc_f90_gga_vxc    (xc_func1,1,rhor_r_shifty(1),sigma2_shifty(1),vxc_dummy(1),vsigma1_shifty(1))
           call xc_f90_gga_vxc    (xc_func1,1,rhor_r_shiftz(1),sigma2_shiftz(1),vxc_dummy(1),vsigma1_shiftz(1))
         else
           exc1=0.0_dp
           vxc1=0.0_dp
           vsigma1=0.0_dp
         endif
         if(dft_xc(2)/=0) then
           call xc_f90_gga_exc_vxc(xc_func2,1,rhor_r(1)       ,sigma2(1)       ,exc2,vxc2(1),vsigma2(1)       )
           call xc_f90_gga_vxc    (xc_func2,1,rhor_r_shiftx(1),sigma2_shiftx(1),vxc_dummy(1),vsigma2_shiftx(1))
           call xc_f90_gga_vxc    (xc_func2,1,rhor_r_shifty(1),sigma2_shifty(1),vxc_dummy(1),vsigma2_shifty(1))
           call xc_f90_gga_vxc    (xc_func2,1,rhor_r_shiftz(1),sigma2_shiftz(1),vxc_dummy(1),vsigma2_shiftz(1))
         else
           exc2=0.0_dp
           vxc2=0.0_dp
           vsigma2=0.0_dp
         endif
       else
         stop'not LDA nor GGA is not implemented'
       endif
  
       exc_xc = exc_xc + weight * fact_becke * ( exc1 + exc2 ) * SUM( rhor_r(:) )
  
       dedd_r(:) = vxc1(:) + vxc2(:)
  
       if(xc_f90_info_family(xc_info1) == XC_FAMILY_GGA) then
         if(nspin==1) then
           dedgd_r       (:,1) = 2.0_dp * ( vsigma1(1)        + vsigma2(1)        ) * grad_rhor(:,1) 
           dedgd_r_shiftx(:,1) = 2.0_dp * ( vsigma1_shiftx(1) + vsigma2_shiftx(1) ) * grad_rhor_shiftx(:,1) 
           dedgd_r_shifty(:,1) = 2.0_dp * ( vsigma1_shifty(1) + vsigma2_shifty(1) ) * grad_rhor_shifty(:,1) 
           dedgd_r_shiftz(:,1) = 2.0_dp * ( vsigma1_shiftz(1) + vsigma2_shiftz(1) ) * grad_rhor_shiftz(:,1) 
         else
           dedgd_r(:,1)        = 2.0_dp * ( vsigma1(1)        + vsigma2(1)        ) * grad_rhor(:,1) &
                                        + ( vsigma1(2)        + vsigma2(2)        ) * grad_rhor(:,2)
           dedgd_r_shiftx(:,1) = 2.0_dp * ( vsigma1_shiftx(1) + vsigma2_shiftx(1) ) * grad_rhor_shiftx(:,1) &
                                        + ( vsigma1_shiftx(2) + vsigma2_shiftx(2) ) * grad_rhor_shiftx(:,2)
           dedgd_r_shifty(:,1) = 2.0_dp * ( vsigma1_shifty(1) + vsigma2_shifty(1) ) * grad_rhor_shifty(:,1) &
                                        + ( vsigma1_shifty(2) + vsigma2_shifty(2) ) * grad_rhor_shifty(:,2)
           dedgd_r_shiftz(:,1) = 2.0_dp * ( vsigma1_shiftz(1) + vsigma2_shiftz(1) ) * grad_rhor_shiftz(:,1) &
                                        + ( vsigma1_shiftz(2) + vsigma2_shiftz(2) ) * grad_rhor_shiftz(:,2)

           dedgd_r(:,2)        = 2.0_dp * ( vsigma1(3)        + vsigma2(3)        ) * grad_rhor(:,2) &
                                        + ( vsigma1(2)        + vsigma2(2)        ) * grad_rhor(:,1)
           dedgd_r_shiftx(:,2) = 2.0_dp * ( vsigma1_shiftx(3) + vsigma2_shiftx(3) ) * grad_rhor_shiftx(:,2) &
                                        + ( vsigma1_shiftx(2) + vsigma2_shiftx(2) ) * grad_rhor_shiftx(:,1)
           dedgd_r_shifty(:,2) = 2.0_dp * ( vsigma1_shifty(3) + vsigma2_shifty(3) ) * grad_rhor_shifty(:,2) &
                                        + ( vsigma1_shifty(2) + vsigma2_shifty(2) ) * grad_rhor_shifty(:,1)
           dedgd_r_shiftz(:,2) = 2.0_dp * ( vsigma1_shiftz(3) + vsigma2_shiftz(3) ) * grad_rhor_shiftz(:,2) &
                                        + ( vsigma1_shiftz(2) + vsigma2_shiftz(2) ) * grad_rhor_shiftz(:,1)

         endif
  
!         call xc_f90_gga_exc_vxc(xc_func1,1,rhor_r(1),sigma2(1),exc1,vxc1(1),vsigma1(1))
!         if(dft_xc(2)/=0) then
!           call xc_f90_gga_exc_vxc(xc_func2,1,rhor_r(1),sigma2(1),exc2,vxc2(1),vsigma2(1))
!         else
!           exc2=0.0_dp
!           vxc2=0.0_dp
!           vsigma2=0.0_dp
!         endif
  
         div(:) = ( dedgd_r_shiftx(1,:) - dedgd_r(1,:) ) / shift &
                + ( dedgd_r_shifty(2,:) - dedgd_r(2,:) ) / shift &
                + ( dedgd_r_shiftz(3,:) - dedgd_r(3,:) ) / shift
!         div(:) = ( dedgd_r(1,:) ) + (  dedgd_r(2,:) ) + (  dedgd_r(3,:) ) 
  
       else
         div(:) = 0.0_dp
       endif

!       if(div(1) > 1000.0) then 
!          write(*,*) '-------------'
!          write(*,*) 'WARNING large div',div(1)
!          write(*,*) 'for point',rr(:)
!          write(*,*) dedgd_r(1,:)
!          write(*,*) dedgd_r_shiftx(1,:)
!          write(*,*) dedgd_r_shifty(1,:)
!          write(*,*) dedgd_r_shiftz(1,:)
!          write(*,*) dedd_r(:)
!          write(*,*)
!       endif
!       stop'ENOUGH'
! div = div / 2.0 * 1.5
! div = 0.0
  
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

contains

 function smooth_step(mu)
 real(dp) :: smooth_step
 real(dp),intent(in) :: mu
!=====

 smooth_step = 0.5_dp * mu * ( 3.0_dp - mu**2 )

 end function smooth_step

end subroutine dft_exc_vxc

