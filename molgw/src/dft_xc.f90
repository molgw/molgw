!=========================================================================
subroutine dft_exc_vxc(basis,p_matrix,ehomo,vxc_ij,exc_xc)
 use m_definitions
 use m_mpi
 use m_timing
 use m_inputparam
 use m_basis_set
 use m_dft_grid
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
#ifdef _OPENMP
 use,intrinsic :: omp_lib
#endif
 use,intrinsic ::  iso_c_binding, only: C_INT,C_DOUBLE
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: p_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in)        :: ehomo(nspin)
 real(dp),intent(out)       :: vxc_ij(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)       :: exc_xc
!=====

 integer  :: idft_xc
 logical  :: require_gradient,require_laplacian
 integer  :: igrid,ibf,jbf,ispin
 real(dp) :: rr(3)
 real(dp) :: normalization(nspin)
 real(dp) :: weight

#ifdef HAVE_LIBXC
 type(xc_f90_pointer_t) :: xc_func(ndft_xc),xc_functest
 type(xc_f90_pointer_t) :: xc_info(ndft_xc),xc_infotest
#endif

 real(dp)             :: basis_function_r(basis%nbf)
 real(dp)             :: basis_function_gradr(3,basis%nbf)
 real(dp)             :: basis_function_laplr(3,basis%nbf)

 real(dp)             :: rhor_r(nspin)
 real(dp)             :: grad_rhor(3,nspin)
 real(dp)             :: sigma(2*nspin-1)
 real(dp)             :: tau(nspin),lapl_rhor(nspin)
 real(dp)             :: vxc_libxc(nspin)
 real(dp)             :: vxc_dummy(nspin)
 real(dp)             :: exc_libxc(1)
 real(dp)             :: vsigma(2*nspin-1)
 real(dp)             :: vlapl_rho(nspin),vtau(nspin)
 real(dp)             :: vxc_av(nspin)
 real(dp)             :: dedd_r(nspin)
 real(dp)             :: dedgd_r(3,nspin)
 real(dp)             :: omega
 character(len=256)   :: string
!=====

 exc_xc = 0.0_dp
 vxc_ij(:,:,:) = 0.0_dp
 if( ndft_xc == 0 ) return

 call start_clock(timing_dft)


#ifdef HAVE_LIBXC

 write(stdout,*) 'Calculate DFT XC potential'
 
 require_gradient =.FALSE.
 require_laplacian=.FALSE.
 do idft_xc=1,ndft_xc

   if( dft_xc_type(idft_xc) < 1000 ) then
     if(nspin==1) then
       call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), dft_xc_type(idft_xc), XC_UNPOLARIZED)
     else
       call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), dft_xc_type(idft_xc), XC_POLARIZED)
     endif
   else if(dft_xc_type(idft_xc) < 2000) then
     write(stdout,*) 'Home-made functional LDA functional'
     ! Fake LIBXC descriptor 
     if(nspin==1) then
       call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), XC_LDA_X, XC_UNPOLARIZED)
     else
       call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), XC_LDA_X, XC_POLARIZED)
     endif
   else
     write(stdout,*) 'Home-made functional GGA functional'
     ! Fake LIBXC descriptor 
     if(nspin==1) then
       call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), XC_GGA_X_PBE, XC_UNPOLARIZED)
     else
       call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), XC_GGA_X_PBE, XC_POLARIZED)
    endif
   endif

   if( dft_xc_type(idft_xc) < 1000 ) then
     call xc_f90_info_name(xc_info(idft_xc),string)
     write(stdout,'(a,i4,a,i6,5x,a)') '   XC functional ',idft_xc,' :  ',xc_f90_info_number(xc_info(idft_xc)),TRIM(string)
   else
     write(stdout,'(a,i4,a,i6,5x,a)') '   XC functional ',idft_xc,' :  ',xc_f90_info_number(xc_info(idft_xc)),'FAKE LIBXC DESCRIPTOR'
   endif

   if(xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_GGA     ) require_gradient  =.TRUE.
   if(xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_HYB_GGA ) require_gradient  =.TRUE.
   if(xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_MGGA    ) require_laplacian =.TRUE.

   if( dft_xc_type(idft_xc) == XC_GGA_X_HJS_PBE ) then
     call xc_f90_gga_x_hjs_set_par(xc_func(idft_xc), REAL(gamma_hybrid,C_DOUBLE) )
   endif
   if( dft_xc_type(idft_xc) == XC_GGA_X_WPBEH ) then
     call xc_f90_gga_x_wpbeh_set_par(xc_func(idft_xc), REAL(gamma_hybrid,C_DOUBLE) )
   endif

 enddo


 !
 ! If it is the first time, then set up the stored arrays
 !
 if( .NOT. ALLOCATED(bfr) )                          call prepare_basis_functions_r(basis)
 if( require_gradient  .AND. .NOT. ALLOCATED(bfgr) ) call prepare_basis_functions_gradr(basis)
 if( require_laplacian .AND. .NOT. ALLOCATED(bfgr) ) call prepare_basis_functions_laplr(basis)

 normalization(:)=0.0_dp


 do igrid=1,ngrid

   rr(:) = rr_grid(:,igrid)
   weight = w_grid(igrid)

   !
   ! Get all the functions and gradients at point rr
   call get_basis_functions_r(basis,igrid,basis_function_r)
   if( require_gradient ) then
     call get_basis_functions_gradr(basis,igrid,basis_function_gradr)
   endif
   if( require_laplacian ) then
     call get_basis_functions_laplr(basis,igrid,basis_function_gradr,basis_function_laplr)
   endif


   !
   ! calculate the density at point r for spin up and spin down
   call calc_density_r(nspin,basis,p_matrix,rr,basis_function_r,rhor_r)
   !
   ! Normalization
   normalization(:) = normalization(:) + rhor_r(:) * weight


   if( require_gradient ) then 
     call calc_density_gradr(nspin,basis%nbf,p_matrix,basis_function_r,basis_function_gradr,grad_rhor)

     sigma(1) = SUM( grad_rhor(:,1)**2 )
     if(nspin==2) then
       sigma(2) = SUM( grad_rhor(:,1) * grad_rhor(:,2) )
       sigma(3) = SUM( grad_rhor(:,2)**2 )
     endif

   endif


   if( require_laplacian ) then

     grad_rhor(:,:)=0.0_dp
     tau(:)        =0.0_dp
     lapl_rhor(:)  =0.0_dp
     do ispin=1,nspin
       do jbf=1,basis%nbf
         do ibf=1,basis%nbf

           grad_rhor(:,ispin) = grad_rhor(:,ispin) + p_matrix(ibf,jbf,ispin) &
                *( basis_function_gradr(:,ibf) * basis_function_r(jbf) &
                 + basis_function_gradr(:,jbf) * basis_function_r(ibf) ) 

           tau(ispin)        = tau(ispin)        + p_matrix(ibf,jbf,ispin) &
                * DOT_PRODUCT( basis_function_gradr(:,ibf) , basis_function_gradr(:,jbf) )

           lapl_rhor(ispin)  = lapl_rhor(ispin)  + p_matrix(ibf,jbf,ispin) &
                              * (  SUM( basis_function_laplr(:,ibf) ) * basis_function_r(jbf)               &
                                 + basis_function_r(ibf) * SUM( basis_function_laplr(:,jbf) )               &
                                 + 2.0_dp * DOT_PRODUCT( basis_function_gradr(:,ibf) , basis_function_gradr(:,jbf) ) )

         enddo
       enddo
     enddo
     sigma(1) = SUM( grad_rhor(:,1)**2 )
     if(nspin==2) then
       sigma(2) = SUM( grad_rhor(:,1) * grad_rhor(:,2) )
       sigma(3) = SUM( grad_rhor(:,2)**2 )
     endif

   endif


   !
   ! LIBXC calls
   !
   dedd_r(:)    = 0.0_dp
   dedgd_r(:,:) = 0.0_dp

   do idft_xc=1,ndft_xc

     select case(xc_f90_info_family(xc_info(idft_xc)))

     case(XC_FAMILY_LDA)
       if( dft_xc_type(idft_xc) < 1000 ) then 
         call xc_f90_lda_exc_vxc(xc_func(idft_xc),1_C_INT,rhor_r(1),exc_libxc(1),vxc_libxc(1))
       else
         call my_lda_exc_vxc(nspin,dft_xc_type(idft_xc),rhor_r,exc_libxc(1),vxc_libxc)
!         call my_lda_exc_vxc_mu(1.00_dp,rhor_r,exc_libxc,vxc_libxc)
       endif

     case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
       if( dft_xc_type(idft_xc) < 2000 ) then 
         !
         ! Remove too small densities to stabilize the computation
         ! especially useful for Becke88
         if( ANY( rhor_r(:) > 1.0e-9_dp ) ) then
           call xc_f90_gga_exc_vxc(xc_func(idft_xc),1_C_INT,rhor_r(1),sigma(1),exc_libxc(1),vxc_libxc(1),vsigma(1))
         else
           exc_libxc(:)     = 0.0_dp
           vxc_libxc(:)     = 0.0_dp
           vsigma(:)        = 0.0_dp
         endif
       else
         call my_gga_exc_vxc_hjs(gamma_hybrid,rhor_r(1),sigma(1),exc_libxc(1),vxc_libxc(1),vsigma(1))
       endif

     case(XC_FAMILY_MGGA)
       call xc_f90_mgga_vxc(xc_func(idft_xc),1_C_INT,rhor_r(1),sigma(1),lapl_rhor(1),tau(1),vxc_libxc(1),vsigma(1),vlapl_rho(1),vtau(1))
       ! no expression for the energy
       exc_libxc(1) = 0.0_dp

     case default
       call die('functional is not LDA nor GGA nor hybrid nor meta-GGA')
     end select

     exc_xc = exc_xc + weight * exc_libxc(1) * SUM( rhor_r(:) ) * dft_xc_coef(idft_xc)

     dedd_r(:) = dedd_r(:) + vxc_libxc(:) * dft_xc_coef(idft_xc)

     !
     ! Set up divergence term if needed (GGA case)
     !
     if( xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_GGA &
        .OR. xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_HYB_GGA ) then
       if(nspin==1) then

         dedgd_r(:,1) = dedgd_r(:,1) + 2.0_dp * vsigma(1) * grad_rhor(:,1) * dft_xc_coef(idft_xc)

       else

         dedgd_r(:,1) = dedgd_r(:,1) + 2.0_dp * vsigma(1) * grad_rhor(:,1) * dft_xc_coef(idft_xc) &
                               + vsigma(2) * grad_rhor(:,2)

         dedgd_r(:,2) = dedgd_r(:,2) + 2.0_dp * vsigma(3) * grad_rhor(:,2) * dft_xc_coef(idft_xc) &
                               + vsigma(2) * grad_rhor(:,1)
       endif

     endif


     !
     ! In the case of the BJ06 meta-GGA functional, a spin-dependent shift is applied
     ! since the potential does not vanish at infinity
     !
     if( dft_xc_type(idft_xc) == XC_MGGA_X_BJ06 ) then
       dedd_r(:) = dedd_r(:) - SQRT( 5.0_dp * ABS(ehomo(:)) / 6.0_dp ) / pi
     endif

   enddo ! loop on the XC functional


   !
   ! Eventually set up the vxc term
   !
   if( .NOT. require_gradient ) then 
     ! LDA
     do ispin=1,nspin
       do jbf=1,basis%nbf
         ! Only the lower part is calculated
         do ibf=1,jbf ! basis%nbf 
           vxc_ij(ibf,jbf,ispin) =  vxc_ij(ibf,jbf,ispin) + weight &
               *  dedd_r(ispin) * basis_function_r(ibf) * basis_function_r(jbf) 
         enddo
       enddo
     enddo

   else 
     ! GGA
     do ispin=1,nspin
       do jbf=1,basis%nbf
         ! Only the lower part is calculated
         do ibf=1,jbf ! basis%nbf 
           vxc_ij(ibf,jbf,ispin) = vxc_ij(ibf,jbf,ispin) + weight  &
               * dedd_r(ispin) * basis_function_r(ibf) * basis_function_r(jbf) 

           vxc_ij(ibf,jbf,ispin) = vxc_ij(ibf,jbf,ispin) + weight &
                    * DOT_PRODUCT( dedgd_r(:,ispin) ,                                     &
                                      basis_function_gradr(:,ibf) * basis_function_r(jbf) &
                                    + basis_function_gradr(:,jbf) * basis_function_r(ibf) )
         enddo
       enddo
     enddo
   endif

 enddo ! loop on the grid point

 ! Symmetrize now
 do ispin=1,nspin
   do jbf=1,basis%nbf
     do ibf=1,jbf-1 ! basis%nbf 
       vxc_ij(jbf,ibf,ispin) = vxc_ij(ibf,jbf,ispin)
     enddo
   enddo
 enddo

 !
 ! Sum up the contributions from all procs only if needed
 if( parallel_grid ) then
   call xsum(normalization)
   call xsum(vxc_ij)
   call xsum(exc_xc)
 endif

 !
 ! Destroy operations
 do idft_xc=1,ndft_xc
   call xc_f90_func_end(xc_func(idft_xc))
 enddo

#else
!TODO write a call to teter to have MOLGW working without LIBXC
!   call teter_lda_vxc_exc(rhor,vxc,exc)
 write(stdout,*) 'XC energy and potential set to zero'
 write(stdout,*) 'LIBXC is not present'
#endif

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization(:)
 write(stdout,'(a,2x,f12.6,/)')    '  DFT xc energy (Ha):',exc_xc

 call stop_clock(timing_dft)

end subroutine dft_exc_vxc


!=========================================================================
subroutine dft_approximate_vhxc(basis,vhxc_ij)
 use m_definitions
 use m_mpi
 use m_timing
 use m_inputparam
 use m_basis_set
 use m_dft_grid
 use m_eri
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
#ifdef _OPENMP
 use omp_lib
#endif
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: vhxc_ij(basis%nbf,basis%nbf)
!=====
 integer              :: idft_xc
 integer              :: igrid,ibf,jbf,ispin
 real(dp)             :: rr(3)
 real(dp)             :: normalization
 real(dp)             :: weight
 real(dp)             :: basis_function_r(basis%nbf)
 real(dp)             :: rhor
 real(dp)             :: vxc,exc
 real(dp)             :: vsigma(2*nspin-1)
 real(dp)             :: vhartree
 real(dp)             :: vhgau(basis%nbf,basis%nbf)
 integer              :: iatom,igau,ngau
 real(dp),allocatable :: alpha(:),coeff(:)
!=====

 vhxc_ij(:,:) = 0.0_dp


 write(stdout,'(/,a)') ' Calculate approximate HXC potential with a superposition of atomic densities'

 do iatom=1,natom
   if( rank /= MODULO(iatom,nproc) ) cycle

   ngau = 4
   allocate(alpha(ngau),coeff(ngau))
   call element_atomicdensity(zatom(iatom),coeff,alpha)


   do igau=1,ngau
     call calculate_eri_approximate_hartree(.FALSE.,basis,basis%nbf,basis%nbf,x(:,iatom),alpha(igau),vhgau)
     vhxc_ij(:,:) = vhxc_ij(:,:) + vhgau(:,:) * coeff(igau) / 2.0_dp**1.25_dp / pi**0.75_dp * alpha(igau)**1.5_dp
   enddo

   deallocate(alpha,coeff)
 enddo

 write(stdout,*) 'Home-made functional LDA functional'
 !
 ! If it is the first time, set up the stored arrays
 !
 if( .NOT. ALLOCATED(bfr) ) call prepare_basis_functions_r(basis)

 normalization = 0.0_dp
 do igrid=1,ngrid

   rr(:) = rr_grid(:,igrid)
   weight = w_grid(igrid)

   !
   ! Get all the functions and gradients at point rr
   call get_basis_functions_r(basis,igrid,basis_function_r)

   !
   ! calculate the density at point r for spin up and spin down
   call setup_atomic_density(rr,rhor,vhartree)

   !
   ! Normalization
   normalization = normalization + rhor * weight

   call teter_lda_vxc_exc(rhor,vxc,exc)
   !
   ! HXC
   do jbf=1,basis%nbf
     do ibf=1,basis%nbf 
       vhxc_ij(ibf,jbf) =  vhxc_ij(ibf,jbf) + weight &
           *  vxc * basis_function_r(ibf) * basis_function_r(jbf)
     enddo
   enddo

 enddo ! loop on the grid point
 !
 ! Sum up the contributions from all procs only if needed
 call xsum(normalization)
 call xsum(vhxc_ij)

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization

end subroutine dft_approximate_vhxc


!=========================================================================
subroutine setup_atomic_density(rr,rhor,vhartree)
 use m_definitions
 use m_atoms
 use m_gaussian
 use m_inputparam
 implicit none

 real(dp),intent(in)  :: rr(3)
 real(dp),intent(out) :: rhor,vhartree
!=====
 real(dp),parameter   :: bondcharge=1.000_dp
 integer              :: iatom,igau,ngau,ibond
 real(dp)             :: dr
 real(dp),allocatable :: alpha(:),coeff(:)
 real(dp)             :: xbond(3)
!=====

 rhor = 0.0_dp
 vhartree = 0.0_dp
 do iatom=1,natom

   ngau = 4
   allocate(alpha(ngau),coeff(ngau))
   call element_atomicdensity(zatom(iatom),coeff,alpha)

   dr=NORM2( rr(:) - x(:,iatom) )

   do igau=1,ngau
     rhor     = rhor     + SQRT(alpha(igau)/pi)**3 * EXP( -alpha(igau)*dr**2) * coeff(igau)
     vhartree = vhartree + ERF(SQRT(alpha(igau))*dr)/dr * coeff(igau)
   enddo

   deallocate(alpha,coeff)
 enddo



end subroutine setup_atomic_density


!=========================================================================
subroutine calc_density_r(nspin,basis,p_matrix,rr,basis_function_r,rhor_r)
 use m_definitions
 use m_mpi
 use m_timing
 use m_basis_set
 use m_dft_grid,only: bf_rad2
 implicit none
 integer,intent(in)         :: nspin
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)  :: p_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in)  :: rr(3),basis_function_r(basis%nbf)
 real(dp),intent(out) :: rhor_r(nspin)
!=====
 integer :: ispin,ibf,jbf
!=====
 !
 ! Calculate the density rho at point r
 rhor_r(:)=0.0_dp
 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO REDUCTION(+:rhor_r) 
   do jbf=1,basis%nbf
     if( SUM( (basis%bf(jbf)%x0(:) - rr(:))**2 ) > bf_rad2(jbf) ) cycle
     ! implementing i <-> j symmetry does not save much time with ifort compiler
     do ibf=1,basis%nbf
       rhor_r(ispin)=rhor_r(ispin)+p_matrix(ibf,jbf,ispin)&
                         * basis_function_r(ibf) &
                         * basis_function_r(jbf)
     enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo


end subroutine calc_density_r


!=========================================================================
subroutine calc_density_gradr(nspin,nbf,p_matrix,basis_function_r,basis_function_gradr,grad_rhor)
 use m_definitions
 use m_mpi
 use m_timing
 implicit none
 integer,intent(in)   :: nspin,nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(in)  :: basis_function_r(nbf)
 real(dp),intent(in)  :: basis_function_gradr(3,nbf)
 real(dp),intent(out) :: grad_rhor(3,nspin)
!=====
 integer :: ispin,ibf,jbf
!=====

 grad_rhor(:,:)=0.0_dp
 do ispin=1,nspin
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO REDUCTION(+:grad_rhor) 
   do jbf=1,nbf
     ! implementing i <-> j symmetry does not save much time with ifort
     ! compiler
     do ibf=1,nbf
       grad_rhor(:,ispin) = grad_rhor(:,ispin) + p_matrix(ibf,jbf,ispin) &
            *( basis_function_gradr(:,ibf) * basis_function_r(jbf) &
             + basis_function_gradr(:,jbf) * basis_function_r(ibf) )
     enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
 enddo

end subroutine calc_density_gradr


!=========================================================================
subroutine teter_lda_vxc_exc(rhor,vxc,exc)
 use m_definitions
 implicit none

 real(dp),intent(in) :: rhor
 real(dp),intent(out) :: vxc,exc
!=====
 !
 ! The usual full LDA parameters of Teter
 real(dp),parameter :: a0p=.4581652932831429_dp
 real(dp),parameter :: a1p=2.217058676663745_dp
 real(dp),parameter :: a2p=0.7405551735357053_dp
 real(dp),parameter :: a3p=0.01968227878617998_dp
 real(dp),parameter :: b1p=1.0_dp
 real(dp),parameter :: b2p=4.504130959426697_dp
 real(dp),parameter :: b3p=1.110667363742916_dp
 real(dp),parameter :: b4p=0.02359291751427506_dp
 real(dp)           :: d1m1
 real(dp)           :: dd1drs,dn1drs,dexcdrs
 real(dp)           :: n1,d1
 real(dp)           :: rs
!=====

 rs = ( 3.0_dp / (4.0_dp*pi*rhor) )**(1.0_dp/3.0_dp) 
 n1 = a0p + rs * (a1p + rs * ( a2p + rs * a3p ) )
 d1 = rs * ( b1p + rs * ( b2p + rs * ( b3p + rs * b4p ) ) )
 d1m1 = 1.0_dp / d1

 ! Firstly, exchange-correlation energy
 exc = -n1 * d1m1

 ! Secondly, exchange-correlation potential
 dn1drs = a1p + rs * ( 2.0_dp * a2p + rs * ( 3.0_dp * a3p ) )
 dd1drs = b1p + rs * ( 2.0_dp * b2p + rs * ( 3.0_dp * b3p + rs * ( 4.0_dp * b4p ) ) )

 ! dexcdrs is d(exc)/d(rs)
 dexcdrs = -( dn1drs + exc * dd1drs ) * d1m1
 vxc = exc - rs * dexcdrs / 3.0_dp


end subroutine teter_lda_vxc_exc


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

! *************************************************************************

 select case(ixc)
 case(1200)
   !
   ! The Slater exchange
   a0p=3.0_dp/8.0_dp*(18.0_dp/pi**2)**(1.0_dp/3.0_dp)
   a1p=0.0_dp
   a2p=0.0_dp
   a3p=0.0_dp
   b1p=1.0_dp
   b2p=0.0_dp
   b3p=0.0_dp
   b4p=0.0_dp
   da0=0.0_dp
   da1=0.0_dp
   da2=0.0_dp
   da3=0.0_dp
   db1=0.0_dp
   db2=0.0_dp
   db3=0.0_dp
   db4=0.0_dp
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
subroutine my_lda_exc_vxc_mu(mu,rhor,exc,vxc)
 use m_definitions
 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: mu
 integer,parameter :: npt=1
 integer,parameter :: order=1
!arrays
 real(dp),intent(in) :: rhor(npt)
 real(dp),intent(out) :: exc(npt),vxc(npt)

!Local variables-------------------------------
!Set value of alpha in "X-alpha" method
!scalars
 integer :: ipt
 real(dp) :: efac,rs,rsm1,vfac
 character(len=500) :: message

 real(dp)           :: biga,kf,fact_mu
 real(dp)           :: rs_step=1.0e-6_dp
 real(dp)           :: rsp,rsm1p
 real(dp)           :: bigap,kfp,fact_mup
 real(dp)           :: omega
 real(dp)           :: rho,aa,f_aa
 real(dp)           :: exp4aa2,rhodaadrho,dfdaa

! *************************************************************************

 efac=0.75_dp*(1.5_dp/pi)**(2.0_dp/3.0_dp)
 vfac=4.0/3.0 * efac

 omega = mu
 rs= (3.0_dp/(4.0_dp*pi*rhor(1)))**(1.0_dp/3.0_dp)
 rho = 3.0 / ( 4.0 * pi * rs**3 )
 kf  = ( 3.0 * pi**2 * rho )**(1.0/3.0)
 aa  = omega / ( 2.0 * kf )

 exp4aa2 = EXP(-1.0/(4.0*aa**2))

 f_aa = 8.0/3.0 * aa &
       * ( SQRT(pi) * ERF(0.5/aa) &
          + (2.0 * aa - 4.0 * aa**3) * exp4aa2 &
          - 3.0 * aa + 4.0 * aa**3 )

 rhodaadrho = - omega/ ( 6.0 * kf )

 dfdaa = f_aa / aa + 8.0 * aa * ( 4.0 * aa**2 * ( 1.0 - exp4aa2 ) - 1.0 )

 exc(1) = -efac/rs * (1.0 - f_aa)

 vxc(1) = -vfac/rs * (1.0 - f_aa) + efac/rs * dfdaa * rhodaadrho


end subroutine my_lda_exc_vxc_mu


!=========================================================================
subroutine my_gga_exc_vxc_hjs(omega,nn,sigma,exc,vxc,vsigma)
 use m_definitions
 implicit none
!=====
 real(dp),intent(in)  :: omega,nn,sigma
 real(dp),intent(out) :: exc,vxc,vsigma
!=====
 real(dp),parameter :: ss0=2.0
 ! HJS parameters
 real(dp),parameter :: aabar= 0.757211
 real(dp),parameter :: bb   =-0.106364
 real(dp),parameter :: cc   =-0.118649
 real(dp),parameter :: dd   = 0.609650
 real(dp),parameter :: ee   =-0.0477963
 ! PBE parameters
 real(dp),parameter :: a2= 0.0159941
 real(dp),parameter :: a3= 0.0852995
 real(dp),parameter :: a4=-0.160368
 real(dp),parameter :: a5= 0.152645
 real(dp),parameter :: a6=-0.0971263
 real(dp),parameter :: a7= 0.0422061
 real(dp),parameter :: b1= 5.33319
 real(dp),parameter :: b2=-12.4780
 real(dp),parameter :: b3=11.0988
 real(dp),parameter :: b4=-5.11013
 real(dp),parameter :: b5= 1.71468
 real(dp),parameter :: b6=-0.610380
 real(dp),parameter :: b7= 0.307555
 real(dp),parameter :: b8=-0.0770547
 real(dp),parameter :: b9= 0.0334840
!=====
 real(dp) :: efac
 real(dp) :: nn_local,sigma_local,rs
 real(dp) :: ss
 real(dp) :: kf
 real(dp) :: nu
 real(dp) :: chi
 real(dp) :: lambda,eta,zeta
 real(dp) :: hh_s
 real(dp) :: ffbar_s
 real(dp) :: ggbar_s
 real(dp) :: factor_w
 real(dp) :: exc_nn,exc_sigma
 real(dp) :: fx,dfxds,dfxdnu
 real(dp) :: dsdsigma,dsdn,dnudn
!=====

 efac=0.75_dp * (1.5_dp/pi)**(2.0_dp/3.0_dp)


!HOME MADE    !
!HOME MADE    ! first calculation
!HOME MADE    nn_local = nn
!HOME MADE    sigma_local = sigma
!HOME MADE    
!HOME MADE    rs = ( 3.0 / (4.0 *pi * nn_local) )**(1./3.)
!HOME MADE    kf = (9.0_dp * pi / 4.0_dp)**(1.0_dp/3.0_dp) / rs
!HOME MADE    nu = omega / kf
!HOME MADE    ss = SQRT(sigma_local) / ( 2.0_dp * kf * nn_local )
!HOME MADE   
!HOME MADE    hh_s = ( a2*ss**2 + a3*ss**3 + a4*ss**4 + a5*ss**5 + a6*ss**6 + a7*ss**7 ) &
!HOME MADE         / ( 1.0_dp + b1*ss + b2*ss**2 + b3*ss**3 + b4*ss**4 + b5*ss**5 + b6*ss**6 + b7*ss**7 + b8*ss**8 + b9*ss**9 )
!HOME MADE   
!HOME MADE    ffbar_s = 1.0_dp - ss**2 / ( 27.0_dp * cc * (1.0_dp + ss**2/ss0**2 ) ) &
!HOME MADE               - ss**2 * hh_s / (2.0_dp * cc )
!HOME MADE   
!HOME MADE   
!HOME MADE    zeta   = ss**2 * hh_s
!HOME MADE    eta    = aabar + ss**2 * hh_s
!HOME MADE    lambda = dd    + ss**2 * hh_s
!HOME MADE    chi = nu / SQRT( lambda + nu**2)
!HOME MADE   
!HOME MADE    ggbar_s = -2./5.*cc*ffbar_s*lambda -4./15.*bb*lambda**2 - 6./5.*aabar*lambda**3 &
!HOME MADE             -4./5.*SQRT(pi)*lambda**(7./2.) &
!HOME MADE             -12./5.*lambda**(7./2.) * ( SQRT(zeta)-SQRT(eta) )
!HOME MADE    ggbar_s = ggbar_s / ee
!HOME MADE   
!HOME MADE   
!HOME MADE    factor_w = aabar - 4./9.*bb/lambda*(1.0-chi) - 4./9.*cc*ffbar_s/lambda**2 * (1.0 - 1.5*chi+0.5*chi**3)  &
!HOME MADE              -8./9.*ee*ggbar_s/lambda**3 * ( 1.0 - 15./8.*chi + 5./4.*chi**3 -3./8.*chi**5 ) &
!HOME MADE              + 2.*nu   * ( SQRT(zeta+nu**2)- SQRT(eta+nu**2) ) &
!HOME MADE              + 2.*zeta * LOG( ( nu + SQRT(zeta+nu**2) ) / ( nu + SQRT(lambda + nu**2) ) ) &
!HOME MADE              - 2.*eta  * LOG( ( nu + SQRT( eta+nu**2) ) / ( nu + SQRT(lambda + nu**2) ) ) 
!HOME MADE   
!HOME MADE   
!HOME MADE    exc = -efac/rs * factor_w
!HOME MADE
!HOME MADE write(stdout,*) 'exc1=',exc

 !
 ! call to the nwchem subroutine
 !
 rs = ( 3.0 / (4.0*pi*nn) )**(1.0/3.0)
 kf = (9.0_dp * pi / 4.0_dp)**(1.0_dp/3.0_dp) / rs
 ss = SQRT(sigma) / ( 2.0_dp * kf * nn )

 call HSE08Fx(omega,1,nn,ss,fx,dfxds,dfxdnu)

 exc = -efac/rs*fx

 dsdsigma= 1.0_dp / ( 4.0_dp * kf * nn * SQRT(sigma) )

 vsigma = -efac/rs * nn * dfxds * dsdsigma

 dsdn  = SQRT(sigma) / (2.0_dp * (3.0*pi**2)**(1./3.) ) * (-4.0/3.0) * nn**(-7.0/3.0)
 dnudn = omega / (3.0*pi**2)**(1./3.) * (-1.0/3.0) * nn**(-4.0/3.0)
 vxc = -efac/rs * nn * ( dfxds * dsdn + dfxdnu * dnudn ) - (4.0/3.0)*efac/rs *fx


end subroutine my_gga_exc_vxc_hjs


!=========================================================================
subroutine HSE08Fx(omega,ipol,rho,s,Fxhse,d10Fxhse,d01Fxhse)

 implicit none
!
!case...start
! #include "case.fh"
!case...end
!
! HSE evaluates the Heyd et al. Screened Coulomb
! Exchange Functional
!
! Calculates the enhancement factor
!
 double precision omega
 integer ipol
 double precision rho,s,Fxhse,d10Fxhse,d01Fxhse

 double precision  A,B,C,D,E
 double precision  ha2,ha3,ha4,ha5,ha6,ha7
 double precision  hb1,hb2,hb3,hb4,hb5,hb6,hb7,hb8,hb9
 double precision  smax,strans,sconst

 double precision  zero,one,two,three,four,five,six,seven,eight
 double precision  nine,ten
 double precision  fifteen,sixteen

 double precision  H,hnum,hden 
 double precision  d1H,d1hnum,d1hden 
 double precision  s2,s3,s4,s5,s6,s7,s8,s9
 double precision  Fs, d1Fs
 double precision  zeta, lambda, eta, kf, nu, chi, lambda2
 double precision  d1zeta,d1lambda,d1eta,d1nu,d1chi,d1lambda2
 double precision  EGs,d1EGs
 double precision  nu2,L2,L3,nu3,nu4,nu5,nu6
 double precision  Js,Ks,Ms,Ns
 double precision  d1Js,d1Ks,d1Ms,d1Ns

 double precision  tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
 double precision  tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15
 double precision  Fxhse1,Fxhse2,Fxhse3,Fxhse4,Fxhse5,Fxhse6
 double precision  d1Fxhse1,d1Fxhse2,d1Fxhse3,d1Fxhse4,d1Fxhse5
 double precision  d1Fxhse6,d1Fxhse7

 double precision  r42,r27,r12,r15,r14,r18,r20,r30,r56,r72
 double precision  r16,r32,r24,r48,r11,r64,r35
 double precision  pi,pi2,srpi,s02
 double precision  f12,f13,f32,f52,f72,f92

!
!Constants for HJS hole
!
 Data A,B,C,D,E  &
     / 7.57211D-1,-1.06364D-1,-1.18649D-1, &
       6.09650D-1,-4.77963D-2 /
!
!Constants for fit of H(s) (PBE hole)
!Taken from JCTC_5_754 (2009)
!
 Data ha2,ha3,ha4,ha5,ha6,ha7 &
     / 1.59941D-2,8.52995D-2,-1.60368D-1,1.52645D-1, &
      -9.71263D-2,4.22061D-2 /

 Data hb1,hb2,hb3,hb4,hb5,hb6,hb7,hb8,hb9 &
      / 5.33319D0,-12.4780D0,11.0988D0,-5.11013D0,&
       1.71468D0,-6.10380D-1,3.07555D-1,-7.70547D-2,&
       3.34840D-2 /

!
!Whole numbers used during evaluation
!
 Data zero,one,two,three,four,five,six,seven,eight,nine,ten &
      / 0D0,1D0,2D0,3D0,4D0,5D0,6D0,7D0,8D0,9D0,10D0 /
  
 Data r11,r12,r14,r15,r16,r18,r20,r24,r27,r30,r32 &
      / 11D0,12D0,14D0,15D0,16D0,18D0,20D0,24D0,27d0,30D0,32D0 /

 Data r35,r42,r48,r56,r64,r72 &
      / 35D0,42D0,48D0,56D0,64D0,72D0 /
!
!Fractions used during evaluation
!
 Data f12     / 0.5D0 /
!
!General constants
!
 f13   = one/three
 f32   = three/two
 f52   = five/two
 f72   = seven/two
 f92   = nine/two
 pi    = ACos(-one)
 pi2   = pi*pi
 srpi = dsqrt(pi)
!
!
!Calculate prelim variables
!
 s2 = s*s
 s02 = s2/four
 s3 = s2*s
 s4 = s3*s
 s5 = s4*s
 s6 = s5*s
 s7 = s6*s
 s8 = s7*s
 s9 = s8*s

!
!Calculate H(s) the model exhange hole
!
 hnum = ha2*s2 + ha3*s3 + ha4*s4 + ha5*s5 + ha6*s6 + ha7*s7 
 hden = one + hb1*s + hb2*s2 + hb3*s3 + hb4*s4 + hb5*s5 + &
        hb6*s6 + hb7*s7 + hb8*s8 + hb9*s9
 H = hnum/hden

!
!Calculate helper variables
!
 zeta = s2*H
 eta = A + zeta
 lambda = D + zeta
 if (ipol.eq.1) then
    kf = (three*pi2*rho)**f13 
 else
    kf = (six*pi2*rho)**f13 
 endif
 nu = omega/kf
 chi = nu/dsqrt(lambda+nu**two)
 lambda2 = (one+chi)*(lambda+nu**two)

!
!Calculate F(H(s)) for the model exhange hole
!
 Fs = one-s2/(r27*C*(one+s02))-zeta/(two*C)

!
!Calculate EG(s) 
!
 EGs = -(two/five)*C*Fs*lambda - (four/r15)*B*lambda**two - &
       (six/five)*A*lambda**three - &
       (four/five)*srpi*lambda**(seven/two) -&
       (r12/five)*(lambda**(seven/two))*(dsqrt(zeta)-dsqrt(eta))
 
!
!Calculate the denominators needed
!

 nu2 = nu*nu
 Js = (dsqrt(zeta+nu2)+dsqrt(eta+nu2))*(dsqrt(zeta+nu2)+nu) 
 Ks = (dsqrt(zeta+nu2)+dsqrt(eta+nu2))*(dsqrt(eta+nu2)+nu) 
 Ms = (dsqrt(zeta+nu2)+dsqrt(lambda+nu2))*(dsqrt(lambda+nu2)+nu) 
 Ns = (dsqrt(eta+nu2)+dsqrt(lambda+nu2))*(dsqrt(lambda+nu2)+nu) 

!
!  The final value for the enhancement factor is
!
 tmp1 = one + f12*chi
 tmp2 = one + (nine/eight)*chi + (three/eight)*chi**two 
 Fxhse1  = A*(zeta/Js + eta/Ks) 
 Fxhse2  = -(four/nine)*B/lambda2
 Fxhse3  = -(four/nine)*C*Fs*tmp1/lambda2**two
 Fxhse4  = -(eight/nine)*EGs*tmp2/lambda2**three
 Fxhse5  = two*zeta*dlog(one -D/Ms)
 Fxhse6  = -two*eta*dlog(one -(D-A)/Ns)

 Fxhse = Fxhse1+Fxhse2+Fxhse3+Fxhse4+Fxhse5+Fxhse6
!
!Calculate the first derivative of H with respect to the
!reduced density gradient, s.
!
 d1hnum = two*ha2*s + three*ha3*s2 + four*ha4*s3 + &
           five*ha5*s4 + six*ha6*s5 + seven*ha7*s6

 d1hden  = hb1 + two*hb2*s +three*hb3*s2 + four*hb4*s3 + &
           five*hb5*s4 + six*hb6*s5 + seven*hb7*s6 +&
           eight*hb8*s7 + nine*hb9*s8 
 d1H =   (hden*d1hnum -hnum*d1hden)/hden**two

!
!calculate first derivative of variables needed with respect to s
!
 d1zeta = two*s*H + s2*d1H
 d1eta  = d1zeta
 d1lambda = d1zeta
 d1chi = -f12*nu*d1zeta/(lambda + nu2)**f32
 d1lambda2 = d1chi*(lambda + nu**two) + (one+chi)*d1lambda
 !d1lambda2 = (d1lambda*(one-chi)+lambda*d1chi)/(one-chi)**two

!
!calculate the first derivative of Fs with respect to s
!
 d1Fs = -two*s/(r27*C*(one+s02)**two) - d1zeta/(two*C)

!
!Calculate the first derivate of EGs with respect to s
!
 d1EGs = -(two/five)*C*(d1Fs*lambda + Fs*d1lambda) -&
         (eight/r15)*B*lambda*d1lambda -&
         (r18/five)*A*lambda*lambda*d1lambda -&
         (r14/five)*srpi*d1lambda*lambda**f52 -&
         (r42/five)*(lambda**f52)*&
         d1lambda*(dsqrt(zeta)-dsqrt(eta))-&
         (six/five)*(lambda**(seven/two))*&
         (d1zeta/dsqrt(zeta)-d1eta/dsqrt(eta))

!
!Calculate the first derivate of denominators needed with respect
!to s
!
 tmp1 = (dsqrt(zeta+nu2)+nu)/(dsqrt(eta+nu2)) 
 tmp2 = (dsqrt(eta+nu2)+nu)/(dsqrt(zeta+nu2))

 d1Js = f12*d1zeta*(two+tmp1+tmp2)
 d1Ks = d1Js

 tmp3 = (dsqrt(zeta+nu2)+nu)/(dsqrt(lambda+nu2))
 tmp4 = (dsqrt(lambda+nu2)+nu)/(dsqrt(zeta+nu2)) 
 d1Ms = f12*d1zeta*(two +tmp3+tmp4)

 tmp5 = (dsqrt(lambda+nu2)+nu)/(dsqrt(eta+nu2))
 tmp6 = (dsqrt(eta+nu2)+nu)/(dsqrt(lambda+nu2))
 d1Ns = f12*d1zeta*(two + tmp5+tmp6)
!
!Calculate the derivative of the 08-Fxhse with respect to s
!
 L2 = lambda2*lambda2
 L3 = lambda2*lambda2*lambda2
 d1Fxhse1  = A*( (Js*d1zeta - zeta*d1Js)/(Js*Js) +&
                 (Ks*d1zeta - eta*d1Ks)/(Ks*Ks) ) 

 d1Fxhse2  = (four/nine)*B*d1lambda2/L2 

 tmp9 = d1lambda2/lambda2
 tmp7 = d1Fs - two*Fs*tmp9
 tmp8 = one + f12*chi
 tmp10 =  f12*Fs*d1chi

 d1Fxhse3 = -(four*C/(nine*L2))*(tmp7*tmp8+tmp10)


   tmp7 = one + (nine/eight)*chi+(three/eight)*chi*chi
   tmp8 = (nine/eight)*d1chi + (six/eight)*chi*d1chi

  d1Fxhse4 = -(eight/(nine*L3))*((d1EGs-three*EGs*tmp9)*tmp7 &
            + EGs*tmp8)
 d1Fxhse5  = two*d1zeta*dlog(one-D/Ms) + &
            two*zeta*D*d1Ms/(Ms*Ms*(one-D/Ms)) 

 d1Fxhse6  = -two*d1eta*dlog(one- (D-A)/Ns) - &
            two*eta*(D-A)*d1Ns/(Ns*Ns*(one-(D-A)/Ns)) 

 d10Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6
!
!Calculate the derivative of 08-Fxhse with respect to nu
!
 nu3 = nu2*nu

 d1Fxhse1 = -((A*(nu + dsqrt(eta + nu2))*zeta)/ &
             (dsqrt(eta + nu2)*dsqrt(nu2 + zeta)*&
             (nu + dsqrt(nu2 + zeta))*&
             (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))))

 d1Fxhse2 = -((A*eta*(nu/dsqrt(eta + nu2) + nu/ &
             dsqrt(nu2 + zeta)))/&
             ((nu + dsqrt(eta + nu2))*&
             (dsqrt(eta + nu2) + dsqrt(nu2 + zeta))**two)) -&
             (A*eta*(one + nu/dsqrt(eta + nu2)))/&
             ((nu + dsqrt(eta + nu2))**two*&
             (dsqrt(eta + nu2) + dsqrt(nu2 + zeta)))

 d1Fxhse3 = (four*B)/(nine*(lambda + nu2)**(f32))

 d1Fxhse4 = (two*C*Fs)/(three*(lambda + nu2)**(f52))

 d1Fxhse5 = (five*EGs*(lambda**two + four*nu3* &
             (nu + dsqrt(lambda + nu2)) +&
             lambda*nu*(five*nu + three*dsqrt(lambda + nu2))))/&
    (three*(lambda + nu2)**four*(nu + dsqrt(lambda + nu2))**three)

 d1Fxhse6 = (two*D*zeta*(nu + dsqrt(nu2 + zeta)))/&
             (dsqrt(lambda + nu2)*dsqrt(nu2 + zeta)*&
             (-D + lambda + (nu + dsqrt(lambda + nu2))*&
             (nu + dsqrt(nu2 + zeta))))

 d1Fxhse7 = (two*(A - D)*eta*(nu + dsqrt(eta + nu2)))/ &
             (dsqrt(eta + nu2)*dsqrt(lambda + nu2)*&
             (A - D + lambda + nu2 + nu*dsqrt(eta + nu2) +&
             nu*dsqrt(lambda + nu2) +&
             dsqrt(eta + nu2)*dsqrt(lambda + nu2)))


 d01Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6+d1Fxhse7
 
end subroutine HSE08Fx

!=========================================================================
