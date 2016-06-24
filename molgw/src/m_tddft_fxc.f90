!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods and data to evaluate the TDDFT kernel f_xc
! and in particular integrals ( i a | f_xc | j b )
!
!=========================================================================
module m_tddft_fxc
 use m_definitions
 use m_warning
 use m_mpi
 use m_memory
 use m_timing
 use m_inputparam
 use m_dft_grid
 use m_basis_set

 integer,private   :: nspin_tddft
 logical,private   :: require_gradient

 !
 ! fxc kernel
 ! Second derivatives with respect to rho, sigma 
 real(dp),allocatable,protected :: v2rho2(:,:)
 real(dp),allocatable,protected :: vsigma(:,:)
 real(dp),allocatable,protected :: v2rhosigma(:,:)
 real(dp),allocatable,protected :: v2sigma2(:,:)

 ! 
 ! Wfn and gradients evaluated on the grid points
 real(dp),allocatable,protected :: wf_r(:,:,:)
 real(dp),allocatable,protected :: wf_gradr(:,:,:,:)
 real(dp),allocatable,protected :: rho_gradr(:,:,:)

 ! Intermediate quantities
 real(dp),allocatable,protected :: grad_ij(:,:,:)
 real(dp),allocatable,protected :: grad_kl(:,:,:)
 real(dp),allocatable,protected :: dot_ij_kl(:,:)
 real(dp),allocatable,protected :: dot_rho_ij(:,:)
 real(dp),allocatable,protected :: dot_rho_kl(:,:)


contains


!=========================================================================
subroutine prepare_tddft(nstate,basis,c_matrix,occupation)
#ifdef HAVE_LIBXC
 use libxc_funcs_m
 use xc_f90_lib_m
 use xc_f90_types_m
#endif
 implicit none

 integer,intent(in)               :: nstate
 type(basis_set),intent(in)       :: basis
 real(dp),intent(in)              :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)              :: occupation(nstate,nspin)
!=====
#ifdef HAVE_LIBXC
 type(xc_f90_pointer_t) :: xc_func(ndft_xc),xc_functest
 type(xc_f90_pointer_t) :: xc_info(ndft_xc),xc_infotest
#endif
 real(dp),parameter   :: kernel_capping=1.0e14_dp
 character(len=256)   :: string
 integer              :: idft_xc,igrid
 integer              :: ispin,ibf,jbf
 real(dp)             :: rr(3)
 real(dp)             :: basis_function_r(basis%nbf)
 real(dp)             :: basis_function_gradr(3,basis%nbf)
 real(dp)             :: rhor_r(nspin)
 real(dp)             :: grad_rhor(3,nspin)
 real(dp)             :: p_matrix(basis%nbf,basis%nbf,nspin)
 real(dp)             :: max_v2sigma2
 real(dp),allocatable :: rho_c(:)
 real(dp),allocatable :: v2rho2_c(:)
 real(dp),allocatable :: sigma_c(:)
 real(dp),allocatable :: vrho_c(:)
 real(dp),allocatable :: vsigma_c(:)
 real(dp),allocatable :: v2rhosigma_c(:)
 real(dp),allocatable :: v2sigma2_c(:)
!=====

 if( is_triplet ) then
   nspin_tddft = 2
 else
   nspin_tddft = nspin
 endif

 call init_dft_grid(grid_level)

 allocate( rho_c(nspin_tddft)            )
 allocate( v2rho2_c(2*nspin_tddft-1)     )
 allocate( sigma_c(2*nspin_tddft-1)      )
 allocate( vrho_c(nspin_tddft)           )
 allocate( vsigma_c(2*nspin_tddft-1)     )
 allocate( v2rhosigma_c(5*nspin_tddft-4) )
 allocate( v2sigma2_c(5*nspin_tddft-4)   )
 !
 ! Prepare DFT kernel calculation with Libxc
 !
 require_gradient  =.FALSE.
 do idft_xc=1,ndft_xc
   if(nspin_tddft==1) then
     call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), dft_xc_type(idft_xc), XC_UNPOLARIZED)
   else
     call xc_f90_func_init(xc_func(idft_xc), xc_info(idft_xc), dft_xc_type(idft_xc), XC_POLARIZED)
   endif
   call xc_f90_info_name(xc_info(idft_xc),string)
   write(stdout,'(a,i4,a,i6,5x,a)') '   XC functional ',idft_xc,' :  ',xc_f90_info_number(xc_info(idft_xc)),&
         TRIM(string)
   if( MODULO(xc_f90_info_flags( xc_info(idft_xc)),XC_FLAGS_HAVE_FXC*2) < XC_FLAGS_HAVE_FXC ) then
     call die('This functional does not have the kernel implemented in Libxc')
   endif
   if(xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_GGA     ) require_gradient  =.TRUE.
   if(xc_f90_info_family(xc_info(idft_xc)) == XC_FAMILY_HYB_GGA ) require_gradient  =.TRUE.
 enddo

 !
 ! calculate rho, grad rho and the kernel
 ! 
 ! Setup the density matrix P from C
! call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)
 do ispin=1,nspin
   do jbf=1,basis%nbf
     do ibf=1,basis%nbf
       p_matrix(ibf,jbf,ispin) = SUM( occupation(:,ispin) * c_matrix(ibf,:,ispin) * c_matrix(jbf,:,ispin) )
     enddo
   enddo
 enddo



 allocate(v2rho2(ngrid,2*nspin_tddft-1),wf_r(ngrid,basis%nbf,nspin))
 v2rho2(:,:) = 0.0_dp

 if(require_gradient) then
   allocate(vsigma(ngrid,2*nspin_tddft-1))
   allocate(v2rhosigma(ngrid,5*nspin_tddft-4))
   allocate(v2sigma2(ngrid,5*nspin_tddft-4))
   allocate(wf_gradr(3,ngrid,basis%nbf,nspin))
   allocate(rho_gradr(3,ngrid,nspin))
   vsigma(:,:)     = 0.0_dp
   v2rhosigma(:,:) = 0.0_dp
   v2sigma2(:,:)   = 0.0_dp
 endif


 max_v2sigma2 = -1.0_dp
 do igrid=1,ngrid

   rr(:) = rr_grid(:,igrid)

   if( .NOT. ALLOCATED(bfr) ) call prepare_basis_functions_r(basis)
   if( require_gradient .AND. .NOT. ALLOCATED(bfgr) ) call prepare_basis_functions_gradr(basis)
   !
   ! Get all the functions and gradients at point rr
   call get_basis_functions_r(basis,igrid,basis_function_r)
   !
   ! store the wavefunction in r
   do ispin=1,nspin
     wf_r(igrid,:,ispin) = MATMUL( basis_function_r(:) , c_matrix(:,:,ispin) )
   enddo

   if( require_gradient ) then
     call get_basis_functions_gradr(basis,igrid,basis_function_gradr)
     !
     ! store the wavefunction in r
     do ispin=1,nspin
       wf_gradr(:,igrid,:,ispin) = MATMUL( basis_function_gradr(:,:) , c_matrix(:,:,ispin) )
     enddo
   endif


   call calc_density_pmatrix(nspin,basis,p_matrix,rr,basis_function_r,rhor_r)
   if( require_gradient ) then
     call calc_density_gradr_pmatrix(nspin,basis%nbf,p_matrix,basis_function_r,basis_function_gradr,grad_rhor)

     rho_gradr(:,igrid,:) = grad_rhor(:,:)

     if( nspin_tddft==1 ) then
       sigma_c(1) = SUM( grad_rhor(:,1)**2 )
     else if( nspin==2 ) then
       sigma_c(2) = SUM( grad_rhor(:,1) * grad_rhor(:,2) )
       sigma_c(3) = SUM( grad_rhor(:,2)**2 )
     else ! triplet excitations from singlet ground-state
       sigma_c(:) = SUM( grad_rhor(:,1)**2 ) * 0.25_dp
     endif

   endif


   if( nspin_tddft==1 ) then
     rho_c(1) = rhor_r(1)
   else if( nspin==2 ) then
     rho_c(:) = rhor_r(:)
   else ! triplet excitations from singlet ground-state
     rho_c(:) = rhor_r(1) * 0.5_dp
   endif

   !
   ! Calculate the kernel
   ! 
   do idft_xc=1,ndft_xc
     select case(xc_f90_info_family(xc_info(idft_xc)))
     case(XC_FAMILY_LDA)
       call xc_f90_lda_fxc(xc_func(idft_xc),1,rho_c(1),v2rho2_c(1))
     case(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
       call xc_f90_gga_vxc(xc_func(idft_xc),1,rho_c(1),sigma_c(1),vrho_c(1),vsigma_c(1))
       call xc_f90_gga_fxc(xc_func(idft_xc),1,rho_c(1),sigma_c(1),v2rho2_c(1),v2rhosigma_c(1),v2sigma2_c(1))
     case default
       call die('Other kernels not yet implemented')
     end select
     !
     ! Remove the too large values for stability
     v2rho2_c(:) = MIN( v2rho2_c(:), kernel_capping )
     if(require_gradient) then
       max_v2sigma2 = MAX(ABS(v2sigma2_c(1)),max_v2sigma2)
       vsigma_c(:)     = MIN( vsigma_c(:), kernel_capping )
       v2rhosigma_c(:) = MIN( v2rhosigma_c(:), kernel_capping )
       v2sigma2_c(:)   = MIN( v2sigma2_c(:), kernel_capping )
     endif

     ! Store the result with the weight
     ! Remove too large values for stability
     v2rho2(igrid,:)     = v2rho2(igrid,:) + v2rho2_c(:) * w_grid(igrid) * dft_xc_coef(idft_xc)
     if(require_gradient) then
       vsigma(igrid,:)     = vsigma(igrid,:)     + vsigma_c(:)     * w_grid(igrid) * dft_xc_coef(idft_xc)
       v2rhosigma(igrid,:) = v2rhosigma(igrid,:) + v2rhosigma_c(:) * w_grid(igrid) * dft_xc_coef(idft_xc)
       v2sigma2(igrid,:)   = v2sigma2(igrid,:)   + v2sigma2_c(:)   * w_grid(igrid) * dft_xc_coef(idft_xc)
     endif

   enddo
 enddo
 if(require_gradient) then
   call xmax(max_v2sigma2)
   write(stdout,'(a,e18.6)') ' Maximum numerical value for fxc: ',max_v2sigma2

   allocate(grad_ij(3,ngrid,nspin))
   allocate(grad_kl(3,ngrid,nspin))
   allocate(dot_ij_kl(ngrid,nspin))
   allocate(dot_rho_ij(ngrid,nspin))
   allocate(dot_rho_kl(ngrid,nspin))

 endif

 deallocate( rho_c )
 deallocate( v2rho2_c )
 deallocate( sigma_c )
 deallocate( vrho_c )
 deallocate( vsigma_c )
 deallocate( v2rhosigma_c )
 deallocate( v2sigma2_c )

end subroutine prepare_tddft


!=========================================================================
function eval_fxc_rks_singlet(istate,jstate,ijspin,kstate,lstate,klspin)
 implicit none
!=====
 real(dp)           :: eval_fxc_rks_singlet
 integer,intent(in) :: istate,jstate,kstate,lstate
 integer,intent(in) :: ijspin,klspin
!=====
 
 eval_fxc_rks_singlet = SUM(  wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) &
                    * wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) &
                    * v2rho2(:,ijspin) * 2.0_dp )

 if(require_gradient) then

   grad_ij(1,:,ijspin) = wf_gradr(1,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(1,:,jstate,ijspin) * wf_r(:,istate,ijspin)
   grad_ij(2,:,ijspin) = wf_gradr(2,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(2,:,jstate,ijspin) * wf_r(:,istate,ijspin)
   grad_ij(3,:,ijspin) = wf_gradr(3,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(3,:,jstate,ijspin) * wf_r(:,istate,ijspin)
   grad_kl(1,:,klspin) = wf_gradr(1,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(1,:,lstate,klspin) * wf_r(:,kstate,klspin)
   grad_kl(2,:,klspin) = wf_gradr(2,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(2,:,lstate,klspin) * wf_r(:,kstate,klspin)
   grad_kl(3,:,klspin) = wf_gradr(3,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(3,:,lstate,klspin) * wf_r(:,kstate,klspin)
   dot_ij_kl(:,ijspin) = grad_ij(1,:,ijspin) * grad_kl(1,:,klspin) + grad_ij(2,:,ijspin) * grad_kl(2,:,klspin) &
                        + grad_ij(3,:,ijspin) * grad_kl(3,:,klspin)
   dot_rho_ij(:,ijspin) = rho_gradr(1,:,1) * grad_ij(1,:,ijspin) + rho_gradr(2,:,1) * grad_ij(2,:,ijspin)  &
                         + rho_gradr(3,:,1) * grad_ij(3,:,ijspin)
   dot_rho_kl(:,klspin) = rho_gradr(1,:,1) * grad_kl(1,:,klspin) + rho_gradr(2,:,1) * grad_kl(2,:,klspin)  &
                         + rho_gradr(3,:,1) * grad_kl(3,:,klspin)

   eval_fxc_rks_singlet = eval_fxc_rks_singlet &
         +  SUM( dot_ij_kl(:,1) * 4.0_dp * vsigma(:,1) ) &
         +  SUM( dot_rho_ij(:,1) * dot_rho_kl(:,1) * 8.0_dp * v2sigma2(:,1) ) &
         +  SUM( ( wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) * dot_rho_kl(:,1)   &
                 + wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) * dot_rho_ij(:,1) ) &
                   * 4.0_dp * v2rhosigma(:,1) )

 endif


end function eval_fxc_rks_singlet


!=========================================================================
function eval_fxc_uks(istate,jstate,ijspin,kstate,lstate,klspin)
 implicit none
!=====
 real(dp)           :: eval_fxc_uks
 integer,intent(in) :: istate,jstate,kstate,lstate
 integer,intent(in) :: ijspin,klspin
!=====

 eval_fxc_uks = SUM(  wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) &
                    * wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) &
                    * ( v2rho2(:,1) + v2rho2(:,2) ) )


 if(require_gradient) then

   grad_ij(1,:,ijspin) = wf_gradr(1,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(1,:,jstate,ijspin) * wf_r(:,istate,ijspin)
   grad_ij(2,:,ijspin) = wf_gradr(2,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(2,:,jstate,ijspin) * wf_r(:,istate,ijspin)
   grad_ij(3,:,ijspin) = wf_gradr(3,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(3,:,jstate,ijspin) * wf_r(:,istate,ijspin)
   grad_kl(1,:,klspin) = wf_gradr(1,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(1,:,lstate,klspin) * wf_r(:,kstate,klspin)
   grad_kl(2,:,klspin) = wf_gradr(2,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(2,:,lstate,klspin) * wf_r(:,kstate,klspin)
   grad_kl(3,:,klspin) = wf_gradr(3,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(3,:,lstate,klspin) * wf_r(:,kstate,klspin)
   dot_ij_kl(:,ijspin) = grad_ij(1,:,ijspin) * grad_kl(1,:,klspin) + grad_ij(2,:,ijspin) * grad_kl(2,:,klspin) &
                        + grad_ij(3,:,ijspin) * grad_kl(3,:,klspin)
   dot_rho_ij(:,ijspin) = rho_gradr(1,:,1) * grad_ij(1,:,ijspin) + rho_gradr(2,:,1) * grad_ij(2,:,ijspin)  &
                         + rho_gradr(3,:,1) * grad_ij(3,:,ijspin)
   dot_rho_kl(:,klspin) = rho_gradr(1,:,1) * grad_kl(1,:,klspin) + rho_gradr(2,:,1) * grad_kl(2,:,klspin)  &
                         + rho_gradr(3,:,1) * grad_kl(3,:,klspin)

   eval_fxc_uks = eval_fxc_uks   &
         +  SUM( dot_ij_kl(:,1) * ( 2.0_dp * vsigma(:,1) + vsigma(:,2) ) ) &
         +  SUM( dot_rho_ij(:,1) * dot_rho_kl(:,1) * ( v2sigma2(:,1) + 0.5_dp * v2sigma2(:,2) +  2.0_dp * v2sigma2(:,3) + v2sigma2(:,4) ) ) &
         +  SUM( ( wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) * dot_rho_kl(:,1)   &
                + wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) * dot_rho_ij(:,1) ) &
                   * ( v2rhosigma(:,1) + v2rhosigma(:,2) + v2rhosigma(:,3) )  )
 endif
    
end function eval_fxc_uks
    

!=========================================================================
function eval_fxc_rks_triplet(istate,jstate,ijspin,kstate,lstate,klspin)
 implicit none
!=====
 real(dp)           :: eval_fxc_rks_triplet
 integer,intent(in) :: istate,jstate,kstate,lstate
 integer,intent(in) :: ijspin,klspin
!=====
 
 eval_fxc_rks_triplet = SUM(  wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) &
                                * wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) &
                                * ( v2rho2(:,1) - v2rho2(:,2) ) )

 if(require_gradient) then

   grad_ij(1,:,ijspin) = wf_gradr(1,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(1,:,jstate,ijspin) * wf_r(:,istate,ijspin)
   grad_ij(2,:,ijspin) = wf_gradr(2,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(2,:,jstate,ijspin) * wf_r(:,istate,ijspin)
   grad_ij(3,:,ijspin) = wf_gradr(3,:,istate,ijspin) * wf_r(:,jstate,ijspin) + wf_gradr(3,:,jstate,ijspin) * wf_r(:,istate,ijspin)
   grad_kl(1,:,klspin) = wf_gradr(1,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(1,:,lstate,klspin) * wf_r(:,kstate,klspin)
   grad_kl(2,:,klspin) = wf_gradr(2,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(2,:,lstate,klspin) * wf_r(:,kstate,klspin)
   grad_kl(3,:,klspin) = wf_gradr(3,:,kstate,klspin) * wf_r(:,lstate,klspin) + wf_gradr(3,:,lstate,klspin) * wf_r(:,kstate,klspin)
   dot_ij_kl(:,ijspin) = grad_ij(1,:,ijspin) * grad_kl(1,:,klspin) + grad_ij(2,:,ijspin) * grad_kl(2,:,klspin) &
                        + grad_ij(3,:,ijspin) * grad_kl(3,:,klspin)
   dot_rho_ij(:,ijspin) = rho_gradr(1,:,1) * grad_ij(1,:,ijspin) + rho_gradr(2,:,1) * grad_ij(2,:,ijspin)  &
                         + rho_gradr(3,:,1) * grad_ij(3,:,ijspin)
   dot_rho_kl(:,klspin) = rho_gradr(1,:,1) * grad_kl(1,:,klspin) + rho_gradr(2,:,1) * grad_kl(2,:,klspin)  &
                         + rho_gradr(3,:,1) * grad_kl(3,:,klspin)

   eval_fxc_rks_triplet = eval_fxc_rks_triplet   &
         +  SUM( dot_ij_kl(:,1) * ( 2.0_dp * vsigma(:,1) - vsigma(:,2) ) ) &
         +  SUM( dot_rho_ij(:,1) * dot_rho_kl(:,1) * ( v2sigma2(:,1) - v2sigma2(:,3) ) ) &   !FIXME 3 or 5 are working, but only one is correct in principle
         +  SUM( ( wf_r(:,istate,ijspin) * wf_r(:,jstate,ijspin) * dot_rho_kl(:,1)   &
                + wf_r(:,kstate,klspin) * wf_r(:,lstate,klspin) * dot_rho_ij(:,1) ) &
                   * ( v2rhosigma(:,1) - v2rhosigma(:,4) )  )   !FIXME 3 and 4 are working, but only one is correct in principle
 endif

end function eval_fxc_rks_triplet
    

!=========================================================================
subroutine destroy_tddft()
 implicit none

 call destroy_dft_grid()

 if( ALLOCATED(v2rho2))     deallocate(v2rho2)
 if( ALLOCATED(vsigma))     deallocate(vsigma)
 if( ALLOCATED(v2rhosigma)) deallocate(v2rhosigma)
 if( ALLOCATED(v2sigma2))   deallocate(v2sigma2)

 if( ALLOCATED(wf_r) )      deallocate(wf_r)
 if( ALLOCATED(wf_gradr) )  deallocate(wf_gradr)
 if( ALLOCATED(rho_gradr) ) deallocate(rho_gradr)

 if( ALLOCATED(grad_ij) )    deallocate(grad_ij)
 if( ALLOCATED(grad_kl) )    deallocate(grad_kl)
 if( ALLOCATED(dot_ij_kl) )  deallocate(dot_ij_kl)
 if( ALLOCATED(dot_rho_ij) ) deallocate(dot_rho_ij)
 if( ALLOCATED(dot_rho_kl) ) deallocate(dot_rho_kl)

end subroutine destroy_tddft


end module m_tddft_fxc


!=========================================================================
