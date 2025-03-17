!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods and data to evaluate the TDDFT kernel f_xc
! and in particular integrals ( i a | f_xc | j b )
!
!=========================================================================
#include "molgw.h"
module m_tddft_fxc
  use m_definitions
  use m_warning
  use m_mpi
  use m_memory
  use m_timing
  use m_inputparam
  use m_dft_grid
  use m_basis_set
  use m_density_tools
  use m_libxc_tools

  integer, private   :: nspin_tddft

  !
  ! fxc kernel
  ! Second derivatives with respect to rho, sigma
  real(dp), allocatable, protected :: v2rho2(:, :)
  real(dp), allocatable, protected :: vsigma(:, :)
  real(dp), allocatable, protected :: v2rhosigma(:, :)
  real(dp), allocatable, protected :: v2sigma2(:, :)

  !
  ! Wfn and gradients evaluated on the grid points
  real(dp), allocatable, protected :: wf_r(:, :, :)
  real(dp), allocatable, protected :: wf_gradr(:, :, :, :)
  real(dp), allocatable, protected :: rho_gradr(:, :, :)

  ! Intermediate quantities
  real(dp), allocatable, protected :: grad_ij(:, :, :)
  real(dp), allocatable, protected :: grad_kl(:, :, :)
  real(dp), allocatable, protected :: dot_ij_kl(:, :)
  real(dp), allocatable, protected :: dot_rho_ij(:, :)
  real(dp), allocatable, protected :: dot_rho_kl(:, :)


contains


!=========================================================================
subroutine prepare_tddft(is_triplet_in, nstate, basis, c_matrix, occupation)
  implicit none

  logical, intent(in)               :: is_triplet_in
  integer, intent(in)               :: nstate
  type(basis_set), intent(in)       :: basis
  real(dp), intent(in)              :: c_matrix(basis%nbf, nstate, nspin)
  real(dp), intent(in)              :: occupation(nstate, nspin)
  !=====
  type(dft_xc_info), allocatable    :: tddft_xc(:)
  real(dp), parameter   :: kernel_capping = 1.0e14_dp ! for numerical stability
  character(len=256)   :: string
  integer              :: ixc, igrid
  integer              :: ispin
  real(dp)             :: basis_function_r(basis%nbf, 1)
  real(dp)             :: bf_gradx(basis%nbf, 1)
  real(dp)             :: bf_grady(basis%nbf, 1)
  real(dp)             :: bf_gradz(basis%nbf, 1)
  real(dp)             :: rhor_r(nspin, 1)
  real(dp)             :: grad_rhor(nspin, 1, 3)
  real(dp)             :: max_v2sigma2
  real(dp), allocatable :: rho_c(:)
  real(dp), allocatable :: v2rho2_c(:)
  real(dp), allocatable :: sigma_c(:)
  real(dp), allocatable :: vrho_c(:)
  real(dp), allocatable :: vsigma_c(:)
  real(dp), allocatable :: v2rhosigma_c(:)
  real(dp), allocatable :: v2sigma2_c(:)
  !=====

#if !defined(NO_LIBXC)

  !
  ! Prepare DFT kernel calculation with Libxc
  !
  nspin_tddft = MERGE(2, nspin, is_triplet_in)
  call copy_libxc_info(dft_xc, tddft_xc)
  do ixc=1, tddft_xc(1)%nxc
    tddft_xc(ixc)%nspin = nspin_tddft
  enddo
  call init_libxc_info(tddft_xc)

  call init_dft_grid(basis, tddft_grid_level, dft_xc(1)%needs_gradient, .FALSE., 1)

  allocate(rho_c(nspin_tddft))
  allocate(v2rho2_c(2*nspin_tddft-1))
  allocate(sigma_c(2*nspin_tddft-1))
  allocate(vrho_c(nspin_tddft))
  allocate(vsigma_c(2*nspin_tddft-1))
  allocate(v2rhosigma_c(5*nspin_tddft-4))
  allocate(v2sigma2_c(5*nspin_tddft-4))

  !
  ! calculate rho, grad rho and the kernel
  !

  allocate(v2rho2(ngrid, 2*nspin_tddft-1), wf_r(ngrid, basis%nbf, nspin))
  v2rho2(:, :) = 0.0_dp

  if( dft_xc(1)%needs_gradient ) then
    allocate(vsigma(ngrid, 2*nspin_tddft-1))
    allocate(v2rhosigma(ngrid, 5*nspin_tddft-4))
    allocate(v2sigma2(ngrid, 5*nspin_tddft-4))
    allocate(wf_gradr(3, ngrid, basis%nbf, nspin))
    allocate(rho_gradr(3, ngrid, nspin))
    vsigma(:, :)     = 0.0_dp
    v2rhosigma(:, :) = 0.0_dp
    v2sigma2(:, :)   = 0.0_dp
  endif


  max_v2sigma2 = -1.0_dp


  do igrid=1, ngrid
    !
    ! Get all the functions and gradients at point rr
    call get_basis_functions_r_batch(basis, igrid, basis_function_r)
    !
    ! store the wavefunction in r
    do ispin=1, nspin
      wf_r(igrid, :, ispin) = MATMUL( basis_function_r(:, 1) , c_matrix(:, :, ispin) )
    enddo

    if( dft_xc(1)%needs_gradient ) then
      call get_basis_functions_gradr_batch(basis, igrid, bf_gradx, bf_grady, bf_gradz)
      !
      ! store the wavefunction in r
      do ispin=1, nspin
        wf_gradr(1, igrid, :, ispin) = MATMUL( bf_gradx(:, 1) , c_matrix(:, :, ispin) )
        wf_gradr(2, igrid, :, ispin) = MATMUL( bf_grady(:, 1) , c_matrix(:, :, ispin) )
        wf_gradr(3, igrid, :, ispin) = MATMUL( bf_gradz(:, 1) , c_matrix(:, :, ispin) )
      enddo
    endif


    if( dft_xc(1)%needs_gradient ) then
      call calc_density_gradr_batch(occupation, c_matrix, basis_function_r, bf_gradx, bf_grady, bf_gradz, rhor_r, grad_rhor)
      rho_gradr(:, igrid, :) = TRANSPOSE(grad_rhor(:, 1, :))

      if( nspin_tddft == 1 ) then
        sigma_c(1) = SUM( grad_rhor(1, 1, :)**2 )
      else if( nspin==2 ) then
        sigma_c(2) = SUM( grad_rhor(1, 1, :) * grad_rhor(2, 1, :) )
        sigma_c(3) = SUM( grad_rhor(2, 1, :)**2 )
      else ! triplet excitations from singlet ground-state
        sigma_c(:) = SUM( grad_rhor(1, 1, :)**2 ) * 0.25_dp
      endif

    else
      call calc_density_r_batch(occupation, c_matrix, basis_function_r, rhor_r)
    endif


    if( nspin_tddft==1 ) then
      rho_c(1) = rhor_r(1, 1)
    else if( nspin==2 ) then
      rho_c(:) = rhor_r(:, 1)
    else ! triplet excitations from singlet ground-state
      rho_c(:) = rhor_r(1, 1) * 0.5_dp
    endif

    !
    ! Calculate the kernel
    !
    do ixc=1, dft_xc(1)%nxc
      if( ABS(dft_xc(ixc)%coeff) < 1.0e-6_dp ) cycle

      select case(dft_xc(ixc)%family)
      case(XC_FAMILY_LDA)
        call xc_lda_fxc(tddft_xc(ixc)%func, 1, rho_c(1), v2rho2_c(1))
      case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        call xc_gga_vxc(tddft_xc(ixc)%func, 1, rho_c(1), sigma_c(1), vrho_c(1), vsigma_c(1))
        call xc_gga_fxc(tddft_xc(ixc)%func, 1, rho_c(1), sigma_c(1), v2rho2_c(1), v2rhosigma_c(1), v2sigma2_c(1))
      case default
        call die('Other kernels not yet implemented')
      end select
      !
      ! Remove the too large values for stability
      v2rho2_c(:) = MIN( v2rho2_c(:), kernel_capping )
      if(dft_xc(1)%needs_gradient) then
        max_v2sigma2 = MAX(ABS(v2sigma2_c(1)), max_v2sigma2)
        vsigma_c(:)     = MIN( vsigma_c(:), kernel_capping )
        v2rhosigma_c(:) = MIN( v2rhosigma_c(:), kernel_capping )
        v2sigma2_c(:)   = MIN( v2sigma2_c(:), kernel_capping )
      endif

      ! Store the result with the weight
      ! Remove too large values for stability
      v2rho2(igrid, :)     = v2rho2(igrid, :) + v2rho2_c(:) * w_grid(igrid) * dft_xc(ixc)%coeff
      if(dft_xc(1)%needs_gradient) then
        vsigma(igrid, :)     = vsigma(igrid, :)     + vsigma_c(:)     * w_grid(igrid) * dft_xc(ixc)%coeff
        v2rhosigma(igrid, :) = v2rhosigma(igrid, :) + v2rhosigma_c(:) * w_grid(igrid) * dft_xc(ixc)%coeff
        v2sigma2(igrid, :)   = v2sigma2(igrid, :)   + v2sigma2_c(:)   * w_grid(igrid) * dft_xc(ixc)%coeff
      endif

    enddo
  enddo
  if(dft_xc(1)%needs_gradient) then
    call world%max(max_v2sigma2)
    write(stdout, '(a,e18.6)') ' Maximum numerical value for fxc: ', max_v2sigma2

    allocate(grad_ij(3, ngrid, nspin))
    allocate(grad_kl(3, ngrid, nspin))
    allocate(dot_ij_kl(ngrid, nspin))
    allocate(dot_rho_ij(ngrid, nspin))
    allocate(dot_rho_kl(ngrid, nspin))

  endif

  deallocate(rho_c)
  deallocate(v2rho2_c)
  deallocate(sigma_c)
  deallocate(vrho_c)
  deallocate(vsigma_c)
  deallocate(v2rhosigma_c)
  deallocate(v2sigma2_c)
  deallocate(tddft_xc)
#else
  call die('prepare_tddft: not available without LIBXC')
#endif

end subroutine prepare_tddft


!=========================================================================
function eval_fxc_rks_singlet(istate, jstate, ijspin, kstate, lstate, klspin)
  implicit none
  !=====
  real(dp)           :: eval_fxc_rks_singlet
  integer, intent(in) :: istate, jstate, kstate, lstate
  integer, intent(in) :: ijspin, klspin
  !=====

  eval_fxc_rks_singlet = SUM(  wf_r(:, istate, ijspin) * wf_r(:, jstate, ijspin) &
                    * wf_r(:, kstate, klspin) * wf_r(:, lstate, klspin) &
                    * v2rho2(:, ijspin) * 2.0_dp )

  if(dft_xc(1)%needs_gradient) then

    grad_ij(1, :, ijspin) = wf_gradr(1, :, istate, ijspin) * wf_r(:, jstate, ijspin) + wf_gradr(1, :, jstate, ijspin) * wf_r(:, istate, ijspin)
    grad_ij(2, :, ijspin) = wf_gradr(2, :, istate, ijspin) * wf_r(:, jstate, ijspin) + wf_gradr(2, :, jstate, ijspin) * wf_r(:, istate, ijspin)
    grad_ij(3, :, ijspin) = wf_gradr(3, :, istate, ijspin) * wf_r(:, jstate, ijspin) + wf_gradr(3, :, jstate, ijspin) * wf_r(:, istate, ijspin)
    grad_kl(1, :, klspin) = wf_gradr(1, :, kstate, klspin) * wf_r(:, lstate, klspin) + wf_gradr(1, :, lstate, klspin) * wf_r(:, kstate, klspin)
    grad_kl(2, :, klspin) = wf_gradr(2, :, kstate, klspin) * wf_r(:, lstate, klspin) + wf_gradr(2, :, lstate, klspin) * wf_r(:, kstate, klspin)
    grad_kl(3, :, klspin) = wf_gradr(3, :, kstate, klspin) * wf_r(:, lstate, klspin) + wf_gradr(3, :, lstate, klspin) * wf_r(:, kstate, klspin)
    dot_ij_kl(:, ijspin) = grad_ij(1, :, ijspin) * grad_kl(1, :, klspin) + grad_ij(2, :, ijspin) * grad_kl(2, :, klspin) &
                        + grad_ij(3, :, ijspin) * grad_kl(3, :, klspin)
    dot_rho_ij(:, ijspin) = rho_gradr(1, :, 1) * grad_ij(1, :, ijspin) + rho_gradr(2, :, 1) * grad_ij(2, :, ijspin)  &
                         + rho_gradr(3, :, 1) * grad_ij(3, :, ijspin)
    dot_rho_kl(:, klspin) = rho_gradr(1, :, 1) * grad_kl(1, :, klspin) + rho_gradr(2, :, 1) * grad_kl(2, :, klspin)  &
                         + rho_gradr(3, :, 1) * grad_kl(3, :, klspin)

    eval_fxc_rks_singlet = eval_fxc_rks_singlet &
         +  SUM( dot_ij_kl(:, 1) * 4.0_dp * vsigma(:, 1) ) &
         +  SUM( dot_rho_ij(:, 1) * dot_rho_kl(:, 1) * 8.0_dp * v2sigma2(:, 1) ) &
         +  SUM( ( wf_r(:, istate, ijspin) * wf_r(:, jstate, ijspin) * dot_rho_kl(:, 1)   &
                 + wf_r(:, kstate, klspin) * wf_r(:, lstate, klspin) * dot_rho_ij(:, 1) ) &
                   * 4.0_dp * v2rhosigma(:, 1) )

  endif


end function eval_fxc_rks_singlet


!=========================================================================
function eval_fxc_uks(istate, jstate, ijspin, kstate, lstate, klspin)
  implicit none
  !=====
  real(dp)           :: eval_fxc_uks
  integer, intent(in) :: istate, jstate, kstate, lstate
  integer, intent(in) :: ijspin, klspin
  !=====

  eval_fxc_uks = SUM(  wf_r(:, istate, ijspin) * wf_r(:, jstate, ijspin) &
                    * wf_r(:, kstate, klspin) * wf_r(:, lstate, klspin) &
                    * ( v2rho2(:, 1) + v2rho2(:, 2) ) )


  if(dft_xc(1)%needs_gradient) then

    grad_ij(1, :, ijspin) = wf_gradr(1, :, istate, ijspin) * wf_r(:, jstate, ijspin) + wf_gradr(1, :, jstate, ijspin) * wf_r(:, istate, ijspin)
    grad_ij(2, :, ijspin) = wf_gradr(2, :, istate, ijspin) * wf_r(:, jstate, ijspin) + wf_gradr(2, :, jstate, ijspin) * wf_r(:, istate, ijspin)
    grad_ij(3, :, ijspin) = wf_gradr(3, :, istate, ijspin) * wf_r(:, jstate, ijspin) + wf_gradr(3, :, jstate, ijspin) * wf_r(:, istate, ijspin)
    grad_kl(1, :, klspin) = wf_gradr(1, :, kstate, klspin) * wf_r(:, lstate, klspin) + wf_gradr(1, :, lstate, klspin) * wf_r(:, kstate, klspin)
    grad_kl(2, :, klspin) = wf_gradr(2, :, kstate, klspin) * wf_r(:, lstate, klspin) + wf_gradr(2, :, lstate, klspin) * wf_r(:, kstate, klspin)
    grad_kl(3, :, klspin) = wf_gradr(3, :, kstate, klspin) * wf_r(:, lstate, klspin) + wf_gradr(3, :, lstate, klspin) * wf_r(:, kstate, klspin)
    dot_ij_kl(:, ijspin) = grad_ij(1, :, ijspin) * grad_kl(1, :, klspin) + grad_ij(2, :, ijspin) * grad_kl(2, :, klspin) &
                        + grad_ij(3, :, ijspin) * grad_kl(3, :, klspin)
    dot_rho_ij(:, ijspin) = rho_gradr(1, :, 1) * grad_ij(1, :, ijspin) + rho_gradr(2, :, 1) * grad_ij(2, :, ijspin)  &
                         + rho_gradr(3, :, 1) * grad_ij(3, :, ijspin)
    dot_rho_kl(:, klspin) = rho_gradr(1, :, 1) * grad_kl(1, :, klspin) + rho_gradr(2, :, 1) * grad_kl(2, :, klspin)  &
                         + rho_gradr(3, :, 1) * grad_kl(3, :, klspin)

    eval_fxc_uks = eval_fxc_uks   &
         +  SUM( dot_ij_kl(:, 1) * ( 2.0_dp * vsigma(:, 1) + vsigma(:, 2) ) ) &
         +  SUM( dot_rho_ij(:, 1) * dot_rho_kl(:, 1)                         &
                 * ( v2sigma2(:, 1) + 0.5_dp * v2sigma2(:, 2) +  2.0_dp * v2sigma2(:, 3) + v2sigma2(:, 4) ) ) &
         +  SUM( ( wf_r(:, istate, ijspin) * wf_r(:, jstate, ijspin) * dot_rho_kl(:, 1)   &
                + wf_r(:, kstate, klspin) * wf_r(:, lstate, klspin) * dot_rho_ij(:, 1) ) &
                   * ( v2rhosigma(:, 1) + v2rhosigma(:, 2) + v2rhosigma(:, 3) )  )
  endif

end function eval_fxc_uks


!=========================================================================
function eval_fxc_rks_triplet(istate, jstate, ijspin, kstate, lstate, klspin)
  implicit none
  !=====
  real(dp)           :: eval_fxc_rks_triplet
  integer, intent(in) :: istate, jstate, kstate, lstate
  integer, intent(in) :: ijspin, klspin
  !=====

  eval_fxc_rks_triplet = SUM(  wf_r(:, istate, ijspin) * wf_r(:, jstate, ijspin) &
                                * wf_r(:, kstate, klspin) * wf_r(:, lstate, klspin) &
                                * ( v2rho2(:, 1) - v2rho2(:, 2) ) )

  if(dft_xc(1)%needs_gradient) then

    grad_ij(1, :, ijspin) = wf_gradr(1, :, istate, ijspin) * wf_r(:, jstate, ijspin) + wf_gradr(1, :, jstate, ijspin) * wf_r(:, istate, ijspin)
    grad_ij(2, :, ijspin) = wf_gradr(2, :, istate, ijspin) * wf_r(:, jstate, ijspin) + wf_gradr(2, :, jstate, ijspin) * wf_r(:, istate, ijspin)
    grad_ij(3, :, ijspin) = wf_gradr(3, :, istate, ijspin) * wf_r(:, jstate, ijspin) + wf_gradr(3, :, jstate, ijspin) * wf_r(:, istate, ijspin)
    grad_kl(1, :, klspin) = wf_gradr(1, :, kstate, klspin) * wf_r(:, lstate, klspin) + wf_gradr(1, :, lstate, klspin) * wf_r(:, kstate, klspin)
    grad_kl(2, :, klspin) = wf_gradr(2, :, kstate, klspin) * wf_r(:, lstate, klspin) + wf_gradr(2, :, lstate, klspin) * wf_r(:, kstate, klspin)
    grad_kl(3, :, klspin) = wf_gradr(3, :, kstate, klspin) * wf_r(:, lstate, klspin) + wf_gradr(3, :, lstate, klspin) * wf_r(:, kstate, klspin)
    dot_ij_kl(:, ijspin) = grad_ij(1, :, ijspin) * grad_kl(1, :, klspin) + grad_ij(2, :, ijspin) * grad_kl(2, :, klspin) &
                        + grad_ij(3, :, ijspin) * grad_kl(3, :, klspin)
    dot_rho_ij(:, ijspin) = rho_gradr(1, :, 1) * grad_ij(1, :, ijspin) + rho_gradr(2, :, 1) * grad_ij(2, :, ijspin)  &
                         + rho_gradr(3, :, 1) * grad_ij(3, :, ijspin)
    dot_rho_kl(:, klspin) = rho_gradr(1, :, 1) * grad_kl(1, :, klspin) + rho_gradr(2, :, 1) * grad_kl(2, :, klspin)  &
                         + rho_gradr(3, :, 1) * grad_kl(3, :, klspin)

    eval_fxc_rks_triplet = eval_fxc_rks_triplet   &
         +  SUM( dot_ij_kl(:, 1) * ( 2.0_dp * vsigma(:, 1) - vsigma(:, 2) ) ) &
         +  SUM( dot_rho_ij(:, 1) * dot_rho_kl(:, 1) * ( v2sigma2(:, 1) - v2sigma2(:, 3) ) ) &   !FIXME 3 or 5 are working, but only one is correct in principle
         +  SUM( ( wf_r(:, istate, ijspin) * wf_r(:, jstate, ijspin) * dot_rho_kl(:, 1)   &
                + wf_r(:, kstate, klspin) * wf_r(:, lstate, klspin) * dot_rho_ij(:, 1) ) &
                   * ( v2rhosigma(:, 1) - v2rhosigma(:, 4) )  )   !FIXME 3 and 4 are working, but only one is correct in principle
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
