!=========================================================================
! This file is part of MOLGW.
! Authors: Fabien Bruneval
!
! This module contains
!
!=========================================================================
#include "molgw.h"
#if !defined(NO_LIBINT)
#include<libint2/libint2_params.h>
#endif

module m_hamiltonian_periodic
  use m_definitions
  use m_timing
  use m_mpi
  use m_scalapack
  use m_warning
  use m_memory
  use m_cart_to_pure
  use m_inputparam
  use m_basis_set
  use m_libint_tools
  use m_libcint_tools
  use m_io
  use m_atoms
  use m_dft_grid
  use m_ecp
  use m_eri
  use m_hamiltonian_tools
  use m_hamiltonian_onebody
#if defined(HAVE_FFTW3)
  use m_fftw3
#endif

  integer, private, parameter :: nx = 1

  integer(C_INT), private :: nfft1
  integer(C_INT), private :: nfft2
  integer(C_INT), private :: nfft3
  integer, private        :: nfft_global
  integer, private        :: nfft_local

  real(dp), allocatable, private :: rgrid(:, :)
  real(dp), allocatable, private :: rhoelecr(:, :)

  real(dp), allocatable, private :: bfr(:, :)

contains


!=========================================================================
! Set FFT grid points so to accommodate variations at Δx scale in real space
!
subroutine set_fft_grid(basis)
  implicit none

  type(basis_set), intent(in) :: basis
  !=====
  integer :: nfftx
  integer :: fft_grids(18) = [ 16, 24, 32, 48, 60, 64, 96, 120, 128, 144, &
                              180, 192, 256, 288, 384, 512, 768, 1024]
  integer :: i1, i2, i3, ifft_global, ifft_local
  integer :: ibf
  real(dp) :: maximum_extension
  !=====

  write(stdout,'(/,1x,a)')    'Setting up the FFT grid'
  write(stdout,'(1x,a,f8.4)') 'spacing fft_delta_x: ', fft_delta_x

  maximum_extension = 0.0_dp
  do ibf=1, basis%nbf
    maximum_extension = MAX(maximum_extension, basis%bff(ibf)%radius)
  enddo
  if( maximum_extension > minimal_image_distance * 0.45_dp) then
    call issue_warning('Basis set too diffuse compared to box dimension')
    write(stdout, '(1x,a,f8.2)') 'Basis extension    (bohr): ', maximum_extension
    write(stdout, '(1x,a,f8.2)') 'Box image distance (bohr): ', minimal_image_distance
  endif

  ! Direction 1
  nfftx = CEILING( NORM2(aprim(:, 1)) / fft_delta_x )
  nfft1 = next_in_grid(nfftx)

  ! Direction 2
  nfftx = CEILING( NORM2(aprim(:, 2)) / fft_delta_x )
  nfft2 = next_in_grid(nfftx)

  ! Direction 3
  nfftx = CEILING( NORM2(aprim(:, 3)) / fft_delta_x )
  nfft3 = next_in_grid(nfftx)

  nfft_global = nfft1 * nfft2 * nfft3

  write(stdout,'(1x,a,1x,i4,a,i4,a,i4)') 'FFT grid:', nfft1, ' x ', nfft2, ' x ', nfft3

  ! Distribute grid points according to dft_grid
  nfft_local = 0
  do ifft_global=1, nfft_global
    if( MODULO(ifft_global - 1, grid%nproc) == grid%rank ) nfft_local = nfft_local + 1
  enddo
  write(stdout,'(1x,a,i9,a,i9)') 'This MPI task treats ', nfft_local , ' grid points out of ', nfft_global

  allocate(rgrid(3, nfft_local))
  ! Create the real-space grid
  ifft_global = 0
  ifft_local  = 0
  do i3=0, nfft3 - 1
    do i2=0, nfft2 - 1
      do i1=0, nfft1 - 1
        ifft_global = ifft_global + 1
        if( MODULO(ifft_global - 1, grid%nproc) == grid%rank ) then
          ifft_local = ifft_local + 1
          rgrid(:, ifft_local) = MATMUL( aprim(:, :), [DBLE(i1)/nfft1, DBLE(i2)/nfft2, DBLE(i3)/nfft3] )
        endif
      enddo
    enddo
  enddo

  allocate(rhoelecr(nfft_local, nspin))
  rhoelecr(:, :) = 0.0_dp

  allocate(bfr(basis%nbf, nfft_local))
  call start_clock(timing_tmp4)
  call calculate_basis_functions_periodic(basis)
  call stop_clock(timing_tmp4)

contains
  ! Find smallest fft_grid greater than fftx
  function next_in_grid(nfftx)
    implicit none

    integer, intent(in) :: nfftx
    integer :: next_in_grid
    !=====
    integer :: i
    !=====

    next_in_grid = -1   ! default value if none found

    do i = 1, SIZE(fft_grids)
      if( fft_grids(i) > nfftx) then
          next_in_grid = fft_grids(i)
          exit
      endif
    enddo

  end function next_in_grid

end subroutine set_fft_grid


!=========================================================================
subroutine setup_overlap_periodic(basis, overlap_ao)
  implicit none

  type(basis_set), intent(in) :: basis
  real(dp), intent(out) :: overlap_ao(:, :)
  !=====
  integer :: i1, i2, i3
  real(dp) :: shift(3)
  real(dp), allocatable :: s_matrix(:, :)
  !=====

  allocate(s_matrix, MOLD=overlap_ao)

  overlap_ao(:, :) = 0.0_dp
  do i3=-nx, nx
    do i2=-nx, nx
      do i1=-nx, nx

        shift(:) = MATMUL( aprim(:, :), [i1, i2, i3] )
        s_matrix(:, :) = 0.0_dp
        call setup_overlap_onecell(basis, shift, s_matrix)
        overlap_ao(:, :) = overlap_ao(:, :) + s_matrix(:, :)

      enddo
    enddo
  enddo

  deallocate(s_matrix)

end subroutine setup_overlap_periodic


!=========================================================================
! Calculate ( \alpha | \beta ) with | beta) that may be shifted in the neighboring cell
!
subroutine setup_overlap_onecell(basis, shift, s_matrix)
  implicit none
  type(basis_set), intent(in) :: basis
  real(dp), intent(in)        :: shift(3)
  real(dp), intent(out)       :: s_matrix(:, :)
  !=====
  integer              :: ishell, jshell
  integer              :: ibf1, ibf2, jbf1, jbf2
  integer              :: ni, nj, ni_cart, nj_cart, li, lj
  real(dp), allocatable :: matrix(:, :)

  real(C_DOUBLE), allocatable :: array_cart(:)
  integer(C_INT)             :: amA, contrdepthA
  real(C_DOUBLE)             :: A(3)
  real(C_DOUBLE), allocatable :: alphaA(:)
  real(C_DOUBLE), allocatable :: cA(:)
  integer(C_INT)             :: amB, contrdepthB
  real(C_DOUBLE)             :: B(3)
  real(C_DOUBLE), allocatable :: alphaB(:)
  real(C_DOUBLE), allocatable :: cB(:)
  !=====
  integer :: i_cart, j_cart, ij
  integer :: ibf_cart, jbf_cart
  !=====

  call start_clock(MERGE(0, timing_overlap, in_rt_tddft))

  do jshell=1, basis%nshell
    lj      = basis%shell(jshell)%am
    nj_cart = number_basis_function_am('CART',lj)
    nj      = number_basis_function_am(basis%gaussian_type, lj)
    jbf1    = basis%shell(jshell)%istart
    jbf2    = basis%shell(jshell)%iend

    call set_libint_shell(basis%shell(jshell), amB, contrdepthB, B, alphaB, cB)
    ! Apply shift here
    B(:) = B(:) + shift(:)

    do ishell=1, basis%nshell
      li      = basis%shell(ishell)%am
      ni_cart = number_basis_function_am('CART',li)
      ni      = number_basis_function_am(basis%gaussian_type, li)
      ibf1    = basis%shell(ishell)%istart
      ibf2    = basis%shell(ishell)%iend

      call set_libint_shell(basis%shell(ishell), amA, contrdepthA, A, alphaA, cA)

      allocate(array_cart(ni_cart*nj_cart))

#if defined(HAVE_LIBCINT)
      call libcint_overlap(amA, contrdepthA, A, alphaA, cA, &
                           amB, contrdepthB, B, alphaB, cB, &
                           array_cart)

      call transform_libint_to_molgw(basis%gaussian_type, li, lj, array_cart, matrix)
#endif

      deallocate(alphaA, cA)


      s_matrix(ibf1:ibf2, jbf1:jbf2) = matrix(:, :)
      s_matrix(jbf1:jbf2, ibf1:ibf2) = TRANSPOSE(matrix(:, :))

      deallocate(array_cart, matrix)

    enddo
    deallocate(alphaB, cB)
  enddo


  call stop_clock(MERGE(0, timing_overlap, in_rt_tddft))


end subroutine setup_overlap_onecell


!=========================================================================
subroutine setup_kinetic_periodic(basis, kin_ao)
  implicit none

  type(basis_set), intent(in) :: basis
  real(dp), intent(out) :: kin_ao(:, :)
  !=====
  integer :: i1, i2, i3
  real(dp) :: shift(3)
  real(dp), allocatable :: hkin(:, :)
  !=====

  allocate(hkin, MOLD=kin_ao)

  kin_ao(:, :) = 0.0_dp
  do i3=-nx, nx
    do i2=-nx, nx
      do i1=-nx, nx

        shift(:) = MATMUL( aprim(:, :), [i1, i2, i3] )
        call setup_kinetic_onecell(basis, shift, hkin)
        kin_ao(:, :) = kin_ao(:, :) + hkin(:, :)

      enddo
    enddo
  enddo

  deallocate(hkin)

end subroutine setup_kinetic_periodic


!=========================================================================
! Calculate ( \alpha | p^2/2 |\beta ) with | beta) that may be shifted in the neighboring cell
!
subroutine setup_kinetic_onecell(basis, shift, kin_ao)
  implicit none
  type(basis_set), intent(in) :: basis
  real(dp), intent(in)        :: shift(3)
  real(dp), intent(out)       :: kin_ao(:, :)
  !=====
  integer              :: ishell, jshell
  integer              :: ibf1, ibf2, jbf1, jbf2
  integer              :: ni, nj, ni_cart, nj_cart, li, lj
  real(dp), allocatable :: matrix(:, :)

  real(C_DOUBLE), allocatable :: array_cart(:)
  integer(C_INT)             :: amA, contrdepthA
  real(C_DOUBLE)             :: A(3)
  real(C_DOUBLE), allocatable :: alphaA(:)
  real(C_DOUBLE), allocatable :: cA(:)
  integer(C_INT)             :: amB, contrdepthB
  real(C_DOUBLE)             :: B(3)
  real(C_DOUBLE), allocatable :: alphaB(:)
  real(C_DOUBLE), allocatable :: cB(:)
  !=====
  integer :: i_cart, j_cart, ij
  integer :: ibf_cart, jbf_cart
  !=====

  call start_clock(MERGE(0, timing_overlap, in_rt_tddft))

  do jshell=1, basis%nshell
    lj      = basis%shell(jshell)%am
    nj_cart = number_basis_function_am('CART',lj)
    nj      = number_basis_function_am(basis%gaussian_type, lj)
    jbf1    = basis%shell(jshell)%istart
    jbf2    = basis%shell(jshell)%iend

    call set_libint_shell(basis%shell(jshell), amB, contrdepthB, B, alphaB, cB)
    ! Apply shift here
    B(:) = B(:) + shift(:)

    do ishell=1, basis%nshell
      li      = basis%shell(ishell)%am
      ni_cart = number_basis_function_am('CART',li)
      ni      = number_basis_function_am(basis%gaussian_type, li)
      ibf1    = basis%shell(ishell)%istart
      ibf2    = basis%shell(ishell)%iend

      call set_libint_shell(basis%shell(ishell), amA, contrdepthA, A, alphaA, cA)


      allocate(array_cart(ni_cart*nj_cart))

#if defined(HAVE_LIBCINT)
      call libcint_kinetic(amA, contrdepthA, A, alphaA, cA, &
                           amB, contrdepthB, B, alphaB, cB, &
                           array_cart)

      call transform_libint_to_molgw(basis%gaussian_type, li, lj, array_cart, matrix)
#endif

      deallocate(alphaA, cA)


      kin_ao(ibf1:ibf2, jbf1:jbf2) = matrix(:, :)
      kin_ao(jbf1:jbf2, ibf1:ibf2) = TRANSPOSE(matrix(:, :))

      deallocate(array_cart, matrix)

    enddo
    deallocate(alphaB, cB)
  enddo


  call stop_clock(MERGE(0, timing_overlap, in_rt_tddft))


end subroutine setup_kinetic_onecell


!=========================================================================
subroutine setup_nucleus_periodic(basis, h_ao, enucnuc)
  implicit none

  type(basis_set), intent(in) :: basis
  real(dp), intent(out) :: h_ao(:, :)
  real(dp), intent(out) :: enucnuc
  !=====
  real(dp) :: nuc_selfenergy
  real(dp) :: rhonuclr(nfft_local)
  real(dp) :: vnuclgrid(nfft_local)
  !=====

#if !defined(HAVE_FFTW3)
  call die('setup_periodic_periodic: requires FFTW at compilation time')
#endif


  if( in_rt_tddft ) then
    call start_clock(timing_tddft_hamiltonian_nuc)
  else
    call start_clock(timing_hamiltonian_nuc)
  endif

  write(stdout, '(/,1x,a)') 'Calculate periodic nucleus term'


  !
  ! Get the electrostatic potential of the nuclei
  !
  call start_clock(timing_tmp1)
!  call prepare_nuclei_density_periodic(rhonuclr, nuc_selfenergy)
  ! Analytic for GTH pseudos only
  call prepare_nuclei_density_analytic_periodic(rhonuclr, nuc_selfenergy)
  call stop_clock(timing_tmp1)

  call poisson_solver_fft(rhonuclr, vnuclgrid) 

  enucnuc = 0.5_dp * SUM( rhonuclr(:) * vnuclgrid(:) ) * volume / REAL(nfft_global, KIND=dp)
  call grid%sum(enucnuc)
  enucnuc = enucnuc - nuc_selfenergy
  write(stdout,'(1x,a,f16.8)') 'Nucleus-nucleus (Ha): ', enucnuc

  call calculate_hao_periodic(basis, vnuclgrid, h_ao)

  call dump_out_matrix(.FALSE., '=== Nuclei contribution (FFTW) ===', h_ao)

  if( in_rt_tddft ) then
    call stop_clock(timing_tddft_hamiltonian_nuc)
  else
    call stop_clock(timing_hamiltonian_nuc)
  endif


end subroutine setup_nucleus_periodic


!=========================================================================
subroutine setup_hartree_periodic(basis, p_matrix, h_ao, ehartree)
  implicit none

  type(basis_set), intent(in) :: basis
  class(*), intent(in)  :: p_matrix(:, :, :)
  real(dp), intent(out) :: h_ao(:, :)
  real(dp), intent(out) :: ehartree
  !=====
  integer               :: timing_xxdft_hartree
  real(dp), allocatable :: p_matrix_local(:, :, :)
  real(dp) :: vhartreegrid(nfft_local, nspin)
  !real(dp) :: dr(3, 3)
  !=====

#if !defined(HAVE_FFTW3)
  call die('setup_hartree_periodic: requires FFTW at compilation time')
#endif

  select type(p_matrix)
  type is(real(dp))
    timing_xxdft_hartree   = timing_hartree
  type is(complex(dp))
    timing_xxdft_hartree   = timing_tddft_hartree
  class default
    call die("setup_hartree_periodic: p_matrix is neither real nor complex")
  end select

  call start_clock(timing_xxdft_hartree)

  write(stdout, '(/,1x,a)') 'Calculate periodic Hartree term'


  allocate(p_matrix_local(basis%nbf, basis%nbf, nspin))
  select type(p_matrix)
  type is(real(dp))
    p_matrix_local(:, :, :) = p_matrix(:, :, :)
  end select


  call start_clock(timing_tmp2)

  call calculate_density_periodic(basis, p_matrix_local, rhoelecr)

  call stop_clock(timing_tmp2)

  
  !
  ! Get the electrostatic potential of the electrons (Hartree)
  !

  !DEBUG
  !dr(:, 1) = aprim(:, 1) / nfft1
  !dr(:, 2) = aprim(:, 2) / nfft2
  !dr(:, 3) = aprim(:, 3) / nfft3
  !call write_cube_file('rhoelecr.cube', nfft1, nfft2, nfft3, dr, rhoelecr(:, 1), comment='test')

  call poisson_solver_fft(rhoelecr(:, 1), vhartreegrid(:, 1)) 

  ehartree = 0.5_dp * SUM( rhoelecr * vhartreegrid ) * volume / REAL(nfft_global, KIND=dp)
  call grid%sum(ehartree)
  write(stdout,'(1x,a,f16.8)') 'Hartree energy (Ha): ', ehartree

  ! Spin channel 1 will contain the total Hartree potential
  vhartreegrid(:, 1) = SUM(vhartreegrid(:, :), DIM=2)

  call calculate_hao_periodic(basis, vhartreegrid(:, 1), h_ao)
  call dump_out_matrix(.FALSE., '=== Hartree contribution (FFTW) ===', h_ao)



  call stop_clock(timing_xxdft_hartree)


end subroutine setup_hartree_periodic


!=========================================================================
subroutine setup_vxc_periodic(basis, vxc_ao, exc_xc, dft_xc_in)
  implicit none

  type(basis_set), intent(in) :: basis
  real(dp), intent(out)       :: vxc_ao(:, :, :)
  real(dp), intent(out)       :: exc_xc
  type(dft_xc_info), intent(in), optional :: dft_xc_in(:)
  !=====
  real(dp), parameter            :: TOL_RHO = 1.0e-9_dp
  type(dft_xc_info), allocatable :: dft_xc_local(:)
  integer              :: nstate
  integer              :: ispin
  integer              :: ixc
  integer              :: igrid_start, igrid_end
  integer              :: timing_xxdft_xc, timing_xxdft_densities, timing_xxdft_libxc, timing_xxdft_vxc
  !real(dp), allocatable :: tmp_batch(:, :)
  !real(dp), allocatable :: basis_function_r_batch(:, :)
  !real(dp), allocatable :: bf_gradx_batch(:, :)
  !real(dp), allocatable :: bf_grady_batch(:, :)
  !real(dp), allocatable :: bf_gradz_batch(:, :)
  !real(dp), allocatable :: dedd_r_batch(:, :)
  !real(dp), allocatable :: grad_rhor_batch(:, :, :)
  !real(dp), allocatable :: dedgd_r_batch(:, :, :)
  real(C_DOUBLE), allocatable :: rhor(:, :)
  !real(C_DOUBLE), allocatable :: sigma_batch(:, :)
  real(C_DOUBLE), allocatable :: vrho(:, :)
  real(C_DOUBLE), allocatable :: excgrid(:)
  !real(C_DOUBLE), allocatable :: vsigma_batch(:, :)
  real(dp), allocatable :: vxcgrid(:, :)
  real(dp) :: vxc_avg
  !=====

  !
  ! create a local copy of dft_xc (from m_inputparam) or dft_xc_in (optional argument)
  !
  if( PRESENT(dft_xc_in) ) then
    call copy_libxc_info(dft_xc_in, dft_xc_local)
  else
    call copy_libxc_info(dft_xc, dft_xc_local)
  endif

  allocate(vxcgrid(nfft_local, nspin))
  exc_xc = 0.0_dp
  vxcgrid(:, :) = 0.0_dp

  if( dft_xc_local(1)%nxc == 0 ) return

  timing_xxdft_xc    = timing_dft_xc
  timing_xxdft_libxc = timing_dft_libxc
  timing_xxdft_vxc   = timing_dft_vxc

  call start_clock(timing_xxdft_xc)
  write(stdout, '(/,1x,a)') 'Calculate DFT XC periodic potential'

  call start_clock(timing_xxdft_libxc)

  allocate(rhor(nspin, nfft_local))
  allocate(vrho(nspin, nfft_local))
  allocate(excgrid(nfft_local))

  ! Need to transpose rho to adapt to libxc convention
  ! nspin x nfft_local
  rhor(:, :) = TRANSPOSE(rhoelecr(:, :))

  do ixc=1, dft_xc_local(1)%nxc
    if( ABS(dft_xc_local(ixc)%coeff) < 1.0e-6_dp ) cycle

    select case(dft_xc_local(ixc)%family)
    case(XC_FAMILY_LDA)
      call xc_lda_exc_vxc(dft_xc_local(ixc)%func, INT(nfft_local, KIND=C_INT), rhor(1, 1), excgrid(1), vrho(1, 1))
    case default
      call die('setup_vxc_periodic: only LDA is implemented right now')
    end select  

    vxcgrid(:, :) = vxcgrid(:, :) + dft_xc_local(ixc)%coeff * TRANSPOSE(vrho(:, :))
    exc_xc = exc_xc + dft_xc_local(ixc)%coeff * SUM(SUM(rhor(:, :), DIM=1) * excgrid(:)) * volume / REAL(nfft_global, dp)

  enddo 

  call grid%sum(exc_xc)

  deallocate(vrho)
  deallocate(excgrid)
  call stop_clock(timing_xxdft_libxc)

  write(stdout, '(1x,a,f14.8)') 'Exc LDA (Ha): ', exc_xc

  vxc_avg = SUM(vxcgrid(:, :)) / REAL(nspin * nfft_global, dp)
  call grid%sum(vxc_avg)
  write(stdout, '(1x,a,f12.5)') 'Average vxc (eV): ', vxc_avg * Ha_eV

  call start_clock(timing_xxdft_vxc)
  call calculate_hao_periodic(basis, vxcgrid, vxc_ao)
  call stop_clock(timing_xxdft_vxc)


  call stop_clock(timing_xxdft_xc)


end subroutine setup_vxc_periodic


!=========================================================================
subroutine calculate_density_periodic(basis, p_matrix, rhor)
  implicit none

  type(basis_set), intent(in) :: basis
  real(dp), intent(in)  :: p_matrix(:, :, :)
  real(dp), intent(out) :: rhor(:, :)
  !=====
  integer :: ispin
  real(dp) :: shift(3)
  real(dp), allocatable :: tmp(:, :)
  !=====

  allocate(tmp(basis%nbf, nfft_local))

  rhor(:, :) = 0.0_dp
  do ispin=1, nspin
    ! phi_α(r) * ( P_αβ * phi_β(r) )
    call DGEMM('N', 'N', basis%nbf, nfft_local, basis%nbf, 1.0d0, p_matrix(:, :, ispin), basis%nbf, &
               bfr(:, :), basis%nbf, 0.0d0, tmp, basis%nbf)
  
    rhor(:, ispin) = rhor(:, ispin) + SUM( bfr(:, :) * tmp(:, :), DIM=1)
  enddo

  deallocate(tmp)

end subroutine calculate_density_periodic


!=========================================================================
subroutine calculate_hao_periodic(basis, vr, h_ao)
  implicit none

  type(basis_set), intent(in) :: basis
  real(dp), intent(in) :: vr(..)
  real(dp), intent(out) :: h_ao(..)
  !=====
  integer :: nx=1
  integer :: ispin, ibf, jbf, ifft_local
  real(dp) :: shift(3)
  !=====



  call start_clock(timing_tmp5)

  select rank(h_ao)
  !
  ! vr and h_ao have a spin index
  !
  rank(3)
    select rank(vr)
      rank(2)
      h_ao(:, :, :) = 0.0_dp
      do ispin=1, nspin
        !do jbf=1, basis%nbf
        !  do ibf=1, basis%nbf
        !    h_ao(ibf, jbf, ispin) = h_ao(ibf, jbf, ispin) + SUM( bfr(ibf, :) * vr(:, ispin) * bfr(jbf, :) )
        !  enddo
        !enddo

        do ifft_local=1, nfft_local
          call DSYR('L', basis%nbf, vr(ifft_local, ispin), bfr(:, ifft_local), 1, h_ao(:, :, ispin), basis%nbf)
        enddo

        !block
        !real(dp), allocatable :: tmp(:, :)
        !allocate(tmp, MOLD=bfr)
        !do concurrent(ifft_local=1:nfft_local)
        !  tmp(:, ifft_local) = bfr(:, ifft_local) * vr(ifft_local, ispin)
        !enddo
        !call DSYR2K('L', 'N', basis%nbf, nfft_local, 0.5d0, bfr, basis%nbf, &
        !            tmp, basis%nbf, 0.0d0, h_ao(:, :, ispin), basis%nbf)
        !deallocate(tmp)
        !end block

        call matrix_lower_to_full(h_ao(:, :, ispin))

      enddo
      h_ao(:, :, :) = h_ao(:, :, :) * volume / REAL(nfft_global, KIND=dp)
      call grid%sum(h_ao)

    rank(1)
      call die('calculate_hao_periodic: rank mismatch') 
    end select

  !
  ! vr and h_ao don't have a spin index
  !
  rank(2)
    select rank(vr)
    rank(1)
      h_ao(:, :) = 0.0_dp
      !do jbf=1, basis%nbf
      !  do ibf=1, basis%nbf
      !    h_ao(ibf, jbf) = h_ao(ibf, jbf) + SUM( bfr(ibf, :) * vr(:) * bfr(jbf, :) )
      !  enddo
      !enddo
      do ifft_local=1, nfft_local
        call DSYR('L', basis%nbf, vr(ifft_local), bfr(:, ifft_local), 1, h_ao(:, :), basis%nbf)
      enddo
      call matrix_lower_to_full(h_ao(:, :))

      h_ao(:, :) = h_ao(:, :) * volume / REAL(nfft_global, KIND=dp)
      call grid%sum(h_ao)

    rank(2)
      call die('calculate_hao_periodic: rank mismatch') 
    end select
  end select

  call stop_clock(timing_tmp5)


end subroutine calculate_hao_periodic


!=========================================================================
! Evaluate the nuclei charge density on the FFT grid
! that induces the periodic local potential
!
! Convention: nuclei have a negative charge (since electrons have positive one)
!
subroutine prepare_nuclei_density_periodic(rhonuclr, selfenergy)
  implicit none

  real(dp), intent(out) :: rhonuclr(:)
  real(dp), intent(out) :: selfenergy
  !=====
  integer :: nrad, irad, irmin, irmax
  real(dp), allocatable :: rrad(:), vloc(:), dvdr(:), d2vdr2(:), rho(:)
  real(dp), allocatable :: rradfinal(:), rhofinal(:)
  real(dp) :: zval, rmax, dr, shift(3), factor
  integer :: i1, i2, i3, icenter, ifft_local
  real(dp) :: agrid, bgrid, irad_float
  !=====


  call read_potential_data('vloc.dat', nrad, rrad, vloc)

  allocate(dvdr, MOLD=vloc)
  call calc_derivative(rrad, vloc, dvdr)
  allocate(d2vdr2, MOLD=vloc)
  call calc_derivative(rrad, dvdr, d2vdr2)


  allocate(rho, MOLD=vloc)
  do irad=1, nrad
    rho(irad) = -( d2vdr2(irad) + 2.0_dp * dvdr(irad) / rrad(irad)) / (4.0_dp * pi)
  enddo
  deallocate(dvdr, d2vdr2)

  ! clean up rho
  ! remove the too low density values at large r
  ! remove the too low r

  do irad=nrad, 1, -1
    if( ABS(rho(irad)) > 1.0e-6_dp ) exit
  enddo
  irmax = irad
  do irad=1, nrad
    if( rrad(irad) > 1.0e-3_dp ) exit
  enddo
  irmin = irad

  nrad = irmax - irmin + 1
  allocate(rradfinal(nrad))
  allocate(rhofinal(nrad))
  rradfinal(1) = 0.0_dp
  rhofinal(1) = SUM(rho(1:irmin)) / REAL(irmin, KIND=dp)
  rradfinal(2:nrad) = rrad(irmin+1:irmin+nrad-1)
  rhofinal(2:nrad) = rho(irmin+1:irmin+nrad-1)
  rmax = rradfinal(nrad)

  zval = calculate_total_charge(rradfinal, rhofinal)
  write(stdout,'(1x,a,f12.8)') 'Integrated nucleus charge: ', zval
  ! Renormalize charge to correct value
  rhofinal(:) = rhofinal(:) / zval * (-1.0_dp)

  call calculate_selfenergy(rradfinal, rhofinal, zval, selfenergy)
  write(stdout,'(1x,a,f12.8)') 'Integrated nucleus charge (after renormalization): ', zval
  write(stdout,'(1x,a,f12.8)') 'Selfenergy per atom (Ha): ', selfenergy
  
  selfenergy = selfenergy * ncenter_nuclei
  write(stdout,'(1x,a,f12.8)') 'Selfenergy (Ha): ', selfenergy


  !DEBUG
  !do irad=1, nrad
  !  write(1001,'(*(es20.10,1x))') rradfinal(irad), rhofinal(irad)
  !enddo

  !
  ! Loop over the unit-cell regular grid
  !

  ! Assume a logarithmic grid (and check it is correct)
  bgrid = rradfinal(2)
  agrid = rradfinal(3) / rradfinal(2)
  if( ABS( rradfinal(101) / rradfinal(100) - agrid ) > 1.0e-6_dp ) then
    call die('prepare_nuclei_density_periodic: vloc grid is not logarithmic. Not implemented')
  endif

  call start_clock(timing_tmp9)
  rhonuclr(:) = 0.0_dp

  !$OMP PARALLEL PRIVATE(shift, dr, irad_float, irad)
  !$OMP DO
  do ifft_local=1, nfft_local
    do i3=-nx, nx
      do i2=-nx, nx
        do i1=-nx, nx
          shift(:) = MATMUL( aprim(:, :), [i1, i2, i3] )
          do icenter=1, ncenter_nuclei
            dr = NORM2( rgrid(:, ifft_local) - xatom(:, icenter) - shift(:) )
            if( dr > rmax ) cycle

            ! Calculate the radial indices [irad-1, irad] that braket dr
            irad_float = ( LOG(dr) - LOG(bgrid) ) / LOG(agrid) + 2.0d0
            irad = CEILING(irad_float)
            if( irad > nrad ) cycle

            !do irad=1, nrad
            !  if( rradfinal(irad) > dr ) exit
            !enddo

            rhonuclr(ifft_local) = rhonuclr(ifft_local) &
                                   + ( rradfinal(irad) - dr   ) / ( rradfinal(irad) - rradfinal(irad-1) ) * rhofinal(irad-1) &
                                   + ( dr - rradfinal(irad-1) ) / ( rradfinal(irad) - rradfinal(irad-1) ) * rhofinal(irad)

          enddo
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  call stop_clock(timing_tmp9)

 
  zval = -SUM(rhonuclr(:)) * volume / REAL(nfft_global, KIND=dp)
  call grid%sum(zval)
  write(stdout, '(1x,a,f12.6)') 'Nuclei charge in the cell: ', zval
  write(stdout, '(1x,a,f12.6)') 'whereas it should be: ', SUM(zvalence(:))

  factor = SUM(zvalence(:)) / zval
 
  ! If the deviation from expected charge is larger than 5%, then stop and ask for a finer fft_delta_x
  if( ABS( factor - 1.0_dp ) > 0.05_dp ) then
    call die('prepare_nuclei_density_periodic: wrong nucleus total charge. ' // &
             'FFT grid is certainly too coarse. Try to decrease fft_delta_x.')
  endif

  rhonuclr(:) = rhonuclr(:) * factor
  zval = -SUM(rhonuclr(:)) * volume / REAL(nfft_global, KIND=dp)
  call grid%sum(zval)
  write(stdout, '(1x,a,f12.6)') 'Nuclei charge in the cell (after renormalization): ', zval


contains

  subroutine read_potential_data(filename, nrad, r, phi)
    implicit none
    character(len=*), intent(in) :: filename
  
    !=====
    integer, intent(out) :: nrad
    real(dp), allocatable, intent(out) :: r(:), phi(:)
    integer :: iostat, irad
    integer :: file_unit
    !=====
    
    ! Count lines
    open(newunit=file_unit, file=filename, status='old', action='read')
    nrad = 0
    do
        read(file_unit, *, iostat=iostat)
        if (iostat /= 0) exit
        nrad = nrad + 1
    enddo
    rewind(file_unit)
    
    ! Allocate and read
    allocate(r(nrad), phi(nrad))
    do irad = 1, nrad
        read(file_unit, *) r(irad), phi(irad)
    enddo
    close(file_unit)
      
  end subroutine read_potential_data
  
  subroutine calc_derivative(r, phi, dphi_dr)
    real(dp), intent(in) :: r(:), phi(:)
    real(dp), intent(out) :: dphi_dr(:)
    real(dp) :: h1, h2, a, b, c
    integer :: i, n
  
    n = SIZE(r)
    
    ! First derivative using non-uniform grid formulas
    ! Forward difference for first point
    h1 = r(2) - r(1)
    h2 = r(3) - r(2)
    a = -(2.0_dp*h1 + h2) / (h1 * (h1 + h2))
    b = (h1 + h2) / (h1 * h2)
    c = -h1 / (h2 * (h1 + h2))
    dphi_dr(1) = a*phi(1) + b*phi(2) + c*phi(3)
    
    ! Central difference for interior points
    do i = 2, n-1
        h1 = r(i) - r(i-1)
        h2 = r(i+1) - r(i)
        a = -h2 / (h1 * (h1 + h2))
        b = (h2 - h1) / (h1 * h2)
        c = h1 / (h2 * (h1 + h2))
        dphi_dr(i) = a*phi(i-1) + b*phi(i) + c*phi(i+1)
    enddo
    
    ! Backward difference for last point
    h1 = r(n-1) - r(n-2)
    h2 = r(n) - r(n-1)
    a = h2 / (h1 * (h1 + h2))
    b = -(h1 + h2) / (h1 * h2)
    c = (2.0_dp*h2 + h1) / (h2 * (h1 + h2))
    dphi_dr(n) = a*phi(n-2) + b*phi(n-1) + c*phi(n)
  
  end subroutine calc_derivative
  
  function calculate_total_charge(r, rho) result(total_charge)
    real(dp), intent(in) :: r(:), rho(:)
    real(dp) :: total_charge
    integer :: irad, nrad
  
    nrad = SIZE(rho)
    ! Integrate 4πr²ρ(r) dr using trapezoidal rule
    total_charge = 0.0_dp
    do irad = 1, nrad - 1
      total_charge = total_charge &
           + 4.0_dp * pi * 0.5_dp * ( r(irad + 1) - r(irad) )   &
               * ( r(irad)**2 * rho(irad) + r(irad + 1)**2 * rho(irad + 1) )
    enddo
  end function calculate_total_charge
  
  subroutine calculate_selfenergy(r, rho, total_charge, selfenergy)
    real(dp), intent(in) :: r(:), rho(:)
    real(dp) :: total_charge, selfenergy
    !=====
    integer  :: irad, nrad
    real(dp), allocatable :: Q(:), S(:), phi(:)
    real(dp) :: dr1
    !=====
  
    nrad = SIZE(rho)
  
    allocate(Q(nrad), S(nrad), phi(nrad))
    
    !------------------------------------------------------------
    ! 1. Compute cumulative charge Q(r)
    !    Q(i) = 4π ∫0^r(i) r'^2 ρ(r') dr'   (trapezoidal, nonuniform)
    !------------------------------------------------------------
    Q(1) = 0.0_dp
    do irad = 2, nrad
       dr1 = r(irad) - r(irad-1)
       Q(irad) = Q(irad-1) + 0.5_dp * 4.0_dp*pi * (r(irad)**2*rho(irad) + r(irad-1)**2*rho(irad-1)) * dr1
    end do
    total_charge = Q(nrad)
    
    !------------------------------------------------------------
    ! 2. Compute reverse cumulative S(r)
    !    S(i) = 4π ∫_r(i)^∞ r' ρ(r') dr'
    !------------------------------------------------------------
    S(nrad) = 0.0_dp
    do irad = nrad-1, 1, -1
       dr1 = r(irad+1) - r(irad)
       S(irad) = S(irad+1) + 0.5_dp * 4.0_dp * pi * ( r(irad) * rho(irad) + r(irad+1) * rho(irad+1) ) * dr1
    end do
    
    !------------------------------------------------------------
    ! 3. Potential φ(r) = Q(r)/r + S(r)
    !------------------------------------------------------------
    phi(1) = S(1)
    phi(2:nrad) = Q(2:nrad) / r(2:nrad) + S(2:nrad)
    
    !------------------------------------------------------------
    ! 4. Energy integral: U = 1/2 * ∫ 4π r^2 ρ(r) φ(r) dr
    !------------------------------------------------------------
    selfenergy = 0.0_dp
    do irad = 2, nrad
       dr1 = r(irad) - r(irad-1)
       selfenergy = selfenergy &
                      + 0.5_dp * 4.0_dp * pi * (  r(irad)**2   * rho(irad)   * phi(irad) &
                                                + r(irad-1)**2 * rho(irad-1) * phi(irad-1) ) * dr1 * 0.5_dp
    end do
    
    deallocate(Q, S, phi)

  end subroutine calculate_selfenergy


end subroutine prepare_nuclei_density_periodic


!=========================================================================
! Evaluate the nuclei charge density on the FFT grid
! that induces the periodic local potential
!
! Convention: nuclei have a negative charge (since electrons have positive one)
!
subroutine prepare_nuclei_density_analytic_periodic(rhonuclr, selfenergy)
  implicit none

  real(dp), intent(out) :: rhonuclr(:)
  real(dp), intent(out) :: selfenergy
  !=====
  real(dp) :: zval, rmax, dr, shift(3), factor
  integer :: i1, i2, i3, icenter, ifft_local, iloc, ie

  real(dp) :: factor_zeff, factor_c1, factor_c2, factor_c3, factor_c4
  real(dp) :: zeff, rloc, ci(4)
  real(dp) :: alphapp, selfenergy_nucleus
  logical :: element_has_ecp
  !=====

  selfenergy = 0.0_dp
  rhonuclr(:) = 0.0_dp

  do icenter=1, ncenter_nuclei

    element_has_ecp = .FALSE.
    do ie=1, nelement_ecp
      if( ABS( element_ecp(ie) - zatom(icenter) ) < 1.0e-5_dp ) then
        element_has_ecp = .TRUE.
        exit
      endif
    enddo
 
    if( ecp(ie)%ecp_format /= ECP_GTH ) call die('prepare_nuclei_density_analytic_periodic: only for GTH pseudos')


    zeff = -zvalence(icenter)
    rloc = ecp(ie)%gth_rpploc
    ci(:) = 0.0_dp
    do iloc=1, ecp(ie)%gth_nloc
      ci(iloc) = ecp(ie)%gth_cipp(iloc)
    enddo
    alphapp = 1.0_dp / SQRT(2.0_dp) / rloc

    !
    ! Loop over the unit-cell regular grid
    !
    call start_clock(timing_tmp9)

    factor_zeff = zeff * ( alphapp / SQRT(pi) )**3
    factor_c1 = ci(1) * alphapp**2 / (2.0_dp * pi)
    factor_c2 = ci(2) * alphapp**2 / pi
    factor_c3 = ci(3) * alphapp**2 * 2.0_dp / pi
    factor_c4 = ci(4) * alphapp**2 * 4.0_dp / pi

    !$OMP PARALLEL PRIVATE(shift, dr, irad_float, irad)
    !$OMP DO
    do ifft_local=1, nfft_local
      do i3=-nx, nx
        do i2=-nx, nx
          do i1=-nx, nx
            shift(:) = MATMUL( aprim(:, :), [i1, i2, i3] )

            dr = NORM2( rgrid(:, ifft_local) - xatom(:, icenter) - shift(:) )

            ! Skip points that are too far, where the density is close to zero
            if( dr > 10.0_dp * rloc ) cycle

            rhonuclr(ifft_local) = rhonuclr(ifft_local) + gth_rhonucl(dr)

          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call stop_clock(timing_tmp9)

    !selfenergy = zeff**2 * alphapp / SQRT(2.0_dp * pi) &
    !            + c1**2 / alphapp * 3.0_dp / 16.0_dp * SQRT(pi / 2.0_dp) &
    !            + c2**2 / alphapp * 33.0_dp / 128.0_dp * SQRT(2.0_dp * pi ) &
    !            + zeff * c1 / ( 2.0_dp * SQRT(2.0_dp)) &
    !            + zeff * c2 * 3.0_dp / 8.0_dp  * SQRT(2.0_dp) &
    !            + c1 * c2 / alphapp * 3.0_dp / 32.0_dp  * SQRT(2.0_dp * pi)

    !write(stdout,'(1x,a,f12.8)') 'Selfenergy per atom (Ha): ', selfenergy

    selfenergy_nucleus = selfenergy_quadrature()
    write(stdout,'(1x,a,i4,a,f12.8)') 'Nucleus selfenergy for center ', icenter, ' (Ha): ', selfenergy_nucleus

    selfenergy = selfenergy + selfenergy_nucleus

  enddo ! icenter
 
  zval = -SUM(rhonuclr(:)) * volume / REAL(nfft_global, KIND=dp)
  call grid%sum(zval)

  write(stdout, '(1x,a,f12.6)') 'Nuclei charge in the cell: ', zval
  write(stdout, '(1x,a,f12.6)') 'whereas it should be: ', SUM(zvalence(:))

  factor = SUM(zvalence(:)) / zval

  ! If the deviation from expected charge is larger than 1%, then stop and ask for a finer fft_delta_x
  if( ABS( factor - 1.0_dp ) > 0.01_dp ) then
    call die('prepare_nuclei_density_periodic: wrong nucleus total charge. ' // &
             'FFT grid is certainly too coarse. Try to decrease fft_delta_x.')
  endif

  rhonuclr(:) = rhonuclr(:) * factor
  zval = -SUM(rhonuclr(:)) * volume / REAL(nfft_global, KIND=dp)
  call grid%sum(zval)
  write(stdout, '(1x,a,f12.6)') 'Nuclei charge in the cell (after renormalization): ', zval

  write(stdout,'(1x,a,f12.8)') 'Total nucleus selfenergy (Ha): ', selfenergy

contains

pure function gth_potnucl(rr) result(potr)
  implicit none

  real(dp), intent(in) :: rr
  real(dp)             :: potr
  !=====
  real(dp) ::  xx
  !=====
 
  xx = (alphapp * rr)**2
  potr = zeff * ERF(alphapp * rr) / rr &
         + EXP(-xx) * (                     ci(1) &
                       + (2.0_dp * xx)    * ci(2) &
                       + (2.0_dp * xx)**2 * ci(3) &
                       + (2.0_dp * xx)**3 * ci(4) )  

end function gth_potnucl

pure function gth_rhonucl(rr) result(rhor)
  implicit none

  real(dp), intent(in) :: rr
  real(dp)             :: rhor
  !=====
  real(dp) ::  xx
  !=====
 
  xx = (alphapp * rr)**2
  rhor = EXP( - xx )  &
           * ( &
                factor_zeff &
              + factor_c1 *         (  3.0_dp -  2.0_dp * xx ) &
              + factor_c2 *         ( -3.0_dp +  7.0_dp * xx - 2.0_dp * xx**2) &
              + factor_c3 * xx    * (-10.0_dp + 11.0_dp * xx - 2.0_dp * xx**2) &
              + factor_c4 * xx**2 * (-21.0_dp + 15.0_dp * xx - 2.0_dp * xx**2) &
             )

end function gth_rhonucl


pure function selfenergy_quadrature() result(integral)
  implicit none

  real(dp) :: integral
  !=====
  real(dp) :: stepr, rr
  integer, parameter :: nrad = 10000
  integer :: irad
  !=====

  stepr = 10.0_dp * rloc / REAL(nrad, KIND=dp)
  
  rr = -stepr / 2.0_dp
  stepr = 0.01_dp
  integral = 0.0_dp
  do irad=1, nrad
    rr = rr + stepr
    integral = integral + 0.5_dp * stepr * 4.0_dp * pi * rr**2  &
                            * gth_rhonucl(rr) * gth_potnucl(rr)
  enddo


end function selfenergy_quadrature


end subroutine prepare_nuclei_density_analytic_periodic


!=========================================================================
! Poisson solver on a FFT grid 
! input:  ρ(r)
! output: vcoul(r)
!
! This routine works on the full FFT grid (not the MPI distributed one).
! Communications are needed and FFTs do not scalable with MPI tasks
!
subroutine poisson_solver_fft(rhor, vcoulr)
  implicit none

  real(dp), intent(in)  :: rhor(:)
  real(dp), intent(out) :: vcoulr(:)
  !=====
#if defined(HAVE_FFTW3)
  type(C_PTR) :: plan
  type(C_PTR) :: pr, pg, pvr, pvg
  real(C_DOUBLE), pointer            :: rhor_fftw(:, :, :), vr_fftw(:, :, :)
  complex(C_DOUBLE_COMPLEX), pointer :: rhog_fftw(:, :, :), vg_fftw(:, :, :)
  integer :: ig1, ig2, ig3, i1, i2, i3, ifft_local, ifft_global
  !=====

  ! Allocate arrays with FFTW functions
  pr = fftw_alloc_real(INT(nfft_global, C_SIZE_T))
  call C_F_POINTER(pr, rhor_fftw, [nfft1, nfft2, nfft3])
  pg = fftw_alloc_complex(INT( (nfft1/2+1) * nfft2 * nfft3, C_SIZE_T))
  call C_F_POINTER(pg, rhog_fftw, [(nfft1/2+1), nfft2, nfft3])

  !
  ! Gather the densities from the different MPI tasks
  ifft_local = 0
  ifft_global = 0
  rhor_fftw(:, :, :) = 0.0_dp
  do i3=1, nfft3
    do i2=1, nfft2
      do i1=1, nfft1
        ifft_global = ifft_global + 1
        if( MODULO(ifft_global - 1, grid%nproc) == grid%rank ) then
          ifft_local = ifft_local + 1
          rhor_fftw(i1, i2, i3) = rhor(ifft_local)
        endif
      enddo
    enddo
  enddo
  call grid%sum(rhor_fftw)

  !
  ! Watch the reversed ordering nfft3, nfft2, nfft1
  plan = fftw_plan_dft_r2c_3d(nfft3, nfft2, nfft1, rhor_fftw, rhog_fftw, FFTW_ESTIMATE)
  call fftw_execute_dft_r2c(plan, rhor_fftw, rhog_fftw)
  call fftw_destroy_plan(plan)

  rhog_fftw(:, :, :) = rhog_fftw(:, :, :) / REAL(nfft_global, KIND=dp)

  write(stdout, '(1x,a,f12.6)') 'Charge of the distribution (real space):       ', &
                                SUM(rhor_fftw) * volume / REAL(nfft_global, KIND=dp)
  write(stdout, '(1x,a,f12.6)') 'Charge of the distribution (reciprocal space): ', &
                                rhog_fftw(1, 1, 1)%re * volume


  pvg = fftw_alloc_complex(int( (nfft1/2+1) * nfft2 * nfft3, C_SIZE_T))
  call C_F_POINTER(pvg, vg_fftw, [(nfft1/2+1), nfft2, nfft3])

  do ig3=1, nfft3
    do ig2=1, nfft2
      do ig1=1, nfft1/2+1
        if( ig1 * ig2 * ig3 > 1 ) then
          vg_fftw(ig1, ig2, ig3) = vcoulg(ig1, ig2, ig3) * rhog_fftw(ig1, ig2, ig3)
        else
          vg_fftw(ig1, ig2, ig3) = 0.0_dp
        endif
      enddo
    enddo
  enddo

  pvr = fftw_alloc_real(INT(nfft_global, C_SIZE_T))
  call C_F_POINTER(pvr, vr_fftw, [nfft1, nfft2, nfft3])

  plan = fftw_plan_dft_c2r_3d(nfft3, nfft2, nfft1, vg_fftw, vr_fftw, FFTW_ESTIMATE)
  call fftw_execute_dft_c2r(plan, vg_fftw, vr_fftw)
  call fftw_destroy_plan(plan)


  !
  ! Keep the potential on the local section of the grid only
  ifft_local = 0
  ifft_global = 0
  rhor_fftw(:, :, :) = 0.0_dp
  do i3=1, nfft3
    do i2=1, nfft2
      do i1=1, nfft1
        ifft_global = ifft_global + 1
        if( MODULO(ifft_global - 1, grid%nproc) == grid%rank ) then
          ifft_local = ifft_local + 1
          vcoulr(ifft_local) = vr_fftw(i1, i2, i3)
        endif
      enddo
    enddo
  enddo

#else
  vcoulr(:) = rhor(:) ! Fake operation to cheat on the unused dummy variable check of the compiler:w
  call die('poisson_solver_fft: requires to be compiled FFTs (-DHAVE_FFTW3)')
#endif


contains

pure function vcoulg(ig1, ig2, ig3) 
  implicit none
  integer, intent(in) :: ig1, ig2, ig3
  real(dp) :: vcoulg
  !=====
  integer  :: ig1m, ig2m, ig3m
  real(dp) :: g2
  !=====
  
  ig1m = MERGE( ig1 - 1, ig1 - 1 - nfft1, ig1 <= nfft1/2)
  ig2m = MERGE( ig2 - 1, ig2 - 1 - nfft2, ig2 <= nfft2/2)
  ig3m = MERGE( ig3 - 1, ig3 - 1 - nfft3, ig3 <= nfft3/2)
  g2 = SUM( MATMUL(bprim(:, :) , [ig1m, ig2m ,ig3m])**2 )

  vcoulg = 4.0_dp * pi / g2  !* (1.0 -  COS( SQRT(g2) * 37.0 ))

end function vcoulg

end subroutine poisson_solver_fft


!=========================================================================
subroutine calculate_basis_functions_periodic(basis)
  implicit none

  type(basis_set), intent(in) :: basis
  !=====
  integer               :: gt
  integer               :: ifft
  integer               :: ishell, ibf1, ibf2, ibf1_cart
  integer               :: i_cart
  integer               :: ni_cart, li
  real(dp), allocatable :: basis_function_r_cart(:, :), dr(:, :)
  !=====

  allocate(dr(3, nfft_local))

  gt = get_gaussian_type_tag(basis%gaussian_type)

  !$OMP PARALLEL PRIVATE(li, ni_cart, ibf1, ibf1_cart, ibf2, basis_function_r_cart, dr)
  !$OMP DO
  do ishell=1, basis%nshell
    li      = basis%shell(ishell)%am
    ni_cart = number_basis_function_am('CART',li)
    ibf1      = basis%shell(ishell)%istart
    ibf1_cart = basis%shell(ishell)%istart_cart
    ibf2      = basis%shell(ishell)%iend

    ! relative position to the shell center
    do concurrent(ifft=1:nfft_local)
      dr(:, ifft) = rgrid(:, ifft) - basis%shell(ishell)%x0(:)
    enddo
    ! transform to reduced coordinates
    dr(:, :) = MATMUL( aprim_inv, dr)
    ! minimum image convention
    do concurrent(ifft=1:nfft_local)
      dr(:, ifft) = dr(:, ifft) - NINT(dr(:, ifft))
    enddo
    ! transform back to cartesian coordinates
    dr(:, :) = MATMUL( aprim, dr) 
    ! transform balck to absolute position
    do concurrent(ifft=1:nfft_local)
      dr(:, ifft) = dr(:, ifft) + basis%shell(ishell)%x0(:)
    enddo

    allocate(basis_function_r_cart(ni_cart, nfft_local))

    do ifft=1, nfft_local
      do i_cart=1, ni_cart
        basis_function_r_cart(i_cart, ifft) = eval_basis_function(basis%bfc(ibf1_cart + i_cart - 1), dr(:, ifft))
      enddo
    enddo

    bfr(ibf1:ibf2, :) = MATMUL(  TRANSPOSE(cart_to_pure(li, gt)%matrix(:, :)) , basis_function_r_cart(:, :) )
    deallocate(basis_function_r_cart)

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(dr)

end subroutine calculate_basis_functions_periodic


!=========================================================================
subroutine destroy_fft_grid()
  implicit none

  !=====
  !=====

  deallocate(rgrid)
  deallocate(rhoelecr)
  deallocate(bfr)

end subroutine destroy_fft_grid

end module m_hamiltonian_periodic


!=========================================================================
