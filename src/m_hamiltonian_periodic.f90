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

  integer(C_INT), parameter :: nfft1 = 120
  integer(C_INT), parameter :: nfft2 = nfft1
  integer(C_INT), parameter :: nfft3 = nfft1
  integer, parameter :: nfft = nfft1 * nfft2 * nfft3

contains


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
subroutine setup_hartree_periodic(basis, p_matrix, h_ao, ehartree, enuc, nucleus_only)
  implicit none

  type(basis_set), intent(in) :: basis
  class(*), intent(in)  :: p_matrix(:, :, :)
  real(dp), intent(out) :: h_ao(:, :)
  real(dp), intent(out) :: ehartree, enuc
  logical, intent(in), optional :: nucleus_only
  !=====
  integer               :: timing_xxdft_hartree
  real(dp), allocatable :: hnucl_ao(:, :)
  real(dp), allocatable :: p_matrix_local(:, :, :)
  integer :: i1, i2, i3, ir
  real(dp) :: rr(3)
  real(dp) :: rgrid(3, nfft), rhonuclr(nfft)
  real(dp), allocatable :: rhoelecr(:, :)
  real(dp) :: vhartreegrid(nfft, nspin), vnuclgrid(nfft)
  real(dp) :: dr(3, 3)
  real(dp), allocatable :: c_matrix(:, :, :), occupation(:, :), s_matrix(:, :)
  logical :: nucleus_only_
  !=====

#if !defined(HAVE_FFTW3)
  call die('setup_hartree_periodic: requires FFTW at compilation time')
#else

  if( PRESENT(nucleus_only) ) then
    nucleus_only_ = nucleus_only
  else
    nucleus_only_ = .FALSE.
  endif


  select type(p_matrix)
  type is(real(dp))
    timing_xxdft_hartree   = timing_hartree
  type is(complex(dp))
    timing_xxdft_hartree   = timing_tddft_hartree
  class default
    call die("setup_hartree_periodic: p_matrix is neither real nor complex")
  end select

  call start_clock(timing_xxdft_hartree)

  if( nucleus_only_ ) then
    write(stdout, '(/,1x,a)') 'Calculate periodic nucleus term'
  else
    write(stdout, '(/,1x,a)') 'Calculate periodic Hartree + nucleus term'
  endif


  ! Create the real-space grid
  ir = 0
  do i3=0, nfft3 - 1
    do i2=0, nfft2 - 1
      do i1=0, nfft1 - 1
        ir = ir + 1
        rgrid(:, ir) = MATMUL( aprim(:, :), [DBLE(i1)/nfft1, DBLE(i2)/nfft2, DBLE(i3)/nfft3] )
      enddo
    enddo
  enddo


  allocate(p_matrix_local(basis%nbf, basis%nbf, nspin))
  select type(p_matrix)
  type is(real(dp))
    p_matrix_local(:, :, :) = p_matrix(:, :, :)
  end select


  if( .NOT. nucleus_only_ ) then
    call start_clock(timing_tmp1)

    allocate(rhoelecr(nfft, nspin))
    call calculate_density_periodic(basis, p_matrix_local, rgrid, rhoelecr)

    call stop_clock(timing_tmp1)
    call output_timing_line('set density on FFT grid' , timing_tmp1 , 1)

    
    !
    ! Get the electrostatic potential of the electrons (Hartree)
    !
    call start_clock(timing_tmp2)

    !write(*,*) 'FBFB charge', SUM(rhoelecr(:, 1))/nfft*volume
    !dr(:, 1) = aprim(:, 1) / nfft1
    !dr(:, 2) = aprim(:, 2) / nfft2
    !dr(:, 3) = aprim(:, 3) / nfft3
    !call write_cube_file('toto.cube', nfft1, nfft2, nfft3, dr, rhoelecr(:, 1), comment='test')

    call poisson_solver_fft(rhoelecr(:, 1), vhartreegrid(:, 1)) 
    call stop_clock(timing_tmp2)
    call output_timing_line('timing 2 FFTs' , timing_tmp2 , 1)

    ehartree = 0.5_dp * SUM( rhoelecr * vhartreegrid ) * volume / REAL(nfft, KIND=dp)
    write(stdout,'(1x,a,f16.8)') 'Hartree energy (Ha): ', ehartree

    call calculate_hao_periodic(basis, rgrid, vhartreegrid, h_ao)
    call dump_out_matrix(.TRUE., '=== Hartree contribution (FFTW) ===', h_ao)

  else
    ehartree = 0.0_dp
    h_ao(:, :) = 0.0_dp
  endif

  !
  ! Get the electrostatic potential of the nuclei
  !
  call prepare_nuclei_density_periodic(rgrid, rhonuclr)

  !dr(:, 1) = aprim(:, 1) / nfft1
  !dr(:, 2) = aprim(:, 2) / nfft2
  !dr(:, 3) = aprim(:, 3) / nfft3
  !call write_cube_file('toto.cube', nfft1, nfft2, nfft3, dr, rhonuclr, comment='test')
  call poisson_solver_fft(rhonuclr, vnuclgrid) 

  if( ALLOCATED(rhoelecr) ) then
    enuc = SUM( rhoelecr(:, 1) * vnuclgrid ) * volume / REAL(nfft, KIND=dp)
    write(stdout,'(1x,a,f16.8)') 'Enucl (Ha): ', enuc
  else
    enuc = 0.0_dp
  endif

  allocate(hnucl_ao, MOLD=h_ao)
  call calculate_hao_periodic(basis, rgrid, RESHAPE(vnuclgrid, [nfft, 1]), hnucl_ao)
  call dump_out_matrix(.TRUE., '=== Nuclei contribution (FFTW) ===', hnucl_ao)

  h_ao(:, :) = h_ao(:, :) + hnucl_ao(:, :)

  deallocate(hnucl_ao)
  if( ALLOCATED(rhoelecr) ) deallocate(rhoelecr)

  call stop_clock(timing_xxdft_hartree)

#endif


end subroutine setup_hartree_periodic


!=========================================================================
subroutine calculate_density_periodic(basis, p_matrix, rr, rhor)
  implicit none

  type(basis_set), intent(in) :: basis
  real(dp), intent(in)  :: p_matrix(:, :, :)
  real(dp), intent(in)  :: rr(:, :)
  real(dp), intent(out) :: rhor(:, :)
  !=====
  integer :: ispin, ir, nr, i1, i2, i3
  real(dp) :: shift(3)
  real(dp), allocatable :: rr_shifted(:, :)
  real(dp), allocatable :: bfr(:, :)
  real(dp), allocatable :: tmp(:, :)
  !=====

  nr = SIZE(rr, DIM=2)

  allocate(rr_shifted, MOLD=rr)
  allocate(bfr(basis%nbf, nr))
  allocate(tmp(basis%nbf, nr))

  rhor(:, :) = 0.0_dp
  call calculate_basis_functions_periodic(basis, rr, bfr)
  do ispin=1, nspin
    ! phi_α(r) * ( P_αβ * phi_β(r) )
    call DGEMM('N', 'N', basis%nbf, nr, basis%nbf, 1.0d0, p_matrix(:, :, ispin), basis%nbf, &
               bfr(:, :), basis%nbf, 0.0d0, tmp, basis%nbf)
  
    rhor(:, ispin) = rhor(:, ispin) + SUM( bfr(:, :) * tmp(:, :), DIM=1)
  enddo

  deallocate(rr_shifted)
  deallocate(bfr, tmp)

  ! write(stdout,*) 'sum_r ρ(r)',SUM(rhor)
  ! write(stdout,*) 'dv', volume /REAL(nfft, KIND=dp)
  ! write(stdout,*) 'FBFB total rho', SUM(rhor) * volume /REAL(nfft, KIND=dp)

end subroutine calculate_density_periodic


!=========================================================================
subroutine calculate_hao_periodic(basis, rr, vr, h_ao)
  implicit none

  type(basis_set), intent(in) :: basis
  real(dp), intent(in)  :: rr(:, :)
  real(dp), intent(in) :: vr(:, :)
  real(dp), intent(out) :: h_ao(:, :)
  !=====
  integer :: nx=1
  integer :: ispin, ir, nr, nstate, i1, i2, i3, ibf, jbf
  real(dp) :: shift(3)
  real(dp), allocatable :: rr_shifted(:, :)
  real(dp), allocatable :: bfr(:, :)
  !=====

  nr = SIZE(rr, DIM=2)

  allocate(rr_shifted, MOLD=rr)
  allocate(bfr(basis%nbf, nr))

  h_ao(:, :) = 0.0_dp
  call calculate_basis_functions_periodic(basis, rr, bfr)
  do ispin=1, nspin
    do jbf=1, basis%nbf
      do ibf=1, basis%nbf
        h_ao(ibf, jbf) = h_ao(ibf, jbf) + SUM( bfr(ibf, :) * vr(:, ispin) * bfr(jbf, :) )
      enddo
    enddo
  enddo

  h_ao(:, :) = h_ao(:, :) * volume / nr

  deallocate(rr_shifted)
  deallocate(bfr)

end subroutine calculate_hao_periodic


!=========================================================================
! Evaluate the nuclei charge density on the FFT grid
! that induces the periodic local potential
!
! Convention: nuclei have a negative charge (since electrons have positive one)
!
subroutine prepare_nuclei_density_periodic(rgrid, rhonuclr)
  implicit none

  real(dp), intent(in)  :: rgrid(:, :)
  real(dp), intent(out) :: rhonuclr(:)
  !=====
  integer :: nr, ir, irmin, irmax, nfft
  real(dp), allocatable :: rrad(:), vloc(:), dvdr(:), d2vdr2(:), rho(:)
  real(dp), allocatable :: rradfinal(:), rhofinal(:)
  real(dp) :: zval, rmax, dr, shift(3), factor
  integer :: i1, i2, i3, icenter, igrid
  !=====

  nfft = SIZE(rgrid, DIM=2)

  call read_potential_data('vloc.dat', nr, rrad, vloc)

  allocate(dvdr, MOLD=vloc)
  call calc_derivative(rrad, vloc, dvdr)
  allocate(d2vdr2, MOLD=vloc)
  call calc_derivative(rrad, dvdr, d2vdr2)


  allocate(rho, MOLD=vloc)
  do ir=1, nr
    rho(ir) = -( d2vdr2(ir) + 2.0_dp * dvdr(ir) / rrad(ir)) / (4.0_dp * pi)
  enddo
  deallocate(dvdr, d2vdr2)

  ! clean up rho
  ! remove the too low density values at large r
  ! remove the too low r

  do ir=nr, 1, -1
    if( ABS(rho(ir)) > 1.0e-6_dp ) exit
  enddo
  irmax = ir
  do ir=1, nr
    if( rrad(ir) > 1.0e-3_dp ) exit
  enddo
  irmin = ir

  nr = irmax-irmin+1
  allocate(rradfinal(nr))
  allocate(rhofinal(nr))
  rradfinal(1) = 0.0_dp
  rhofinal(1) = SUM(rho(1:irmin)) / REAL(irmin, KIND=dp)
  rradfinal(2:nr) = rrad(irmin+1:irmin+nr-1)
  rhofinal(2:nr) = rho(irmin+1:irmin+nr-1)
  rmax = rradfinal(nr)

  ! FBFB debug
  !do ir=1, nr
  !  write(1001,'(*(es20.10,1x))') rradfinal(ir), rhofinal(ir)
  !enddo
  zval = calculate_total_charge(rradfinal, rhofinal)
  write(stdout,'(1x,a,f12.8)') 'Integrated nucleus charge: ', zval
  rhofinal(:) = rhofinal(:) / zval * (-1.0_dp)
  zval = calculate_total_charge(rradfinal, rhofinal)
  write(stdout,'(1x,a,f12.8)') 'Integrated nucleus charge (after renormalization): ', zval


  !
  ! Loop over the unitcell regular grid
  !
  rhonuclr(:) = 0.0_dp
  do igrid=1, nfft
    do i3=-nx, nx
      do i2=-nx, nx
        do i1=-nx, nx
          shift(:) = MATMUL( aprim(:, :), [i1, i2, i3] )
          do icenter=1, ncenter_nuclei
            dr = NORM2( rgrid(:, igrid) - xatom(:, icenter) - shift(:) )
            if( dr > rmax ) cycle
            do ir=1, nr
              if( rradfinal(ir) > dr ) exit
            enddo
            rhonuclr(igrid) = rhonuclr(igrid) + ( rradfinal(ir) - dr   ) / ( rradfinal(ir) - rradfinal(ir-1) ) * rhofinal(ir-1) &
                                              + ( dr - rradfinal(ir-1) ) / ( rradfinal(ir) - rradfinal(ir-1) ) * rhofinal(ir)

          enddo
        enddo
      enddo
    enddo
  enddo
 
  zval = -SUM(rhonuclr(:)) * volume / REAL(nfft, KIND=dp)
  write(stdout, '(1x,a,f12.6)') 'Nuclei charge in the cell: ', zval
  write(stdout, '(1x,a,f12.6)') 'whereas it should be: ', SUM(zvalence(:))

  factor = SUM(zvalence(:)) / zval
  rhonuclr(:) = rhonuclr(:) * factor

  if( ABS( factor - 1.0_dp ) > 0.05_dp ) then
    call die('prepare_nuclei_density_periodic: wrong nucleus total charge. FFT grid is certainly too coarse')
  endif

  zval = -SUM(rhonuclr(:)) * volume / REAL(nfft, KIND=dp)
  write(stdout, '(1x,a,f12.6)') 'Nuclei charge in the cell (after renormalization): ', zval


contains

  subroutine read_potential_data(filename, nr, r, phi)
    implicit none
    character(len=*), intent(in) :: filename

    !=====
    integer, intent(out) :: nr
    real(dp), allocatable, intent(out) :: r(:), phi(:)
    integer :: iostat, ir
    integer :: file_unit
    !=====
    
    ! Count lines
    open(newunit=file_unit, file=filename, status='old', action='read')
    nr = 0
    do
        read(file_unit, *, iostat=iostat)
        if (iostat /= 0) exit
        nr = nr + 1
    enddo
    rewind(file_unit)
    
    ! Allocate and read
    allocate(r(nr), phi(nr))
    do ir = 1, nr
        read(file_unit, *) r(ir), phi(ir)
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
    integer :: ir, nr

    nr = SIZE(rho)
    ! Integrate 4πr²ρ(r) dr using trapezoidal rule
    total_charge = 0.0_dp
    do ir = 1, nr - 1
      total_charge = total_charge &
           + 4.0_dp * pi * 0.5_dp * ( r(ir + 1) - r(ir) )   &
               * ( r(ir)**2 * rho(ir) + r(ir + 1)**2 * rho(ir + 1) )
    enddo
  end function calculate_total_charge


end subroutine prepare_nuclei_density_periodic


!=========================================================================
! Poisson solver on a FFT grid 
! input:  ρ(r)
! output: vcoul(r)
subroutine poisson_solver_fft(rhor, vcoulr)
  implicit none

  real(dp), intent(in)  :: rhor(:)
  real(dp), intent(out) :: vcoulr(:)
  !=====
  type(C_PTR) :: plan
  type(C_PTR) :: pr, pg, pvr, pvg
  real(C_DOUBLE), pointer :: rhor_fftw(:, :, :), vr_fftw(:, :, :)
  complex(C_DOUBLE_COMPLEX), pointer ::  rhog_fftw(:, :, :), vg_fftw(:, :, :)
  integer :: ig1, ig2, ig3
  !=====

  ! Allocate arrays with FFTW functions
  pr = fftw_alloc_real(int(nfft1 * nfft2 * nfft3, C_SIZE_T))
  call C_F_POINTER(pr, rhor_fftw, [nfft1, nfft2, nfft3])
  pg = fftw_alloc_complex(int( (nfft1/2+1) * nfft2 * nfft3, C_SIZE_T))
  call C_F_POINTER(pg, rhog_fftw, [(nfft1/2+1), nfft2, nfft3])

  rhor_fftw(:, :, :) = RESHAPE(rhor(:), [nfft1, nfft2, nfft3])
  ! reversed order nfft3, nfft2, nfft1
  plan = fftw_plan_dft_r2c_3d(nfft3, nfft2, nfft1, rhor_fftw, rhog_fftw, FFTW_ESTIMATE)
  call fftw_execute_dft_r2c(plan, rhor_fftw, rhog_fftw)
  call fftw_destroy_plan(plan)

  rhog_fftw(:, :, :) = rhog_fftw(:, :, :) / REAL(nfft, KIND=dp)

  write(stdout, '(1x,a,f12.6)') 'Charge of the distribution (real space):       ', SUM(rhor_fftw) * volume / REAL(nfft, KIND=dp)
  write(stdout, '(1x,a,f12.6)') 'Charge of the distribution (reciprocal space): ', rhog_fftw(1, 1, 1)%re * volume


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

  pvr = fftw_alloc_real(int(nfft1 * nfft2 * nfft3, C_SIZE_T))
  call C_F_POINTER(pvr, vr_fftw, [nfft1, nfft2, nfft3])

  plan = fftw_plan_dft_c2r_3d(nfft3, nfft2, nfft1, vg_fftw, vr_fftw, FFTW_ESTIMATE)
  call fftw_execute_dft_c2r(plan, vg_fftw, vr_fftw)
  call fftw_destroy_plan(plan)

  vcoulr(:) = RESHAPE(vr_fftw, [nfft])

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
subroutine calculate_basis_functions_periodic(basis, rr, basis_function_r)
  implicit none

  type(basis_set), intent(in) :: basis
  real(dp), intent(in)        :: rr(:, :)
  real(dp), intent(out)       :: basis_function_r(:, :)
  !=====
  integer               :: gt
  integer               :: ir, nr
  integer               :: ishell, ibf1, ibf2, ibf1_cart
  integer               :: i_cart
  integer               :: ni_cart, li
  real(dp), allocatable :: basis_function_r_cart(:, :), dr(:, :)
  !=====

  nr = SIZE(basis_function_r(:, :), DIM=2)
  allocate(dr(3, nr))

  gt = get_gaussian_type_tag(basis%gaussian_type)

  !$OMP PARALLEL PRIVATE(li, ni_cart, ibf1, ibf1_cart, ibf2, basis_function_r_cart)
  !$OMP DO
  do ishell=1, basis%nshell
    li      = basis%shell(ishell)%am
    ni_cart = number_basis_function_am('CART',li)
    ibf1      = basis%shell(ishell)%istart
    ibf1_cart = basis%shell(ishell)%istart_cart
    ibf2      = basis%shell(ishell)%iend

    ! relative position to the shell center
    do concurrent(ir=1:nr)
      dr(:, ir) = rr(:, ir) - basis%shell(ishell)%x0(:)
    enddo
    ! transform to reduced coordinates
    dr(:, :) = MATMUL( aprim_inv, dr)
    ! minimum image convention
    do concurrent(ir=1:nr)
      dr(:, ir) = dr(:, ir) - NINT(dr(:, ir))
    enddo
    ! transform back to cartesian coordinates
    dr(:, :) = MATMUL( aprim, dr) 
    ! transform balck to absolute position
    do concurrent(ir=1:nr)
      dr(:, ir) = dr(:, ir) + basis%shell(ishell)%x0(:)
    enddo

    allocate(basis_function_r_cart(ni_cart, nr))

    do ir=1, nr
      do i_cart=1, ni_cart
        basis_function_r_cart(i_cart, ir) = eval_basis_function(basis%bfc(ibf1_cart + i_cart - 1), dr(:, ir))
      enddo
    enddo

    basis_function_r(ibf1:ibf2, :) = MATMUL(  TRANSPOSE(cart_to_pure(li, gt)%matrix(:, :)) , basis_function_r_cart(:, :) )
    deallocate(basis_function_r_cart)

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(dr)

end subroutine calculate_basis_functions_periodic


end module m_hamiltonian_periodic


!=========================================================================
