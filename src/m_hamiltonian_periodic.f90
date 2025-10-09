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

  integer, private, parameter :: nx=1

  integer(C_INT), parameter :: nfft1=60
  integer(C_INT), parameter :: nfft2=60
  integer(C_INT), parameter :: nfft3=60
  integer :: nfft = nfft1 * nfft2 * nfft3

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

      deallocate(array_cart, matrix)

    enddo
    deallocate(alphaB, cB)
  enddo


  call stop_clock(MERGE(0, timing_overlap, in_rt_tddft))


end subroutine setup_kinetic_onecell


!=========================================================================
subroutine setup_hartree_periodic(basis, p_matrix, hartree_ao, ehartree)
  implicit none

  type(basis_set), intent(in) :: basis
  class(*), intent(in)  :: p_matrix(:, :, :)
  real(dp), intent(out) :: hartree_ao(:, :)
  real(dp), intent(out) :: ehartree
  !=====
  integer              :: nbf
  integer              :: ibf, kbf, lbf
  integer              :: ipair
  real(dp), allocatable :: x_vector(:)
  integer              :: timing_xxdft_hartree
  real(dp), allocatable :: pmat(:)
  real(dp), allocatable :: hnucl_ao(:, :)

  integer :: i1, i2, i3, ir
  real(dp) :: rr(3)
  real(dp) :: rgrid(3, nfft), rhogrid(nfft, nspin), rhonuclr(nfft)
  real(dp) :: vhartreegrid(nfft, nspin), vnuclgrid(nfft)
  real(dp) :: dr(3, 3)
  real(dp), allocatable :: c_matrix(:, :, :), occupation(:, :)
  !=====

#if !defined(HAVE_FFTW3)
  call die('setup_hartree_periodic: requires FFTW at compilation time')
#else

  nbf = SIZE(hartree_ao(:, :), DIM=1)

  select type(p_matrix)
  type is(real(dp))
    timing_xxdft_hartree   = timing_hartree
  type is(complex(dp))
    timing_xxdft_hartree   = timing_tddft_hartree
  class default
    call die("setup_hartree_periodic: p_matrix is neither real nor complex")
  end select

  call start_clock(timing_xxdft_hartree)

  write(stdout, *) 'Calculate periodic Hartree term'

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


  select type(p_matrix)
  type is(real(dp))
    call get_c_matrix_from_p_matrix(p_matrix, c_matrix, occupation)
  end select
  call start_clock(timing_tmp1)
  call calculate_density_periodic(basis, c_matrix, occupation, rgrid, rhogrid)
  call stop_clock(timing_tmp1)
  call output_timing_line('tmp1' , timing_tmp1 , 0)

  
  ! FBFB testing rhonuclr
  call prepare_nuclei_density_periodic(rgrid, rhonuclr)

  dr(:, 1) = aprim(:, 1) / nfft1
  dr(:, 2) = aprim(:, 2) / nfft2
  dr(:, 3) = aprim(:, 3) / nfft3
  call write_cube_file('toto.cube', nfft1, nfft2, nfft3, dr, rhonuclr, comment='test')

  
  do ir=1, nfft 
    write(1001, '(1x,i8,1x,3(f12.6,1x),4x,f16.6)') ir, rgrid(:, ir), rhogrid(ir, 1)
  enddo

  !
  ! Get the electrostatic potential of the electrons (Hartree)
  !
  call poisson_solver_fft(rhogrid(:, 1), vhartreegrid(:, 1)) 

  write(stdout,*) 'FBFB Eh', 0.5_dp * SUM( rhogrid * vhartreegrid ) * volume / REAL(nfft, KIND=dp)

  call calculate_hao_periodic(basis, rgrid, vhartreegrid, hartree_ao)
  call dump_out_matrix(.TRUE., '=== Hartree contribution (FFTW) ===', hartree_ao)

  !
  ! Get the electrostatic potential of the nuclei
  !
  call poisson_solver_fft(rhonuclr, vnuclgrid) 

  write(stdout,*) 'FBFB Enucl', SUM( rhogrid(:, 1) * vnuclgrid ) * volume / REAL(nfft, KIND=dp)
  allocate(hnucl_ao, MOLD=hartree_ao)
  call calculate_hao_periodic(basis, rgrid, RESHAPE(vnuclgrid, [nfft, 1]), hnucl_ao)
  call dump_out_matrix(.TRUE., '=== Nuclei contribution (FFTW) ===', hnucl_ao)

  deallocate(hnucl_ao)


  hartree_ao(:, :) = 0.0_dp

  allocate(pmat(npair))
  allocate(x_vector(nauxil_local))

  select type(p_matrix)
  type is(real(dp))
    do ipair=1, npair
      kbf = index_basis(1, ipair)
      lbf = index_basis(2, ipair)
      pmat(ipair) = SUM(p_matrix(kbf, lbf, :)) * 2.0_dp
    enddo
  type is(complex(dp))
    do ipair=1, npair
      kbf = index_basis(1, ipair)
      lbf = index_basis(2, ipair)
      ! As all pairs contribute twice for (k, l) and (l, k) and as P is hermitian,
      ! only the real part survives
      pmat(ipair) = SUM(p_matrix(kbf, lbf, :)%re) * 2.0_dp
    enddo
  end select

  ! X_P = \sum_{\alpha \beta} P_{\alpha \beta} * ( \alpha \beta | P )
  call DGEMV('T', npair, nauxil_local, 1.0d0, eri_3center, npair, pmat, 1, 0.0d0, x_vector, 1)
  ! v_H_{alpha beta} = \sum_P ( alpha beta | P ) * X_P
  call DGEMV('N', npair, nauxil_local, 1.0d0, eri_3center, npair, x_vector, 1, 0.0d0, pmat, 1)

  !$OMP PARALLEL PRIVATE(kbf, lbf)
  !$OMP DO
  do ipair=1, npair
    kbf = index_basis(1, ipair)
    lbf = index_basis(2, ipair)
    hartree_ao(kbf, lbf) = pmat(ipair)
    hartree_ao(lbf, kbf) = pmat(ipair)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! Do not forget that the eri_3center(ibf, ibf | P ) included a factor 0.50
  do ibf=1, nbf
    hartree_ao(ibf, ibf) = hartree_ao(ibf, ibf) * 2.0_dp
  enddo
  deallocate(x_vector, pmat)

  !
  ! Sum up the different contribution from different procs
  call world%sum(hartree_ao)
  hartree_ao(:, :) = hartree_ao(:, :) / REAL(poorman%nproc, dp)

  call dump_out_matrix(.TRUE., '=== Hartree contribution (ref) ===', hartree_ao)

  select type(p_matrix)
  type is(real(dp))
    ehartree = 0.5_dp * SUM( hartree_ao(:, :) * SUM(p_matrix(:, :, :), DIM=3) )
  type is(complex(dp))
    ehartree = 0.5_dp * SUM( hartree_ao(:, :) * SUM(REAL(p_matrix(:, :, :), dp), DIM=3) )
  end select

  write(stdout,*) 'Eh 2 gaussian', ehartree

  call stop_clock(timing_xxdft_hartree)
#endif

  call die("STOP: ENOUGH")

end subroutine setup_hartree_periodic


!=========================================================================
subroutine calculate_density_periodic(basis, c_matrix, occupation, rr, rhor)
  implicit none

  type(basis_set), intent(in) :: basis
  real(dp), intent(in)  :: c_matrix(:, :, :)
  real(dp), intent(in)  :: occupation(:, :)
  real(dp), intent(in)  :: rr(:, :)
  real(dp), intent(out) :: rhor(:, :)
  !=====
  integer :: ispin, ir, nr, nstate, i1, i2, i3
  real(dp) :: shift(3)
  real(dp), allocatable :: rr_shifted(:, :)
  real(dp), allocatable :: bfr(:, :)
  real(dp), allocatable :: phir(:, :)
  !=====

  nr = SIZE(rr, DIM=2)
  nstate = SIZE(c_matrix, DIM=2)

  allocate(rr_shifted, MOLD=rr)
  allocate(bfr(basis%nbf,nr))
  allocate(phir(nr, nstate))

  rhor(:, :) = 0.0_dp
  do i3=-nx, nx
    do i2=-nx ,nx
      do i1=-nx, nx
        shift(:) = MATMUL( aprim(:, :), [i1, i2, i3] )

        do ir=1, nr
          rr_shifted(:, ir) = rr(:, ir) - shift(:)
          call calculate_basis_functions_r(basis, rr_shifted(:, ir), bfr(:, ir))
        enddo

        do ispin=1, nspin
          !phir(:, :) = MATMUL( TRANPOSE(bfr(:, :)) , c_matrix(:, :, ispin) )
          call DGEMM('T', 'N', nr, nstate, basis%nbf, 1.0d0, bfr, basis%nbf, c_matrix(:,:,ispin), basis%nbf, 0.0d0, phir, nr)

          rhor(:, :) = rhor(:, :) + MATMUL( phir(:, :)**2, occupation(:, :) )
        enddo

      enddo
    enddo
  enddo

  deallocate(rr_shifted)
  deallocate(bfr, phir)

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
  allocate(bfr(basis%nbf,nr))

  h_ao(:, :) = 0.0_dp
  do i3=-nx, nx
    do i2=-nx, nx
      do i1=-nx, nx
        shift(:) = MATMUL( aprim(:, :), [i1, i2, i3] )

        do ir=1, nr
          rr_shifted(:, ir) = rr(:, ir) - shift(:)
          call calculate_basis_functions_r(basis, rr_shifted(:, ir), bfr(:, ir))
        enddo

        do ispin=1, nspin
          do jbf=1, basis%nbf
            do ibf=1, basis%nbf
              h_ao(ibf, jbf) = h_ao(ibf, jbf) + SUM( bfr(ibf, :) * vr(:, ispin) * bfr(jbf, :) )
            enddo
          enddo
        enddo

      enddo
    enddo
  enddo
  h_ao(:, :) = h_ao(:, :) * volume / nr

  deallocate(rr_shifted)
  deallocate(bfr)

end subroutine calculate_hao_periodic


!=========================================================================
subroutine prepare_nuclei_density_periodic(rgrid, rhonuclr)
  implicit none

  real(dp), intent(in)  :: rgrid(:, :)
  real(dp), intent(out) :: rhonuclr(:)
  !=====
  integer :: nr, ir, irmin, irmax, ngrid
  real(dp), allocatable :: rrad(:), vloc(:), dvdr(:), d2vdr2(:), rho(:)
  real(dp), allocatable :: rradfinal(:), rhofinal(:)
  real(dp) :: zval, rmax, dr, shift(3)
  integer :: i1, i2, i3, icenter, igrid
  !=====

  ngrid = SIZE(rgrid, DIM=2)

  call read_potential_data('vloc.dat', nr, rrad, vloc)

  allocate(dvdr, MOLD=vloc)
  call calc_derivative(rrad, vloc, dvdr)
  allocate(d2vdr2, MOLD=vloc)
  call calc_derivative(rrad, dvdr, d2vdr2)


  allocate(rho, MOLD=vloc)
  do ir=1, nr
    rho(ir) = -( d2vdr2(ir) + 2.0_dp * dvdr(ir) / rrad(ir)) / (4.0_dp * pi)
    !write(1001,'(*(es20.10,1x))') rrad(ir), vloc(ir), dvdr(ir), d2vdr2(ir), rho(ir)
  enddo

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

  do ir=1, nr
    write(1001,'(*(es20.10,1x))') rradfinal(ir), rhofinal(ir)
  enddo
  zval = calculate_total_charge(rradfinal, rhofinal)
  write(stdout,'(1x,a,f12.8)') 'Integrated nucleus charge: ', zval
  rhofinal(:) = rhofinal(:) / zval * (-1.0_dp)
  zval = calculate_total_charge(rradfinal, rhofinal)
  write(stdout,'(1x,a,f12.8)') 'Integrated nucleus charge (after renormalization): ', zval



  !
  ! Loop over the unitcell regular grid
  !
  rhonuclr(:) = 0.0_dp
  do igrid=1, ngrid
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
            rhonuclr(igrid) = rhonuclr(igrid) + ( dr - rradfinal(ir)   ) / ( rradfinal(ir) - rradfinal(ir-1) ) * rhofinal(ir-1) &
                                              + ( rradfinal(ir-1) - dr ) / ( rradfinal(ir) - rradfinal(ir-1) ) * rhofinal(ir)

          enddo
        enddo
      enddo
    enddo
  enddo
 

  deallocate(dvdr, d2vdr2)

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




end module m_hamiltonian_periodic


!=========================================================================
