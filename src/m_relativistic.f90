!==================================================================
! This file is part of MOLGW.
! Author: M. Rodriguez-Mayorga
!
! This file contains
! - The construction of the Relativistic Hamiltonians (4C and X2C) 
!=========================================================================
#include "molgw.h"
module m_relativistic
  use m_definitions
  use m_warning
  use m_timing
  use m_mpi
  use m_elements
  use m_string_tools, only: orbital_momentum_name, append_to_list
  use m_atoms
  use m_ecp
  use m_gaussian
  use m_cart_to_pure
  use m_basis_set
  use m_eri_ao_mo
  use m_inputparam
  use m_hamiltonian_onebody
  use m_hamiltonian_tools
  use m_hamiltonian_wrapper
  use m_libcint_tools

  implicit none

  ! Tolerance on the (Jacobi scaled) small component overlaps: tighter than min_overlap used for the large component
  real(dp), parameter :: TOL_SMALL_OVERLAP = 1.0e-8_dp

 public::relativistic_init

 private::TOL_SMALL_OVERLAP, build_rkb_from_large, print_spinor_energies, time_reversal_deviation, &
   & hermitian_canonical_x, hermitian_inverse_sqrt, invert_general

contains


!=========================================================================
! Relativistic Hamiltonians (4C and X2C) from the large component basis only:
! the restricted-kinetic-balance small component sigma.p |large> needs only the integrals
! S, T, V and sigma.p V sigma.p over the (canonical, spherical) large component basis.
subroutine relativistic_init(is_x2c, electrons_in, nstate, c_matrix, s_matrix, x_matrix, &
    &   H_rel_rkb_mat, energy_rel)

  logical, intent(in)             :: is_x2c
  real(dp), intent(in)            :: electrons_in
  integer, intent(out)            :: nstate
  real(dp), allocatable, dimension(:), intent(inout)     ::energy_rel
  complex(dp), allocatable, dimension(:, :), intent(inout)::c_matrix
  complex(dp), allocatable, dimension(:, :), intent(inout)::s_matrix
  complex(dp), allocatable, dimension(:, :), intent(inout)::x_matrix
  complex(dp), allocatable, dimension(:, :), intent(inout)::H_rel_rkb_mat
  !====
  type(basis_set)                :: basis_nrel
  character(len=100), allocatable :: basis_name_nrel(:)
  integer                        :: gt, nbf_large, nbf_pure, nstate_large, nstate_c2, nstate_rkb
  integer                        :: ibf, jbf, ishell, icart, ipure, ncart_shell, npure_shell, am_shell, kcomp
  integer                        :: ielectrons, nkept
  real(dp)                       :: eext, max_dev
  real(dp), allocatable           :: W(:), W_x2c(:)
  real(dp), allocatable           :: s_matrix_large(:, :), x_matrix_large(:, :)
  real(dp), allocatable           :: t_sph(:, :), s_pure(:, :), x_cart(:, :)
  real(dp), allocatable           :: t_kin(:, :), v_nuc(:, :), w_sp(:, :, :)
  real(dp), allocatable           :: kin_c(:, :), v_c(:, :), w_c(:, :, :)
  complex(dp), allocatable        :: t_pure(:, :)           ! Canonical large spin-orbitals in terms of the SCF AO^L
  complex(dp), allocatable        :: s_nrel(:, :)           ! Spin-duplicated overlap of the SCF AO^L
  complex(dp), allocatable        :: h_rkb(:, :), h_rkb_ortho(:, :)
  complex(dp), allocatable        :: s_full(:, :), x_full(:, :)
  complex(dp), allocatable        :: u_mat(:, :), c_full(:, :)
  complex(dp), allocatable        :: c_l_inv(:, :), r_mat(:, :), lambda_mat(:, :), lambda_inv_sqrt(:, :)
  complex(dp), allocatable        :: h_x2c(:, :), h_x2c_ortho(:, :), tmp_matrix(:, :), s_t_large(:, :)
  !=====

  call timer_relativistic%start()


#if defined(HAVE_LIBCINT)

  write(stdout, '(/,a)') ' Relativitistic Hamiltonian construction'
  write(stdout, '(a,/)') ' ======================================='

  !! Large component basis: cartesian gaussians, transformed below to spherical (X2C is always spherical)
  write(stdout, '(/,a,/)') ' Computing large component reference basis'
  allocate(basis_name_nrel(ncenter_basis))
  basis_name_nrel(:) = basis_name(1)
  call init_basis_set(basis_path, basis_name_nrel, ecp_basis_name, 'CART', &
                        even_tempered_alpha, even_tempered_beta, even_tempered_n_list, basis_nrel)
  deallocate(basis_name_nrel)
  nbf_large=basis_nrel%nbf

  gt=get_gaussian_type_tag('PURE')
  nbf_pure=0
  do ishell=1, basis_nrel%nshell
    nbf_pure=nbf_pure+number_basis_function_am('PURE', basis_nrel%shell(ishell)%am)
  end do
  allocate(t_sph(nbf_large, nbf_pure))
  t_sph=0.0_dp
  icart=0; ipure=0;
  do ishell=1, basis_nrel%nshell
    am_shell=basis_nrel%shell(ishell)%am
    ncart_shell=number_basis_function_am('CART', am_shell)
    npure_shell=number_basis_function_am('PURE', am_shell)
    t_sph(icart+1:icart+ncart_shell, ipure+1:ipure+npure_shell)=cart_to_pure(am_shell, gt)%matrix(:, :)
    icart=icart+ncart_shell
    ipure=ipure+npure_shell
  end do

  !! Integrals over the cartesian large component basis
  call clean_allocate('Large overlap matrix S', s_matrix_large, nbf_large, nbf_large)
  allocate(t_kin(nbf_large, nbf_large), v_nuc(nbf_large, nbf_large))
  call init_libcint(basis_nrel)
  call setup_overlap(basis_nrel, s_matrix_large)
  call setup_kinetic(basis_nrel, t_kin)
  call setup_nucleus(basis_nrel, v_nuc)
  call setup_electric_field(basis_nrel, eext, v_nuc)
  call setup_spnucsp_cart(basis_nrel, w_sp)
  call destroy_libcint(basis_nrel)
  call destroy_basis_set(basis_nrel)

  !! Canonical large component basis: linearly dependent directions are dropped (nstate_large <= nbf_pure)
  write(stdout, '(a)') ' Building the canonical large component basis'
  allocate(s_pure(nbf_pure, nbf_pure))
  s_pure=matmul(transpose(t_sph), matmul(s_matrix_large, t_sph))
  call clean_deallocate('Large overlap matrix S', s_matrix_large)
  call setup_x_matrix(min_overlap, s_pure, nstate_large, x_matrix_large)
  nstate_c2=2*nstate_large
  nstate_rkb=2*nstate_c2
  write(stdout, '(a,i10)') ' Canonical large spin-orbitals', nstate_c2

  allocate(x_cart(nbf_large, nstate_large))
  x_cart=matmul(t_sph, x_matrix_large)
  allocate(kin_c(nstate_large, nstate_large), v_c(nstate_large, nstate_large), w_c(nstate_large, nstate_large, 4))
  kin_c=matmul(transpose(x_cart), matmul(t_kin, x_cart))
  v_c=matmul(transpose(x_cart), matmul(v_nuc, x_cart))
  do kcomp=1, 4
    w_c(:, :, kcomp)=matmul(transpose(x_cart), matmul(w_sp(:, :, kcomp), x_cart))
  end do
  deallocate(t_kin, v_nuc, w_sp, x_cart, t_sph)

  ! Canonical large spin-orbitals in terms of the spherical AO^L, and the spin-duplicated spherical AO^L overlap
  allocate(t_pure(2*nbf_pure, nstate_c2), s_nrel(2*nbf_pure, 2*nbf_pure))
  t_pure=complex_zero; s_nrel=complex_zero;
  do ibf=1, nbf_pure
    do jbf=1, nbf_pure
      s_nrel(2*ibf-1, 2*jbf-1)=s_pure(ibf, jbf)
      s_nrel(2*ibf, 2*jbf)=s_pure(ibf, jbf)
    end do
    do jbf=1, nstate_large
      t_pure(2*ibf-1, 2*jbf-1)=x_matrix_large(ibf, jbf)
      t_pure(2*ibf, 2*jbf)=x_matrix_large(ibf, jbf)
    end do
  end do
  deallocate(s_pure)
  call clean_deallocate('Overlap X * X**H = S**-1', x_matrix_large)

  !! RKB Hamiltonian and metric: H_RKB, S_full = diag(1, 2T) and X_full = diag(1, (2T)^-1/2)
  write(stdout, '(/,a)') ' Building the H_rel in RKB from the large component integrals'
  write(stdout, '(a,i10)') ' RKB Nbasis ', nstate_rkb
  call build_rkb_from_large(nstate_large, kin_c, v_c, w_c, h_rkb, s_full, x_full)
  deallocate(kin_c, v_c, w_c)
  allocate(tmp_matrix(nstate_rkb, nstate_rkb))
  tmp_matrix=matmul(conjg(transpose(x_full)), matmul(s_full, x_full))
  do ibf=1, nstate_rkb
    tmp_matrix(ibf, ibf)=tmp_matrix(ibf, ibf)-(1.0_dp, 0.0_dp)
  end do
  max_dev=maxval(abs(tmp_matrix))
  write(stdout, '(a,es12.4)') ' Max |X^dagger S X - I| in RKB ', max_dev
  deallocate(tmp_matrix)
  write(stdout, '(a,/)') ' Completed the H_rel in RKB'

  !! Diagonalize the H_rel in RKB: H^RKB C = S C e
  write(stdout, '(/,a)') ' Diagonalizing the H_rel in RKB'
  allocate(h_rkb_ortho(nstate_rkb, nstate_rkb), u_mat(nstate_rkb, nstate_rkb), c_full(nstate_rkb, nstate_rkb))
  allocate(W(nstate_rkb))
  h_rkb_ortho=matmul(conjg(transpose(x_full)), matmul(h_rkb, x_full))
  call diagonalize(' ', h_rkb_ortho, W, u_mat)
  deallocate(h_rkb_ortho)
  c_full=matmul(x_full, u_mat)
  deallocate(u_mat)
  ielectrons=nint(electrons_in)
  call print_spinor_energies(W, nstate_c2/2)
  write(stdout, '(a,/)') ' Diagonalized the H_rel in RKB'
  write(stdout, '(/,/,a25,1x,f19.10,/)') 'Rel Hcore Energy (Ha):', sum(W(nstate_c2+1:nstate_c2+ielectrons))

  if(.not.is_x2c) then ! 4c-calculations

   ! Set the crucial nstate
   nstate=nstate_rkb

   ! TODO these arrays should be passed to molgw.f90 for 4C calcs.
   call clean_allocate('Full RKB wavefunctions C', c_matrix, nstate_rkb, nstate_rkb)
   call clean_allocate('Full RKB S matrix', s_matrix, nstate_rkb, nstate_rkb)
   call clean_allocate('Full RKB X matrix', x_matrix, nstate_rkb, nstate_rkb)
   call clean_allocate('H_rel in RKB', H_rel_rkb_mat, nstate_rkb, nstate_rkb)
   c_matrix(:, :)=c_full(:, :)
   s_matrix(:, :)=s_full(:, :)
   x_matrix(:, :)=x_full(:, :)
   H_rel_rkb_mat(:, :)=h_rkb(:, :)
   deallocate(c_full, s_full, x_full, h_rkb, W)

   write(stdout, '(/,a)') ' Completed Relativistic Hamiltonian construction'
   write(stdout, '(a,/)') ' ==============================================='

  else ! X2C

   !! Decoupling R = C_S(+) C_L(+)^-1 from the (+) energy states (columns nstate_c2+1:nstate_rkb)
   write(stdout, '(/,a)') ' Computing X2C decoupling matrix R'
   allocate(c_l_inv(nstate_c2, nstate_c2), r_mat(nstate_c2, nstate_c2), tmp_matrix(nstate_c2, nstate_c2))
   c_l_inv(:, :)=c_full(1:nstate_c2, nstate_c2+1:)
   call invert_general(nstate_c2, c_l_inv)
   r_mat=matmul(c_full(nstate_c2+1:, nstate_c2+1:), c_l_inv)
   deallocate(c_l_inv)
   ! Check R C_L(+) = C_S(+)
   tmp_matrix=matmul(r_mat, c_full(1:nstate_c2, nstate_c2+1:))-c_full(nstate_c2+1:, nstate_c2+1:)
   max_dev=maxval(abs(tmp_matrix))
   write(stdout, '(a,es12.4)') ' Max |R C_L - C_S| ', max_dev
   deallocate(c_full)

   !! Lambda = S_LL + R^dagger S_SS R and H^X2C = H_LL + H_LS R + (H_LS R)^dagger + R^dagger H_SS R
   write(stdout, '(a)') ' Building H^X2C'
   allocate(lambda_mat(nstate_c2, nstate_c2), h_x2c(nstate_c2, nstate_c2))
   lambda_mat=s_full(1:nstate_c2, 1:nstate_c2) &
     &   + matmul(conjg(transpose(r_mat)), matmul(s_full(nstate_c2+1:, nstate_c2+1:), r_mat))
   tmp_matrix=matmul(h_rkb(1:nstate_c2, nstate_c2+1:), r_mat)
   h_x2c=h_rkb(1:nstate_c2, 1:nstate_c2) + tmp_matrix + conjg(transpose(tmp_matrix)) &
     &   + matmul(conjg(transpose(r_mat)), matmul(h_rkb(nstate_c2+1:, nstate_c2+1:), r_mat))
   tmp_matrix=h_x2c-conjg(transpose(h_x2c))
   max_dev=maxval(abs(tmp_matrix))
   write(stdout, '(a,es12.4)') ' Max |H^X2C - (H^X2C)^dagger| ', max_dev
   max_dev=time_reversal_deviation(h_x2c)
   write(stdout, '(a,es12.4)') ' Max time-reversal violation of H^X2C ', max_dev
   if( max_dev > 1.0e-8_dp ) call issue_warning('Relativistic: H^X2C does not preserve time-reversal symmetry')
   deallocate(tmp_matrix, r_mat, s_full, x_full, h_rkb)

   !! Renormalization R^dagger Lambda R = S_LL = 1 (canonical basis), R = Lambda^-1/2: H^FW = R^dagger H^X2C R
   !! Its spectrum must reproduce the (+) energy states
   write(stdout, '(a)') ' Renormalizing H^X2C with Lambda^-1/2 and diagonalizing H^FW'
   call hermitian_inverse_sqrt(lambda_mat, lambda_inv_sqrt)
   deallocate(lambda_mat)
   allocate(h_x2c_ortho(nstate_c2, nstate_c2), u_mat(nstate_c2, nstate_c2), W_x2c(nstate_c2))
   h_x2c_ortho=matmul(conjg(transpose(lambda_inv_sqrt)), matmul(h_x2c, lambda_inv_sqrt))
   call diagonalize(' ', h_x2c_ortho, W_x2c, u_mat)
   h_x2c(:, :)=h_x2c_ortho(:, :)
   deallocate(h_x2c_ortho, u_mat, lambda_inv_sqrt)
   max_dev=maxval(abs(W_x2c(:)-W(nstate_c2+1:nstate_rkb)))
   write(stdout, '(a,es12.4)') ' Max |E(X2C) - E(RKB, +)| ', max_dev
   if( max_dev > 1.0e-8_dp ) call issue_warning('Relativistic: X2C energies differ from the (+) energies of the 4C Hamiltonian')
   allocate(energy_rel(nstate_c2))
   energy_rel(:)=W_x2c(:)
   deallocate(W_x2c, W)

   !! Back to the AO^L representation: H_AO = S T H^X2C T^dagger S, so that T^dagger H_AO T = H^X2C
   write(stdout, '(a)') ' Transforming H^FW to the AO^L basis'
   allocate(s_t_large(2*nbf_pure, nstate_c2))
   s_t_large=matmul(s_nrel, t_pure)
   call clean_allocate('H_rel in RKB', H_rel_rkb_mat, 2*nbf_pure, 2*nbf_pure)
   H_rel_rkb_mat=matmul(s_t_large, matmul(h_x2c, conjg(transpose(s_t_large))))
   deallocate(s_t_large)
   call clean_allocate('Full RKB X matrix', x_matrix, 2*nbf_pure, nstate_c2)
   call clean_allocate('Full RKB S matrix', s_matrix, 2*nbf_pure, 2*nbf_pure)
   x_matrix(:, :)=t_pure(:, :)
   s_matrix(:, :)=s_nrel(:, :)

   !! Initial guess: H^FW orthonormalized with the large component only (T^dagger S T = I)
   allocate(u_mat(nstate_c2, nstate_c2), W_x2c(nstate_c2))
   call diagonalize(' ', h_x2c, W_x2c, u_mat)
   call clean_allocate('Full RKB wavefunctions C', c_matrix, 2*nbf_pure, nstate_c2)
   c_matrix=matmul(t_pure, u_mat)
   write(stdout, '(a)') ' Energies of the H^FW orthonormalized with the large component only'
   call print_spinor_energies(W_x2c, ielectrons/2)
   deallocate(u_mat, W_x2c, h_x2c)

   ! Set the crucial nstate
   nstate=nstate_c2

   write(stdout, '(/,a)') ' Completed X2C Hamiltonian construction'
   write(stdout, '(a,/)') ' ======================================'

  end if

  deallocate(t_pure, s_nrel)

#endif

  call timer_relativistic%stop()

end subroutine relativistic_init


!==================================================================
! RKB Dirac Hamiltonian from the large component integrals only, in the basis ordered as
! (large, small) with (alpha, beta) pairs: the small functions are sigma.p |large> so that
!   H_LL = V,  H_LS = c sigma.p sigma.p = 2 c T,  H_SS = sigma.p (V - 2c^2) sigma.p,  S_SS = 2 T
! with sigma.p V sigma.p = W_4 + i sigma.W_(1:3) (libcint components), the rest energy c^2 being removed.
subroutine build_rkb_from_large(nc, kin_c, v_c, w_c, h_rkb, s_full, x_full)
  integer, intent(in)   :: nc
  real(dp), intent(in)  :: kin_c(nc, nc), v_c(nc, nc), w_c(nc, nc, 4)
  complex(dp), allocatable, intent(out) :: h_rkb(:, :), s_full(:, :), x_full(:, :)
  !=====
  integer                 :: nkept, imo, jmo, nlarge
  complex(dp), allocatable :: x_small(:, :)
  real(dp)                :: w_scalar
  !=====

  nlarge=2*nc
  allocate(h_rkb(4*nc, 4*nc), s_full(4*nc, 4*nc), x_full(4*nc, 4*nc))
  h_rkb=(0.0_dp, 0.0_dp); s_full=(0.0_dp, 0.0_dp); x_full=(0.0_dp, 0.0_dp)
  do jmo=1, nc
    do imo=1, nc
      ! LL
      h_rkb(2*imo-1, 2*jmo-1)=v_c(imo, jmo)
      h_rkb(2*imo  , 2*jmo  )=v_c(imo, jmo)
      ! LS and SL
      h_rkb(2*imo-1, nlarge+2*jmo-1)=2.0_dp * c_speedlight * kin_c(imo, jmo)
      h_rkb(2*imo  , nlarge+2*jmo  )=2.0_dp * c_speedlight * kin_c(imo, jmo)
      h_rkb(nlarge+2*imo-1, 2*jmo-1)=2.0_dp * c_speedlight * kin_c(imo, jmo)
      h_rkb(nlarge+2*imo  , 2*jmo  )=2.0_dp * c_speedlight * kin_c(imo, jmo)
      ! SS
      w_scalar = w_c(imo, jmo, 4) - 4.0_dp * c_speedlight**2 * kin_c(imo, jmo)
      h_rkb(nlarge+2*imo-1, nlarge+2*jmo-1)=w_scalar + im * w_c(imo, jmo, 3)
      h_rkb(nlarge+2*imo  , nlarge+2*jmo  )=w_scalar - im * w_c(imo, jmo, 3)
      h_rkb(nlarge+2*imo-1, nlarge+2*jmo  )=im * w_c(imo, jmo, 1) + w_c(imo, jmo, 2)
      h_rkb(nlarge+2*imo  , nlarge+2*jmo-1)=im * w_c(imo, jmo, 1) - w_c(imo, jmo, 2)
      ! Metric
      s_full(nlarge+2*imo-1, nlarge+2*jmo-1)=2.0_dp * kin_c(imo, jmo)
      s_full(nlarge+2*imo  , nlarge+2*jmo  )=2.0_dp * kin_c(imo, jmo)
    end do
  end do
  do imo=1, nlarge
    s_full(imo, imo)=(1.0_dp, 0.0_dp)
    x_full(imo, imo)=(1.0_dp, 0.0_dp)
  end do
  call hermitian_canonical_x(TOL_SMALL_OVERLAP, s_full(nlarge+1:, nlarge+1:), x_small, nkept)
  write(stdout, '(a,es12.4)') ' Tolerance on RKB S_SS eigenvalues ', TOL_SMALL_OVERLAP
  if( nkept /= nlarge ) then
    write(stdout, '(/,a,i10,i10,/)') ' RKB S_SS kept and expected directions ', nkept, nlarge
    call die("Relativistic: the RKB small component overlap is linearly dependent")
  end if
  x_full(nlarge+1:, nlarge+1:)=x_small(:, :)
  deallocate(x_small)

end subroutine build_rkb_from_large


!==================================================================
subroutine print_spinor_energies(eigval, ipair_marker)
  real(dp), intent(in)  :: eigval(:)
  integer, intent(in)   :: ipair_marker
  !=====
  integer :: ipair
  !=====

  write(stdout, '(1x,a)') '=== Energies ==='
  write(stdout, '(a)') '   #                               (Ha)                       &
  &                         (eV)      '
  write(stdout, '(a)') '                         bar                       ubar       &
  &               bar                       ubar '
  do ipair=1, SIZE(eigval)/2
   write(stdout, '(1x,i5,2(2(1x,f25.5)))') ipair, eigval(2*ipair-1), eigval(2*ipair), &
     &  eigval(2*ipair-1)*Ha_eV, eigval(2*ipair)*Ha_eV
   if(ipair==ipair_marker) write(stdout, '(a)') '  --------------------------------------------------------------'
  end do

end subroutine print_spinor_energies


!==================================================================
! Max violation of time-reversal symmetry (T = -i sigma_y K) of a one-body matrix in a basis ordered
! as (n alpha, n beta) pairs: h_(n b, m b) = conj(h_(n a, m a)) and h_(n a, m b) = - conj(h_(n b, m a))
function time_reversal_deviation(hmat) result(dev)
  complex(dp), intent(in) :: hmat(:, :)
  real(dp)                :: dev
  !=====
  integer :: ipair, jpair
  !=====

  dev = 0.0_dp
  do jpair=1, SIZE(hmat, DIM=2)/2
    do ipair=1, SIZE(hmat, DIM=1)/2
      dev = MAX(dev, ABS(hmat(2*ipair, 2*jpair) - CONJG(hmat(2*ipair-1, 2*jpair-1))), &
                     ABS(hmat(2*ipair-1, 2*jpair) + CONJG(hmat(2*ipair, 2*jpair-1))))
    end do
  end do

end function time_reversal_deviation


!==================================================================
! Canonical orthonormalization of a Hermitian positive semidefinite matrix S (X^dagger S X = I).
! S is Jacobi scaled first, so that the tolerance acts on the linear dependence and not on the scale.
! The columns of X belonging to scaled eigenvalues below TOL_OVERLAP are dropped.
subroutine hermitian_canonical_x(tol_overlap, smat, xmat, nkept)
  real(dp), intent(in)                 :: tol_overlap
  complex(dp), intent(in)              :: smat(:, :)
  complex(dp), allocatable, intent(out) :: xmat(:, :)
  integer, intent(out)                 :: nkept
  !=====
  integer                 :: nmat, ibf, jbf, ikept
  real(dp), allocatable    :: scale_vec(:), eigval(:)
  complex(dp), allocatable :: scaled(:, :), eigvec(:, :)
  !=====

  nmat=SIZE(smat, DIM=1)
  allocate(scale_vec(nmat), eigval(nmat), scaled(nmat, nmat), eigvec(nmat, nmat))
  do ibf=1, nmat
    scale_vec(ibf)=sqrt(real(smat(ibf, ibf), dp))
  end do
  do jbf=1, nmat
    do ibf=1, nmat
      scaled(ibf, jbf)=smat(ibf, jbf)/(scale_vec(ibf)*scale_vec(jbf))
    end do
  end do
  call diagonalize(' ', scaled, eigval, eigvec)
  nkept=COUNT(eigval(:) > tol_overlap)
  write(stdout, '(a,es12.4)') ' Lowest scaled overlap eigenvalue ', MINVAL(eigval)
  allocate(xmat(nmat, nkept))
  ikept=0
  do jbf=1, nmat
    if( eigval(jbf) > tol_overlap ) then
      ikept=ikept+1
      do ibf=1, nmat
        xmat(ibf, ikept)=eigvec(ibf, jbf)/sqrt(eigval(jbf))/scale_vec(ibf)
      end do
    end if
  end do
  deallocate(scale_vec, eigval, scaled, eigvec)

end subroutine hermitian_canonical_x


!==================================================================
! X = S^-1/2 (symmetric) of a Hermitian positive definite matrix
subroutine hermitian_inverse_sqrt(smat, xmat)
  complex(dp), intent(in)              :: smat(:, :)
  complex(dp), allocatable, intent(out) :: xmat(:, :)
  !=====
  integer                 :: nmat, ibf
  real(dp), allocatable    :: eigval(:)
  complex(dp), allocatable :: eigvec(:, :), scaled(:, :)
  !=====

  nmat=SIZE(smat, DIM=1)
  allocate(eigval(nmat), eigvec(nmat, nmat), scaled(nmat, nmat), xmat(nmat, nmat))
  call diagonalize(' ', smat, eigval, eigvec)
  if( MINVAL(eigval) < 1.0e-10_dp ) call die('hermitian_inverse_sqrt: matrix is not positive definite')
  do ibf=1, nmat
    scaled(:, ibf)=eigvec(:, ibf)/sqrt(eigval(ibf))
  end do
  xmat=matmul(scaled, conjg(transpose(eigvec)))
  deallocate(eigval, eigvec, scaled)

end subroutine hermitian_inverse_sqrt


!==================================================================
! In-place inverse of a general complex matrix (LU)
subroutine invert_general(nmat, amat)
  integer, intent(in)         :: nmat
  complex(dp), intent(inout)  :: amat(:, :)
  !=====
  integer                 :: info, lwork
  integer, allocatable     :: ipiv(:)
  complex(dp)             :: work_query(1)
  complex(dp), allocatable :: work(:)
  !=====

  allocate(ipiv(nmat))
  call zgetrf(nmat, nmat, amat, nmat, ipiv, info)
  if(info/=0) call die("invert_general: error in zgetrf")
  lwork=-1
  call zgetri(nmat, amat, nmat, ipiv, work_query, lwork, info)
  lwork=nint(real(work_query(1)))
  allocate(work(lwork))
  call zgetri(nmat, amat, nmat, ipiv, work, lwork, info)
  if(info/=0) call die("invert_general: error in zgetri")
  deallocate(ipiv, work)

end subroutine invert_general

!==================================================================

end module m_relativistic
!==================================================================

