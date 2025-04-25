!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the spectral decomposition of the screened Coulomb interaction W
!
!=========================================================================
#include "molgw.h"
module m_spectral_function
  use m_definitions
  use m_mpi
  use m_timing
  use m_warning
  use m_memory
  use m_scalapack
  use m_inputparam
  use m_eri, only: iproc_ibf_auxil, ibf_auxil_l
  use m_eri_calculate, only: nauxil_global, nauxil_local
  use m_numerical_tools, only: coeffs_gausslegint
  use m_basis_set
  use m_eri_ao_mo

  integer, parameter :: NO_GRID        = 0
  integer, parameter :: IMAGINARY_QUAD = 1
  integer, parameter :: REAL_LINEAR    = 2
  integer, parameter :: STATIC         = 3
  integer, parameter :: MANUAL         = 4

  !
  ! General form of any spectral function
  ! z complex number
  ! i, j running on the basis set
  ! sf_ij(z) = \sum_n L_n(i) L_n(j) / ( z - w_n )
  !
  ! or on a grid of frequencies
  !

  type chi_type
    real(dp), allocatable :: eigvec(:, :)
    real(dp), allocatable :: eigval(:)
  contains
    procedure :: destroy => ct_destroy
  end type chi_type

  type spectral_function

    integer              :: type                   ! type encodes the storage method

    integer              :: nprodbasis_total       ! total over all procs
    integer              :: nprodbasis             ! for this proc

    !
    ! Information about all the poles
    !
    integer              :: npole_reso
    integer, allocatable  :: transition_table(:, :)  ! correspondance table from
                                                   ! transition index to state pair indexes

    real(dp), allocatable :: pole(:)
    real(dp), allocatable :: residue_left(:, :)       ! first index runs on n, second index on i

    !
    ! Static or Dynamic W might be stored directly in the auxiliary basis
    real(dp), allocatable    :: chi(:, :, :)
    integer                 :: mchi, nchi
    integer                 :: desc_chi(NDEL)
    integer                 :: nomega           ! Number of quadrature points
    real(dp), allocatable    :: weight_quad(:)   ! quadrature weight for numerical integration
    complex(dp), allocatable :: omega(:)         ! frequency grid (can real or imaginary)

    type(chi_type), allocatable :: vchiv_sqrt(:)

  contains

    generic :: evaluate => sf_evaluate_several_omegas
    generic :: evaluate => sf_evaluate_one_omega
    generic :: evaluate => sf_evaluate_one_real_omega
    procedure :: sf_evaluate_several_omegas
    procedure :: sf_evaluate_one_omega
    procedure :: sf_evaluate_one_real_omega

    procedure :: vsqrt_chi_vsqrt_rpa => sf_vsqrt_chi_vsqrt_rpa
    procedure :: interpolate_vsqrt_chi_vsqrt => sf_interpolate_vsqrt_chi_vsqrt
    procedure :: init    => sf_init
    procedure :: destroy => sf_destroy

  end type spectral_function

  !
  ! Highest occupied state
  integer, protected :: nhomo_W

  !
  ! Lowest unoccupied state
  integer, protected :: nlumo_W

  !
  ! frozen core approximation parameter
  integer, protected :: ncore_W

  !
  ! frozen virtual approximation parameter
  integer, protected :: nvirtual_W


contains


!=========================================================================
pure function index_prodstate(istate, jstate)
  implicit none
  integer, intent(in)  :: istate, jstate
  integer             :: index_prodstate
  !=====
  integer             :: imin, imax
  !=====

  ! Index (i,j) transformed into (I)
  ! with this ordering:
  !   (  1  2  4  7  ... )
  !   (  2  3  5  8  ... )
  !   (  4  5  6  9  ... )
  !   (  7  8  9 10  ... )
  !   ( ................ )

  imin = MIN(istate, jstate)
  imax = MAX(istate, jstate)

  index_prodstate = ( imax * (imax - 1) ) / 2 + imin


end function index_prodstate


!=========================================================================
subroutine sf_init(sf, nstate, occupation, nomega_in, grid_type, omega_max, verbose)
  implicit none
  class(spectral_function), intent(out)  :: sf
  integer, intent(in)                    :: nstate
  real(dp), intent(in)                   :: occupation(:, :)
  integer, intent(in)                    :: nomega_in
  integer, optional, intent(in)           :: grid_type
  real(dp), optional, intent(in)          :: omega_max
  logical, optional, intent(in)           :: verbose
  !=====
  integer                               :: grid_
  real(dp)                              :: omega_max_
  integer                               :: stdout_
  integer                               :: ijspin, istate, jstate, itrans
  integer                               :: iomega
  integer                               :: nlumo_W_spin(nspin)
  integer                               :: nhomo_W_spin(nspin)
  real(dp), parameter                    :: alpha=1.0_dp ! 0.50_dp
  real(dp), parameter                    :: beta=2.0_dp ! 6.0_dp
  real(dp), parameter                    :: omega_0=2.0_dp ! tweaked value from Arno Foerster
  real(dp), allocatable                  :: omega_quad(:)
  !=====

  if( nstate > SIZE( occupation(:, :) , DIM=1 ) ) then
    call die('sf_init: nstate is too large')
  endif
  if( PRESENT(grid_type) ) then
    grid_ = grid_type
  else
    grid_ = NO_GRID
  endif
  if( PRESENT(omega_max) ) then
    omega_max_ = omega_max
  else
    omega_max_ = 1.0_dp
  endif
  stdout_ = stdout
  if( PRESENT(verbose) ) then
    if( .NOT. verbose ) then
      open(newunit=stdout_, file='/dev/null')
    endif
  endif

  sf%type = grid_

  ncore_W      = ncorew
  nvirtual_W   = MIN(nvirtualw, nstate+1)

  if(is_frozencore) then
    if( ncore_W == 0) ncore_W = atoms_core_states()
  endif
  if( ncore_W > 0 ) then
    write(stdout_, '(1x,a,i4,2x,i4)') 'frozen core approximation in W switched on up to state = ', ncore_W
  endif

  if( nvirtual_W <= nstate ) then
    write(stdout_, '(1x,a,i4,2x,i4)') 'frozen virtual approximation in W switched on starting with state = ', nvirtual_W
  endif

  !
  ! Find the highest occupied state
  nhomo_W         = 0
  nhomo_W_spin(:) = 0
  do ijspin=1, nspin
    do istate=1, nstate
      if( occupation(istate, ijspin) / spin_fact < completely_empty ) cycle
      nhomo_W              = MAX(nhomo_W, istate)
      nhomo_W_spin(ijspin) = MAX(nhomo_W_spin(ijspin), istate)
    enddo
  enddo

  !
  ! Find the lowest occupied state
  nlumo_W         = 100000
  nlumo_W_spin(:) = 100000
  do ijspin=1, nspin
    do istate=1, nstate
      if( (spin_fact - occupation(istate, ijspin)) / spin_fact < completely_empty) cycle
      nlumo_W              = MIN(nlumo_W, istate)
      nlumo_W_spin(ijspin) = MIN(nlumo_W_spin(ijspin), istate)
    enddo
  enddo

  write(stdout_, '(/,1x,a)') 'Prepare a polarizability spectral function with'
  if( nspin == 1 ) then
    write(stdout_, '(30x,a,i8)') ' Occupied states: ', nhomo_W-ncore_W
    write(stdout_, '(30x,a,i8)') '  Virtual states: ', nvirtual_W-nlumo_W
    write(stdout_, '(30x,a,i8)') 'Transition space: ', (nvirtual_W-nlumo_W)*(nhomo_W-ncore_W)
  else
    write(stdout_, '(30x,a,i8,2x,i8)') ' Occupied states: ', nhomo_W_spin(:)-ncore_W
    write(stdout_, '(30x,a,i8,2x,i8)') '  Virtual states: ', nvirtual_W-nlumo_W_spin(:)
    write(stdout_, '(30x,a,i8)')       'Transition space: ', (nvirtual_W-nlumo_W_spin(1))*(nhomo_W_spin(1)-ncore_W) &
                                                         + (nvirtual_W-nlumo_W_spin(nspin))*(nhomo_W_spin(nspin)-ncore_W)
  endif

  !
  ! First, count the number of resonant transitions
  itrans=0
  do ijspin=1, nspin
    do istate=1, nstate
      do jstate=1, nstate
        if( skip_transition(jstate, istate, occupation(jstate, ijspin), occupation(istate, ijspin)) ) cycle
        if( occupation(jstate, ijspin) - occupation(istate, ijspin) > 0.0_dp ) cycle
        itrans = itrans + 1
      enddo
    enddo
  enddo

  sf%npole_reso = itrans
  allocate(sf%transition_table(3, sf%npole_reso))
  ! Set the transition_table
  itrans=0
  do ijspin=1, nspin
    do istate=1, nstate
      do jstate=1, nstate
        if( skip_transition(jstate, istate, occupation(jstate, ijspin), occupation(istate, ijspin)) ) cycle
        if( occupation(jstate, ijspin) - occupation(istate, ijspin) > 0.0_dp ) cycle
        itrans = itrans + 1
        ! Set the resonant transition table
        sf%transition_table(1, itrans) = istate
        sf%transition_table(2, itrans) = jstate
        sf%transition_table(3, itrans) = ijspin
      enddo
    enddo
  enddo

  if( has_auxil_basis ) then
    sf%nprodbasis_total = nauxil_global
  else
    sf%nprodbasis_total = index_prodstate(nvirtual_W-1, nvirtual_W-1) * nspin
  endif


  !
  ! Set the sampling points for Chi
  sf%nomega    = nomega_in

  select case(grid_)
  case(STATIC)
    if( nomega_in /= 1 ) call die('sf_init: static chi requested, number of frequencies should be 1')
    allocate(sf%weight_quad(sf%nomega))
    allocate(sf%omega(sf%nomega))
    sf%weight_quad(1) = 1.0_dp
    sf%omega(1)       = (0.0_dp, 0.0_dp)
  case(IMAGINARY_QUAD)
    if( nomega_in < 1 ) call die('sf_init: grid points is zero whereas a grid is requested')
    allocate(sf%weight_quad(sf%nomega))
    allocate(omega_quad(sf%nomega))
    allocate(sf%omega(sf%nomega))

    if( .TRUE. ) then
      call coeffs_gausslegint(0.0_dp, 1.0_dp, omega_quad, sf%weight_quad, sf%nomega)

      write(stdout_, '(/,1x,a)') 'Numerical integration on a grid along the imaginary axis'
      ! Variable change [0,1] -> [0,+\inf[
      write(stdout_, '(a)') '    #    Frequencies (eV)    Quadrature weights'
      do iomega=1, sf%nomega
        sf%weight_quad(iomega) = sf%weight_quad(iomega) / ( 2.0_dp**alpha - 1.0_dp ) * alpha &
                                * (1.0_dp -  omega_quad(iomega))**(-alpha-1.0_dp) * beta
        omega_quad(iomega)  =  1.0_dp / ( 2.0_dp**alpha - 1.0_dp ) &
                                  * ( 1.0_dp / (1.0_dp-omega_quad(iomega))**alpha - 1.0_dp ) * beta
        sf%omega(iomega)       =  omega_quad(iomega) * im
        write(stdout_, '(i5,2(2x,f14.6))') iomega, sf%omega(iomega)%im*Ha_eV, sf%weight_quad(iomega)
      enddo
    else
      ! Grid from Erhard, Goerling JChemPhys 157, 114105 (2022)
      call issue_warning("FBFB Goerling hack omega freqs")
      call coeffs_gausslegint(-1.0_dp, 1.0_dp, omega_quad, sf%weight_quad, sf%nomega)

      write(stdout_, '(/,1x,a)') 'Numerical integration on a grid along the imaginary axis'
      ! Variable change [-1,1] -> [0,+\inf[
      write(stdout_, '(a)') '    #    Frequencies (eV)    Quadrature weights'
      do iomega=1, sf%nomega
        sf%weight_quad(iomega) = sf%weight_quad(iomega) * 2.0_dp * omega_0 / ( 1.0_dp - omega_quad(iomega) )**2
        sf%omega(iomega)       = im * omega_0 * ( 1.0_dp + omega_quad(iomega) ) / ( 1.0_dp - omega_quad(iomega) )
        write(stdout_, '(i5,2(2x,f14.6))') iomega, sf%omega(iomega)%im*Ha_eV, sf%weight_quad(iomega)
      enddo
    endif
    deallocate(omega_quad)

  case(REAL_LINEAR)
    if( nomega_in < 1 ) call die('sf_init: grid points is zero whereas a grid is requested')
    allocate(sf%weight_quad(sf%nomega))
    allocate(sf%omega(sf%nomega))
    write(stdout_, '(a)') '    #    Frequencies (eV)   real    / imaginary '
    do iomega=1, nomega_in
      sf%omega(iomega) = REAL(iomega-1, dp)/REAL(nomega_in-1, dp) * omega_max_
      write(stdout_, '(i5,2(2x,f14.6))') iomega, sf%omega(iomega)%re*Ha_eV, sf%omega(iomega)%im*Ha_eV
    enddo
  case(MANUAL)
    if( nomega_in < 1 ) call die('sf_init: grid points is zero whereas a grid is requested')
    ! Just allocate and fill with default values
    allocate(sf%weight_quad(sf%nomega))
    allocate(sf%omega(sf%nomega))
    sf%weight_quad(1) = 1.0_dp
    sf%omega(1)       = (0.0_dp, 0.0_dp)
  case(NO_GRID)
    ! No grid, nothing to do
  case default
    call die("sf_init: grid choice not valid")
  end select


  if( PRESENT(verbose) ) then
    if( .NOT. verbose ) then
      close(stdout_)
    endif
  endif

end subroutine sf_init


!=========================================================================
subroutine allocate_spectral_function(nprodbasis, sf)
  implicit none
  integer, intent(in)                    :: nprodbasis
  type(spectral_function), intent(inout) :: sf
  !=====

  sf%nprodbasis = nprodbasis

  write(stdout, '(/,a,i8)') ' Spectral function initialized with Coulomb basis functions: ', sf%nprodbasis
  write(stdout, '(a,i8)')   ' Spectral function initialized with resonant poles         : ', sf%npole_reso

  allocate(sf%pole(sf%npole_reso))
  call clean_allocate('Left residue', sf%residue_left, sf%nprodbasis, sf%npole_reso)


end subroutine allocate_spectral_function


!=========================================================================
pure function skip_transition(ib1, ib2, occ1, occ2)
  implicit none
  logical             :: skip_transition
  integer, intent(in)  :: ib1, ib2
  real(dp), intent(in) :: occ1, occ2
  !=====

  skip_transition = .FALSE.

  !
  ! skip the core states if asked for a frozen-core calculation
  if( ib1 <= ncore_W .OR. ib2 <= ncore_W ) skip_transition = .TRUE.
  !
  ! skip the virtual states if asked for a frozen-virtual calculation
  if( ib1 >= nvirtual_W .OR. ib2 >= nvirtual_W ) skip_transition = .TRUE.

  if( occ1 < completely_empty             .AND. occ2 < completely_empty )             skip_transition = .TRUE.
  if( occ1 > spin_fact - completely_empty .AND. occ2 > spin_fact - completely_empty ) skip_transition = .TRUE.


end function skip_transition


!=========================================================================
subroutine sf_destroy(sf, verbose)
  implicit none
  class(spectral_function), intent(inout) :: sf
  logical, optional, intent(in)            :: verbose
  !=====
  logical :: verbose_
  !=====

  verbose_ = .TRUE.
  if( PRESENT(verbose) ) then
    verbose_ = verbose
  endif

  if(ALLOCATED(sf%transition_table)) deallocate(sf%transition_table)
  if(ALLOCATED(sf%pole))             deallocate(sf%pole)
  if(ALLOCATED(sf%residue_left)) then
    call clean_deallocate('Left residue', sf%residue_left, verbose=verbose_)
  endif
  if(ALLOCATED(sf%chi)) then
    call clean_deallocate('Chi', sf%chi, verbose=verbose_)
  endif
  if(ALLOCATED(sf%weight_quad)) deallocate(sf%weight_quad)
  if(ALLOCATED(sf%omega))       deallocate(sf%omega)
  if(ALLOCATED(sf%vchiv_sqrt))  deallocate(sf%vchiv_sqrt)

  if( verbose_ ) then
    write(stdout, '(/,a)') ' Spectral function destroyed'
  endif

end subroutine sf_destroy


!=========================================================================
subroutine write_spectral_function(sf)
  implicit none
  type(spectral_function), intent(in) :: sf
  !=====
  integer              :: wfile
#if defined(HAVE_MPI)
  integer              :: iprodbasis
  real(dp), allocatable :: buffer(:)
  integer              :: ierr
  integer              :: ibf_auxil
  integer(kind=MPI_OFFSET_KIND) :: disp
#endif
  !=====
  integer :: ipole
  !=====

  write(stdout, '(/,a,/)') ' Writing the spectral function on file: SCREENED_COULOMB'

  write(stdout, *) 'Number of poles written down:', sf%npole_reso


#if !defined(HAVE_MPI)
  if( is_iomaster ) then
    open(newunit=wfile, file='SCREENED_COULOMB', form='unformatted')

    write(wfile) calc_type%postscf_name
    write(wfile) sf%nprodbasis_total
    write(wfile) sf%npole_reso
    write(wfile) sf%pole(:)
    do ipole=1, sf%npole_reso
      write(wfile) sf%residue_left(:, ipole)
    enddo

    close(wfile)
  endif

#else

  call MPI_FILE_OPEN(MPI_COMM_WORLD, 'SCREENED_COULOMB', &
                    MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                    MPI_INFO_NULL, wfile, ierr)

  ! Only one proc has to write the poles
  disp  = 0
  if(is_iomaster) &
   call MPI_FILE_WRITE_AT(wfile, disp, calc_type%postscf_name, LEN(calc_type%postscf_name), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
  disp = disp + STORAGE_SIZE(calc_type%postscf_name)

  if(is_iomaster) &
   call MPI_FILE_WRITE_AT(wfile, disp, sf%nprodbasis_total, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
  disp = disp + STORAGE_SIZE(sf%nprodbasis_total)

  if(is_iomaster) &
   call MPI_FILE_WRITE_AT(wfile, disp, sf%npole_reso, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
  disp = disp + STORAGE_SIZE(sf%npole_reso)

  if(is_iomaster) &
   call MPI_FILE_WRITE_AT(wfile, disp, sf%pole, sf%npole_reso, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
  disp = disp + sf%npole_reso * STORAGE_SIZE(sf%pole(1))


  if( has_auxil_basis ) then
    !
    ! Write the residue in "the" universal ordering that does not depend on the
    ! data distribution
    allocate(buffer(sf%npole_reso))
    do ibf_auxil=1, sf%nprodbasis_total
      if( auxil%rank == iproc_ibf_auxil(ibf_auxil) ) then

        buffer(:) = sf%residue_left(ibf_auxil_l(ibf_auxil), :)
        call MPI_FILE_WRITE_AT(wfile, disp, buffer, sf%npole_reso, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)

      endif
      disp = disp + sf%npole_reso * STORAGE_SIZE(sf%residue_left(1, 1))
    enddo
    deallocate(buffer)
  else
    if(is_iomaster) then
      do iprodbasis=1, sf%nprodbasis_total
        call MPI_FILE_WRITE_AT(wfile, disp, sf%residue_left(iprodbasis, :), sf%npole_reso, &
                               MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
        disp = disp + sf%npole_reso * STORAGE_SIZE(sf%residue_left(1, 1))
      enddo
    endif
  endif

  call MPI_FILE_CLOSE(wfile, ierr)

#endif


end subroutine write_spectral_function


!=========================================================================
subroutine read_spectral_function(sf, reading_status)
  implicit none
  type(spectral_function), intent(inout) :: sf
  integer, intent(out)                   :: reading_status
  !=====
  integer            :: wfile
  character(len=100) :: postscf_name_read
  logical            :: file_exists
  integer            :: npole_read, nprodbasis_read
#if defined(HAVE_MPI)
  integer            :: ibf_auxil
  integer            :: iprodbasis
  integer            :: ierr
  integer(kind=MPI_OFFSET_KIND) :: disp
  real(dp), allocatable :: buffer(:)
#else
  integer :: ipole_read
#endif
  !=====

  write(stdout, '(/,a)') ' Try to read spectral function from file SCREENED_COULOMB'

  inquire(file='SCREENED_COULOMB', exist=file_exists)
  if( .NOT. file_exists ) then
    write(stdout, '(a,/)') ' File does not exist'
    reading_status=1
    return
  endif

#if !defined(HAVE_MPI)
  open(newunit=wfile, file='SCREENED_COULOMB', status='old', form='unformatted')

  read(wfile) postscf_name_read
  read(wfile) nprodbasis_read
  read(wfile) npole_read

  sf%npole_reso = npole_read
  call allocate_spectral_function(nprodbasis_read, sf)


  read(wfile) sf%pole(:)
  do ipole_read=1, npole_read
    read(wfile) sf%residue_left(:, ipole_read)
  enddo

  reading_status=0
  msg='reading spectral function from SCREENED_COULOMB obtained from '//TRIM(postscf_name_read)
  call issue_warning(msg)

  close(wfile)
#else

  write(stdout, *) 'Reading file SCREENED_COULOMB'
  call MPI_FILE_OPEN(MPI_COMM_WORLD, 'SCREENED_COULOMB', &
                    MPI_MODE_RDONLY, &
                    MPI_INFO_NULL, wfile, ierr)

  disp = 0
  call MPI_FILE_READ_AT(wfile, disp, postscf_name_read, LEN(postscf_name_read), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
  disp = disp + STORAGE_SIZE(postscf_name_read)

  call MPI_FILE_READ_AT(wfile, disp, nprodbasis_read, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
  disp = disp + STORAGE_SIZE(nprodbasis_read)

  sf%nprodbasis_total = nprodbasis_read

  call MPI_FILE_READ_AT(wfile, disp, npole_read, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
  disp = disp + STORAGE_SIZE(npole_read)

  sf%npole_reso = npole_read

  if( has_auxil_basis ) then
    call allocate_spectral_function(nauxil_local, sf)
  else
    call allocate_spectral_function(nprodbasis_read, sf)
  endif

  call MPI_FILE_READ_AT(wfile, disp, sf%pole, sf%npole_reso, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
  disp = disp + sf%npole_reso * STORAGE_SIZE(sf%pole(1))


  if( has_auxil_basis ) then
    !
    ! Read the residue from "the" universal ordering that does not depend on the
    ! data distribution
    allocate(buffer(sf%npole_reso))
    do ibf_auxil=1, nauxil_global
      if( auxil%rank == iproc_ibf_auxil(ibf_auxil) ) then
        call MPI_FILE_READ_AT(wfile, disp, buffer, sf%npole_reso, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
        sf%residue_left(ibf_auxil_l(ibf_auxil), :) = buffer(:)
      endif
      disp = disp + sf%npole_reso * STORAGE_SIZE(sf%residue_left(1, 1))
    enddo
    deallocate(buffer)
  else
    do iprodbasis=1, sf%nprodbasis
      call MPI_FILE_READ_AT(wfile, disp, sf%residue_left(iprodbasis, :), sf%npole_reso, &
                            MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
      disp = disp + sf%npole_reso * STORAGE_SIZE(sf%residue_left(1, 1))
    enddo
  endif



  call MPI_FILE_CLOSE(wfile, ierr)

  msg='reading spectral function from SCREENED_COULOMB obtained from '//TRIM(postscf_name_read)
  call issue_warning(msg)
  reading_status=0

#endif


end subroutine read_spectral_function


!=========================================================================
subroutine sf_evaluate_several_omegas(sf, omega_cmplx, chi)
  implicit none

  class(spectral_function), intent(in) :: sf
  complex(dp), intent(in) :: omega_cmplx(:)
  real(dp), intent(out) :: chi(:, :, :)
  !=====
  integer :: nomega, iomega, ipole
  integer :: jauxil, iauxil
  real(dp), allocatable :: tmp(:, :)
  !=====
  if( nauxil_global /= nauxil_local ) call die('sf_evaluate_several_omegas: not implemented with distributed auxiliary basis')
  if( .NOT. ALLOCATED(sf%residue_left) ) call die('sf_evaluate_severals_omegas: should have sf%residue_left available')

  nomega = SIZE(omega_cmplx)

  chi(:, :, :) = 0.0_dp

  !
  ! for real frequencies, use a naive implementation
  !
  if( ANY( ABS(omega_cmplx(:)%re) > 1.0e-6_dp ) ) then
    do iomega=1, nomega
      do ipole=1, sf%npole_reso
        do jauxil=1, sf%nprodbasis
          chi(:, jauxil, iomega) = chi(:, jauxil, iomega) &
                 + sf%residue_left(:, ipole) * sf%residue_left(jauxil, ipole) &
                       * REAL( 1.0_dp / ( omega_cmplx(iomega) - sf%pole(ipole) + ieta ) &
                              -1.0_dp / ( omega_cmplx(iomega) + sf%pole(ipole) - ieta ) )
        enddo
      enddo
    enddo

  else

    allocate(tmp, MOLD=sf%residue_left)

    do iomega=1, nomega
      tmp(:, :) = sf%residue_left(:, :)
      do ipole=1, sf%npole_reso
        tmp(:, ipole) = tmp(:, ipole) * SQRT( 2.0_dp * sf%pole(ipole) / ( ABS(omega_cmplx(iomega))**2 + sf%pole(ipole)**2 ) )
      enddo

      call DSYRK('L', 'N', nauxil_global, sf%npole_reso, -1.0d0, tmp, nauxil_global, 0.0d0, chi(:, :, iomega),nauxil_global)
      do jauxil=1, nauxil_global
        do iauxil=jauxil+1, nauxil_global
          chi(jauxil, iauxil, iomega) = chi(iauxil, jauxil, iomega)
        enddo
      enddo
    enddo

    deallocate(tmp)
  endif


end subroutine sf_evaluate_several_omegas


!=========================================================================
subroutine sf_evaluate_one_real_omega(sf, omega_real, chi)
  implicit none

  class(spectral_function), intent(in) :: sf
  real(dp), intent(in) :: omega_real
  complex(dp), intent(out) :: chi(:, :)
  !=====
  integer :: ipole
  integer :: jauxil
  !=====
  if( nauxil_global /= nauxil_local ) call die('sf_evaluate_one_omega: not implemented with distributed auxiliary basis')
  if( .NOT. ALLOCATED(sf%residue_left) ) call die('sf_evaluate_one_omega: should have sf%residue_left available')


  chi(:, :) = 0.0_dp

  !
  ! for real frequencies, use a naive implementation
  !
  do ipole=1, sf%npole_reso
    do jauxil=1, sf%nprodbasis
      chi(:, jauxil) = chi(:, jauxil) &
             + sf%residue_left(:, ipole) * sf%residue_left(jauxil, ipole) &
                   * ( 1.0_dp / ( omega_real - sf%pole(ipole) + ieta ) &
                      -1.0_dp / ( omega_real + sf%pole(ipole) - ieta ) )
    enddo
  enddo



end subroutine sf_evaluate_one_real_omega


!=========================================================================
subroutine sf_evaluate_one_omega(sf, omega_cmplx, chi)
  implicit none

  class(spectral_function), intent(in) :: sf
  complex(dp), intent(in) :: omega_cmplx
  real(dp), intent(out) :: chi(:, :)
  !=====
  integer :: ipole
  integer :: jauxil, iauxil
  real(dp), allocatable :: tmp(:, :)
  !=====
  if( nauxil_global /= nauxil_local ) call die('sf_evaluate_one_omega: not implemented with distributed auxiliary basis')
  if( .NOT. ALLOCATED(sf%residue_left) ) call die('sf_evaluate_one_omega: should have sf%residue_left available')


  chi(:, :) = 0.0_dp

  !
  ! for real frequencies, use a naive implementation
  !
  if( ABS(omega_cmplx%re) > 1.0e-6_dp  ) then
      do ipole=1, sf%npole_reso
        do jauxil=1, sf%nprodbasis
          chi(:, jauxil) = chi(:, jauxil) &
                 + sf%residue_left(:, ipole) * sf%residue_left(jauxil, ipole) &
                       * REAL( 1.0_dp / ( omega_cmplx - sf%pole(ipole) + ieta ) &
                              -1.0_dp / ( omega_cmplx + sf%pole(ipole) - ieta ) )
        enddo
      enddo

  else

    allocate(tmp, MOLD=sf%residue_left)

    tmp(:, :) = sf%residue_left(:, :)
    do ipole=1, sf%npole_reso
      tmp(:, ipole) = tmp(:, ipole) * SQRT( 2.0_dp * sf%pole(ipole) / ( ABS(omega_cmplx)**2 + sf%pole(ipole)**2 ) )
    enddo

    call DSYRK('L', 'N', nauxil_global, sf%npole_reso, -1.0d0, tmp, nauxil_global, 0.0d0, chi(:,:),nauxil_global)
    do jauxil=1, nauxil_global
      do iauxil=jauxil+1, nauxil_global
        chi(jauxil, iauxil) = chi(iauxil, jauxil)
      enddo
    enddo

    deallocate(tmp)
  endif


end subroutine sf_evaluate_one_omega


!=========================================================================
subroutine sf_vsqrt_chi_vsqrt_rpa(sf, occupation, energy, c_matrix, low_rank, verbose)
  implicit none

  class(spectral_function), intent(inout) :: sf
  real(dp), intent(in)                    :: occupation(:, :)
  real(dp), intent(in)                    :: energy(:, :)
  real(dp), intent(in)                    :: c_matrix(:, :, :)
  logical, optional, intent(in)            :: low_rank, verbose
  !=====
  real(dp), parameter   :: TOL_LOW_EIGVAL = 1.0e-3_dp
  logical              :: verbose_, low_rank_, eri_3center_mo_available
  integer              :: stdout_
  integer              :: nstate, non_negligible
  integer              :: iomega, ieig, jeig, jauxil
  integer              :: t_ia
  integer              :: istate, astate, iaspin
  integer              :: info
  real(dp)             :: docc, de, factor
  real(dp), allocatable :: eri3_t1(:, :), eri3_t2(:, :)
  real(dp), allocatable :: chi0(:, :)
  real(dp), allocatable :: chi0tmp(:, :)
  real(dp)             :: eigval(nauxil_global)
  !=====

  call start_clock(timing_rpa_dynamic)

  stdout_  = stdout
  verbose_ = .TRUE.
  if( PRESENT(verbose) ) then
    verbose_ = verbose
    if( .NOT. verbose ) then
      open(newunit=stdout_, file='/dev/null')
    endif
  endif
  low_rank_ = .FALSE.
  if( PRESENT(low_rank) ) then
    low_rank_ = low_rank
  endif

  nstate = SIZE(occupation, DIM=1)

  write(stdout_, '(/,1x,a,i5,a)') 'Calculation of RPA polarizability on complex grid of ', sf%nomega, ' frequencies'


  if( sf%nomega < 1 ) call die('sf_vsqrt_chi_vsqrt_rpa: sf contains no frequency (nomega==0)')
  if( .NOT. has_auxil_basis ) call die('sf_vsqrt_chi_vsqrt_rpa: only works with an auxiliary basis')
  if( nauxil_global /= nauxil_local ) call die('sf_vsqrt_chi_vsqrt_rpa: not implemented with distributed auxiliary basis')

  sf%nprodbasis = nauxil_local
  if( low_rank_ ) then
    allocate(sf%vchiv_sqrt(sf%nomega))
  endif

  ! Check if (I | p q) integrals are already available
  !   if not, then calculate them
  eri_3center_mo_available = ALLOCATED(eri_3center_eigen)
  if( .NOT. eri_3center_mo_available ) then
    call calculate_eri_3center_eigen(c_matrix, ncore_W+1, nhomo_W, nlumo_W, nvirtual_W-1, timing=timing_aomo_pola)
  else
    ! eri_3center_eigen is already available
    ! check if it has the correct dimensions
    if(    LBOUND(eri_3center_eigen, DIM=2) > ncore_W+1      &
      .OR. UBOUND(eri_3center_eigen, DIM=2) < nhomo_W        &
      .OR. LBOUND(eri_3center_eigen, DIM=3) > nlumo_W        &
      .OR. UBOUND(eri_3center_eigen, DIM=3) > nvirtual_W-1 ) then
      call die('sf_vsqrt_chi_vsqrt_rpa: eri_3center_eigen does not contain all the needed states')
    endif
  endif

  sf%mchi = nauxil_global
  sf%nchi = nauxil_global
  call DESCINIT(sf%desc_chi, nauxil_global, nauxil_global, block_row, block_col, first_row, first_col, &
                cntxt_sd, MAX(1, sf%mchi), info)
  call clean_allocate('Chi', sf%chi, sf%mchi, sf%nchi, sf%nomega, verbose=verbose_)
  write(stdout_, '(1x,a,i7,a,i7)') 'Matrix sizes   ', nauxil_global, ' x ', nauxil_global
  write(stdout_, '(1x,a,i7,a,i7)') 'Distributed in ', sf%mchi, ' x ', sf%nchi

  call clean_allocate('TMP 3-center MO integrals', eri3_t1, nauxil_local, sf%npole_reso, verbose=verbose_)
  call clean_allocate('TMP 3-center MO integrals', eri3_t2, nauxil_local, sf%npole_reso, verbose=verbose_)
  call clean_allocate('Chi0', chi0, sf%mchi, sf%nchi, verbose=verbose_)


  do iomega=1, sf%nomega

    write(stdout_, '(1x,a,i4,a,i4)') 'Loop on frequencies: ', iomega, ' / ', sf%nomega

    !
    ! First evaluate v^{1/2} \chi_0 v^{1/2}
    !
    ! Loop over resonant transitions
    do t_ia=1, sf%npole_reso
      istate = sf%transition_table(1, t_ia)
      astate = sf%transition_table(2, t_ia)
      iaspin = sf%transition_table(3, t_ia)

      docc = occupation(istate, iaspin) - occupation(astate, iaspin)
      de   = energy(astate, iaspin)     - energy(istate, iaspin)
      factor = REAL( 2.0_dp * docc * de / ( sf%omega(iomega)**2 - de**2 ) )

      eri3_t1(:, t_ia) = eri_3center_eigen(:, istate, astate, iaspin) * factor
      eri3_t2(:, t_ia) = eri_3center_eigen(:, istate, astate, iaspin)

    enddo

#if defined(HAVE_MKL)
    call DGEMMT('L', 'N', 'T', nauxil_global, sf%npole_reso, 1.0d0, eri3_t1, nauxil_global,eri3_t2,nauxil_global, &
                0.0d0, chi0, nauxil_global)
    call matrix_lower_to_full(chi0)
#else
    call DGEMM('N', 'T', nauxil_global, nauxil_global, sf%npole_reso, 1.0d0, eri3_t1, nauxil_global, eri3_t2,nauxil_global, &
               0.0d0, chi0, nauxil_global)
#endif




    allocate(chi0tmp, SOURCE=chi0)

    ! Second
    !
    !   Diagonalize chi0 (low_rank route)
    ! or 
    !   invert (1 - chi0)  (standard route)
    !
    if( low_rank_ ) then
      call diagonalize(postscf_diago_flavor, chi0tmp, eigval)
      !
      ! Transform vsqrt * chi0 * vsqrt in to vsqrt * chi * vsqrt
      eigval(:) = eigval(:) / (1.0_dp - eigval(:))
      non_negligible = COUNT( ABS(eigval(:)) > TOL_LOW_EIGVAL )
      write(stdout_, *) 'Number of non-negligible eigenvalues:', non_negligible, nauxil_global

      allocate(sf%vchiv_sqrt(iomega)%eigvec(nauxil_global, non_negligible))
      allocate(sf%vchiv_sqrt(iomega)%eigval(non_negligible))
      !
      ! Store the eigenelements of vsqrt_chi_vsqrt
      !
      ieig = 0
      do jeig=1, nauxil_global
        if( ABS(eigval(jeig)) > TOL_LOW_EIGVAL ) then
          ieig = ieig + 1
          sf%vchiv_sqrt(iomega)%eigvec(:, ieig) = chi0tmp(:, jeig)
          sf%vchiv_sqrt(iomega)%eigval(ieig)   = eigval(jeig)
        endif
      enddo

    else
      !
      ! Calculate and store the full matrix vsqrt_chi_vsqrt
      !
      chi0tmp(:, :) = -chi0tmp(:, :)
      do jauxil=1, nauxil_global
        chi0tmp(jauxil, jauxil) = 1.0_dp + chi0tmp(jauxil, jauxil)
      enddo
      call invert_symmetric(chi0tmp)
      !sf%chi(:,:,iomega) = MATMUL( chi0tmp, chi0 )
      call DSYMM('L', 'L', nauxil_global, nauxil_global, 1.0d0, chi0tmp, nauxil_global, chi0,nauxil_global, &
                 0.0d0, sf%chi(:, :, iomega), nauxil_global)

    endif

    deallocate(chi0tmp)


  enddo

  call clean_deallocate('TMP 3-center MO integrals', eri3_t1, verbose=verbose_)
  call clean_deallocate('TMP 3-center MO integrals', eri3_t2, verbose=verbose_)
  call clean_deallocate('Chi0', chi0, verbose=verbose_)

  if( .NOT. eri_3center_mo_available ) then
    call destroy_eri_3center_eigen()
  endif

  if( PRESENT(verbose) ) then
    if( .NOT. verbose ) then
      close(stdout_)
    endif
  endif

  call stop_clock(timing_rpa_dynamic)

end subroutine sf_vsqrt_chi_vsqrt_rpa


!=========================================================================
subroutine sf_interpolate_vsqrt_chi_vsqrt(sf, omega, vchiv_sqrt_omega)
  implicit none

  class(spectral_function), intent(in) :: sf
  real(dp), intent(in)                 :: omega
  type(chi_type), intent(out)          :: vchiv_sqrt_omega
  !=====
  integer  :: jomega
  !=====


  if( ANY( ABS(sf%omega(:)%im) > 1.0e-6_dp ) ) then
    write(stdout, *) sf%omega(:)
    call die('sf_interpolate_vsqrt_chi_vsqrt_rpa: for real frequencies only')
  endif

  if( omega > MAXVAL(sf%omega(:)%re) .OR. omega < MINVAL(sf%omega(:)%re) ) then
    write(stdout, *) omega
    write(stdout, *) sf%omega(1)%re
    write(stdout, *) sf%omega(SIZE(sf%omega))%re
    call die('sf_interpolate_vsqrt_chi_vsqrt_rpa: requested frequency out of range')
  endif

  ! no interpolation, but just take the closest omega
  jomega = MINLOC( ABS(sf%omega(:)%re-omega) , DIM=1 )

  allocate(vchiv_sqrt_omega%eigvec, SOURCE=sf%vchiv_sqrt(jomega)%eigvec)
  allocate(vchiv_sqrt_omega%eigval, SOURCE=sf%vchiv_sqrt(jomega)%eigval)

end subroutine sf_interpolate_vsqrt_chi_vsqrt


!=========================================================================
subroutine ct_destroy(chi)
  implicit none
  !=====
  class(chi_type), intent(inout) :: chi
  !=====
  if( ALLOCATED(chi%eigvec) ) deallocate(chi%eigvec)
  if( ALLOCATED(chi%eigval) ) deallocate(chi%eigval)

end subroutine ct_destroy


end module m_spectral_function
!=========================================================================
