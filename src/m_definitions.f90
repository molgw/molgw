!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the most basic definitions. It should be "used" in all routines.
!
!=========================================================================
module m_definitions
  use,intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT,ERROR_UNIT
  use,intrinsic :: ISO_C_BINDING, only: C_INT,C_DOUBLE,C_BOOL,C_PTR,C_CHAR,C_NULL_PTR,C_F_POINTER
#if defined(_OPENMP)
  use,intrinsic :: OMP_LIB, only: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,OMP_GET_MAX_THREADS
#endif


  integer,parameter  :: dp = KIND(0.d0)
  integer,parameter  :: sp = KIND(0.0)

  integer,parameter  :: int8 = 8

  integer,protected  :: stdout = OUTPUT_UNIT
  integer,parameter  :: stderr = ERROR_UNIT

  integer,protected  :: MOLGW_LMAX

  !
  ! Physical constants
  ! Values from NIST CODATA 2010
  real(dp),parameter    :: au_debye     = 2.54174623105_dp
  real(dp),parameter    :: Ha_eV        = 27.21138505_dp
  real(dp),parameter    :: bohr_A       = 0.52917721092_dp
  real(dp),parameter    :: c_speedlight = 137.035999074_dp
  real(dp),parameter    :: Ha_K         = 315775.0431159572_dp


  !
  ! Mathematical constants
  real(dp),parameter     :: pi    = 3.14159265358979323_dp
  real(dp),parameter     :: pi2   = pi**2
  complex(dp),parameter  :: im    = (0.0_dp,1.0_dp)
  complex(dp),parameter  :: COMPLEX_ONE  = (1.0_dp,0.0_dp)
  complex(dp),parameter  :: COMPLEX_ZERO = (0.0_dp,0.0_dp)


  !
  ! Thresholds
  real(dp),parameter :: completely_empty = 1.0e-5_dp


  !
  ! Quality levels used of grids and Coulomb integrals
  integer,parameter :: low       = 10
  integer,parameter :: medium    = 20
  integer,parameter :: high      = 30
  integer,parameter :: very_high = 40
  integer,parameter :: insane    = 50


contains


!=========================================================================
subroutine set_molgw_lmax(lmax)
  implicit none
  integer,intent(in) :: lmax
  !=====
  !=====

  MOLGW_LMAX = lmax

end subroutine set_molgw_lmax


!=========================================================================
subroutine set_standard_output(unit_stdout)
  implicit none
  integer,intent(in) :: unit_stdout
  !=====
  !=====

  if( unit_stdout /= OUTPUT_UNIT ) then
    close(OUTPUT_UNIT)
    stdout = unit_stdout
    open(unit=stdout)
  endif

end subroutine set_standard_output


!=========================================================================
end module m_definitions
!=========================================================================
