!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the most basic definitions. It should be "used" in all routines.
!
!=========================================================================
#include "molgw.h"
module m_definitions
  use,intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT,ERROR_UNIT
  use,intrinsic :: ISO_C_BINDING, only: C_INT,C_LONG,C_DOUBLE,C_BOOL,C_PTR,C_CHAR,C_NULL_PTR,C_F_POINTER,C_NULL_CHAR
#if defined(_OPENMP)
  use OMP_LIB, only: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,OMP_GET_MAX_THREADS
#endif

#if defined(DEBUG)
  logical,parameter :: debug=.TRUE.
#else
  logical,parameter :: debug=.FALSE.
#endif

  integer,parameter  :: dp = KIND(0.0d0)
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
  real(dp), parameter :: zero   = 0.0_dp
  real(dp), parameter :: half   = 0.5_dp
  real(dp), parameter :: one    = 1.0_dp
  real(dp), parameter :: two    = 2.0_dp
  real(dp), parameter :: four   = 4.0_dp
  real(dp), parameter :: eight  = 8.0_dp
  real(dp), parameter :: ten    = 1.0d1
  real(dp), parameter :: twelve = 1.2d1
  complex(dp),parameter  :: im    = (0.0_dp,1.0_dp)
  complex(dp),parameter  :: COMPLEX_ONE  = (1.0_dp,0.0_dp)
  complex(dp),parameter  :: COMPLEX_ZERO = (0.0_dp,0.0_dp)

  !
  ! Thresholds
  real(dp),parameter :: completely_empty = 1.0e-5_dp
  real(dp), parameter :: tol1 = 1.0e-1_dp
  real(dp), parameter :: tol3 = 1.0e-3_dp
  real(dp), parameter :: tol5 = 1.0e-5_dp
  real(dp), parameter :: tol6 = 1.0e-6_dp
  real(dp), parameter :: tol8 = 1.0e-8_dp
  real(dp), parameter :: tol9 = 1.0e-9_dp
  real(dp), parameter :: tol16 = 1.0e-16_dp
  real(dp), parameter :: tol20 = 1.0e-20_dp
  real(dp), parameter :: thousand = 1.0e3_dp


  !
  ! Quality levels used of grids and Coulomb integrals
  integer,parameter :: low       = 10
  integer,parameter :: medium    = 20
  integer,parameter :: high      = 30
  integer,parameter :: very_high = 40
  integer,parameter :: insane    = 50

  !
  ! Name of the output file used by NOFT calcs.
  character(len=140)::output_name

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
