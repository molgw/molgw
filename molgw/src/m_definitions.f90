!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the most basic definitions. It should be "used" in all routines.
!
!=========================================================================
module m_definitions
 use,intrinsic :: iso_fortran_env, only: OUTPUT_UNIT,ERROR_UNIT

 integer,parameter  :: dp=KIND(0.d0)
 integer,parameter  :: sp=KIND(0.0)
 integer,parameter  :: dpc=KIND((0.d0,0.d0))

 integer,parameter  :: prec_eri=dp

 integer,parameter  :: stdout = OUTPUT_UNIT
 integer,parameter  :: stderr = ERROR_UNIT


 !
 ! Physical constants
 ! Values from NIST CODATA 2010
 real(dp),parameter    :: Ha_eV        = 27.21138505_dp
 real(dp),parameter    :: bohr_A       = 0.52917721092_dp
 real(dp),parameter    :: c_speedlight = 137.035999074_dp


 !
 ! Mathematical constants
 real(dp),parameter     :: pi    =3.14159265358979323_dp
 real(dp),parameter     :: pi2   =pi**2
 complex(dpc),parameter :: im    =(0.0_dp,1.0_dp)


 ! 
 ! Thresholds
 real(dp),parameter :: completely_empty=1.0e-5_dp


 !
 ! Quality levels used of grids and Coulomb integrals
 integer,parameter :: low       = 10
 integer,parameter :: medium    = 20
 integer,parameter :: high      = 30
 integer,parameter :: very_high = 40
 integer,parameter :: insane    = 50

 !
 ! Restart file types
 integer,parameter ::           NO_RESTART = 0
 integer,parameter ::        SMALL_RESTART = 1
 integer,parameter ::          BIG_RESTART = 2
 integer,parameter ::        BASIS_RESTART = 3
 integer,parameter :: EMPTY_STATES_RESTART = 4


end module m_definitions


!=========================================================================
