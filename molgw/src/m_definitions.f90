!=========================================================================
module m_definitions
 use,intrinsic :: iso_fortran_env, only: OUTPUT_UNIT,ERROR_UNIT

 integer,parameter  :: dp=kind(0.d0)
 integer,parameter  :: sp=kind(0.0)

 integer,parameter  :: prec_eri=dp
#ifndef TD_SP
 integer,parameter  :: prec_td=dp
#else
 integer,parameter  :: prec_td=sp
#endif

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
 real(dp),parameter    :: pi    =3.14159265358979323_dp
 real(dp),parameter    :: pi2   =pi**2
 complex(dp),parameter :: im    =(0.0_dp,1.0_dp)


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

end module m_definitions
!=========================================================================
