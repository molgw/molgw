!=========================================================================
#include "macros.h"
!=========================================================================
module m_definitions

 integer,parameter  :: dp=kind(0.d0)
 integer,parameter  :: sp=kind(0.)
 integer,parameter  :: prec_eri=dp

 !
 ! Physical constants
 real(dp),parameter :: Ha_eV =27.2113961_dp
 real(dp),parameter :: bohr_A=0.529177249_dp
 real(dp),parameter :: pi    =3.14159265358979323_dp
 real(dp),parameter :: pi2   =pi**2
 complex(dp)        :: im=(0.0_dp,1.0_dp)

 real(dp),parameter :: completely_empty=1.0e-5_dp

end module m_definitions

