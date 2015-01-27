!=========================================================================
#include "macros.h"
!=========================================================================
module m_elements
 use m_definitions
 use m_mpi








contains

!=========================================================================
function element_covalent_radius(zatom)
 implicit none
 real(dp),intent(in) :: zatom
 real(dp)            :: element_covalent_radius
!=====

 !
 ! Data from Cambridge Structural Database
 ! http://en.wikipedia.org/wiki/Covalent_radius
 !
 ! Values are first given in picometer
 ! They will be converted in bohr just after
 select case(NINT(zatom))
 case( 1)
   element_covalent_radius =  31
 case( 2)
   element_covalent_radius =  28
 case( 3)
   element_covalent_radius = 128
 case( 4)
   element_covalent_radius =  96
 case( 5)
   element_covalent_radius =  84
 case( 6)
   element_covalent_radius =  73
 case( 7)
   element_covalent_radius =  71
 case( 8)
   element_covalent_radius =  66
 case( 9)
   element_covalent_radius =  57
 case(10) ! Ne
   element_covalent_radius =  58
 case(11)
   element_covalent_radius = 166
 case(12)
   element_covalent_radius = 141
 case(13)
   element_covalent_radius = 121
 case(14)
   element_covalent_radius = 111
 case(15)
   element_covalent_radius = 107
 case(16)
   element_covalent_radius = 105
 case(17)
   element_covalent_radius = 102
 case(18) ! Ar
   element_covalent_radius = 106
 case(19)
   element_covalent_radius = 203
 case(20)
   element_covalent_radius = 176
 case(21)
   element_covalent_radius = 170
 case(22)
   element_covalent_radius = 160
 case(23)
   element_covalent_radius = 153
 case(24)
   element_covalent_radius = 139
 case(25)
   element_covalent_radius = 145
 case(26)
   element_covalent_radius = 145
 case(27)
   element_covalent_radius = 140
 case(28)
   element_covalent_radius = 124
 case(29)
   element_covalent_radius = 132
 case(30)
   element_covalent_radius = 122
 case(31)
   element_covalent_radius = 120
 case(32)
   element_covalent_radius = 119
 case(34)
   element_covalent_radius = 120
 case(35)
   element_covalent_radius = 120
 case(36)
   element_covalent_radius = 116
 case default
   stop'radius not available'
 end select
  
 !
 ! pm to bohr conversion
 element_covalent_radius = element_covalent_radius / bohr_A * 0.01_dp

end function element_covalent_radius


!=========================================================================
function element_name(zatom)
 implicit none
 real(dp),intent(in) :: zatom
 character(len=2)  :: element_name
!=====

 select case(NINT(zatom))
 case( 0)
   element_name='X'
 case( 1)
   element_name='H'
 case( 2)
   element_name='He'
 case( 3)
   element_name='Li'
 case( 4)
   element_name='Be'
 case( 5)
   element_name='B'
 case( 6)
   element_name='C'
 case( 7)
   element_name='N'
 case( 8)
   element_name='O'
 case( 9)
   element_name='F'
 case(10)
   element_name='Ne'
 case(11)
   element_name='Na'
 case(12)
   element_name='Mg'
 case(13)
   element_name='Al'
 case(14)
   element_name='Si'
 case(15)
   element_name='P'
 case(16)
   element_name='S'
 case(17)
   element_name='Cl'
 case(18)
   element_name='Ar'
 case(19)
   element_name='K'
 case(20)
   element_name='Ca'
 case(21)
   element_name='Sc'
 case(22)
   element_name='Ti'
 case(23)
   element_name='V'
 case(24)
   element_name='Cr'
 case(25)
   element_name='Mn'
 case(26)
   element_name='Fe'
 case(27)
   element_name='Co'
 case(28)
   element_name='Ni'
 case(29)
   element_name='Cu'
 case(30)
   element_name='Zn'
 case(31)
   element_name='Ga'
 case(32)
   element_name='Ge'
 case(33)
   element_name='As'
 case(34)
   element_name='Se'
 case(35)
   element_name='Br'
 case(36)
   element_name='Kr'
 case(47)
   element_name='Ag'
 case default
   WRITE_MASTER(*,*) 'Z too large: element not yet coded'
   WRITE_MASTER(*,*) 'Give a generic element name: X'
   element_name='X'
 end select

end function element_name










end module m_elements


!=========================================================================
