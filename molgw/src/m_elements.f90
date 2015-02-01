!=========================================================================
#include "macros.h"
!=========================================================================
module m_elements
 use m_definitions
 use m_mpi


 integer,parameter,private :: nelement_max = 54
 character(len=2),parameter,private :: element_list(nelement_max) =                            &
  (/' H',                                                                                'He', &  !  2
    'Li','Be',                                                  ' B',' C',' N',' O',' F','Ne', &  ! 10
    'Na','Mg',                                                  'Al','Si',' P',' S','Cl','Ar', &  ! 18
    ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &  ! 36
    'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe'  /) ! 54




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
   element_covalent_radius =  31.
 case( 2)
   element_covalent_radius =  28.
 case( 3)
   element_covalent_radius = 128.
 case( 4)
   element_covalent_radius =  96.
 case( 5)
   element_covalent_radius =  84.
 case( 6)
   element_covalent_radius =  73.
 case( 7)
   element_covalent_radius =  71.
 case( 8)
   element_covalent_radius =  66.
 case( 9)
   element_covalent_radius =  57.
 case(10) ! Ne.
   element_covalent_radius =  58.
 case(11)
   element_covalent_radius = 166.
 case(12)
   element_covalent_radius = 141.
 case(13)
   element_covalent_radius = 121.
 case(14)
   element_covalent_radius = 111.
 case(15)
   element_covalent_radius = 107.
 case(16)
   element_covalent_radius = 105.
 case(17)
   element_covalent_radius = 102.
 case(18) ! Ar.
   element_covalent_radius = 106.
 case(19)
   element_covalent_radius = 203.
 case(20)
   element_covalent_radius = 176.
 case(21)
   element_covalent_radius = 170.
 case(22)
   element_covalent_radius = 160.
 case(23)
   element_covalent_radius = 153.
 case(24)
   element_covalent_radius = 139.
 case(25)
   element_covalent_radius = 145.
 case(26)
   element_covalent_radius = 145.
 case(27)
   element_covalent_radius = 140.
 case(28)
   element_covalent_radius = 124.
 case(29)
   element_covalent_radius = 132.
 case(30)
   element_covalent_radius = 122.
 case(31)
   element_covalent_radius = 120.
 case(32)
   element_covalent_radius = 119.
 case(34)
   element_covalent_radius = 120.
 case(35)
   element_covalent_radius = 120.
 case(36) ! Kr.
   element_covalent_radius = 116.
 case default
   stop'radius not available'
 end select
  
 !
 ! pm to bohr conversion
 element_covalent_radius = element_covalent_radius / bohr_A * 0.01_dp

end function element_covalent_radius


!=========================================================================
function element_number(element_name)
 implicit none
 character(len=2),intent(in) :: element_name
 integer                     :: element_number
!=====
 integer :: ielement
!=====
 ielement=1
 do while( ADJUSTL(element_name) /= ADJUSTL(element_list(ielement)) )
   if( ielement == nelement_max ) then
     WRITE_MASTER(*,'(a,a)')    ' Input symbol ',element_name
     WRITE_MASTER(*,'(a,i3,a)') ' Element symbol is not one of first ',nelement_max,' elements'
     stop'element symbol not understood'
   endif
   ielement = ielement + 1
 enddo

 element_number = ielement
 

end function element_number


!=========================================================================
function element_name(zatom)
 implicit none
 real(dp),intent(in) :: zatom
 character(len=2)  :: element_name
!=====
 if( NINT(zatom) == 0 ) then
   element_name='X'
   return
 endif
 if( NINT(zatom) > nelement_max ) then
   WRITE_MASTER(*,'(a,i3,a)') 'Element symbol is not one of first ',nelement_max,' elements'
   stop'element symbol not understood'
 endif

 element_name = element_list(NINT(zatom))


end function element_name



end module m_elements


!=========================================================================
