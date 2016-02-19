!=========================================================================
! This file is part of MOLGW.
!=========================================================================
module m_elements
 use m_definitions
 use m_warning,only: die


 integer,parameter,private :: nelement_max = 54
 character(len=2),parameter,private :: element_list(nelement_max) =                            &
  (/' H',                                                                                'He', &  !  2
    'Li','Be',                                                  ' B',' C',' N',' O',' F','Ne', &  ! 10
    'Na','Mg',                                                  'Al','Si',' P',' S','Cl','Ar', &  ! 18
    ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &  ! 36
    'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe'  /) ! 54




contains


!=========================================================================
function element_core(zatom)
 implicit none
 real(dp),intent(in) :: zatom
 real(dp)            :: element_core
!=====

 if( zatom <= 4.00001 ) then  ! up to Be
  element_core = 0
 else if( zatom <= 12.00001 ) then  ! up to Mg
  element_core = 1
 else if( zatom <= 30.00001 ) then  ! up to Ca
  element_core = 5
 else if( zatom <= 48.00001 ) then  ! up to Sr
  element_core = 9
 else
   call die('not imlemented in element_core')
 endif


end function element_core


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
   call die('radius not available')
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
     write(stdout,'(a,a)')    ' Input symbol ',element_name
     write(stdout,'(a,i3,a)') ' Element symbol is not one of first ',nelement_max,' elements'
     call die('element symbol not understood')
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
   write(stdout,'(a,i3,a)') 'Element symbol is not one of first ',nelement_max,' elements'
   call die('element symbol not understood')
 endif

 element_name = element_list(NINT(zatom))


end function element_name


!=========================================================================
subroutine element_atomicdensity(zatom,coeff,alpha)
 implicit none
 real(dp),intent(in)  :: zatom
 real(dp),intent(out) :: coeff(4)
 real(dp),intent(out) :: alpha(4)
!=====

 select case(NINT(zatom))
 case( 1)
   alpha(1) = 35.7662 
   alpha(2) = 5.66045 
   alpha(3) = 1.34388 
   alpha(4) = 0.340268
   coeff(2) = 0.038544
   coeff(3) = 0.317984
   coeff(4) = 0.641901 
   coeff(1) = 1.0_dp - SUM(coeff(2:4))
!  alpha(1) = 0.191852  
!  alpha(2) = 5.6143    
!  alpha(3) = 0.510303  
!  alpha(4) = 1.47022   
!  coeff(1) = 0.208981  
!  coeff(2) = 0.0287008 
!  coeff(3) = 0.358982  
!  coeff(4) = 0.168689  
 case( 6)
   alpha(1) = 167.342
   alpha(2) =  28.2192
   alpha(3) =   0.96127
   alpha(4) = 1832.38 
   coeff(2) =  1.3663 
   coeff(3) =  4.49642
   coeff(4) =  0.00210702
   coeff(1) = 6.0_dp - SUM(coeff(2:4))
!   alpha(1) = 0.294412 
!   alpha(2) = 68.2412  
!   alpha(3) = 16.7565  
!   alpha(4) = 0.852332 
!   coeff(1) = 1.90506  
!   coeff(2) = 0.485794 
!   coeff(3) = 1.2976   
!   coeff(4) = 2.24937  
! case( 7)
!   alpha(1) = 0.404096 
!   alpha(2) = 93.8331  
!   alpha(3) = 1.24113  
!   alpha(4) = 23.2099  
!   coeff(1) = 5.34772  
!   coeff(2) = 1.10327  
!   coeff(3) = 6.26128  
!   coeff(4) = 2.90708  
! case( 8)
!   alpha(1) = 122.833 
!   alpha(2) = 30.5529 
!   alpha(3) = 0.536011
!   alpha(4) = 1.71315 
!   coeff(1) = 0.495436
!   coeff(2) = 1.27809 
!   coeff(3) = 2.85614 
!   coeff(4) = 3.27719 
 case default
   alpha(1) = 0.3_dp
   coeff(1) = 2.0_dp
   alpha(2) = 2.6_dp
   coeff(2) = zatom - coeff(1) 
   coeff(3) = 0.0_dp
   coeff(4) = 0.0_dp
   alpha(3) = 1.0_dp
   alpha(4) = 1.0_dp
 end select

 ! Ensure the proper normalization
 coeff(:) = zatom / SUM(coeff(:))

end subroutine element_atomicdensity


!=========================================================================
end module m_elements
