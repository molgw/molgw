!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods and data concerning the periodic table elements
!
!=========================================================================
module m_elements
 use m_definitions
 use m_warning,only: die


 integer,parameter          :: nelement_max = 103
 character(len=2),parameter :: element_list(nelement_max) =                            &
  (/' H',                                                                                'He', &  !   2
    'Li','Be',                                                  ' B',' C',' N',' O',' F','Ne', &  !  10
    'Na','Mg',                                                  'Al','Si',' P',' S','Cl','Ar', &  !  18
    ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &  !  36
    'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe', &  !  54
    'Cs','Ba',                                                                                 &  !  56
    'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',                     &  !  70
              'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &  !  86
    'Fr','Ra',                                                                                 &  !  88
              'Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',           &  ! 102
    'Lr'                                                                                       &  ! 103
  /)





contains


!=========================================================================
function element_core(zval,zatom)
 implicit none
 real(dp),intent(in) :: zval
 real(dp),intent(in) :: zatom
 integer             :: element_core
!=====

 !
 ! If zval /= zatom, this is certainly an effective core potential
 ! and no core states should be frozen.
 if( ABS(zval - zatom) > 1.0e-3_dp ) then
   element_core = 0
 else

   if( zval <= 4.00001 ) then  ! up to Be
     element_core = 0
   else if( zval <= 12.00001 ) then  ! up to Mg
     element_core = 1
   else if( zval <= 30.00001 ) then  ! up to Ca
     element_core = 5
   else if( zval <= 48.00001 ) then  ! up to Sr
     element_core = 9
   else
     call die('not imlemented in element_core')
   endif

 endif


end function element_core


!=========================================================================
function element_covalent_radius(zatom)
 implicit none
 integer,intent(in) :: zatom
 real(dp)           :: element_covalent_radius
!=====

 !
 ! Data from Cambridge Structural Database
 ! http://en.wikipedia.org/wiki/Covalent_radius
 !
 ! Values are first given in picometer
 ! They will be converted in bohr just after
 select case(zatom)
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
   ! If not available, then invent one!
   element_covalent_radius = 150.
 end select

 !
 ! pm to bohr conversion
 element_covalent_radius = element_covalent_radius / bohr_A * 0.01_dp

end function element_covalent_radius


!=========================================================================
function element_number(element_name)
 implicit none
 character(len=*),intent(in) :: element_name
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

 element_name = element_list(NINT(ABS(zatom)))


end function element_name


!=========================================================================
function element_name_long(zatom)
 implicit none
 real(dp),intent(in) :: zatom
 character(len=8)  :: element_name_long
!=====

 if( NINT(zatom) > nelement_max ) then
   write(stdout,'(a,i3,a)') 'Element symbol is not one of first ',nelement_max,' elements'
   call die('element symbol not understood')
 endif

 if( zatom > 0.0 .AND. ABS(NINT(zatom)-zatom) < 1.0e-9_dp ) then
   element_name_long = element_list(NINT(ABS(zatom)))
 else
   write(element_name_long,'(f8.4)') zatom
 endif


end function element_name_long


!=========================================================================
subroutine element_atomicdensity(zval,zatom,coeff,alpha)
 implicit none
 real(dp),intent(in)  :: zval
 real(dp),intent(in)   :: zatom
 real(dp),intent(out) :: coeff(4)
 real(dp),intent(out) :: alpha(4)
!=====

 select case(NINT(zatom))
 case( 1)
  alpha(1) = 0.191852
  alpha(2) = 5.6143
  alpha(3) = 0.510303
  alpha(4) = 1.47022
  coeff(1) = 0.208981
  coeff(2) = 0.0287008
  coeff(3) = 0.358982
  coeff(4) = 0.168689
 case( 6)
   alpha(1) = 0.294412
   alpha(2) = 68.2412
   alpha(3) = 16.7565
   alpha(4) = 0.852332
   coeff(1) = 1.90506
   coeff(2) = 0.485794
   coeff(3) = 1.2976
   coeff(4) = 2.24937
 case( 7)
   alpha(1) = 0.404121
   alpha(2) = 23.2108
   alpha(3) = 1.24119
   alpha(4) = 93.8394
   coeff(1) = 2.37996
   coeff(2) = 1.29369
   coeff(3) = 2.7861
   coeff(4) = 0.49092
 case( 8)
   alpha(1) = 122.833
   alpha(2) = 30.5529
   alpha(3) = 0.536011
   alpha(4) = 1.71315
   coeff(1) = 0.495436
   coeff(2) = 1.27809
   coeff(3) = 2.85614
   coeff(4) = 3.27719
 case(79) ! Au
   if( NINT(zval) == 19 ) then   ! Au ECP
     coeff(1) = zval -2.0_dp + 4.1_dp
     alpha(1) = 1.0_dp
     coeff(2) = 0.2_dp
     alpha(2) = 2.0_dp
     coeff(3) =-4.1_dp
     alpha(3) =  3.0_dp
     coeff(4) = 0.0_dp
     alpha(4) = 2.0_dp
   else
     alpha(1) = 0.3_dp
     coeff(1) = 2.0_dp
     alpha(2) = 2.6_dp
     coeff(2) = zval - coeff(1)
     coeff(3) = 0.0_dp
     coeff(4) = 0.0_dp
     alpha(3) = 1.0_dp
     alpha(4) = 1.0_dp
   endif
 case default
   alpha(1) = 0.3_dp
   coeff(1) = 2.0_dp
   alpha(2) = 2.6_dp
   coeff(2) = zval - coeff(1)
   coeff(3) = 0.0_dp
   coeff(4) = 0.0_dp
   alpha(3) = 1.0_dp
   alpha(4) = 1.0_dp
 end select

 ! Ensure the proper normalization
 coeff(:) = zval / SUM(coeff(:)) * coeff(:)

end subroutine element_atomicdensity


!=========================================================================
end module m_elements
