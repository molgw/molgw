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
 use,intrinsic :: iso_c_binding, only: C_INT,C_DOUBLE,C_BOOL

 integer,parameter  :: dp = KIND(0.d0)
 integer,parameter  :: sp = KIND(0.0)

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

 !
 ! Restart file types
 integer,parameter ::           NO_RESTART = 0
 integer,parameter ::        SMALL_RESTART = 1
 integer,parameter ::          BIG_RESTART = 2
 integer,parameter ::        BASIS_RESTART = 3
 integer,parameter :: EMPTY_STATES_RESTART = 4

#ifdef NEED_NORM2
 interface norm2
   module procedure norm2_ra1
   module procedure norm2_ra2
   module procedure norm2_ra2_dim
   module procedure norm2_ra3
   module procedure norm2_ra4
 end interface
#endif

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


#ifdef NEED_NORM2
!=========================================================================
function norm2_ra1(array) RESULT(norm2)
 implicit none
 real(dp),intent(in) :: array(:)
 real(dp)            :: norm2
!=====
!=====

 norm2 = SQRT( SUM( array**2 ) )

end function norm2_ra1


!=========================================================================
function norm2_ra2(array) RESULT(norm2)
 implicit none
 real(dp),intent(in)          :: array(:,:)
 real(dp)                     :: norm2
!=====
!=====

 norm2 = SQRT( SUM( array**2 ) )

end function norm2_ra2


!=========================================================================
function norm2_ra2_dim(array,dim) RESULT(norm2)
 implicit none
 real(dp),intent(in) :: array(:,:)
 integer,intent(in)  :: dim
 real(dp)            :: norm2(1:SIZE(array,DIM=2))
!=====
!=====

 if( dim /= 1 ) stop 'NORM2 not correctly implemented. Use a newer compiler that implemented NORM2'

 norm2(:) = SQRT( SUM( array(:,:)**2 , DIM=1) )

end function norm2_ra2_dim


!=========================================================================
function norm2_ra3(array) RESULT(norm2)
 implicit none
 real(dp),intent(in) :: array(:,:,:)
 real(dp)            :: norm2
!=====
!=====

 norm2 = SQRT( SUM( array**2 ) )

end function norm2_ra3


!=========================================================================
function norm2_ra4(array) RESULT(norm2)
 implicit none
 real(dp),intent(in) :: array(:,:,:,:)
 real(dp)            :: norm2
!=====
!=====

 norm2 = SQRT( SUM( array**2 ) )

end function norm2_ra4
#endif


!=========================================================================
end module m_definitions
!=========================================================================
