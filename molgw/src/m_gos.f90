!=========================================================================
! This module contains
! the procedures to calculate the Generalized Oscillator Strengths (GOS)
!=========================================================================
module m_gos
 use m_definitions

 ! Maximum angular momentum hard-coded as of today
 integer,parameter :: gos_lmax=2

 ! Table for the fast calculation of the GOS(l,l')
 type gos_llp
   integer             :: np
   integer,allocatable :: mu(:)
   integer,allocatable :: alpha(:)
   integer,allocatable :: beta(:)
   integer,allocatable :: gamma(:)
   integer,allocatable :: delta(:)
   integer,allocatable :: epsilon(:)
 end type gos_llp

 ! Table containing all the GOS for all l, all l'
 type(gos_llp),allocatable :: gos(:,:)


contains


!=========================================================================
subroutine setup_gos_llp()
 implicit none
!=====
 integer :: l1,l2
!=====

 allocate(gos(0:gos_lmax,0:gos_lmax))

 ! new
 l1 = 0 ; l2 = 0
 call gos_allocate(l1,l2,1)
 gos(l1,l2)%mu     (1) = 1
 gos(l1,l2)%alpha  (1) = 0
 gos(l1,l2)%beta   (1) = 0
 gos(l1,l2)%gamma  (1) = 0
 gos(l1,l2)%delta  (1) = 0
 gos(l1,l2)%epsilon(1) = 0

 ! new
 l1 = 0 ; l2 = 1
 call gos_allocate(l1,l2,1)
 gos(l1,l2)%mu     (1) = 1
 gos(l1,l2)%alpha  (1) = 0
 gos(l1,l2)%beta   (1) = 1
 gos(l1,l2)%gamma  (1) = 0
 gos(l1,l2)%delta  (1) = 1
 gos(l1,l2)%epsilon(1) = 0

 ! new
 l1 = 0 ; l2 = 2
 call gos_allocate(l1,l2,3)
 gos(l1,l2)%mu     (1) = 1
 gos(l1,l2)%alpha  (1) = 0
 gos(l1,l2)%beta   (1) = 2
 gos(l1,l2)%gamma  (1) = 0
 gos(l1,l2)%delta  (1) = 2
 gos(l1,l2)%epsilon(1) = 0
 gos(l1,l2)%mu     (2) =-1
 gos(l1,l2)%alpha  (2) = 0
 gos(l1,l2)%beta   (2) = 0
 gos(l1,l2)%gamma  (2) = 0
 gos(l1,l2)%delta  (2) = 2
 gos(l1,l2)%epsilon(2) = 1
 gos(l1,l2)%mu     (3) = 1
 gos(l1,l2)%alpha  (3) = 0
 gos(l1,l2)%beta   (3) = 0
 gos(l1,l2)%gamma  (3) = 0
 gos(l1,l2)%delta  (3) = 1
 gos(l1,l2)%epsilon(3) = 0

 ! new
 l1 = 1 ; l2 = 0
 call gos_allocate(l1,l2,1)
 gos(l1,l2)%mu     (1) = 1
 gos(l1,l2)%alpha  (1) = 1
 gos(l1,l2)%beta   (1) = 0
 gos(l1,l2)%gamma  (1) = 1
 gos(l1,l2)%delta  (1) = 0
 gos(l1,l2)%epsilon(1) = 0
 
 ! new
 l1 = 1 ; l2 = 1
 call gos_allocate(l1,l2,2)
 gos(l1,l2)%mu     (1) = 1
 gos(l1,l2)%alpha  (1) = 1
 gos(l1,l2)%beta   (1) = 1
 gos(l1,l2)%gamma  (1) = 1
 gos(l1,l2)%delta  (1) = 1
 gos(l1,l2)%epsilon(1) = 0
 gos(l1,l2)%mu     (2) = 1
 gos(l1,l2)%alpha  (2) = 0
 gos(l1,l2)%beta   (2) = 0
 gos(l1,l2)%gamma  (2) = 1
 gos(l1,l2)%delta  (2) = 1
 gos(l1,l2)%epsilon(2) = 1

 ! new
 l1 = 1 ; l2 = 2
 call gos_allocate(l1,l2,4)
 gos(l1,l2)%mu     (1) = 1
 gos(l1,l2)%alpha  (1) = 1
 gos(l1,l2)%beta   (1) = 2
 gos(l1,l2)%gamma  (1) = 1
 gos(l1,l2)%delta  (1) = 2
 gos(l1,l2)%epsilon(1) = 0
 gos(l1,l2)%mu     (2) =-1
 gos(l1,l2)%alpha  (2) = 1
 gos(l1,l2)%beta   (2) = 0
 gos(l1,l2)%gamma  (2) = 1
 gos(l1,l2)%delta  (2) = 2
 gos(l1,l2)%epsilon(2) = 1
 gos(l1,l2)%mu     (3) = 1
 gos(l1,l2)%alpha  (3) = 1
 gos(l1,l2)%beta   (3) = 0
 gos(l1,l2)%gamma  (3) = 1
 gos(l1,l2)%delta  (3) = 1
 gos(l1,l2)%epsilon(3) = 0
 gos(l1,l2)%mu     (4) = 2
 gos(l1,l2)%alpha  (4) = 0
 gos(l1,l2)%beta   (4) = 1
 gos(l1,l2)%gamma  (4) = 1
 gos(l1,l2)%delta  (4) = 2
 gos(l1,l2)%epsilon(4) = 1

 ! new
 l1 = 2 ; l2 = 0
 call gos_allocate(l1,l2,3)
 gos(l1,l2)%mu     (1) = 1
 gos(l1,l2)%alpha  (1) = 2
 gos(l1,l2)%beta   (1) = 0
 gos(l1,l2)%gamma  (1) = 2
 gos(l1,l2)%delta  (1) = 0
 gos(l1,l2)%epsilon(1) = 0
 gos(l1,l2)%mu     (2) =-1
 gos(l1,l2)%alpha  (2) = 0
 gos(l1,l2)%beta   (2) = 0
 gos(l1,l2)%gamma  (2) = 2
 gos(l1,l2)%delta  (2) = 0
 gos(l1,l2)%epsilon(2) = 1
 gos(l1,l2)%mu     (3) = 1
 gos(l1,l2)%alpha  (3) = 0
 gos(l1,l2)%beta   (3) = 0
 gos(l1,l2)%gamma  (3) = 1
 gos(l1,l2)%delta  (3) = 0
 gos(l1,l2)%epsilon(3) = 0

 ! new
 l1 = 2 ; l2 = 1
 call gos_allocate(l1,l2,4)
 gos(l1,l2)%mu     (1) = 1
 gos(l1,l2)%alpha  (1) = 2
 gos(l1,l2)%beta   (1) = 1
 gos(l1,l2)%gamma  (1) = 2
 gos(l1,l2)%delta  (1) = 1
 gos(l1,l2)%epsilon(1) = 0
 gos(l1,l2)%mu     (2) = 2
 gos(l1,l2)%alpha  (2) = 1
 gos(l1,l2)%beta   (2) = 0
 gos(l1,l2)%gamma  (2) = 2
 gos(l1,l2)%delta  (2) = 1
 gos(l1,l2)%epsilon(2) = 1
 gos(l1,l2)%mu     (3) =-1
 gos(l1,l2)%alpha  (3) = 0
 gos(l1,l2)%beta   (3) = 1
 gos(l1,l2)%gamma  (3) = 2
 gos(l1,l2)%delta  (3) = 1
 gos(l1,l2)%epsilon(3) = 1
 gos(l1,l2)%mu     (4) = 1
 gos(l1,l2)%alpha  (4) = 0
 gos(l1,l2)%beta   (4) = 1
 gos(l1,l2)%gamma  (4) = 1
 gos(l1,l2)%delta  (4) = 1
 gos(l1,l2)%epsilon(4) = 0

 ! new
 l1 = 2 ; l2 = 2
 call gos_allocate(l1,l2,10)
 gos(l1,l2)%mu     ( 1) = 1
 gos(l1,l2)%alpha  ( 1) = 2
 gos(l1,l2)%beta   ( 1) = 2
 gos(l1,l2)%gamma  ( 1) = 2
 gos(l1,l2)%delta  ( 1) = 2
 gos(l1,l2)%epsilon( 1) = 0
 gos(l1,l2)%mu     ( 2) =-1
 gos(l1,l2)%alpha  ( 2) = 2
 gos(l1,l2)%beta   ( 2) = 0
 gos(l1,l2)%gamma  ( 2) = 2
 gos(l1,l2)%delta  ( 2) = 2
 gos(l1,l2)%epsilon( 2) = 1
 gos(l1,l2)%mu     ( 3) = 1
 gos(l1,l2)%alpha  ( 3) = 2
 gos(l1,l2)%beta   ( 3) = 0
 gos(l1,l2)%gamma  ( 3) = 2
 gos(l1,l2)%delta  ( 3) = 1
 gos(l1,l2)%epsilon( 3) = 0
 gos(l1,l2)%mu     ( 4) = 4
 gos(l1,l2)%alpha  ( 4) = 1
 gos(l1,l2)%beta   ( 4) = 1
 gos(l1,l2)%gamma  ( 4) = 2
 gos(l1,l2)%delta  ( 4) = 2
 gos(l1,l2)%epsilon( 4) = 1
 gos(l1,l2)%mu     ( 5) =-1
 gos(l1,l2)%alpha  ( 5) = 0
 gos(l1,l2)%beta   ( 5) = 2
 gos(l1,l2)%gamma  ( 5) = 2
 gos(l1,l2)%delta  ( 5) = 2
 gos(l1,l2)%epsilon( 5) = 1
 gos(l1,l2)%mu     ( 6) = 1
 gos(l1,l2)%alpha  ( 6) = 0
 gos(l1,l2)%beta   ( 6) = 2
 gos(l1,l2)%gamma  ( 6) = 1
 gos(l1,l2)%delta  ( 6) = 2
 gos(l1,l2)%epsilon( 6) = 0
 gos(l1,l2)%mu     ( 7) =-1
 gos(l1,l2)%alpha  ( 7) = 0
 gos(l1,l2)%beta   ( 7) = 0
 gos(l1,l2)%gamma  ( 7) = 2
 gos(l1,l2)%delta  ( 7) = 1
 gos(l1,l2)%epsilon( 7) = 1
 gos(l1,l2)%mu     ( 8) = 3
 gos(l1,l2)%alpha  ( 8) = 0
 gos(l1,l2)%beta   ( 8) = 0
 gos(l1,l2)%gamma  ( 8) = 2
 gos(l1,l2)%delta  ( 8) = 2
 gos(l1,l2)%epsilon( 8) = 2
 gos(l1,l2)%mu     ( 9) =-1
 gos(l1,l2)%alpha  ( 9) = 0
 gos(l1,l2)%beta   ( 9) = 0
 gos(l1,l2)%gamma  ( 9) = 1
 gos(l1,l2)%delta  ( 9) = 2
 gos(l1,l2)%epsilon( 9) = 1
 gos(l1,l2)%mu     (10) = 1
 gos(l1,l2)%alpha  (10) = 0
 gos(l1,l2)%beta   (10) = 0
 gos(l1,l2)%gamma  (10) = 1
 gos(l1,l2)%delta  (10) = 1
 gos(l1,l2)%epsilon(10) = 0

end subroutine setup_gos_llp


!=========================================================================
subroutine gos_allocate(l1,l2,np)
 implicit none
 integer,intent(in) :: l1,l2,np
!=====
 gos(l1,l2)%np = np
 allocate(gos(l1,l2)%mu     (np))
 allocate(gos(l1,l2)%alpha  (np))
 allocate(gos(l1,l2)%beta   (np))
 allocate(gos(l1,l2)%gamma  (np))
 allocate(gos(l1,l2)%delta  (np))
 allocate(gos(l1,l2)%epsilon(np))
end subroutine gos_allocate


!=========================================================================
end module m_gos
