!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the high-level MPI routines for parallelization over all procs (=world)
!
!=========================================================================
module m_mpi_tools
  use m_definitions
  use m_warning,only: die
#if defined(HAVE_MPI)
!!!#include<mpif.h>
  use mpi
#endif


  type mpi_communicator
    integer    :: comm       ! MPI communicator
    integer    :: nproc      ! number of procs in the communicator comm
    integer    :: rank       ! index           in the communicator comm
  contains
    procedure :: init => mpic_init
    procedure :: barrier => mpic_barrier
    ! sum
    generic :: sum  => mpic_sum_dp
    generic :: sum  => mpic_sum_cdp
    generic :: sum  => mpic_sum_i8
    procedure :: mpic_sum_dp
    procedure :: mpic_sum_cdp
    procedure :: mpic_sum_i8
    ! broadcast
    generic :: bcast  => mpic_bcast_dp
    procedure :: mpic_bcast_dp
  end type mpi_communicator


#if 0
  !
 ! Interfaces for high-level MPI reduce operations
 ! "world" series
 !
 interface xmin_world
   module procedure xmin_world_r
 end interface

 interface xmax_world
   module procedure xmax_world_i
   module procedure xmax_world_r
   module procedure xmax_world_ia2d
   module procedure xmax_world_ra1d
 end interface

 interface xbcast_world
   module procedure xbcast_world_ra1d
   module procedure xbcast_world_ra2d
   module procedure xbcast_world_ra3d
   module procedure xbcast_world_ca2d
 end interface

 interface xand_world
   module procedure xand_world_l
   module procedure xand_world_la1d
   module procedure xand_world_la2d
 end interface
#endif


contains


!=========================================================================
subroutine mpic_init(mpic,comm_in)
  implicit none

  class(mpi_communicator),intent(inout) :: mpic
  integer,intent(in)                    :: comm_in
  !=====
  integer :: ierror
  !=====

#if defined(HAVE_MPI)
  mpic%comm = comm_in
  call MPI_COMM_SIZE(mpic%comm,mpic%nproc,ierror)
  call MPI_COMM_RANK(mpic%comm,mpic%rank,ierror)
#else
  mpic%comm  = 1
  mpic%nproc = 1
  mpic%rank  = 0
#endif

end subroutine mpic_init


!=========================================================================
subroutine mpic_barrier(mpic)
  implicit none

  class(mpi_communicator),intent(in) :: mpic
  !=====
  integer :: ierror
  !=====

#if defined(HAVE_MPI)
 call MPI_BARRIER(mpic%comm,ierror)
#endif

end subroutine mpic_barrier


!=========================================================================
subroutine sum_dp1d(mpic,array)
  implicit none
  class(mpi_communicator),intent(in) :: mpic
  real(dp),intent(inout) :: array(:)
  !=====
  integer :: n
  integer :: ierror=0
  !=====

  if( mpic%nproc == 1 ) return

  n = SIZE(array)

#if defined(HAVE_MPI)
  call MPI_ALLREDUCE( MPI_IN_PLACE, array, n, MPI_DOUBLE_PRECISION, MPI_SUM, mpic%comm, ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_ALLREDUCE'
  endif

end subroutine sum_dp1d


!=========================================================================
subroutine mpic_sum_dp(mpic,array)
  implicit none
  class(mpi_communicator),intent(in) :: mpic
  real(dp),intent(inout) :: array(..)
  !=====
  integer :: nsize
  integer :: ierror=0
  !=====

  if( mpic%nproc == 1 ) return

  nsize = SIZE(array)

#if defined(HAVE_MPI)
  call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_PRECISION, MPI_SUM, mpic%comm, ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_ALLREDUCE'
  endif

end subroutine mpic_sum_dp


!=========================================================================
subroutine mpic_sum_cdp(mpic,array)
  implicit none
  class(mpi_communicator),intent(in) :: mpic
  complex(dp),intent(inout) :: array(..)
  !=====
  integer :: nsize
  integer :: ierror=0
  !=====

  if( mpic%nproc == 1 ) return

  nsize = SIZE(array)

#if defined(HAVE_MPI)
  call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_COMPLEX, MPI_SUM, mpic%comm, ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_ALLREDUCE'
  endif

end subroutine mpic_sum_cdp


!=========================================================================
subroutine mpic_sum_i8(mpic,array)
  implicit none
  class(mpi_communicator),intent(in) :: mpic
  complex(dp),intent(inout) :: array(..)
  !=====
  integer :: nsize
  integer :: ierror=0
  !=====

  if( mpic%nproc == 1 ) return

  nsize = SIZE(array)

#if defined(HAVE_MPI)
  call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_INTEGER8, MPI_SUM, mpic%comm, ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_ALLREDUCE'
  endif

end subroutine mpic_sum_i8


!=========================================================================
subroutine mpic_bcast_dp(mpic,rank,array)
  implicit none
  class(mpi_communicator),intent(in) :: mpic
  integer,intent(in)     :: rank
  real(dp),intent(inout) :: array(..)
  !=====
  integer :: nsize
  integer :: ierror=0
  !=====

  if( mpic%nproc == 1 ) return

  nsize = SIZE(array)

#if defined(HAVE_MPI)
  call MPI_BCAST(array,n1,MPI_DOUBLE_PRECISION,rank,mpic%comm,ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_BCAST'
  endif

end subroutine mpic_bcast_dp





end module m_mpi_tools


!=========================================================================
