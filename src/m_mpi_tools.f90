!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the high-level MPI routines for parallelization over all procs (=world)
!
!=========================================================================
#include "molgw.h"
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
    ! min
    generic :: min  => mpic_min_dp
    generic :: min  => mpic_min_i
    procedure :: mpic_min_dp
    procedure :: mpic_min_i
    ! max
    generic :: max  => mpic_max_dp
    generic :: max  => mpic_max_i
    procedure :: mpic_max_dp
    procedure :: mpic_max_i
    ! and
    procedure :: and => mpic_and_l
    ! broadcast
    generic :: bcast  => mpic_bcast_dp
    generic :: bcast  => mpic_bcast_cdp
    generic :: bcast  => mpic_bcast_i
    procedure :: mpic_bcast_dp
    procedure :: mpic_bcast_cdp
    procedure :: mpic_bcast_i
  end type mpi_communicator



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
  integer(kind=int8),intent(inout) :: array(..)
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
subroutine mpic_max_dp(mpic,array)
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
  call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_PRECISION, MPI_MAX, mpic%comm, ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_ALLREDUCE'
  endif

end subroutine mpic_max_dp


!=========================================================================
subroutine mpic_max_i(mpic,array)
  implicit none
  class(mpi_communicator),intent(in) :: mpic
  integer,intent(inout) :: array(..)
  !=====
  integer :: nsize
  integer :: ierror=0
  !=====

  if( mpic%nproc == 1 ) return

  nsize = SIZE(array)

#if defined(HAVE_MPI)
  call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_INTEGER, MPI_MAX, mpic%comm, ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_ALLREDUCE'
  endif

end subroutine mpic_max_i


!=========================================================================
subroutine mpic_min_dp(mpic,array)
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
  call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_PRECISION, MPI_MIN, mpic%comm, ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_ALLREDUCE'
  endif

end subroutine mpic_min_dp


!=========================================================================
subroutine mpic_min_i(mpic,array)
  implicit none
  class(mpi_communicator),intent(in) :: mpic
  integer,intent(inout) :: array(..)
  !=====
  integer :: nsize
  integer :: ierror=0
  !=====

  if( mpic%nproc == 1 ) return

  nsize = SIZE(array)

#if defined(HAVE_MPI)
  call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_INTEGER, MPI_MIN, mpic%comm, ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_ALLREDUCE'
  endif

end subroutine mpic_min_i


!=========================================================================
subroutine mpic_and_l(mpic,array)
  implicit none
  class(mpi_communicator),intent(in) :: mpic
  logical,intent(inout) :: array(..)
  !=====
  integer :: nsize
  integer :: ierror=0
  !=====

  if( mpic%nproc == 1 ) return

  nsize = SIZE(array)

#if defined(HAVE_MPI)
  call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_LOGICAL, MPI_LAND, mpic%comm, ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_ALLREDUCE'
  endif

end subroutine mpic_and_l


!=========================================================================
subroutine mpic_bcast_i(mpic,rank,array)
  implicit none
  class(mpi_communicator),intent(in) :: mpic
  integer,intent(in)    :: rank
  integer,intent(inout) :: array(:)
  !=====
  integer :: nsize
  integer :: ierror=0
  !=====

  if( mpic%nproc == 1 ) return

  nsize = SIZE(array)

#if defined(HAVE_MPI)
  call MPI_BCAST(array,nsize,MPI_INTEGER,rank,mpic%comm,ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_BCAST'
  endif

end subroutine mpic_bcast_i


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
  call MPI_BCAST(array,nsize,MPI_DOUBLE_PRECISION,rank,mpic%comm,ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_BCAST'
  endif

end subroutine mpic_bcast_dp


!=========================================================================
subroutine mpic_bcast_cdp(mpic,rank,array)
  implicit none
  class(mpi_communicator),intent(in) :: mpic
  integer,intent(in)     :: rank
  complex(dp),intent(inout) :: array(..)
  !=====
  integer :: nsize
  integer :: ierror=0
  !=====

  if( mpic%nproc == 1 ) return

  nsize = SIZE(array)

#if defined(HAVE_MPI)
  call MPI_BCAST(array,nsize,MPI_DOUBLE_COMPLEX,rank,mpic%comm,ierror)
#endif
  if( ierror /= 0 ) then
    write(stdout,*) 'error in MPI_BCAST'
  endif

end subroutine mpic_bcast_cdp



end module m_mpi_tools


!=========================================================================
