!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This modules contains
! the MPI information and MPI basic functions
!
!=========================================================================
module m_mpi
 use m_definitions
 use m_warning,only: die
 use m_mpi_world
 use m_mpi_auxil
 use m_mpi_grid
 use m_mpi_ortho
 use m_mpi_local
 use m_mpi_trans
#ifdef HAVE_MPI
 use mpi
#endif


 logical,parameter :: parallel_grid      = .TRUE.
 logical,parameter :: parallel_auxil     = .TRUE.

 logical,protected :: parallel_buffer    = .TRUE.

#ifdef HAVE_SCALAPACK
 logical,parameter :: parallel_scalapack = .TRUE.
#else
 logical,parameter :: parallel_scalapack = .FALSE.
#endif

!===================================================
! MPI distribution
!  Example: nproc_ortho = 2 x  nproc_auxil = 8  = nproc_world = 16
!
! comm_world
!
! rank_auxil         0 |  1 |  2 |     |  7
! rank_ortho       ---------------------------
!      0             0 |  2 |  4 | ... | 14 |-> comm_auxil
!      1             1 |  3 |  5 | ... | 15 |-> comm_auxil
!                  ---------------------------
!                    | |    |    | ... |  | |
!                    v                    v
!                 comm_ortho           comm_ortho
!===================================================


 integer,protected :: iomaster = 0
 logical,protected :: is_iomaster = .TRUE.

 integer,private :: ngrid_mpi

 integer,allocatable,private :: task_proc(:)
 integer,allocatable,private :: ntask_proc(:)
 integer,allocatable,private :: task_number(:)

 integer,allocatable,private :: task_grid_proc(:)    ! index of the processor working for this grid point
 integer,allocatable,private :: ntask_grid_proc(:)   ! number of grid points for each procressor
 integer,allocatable,private :: task_grid_number(:)  ! local index of the grid point



contains


!=========================================================================
subroutine init_mpi_world()
 implicit none

!=====
 integer :: ier
!=====

#ifdef HAVE_MPI
 call MPI_INIT(ier)

 comm_world = MPI_COMM_WORLD
 call MPI_COMM_SIZE(comm_world,nproc_world,ier)
 call MPI_COMM_RANK(comm_world,rank_world,ier)

#else
 nproc_world = 1
 rank_world  = 0
#endif


#ifndef DEBUG
 if( rank_world /= iomaster ) then
   is_iomaster = .FALSE.
   close(stdout)
   open(unit=stdout,file='/dev/null')
 endif
#else
 if( rank_world /= iomaster ) then
   is_iomaster = .FALSE.
   call set_standard_output(2000+rank_world)
 endif
#endif


end subroutine init_mpi_world


!=========================================================================
subroutine init_mpi_other_communicators(nproc_ortho_in)
 implicit none

 integer,intent(in) :: nproc_ortho_in
!=====
 integer :: color
 integer :: ier
!=====

#ifdef HAVE_MPI

 nproc_ortho = nproc_ortho_in

 !
 ! Set up grid communicator
 !
! nproc_grid = nproc_world / nproc_ortho
!
! color = MODULO( rank_world , nproc_ortho )
! call MPI_COMM_SPLIT(comm_world,color,rank_world,comm_grid,ier);
!
! call MPI_COMM_SIZE(comm_grid,nproc_grid,ier)
! call MPI_COMM_RANK(comm_grid,rank_grid,ier)
! if( nproc_grid /= nproc_world / nproc_ortho ) then
!   write(stdout,*) rank_world,color,nproc_grid,nproc_world,nproc_ortho
!   call die('Problem in init_mpi')
! endif

 nproc_grid = nproc_world
 comm_grid = comm_world
 rank_grid = rank_world

 !
 ! Set up auxil communicator
 !
 nproc_auxil = nproc_world / nproc_ortho

 color = MODULO( rank_world , nproc_ortho )
 call MPI_COMM_SPLIT(comm_world,color,rank_world,comm_auxil,ier);

 call MPI_COMM_SIZE(comm_auxil,nproc_auxil,ier)
 call MPI_COMM_RANK(comm_auxil,rank_auxil,ier)
 if( nproc_auxil /= nproc_world / nproc_ortho ) then
   write(stdout,*) rank_world,color,nproc_auxil,nproc_world,nproc_ortho
   call die('Problem in init_mpi')
 endif

 !
 ! Set up ortho communicator
 !
 nproc_ortho = nproc_world / nproc_auxil

 color = rank_world / nproc_ortho
 call MPI_COMM_SPLIT(comm_world,color,rank_world,comm_ortho,ier);
 call MPI_COMM_RANK(comm_ortho,rank_ortho,ier)




#else
 nproc_ortho = 1
 rank_ortho  = 0
 nproc_auxil = 1
 rank_auxil  = 0
 nproc_grid  = 1
 rank_grid   = 0
#endif


#ifdef HAVE_MPI
  write(stdout,'(/,a)')       ' ==== MPI info'
  write(stdout,'(a50,1x,i6)')  'Number of proc:',nproc_world
  write(stdout,'(a50,1x,i6)')  'nproc_grid:    ',nproc_grid
  write(stdout,'(a50,1x,i6)')  'nproc_auxil:   ',nproc_auxil
  write(stdout,'(a50,1x,i6)')  'nproc_ortho:   ',nproc_ortho
  write(stdout,'(a50,1x,i6)')  'Master proc is:',iomaster
  write(stdout,'(a50,6x,l1)') 'Parallelize auxiliary basis:',parallel_auxil
  write(stdout,'(a50,6x,l1)')  'Parallelize XC grid points:',parallel_grid
  write(stdout,'(a50,6x,l1)')               'Use SCALAPACK:',parallel_scalapack
  write(stdout,'(a50,6x,l1)')    'Parallel using a buffer :',parallel_buffer
  write(stdout,'(/)')
#endif

end subroutine init_mpi_other_communicators


!=========================================================================
subroutine finish_mpi()
 implicit none
 integer :: ier
!=====

 if(ALLOCATED(task_proc))   deallocate(task_proc)
 if(ALLOCATED(ntask_proc))  deallocate(ntask_proc)
 if(ALLOCATED(task_number)) deallocate(task_number)

 if(ALLOCATED(task_grid_proc))   deallocate(task_grid_proc)
 if(ALLOCATED(ntask_grid_proc))  deallocate(ntask_grid_proc)
 if(ALLOCATED(task_grid_number)) deallocate(task_grid_number)

#ifdef HAVE_MPI
 call MPI_FINALIZE(ier)
#endif

end subroutine finish_mpi


!=========================================================================
subroutine init_dft_grid_distribution(ngrid)
 implicit none
 integer,intent(inout) :: ngrid
!=====

 ngrid_mpi = ngrid

 if( parallel_buffer ) then
   if( nproc_grid > 1 .AND. parallel_grid ) then
     write(stdout,'(/,a)') ' Initializing the distribution of the quadrature grid points'
   endif
 else
   if( nproc_local > 1 .AND. parallel_grid ) then
     write(stdout,'(/,a)') ' Initializing the distribution of the quadrature grid points'
   endif
 endif

 call distribute_grid_workload()

 if( parallel_buffer ) then
   ngrid = ntask_grid_proc(rank_grid)
 else
   ngrid = ntask_grid_proc(rank_local)
 endif

end subroutine init_dft_grid_distribution


!=========================================================================
subroutine destroy_dft_grid_distribution()
 implicit none
!=====

 if( ALLOCATED(task_grid_proc) )   deallocate(task_grid_proc)
 if( ALLOCATED(ntask_grid_proc) )  deallocate(ntask_grid_proc)
 if( ALLOCATED(task_grid_number) ) deallocate(task_grid_number)

end subroutine destroy_dft_grid_distribution


!=========================================================================
function is_my_grid_task(igrid)
 implicit none
 integer,intent(in) :: igrid
 logical            :: is_my_grid_task
!=====

 if( parallel_buffer ) then
   is_my_grid_task = ( rank_grid  == task_grid_proc(igrid) )
 else
   is_my_grid_task = ( rank_local == task_grid_proc(igrid) )
 endif

end function is_my_grid_task


!=========================================================================
subroutine distribute_grid_workload()
 implicit none
!=====
 integer            :: igrid,iproc_local
 integer            :: igrid_current
 integer            :: max_grid_per_proc
!=====


 if( .NOT. parallel_buffer ) then

   allocate(task_grid_proc(ngrid_mpi))
   allocate(ntask_grid_proc(0:nproc_local-1))
   allocate(task_grid_number(ngrid_mpi))

   if( parallel_grid) then

     write(stdout,'(/,a)') ' Distributing the grid among procs'

     ntask_grid_proc(:) = 0
     max_grid_per_proc = CEILING( DBLE(ngrid_mpi)/DBLE(nproc_local) )
     write(stdout,*) 'Maximum number of grid points for a single proc',max_grid_per_proc

     iproc_local=0
     do igrid=1,ngrid_mpi

       iproc_local = MODULO(igrid-1,nproc_local)

       !
       ! A simple check to avoid unexpected surprises
       if( iproc_local < 0 .OR. iproc_local >= nproc_local ) then
         call die('error in the distribution')
       endif

       task_grid_proc(igrid)        = iproc_local
       ntask_grid_proc(iproc_local) = ntask_grid_proc(iproc_local) + 1

     enddo

     task_grid_number(:)=0
     igrid_current=0
     do igrid=1,ngrid_mpi
       if( rank_local == task_grid_proc(igrid) ) then
         igrid_current = igrid_current + 1
         task_grid_number(igrid) = igrid_current
       endif
     enddo

   else
     !
     ! if parallel_grid is false,
     ! faking the code with trivial values
     ntask_grid_proc(:) = ngrid_mpi
     task_grid_proc(:)  = rank_local
     do igrid=1,ngrid_mpi
       task_grid_number(igrid) = igrid
     enddo

   endif

 else

   allocate(task_grid_proc(ngrid_mpi))
   allocate(ntask_grid_proc(0:nproc_grid-1))
   allocate(task_grid_number(ngrid_mpi))

   if( parallel_grid) then

     write(stdout,'(/,a)') ' Distributing the grid among procs'

     ntask_grid_proc(:) = 0
     max_grid_per_proc = CEILING( DBLE(ngrid_mpi)/DBLE(nproc_grid) )
     write(stdout,*) 'Maximum number of grid points for a single proc',max_grid_per_proc

     iproc_local=0
     do igrid=1,ngrid_mpi

       iproc_local = MODULO(igrid-1,nproc_grid)

       !
       ! A simple check to avoid unexpected surprises
       if( iproc_local < 0 .OR. iproc_local >= nproc_grid ) then
         call die('error in the distribution')
       endif

       task_grid_proc(igrid)        = iproc_local
       ntask_grid_proc(iproc_local) = ntask_grid_proc(iproc_local) + 1

     enddo

     task_grid_number(:)=0
     igrid_current=0
     do igrid=1,ngrid_mpi
       if( rank_grid == task_grid_proc(igrid) ) then
         igrid_current = igrid_current + 1
         task_grid_number(igrid) = igrid_current
       endif
     enddo

!     if( nproc_grid > 1 ) then
!       write(stdout,'(/,a)') ' Distribute work load among procs'
!       write(stdout,'(a,1x,f8.2)') ' Avg. tasks per cpu:',REAL(ngrid_mpi,dp) / REAL(nproc_grid,dp)
!       write(stdout,'(a,i6,a,i10)') ' proc # , grid points',rank_grid,' , ',ntask_grid_proc(rank_grid)
!     endif

   else
     !
     ! if parallel_grid is false,
     ! faking the code with trivial values
     ntask_grid_proc(:) = ngrid_mpi
     task_grid_proc(:)  = rank_grid
     do igrid=1,ngrid_mpi
       task_grid_number(igrid) = igrid
     enddo

   endif
 endif


end subroutine distribute_grid_workload




!=========================================================================
end module m_mpi
!=========================================================================
