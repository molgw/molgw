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
  use m_mpi_tools
#if defined(HAVE_MPI)
  use mpi
#endif


  logical,parameter :: parallel_grid      = .TRUE.
  logical,parameter :: parallel_auxil     = .TRUE.

!===================================================
! MPI distribution
!  Example: ortho%nproc = 2 x  auxil%nproc = 8  = world%nproc = 16
!
! world%comm
!
! auxil%rank         0 |  1 |  2 |     |  7
! ortho%rank       ---------------------------
!      0             0 |  2 |  4 | ... | 14 |-> auxil%comm
!      1             1 |  3 |  5 | ... | 15 |-> auxil%comm
!                  ---------------------------
!                    | |    |    | ... |  | |
!                    v                    v
!                 ortho%comm           ortho%comm
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

  ! All the MPI communicators are here:
  type(mpi_communicator),protected :: world
  type(mpi_communicator),protected :: auxil
  type(mpi_communicator),protected :: ortho
  type(mpi_communicator),protected :: grid


contains


!=========================================================================
subroutine init_mpi_world()
  implicit none

  !=====
  integer :: ierror
  !=====

#if defined(HAVE_MPI)
  call MPI_INIT(ierror)
#endif

  call world%init(MPI_COMM_WORLD)

  if( world%rank /= iomaster ) then
    is_iomaster = .FALSE.
#if defined(DEBUG)
    call set_standard_output(2000+world%rank)
#else
    close(stdout)
    open(unit=stdout,file='/dev/null')
#endif
  endif

end subroutine init_mpi_world


!=========================================================================
subroutine init_mpi_other_communicators(nproc_ortho_in)
  implicit none

  integer,intent(in) :: nproc_ortho_in
 !=====
  integer :: color
  integer :: ier
  !=====

#if defined(HAVE_MPI)

  ortho%nproc = nproc_ortho_in

  !
  ! Set up grid communicator
  !
  call grid%init(world%comm)

  !
  ! Set up auxil communicator
  !
  auxil%nproc = world%nproc / ortho%nproc

  color = MODULO( world%rank , ortho%nproc )
  call MPI_COMM_SPLIT(world%comm,color,world%rank,auxil%comm,ier);

  call auxil%init(auxil%comm)

  if( auxil%nproc /= world%nproc / ortho%nproc ) then
    write(stdout,*) world%rank,color,auxil%nproc,world%nproc,ortho%nproc
    call die('Problem in init_mpi')
  endif

  !
  ! Set up ortho communicator
  !
  ortho%nproc = world%nproc / auxil%nproc

  color = world%rank / ortho%nproc
  call MPI_COMM_SPLIT(world%comm,color,world%rank,ortho%comm,ier);

  call ortho%init(ortho%comm)


#else
  ortho%rank  = 0
  ortho%nproc = 1
  auxil%nproc = 1
  auxil%rank  = 0
  grid%nproc  = 1
  grid%rank   = 0
#endif


#if defined(HAVE_MPI)
  write(stdout,'(/,a)')       ' ==== MPI info'
  write(stdout,'(a50,1x,i6)')  'Number of proc:',world%nproc
  write(stdout,'(a50,1x,i6)')  'grid%nproc:    ',grid%nproc
  write(stdout,'(a50,1x,i6)')  'auxil%nproc:   ',auxil%nproc
  write(stdout,'(a50,1x,i6)')  'ortho%nproc:   ',ortho%nproc
  write(stdout,'(a50,1x,i6)')  'Master proc is:',iomaster
  write(stdout,'(a50,6x,l1)') 'Parallelize auxiliary basis:',parallel_auxil
  write(stdout,'(a50,6x,l1)')  'Parallelize XC grid points:',parallel_grid
#if defined(HAVE_SCALAPACK)
  write(stdout,'(a50,6x,l1)')               'Use SCALAPACK:',.TRUE.
#else
  write(stdout,'(a50,6x,l1)')               'Use SCALAPACK:',.FALSE.
#endif
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

#if defined(HAVE_MPI)
 call MPI_FINALIZE(ier)
#endif

end subroutine finish_mpi


!=========================================================================
subroutine init_dft_grid_distribution(ngrid)
 implicit none
 integer,intent(inout) :: ngrid
!=====

 ngrid_mpi = ngrid

 if( grid%nproc > 1 .AND. parallel_grid ) then
   write(stdout,'(/,a)') ' Initializing the distribution of the quadrature grid points'
 endif

 call distribute_grid_workload()

 ngrid = ntask_grid_proc(grid%rank)

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

 is_my_grid_task = ( grid%rank  == task_grid_proc(igrid) )

end function is_my_grid_task


!=========================================================================
subroutine distribute_grid_workload()
 implicit none
!=====
 integer            :: igrid,iproc_local
 integer            :: igrid_current
 integer            :: max_grid_per_proc
!=====


 allocate(task_grid_proc(ngrid_mpi))
 allocate(ntask_grid_proc(0:grid%nproc-1))
 allocate(task_grid_number(ngrid_mpi))

 if( parallel_grid) then

   write(stdout,'(/,a)') ' Distributing the grid among procs'

   ntask_grid_proc(:) = 0
   max_grid_per_proc = CEILING( DBLE(ngrid_mpi)/DBLE(grid%nproc) )
   write(stdout,*) 'Maximum number of grid points for a single proc',max_grid_per_proc

   iproc_local=0
   do igrid=1,ngrid_mpi

     iproc_local = MODULO(igrid-1,grid%nproc)

     !
     ! A simple check to avoid unexpected surprises
     if( iproc_local < 0 .OR. iproc_local >= grid%nproc ) then
       call die('error in the distribution')
     endif

     task_grid_proc(igrid)        = iproc_local
     ntask_grid_proc(iproc_local) = ntask_grid_proc(iproc_local) + 1

   enddo

   task_grid_number(:)=0
   igrid_current=0
   do igrid=1,ngrid_mpi
     if( grid%rank == task_grid_proc(igrid) ) then
       igrid_current = igrid_current + 1
       task_grid_number(igrid) = igrid_current
     endif
   enddo

 else
   !
   ! if parallel_grid is false,
   ! faking the code with trivial values
   ntask_grid_proc(:) = ngrid_mpi
   task_grid_proc(:)  = grid%rank
   do igrid=1,ngrid_mpi
     task_grid_number(igrid) = igrid
   enddo

 endif


end subroutine distribute_grid_workload



!=========================================================================
end module m_mpi
!=========================================================================
