!=========================================================================
#include "macros.h"
!=========================================================================
module m_mpi
 use m_definitions
#ifdef MPI
 use mpi
#endif

 private
 public  :: init_mpi,finish_mpi,init_distribution,init_grid_distribution,get_ntask,get_task_number,is_iomaster
 public  :: xsum
 public  :: is_my_task
 public  :: is_my_grid_task
 public  :: parallel_grid,parallel_integral

 logical,parameter :: parallel_grid      = .TRUE.
 logical,parameter :: parallel_integral  = .FALSE.


 integer :: nproc  = 1
 integer :: rank   = 0
 integer :: ioproc = 0

 integer :: mpi_comm

 integer :: nbf_mpi
 integer :: ngrid_mpi

 integer,allocatable :: task_proc(:)
 integer,allocatable :: ntask_proc(:)
 integer,allocatable :: task_number(:)

 integer,allocatable :: task_grid_proc(:)    ! index of the processor working for this grid point
 integer,allocatable :: ntask_grid_proc(:)   ! number of grid points for each procressor
 integer,allocatable :: task_grid_number(:)  ! local index of the grid point

 interface xsum
   module procedure xsum_r
   module procedure xsum_ra1d
   module procedure xsum_ra3d
 end interface


contains


!=========================================================================
subroutine init_mpi()
 implicit none
 integer :: ier
!=====

#ifdef MPI
 call MPI_INIT(ier)
 mpi_comm = MPI_COMM_WORLD
#endif

 call get_size()
 call get_rank()

#ifdef MPI
  WRITE_MASTER(*,'(/,a)')      ' ==== MPI info'
  WRITE_MASTER(*,'(a50,x,i8)') 'Number of proc:',nproc
  WRITE_MASTER(*,'(a50,x,i8)') 'Master proc is:',ioproc
  WRITE_MASTER(*,'(a50,x,l1)') 'Parallelize Coulomb integrals:',parallel_integral
  WRITE_MASTER(*,'(a50,x,l1)') 'Parallelize XC grid points   :',parallel_grid
  WRITE_MASTER(*,'(/)')
#endif

end subroutine init_mpi


!=========================================================================
subroutine finish_mpi()
 implicit none
 integer :: ier
!=====

 if(allocated(task_proc))   deallocate(task_proc)
 if(allocated(ntask_proc))  deallocate(ntask_proc)
 if(allocated(task_number)) deallocate(task_number)

 if(allocated(task_grid_proc))   deallocate(task_grid_proc)
 if(allocated(ntask_grid_proc))  deallocate(ntask_grid_proc)
 if(allocated(task_grid_number)) deallocate(task_grid_number)

#ifdef MPI
 call MPI_FINALIZE(ier)
#endif

end subroutine finish_mpi


!=========================================================================
subroutine get_size
 implicit none
 integer :: ier=0
!=====

#ifdef MPI
 call MPI_COMM_SIZE(mpi_comm,nproc,ier)
#endif
 if(ier/=0) then
   WRITE_ME(*,*) 'error in get_size'
 endif

end subroutine get_size


!=========================================================================
subroutine get_rank
 implicit none
 integer :: ier=0
!=====

#ifdef MPI
 call MPI_COMM_RANK(mpi_comm,rank,ier)
#endif
 if(ier/=0) then
   WRITE_ME(*,*) 'error in get_rank'
 endif

end subroutine get_rank


!=========================================================================
subroutine init_distribution(nbf)
 implicit none
 integer,intent(in) :: nbf
!=====
 integer            :: ntask
!=====

 nbf_mpi = nbf

 if( nproc>1 .AND. parallel_integral ) then
   WRITE_MASTER(*,'(/,a)') ' Initializing distribution: 2-index distribution'
 endif

 ntask = index_prod(nbf,nbf)
 call distribute_workload(ntask)

end subroutine init_distribution


!=========================================================================
subroutine init_grid_distribution(ngrid)
 implicit none
 integer,intent(inout) :: ngrid
!=====

 ngrid_mpi = ngrid

 if( nproc>1 .AND. parallel_grid ) then
   WRITE_MASTER(*,'(/,a)') ' Initializing the distribution of the quadrature grid points'
 endif

 call distribute_grid_workload()

 ngrid = ntask_grid_proc(rank)

end subroutine init_grid_distribution


!=========================================================================
function is_iomaster()
 implicit none
 logical            :: is_iomaster
!=====
 
 is_iomaster = ( rank == ioproc )
 
end function is_iomaster


!=========================================================================
function is_my_task(ibf,jbf)
 implicit none
 integer,intent(in) :: ibf,jbf
 logical            :: is_my_task
!=====
 integer            :: task_number
!=====
 
 !
 ! Distribution among procs is performed on the third and fourth indexes of the ERIs
 !
 ! (kl|ij) = (K|I)       where I and K are composite index number that take into account
 !                       the symmetry I=(ij) = (ji)
 !
 ! The distribution is then performed on index I
 !

 task_number = index_prod(ibf,jbf)

 is_my_task = ( rank == task_proc(task_number) )
 
end function is_my_task


!=========================================================================
function is_my_grid_task(igrid)
 implicit none
 integer,intent(in) :: igrid
 logical            :: is_my_grid_task
!=====
 
 is_my_grid_task = ( rank == task_grid_proc(igrid) )
 
end function is_my_grid_task


!=========================================================================
subroutine distribute_grid_workload()
 implicit none
!=====
 integer            :: igrid,iproc
 integer            :: igrid_current
 integer            :: max_grid_per_proc
!=====

 allocate(task_grid_proc(ngrid_mpi))
 allocate(ntask_grid_proc(0:nproc-1))
 allocate(task_grid_number(ngrid_mpi))

 if( parallel_grid) then

   WRITE_MASTER(*,*) 
   WRITE_MASTER(*,*) 'Distributing the grid among procs'
   
   ntask_grid_proc(:)=0
   max_grid_per_proc = CEILING( DBLE(ngrid_mpi)/DBLE(nproc) )
   WRITE_MASTER(*,*) 'Maximum number of grid points for a single proc',max_grid_per_proc

   iproc=0
   do igrid=1,ngrid_mpi

     iproc = MODULO(igrid-1,nproc)

     !
     ! A simple check to avoid unexpected surprises
     if(iproc < 0 .OR. iproc >= nproc) then
       WRITE_MASTER(*,*) 'error in the distribution'
       stop'STOP'
     endif

     task_grid_proc(igrid)  = iproc
     ntask_grid_proc(iproc) = ntask_grid_proc(iproc) + 1 

   enddo

   task_grid_number(:)=0
   igrid_current=0
   do igrid=1,ngrid_mpi
     if( rank == task_grid_proc(igrid) ) then
       igrid_current = igrid_current + 1 
       task_grid_number(igrid) = igrid_current
     endif
   enddo

   if(nproc>1) then
     WRITE_MASTER(*,'(/,a)') ' Distribute work load among procs'
     WRITE_MASTER(*,'(a,x,f8.2)') ' Avg. tasks per cpu:',REAL(ngrid_mpi,dp)/REAL(nproc,dp)
     do iproc=0,nproc-1
       WRITE_MASTER(*,'(a,i6,a,i10)') ' proc # , grid points',iproc,' , ',ntask_grid_proc(iproc)
     enddo
   endif

 else
   !
   ! if parallel_grid is false,
   ! faking the code with trivial values
   ntask_grid_proc(:) = ngrid_mpi
   task_grid_proc(:)  = rank
   do igrid=1,ngrid_mpi
     task_grid_number(igrid) = igrid
   enddo

 endif

end subroutine distribute_grid_workload


!=========================================================================
subroutine distribute_workload(ntask)
 implicit none
 integer,intent(in) :: ntask
!=====
 integer            :: itask,iproc
 integer            :: itask_current
 integer            :: max_task_per_proc
!=====

 allocate(task_proc(ntask))
 allocate(ntask_proc(0:nproc-1))
 allocate(task_number(ntask))

 if( parallel_integral) then
   WRITE_MASTER(*,*) 
   WRITE_MASTER(*,*) 'Distributing the work load among procs'
   
   ntask_proc(:)=0
   max_task_per_proc = CEILING( DBLE(ntask)/DBLE(nproc) )
   WRITE_MASTER(*,*) 'Maximum number of tasks for a single proc',max_task_per_proc
   iproc=0
   do itask=1,ntask

!     iproc = MODULO(itask,nproc)
!     iproc = FLOOR( itask / ( DBLE(ntask)/DBLE(nproc) ) )  ! This distribution should better preserve the shell structures

     !
     ! A simple check to avoid unexpected surprises
     if(iproc < 0 .OR. iproc >= nproc) then
       WRITE_MASTER(*,*) 'error in the distribution'
       stop'STOP'
     endif

     task_proc(itask)  = iproc
     ntask_proc(iproc) = ntask_proc(iproc) + 1 
     if( ntask_proc(iproc) == max_task_per_proc ) then
       iproc = iproc + 1
     endif

   enddo

   task_number(:)=0
   itask_current=0
   do itask=1,ntask
     if( rank == task_proc(itask) ) then
       itask_current = itask_current + 1 
       task_number(itask) = itask_current
     endif
   enddo

   if(nproc>1) then
     WRITE_MASTER(*,'(/,a)') ' Distribute work load among procs'
     WRITE_MASTER(*,'(a,x,f8.2)') ' Avg. tasks per cpu:',REAL(ntask,dp)/REAL(nproc,dp)
     do iproc=0,nproc-1
       WRITE_MASTER(*,'(a,i6,a,i10)') ' proc # , tasks',iproc,' , ',ntask_proc(iproc)
     enddo
   endif

 else
   !
   ! if parallel_integral is false,
   ! faking the code with trivial values
   ntask_proc(:) = ntask
   task_proc(:)  = rank
   do itask=1,ntask
     task_number(itask) = itask
   enddo

 endif

end subroutine distribute_workload


!=========================================================================
function get_ntask()
 implicit none
 integer :: get_ntask
!=====

 get_ntask = ntask_proc(rank)

end function get_ntask


!=========================================================================
function get_task_number(ibf,jbf)
 implicit none
 integer,intent(in) :: ibf,jbf
!=====
 integer            :: itask
 integer            :: get_task_number
!=====

 itask = index_prod(ibf,jbf)
 get_task_number = task_number(itask)
 
 !
 ! Check
 !
 if(get_task_number == 0) then
   WRITE_ME(*,*) '=======',rank
   WRITE_ME(*,*) ibf,jbf,itask
   WRITE_ME(*,*) task_proc(itask)
   stop' *** That should not happen ***'
 endif

end function get_task_number


!=========================================================================
subroutine xsum_r(real_number)
 implicit none
 real(dp),intent(inout) :: real_number
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, real_number, n1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ier)
#endif
 if(ier/=0) then
   WRITE_ME(*,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_r


!=========================================================================
subroutine xsum_ra1d(array)
 implicit none
 real(dp),intent(inout) :: array(:)
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )

#ifdef MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ier)
#endif
 if(ier/=0) then
   WRITE_ME(*,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_ra1d


!=========================================================================
subroutine xsum_ra3d(array)
 implicit none
 real(dp),intent(inout) :: array(:,:,:)
!=====
 integer :: n1,n2,n3
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )
 n3 = SIZE( array, DIM=3 )

#ifdef MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ier)
#endif
 if(ier/=0) then
   WRITE_ME(*,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_ra3d


!=========================================================================
function index_prod(ibf,jbf)
 implicit none
 integer,intent(in) :: ibf,jbf
 integer            :: index_prod
!=====
 integer            :: jmin,imax
!=====

 imax=MAX(ibf,jbf)
 jmin=MIN(ibf,jbf)
 index_prod = (jmin-1)*nbf_mpi - (jmin-1)*(jmin-2)/2 + imax-jmin+1

end function index_prod


!=========================================================================
end module m_mpi
!=========================================================================
