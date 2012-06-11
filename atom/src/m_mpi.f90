!=========================================================================
#include "macros.h"
!=========================================================================
module m_mpi
 use m_definitions
#ifdef MPI
 use mpi
#endif

 private
! public  :: nproc,rank,ioproc
 public  :: init_mpi,finish_mpi,init_distribution,is_my_task,get_ntask,get_task_number,is_iomaster
 public  :: xsum

 integer :: nproc  = 1
 integer :: rank   = 0
 integer :: ioproc = 0

 integer :: mpi_comm

 integer :: nbf_mpi

 integer,allocatable :: task_proc(:)
 integer,allocatable :: ntask_proc(:)
 integer,allocatable :: task_number(:)

 interface xsum
   module procedure xsum_rrr
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

 if(nproc>1) then
   WRITE_MASTER(*,'(/,a)') ' Distribution initialized: 1 index distribution'
 endif

 ntask = index_prod(nbf,nbf)
 call distribute_workload(ntask)

end subroutine init_distribution

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
subroutine distribute_workload(ntask)
 implicit none
 integer,intent(in) :: ntask
!=====
 integer            :: itask,iproc
 integer            :: itask_current
!=====

 allocate(task_proc(ntask))
 allocate(ntask_proc(0:nproc-1))
 allocate(task_number(ntask))
 
 ntask_proc(:)=0
 do itask=1,ntask
!   iproc = MODULO(itask,nproc)
   iproc = FLOOR( itask / ( DBLE(ntask)/DBLE(nproc) ) -0.0001 )  ! This distribution should better preserve the shell structures
   task_proc(itask)  = iproc
   ntask_proc(iproc) = ntask_proc(iproc) + 1 
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
subroutine xsum_rrr(array)
 implicit none
 real(dp),intent(inout) :: array(:,:,:)
!=====
 integer :: n1,n2,n3
 integer :: ier
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

end subroutine xsum_rrr

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
