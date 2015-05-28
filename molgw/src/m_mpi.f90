!=========================================================================
module m_mpi
 use m_definitions
#ifdef HAVE_MPI
 use mpi
#endif


 logical,parameter :: parallel_grid      = .TRUE.
 logical,parameter :: parallel_auxil     = .TRUE.
 logical,parameter :: parallel_integral  = .FALSE.   !FIXME remove this parameter
#ifdef HAVE_SCALAPACK
 logical,parameter :: parallel_scalapack = .TRUE.
#else
 logical,parameter :: parallel_scalapack = .FALSE.
#endif


 integer,protected :: nproc  = 1
 integer,protected :: rank   = 0
 integer,private   :: iomaster = 0
 logical,protected :: is_iomaster = .TRUE.

 integer,private :: mpi_comm

 integer,private :: nbf_mpi
 integer,private :: ngrid_mpi
 integer,private :: nocc_mp

! Parallelization on the auxiliary basis
 integer,allocatable,public :: iproc_ibf_auxil(:)
 integer,allocatable,public :: ibf_auxil_g(:)
 integer,allocatable,public :: ibf_auxil_l(:)
 integer,allocatable,public :: nbf_local_iproc(:)


 integer,allocatable,private :: task_proc(:)
 integer,allocatable,private :: ntask_proc(:)
 integer,allocatable,private :: task_number(:)

 integer,allocatable,private :: task_grid_proc(:)    ! index of the processor working for this grid point
 integer,allocatable,private :: ntask_grid_proc(:)   ! number of grid points for each procressor
 integer,allocatable,private :: task_grid_number(:)  ! local index of the grid point

 integer,allocatable,private :: task_fast_proc(:)    ! index of the processor working for this grid point
! integer,allocatable :: ntask_fast_proc(:)   ! number of grid points for each procressor
! integer,allocatable :: task_fast_number(:)  ! local index of the grid point

 interface xmin
   module procedure xmin_i
 end interface

 interface xmax
   module procedure xmax_i
 end interface

 interface xsum
   module procedure xsum_r
   module procedure xsum_ra1d
   module procedure xsum_ra2d
   module procedure xsum_ra3d
   module procedure xsum_ra4d
   module procedure xsum_ca4d
 end interface

 interface colindex_local_to_global
   module procedure colindex_local_to_global_distrib
   module procedure colindex_local_to_global_procindex
 end interface

 interface rowindex_local_to_global
   module procedure rowindex_local_to_global_distrib
   module procedure rowindex_local_to_global_procindex
 end interface


 !
 ! SCALAPACK variables
 !
 integer,parameter :: ndel=9
 integer,parameter :: block_col = 64
 integer,parameter :: block_row = 64
 integer,parameter :: first_row = 0
 integer,parameter :: first_col = 0

 integer,protected :: nproc_sca = 1
 integer,protected :: iproc_sca = 0

 ! SCALAPACK grid: square distribution
 integer,protected :: cntxt_sd
 integer,protected :: nprow_sd
 integer,protected :: npcol_sd
 integer,protected :: iprow_sd
 integer,protected :: ipcol_sd

 ! SCALAPACK grid: row distribution
 integer,protected :: cntxt_rd
 integer,protected :: nprow_rd
 integer,protected :: npcol_rd
 integer,protected :: iprow_rd
 integer,protected :: ipcol_rd

 ! SCALAPACK grid: column distribution
 integer,protected :: cntxt_cd
 integer,protected :: nprow_cd
 integer,protected :: npcol_cd
 integer,protected :: iprow_cd
 integer,protected :: ipcol_cd

 ! SCALAPACK grid: no distribution
 integer,protected :: cntxt_nd
 integer,protected :: nprow_nd
 integer,protected :: npcol_nd
 integer,protected :: iprow_nd
 integer,protected :: ipcol_nd

contains


!=========================================================================
subroutine init_mpi()
 implicit none
 integer :: ier
!=====

#ifdef HAVE_MPI
 call MPI_INIT(ier)
 mpi_comm = MPI_COMM_WORLD
#endif

 call get_size()
 call get_rank()

 if( rank /= iomaster ) then
   is_iomaster = .FALSE.
   close(stdout)
   open(unit=stdout,file='/dev/null')
 endif

#ifdef HAVE_MPI
  write(stdout,'(/,a)')       ' ==== MPI info'
  write(stdout,'(a50,x,i6)')  'Number of proc:',nproc
  write(stdout,'(a50,x,i6)')  'Master proc is:',iomaster
  write(stdout,'(a50,6x,l1)') 'Parallelize auxiliary basis:',parallel_auxil
  write(stdout,'(a50,6x,l1)')  'Parallelize XC grid points:',parallel_grid
  write(stdout,'(a50,6x,l1)')               'Use SCALAPACK:',parallel_scalapack
  write(stdout,'(/)')
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

#ifdef HAVE_MPI
 call MPI_FINALIZE(ier)
#endif

end subroutine finish_mpi


!=========================================================================
subroutine barrier()
 implicit none
 integer :: ier
!=====

#ifdef HAVE_MPI
 call MPI_BARRIER(mpi_comm,ier)
#endif

end subroutine barrier


!=========================================================================
subroutine get_size()
 implicit none
 integer :: ier=0
!=====

#ifdef HAVE_MPI
 call MPI_COMM_SIZE(mpi_comm,nproc,ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in get_size'
 endif

end subroutine get_size


!=========================================================================
subroutine get_rank()
 implicit none
 integer :: ier=0
!=====

#ifdef HAVE_MPI
 call MPI_COMM_RANK(mpi_comm,rank,ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in get_rank'
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
   write(stdout,'(/,a)') ' Initializing distribution: 2-index distribution'
 endif

 ntask = index_prod_mpi(nbf,nbf)
 call distribute_workload(ntask)

end subroutine init_distribution


!=========================================================================
subroutine init_grid_distribution(ngrid)
 implicit none
 integer,intent(inout) :: ngrid
!=====

 ngrid_mpi = ngrid

 if( nproc>1 .AND. parallel_grid ) then
   write(stdout,'(/,a)') ' Initializing the distribution of the quadrature grid points'
 endif

 call distribute_grid_workload()

 ngrid = ntask_grid_proc(rank)

end subroutine init_grid_distribution


!=========================================================================
subroutine destroy_grid_distribution()
 implicit none
!=====

 if( allocated(task_grid_proc) )   deallocate(task_grid_proc)
 if( allocated(ntask_grid_proc) )  deallocate(ntask_grid_proc)
 if( allocated(task_grid_number) ) deallocate(task_grid_number)

end subroutine destroy_grid_distribution


!=========================================================================
subroutine init_fast_distribution(ntask)
 implicit none
 integer,intent(in) :: ntask
!=====
 integer            :: itask
!=====


 write(stdout,'(/,a)') ' Initializing the distribution of the lowest level of MPI parallelization'

 allocate(task_fast_proc(ntask))
 do itask=1,ntask
   task_fast_proc(itask) = MODULO(itask-1,nproc)
 enddo


end subroutine init_fast_distribution


!=========================================================================
subroutine destroy_fast_distribution()
 implicit none
!=====

 write(stdout,'(/,a)') ' End of the lowest level of MPI parallelization'

 deallocate(task_fast_proc)

end subroutine destroy_fast_distribution


!=========================================================================
function is_my_fast_task(itask)
 implicit none
 integer,intent(in) :: itask
 logical            :: is_my_fast_task
!=====
 
 is_my_fast_task = ( rank == task_fast_proc(itask) )
 
end function is_my_fast_task


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

 task_number = index_prod_mpi(ibf,jbf)

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

   write(stdout,'(/,a)') ' Distributing the grid among procs'
   
   ntask_grid_proc(:)=0
   max_grid_per_proc = CEILING( DBLE(ngrid_mpi)/DBLE(nproc) )
   write(stdout,*) 'Maximum number of grid points for a single proc',max_grid_per_proc

   iproc=0
   do igrid=1,ngrid_mpi

     iproc = MODULO(igrid-1,nproc)

     !
     ! A simple check to avoid unexpected surprises
     if(iproc < 0 .OR. iproc >= nproc) then
       write(stdout,*) 'error in the distribution'
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
     write(stdout,'(/,a)') ' Distribute work load among procs'
     write(stdout,'(a,x,f8.2)') ' Avg. tasks per cpu:',REAL(ngrid_mpi,dp)/REAL(nproc,dp)
     do iproc=0,nproc-1
       write(stdout,'(a,i6,a,i10)') ' proc # , grid points',iproc,' , ',ntask_grid_proc(iproc)
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

   write(stdout,'(/,a)') ' Distributing the work load among procs'
   
   ntask_proc(:)=0
   max_task_per_proc = CEILING( DBLE(ntask)/DBLE(nproc) )
   write(stdout,*) 'Maximum number of tasks for a single proc',max_task_per_proc
   iproc=0
   do itask=1,ntask

!     iproc = MODULO(itask,nproc)
!     iproc = FLOOR( itask / ( DBLE(ntask)/DBLE(nproc) ) )  ! This distribution should better preserve the shell structures

     !
     ! A simple check to avoid unexpected surprises
     if(iproc < 0 .OR. iproc >= nproc) then
       write(stdout,*) 'error in the distribution'
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
     write(stdout,'(/,a)') ' Distribute work load among procs'
     write(stdout,'(a,x,f8.2)') ' Avg. tasks per cpu:',REAL(ntask,dp)/REAL(nproc,dp)
     do iproc=0,nproc-1
       write(stdout,'(a,i6,a,i10)') ' proc # , tasks',iproc,' , ',ntask_proc(iproc)
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

 itask = index_prod_mpi(ibf,jbf)
 get_task_number = task_number(itask)
 
 !
 ! Check
 !
 if(get_task_number == 0) then
   write(stdout,*) '=======',rank
   write(stdout,*) ibf,jbf,itask
   write(stdout,*) task_proc(itask)
   stop' *** That should not happen ***'
 endif

end function get_task_number


!=========================================================================
subroutine xmin_i(integer_number)
 implicit none
 integer,intent(inout) :: integer_number
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, integer_number, n1, MPI_INTEGER, MPI_MIN, mpi_comm, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xmin_i


!=========================================================================
subroutine xmax_i(integer_number)
 implicit none
 integer,intent(inout) :: integer_number
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, integer_number, n1, MPI_INTEGER, MPI_MAX, mpi_comm, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xmax_i


!=========================================================================
subroutine xsum_r(real_number)
 implicit none
 real(dp),intent(inout) :: real_number
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, real_number, n1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
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

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_ra1d


!=========================================================================
subroutine xsum_ra2d(array)
 implicit none
 real(dp),intent(inout) :: array(:,:)
!=====
 integer :: n1,n2
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_ra2d


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

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_ra3d


!=========================================================================
subroutine xsum_ra4d(array)
 implicit none
 real(dp),intent(inout) :: array(:,:,:,:)
!=====
 integer :: n1,n2,n3,n4
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )
 n3 = SIZE( array, DIM=3 )
 n4 = SIZE( array, DIM=4 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3*n4, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_ra4d


!=========================================================================
subroutine xsum_ca4d(array)
 implicit none
 complex(dp),intent(inout) :: array(:,:,:,:)
!=====
 integer :: n1,n2,n3,n4
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )
 n3 = SIZE( array, DIM=3 )
 n4 = SIZE( array, DIM=4 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3*n4, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_ca4d



!=========================================================================
function index_prod_mpi(ibf,jbf)
 implicit none
 integer,intent(in) :: ibf,jbf
 integer            :: index_prod_mpi
!=====
 integer            :: jmin,imax
!=====

 imax=MAX(ibf,jbf)
 jmin=MIN(ibf,jbf)
 index_prod_mpi = (jmin-1)*nbf_mpi - (jmin-1)*(jmin-2)/2 + imax-jmin+1

end function index_prod_mpi


!=========================================================================
subroutine init_scalapack()
 implicit none
!=====

#ifdef HAVE_SCALAPACK

 ! Get iproc_sca and nproc_sca
 call BLACS_PINFO( iproc_sca, nproc_sca )
 
 ! Set nprow, npcol

 ! Squared division of tasks
 nprow_sd = INT(SQRT(REAL(nproc_sca,dp)))
 npcol_sd = nproc_sca / nprow_sd
 if( nprow_sd * npcol_sd /= nproc_sca ) then
   write(stdout,'(a)') ' Some procs will be idling in the SCALAPACK distribution'
   write(stdout,'(a)') ' This is a waste and it is not yet coded anyway!'
   write(stdout,'(a)') ' => select a number of procs that is not a prime number'
   stop'not permitted'
 endif

 call BLACS_GET( -1, 0, cntxt_sd )
 call BLACS_GRIDINIT( cntxt_sd, 'R', nprow_sd, npcol_sd )
 call BLACS_GRIDINFO( cntxt_sd, nprow_sd, npcol_sd, iprow_sd, ipcol_sd )

 ! Row-only division of tasks
 nprow_rd = nproc_sca
 npcol_rd = 1
 call BLACS_GET( -1, 0, cntxt_rd )
 call BLACS_GRIDINIT( cntxt_rd, 'R', nprow_rd, npcol_rd )
 call BLACS_GRIDINFO( cntxt_rd, nprow_rd, npcol_rd, iprow_rd, ipcol_rd )

 ! Column-only division of tasks
 nprow_cd = 1
 npcol_cd = nproc_sca 
 call BLACS_GET( -1, 0, cntxt_cd )
 call BLACS_GRIDINIT( cntxt_cd, 'R', nprow_cd, npcol_cd )
 call BLACS_GRIDINFO( cntxt_cd, nprow_cd, npcol_cd, iprow_cd, ipcol_cd )

 ! No division of tasks
 nprow_nd = 1
 npcol_nd = 1
 call BLACS_GET( -1, 0, cntxt_nd )
 call BLACS_GRIDINIT( cntxt_nd, 'R', nprow_nd, npcol_nd )
 call BLACS_GRIDINFO( cntxt_nd, nprow_nd, npcol_nd, iprow_nd, ipcol_nd )

 write(stdout,'(/,a)')           ' ==== SCALAPACK info'
 write(stdout,'(a)')             '   Squared distribution'
 write(stdout,'(a50,x,i8)')      'Number of proc:',nprow_sd*npcol_sd
 write(stdout,'(a50,x,i8,x,i8)') 'Grid of procs:',nprow_sd,npcol_sd
! write(stdout,'(a)')             '       Row distribution'
! write(stdout,'(a50,x,i8)')      'Number of proc:',nprow_rd*npcol_rd
! write(stdout,'(a50,x,i8,x,i8)') 'Grid of procs:',nprow_rd,npcol_rd
! write(stdout,'(a)')             '    Column distribution'
! write(stdout,'(a50,x,i8)')      'Number of proc:',nprow_cd*npcol_cd
! write(stdout,'(a50,x,i8,x,i8)') 'Grid of procs:',nprow_cd,npcol_cd
 write(stdout,'(/)')
 
#else
 nprow_sd = 1
 npcol_sd = 1
 iprow_sd = 0
 ipcol_sd = 0
#endif

end subroutine init_scalapack


!=========================================================================
function row_block_size(mglobal,iprow,nprow)
 implicit none

 integer,intent(in) :: mglobal,iprow,nprow
 integer            :: row_block_size
#ifdef HAVE_SCALAPACK
 integer,external :: NUMROC 
#endif
!=====

#ifdef HAVE_SCALAPACK
 row_block_size = NUMROC(mglobal,block_row,iprow,first_row,nprow)
#else
 row_block_size = mglobal
#endif

end function row_block_size


!=========================================================================
function col_block_size(nglobal,ipcol,npcol)
 implicit none

 integer,intent(in) :: nglobal,ipcol,npcol
 integer            :: col_block_size
#ifdef HAVE_SCALAPACK
 integer,external :: NUMROC 
#endif
!=====

#ifdef HAVE_SCALAPACK
 col_block_size = NUMROC(nglobal,block_col,ipcol,first_col,npcol)
#else
 col_block_size = nglobal
#endif

end function col_block_size


!=========================================================================
subroutine init_desc(distribution,mglobal,nglobal,desc,mlocal,nlocal)
 implicit none
 character(len=1),intent(in) :: distribution
 integer,intent(in)          :: mglobal,nglobal
 integer,intent(out)         :: desc(ndel),mlocal,nlocal
!=====
 integer :: info
#ifdef HAVE_SCALAPACK
 integer,external :: NUMROC 
#endif
!=====

#ifdef HAVE_SCALAPACK
 select case(distribution)
 case('S')
   mlocal = NUMROC(mglobal,block_row,iprow_sd,first_row,nprow_sd)
   nlocal = NUMROC(nglobal,block_col,ipcol_sd,first_col,npcol_sd)
   call DESCINIT(desc,mglobal,nglobal,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,mlocal),info)
 case('R')
   mlocal = NUMROC(mglobal,block_row,iprow_rd,first_row,nprow_rd)
   nlocal = NUMROC(nglobal,block_col,ipcol_rd,first_col,npcol_rd)
   call DESCINIT(desc,mglobal,nglobal,block_row,block_col,first_row,first_col,cntxt_rd,MAX(1,mlocal),info)
 case('C')
   mlocal = NUMROC(mglobal,block_row,iprow_cd,first_row,nprow_cd)
   nlocal = NUMROC(nglobal,block_col,ipcol_cd,first_col,npcol_cd)
   call DESCINIT(desc,mglobal,nglobal,block_row,block_col,first_row,first_col,cntxt_cd,MAX(1,mlocal),info)
 case('N')
   mlocal = NUMROC(mglobal,block_row,iprow_nd,first_row,nprow_nd)
   nlocal = NUMROC(nglobal,block_col,ipcol_nd,first_col,npcol_nd)
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   stop'BUG'
 end select

 write(stdout,'(/,a,i6,a,i6,4x,i6)') ' SCALAPACK info: size of the local matrix for proc #', mlocal,' x ',nlocal,iproc_sca

#else
 desc(:)= 0
 desc(3)= mglobal
 desc(4)= nglobal
 desc(9)= mglobal
 mlocal = mglobal
 nlocal = nglobal
#endif

end subroutine init_desc



!=========================================================================
subroutine sum_sca_sd(real_number)
 implicit none
 real(dp),intent(inout) :: real_number
!=====

#ifdef HAVE_SCALAPACK
 call dgsum2d(cntxt_sd,'All',' ',1,1,real_number,1,-1,-1)
#endif

end subroutine sum_sca_sd


!=========================================================================
subroutine max_matrix_sca_sd(m,n,real_matrix)
 implicit none
 integer, intent(in)    :: m,n
 real(dp),intent(inout) :: real_matrix(m,n)
!=====
 integer :: idum1(1),idum2(1)
!=====

#ifdef HAVE_SCALAPACK
 call dgamx2d(cntxt_sd,'All',' ',m,n,real_matrix,m,-1,idum1,idum2,-1,-1)
#endif


end subroutine max_matrix_sca_sd


!=========================================================================
subroutine diagonalize_sca(desc,nglobal,mlocal,nlocal,matrix,eigval)
 implicit none
 integer,intent(in)     :: desc(ndel),nglobal,mlocal,nlocal
 real(dp),intent(inout) :: matrix(mlocal,nlocal)
 real(dp),intent(out)   :: eigval(nglobal)
!=====
 integer              :: desc_tmp(ndel)
 integer              :: lwork,info
 real(dp),allocatable :: work(:)
 real(dp)             :: eigvec(mlocal,nlocal)
!=====

#ifdef HAVE_SCALAPACK
 ! fake descriptor ! Why do I need this?
 desc_tmp(:) = desc(:)

 !
 ! First call to get the dimension of the array work
 lwork = -1
 allocate(work(1))
 call PDSYEV('N','U',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_tmp,work,lwork,info)
 lwork = NINT(work(1))
 deallocate(work)
 !
 ! Second call indeed perform the diago
 allocate(work(lwork))
 call PDSYEV('N','U',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_tmp,work,lwork,info)
 deallocate(work)

#else
 eigval(:) = 0.0_dp
#endif


end subroutine diagonalize_sca


!=========================================================================
function rowindex_global_to_local(distribution,iglobal)
 implicit none
 character(len=1),intent(in) :: distribution
 integer,intent(in)          :: iglobal
 integer                     :: rowindex_global_to_local
!=====
#ifdef HAVE_SCALAPACK
 integer,external :: INDXG2P,INDXG2L
#endif
!=====
 !
 ! returns the local index if this is proc in charge
 ! else returns 0 

#ifdef HAVE_SCALAPACK
 select case(distribution)
 case('S')
   if( iprow_sd == INDXG2P(iglobal,block_row,0,first_row,nprow_sd) ) then
     rowindex_global_to_local = INDXG2L(iglobal,block_row,0,first_row,nprow_sd)
   else
     rowindex_global_to_local = 0
   endif
 case('R')
   if( iprow_rd == INDXG2P(iglobal,block_row,0,first_row,nprow_rd) ) then
     rowindex_global_to_local = INDXG2L(iglobal,block_row,0,first_row,nprow_rd)
   else
     rowindex_global_to_local = 0
   endif
 case('C')
   if( iprow_cd == INDXG2P(iglobal,block_row,0,first_row,nprow_cd) ) then
     rowindex_global_to_local = INDXG2L(iglobal,block_row,0,first_row,nprow_cd)
   else
     rowindex_global_to_local = 0
   endif
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   stop'BUG'
 end select
#else
 rowindex_global_to_local = iglobal
#endif

end function rowindex_global_to_local


!=========================================================================
function colindex_global_to_local(distribution,iglobal)
 implicit none
 character(len=1),intent(in) :: distribution
 integer,intent(in)          :: iglobal
 integer                     :: colindex_global_to_local
!=====
#ifdef HAVE_SCALAPACK
 integer,external :: INDXG2P,INDXG2L
#endif
!=====
 !
 ! returns the local index if this is proc in charge
 ! else returns 0 

#ifdef HAVE_SCALAPACK
 select case(distribution)
 case('S')
   if( ipcol_sd == INDXG2P(iglobal,block_col,0,first_col,npcol_sd) ) then
     colindex_global_to_local = INDXG2L(iglobal,block_col,0,first_col,npcol_sd)
   else
     colindex_global_to_local = 0
   endif
 case('R')
   if( ipcol_rd == INDXG2P(iglobal,block_col,0,first_col,npcol_rd) ) then
     colindex_global_to_local = INDXG2L(iglobal,block_col,0,first_col,npcol_rd)
   else
     colindex_global_to_local = 0
   endif
 case('C')
   if( ipcol_cd == INDXG2P(iglobal,block_col,0,first_col,npcol_cd) ) then
     colindex_global_to_local = INDXG2L(iglobal,block_col,0,first_col,npcol_cd)
   else
     colindex_global_to_local = 0
   endif
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   stop'BUG'
 end select
#else
 colindex_global_to_local = iglobal
#endif

end function colindex_global_to_local


!=========================================================================
function rowindex_local_to_global_distrib(distribution,ilocal)
 implicit none
 character(len=1),intent(in) :: distribution
 integer,intent(in)          :: ilocal
 integer                     :: rowindex_local_to_global_distrib
!=====
#ifdef HAVE_SCALAPACK
 integer,external :: INDXL2G
#endif
!=====

#ifdef HAVE_SCALAPACK
 select case(distribution)
 case('S')
   rowindex_local_to_global_distrib = INDXL2G(ilocal,block_row,iprow_sd,first_row,nprow_sd)
 case('R')
   rowindex_local_to_global_distrib = INDXL2G(ilocal,block_row,iprow_rd,first_row,nprow_rd)
 case('C')
   rowindex_local_to_global_distrib = INDXL2G(ilocal,block_row,iprow_cd,first_row,nprow_cd)
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   stop'BUG'
 end select
#else
 rowindex_local_to_global_distrib = ilocal
#endif

end function rowindex_local_to_global_distrib


!=========================================================================
function rowindex_local_to_global_procindex(iprow,nprow,ilocal)
 implicit none
 integer,intent(in)          :: iprow,nprow,ilocal
 integer                     :: rowindex_local_to_global_procindex
!=====
#ifdef HAVE_SCALAPACK
 integer,external :: INDXL2G
#endif
!=====

#ifdef HAVE_SCALAPACK
 rowindex_local_to_global_procindex = INDXL2G(ilocal,block_row,iprow,first_row,nprow)
#else
 rowindex_local_to_global_procindex = ilocal
#endif

end function rowindex_local_to_global_procindex


!=========================================================================
function colindex_local_to_global_distrib(distribution,ilocal)
 implicit none
 character(len=1),intent(in) :: distribution
 integer,intent(in)          :: ilocal
 integer                     :: colindex_local_to_global_distrib
!=====
#ifdef HAVE_SCALAPACK
 integer,external :: INDXL2G
#endif
!=====

#ifdef HAVE_SCALAPACK
 select case(distribution)
 case('S')
   colindex_local_to_global_distrib = INDXL2G(ilocal,block_col,ipcol_sd,first_col,npcol_sd)
 case('R')
   colindex_local_to_global_distrib = INDXL2G(ilocal,block_col,ipcol_rd,first_col,npcol_rd)
 case('C')
   colindex_local_to_global_distrib = INDXL2G(ilocal,block_col,ipcol_cd,first_col,npcol_cd)
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   stop'BUG'
 end select
#else
 colindex_local_to_global_distrib = ilocal
#endif

end function colindex_local_to_global_distrib


!=========================================================================
function colindex_local_to_global_procindex(ipcol,npcol,ilocal)
 implicit none
 integer,intent(in)          :: ipcol,npcol,ilocal
 integer                     :: colindex_local_to_global_procindex
!=====
#ifdef HAVE_SCALAPACK
 integer,external :: INDXL2G
#endif
!=====

#ifdef HAVE_SCALAPACK
 colindex_local_to_global_procindex = INDXL2G(ilocal,block_col,ipcol,first_col,npcol)
#else
 colindex_local_to_global_procindex = ilocal
#endif

end function colindex_local_to_global_procindex



!=========================================================================
subroutine symmetrize_matrix(desc,m_mat,n_mat,mat)
 implicit none
 integer,intent(in)          :: desc(ndel)
 integer,intent(in)          :: m_mat,n_mat
 real(dp),intent(inout)      :: mat(m_mat,n_mat)
!=====
 integer                     :: imat,jmat,nmat
 real(dp)                    :: alpha
!=====
 
 ! On input, mat is lower triangular
 ! On output, mat is real symmetric

#ifdef HAVE_SCALAPACK
 if( desc(3)/=desc(4) ) stop'matrix has to be square'
 nmat = desc(3)
 do jmat=1,nmat
   do imat=jmat+1,nmat
     call PDELGET('All',' ',alpha,mat,imat,jmat,desc)
     call PDELSET(mat,jmat,imat,desc,alpha)
   enddo
 enddo
#else
 nmat=m_mat
 if (m_mat/=n_mat) stop'matrix should be squared if SCALAPACK is not used'
 if (m_mat/=desc(3)) then
  write(stderr,*) desc(3),m_mat
  stop'descriptor inconsistent row'
 endif
 if (n_mat/=desc(4)) stop'descriptor inconsistent col'
 do jmat=1,nmat
   do imat=jmat+1,nmat
     mat(jmat,imat) = mat(imat,jmat)
   enddo
 enddo
#endif

end subroutine symmetrize_matrix


!=========================================================================
subroutine finish_scalapack()
 implicit none
!=====

#ifdef HAVE_SCALAPACK
 call BLACS_GRIDEXIT( cntxt_sd )
 call BLACS_GRIDEXIT( cntxt_cd )
 call BLACS_GRIDEXIT( cntxt_rd )
 call BLACS_EXIT( 0 )
#endif

end subroutine finish_scalapack


!=========================================================================
end module m_mpi
!=========================================================================
