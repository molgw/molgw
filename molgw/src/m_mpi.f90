!=========================================================================
! This file is part of MOLGW.
!=========================================================================
module m_mpi
 use m_definitions
 use m_warning,only: die
#ifdef HAVE_MPI
 use mpi
#endif


 logical,parameter :: parallel_grid      = .TRUE.
 logical,parameter :: parallel_auxil     = .TRUE.

 logical,protected :: parallel_ham       = .FALSE.

#ifdef HAVE_SCALAPACK
 logical,parameter :: parallel_scalapack = .TRUE.
#else
 logical,parameter :: parallel_scalapack = .FALSE.
#endif

 integer,parameter :: SCALAPACK_MIN = 400  ! TODO Better tune this parameter in the future

 integer,private   :: comm_world
 integer,protected :: nproc  = 1
 integer,protected :: rank   = 0
 integer,private   :: iomaster = 0
 logical,protected :: is_iomaster = .TRUE.

 integer,private   :: comm_local
 integer,protected :: nproc_local  = 1
 integer,protected :: rank_local   = 0

 integer,protected :: comm_trans
 integer,protected :: nproc_trans  = 1
 integer,protected :: rank_trans   = 0

 integer,private :: nbf_mpi
 integer,private :: ngrid_mpi
 integer,private :: nocc_mp

! Parallelization on the auxiliary basis
 integer,allocatable,protected :: iproc_ibf_auxil(:)
 integer,allocatable,protected :: ibf_auxil_g(:)
 integer,allocatable,protected :: ibf_auxil_l(:)
 integer,allocatable,protected :: nbf_local_iproc(:)

! Parallelization on the auxiliary basis (LR part)
 integer,allocatable,public :: iproc_ibf_auxil_lr(:)
 integer,allocatable,public :: ibf_auxil_g_lr(:)
 integer,allocatable,public :: ibf_auxil_l_lr(:)
 integer,allocatable,public :: nbf_local_iproc_lr(:)

 integer,allocatable,protected :: rank_sca_to_mpi(:,:)
 integer,allocatable,protected :: rank_ham_sca_to_mpi(:,:)

 integer,allocatable,private :: task_proc(:)
 integer,allocatable,private :: ntask_proc(:)
 integer,allocatable,private :: task_number(:)

 integer,allocatable,private :: task_grid_proc(:)    ! index of the processor working for this grid point
 integer,allocatable,private :: ntask_grid_proc(:)   ! number of grid points for each procressor
 integer,allocatable,private :: task_grid_number(:)  ! local index of the grid point

 integer,allocatable,private :: task_fast_proc(:)    ! index of the processor working for this grid point
! integer,allocatable :: ntask_fast_proc(:)   ! number of grid points for each procressor
! integer,allocatable :: task_fast_number(:)  ! local index of the grid point

 interface xbcast
   module procedure xbcast_ra2d
 end interface

 interface xand
   module procedure xand_l
   module procedure xand_la1d
   module procedure xand_la2d
 end interface

 interface xmin
   module procedure xmin_i
 end interface

 interface xmax
   module procedure xmax_i
   module procedure xmax_r
   module procedure xmax_ia2d
   module procedure xmax_ra1d
 end interface

 interface xlocal_max
   module procedure xlocal_max_i
 end interface

 interface xtrans_max
   module procedure xtrans_max_ia2d
 end interface

 interface xsum
   module procedure xsum_r
   module procedure xsum_ra1d
   module procedure xsum_ra2d
   module procedure xsum_ra3d
   module procedure xsum_ra4d
   module procedure xsum_ca1d
   module procedure xsum_ca2d
   module procedure xsum_ca4d
   module procedure xsum_procindex_ra2d
 end interface

 interface xlocal_sum
   module procedure xlocal_sum_ra2d
   module procedure xlocal_sum_ra3d
 end interface

 interface xtrans_sum
   module procedure xtrans_sum_r
   module procedure xtrans_sum_ra1d
   module procedure xtrans_sum_ra2d
   module procedure xtrans_sum_ra3d
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
 integer,parameter :: block_col = 32
 integer,parameter :: block_row = 32
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

 ! SCALAPACK grid: hamiltonian
 integer,protected :: cntxt_ham
 integer,protected :: nprow_ham
 integer,protected :: npcol_ham
 integer,protected :: iprow_ham
 integer,protected :: ipcol_ham
 integer,protected :: desc_ham(ndel)
 integer,public    :: desc_c(ndel)
 integer,public    :: desc_small(ndel)

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
 comm_world = MPI_COMM_WORLD

 call MPI_COMM_SIZE(comm_world,nproc,ier)
 call MPI_COMM_RANK(comm_world,rank,ier)
#endif

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
subroutine barrier()
 implicit none
 integer :: ier
!=====

#ifdef HAVE_MPI
 call MPI_BARRIER(comm_world,ier)
#endif

end subroutine barrier


!=========================================================================
subroutine distribute_auxil_basis(nbf_auxil_basis,nbf_auxil_basis_local)
 implicit none

 integer,intent(in)  :: nbf_auxil_basis
 integer,intent(out) :: nbf_auxil_basis_local
!=====
 integer :: ibf
 integer :: ibf_local
 integer :: iproc
!=====

#ifdef TODAY
 ! Use nproc instead nproc_local

 allocate(iproc_ibf_auxil(nbf_auxil_basis))
 allocate(nbf_local_iproc(0:nproc-1))

 iproc              = nproc - 1
 nbf_local_iproc(:) = 0
 do ibf=1,nbf_auxil_basis

   iproc = MODULO(iproc+1,nproc)

   iproc_ibf_auxil(ibf) = iproc

   nbf_local_iproc(iproc) = nbf_local_iproc(iproc) + 1

 enddo

 nbf_auxil_basis_local = nbf_local_iproc(rank)

 allocate(ibf_auxil_g(nbf_auxil_basis_local))
 allocate(ibf_auxil_l(nbf_auxil_basis))
 ibf_local = 0
 do ibf=1,nbf_auxil_basis
   if( rank == iproc_ibf_auxil(ibf) ) then
     ibf_local = ibf_local + 1
     ibf_auxil_g(ibf_local) = ibf
     ibf_auxil_l(ibf)       = ibf_local
   endif
 enddo

 write(stdout,'(/,a)') ' Distribute auxiliary basis functions among processors'
 do iproc=0,0
   write(stdout,'(a,i4,a,i6,a)')   ' Proc: ',iproc,' treats ',nbf_local_iproc(iproc),' auxiliary basis functions'
 enddo

#else

 allocate(iproc_ibf_auxil(nbf_auxil_basis))
 allocate(nbf_local_iproc(0:nproc_local-1))

 iproc              = nproc_local-1
 nbf_local_iproc(:) = 0
 do ibf=1,nbf_auxil_basis

   iproc = MODULO(iproc+1,nproc_local)

   iproc_ibf_auxil(ibf) = iproc

   nbf_local_iproc(iproc) = nbf_local_iproc(iproc) + 1

 enddo

 nbf_auxil_basis_local = nbf_local_iproc(rank_local)

 allocate(ibf_auxil_g(nbf_auxil_basis_local))
 allocate(ibf_auxil_l(nbf_auxil_basis))
 ibf_local = 0
 do ibf=1,nbf_auxil_basis
   if( rank_local == iproc_ibf_auxil(ibf) ) then
     ibf_local = ibf_local + 1
     ibf_auxil_g(ibf_local) = ibf
     ibf_auxil_l(ibf)       = ibf_local
   endif
 enddo

 write(stdout,'(/,a)') ' Distribute auxiliary basis functions among processors'
 do iproc=0,0
   write(stdout,'(a,i4,a,i6,a)')   ' Proc: ',iproc,' treats ',nbf_local_iproc(iproc),' auxiliary basis functions'
 enddo

#endif


end subroutine distribute_auxil_basis


!=========================================================================
subroutine distribute_auxil_basis_lr(nbf_auxil_basis,nbf_auxil_basis_local)
 implicit none

 integer,intent(in)  :: nbf_auxil_basis
 integer,intent(out) :: nbf_auxil_basis_local
!=====
 integer :: ibf
 integer :: ibf_local
 integer :: iproc
!=====

 allocate(iproc_ibf_auxil_lr(nbf_auxil_basis))
 allocate(nbf_local_iproc_lr(0:nproc-1))

 iproc = nproc-1
 nbf_local_iproc_lr(:) = 0
 do ibf=1,nbf_auxil_basis

   iproc = MODULO(iproc+1,nproc)

   iproc_ibf_auxil_lr(ibf) = iproc

   nbf_local_iproc_lr(iproc) = nbf_local_iproc_lr(iproc) + 1

 enddo

 nbf_auxil_basis_local = nbf_local_iproc_lr(rank)

 allocate(ibf_auxil_g_lr(nbf_auxil_basis_local))
 allocate(ibf_auxil_l_lr(nbf_auxil_basis))
 ibf_local = 0
 do ibf=1,nbf_auxil_basis
   if( rank == iproc_ibf_auxil_lr(ibf) ) then
     ibf_local = ibf_local + 1
     ibf_auxil_g_lr(ibf_local) = ibf
     ibf_auxil_l_lr(ibf)       = ibf_local
   endif
 enddo

 write(stdout,'(/,a)') ' Distribute LR auxiliary basis functions among processors'
 do iproc=0,0
   write(stdout,'(a,i4,a,i6,a)')   ' Proc: ',iproc,' treats ',nbf_local_iproc_lr(iproc),' auxiliary basis functions'
 enddo


end subroutine distribute_auxil_basis_lr


!=========================================================================
subroutine init_grid_distribution(ngrid)
 implicit none
 integer,intent(inout) :: ngrid
!=====

 ngrid_mpi = ngrid

#ifndef TODAY
 if( nproc_local > 1 .AND. parallel_grid ) then
   write(stdout,'(/,a)') ' Initializing the distribution of the quadrature grid points'
 endif
#else
 if( nproc > 1 .AND. parallel_grid ) then
   write(stdout,'(/,a)') ' Initializing the distribution of the quadrature grid points'
 endif
#endif

 call distribute_grid_workload()

#ifndef TODAY
 ngrid = ntask_grid_proc(rank_local)
#else
 ngrid = ntask_grid_proc(rank)
#endif

end subroutine init_grid_distribution


!=========================================================================
subroutine destroy_grid_distribution()
 implicit none
!=====

 if( ALLOCATED(task_grid_proc) )   deallocate(task_grid_proc)
 if( ALLOCATED(ntask_grid_proc) )  deallocate(ntask_grid_proc)
 if( ALLOCATED(task_grid_number) ) deallocate(task_grid_number)

end subroutine destroy_grid_distribution


!=========================================================================
function is_my_grid_task(igrid)
 implicit none
 integer,intent(in) :: igrid
 logical            :: is_my_grid_task
!=====
 
 is_my_grid_task = ( rank_local == task_grid_proc(igrid) )
 
end function is_my_grid_task


!=========================================================================
subroutine distribute_grid_workload()
 implicit none
!=====
 integer            :: igrid,iproc_local
 integer            :: igrid_current
 integer            :: max_grid_per_proc
!=====

#ifndef TODAY
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
#else
 allocate(task_grid_proc(ngrid_mpi))
 allocate(ntask_grid_proc(0:nproc-1))
 allocate(task_grid_number(ngrid_mpi))

 if( parallel_grid) then

   write(stdout,'(/,a)') ' Distributing the grid among procs'
   
   ntask_grid_proc(:) = 0
   max_grid_per_proc = CEILING( DBLE(ngrid_mpi)/DBLE(nproc) )
   write(stdout,*) 'Maximum number of grid points for a single proc',max_grid_per_proc

   iproc_local=0
   do igrid=1,ngrid_mpi

     iproc_local = MODULO(igrid-1,nproc)

     !
     ! A simple check to avoid unexpected surprises
     if( iproc_local < 0 .OR. iproc_local >= nproc ) then
       call die('error in the distribution')
     endif

     task_grid_proc(igrid)        = iproc_local
     ntask_grid_proc(iproc_local) = ntask_grid_proc(iproc_local) + 1 

   enddo

   task_grid_number(:)=0
   igrid_current=0
   do igrid=1,ngrid_mpi
     if( rank == task_grid_proc(igrid) ) then
       igrid_current = igrid_current + 1 
       task_grid_number(igrid) = igrid_current
     endif
   enddo

   if( nproc > 1 ) then
     write(stdout,'(/,a)') ' Distribute work load among procs'
     write(stdout,'(a,x,f8.2)') ' Avg. tasks per cpu:',REAL(ngrid_mpi,dp) / REAL(nproc,dp)
     write(stdout,'(a,i6,a,i10)') ' proc # , grid points',rank,' , ',ntask_grid_proc(rank)
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
#endif


end subroutine distribute_grid_workload


!=========================================================================
function get_ntask()
 implicit none
 integer :: get_ntask
!=====

 get_ntask = ntask_proc(rank)

end function get_ntask


!=========================================================================
subroutine xbcast_ra2d(iproc,array)
 implicit none
 integer,intent(in)     :: iproc
 real(dp),intent(inout) :: array(:,:)
!=====
 integer :: n1,n2
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )

#ifdef HAVE_MPI
 call MPI_BCAST(array,n1*n2,MPI_DOUBLE_PRECISION,iproc,comm_world,ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xbcast_ra2d


!=========================================================================
subroutine xand_l(logical_variable)
 implicit none
 logical,intent(inout) :: logical_variable
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, logical_variable, n1, MPI_LOGICAL, MPI_LAND, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xand_l


!=========================================================================
subroutine xand_la1d(logical_array)
 implicit none
 logical,intent(inout) :: logical_array(:)
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = SIZE(logical_array,DIM=1)

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, logical_array, n1, MPI_LOGICAL, MPI_LAND, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xand_la1d


!=========================================================================
subroutine xand_la2d(logical_array)
 implicit none
 logical,intent(inout) :: logical_array(:,:)
!=====
 integer :: n1,n2
 integer :: ier=0
!=====

 n1 = SIZE(logical_array,DIM=1)
 n2 = SIZE(logical_array,DIM=2)

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, logical_array, n1*n2, MPI_LOGICAL, MPI_LAND, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xand_la2d


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
 call MPI_ALLREDUCE( MPI_IN_PLACE, integer_number, n1, MPI_INTEGER, MPI_MIN, comm_world, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, integer_number, n1, MPI_INTEGER, MPI_MAX, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xmax_i


!=========================================================================
subroutine xmax_r(real_number)
 implicit none
 real(dp),intent(inout) :: real_number
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, real_number, n1, MPI_DOUBLE, MPI_MAX, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xmax_r


!=========================================================================
subroutine xmax_ia2d(array)
 implicit none
 integer,intent(inout) :: array(:,:)
!=====
 integer :: n1,n2
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_INTEGER, MPI_MAX, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xmax_ia2d


!=========================================================================
subroutine xmax_ra1d(array)
 implicit none
 real(dp),intent(inout) :: array(:)
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1, MPI_DOUBLE, MPI_MAX, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xmax_ra1d


!=========================================================================
subroutine xlocal_max_i(integer_number)
 implicit none
 integer,intent(inout) :: integer_number
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, integer_number, n1, MPI_INTEGER, MPI_MAX, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xlocal_max_i


!=========================================================================
subroutine xtrans_max_ia2d(array)
 implicit none
 integer,intent(inout) :: array(:,:)
!=====
 integer :: n1,n2
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_INTEGER, MPI_MAX, comm_trans, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xtrans_max_ia2d


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
 call MPI_ALLREDUCE( MPI_IN_PLACE, real_number, n1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3*n4, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_ra4d


!=========================================================================
subroutine xsum_ca1d(array)
 implicit none
 complex(dpc),intent(inout) :: array(:)
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_ca1d


!=========================================================================
subroutine xsum_ca2d(array)
 implicit none
 complex(dpc),intent(inout) :: array(:,:)
!=====
 integer :: n1,n2
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_ca2d


!=========================================================================
subroutine xsum_ca4d(array)
 implicit none
 complex(dpc),intent(inout) :: array(:,:,:,:)
!=====
 integer :: n1,n2,n3,n4
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )
 n3 = SIZE( array, DIM=3 )
 n4 = SIZE( array, DIM=4 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3*n4, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_ca4d


!=========================================================================
subroutine xsum_procindex_ra2d(iproc,array)
 implicit none
 integer,intent(in)     :: iproc
 real(dp),intent(inout) :: array(:,:)
!=====
 integer :: n1,n2
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )

#ifdef HAVE_MPI
 if( rank == iproc ) then
   call MPI_REDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, iproc, comm_world, ier)
 else
   call MPI_REDUCE( array, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, iproc, comm_world, ier)
 endif
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_procindex_ra2d


!=========================================================================
subroutine xlocal_sum_ra2d(array)
 implicit none
 real(dp),intent(inout) :: array(:,:)
!=====
 integer :: n1,n2
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xlocal_sum_ra2d


!=========================================================================
subroutine xlocal_sum_ra3d(array)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3, MPI_DOUBLE_PRECISION, MPI_SUM, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xlocal_sum_ra3d


!=========================================================================
subroutine xtrans_sum_r(real_number)
 implicit none
 real(dp),intent(inout) :: real_number
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, real_number, n1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_trans, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xtrans_sum_r


!=========================================================================
subroutine xtrans_sum_ra1d(array)
 implicit none
 real(dp),intent(inout) :: array(:)
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_trans, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xtrans_sum_ra1d


!=========================================================================
subroutine xtrans_sum_ra2d(array)
 implicit none
 real(dp),intent(inout) :: array(:,:)
!=====
 integer :: n1,n2
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, comm_trans, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xtrans_sum_ra2d


!=========================================================================
subroutine xtrans_sum_ra3d(array)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3, MPI_DOUBLE_PRECISION, MPI_SUM, comm_trans, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xtrans_sum_ra3d


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

   ! Since the attempted distribution does not fit the total number of CPUs
   ! make a new attempt.
   nprow_sd = MAX( nprow_sd - 1 , 1 )
   npcol_sd = nproc_sca / nprow_sd
   if( nprow_sd * npcol_sd /= nproc_sca ) then
     ! Since the attempted distribution does not fit the total number of CPUs
     ! make a new attempt.
     nprow_sd = MAX( nprow_sd - 1 , 1 )
     npcol_sd = nproc_sca / nprow_sd

     if( nprow_sd * npcol_sd /= nproc_sca ) then
       ! Since the attempted distribution does not fit the total number of CPUs
       ! make a new attempt.
       nprow_sd = MAX( nprow_sd - 1 , 1 )
       npcol_sd = nproc_sca / nprow_sd
       if( nprow_sd * npcol_sd /= nproc_sca ) then
         write(stdout,'(a)') ' Some procs will be idling in the SCALAPACK distribution'
         write(stdout,'(a)') ' This is a waste and it is not yet coded anyway!'
         write(stdout,'(a)') ' => select a number of procs that is not a prime number'
         call die('proc distribution not permitted')
       endif
     endif
   endif
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

! write(stdout,'(/,a)')           ' ==== SCALAPACK info'
! write(stdout,'(a)')             '   Squared distribution'
! write(stdout,'(a50,x,i8)')      'Number of proc:',nprow_sd*npcol_sd
! write(stdout,'(a50,x,i8,x,i8)') 'Grid of procs:',nprow_sd,npcol_sd
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
subroutine init_scalapack_ham(nbf,m_ham,n_ham)
 implicit none

 integer,intent(in)  :: nbf
 integer,intent(out) :: m_ham,n_ham
!=====
 integer :: ier=0
 integer :: color
#ifdef HAVE_SCALAPACK
 integer,external :: NUMROC 
#endif
 integer :: unitfile
!=====

#ifdef HAVE_SCALAPACK
 inquire(file='SCALAPACK_GRID',exist=parallel_ham)

 if( parallel_ham ) then

   open(newunit=unitfile,file='SCALAPACK_GRID',status='old')
   read(unitfile,*) nprow_ham,npcol_ham
   if( nprow_ham * npcol_ham > nproc ) call die('SCALAPACK manual distribution asks for too many processors')
   close(unitfile)

   call BLACS_GET( -1, 0, cntxt_ham )
   call BLACS_GRIDINIT( cntxt_ham, 'R', nprow_ham, npcol_ham )
   call BLACS_GRIDINFO( cntxt_ham, nprow_ham, npcol_ham, iprow_ham, ipcol_ham )

   ! Propagate the scalapack grid to all processors
   call xmax(nprow_ham)
   call xmax(npcol_ham)


   if( cntxt_ham > 0 ) then
     call init_desc('H',nbf,nbf,desc_ham,m_ham,n_ham)
   else
     m_ham = 0
     n_ham = 0
   endif

   allocate(rank_ham_sca_to_mpi(0:nprow_ham-1,0:npcol_ham-1))
   rank_ham_sca_to_mpi(:,:) = -1
   if( iprow_ham >= 0 .AND. ipcol_ham >= 0 ) &
     rank_ham_sca_to_mpi(iprow_ham,ipcol_ham) = rank
   call xmax(rank_ham_sca_to_mpi)

   write(stdout,'(/,a)')           ' ==== SCALAPACK Hamiltonian'
   write(stdout,'(a50,x,i8)')      'Number of dedicated processors:',nprow_ham * npcol_ham
   write(stdout,'(a50,x,i8,x,i8)')   'Grid of dedicated processors:',nprow_ham,npcol_ham

   ! Distribute the remaing procs for auxiliary basis and grid points
   color = MODULO( rank , nprow_ham * npcol_ham )

   call MPI_COMM_SPLIT(comm_world,color,rank,comm_local,ier);
   call MPI_COMM_SIZE(comm_local,nproc_local,ier)
   call MPI_COMM_RANK(comm_local,rank_local,ier)

   write(stdout,'(a50,x,i8)')      'Number of local processors:',nproc_local

   call xlocal_max(m_ham)
   call xlocal_max(n_ham)
   call xlocal_max(iprow_ham)
   call xlocal_max(ipcol_ham)


   ! Define the transversal communicator
   color = rank / ( nprow_ham * npcol_ham )

   call MPI_COMM_SPLIT(comm_world,color,rank,comm_trans,ier);
   call MPI_COMM_SIZE(comm_trans,nproc_trans,ier)
   call MPI_COMM_RANK(comm_trans,rank_trans,ier)

   allocate(rank_sca_to_mpi(0:nprow_ham-1,0:npcol_ham-1))
   rank_sca_to_mpi(:,:) = -1
   rank_sca_to_mpi(iprow_ham,ipcol_ham) = rank_trans
   call xtrans_max(rank_sca_to_mpi)


 else

   if( rank == 0 ) then
     cntxt_ham = 1
   else
     cntxt_ham = -1
   endif
   nprow_ham = 1
   npcol_ham = 1
   iprow_ham = 0
   ipcol_ham = 0
   m_ham = nbf
   n_ham = nbf

   comm_local  = comm_world
   nproc_local = nproc
   rank_local  = rank

   comm_trans  = MPI_COMM_SELF
   nproc_trans = 1
   rank_trans  = 0

 endif

#else

 ! Fake values
 cntxt_ham = 1
 nprow_ham = 1
 npcol_ham = 1
 iprow_ham = 0
 ipcol_ham = 0
 m_ham = nbf
 n_ham = nbf

 nproc_local = 1
 rank_local  = 0

 nproc_trans = 1
 rank_trans  = 0

#endif


end subroutine init_scalapack_ham


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
 case('H')
   mlocal = NUMROC(mglobal,block_row,iprow_ham,first_row,nprow_ham)
   nlocal = NUMROC(nglobal,block_col,ipcol_ham,first_col,npcol_ham)
   call DESCINIT(desc,mglobal,nglobal,block_row,block_col,first_row,first_col,cntxt_ham,MAX(1,mlocal),info)
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
   call die('BUG')
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
 integer              :: desc_eigvec(ndel)
 integer              :: lwork,info
 real(dp),allocatable :: work(:)
 real(dp)             :: eigvec(mlocal,nlocal)
 integer              :: neigval,neigvec
 integer,allocatable  :: iwork(:)
 integer              :: liwork
#ifndef SELECT_PDSYEVR
 real(dp)             :: ABSTOL
 integer              :: iclustr(2*nprow_sd*npcol_sd)
 real(dp)             :: gap(nprow_sd*npcol_sd)
 integer              :: ifail(nglobal)
 real(dp),external    :: PDLAMCH
#endif
!=====

#ifdef HAVE_SCALAPACK
 ! fake descriptor ! Why do I need this?
 desc_eigvec(:) = desc(:)

 

 !
 ! First call to get the dimension of the array work
 lwork = -1
 liwork = -1
 allocate(work(1))
 allocate(iwork(1))
#ifndef SELECT_PDSYEVR
 ABSTOL = PDLAMCH(desc(2), 'U')
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
#else
 call PDSYEVR('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              neigval,neigvec,eigval,                                &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              info)
#endif


 !
 ! Second call to actually perform the diago
 lwork = NINT(work(1))
 liwork = iwork(1)

 deallocate(work)
 deallocate(iwork)
 allocate(work(lwork))
 allocate(iwork(liwork))
#ifndef SELECT_PDSYEVR
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
#else
 call PDSYEVR('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              neigval,neigvec,eigval,                                &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              info)
#endif
 deallocate(work)
 deallocate(iwork)


 matrix(:,:) = eigvec(:,:)

#else
 eigval(:) = 0.0_dp
 eigvec(:,:) = 0.0_dp
#endif


end subroutine diagonalize_sca


!=========================================================================
subroutine diagonalize_sca_outofplace(desc,nglobal,mlocal,nlocal,matrix,eigval, &
             desc_eigvec,m_eigvec,n_eigvec,eigvec)
 implicit none
 integer,intent(in)     :: desc(ndel),nglobal,mlocal,nlocal
 real(dp),intent(inout) :: matrix(mlocal,nlocal)
 real(dp),intent(out)   :: eigval(nglobal)
 integer,intent(in)     :: desc_eigvec(ndel),m_eigvec,n_eigvec
 real(dp),intent(out)   :: eigvec(m_eigvec,n_eigvec)
!=====
 integer              :: lwork,info
 real(dp),allocatable :: work(:)
 integer              :: neigval,neigvec
 integer,allocatable  :: iwork(:)
 integer              :: liwork
#ifndef SELECT_PDSYEVR
 real(dp)             :: ABSTOL
 integer              :: iclustr(2*nprow_sd*npcol_sd)
 real(dp)             :: gap(nprow_sd*npcol_sd)
 integer              :: ifail(nglobal)
 real(dp),external    :: PDLAMCH
#endif
!=====

#ifdef HAVE_SCALAPACK


 !
 ! First call to get the dimension of the array work
 lwork = -1
 liwork = -1
 allocate(work(1))
 allocate(iwork(1))
#ifndef SELECT_PDSYEVR
 ABSTOL = PDLAMCH(desc(2), 'U')
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
#else
 call PDSYEVR('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              neigval,neigvec,eigval,                                &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              info)
#endif


 !
 ! Second call to actually perform the diago
 lwork = NINT(work(1))
 liwork = iwork(1)

 deallocate(work)
 deallocate(iwork)
 allocate(work(lwork))
 allocate(iwork(liwork))
#ifndef SELECT_PDSYEVR
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
#else
 call PDSYEVR('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              neigval,neigvec,eigval,                                &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              info)
#endif
 deallocate(work)
 deallocate(iwork)


#else
 eigval(:) = 0.0_dp
 eigvec(:,:) = 0.0_dp
#endif


end subroutine diagonalize_sca_outofplace


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
 case('H')
   if( iprow_ham == INDXG2P(iglobal,block_row,0,first_row,nprow_ham) ) then
     rowindex_global_to_local = INDXG2L(iglobal,block_row,0,first_row,nprow_ham)
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
   call die('BUG')
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
 case('H')
   if( ipcol_ham == INDXG2P(iglobal,block_col,0,first_col,npcol_ham) ) then
     colindex_global_to_local = INDXG2L(iglobal,block_col,0,first_col,npcol_ham)
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
   call die('BUG')
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
 case('H')
   rowindex_local_to_global_distrib = INDXL2G(ilocal,block_row,iprow_ham,first_row,nprow_ham)
 case('R')
   rowindex_local_to_global_distrib = INDXL2G(ilocal,block_row,iprow_rd,first_row,nprow_rd)
 case('C')
   rowindex_local_to_global_distrib = INDXL2G(ilocal,block_row,iprow_cd,first_row,nprow_cd)
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   call die('BUG')
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
 case('H')
   colindex_local_to_global_distrib = INDXL2G(ilocal,block_col,ipcol_ham,first_col,npcol_ham)
 case('R')
   colindex_local_to_global_distrib = INDXL2G(ilocal,block_col,ipcol_rd,first_col,npcol_rd)
 case('C')
   colindex_local_to_global_distrib = INDXL2G(ilocal,block_col,ipcol_cd,first_col,npcol_cd)
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   call die('BUG')
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


#ifdef HAVE_SCALAPACK
!=========================================================================
subroutine diagonalize_scalapack(nmat,matrix,eigval)
 implicit none
 integer,intent(in)     :: nmat
 real(dp),intent(inout) :: matrix(nmat,nmat)
 real(dp),intent(out)   :: eigval(nmat)
!=====
 integer :: cntxt
 integer :: mlocal,nlocal
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
 integer :: imat,jmat,imat_local,jmat_local
 integer :: descm(ndel),descz(ndel)
 real(dp) :: alpha
 real(dp),allocatable :: matrix_local(:,:)
 real(dp),allocatable :: work(:)
 integer :: lwork
 integer :: rank_sca,nprocs_sca
 integer,external :: NUMROC,INDXG2L,INDXG2P
!=====

 nprow = MIN(nprow_sd,nmat/SCALAPACK_MIN)
 npcol = MIN(npcol_sd,nmat/SCALAPACK_MIN)
 nprow = MAX(nprow,1)
 npcol = MAX(npcol,1)

 call BLACS_PINFO( rank_sca, nprocs_sca )

 call BLACS_GET( -1, 0, cntxt )
 call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
 write(stdout,'(a,i4,a,i4)') ' Diagonalization using SCALAPACK with a grid',nprow,' x ',npcol

 !
 ! Participate to the diagonalization only if the CPU has been selected 
 ! in the grid
 if( cntxt >= 0  ) then

   mlocal = NUMROC(nmat,block_row,iprow,first_row,nprow)
   nlocal = NUMROC(nmat,block_col,ipcol,first_col,npcol)
  
   allocate(matrix_local(mlocal,nlocal))
    
   call DESCINIT(descm,nmat,nmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mlocal),info)

  
   ! set up the local copy of the matrix
   do jmat=1,nmat
     if( INDXG2P(jmat,block_col,0,first_col,npcol) /= ipcol ) cycle
     do imat=1,nmat
       if( INDXG2P(imat,block_row,0,first_row,nprow) /= iprow ) cycle
       imat_local = INDXG2L(imat,block_row,0,first_row,nprow)
       jmat_local = INDXG2L(jmat,block_col,0,first_col,npcol)

       matrix_local(imat_local,jmat_local) = matrix(imat,jmat)

     enddo
   enddo
  

   call diagonalize_sca(descm,nmat,mlocal,nlocal,matrix_local,eigval)
  
   ! Nullify the eigval array for all CPUs but one, so that the all_reduce
   ! operation in the end yields the correct value
   ! Of course, using a broadcast instead would be a better solution, but I'm so lazy
   if(rank_sca /= 0 ) eigval(:) = 0.0_dp
  
   matrix(:,:) = 0.0_dp
   do jmat=1,nmat
     if( INDXG2P(jmat,block_col,0,first_col,npcol) /= ipcol ) cycle
     do imat=1,nmat
       if( INDXG2P(imat,block_row,0,first_row,nprow) /= iprow ) cycle
       imat_local = INDXG2L(imat,block_row,0,first_row,nprow)
       jmat_local = INDXG2L(jmat,block_col,0,first_col,npcol)
  
       matrix(imat,jmat) = matrix_local(imat_local,jmat_local)
  
     enddo
   enddo

   deallocate(matrix_local)
   call BLACS_GRIDEXIT( cntxt )

 else
   ! Those CPUs that are not used in the SCALAPACK process
   matrix(:,:) = 0.0_dp
   eigval(:) = 0.0_dp
 endif

 call xsum(matrix)
 call xsum(eigval)


end subroutine diagonalize_scalapack
#endif


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
