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
 use m_mpi_auxil_grid
 use m_mpi_ortho
#ifdef HAVE_MPI
 use mpi
#endif


 logical,parameter :: parallel_grid      = .TRUE.
 logical,parameter :: parallel_auxil     = .TRUE.

 logical,protected :: parallel_ham       = .FALSE.
 logical,protected :: parallel_buffer    = .TRUE.

#ifdef HAVE_SCALAPACK
 logical,parameter :: parallel_scalapack = .TRUE.
#else
 logical,parameter :: parallel_scalapack = .FALSE.
#endif

!===================================================
! MPI distribution 
!  Example: nproc_ortho = 2 x  nproc_auxil_grid = 8  = nproc_world = 16
!
! comm_world
!                                        
! rank_auxil_grid    0 |  1 |  2 |     |  7 
! rank_ortho       ---------------------------
!      0             0 |  2 |  4 | ... | 14 |-> comm_auxil_grid
!      1             1 |  3 |  5 | ... | 15 |-> comm_auxil_grid
!                  ---------------------------
!                    | |    |    | ... |  | |
!                    v                    v
!                 comm_ortho           comm_ortho
!===================================================


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


 !
 ! Interfaces for high-level MPI reduce operations
 ! auxil or grid series
 !
 interface xlocal_max
   module procedure xlocal_max_i
 end interface

 interface xtrans_max
   module procedure xtrans_max_ia2d
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
   module procedure colindex_local_to_global_descriptor
 end interface

 interface rowindex_local_to_global
   module procedure rowindex_local_to_global_distrib
   module procedure rowindex_local_to_global_procindex
   module procedure rowindex_local_to_global_descriptor
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
subroutine init_mpi_world()
 implicit none

!=====
 integer :: color
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
subroutine init_mpi_details(nproc_ortho_in)
 implicit none

 integer,intent(in) :: nproc_ortho_in
!=====
 integer :: color
 integer :: ier
!=====

#ifdef HAVE_MPI
 !
 ! Set up auxil or grid communicator
 !
 nproc_ortho = nproc_ortho_in
 
 nproc_auxil_grid = nproc_world / nproc_ortho

 color = MODULO( rank_world , nproc_ortho )
 call MPI_COMM_SPLIT(comm_world,color,rank_world,comm_auxil_grid,ier);

 call MPI_COMM_SIZE(comm_auxil_grid,nproc_auxil_grid,ier)
 call MPI_COMM_RANK(comm_auxil_grid,rank_auxil_grid,ier)
 if( nproc_auxil_grid /= nproc_world / nproc_ortho ) then
   write(stdout,*) rank_world,color,nproc_auxil_grid,nproc_world,nproc_ortho
   call die('Problem in init_mpi')
 endif

 !
 ! Set up ortho communicator
 !
 nproc_ortho = nproc_world / nproc_auxil_grid

 color = rank_world / nproc_ortho
 call MPI_COMM_SPLIT(comm_world,color,rank_world,comm_ortho,ier);
 call MPI_COMM_RANK(comm_ortho,rank_ortho,ier)


#else
 nproc_ortho = 1
 rank_ortho  = 0
 nproc_auxil_grid = 1
 rank_auxil_grid  = 0
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

#ifdef HAVE_MPI
  write(stdout,'(/,a)')       ' ==== MPI info'
  write(stdout,'(a50,x,i6)')  'Number of proc:',nproc_world
  write(stdout,'(a50,x,i6)')  'nproc_auxil or grid',nproc_auxil_grid
  write(stdout,'(a50,x,i6)')  'nproc_ortho        ',nproc_ortho
  write(stdout,'(a50,x,i6)')  'Master proc is:',iomaster
  write(stdout,'(a50,6x,l1)') 'Parallelize auxiliary basis:',parallel_auxil
  write(stdout,'(a50,6x,l1)')  'Parallelize XC grid points:',parallel_grid
  write(stdout,'(a50,6x,l1)')               'Use SCALAPACK:',parallel_scalapack
  write(stdout,'(a50,6x,l1)')    'Parallel using a buffer :',parallel_buffer
  write(stdout,'(/)')
#endif

end subroutine init_mpi_details


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
subroutine distribute_auxil_basis(nbf_auxil_basis,nbf_auxil_basis_local)
 implicit none

 integer,intent(in)  :: nbf_auxil_basis
 integer,intent(out) :: nbf_auxil_basis_local
!=====
 integer :: ibf
 integer :: ibf_local
 integer :: iproc
!=====

 if( parallel_buffer ) then
   ! Use nproc instead nproc_local
  
   allocate(iproc_ibf_auxil(nbf_auxil_basis))
   allocate(nbf_local_iproc(0:nproc_auxil_grid-1))
  
   iproc              = nproc_auxil_grid - 1
   nbf_local_iproc(:) = 0
   do ibf=1,nbf_auxil_basis
  
     iproc = MODULO(iproc+1,nproc_auxil_grid)
  
     iproc_ibf_auxil(ibf) = iproc
  
     nbf_local_iproc(iproc) = nbf_local_iproc(iproc) + 1
  
   enddo
  
   nbf_auxil_basis_local = nbf_local_iproc(rank_auxil_grid)
  
   allocate(ibf_auxil_g(nbf_auxil_basis_local))
   allocate(ibf_auxil_l(nbf_auxil_basis))
   ibf_local = 0
   do ibf=1,nbf_auxil_basis
     if( rank_auxil_grid == iproc_ibf_auxil(ibf) ) then
       ibf_local = ibf_local + 1
       ibf_auxil_g(ibf_local) = ibf
       ibf_auxil_l(ibf)       = ibf_local
     endif
   enddo
  
 else

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
  
 endif

 write(stdout,'(/,a)') ' Distribute auxiliary basis functions among processors'
 do iproc=0,0
   write(stdout,'(a,i4,a,i6,a)')   ' Proc: ',iproc,' treats ',nbf_local_iproc(iproc),' auxiliary basis functions'
 enddo



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
 allocate(nbf_local_iproc_lr(0:nproc_auxil_grid-1))

 iproc = nproc_auxil_grid-1
 nbf_local_iproc_lr(:) = 0
 do ibf=1,nbf_auxil_basis

   iproc = MODULO(iproc+1,nproc_auxil_grid)

   iproc_ibf_auxil_lr(ibf) = iproc

   nbf_local_iproc_lr(iproc) = nbf_local_iproc_lr(iproc) + 1

 enddo

 nbf_auxil_basis_local = nbf_local_iproc_lr(rank_auxil_grid)

 allocate(ibf_auxil_g_lr(nbf_auxil_basis_local))
 allocate(ibf_auxil_l_lr(nbf_auxil_basis))
 ibf_local = 0
 do ibf=1,nbf_auxil_basis
   if( rank_auxil_grid == iproc_ibf_auxil_lr(ibf) ) then
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
subroutine init_dft_grid_distribution(ngrid)
 implicit none
 integer,intent(inout) :: ngrid
!=====

 ngrid_mpi = ngrid

 if( parallel_buffer ) then
   if( nproc_auxil_grid > 1 .AND. parallel_grid ) then
     write(stdout,'(/,a)') ' Initializing the distribution of the quadrature grid points'
   endif
 else
   if( nproc_local > 1 .AND. parallel_grid ) then
     write(stdout,'(/,a)') ' Initializing the distribution of the quadrature grid points'
   endif
 endif

 call distribute_grid_workload()

 if( parallel_buffer ) then
   ngrid = ntask_grid_proc(rank_auxil_grid)
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
   is_my_grid_task = ( rank_auxil_grid == task_grid_proc(igrid) )
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
   allocate(ntask_grid_proc(0:nproc_auxil_grid-1))
   allocate(task_grid_number(ngrid_mpi))
  
   if( parallel_grid) then
  
     write(stdout,'(/,a)') ' Distributing the grid among procs'
     
     ntask_grid_proc(:) = 0
     max_grid_per_proc = CEILING( DBLE(ngrid_mpi)/DBLE(nproc_auxil_grid) )
     write(stdout,*) 'Maximum number of grid points for a single proc',max_grid_per_proc
  
     iproc_local=0
     do igrid=1,ngrid_mpi
  
       iproc_local = MODULO(igrid-1,nproc_auxil_grid)
  
       !
       ! A simple check to avoid unexpected surprises
       if( iproc_local < 0 .OR. iproc_local >= nproc_auxil_grid ) then
         call die('error in the distribution')
       endif
  
       task_grid_proc(igrid)        = iproc_local
       ntask_grid_proc(iproc_local) = ntask_grid_proc(iproc_local) + 1 
  
     enddo
  
     task_grid_number(:)=0
     igrid_current=0
     do igrid=1,ngrid_mpi
       if( rank_auxil_grid == task_grid_proc(igrid) ) then
         igrid_current = igrid_current + 1 
         task_grid_number(igrid) = igrid_current
       endif
     enddo
  
     if( nproc_auxil_grid > 1 ) then
       write(stdout,'(/,a)') ' Distribute work load among procs'
       write(stdout,'(a,x,f8.2)') ' Avg. tasks per cpu:',REAL(ngrid_mpi,dp) / REAL(nproc_auxil_grid,dp)
       write(stdout,'(a,i6,a,i10)') ' proc # , grid points',rank_auxil_grid,' , ',ntask_grid_proc(rank_auxil_grid)
     endif
  
   else
     !
     ! if parallel_grid is false,
     ! faking the code with trivial values
     ntask_grid_proc(:) = ngrid_mpi
     task_grid_proc(:)  = rank_auxil_grid
     do igrid=1,ngrid_mpi
       task_grid_number(igrid) = igrid
     enddo
  
   endif
 endif


end subroutine distribute_grid_workload


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
subroutine init_scalapack_ham(nbf,scalapack_nprow,scalapack_npcol,m_ham,n_ham)
 implicit none

 integer,intent(in)  :: nbf
 integer,intent(in)  :: scalapack_nprow,scalapack_npcol
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
 parallel_ham = scalapack_nprow * scalapack_npcol > 1

 if( parallel_ham ) then

   nprow_ham = scalapack_nprow
   npcol_ham = scalapack_npcol
   if( nprow_ham * npcol_ham > nproc_auxil_grid ) call die('SCALAPACK manual distribution asks for too many processors')

   call BLACS_GET( -1, 0, cntxt_ham )
   call BLACS_GRIDINIT( cntxt_ham, 'R', nprow_ham, npcol_ham )
   call BLACS_GRIDINFO( cntxt_ham, nprow_ham, npcol_ham, iprow_ham, ipcol_ham )

   ! Propagate the scalapack grid to all processors
   call xmax_world(nprow_ham)
   call xmax_world(npcol_ham)


   if( cntxt_ham > 0 ) then
     call init_desc('H',nbf,nbf,desc_ham,m_ham,n_ham)
   else
     m_ham = 0
     n_ham = 0
   endif

   allocate(rank_ham_sca_to_mpi(0:nprow_ham-1,0:npcol_ham-1))
   rank_ham_sca_to_mpi(:,:) = -1
   if( iprow_ham >= 0 .AND. ipcol_ham >= 0 ) &
     rank_ham_sca_to_mpi(iprow_ham,ipcol_ham) = rank_auxil_grid
   call xmax_world(rank_ham_sca_to_mpi)

   write(stdout,'(/,a)')           ' ==== SCALAPACK Hamiltonian'
   write(stdout,'(a50,x,i8)')      'Number of dedicated processors:',nprow_ham * npcol_ham
   write(stdout,'(a50,x,i8,x,i8)')   'Grid of dedicated processors:',nprow_ham,npcol_ham

   ! Distribute the remaing procs for auxiliary basis and grid points
   color = MODULO( rank_auxil_grid , nprow_ham * npcol_ham )
   call MPI_COMM_SPLIT(comm_auxil_grid,color,rank_auxil_grid,comm_local,ier);
   call MPI_COMM_SIZE(comm_local,nproc_local,ier)
   call MPI_COMM_RANK(comm_local,rank_local,ier)

   write(stdout,'(a50,x,i8)')      'Number of local processors:',nproc_local

   call xlocal_max(m_ham)
   call xlocal_max(n_ham)
   call xlocal_max(iprow_ham)
   call xlocal_max(ipcol_ham)


   ! Define the transversal communicator
   color = rank_auxil_grid / ( nprow_ham * npcol_ham )

   call MPI_COMM_SPLIT(comm_auxil_grid,color,rank_auxil_grid,comm_trans,ier);
   call MPI_COMM_SIZE(comm_trans,nproc_trans,ier)
   call MPI_COMM_RANK(comm_trans,rank_trans,ier)

   allocate(rank_sca_to_mpi(0:nprow_ham-1,0:npcol_ham-1))
   rank_sca_to_mpi(:,:) = -1
   rank_sca_to_mpi(iprow_ham,ipcol_ham) = rank_trans
   call xtrans_max(rank_sca_to_mpi)


 else

   !
   ! Hamiltonian and C matrix are NOT distributed with SCALAPACK
   !
   cntxt_ham = 1
   nprow_ham = 1
   npcol_ham = 1
   iprow_ham = 0
   ipcol_ham = 0
   m_ham = nbf
   n_ham = nbf

   comm_local  = comm_world
   nproc_local = nproc_world
   rank_local  = rank_world

   comm_trans  = MPI_COMM_SELF
   nproc_trans = 1
   rank_trans  = 0

 endif

#else

 ! Fake values
 if( rank_world == 0 ) then
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
function rowindex_local_to_global_descriptor(desc,ilocal)
 implicit none
 integer,intent(in)          :: desc(ndel),ilocal
 integer                     :: rowindex_local_to_global_descriptor
!=====
#ifdef HAVE_SCALAPACK
 integer          :: iprow,ipcol,nprow,npcol
 integer,external :: INDXL2G
#endif
!=====

#ifdef HAVE_SCALAPACK
 call BLACS_GRIDINFO(desc(2),nprow,npcol,iprow,ipcol)
 rowindex_local_to_global_descriptor = INDXL2G(ilocal,block_row,iprow,first_row,nprow)
#else
 rowindex_local_to_global_descriptor = ilocal
#endif

end function rowindex_local_to_global_descriptor


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


!=========================================================================
function colindex_local_to_global_descriptor(desc,ilocal)
 implicit none
 integer,intent(in)          :: desc(ndel),ilocal
 integer                     :: colindex_local_to_global_descriptor
!=====
#ifdef HAVE_SCALAPACK
 integer          :: iprow,ipcol,nprow,npcol
 integer,external :: INDXL2G
#endif
!=====

#ifdef HAVE_SCALAPACK
 call BLACS_GRIDINFO(desc(2),nprow,npcol,iprow,ipcol)
 colindex_local_to_global_descriptor = INDXL2G(ilocal,block_col,ipcol,first_col,npcol)
#else
 colindex_local_to_global_descriptor = ilocal
#endif

end function colindex_local_to_global_descriptor


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
