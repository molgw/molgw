!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This modules contains
! the MPI information for auxil/grid parallelization
!
!=========================================================================
module m_mpi_auxil_grid
 use m_definitions
 use m_warning,only: die
#ifdef HAVE_MPI
 use mpi
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

 !
 ! "auxil/grid" communicator
 ! 
 integer,public    :: comm_auxil_grid        ! communicator over auxiliary basis functions OR DFT grid points
 integer,public    :: nproc_auxil_grid  = 1  ! number of procs in the auxil_grid communicator
 integer,public    :: rank_auxil_grid   = 0  ! index           in the auxil_grid communicator


 !
 ! Interfaces for high-level MPI reduce operations
 ! "auxil/grid" series
 !
 interface xbcast
   module procedure xbcast_ra1d
   module procedure xbcast_ra2d
 end interface

 interface xand
   module procedure xand_l
   module procedure xand_la1d
   module procedure xand_la2d
 end interface

 interface xmin_auxil
   module procedure xmin_i
 end interface

 interface xmin_grid
   module procedure xmin_i
 end interface

 interface xmax_auxil
   module procedure xmax_i
   module procedure xmax_r
   module procedure xmax_ia2d
   module procedure xmax_ra1d
 end interface

 interface xmax_grid
   module procedure xmax_i
   module procedure xmax_r
   module procedure xmax_ia2d
   module procedure xmax_ra1d
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

 interface xsum_auxil
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

 interface xsum_grid
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



contains


!=========================================================================
subroutine xbcast_ra1d(iproc,array)
 implicit none
 integer,intent(in)     :: iproc
 real(dp),intent(inout) :: array(:)
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )

#ifdef HAVE_MPI
 call MPI_BCAST(array,n1,MPI_DOUBLE_PRECISION,iproc,comm_auxil_grid,ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xbcast_ra1d


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
 call MPI_BCAST(array,n1*n2,MPI_DOUBLE_PRECISION,iproc,comm_auxil_grid,ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, logical_variable, n1, MPI_LOGICAL, MPI_LAND, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, logical_array, n1, MPI_LOGICAL, MPI_LAND, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, logical_array, n1*n2, MPI_LOGICAL, MPI_LAND, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, integer_number, n1, MPI_INTEGER, MPI_MIN, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, integer_number, n1, MPI_INTEGER, MPI_MAX, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, real_number, n1, MPI_DOUBLE, MPI_MAX, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_INTEGER, MPI_MAX, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1, MPI_DOUBLE, MPI_MAX, comm_auxil_grid, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xmax_ra1d


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
 call MPI_ALLREDUCE( MPI_IN_PLACE, real_number, n1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3, MPI_DOUBLE_PRECISION, MPI_SUM, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3*n4, MPI_DOUBLE_PRECISION, MPI_SUM, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_auxil_grid, ier)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3*n4, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_auxil_grid, ier)
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
 if( rank_auxil_grid == iproc ) then
   call MPI_REDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, iproc, comm_auxil_grid, ier)
 else
   call MPI_REDUCE( array, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, iproc, comm_auxil_grid, ier)
 endif
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_reduce'
 endif

end subroutine xsum_procindex_ra2d


!=========================================================================
end module m_mpi_auxil_grid
!=========================================================================
