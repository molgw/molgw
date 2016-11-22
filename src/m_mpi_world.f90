!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This modules contains
! the MPI information for parallelization over all procs (=world)
!
!=========================================================================
module m_mpi_world
 use m_definitions
 use m_warning,only: die
#ifdef HAVE_MPI
 use mpi
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

 !
 ! "world" communicator
 ! 
 integer,public    :: comm_world             ! communicator over all the procs (=world)
 integer,public    :: nproc_world  = 1       ! number of procs in the world communicator
 integer,public    :: rank_world   = 0       ! index           in the world communicator


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
 
 interface xsum_world
   module procedure xsum_world_r
   module procedure xsum_world_ra1d
   module procedure xsum_world_ra2d
   module procedure xsum_world_ra3d
   module procedure xsum_world_ra4d
   module procedure xsum_world_ca1d
   module procedure xsum_world_ca2d
   module procedure xsum_world_ca3d
   module procedure xsum_world_ca4d
 end interface

 interface xbcast_world
   module procedure xbcast_world_ra1d
   module procedure xbcast_world_ra2d
   module procedure xbcast_world_ra3d
 end interface

 interface xand_world
   module procedure xand_world_l
   module procedure xand_world_la1d
   module procedure xand_world_la2d
 end interface


contains


!=========================================================================
subroutine barrier_world()
 implicit none
 integer :: ier
!=====

#ifdef HAVE_MPI
 call MPI_BARRIER(MPI_COMM_WORLD,ier)
#endif

end subroutine barrier_world


!=========================================================================
subroutine xbcast_world_ra1d(iproc,array)
 implicit none
 integer,intent(in)     :: iproc
 real(dp),intent(inout) :: array(:)
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )

#ifdef HAVE_MPI
 call MPI_BCAST(array,n1,MPI_DOUBLE_PRECISION,iproc,comm_world,ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xbcast_world_ra1d


!=========================================================================
subroutine xbcast_world_ra2d(iproc,array)
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

end subroutine xbcast_world_ra2d


!=========================================================================
subroutine xbcast_world_ra3d(iproc,array)
 implicit none
 integer,intent(in)     :: iproc
 real(dp),intent(inout) :: array(:,:,:)
!=====
 integer :: n1,n2,n3
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )
 n3 = SIZE( array, DIM=3 )

#ifdef HAVE_MPI
 call MPI_BCAST(array,n1*n2*n3,MPI_DOUBLE_PRECISION,iproc,comm_world,ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xbcast_world_ra3d


!=========================================================================
subroutine xand_world_l(logical_variable)
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

end subroutine xand_world_l


!=========================================================================
subroutine xand_world_la1d(logical_array)
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

end subroutine xand_world_la1d


!=========================================================================
subroutine xand_world_la2d(logical_array)
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

end subroutine xand_world_la2d


!=========================================================================
subroutine xmin_world_r(real_number)
 implicit none
 real(dp),intent(inout) :: real_number
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, real_number, n1, MPI_DOUBLE, MPI_MIN, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xmin_world_r


!=========================================================================
subroutine xmax_world_i(integer_number)
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

end subroutine xmax_world_i


!=========================================================================
subroutine xmax_world_r(real_number)
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

end subroutine xmax_world_r


!=========================================================================
subroutine xmax_world_ia2d(array)
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

end subroutine xmax_world_ia2d


!=========================================================================
subroutine xmax_world_ra1d(array)
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

end subroutine xmax_world_ra1d


!=========================================================================
subroutine xsum_world_r(real_number)
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

end subroutine xsum_world_r


!=========================================================================
subroutine xsum_world_ra1d(array)
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

end subroutine xsum_world_ra1d


!=========================================================================
subroutine xsum_world_ra2d(array)
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

end subroutine xsum_world_ra2d


!=========================================================================
subroutine xsum_world_ra3d(array)
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

end subroutine xsum_world_ra3d


!=========================================================================
subroutine xsum_world_ra4d(array)
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

end subroutine xsum_world_ra4d


!=========================================================================
subroutine xsum_world_ca1d(array)
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

end subroutine xsum_world_ca1d


!=========================================================================
subroutine xsum_world_ca2d(array)
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

end subroutine xsum_world_ca2d


!=========================================================================
subroutine xsum_world_ca3d(array)
 implicit none
 complex(dpc),intent(inout) :: array(:,:,:)
!=====
 integer :: n1,n2,n3
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )
 n3 = SIZE( array, DIM=3 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_world, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_world_ca3d


!=========================================================================
subroutine xsum_world_ca4d(array)
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

end subroutine xsum_world_ca4d


!=========================================================================
end module m_mpi_world
!=========================================================================
