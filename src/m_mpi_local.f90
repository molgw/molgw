!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This modules contains
! the MPI information for local parallelization
!
!=========================================================================
module m_mpi_local
 use m_definitions
 use m_warning,only: die
#ifdef HAVE_MPI
 use mpi
#endif


 !
 ! "local" communicator
 !
 integer,public    :: comm_local        ! communicator over local processors
 integer,public    :: nproc_local  = 1  ! number of procs in the local communicator
 integer,public    :: rank_local   = 0  ! index           in the local communicator


 !
 ! Interfaces for high-level MPI reduce operations
 ! "local" series
 !
 interface xbcast_local
   module procedure xbcast_local_ra1d
   module procedure xbcast_local_ra2d
 end interface

 interface xand_local
   module procedure xand_local_l
   module procedure xand_local_la1d
   module procedure xand_local_la2d
 end interface

 interface xmin_local
   module procedure xmin_local_i
 end interface

 interface xmax_local
   module procedure xmax_local_i
   module procedure xmax_local_r
   module procedure xmax_local_ia2d
   module procedure xmax_local_ra1d
 end interface


 interface xsum_local
   module procedure xsum_local_r
   module procedure xsum_local_ra1d
   module procedure xsum_local_ra2d
   module procedure xsum_local_ra3d
   module procedure xsum_local_ra4d
   module procedure xsum_local_ca1d
   module procedure xsum_local_ca2d
   module procedure xsum_local_ca4d
   module procedure xsum_local_procindex_ra2d
 end interface



contains


!=========================================================================
subroutine xbcast_local_ra1d(iproc,array)
 implicit none
 integer,intent(in)     :: iproc
 real(dp),intent(inout) :: array(:)
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )

#ifdef HAVE_MPI
 call MPI_BCAST(array,n1,MPI_DOUBLE_PRECISION,iproc,comm_local,ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xbcast_local_ra1d


!=========================================================================
subroutine xbcast_local_ra2d(iproc,array)
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
 call MPI_BCAST(array,n1*n2,MPI_DOUBLE_PRECISION,iproc,comm_local,ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xbcast_local_ra2d


!=========================================================================
subroutine xand_local_l(logical_variable)
 implicit none
 logical,intent(inout) :: logical_variable
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, logical_variable, n1, MPI_LOGICAL, MPI_LAND, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xand_local_l


!=========================================================================
subroutine xand_local_la1d(logical_array)
 implicit none
 logical,intent(inout) :: logical_array(:)
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = SIZE(logical_array,DIM=1)

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, logical_array, n1, MPI_LOGICAL, MPI_LAND, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xand_local_la1d


!=========================================================================
subroutine xand_local_la2d(logical_array)
 implicit none
 logical,intent(inout) :: logical_array(:,:)
!=====
 integer :: n1,n2
 integer :: ier=0
!=====

 n1 = SIZE(logical_array,DIM=1)
 n2 = SIZE(logical_array,DIM=2)

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, logical_array, n1*n2, MPI_LOGICAL, MPI_LAND, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xand_local_la2d


!=========================================================================
subroutine xmin_local_i(integer_number)
 implicit none
 integer,intent(inout) :: integer_number
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, integer_number, n1, MPI_INTEGER, MPI_MIN, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xmin_local_i


!=========================================================================
subroutine xmax_local_i(integer_number)
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

end subroutine xmax_local_i


!=========================================================================
subroutine xmax_local_r(real_number)
 implicit none
 real(dp),intent(inout) :: real_number
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, real_number, n1, MPI_DOUBLE, MPI_MAX, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xmax_local_r


!=========================================================================
subroutine xmax_local_ia2d(array)
 implicit none
 integer,intent(inout) :: array(:,:)
!=====
 integer :: n1,n2
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_INTEGER, MPI_MAX, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xmax_local_ia2d


!=========================================================================
subroutine xmax_local_ra1d(array)
 implicit none
 real(dp),intent(inout) :: array(:)
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1, MPI_DOUBLE, MPI_MAX, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xmax_local_ra1d


!=========================================================================
subroutine xsum_local_r(real_number)
 implicit none
 real(dp),intent(inout) :: real_number
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = 1

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, real_number, n1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_local_r


!=========================================================================
subroutine xsum_local_ra1d(array)
 implicit none
 real(dp),intent(inout) :: array(:)
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_local_ra1d


!=========================================================================
subroutine xsum_local_ra2d(array)
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

end subroutine xsum_local_ra2d


!=========================================================================
subroutine xsum_local_ra3d(array)
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

end subroutine xsum_local_ra3d


!=========================================================================
subroutine xsum_local_ra4d(array)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3*n4, MPI_DOUBLE_PRECISION, MPI_SUM, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_local_ra4d


!=========================================================================
subroutine xsum_local_ca1d(array)
 implicit none
 complex(dp),intent(inout) :: array(:)
!=====
 integer :: n1
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_local_ca1d


!=========================================================================
subroutine xsum_local_ca2d(array)
 implicit none
 complex(dp),intent(inout) :: array(:,:)
!=====
 integer :: n1,n2
 integer :: ier=0
!=====

 n1 = SIZE( array, DIM=1 )
 n2 = SIZE( array, DIM=2 )

#ifdef HAVE_MPI
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_local_ca2d


!=========================================================================
subroutine xsum_local_ca4d(array)
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
 call MPI_ALLREDUCE( MPI_IN_PLACE, array, n1*n2*n3*n4, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_local, ier)
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_allreduce'
 endif

end subroutine xsum_local_ca4d


!=========================================================================
subroutine xsum_local_procindex_ra2d(iproc,array)
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
 if( rank_local == iproc ) then
   call MPI_REDUCE( MPI_IN_PLACE, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, iproc, comm_local, ier)
 else
   call MPI_REDUCE( array, array, n1*n2, MPI_DOUBLE_PRECISION, MPI_SUM, iproc, comm_local, ier)
 endif
#endif
 if(ier/=0) then
   write(stdout,*) 'error in mpi_reduce'
 endif

end subroutine xsum_local_procindex_ra2d


!=========================================================================
end module m_mpi_local
!=========================================================================
