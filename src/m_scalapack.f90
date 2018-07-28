!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! a module with some SCALAPACK wrappers
! These wrappers are meant to be independant from the global variables
! defined in module m_mpi:
! All the required SCALAPACK proc grid information, block size are contained in a descriptor
!
! All the subroutine are meant to be valid with OR without SCALAPACK
! Do not forget the preprocessor flags!
!
!=========================================================================
module m_scalapack
 use m_definitions
 use m_warning
 use m_tools
 use m_mpi
#ifdef HAVE_MPI
 use mpi
#endif

 ! Indexes in the BLACS descriptor
 integer,parameter :: DTYPE_ = 1     ! Dense or Sparse (always dense here)
 integer,parameter :: CTXT_  = 2     ! Context
 integer,parameter :: M_     = 3     ! Number of rows
 integer,parameter :: N_     = 4     ! Number of cols
 integer,parameter :: MB_    = 5     ! Blocking factor for rows
 integer,parameter :: NB_    = 6     ! Blocking factor for cols
 integer,parameter :: RSRC_  = 7     ! First row
 integer,parameter :: CSRC_  = 8     ! First col
 integer,parameter :: LLD_   = 9     ! Number of rows in the local sub-matrix

 !
 ! SCALAPACK variables
 !
 integer,parameter :: NDEL=9
 integer,parameter :: block_row = 32
 integer,parameter :: block_col = 32
 integer,parameter :: first_row = 0
 integer,parameter :: first_col = 0

 ! Specific values for the MPI / SCALAPACK transposition
 integer,protected :: MBLOCK_AUXIL = 1
 integer,parameter :: NBLOCK_AUXIL = 1

 integer,protected :: nproc_sca = 1
 integer,protected :: iproc_sca = 0


 ! SCALAPACK grid: auxil distribution (so to mimic MPI distribution on auxiliary functions)
 integer,protected :: cntxt_auxil
 integer,protected :: nprow_auxil,npcol_auxil,iprow_auxil,ipcol_auxil

 ! SCALAPACK grid: square distribution
 integer,protected :: cntxt_3center
 integer,protected :: nprow_3center
 integer,protected :: npcol_3center
 integer,protected :: iprow_3center
 integer,protected :: ipcol_3center

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
 integer,protected :: desc_ham(NDEL)
 integer,public    :: desc_c(NDEL)
 integer,public    :: desc_small(NDEL)

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

 logical,protected :: parallel_ham       = .FALSE.

 ! Correspondence between SCALAPACK ham and MPI
 integer,allocatable,protected :: rank_sca_to_mpi(:,:)

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

 interface rowindex_global_to_local
   module procedure rowindex_global_to_local_distrib
   module procedure rowindex_global_to_local_descriptor
 end interface

 interface colindex_global_to_local
   module procedure colindex_global_to_local_distrib
   module procedure colindex_global_to_local_descriptor
 end interface

 interface diagonalize_sca
   module procedure diagonalize_inplace_sca_dp
   module procedure diagonalize_outofplace_sca_dp
   module procedure diagonalize_inplace_sca_cdp
 end interface

 interface diagonalize_scalapack
   module procedure diagonalize_scalapack_dp
   module procedure diagonalize_scalapack_cdp
 end interface

 interface create_distributed_copy
   module procedure create_distributed_copy_nospin_dp
   module procedure create_distributed_copy_spin_dp
   module procedure create_distributed_copy_nospin_cdp
 end interface create_distributed_copy

 interface gather_distributed_copy
   module procedure gather_distributed_copy_nospin_dp
   module procedure gather_distributed_copy_spin_dp
   module procedure gather_distributed_copy_nospin_cdp
 end interface gather_distributed_copy

 interface matmul_ab_scalapack
   module procedure matmul_ab_scalapack_dp
   module procedure matmul_ab_scalapack_cdp
 end interface matmul_ab_scalapack

 interface matmul_abc_scalapack
   module procedure matmul_abc_scalapack_dp
   module procedure matmul_abc_scalapack_cdp
 end interface matmul_abc_scalapack

 interface matmul_transaba_scalapack
   module procedure matmul_transaba_scalapack_dp
   module procedure matmul_transaba_scalapack_cdp
 end interface

#ifdef HAVE_SCALAPACK
 integer,external :: NUMROC,INDXL2G,INDXG2L,INDXG2P,PDLATRA,PDLAMCH
#endif


contains


#ifndef HAVE_SCALAPACK
!=========================================================================
! Fake SCALAPACK subroutines to be able to run without SCALAPACK
!
!=========================================================================
function NUMROC(n_in,idum1,idum2,idum3,idum4)
 implicit none
 integer,intent(in) :: n_in
 integer,intent(in) :: idum1,idum2,idum3,idum4
 integer            :: NUMROC
!=====

 NUMROC = n_in

end function NUMROC


!=========================================================================
function INDXL2G(indxloc, nb, iproc, isrcproc, nprocs )
 implicit none
 integer,intent(in)  :: indxloc, iproc, isrcproc, nb, nprocs
 integer             :: INDXL2G
!=====

  INDXL2G = nprocs * nb *( ( indxloc - 1 ) / nb ) + MOD( indxloc - 1 , nb ) &
             + MOD( nprocs + iproc - isrcproc , nprocs ) * nb + 1

end function INDXL2G


!=========================================================================
subroutine DESCINIT(desc,mmat,nmat,idum3,idum4,idum5,idum6,cntxtdum,idum7,info)
 implicit none
 integer,intent(inout) :: desc(NDEL)
 integer,intent(in)    :: mmat,nmat,idum3,idum4,idum5,idum6,idum7
 integer,intent(in)    :: cntxtdum
 integer,intent(out)   :: info
!=====

 desc(M_) = mmat
 desc(N_) = nmat

 info = 0

end subroutine DESCINIT


!=========================================================================
function INDXG2L(iglobal,idum1,idum2,idum3,idum4)
 implicit none
 integer,intent(in)  :: iglobal
 integer,intent(in)  :: idum1,idum2,idum3,idum4
 integer             :: INDXG2L
!=====

 INDXG2L = iglobal

end function INDXG2L


!=========================================================================
function INDXG2P(iglobal,idum1,idum2,idum3,idum4)
 implicit none
 integer,intent(in)  :: iglobal
 integer,intent(in)  :: idum1,idum2,idum3,idum4
 integer             :: INDXG2P
!=====

 INDXG2P = 0

end function INDXG2P


!=========================================================================
subroutine BLACS_GRIDINFO(icntxt,nprow,npcol,iprow,ipcol)
 implicit none
 integer,intent(in)  :: icntxt
 integer,intent(out) :: nprow,npcol,iprow,ipcol
!=====
!=====
 nprow = 1
 npcol = 1
 iprow = 0
 ipcol = 0

end subroutine BLACS_GRIDINFO


#endif


!=========================================================================
!
!  This is a pure function copy of the regular BLACS INDXL2G function
!  INDXL2G computes the global index of a distributed matrix entry
!  pointed to by the local index INDXLOC of the process indicated by
!  IPROC.
!
!=========================================================================
pure function indxl2g_pure(indxloc, nb, iproc, isrcproc, nprocs ) RESULT(indxl2g)
 implicit none
 integer,intent(in)  :: indxloc, iproc, isrcproc, nb, nprocs
 integer             :: indxl2g
!=====

  indxl2g = nprocs * nb *( ( indxloc - 1 ) / nb ) + MOD( indxloc - 1 , nb ) &
             + MOD( nprocs + iproc - isrcproc , nprocs ) * nb + 1

end function indxl2g_pure


!=========================================================================
! Create a distributed copy of a global matrix owned by each core
!
!=========================================================================
subroutine create_distributed_copy_nospin_dp(matrix_global,desc,matrix)
 implicit none
 integer,intent(in)   :: desc(NDEL)
 real(dp),intent(in)  :: matrix_global(:,:)
 real(dp),intent(out) :: matrix(:,:)
!=====
 integer              :: mlocal,nlocal
 integer              :: ilocal,jlocal,iglobal,jglobal
!=====

 mlocal = SIZE( matrix , DIM=1 )
 nlocal = SIZE( matrix , DIM=2 )

 if( desc(CTXT_) > 0 ) then
   do jlocal=1,nlocal
     jglobal = colindex_local_to_global(desc,jlocal)
     do ilocal=1,mlocal
       iglobal = rowindex_local_to_global(desc,ilocal)
       matrix(ilocal,jlocal) = matrix_global(iglobal,jglobal)
     enddo
   enddo
 endif


end subroutine create_distributed_copy_nospin_dp


!=========================================================================
! Create a distributed copy of a global matrix owned by each core
! with spin
!=========================================================================
subroutine create_distributed_copy_spin_dp(matrix_global,desc,matrix)
 implicit none
 integer,intent(in)   :: desc(NDEL)
 real(dp),intent(in)  :: matrix_global(:,:,:)
 real(dp),intent(out) :: matrix(:,:,:)
!=====
 integer              :: idim3,ndim3
 integer              :: mlocal,nlocal
 integer              :: ilocal,jlocal,iglobal,jglobal
!=====

 mlocal = SIZE( matrix , DIM=1 )
 nlocal = SIZE( matrix , DIM=2 )
 ndim3  = SIZE( matrix , DIM=3 )

 if( desc(CTXT_) > 0 ) then
   do idim3=1,ndim3
     do jlocal=1,nlocal
       jglobal = colindex_local_to_global(desc,jlocal)
       do ilocal=1,mlocal
         iglobal = rowindex_local_to_global(desc,ilocal)
         matrix(ilocal,jlocal,idim3) = matrix_global(iglobal,jglobal,idim3)
       enddo
     enddo
   enddo
 endif


end subroutine create_distributed_copy_spin_dp


!=========================================================================
! Create a distributed copy of a global matrix owned by each core (complex)
!
!=========================================================================
subroutine create_distributed_copy_nospin_cdp(matrix_global,desc,matrix)
 implicit none
 integer,intent(in)      :: desc(NDEL)
 complex(dp),intent(in)  :: matrix_global(:,:)
 complex(dp),intent(out) :: matrix(:,:)
!=====
 integer              :: mlocal,nlocal
 integer              :: ilocal,jlocal,iglobal,jglobal
!=====

 mlocal = SIZE( matrix , DIM=1 )
 nlocal = SIZE( matrix , DIM=2 )

 if( desc(CTXT_) > 0 ) then
   do jlocal=1,nlocal
     jglobal = colindex_local_to_global(desc,jlocal)
     do ilocal=1,mlocal
       iglobal = rowindex_local_to_global(desc,ilocal)
       matrix(ilocal,jlocal) = matrix_global(iglobal,jglobal)
     enddo
   enddo
 endif


end subroutine create_distributed_copy_nospin_cdp


!=========================================================================
! Gather a distributed matrix into a global matrix owned by each core
!
!=========================================================================
subroutine gather_distributed_copy_nospin_dp(desc,matrix,matrix_global)
 implicit none
 integer,intent(in)               :: desc(NDEL)
 real(dp),allocatable,intent(in)  :: matrix(:,:)
 real(dp),intent(out)             :: matrix_global(:,:)
!=====
 integer              :: cntxt
 integer              :: mlocal,nlocal,mglobal,nglobal
 integer              :: ilocal,jlocal,iglobal,jglobal
 integer              :: rank_master,iprow,ipcol,nprow,npcol
!=====

#ifdef HAVE_SCALAPACK

 cntxt = desc(CTXT_)
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )

 ! Find the master
 if( iprow == 0 .AND. ipcol == 0 ) then
   rank_master = rank_world
 else
   rank_master = -1
 endif
 call xmax_world(rank_master)

 if( cntxt > 0 ) then

   mlocal  = SIZE( matrix , DIM=1 )
   nlocal  = SIZE( matrix , DIM=2 )
   mglobal = SIZE( matrix_global , DIM=1 )
   nglobal = SIZE( matrix_global , DIM=2 )

   matrix_global(:,:) = 0.0_dp
   do jlocal=1,nlocal
     jglobal = colindex_local_to_global(desc,jlocal)
     do ilocal=1,mlocal
       iglobal = rowindex_local_to_global(desc,ilocal)
       matrix_global(iglobal,jglobal) = matrix(ilocal,jlocal)
     enddo
   enddo

   ! Only the master proc (0,0) gets the complete information
   call DGSUM2D(cntxt,'A',' ',mglobal,nglobal,matrix_global,nglobal,0,0)

 endif

 call xbcast_world(rank_master,matrix_global)

#else

 matrix_global(:,:) = matrix(:,:)

#endif

end subroutine gather_distributed_copy_nospin_dp


!=========================================================================
! Gather a distributed matrix into a global matrix owned by each core
!
!=========================================================================
subroutine gather_distributed_copy_spin_dp(desc,matrix,matrix_global)
 implicit none
 integer,intent(in)               :: desc(NDEL)
 real(dp),allocatable,intent(in)  :: matrix(:,:,:)
 real(dp),intent(out)             :: matrix_global(:,:,:)
!=====
 integer              :: cntxt
 integer              :: idim3,ndim3
 integer              :: mlocal,nlocal,mglobal,nglobal
 integer              :: ilocal,jlocal,iglobal,jglobal
 integer              :: rank_master,iprow,ipcol,nprow,npcol
!=====

#ifdef HAVE_SCALAPACK

 cntxt = desc(CTXT_)
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )

 ! Find the master
 if( iprow == 0 .AND. ipcol == 0 ) then
   rank_master = rank_world
 else
   rank_master = -1
 endif
 call xmax_world(rank_master)

 if( cntxt > 0 ) then
   mlocal  = SIZE( matrix , DIM=1 )
   nlocal  = SIZE( matrix , DIM=2 )
   mglobal = SIZE( matrix_global , DIM=1 )
   nglobal = SIZE( matrix_global , DIM=2 )
   ndim3   = SIZE( matrix , DIM=3 )

   matrix_global(:,:,:) = 0.0_dp
   do idim3=1,ndim3
     do jlocal=1,nlocal
       jglobal = colindex_local_to_global(desc,jlocal)
       do ilocal=1,mlocal
         iglobal = rowindex_local_to_global(desc,ilocal)
         matrix_global(iglobal,jglobal,idim3) = matrix(ilocal,jlocal,idim3)
       enddo
     enddo

     ! Only the master proc (0,0) gets the complete information
     call DGSUM2D(cntxt,'A',' ',mglobal,nglobal,matrix_global(1,1,idim3),nglobal,0,0)
   enddo

 endif
 call xbcast_world(rank_master,matrix_global)

#else
  matrix_global(:,:,:) = matrix(:,:,:)
#endif

end subroutine gather_distributed_copy_spin_dp


!=========================================================================
! Gather a distributed matrix into a global matrix owned by each core (complex)
!
!=========================================================================
subroutine gather_distributed_copy_nospin_cdp(desc,matrix,matrix_global)
 implicit none
 integer,intent(in)                  :: desc(NDEL)
 complex(dp),allocatable,intent(in)  :: matrix(:,:)
 complex(dp),intent(out)             :: matrix_global(:,:)
!=====
 integer              :: cntxt
 integer              :: mlocal,nlocal,mglobal,nglobal
 integer              :: ilocal,jlocal,iglobal,jglobal
 integer              :: rank_master,iprow,ipcol,nprow,npcol
!=====

#ifdef HAVE_SCALAPACK

 cntxt = desc(CTXT_)
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )

 ! Find the master
 if( iprow == 0 .AND. ipcol == 0 ) then
   rank_master = rank_world
 else
   rank_master = -1
 endif
 call xmax_world(rank_master)

 if( cntxt > 0 ) then

   mlocal  = SIZE( matrix , DIM=1 )
   nlocal  = SIZE( matrix , DIM=2 )
   mglobal = SIZE( matrix_global , DIM=1 )
   nglobal = SIZE( matrix_global , DIM=2 )

   matrix_global(:,:) = (0.0_dp, 0.0_dp)
   do jlocal=1,nlocal
     jglobal = colindex_local_to_global(desc,jlocal)
     do ilocal=1,mlocal
       iglobal = rowindex_local_to_global(desc,ilocal)
       matrix_global(iglobal,jglobal) = matrix(ilocal,jlocal)
     enddo
   enddo

   ! Only the master proc (0,0) gets the complete information
   call ZGSUM2D(cntxt,'A',' ',mglobal,nglobal,matrix_global,nglobal,0,0)

 endif

 call xbcast_world(rank_master,matrix_global)

#else

 matrix_global(:,:) = matrix(:,:)

#endif

end subroutine gather_distributed_copy_nospin_cdp


!=========================================================================
! Multiply on the left or on the right a distributed square matrix
! with a non-distributed diagonal
!=========================================================================
subroutine matmul_diag_sca(side,diag,desc,matrix)
 implicit none
 character(len=1),intent(in) :: side
 integer,intent(in)          :: desc(NDEL)
 real(dp),intent(in)         :: diag(:)
 real(dp),intent(inout)      :: matrix(:,:)
!=====
 integer                     :: nglobal
 integer                     :: mlocal,nlocal
 integer                     :: iglobal,jglobal,ilocal,jlocal
!=====

 nglobal = SIZE( diag(:) )
 mlocal  = SIZE( matrix , DIM=1 )
 nlocal  = SIZE( matrix , DIM=2 )


 select case(side)
 case('L')
#ifndef HAVE_SCALAPACK
   forall(iglobal=1:nglobal)
     matrix(iglobal,:) = matrix(iglobal,:) * diag(iglobal)
   end forall
#else
! Alternative coding
!   do iglobal=1,nglobal
!     call PDSCAL(nglobal,diag(iglobal),matrix,iglobal,1,desc,nglobal)
!   enddo
   do ilocal=1,mlocal
     iglobal = rowindex_local_to_global(desc,ilocal)
     matrix(ilocal,:) = matrix(ilocal,:) * diag(iglobal)
   enddo
#endif


 case('R')
#ifndef HAVE_SCALAPACK
   forall(jglobal=1:nglobal)
     matrix(:,jglobal) = matrix(:,jglobal) * diag(jglobal)
   end forall
#else
! Alternative coding
!   do jglobal=1,nglobal
!     call PDSCAL(nglobal,diag(jglobal),matrix,1,jglobal,desc,1)
!   enddo
   do jlocal=1,nlocal
     jglobal = colindex_local_to_global(desc,jlocal)
     matrix(:,jlocal) = matrix(:,jlocal) * diag(jglobal)
   enddo
#endif


 case default
   if( side /= 'L' .AND. side /= 'R' ) then
     call die('matmul_diag_sca: argument side should be L or R')
   endif
 end select


end subroutine matmul_diag_sca


!=========================================================================
! Diagonalize a REAL SYMMETRIC distributed matrix to get eigenvalues only
!
!=========================================================================
subroutine diagonalize_eigval_sca(nglobal,desc,matrix,eigval)
 implicit none
 integer,intent(in)     :: desc(NDEL),nglobal
 real(dp),intent(inout) :: matrix(:,:)
 real(dp),intent(out)   :: eigval(nglobal)
!=====
 integer              :: desc_eigvec(NDEL)
 integer              :: mlocal,nlocal
 integer              :: lwork,info
 real(dp),allocatable :: work(:)
 real(dp),allocatable :: eigvec(:,:)
 integer              :: neigval,neigvec
 integer,allocatable  :: iwork(:)
 integer              :: liwork
!=====

#ifdef HAVE_SCALAPACK
 desc_eigvec = desc

 mlocal = SIZE( matrix , DIM=1 )
 nlocal = SIZE( matrix , DIM=2 )

 allocate(eigvec(1,1))

 !
 ! First call to get the dimension of the array work
 lwork = -1
 allocate(work(1))
 call PDSYEV('N','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)


 !
 ! Second call to actually perform the diago
 lwork = NINT(work(1))

 deallocate(work)
 allocate(work(lwork))
 call PDSYEV('N','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)

 deallocate(work)


 deallocate(eigvec)

#else

 call diagonalize_wo_vectors(matrix,eigval)

#endif


end subroutine diagonalize_eigval_sca


!=========================================================================
! Diagonalize a REAL SYMMETRIC distributed matrix
!
!=========================================================================
subroutine diagonalize_inplace_sca_dp(nglobal,desc,matrix,eigval)
 implicit none
 integer,intent(in)     :: desc(NDEL),nglobal
 real(dp),intent(inout) :: matrix(:,:)
 real(dp),intent(out)   :: eigval(nglobal)
!=====
 integer              :: desc_eigvec(NDEL)
 integer              :: mlocal,nlocal
 integer              :: lwork,info
 real(dp),allocatable :: work(:)
 real(dp),allocatable :: eigvec(:,:)
 integer              :: neigval,neigvec
 integer,allocatable  :: iwork(:)
 integer              :: liwork
!=====

#ifdef HAVE_SCALAPACK
 desc_eigvec = desc

 mlocal = SIZE( matrix , DIM=1 )
 nlocal = SIZE( matrix , DIM=2 )

 allocate(eigvec(mlocal,nlocal))

 !
 ! First call to get the dimension of the array work
 lwork = -1
 allocate(work(1))
 call PDSYEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)


 !
 ! Second call to actually perform the diago
 lwork = NINT(work(1))

 deallocate(work)
 allocate(work(lwork))
 call PDSYEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)

 deallocate(work)


 matrix(:,:) = eigvec(:,:)

 deallocate(eigvec)

#else

 call diagonalize(matrix,eigval)

#endif


end subroutine diagonalize_inplace_sca_dp


!=========================================================================
! Diagonalize a COMPLEX HERMITIAN distributed matrix
!
!=========================================================================
subroutine diagonalize_inplace_sca_cdp(nglobal,desc,matrix,eigval)
 implicit none
 integer,intent(in)        :: desc(NDEL),nglobal
 complex(dp),intent(inout) :: matrix(:,:)
 real(dp),intent(out)      :: eigval(nglobal)
!=====
 integer              :: desc_eigvec(NDEL)
 integer              :: mlocal,nlocal
 integer              :: lwork,lrwork,info
 complex(dp),allocatable :: work(:),rwork(:)
 complex(dp),allocatable :: eigvec(:,:)
!=====

#ifdef HAVE_SCALAPACK
 desc_eigvec = desc

 mlocal = SIZE( matrix , DIM=1 )
 nlocal = SIZE( matrix , DIM=2 )

 allocate(eigvec(mlocal,nlocal))

 !
 ! First call to get the dimension of the array work
 lwork  = -1
 lrwork = -1
 allocate(work(1))
 allocate(rwork(1))
 call PZHEEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,rwork,lrwork,info)


 !
 ! Second call to actually perform the diago
 lwork  = NINT(REAL(work(1),dp))
 lrwork = NINT(REAL(rwork(1),dp))

 deallocate(work)
 allocate(work(lwork))
 call PZHEEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,rwork,lrwork,info)

 deallocate(work)


 matrix(:,:) = eigvec(:,:)

 deallocate(eigvec)

#else

 call diagonalize(matrix,eigval)

#endif


end subroutine diagonalize_inplace_sca_cdp


!=========================================================================
! Diagonalize a REAL SYMMETRIC distributed matrix
!
!=========================================================================
subroutine diagonalize_outofplace_sca_dp(nglobal,desc,matrix,eigval,desc_eigvec,eigvec)
 implicit none
 integer,intent(in)     :: nglobal
 integer,intent(in)     :: desc(NDEL)
 integer,intent(in)     :: desc_eigvec(NDEL)
 real(dp),intent(inout) :: matrix(:,:)
 real(dp),intent(out)   :: eigval(nglobal)
 real(dp),intent(out)   :: eigvec(:,:)
!=====
 integer              :: lwork,info
 real(dp),allocatable :: work(:)
 integer              :: neigval,neigvec
 integer,allocatable  :: iwork(:)
 integer              :: liwork
!=====

#ifdef HAVE_SCALAPACK

 !
 ! First call to get the dimension of the array work
 lwork = -1
 allocate(work(1))
 call PDSYEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)


 !
 ! Second call to actually perform the diago
 lwork = NINT(work(1))

 deallocate(work)
 allocate(work(lwork))
 call PDSYEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)

 deallocate(work)


#else

 call diagonalize(matrix,eigval,eigvec)

#endif


end subroutine diagonalize_outofplace_sca_dp


!=========================================================================
! Diagonalize a distributed matrix with PDSYEVR
!
!=========================================================================
subroutine diagonalize_sca_pdsyevr(nglobal,desc,matrix,eigval,desc_eigvec,eigvec)
 implicit none
 integer,intent(in)     :: nglobal
 integer,intent(in)     :: desc(NDEL)
 integer,intent(in)     :: desc_eigvec(NDEL)
 real(dp),intent(inout) :: matrix(:,:)
 real(dp),intent(out)   :: eigval(nglobal)
 real(dp),intent(out)   :: eigvec(:,:)
!=====
 integer              :: lwork,info
 real(dp),allocatable :: work(:)
 integer              :: neigval,neigvec
 integer,allocatable  :: iwork(:)
 integer              :: liwork
#ifdef SELECT_PDSYEVX
 real(dp)             :: ABSTOL
 integer              :: iclustr(2*nprow_sd*npcol_sd)
 real(dp)             :: gap(nprow_sd*npcol_sd)
 integer              :: ifail(nglobal)
#endif
!=====

#ifdef HAVE_SCALAPACK

 !
 ! First call to get the dimension of the array work
 lwork = -1
 allocate(work(3))
 liwork = -1
 allocate(iwork(1))
#ifdef SELECT_PDSYEVX
 ABSTOL = PDLAMCH(desc(CTXT_), 'U')
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
#else
 call PDSYEVR('V','A','L',nglobal,matrix,1,1,desc,0.0d0,0.0d0,0,0,neigval,neigvec,eigval,eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,info)
#endif


 !
 ! Second call to actually perform the diago
 lwork = NINT(work(1))

 deallocate(work)
 allocate(work(lwork))
 liwork = iwork(1)
 deallocate(iwork)
 allocate(iwork(liwork))
#ifdef SELECT_PDSYEVX
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
#else
 call PDSYEVR('V','A','L',nglobal,matrix,1,1,desc,0.0d0,0.0d0,0,0,neigval,neigvec,eigval,eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,info)
#endif

 deallocate(work)
 deallocate(iwork)


#else

 call diagonalize(matrix,eigval,eigvec)

#endif


end subroutine diagonalize_sca_pdsyevr


!=========================================================================
! Partially diagonalize a distributed matrix with PDSYEVR
!=========================================================================
subroutine diagonalize_sca_pdsyevr_partial(nglobal,neig,desc,matrix,eigval,desc_eigvec,eigvec)
 implicit none
 integer,intent(in)     :: nglobal,neig
 integer,intent(in)     :: desc(NDEL)
 integer,intent(in)     :: desc_eigvec(NDEL)
 real(dp),intent(inout) :: matrix(:,:)
 real(dp),intent(out)   :: eigval(nglobal)
 real(dp),intent(out)   :: eigvec(:,:)
!=====
 integer              :: lwork,info
 real(dp),allocatable :: work(:)
 integer              :: neigval,neigvec
 integer,allocatable  :: iwork(:)
 integer              :: liwork
 integer              :: ifail(nglobal)
 real(dp)             :: ABSTOL
#ifdef HAVE_SCALAPACK
 integer              :: iclustr(2*nprow_sd*npcol_sd)
 real(dp)             :: gap(nprow_sd*npcol_sd)
#endif
 real(dp),external    :: DLAMCH
!=====

#ifdef HAVE_SCALAPACK

 !
 ! First call to get the dimension of the array work
 lwork = -1
 allocate(work(3))
 liwork = -1
 allocate(iwork(1))
! call PDSYEVR('V','I','L',nglobal,matrix,1,1,desc,0.0d0,0.0d0,1,neig,neigval,neigvec,eigval,eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,info)
 ABSTOL = PDLAMCH(desc(CTXT_), 'U')
 call PDSYEVX('V','I','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,1,neig, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                     &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,           &
              ifail,iclustr,gap,info)


 !
 ! Second call to actually perform the diago
 lwork = NINT(work(1))

 deallocate(work)
 allocate(work(lwork))
 liwork = iwork(1)
 deallocate(iwork)
 allocate(iwork(liwork))
! call PDSYEVR('V','I','L',nglobal,matrix,1,1,desc,0.0d0,0.0d0,1,neig,neigval,neigvec,eigval,eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,info)
 call PDSYEVX('V','I','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,1,neig, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                     &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,           &
              ifail,iclustr,gap,info)

 deallocate(work)
 deallocate(iwork)

#else

 ABSTOL = 2.0_dp * DLAMCH('S')
 allocate(iwork(5*nglobal))
 lwork=-1
 allocate(work(1))
 call DSYEVX('V','I','L',nglobal,matrix,nglobal,0.d0,0.d0,1,neig,ABSTOL,neigval,eigval,eigvec,nglobal,work,lwork,iwork,ifail,info)

 lwork = NINT(work(1))
 deallocate(work)
 allocate(work(lwork))
 call DSYEVX('V','I','L',nglobal,matrix,nglobal,0.d0,0.d0,1,neig,ABSTOL,neigval,eigval,eigvec,nglobal,work,lwork,iwork,ifail,info)
 deallocate(work)
 deallocate(iwork)


#endif


end subroutine diagonalize_sca_pdsyevr_partial


!=========================================================================
! Diagonalize a non-distributed matrix
!
!=========================================================================
subroutine diagonalize_scalapack_dp(scalapack_block_min,nmat,matrix_global,eigval)
 implicit none
 integer,intent(in)     :: scalapack_block_min
 integer,intent(in)     :: nmat
 real(dp),intent(inout) :: matrix_global(nmat,nmat)
 real(dp),intent(out)   :: eigval(nmat)
!=====
 integer :: cntxt
 integer :: mlocal,nlocal
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
 integer :: iglobal,jglobal,ilocal,jlocal
 integer :: descm(NDEL),descz(NDEL)
 real(dp),allocatable :: matrix(:,:)
 integer :: rank_master
!=====

#ifdef HAVE_SCALAPACK
 call select_nprow_npcol(scalapack_block_min,nmat,nmat,nprow,npcol)

 if( nprow /= 1 .OR. npcol /= 1 ) then

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
   write(stdout,'(a,i4,a,i4)') ' Diagonalization using SCALAPACK with a grid',nprow,' x ',npcol

   ! Find the master
   if( iprow == 0 .AND. ipcol == 0 ) then
     rank_master = rank_world
   else
     rank_master = -1
   endif
   call xmax_world(rank_master)

   !
   ! Participate to the diagonalization only if the CPU has been selected
   ! in the grid
   if( cntxt > 0 ) then
     mlocal = NUMROC(nmat,block_row,iprow,first_row,nprow)
     nlocal = NUMROC(nmat,block_col,ipcol,first_col,npcol)

     allocate(matrix(mlocal,nlocal))

     call DESCINIT(descm,nmat,nmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mlocal),info)

     ! Set up the local copy of the matrix_global
     call create_distributed_copy(matrix_global,descm,matrix)

     call diagonalize_sca(nmat,descm,matrix,eigval)

   endif

   call gather_distributed_copy(descm,matrix,matrix_global)

   if( cntxt > 0 ) then
     deallocate(matrix)
     call BLACS_GRIDEXIT( cntxt )
   endif

   ! Then the master proc (0,0) broadcasts to all the others
   call xbcast_world(rank_master,eigval)


 else ! Only one SCALAPACK proc

   call diagonalize(matrix_global,eigval)

 endif

#else

 call diagonalize(matrix_global,eigval)

#endif


end subroutine diagonalize_scalapack_dp


!=========================================================================
! Diagonalize a non-distributed matrix
!
!=========================================================================
subroutine diagonalize_scalapack_cdp(scalapack_block_min,nmat,matrix_global,eigval)
 implicit none
 integer,intent(in)     :: scalapack_block_min
 integer,intent(in)     :: nmat
 complex(dp),intent(inout) :: matrix_global(nmat,nmat)
 real(dp),intent(out)   :: eigval(nmat)
!=====
 integer :: cntxt
 integer :: mlocal,nlocal
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
 integer :: iglobal,jglobal,ilocal,jlocal
 integer :: descm(NDEL),descz(NDEL)
 complex(dp),allocatable :: matrix(:,:)
 integer :: rank_master
!=====

#ifdef HAVE_SCALAPACK
 call select_nprow_npcol(scalapack_block_min,nmat,nmat,nprow,npcol)

 if( nprow /= 1 .OR. npcol /= 1 ) then

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
   write(stdout,'(a,i4,a,i4)') ' Diagonalization using SCALAPACK with a grid',nprow,' x ',npcol

   ! Find the master
   if( iprow == 0 .AND. ipcol == 0 ) then
     rank_master = rank_world
   else
     rank_master = -1
   endif
   call xmax_world(rank_master)

   !
   ! Participate to the diagonalization only if the CPU has been selected
   ! in the grid
   if( cntxt > 0 ) then
     mlocal = NUMROC(nmat,block_row,iprow,first_row,nprow)
     nlocal = NUMROC(nmat,block_col,ipcol,first_col,npcol)

     allocate(matrix(mlocal,nlocal))

     call DESCINIT(descm,nmat,nmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mlocal),info)

     ! Set up the local copy of the matrix_global
     call create_distributed_copy(matrix_global,descm,matrix)

     call diagonalize_sca(nmat,descm,matrix,eigval)

   endif

   call gather_distributed_copy(descm,matrix,matrix_global)

   if( cntxt > 0 ) then
     deallocate(matrix)
     call BLACS_GRIDEXIT( cntxt )
   endif

   ! Then the master proc (0,0) broadcasts to all the others
   call xbcast_world(rank_master,eigval)


 else ! Only one SCALAPACK proc

   call diagonalize(matrix_global,eigval)

 endif

#else

 call diagonalize(matrix_global,eigval)

#endif


end subroutine diagonalize_scalapack_cdp


!=========================================================================
! Calculate C = A * B for non-distributed REAL matrices
!
!=========================================================================
subroutine matmul_ab_scalapack_dp(scalapack_block_min,a_matrix,b_matrix,c_matrix)
 implicit none
 integer,intent(in)     :: scalapack_block_min
 real(dp),intent(in)    :: a_matrix(:,:)
 real(dp),intent(in)    :: b_matrix(:,:)
 real(dp),intent(out)   :: c_matrix(:,:)
!=====
 integer                :: mmat,nmat,kmat,lmat
 integer                :: kmat1,lmat1
 real(dp),allocatable   :: a_matrix_local(:,:)
 real(dp),allocatable   :: b_matrix_local(:,:)
 real(dp),allocatable   :: c_matrix_local(:,:)
 integer :: cntxt
 integer :: ma,na,mb,nb,mc,nc
 integer :: desca(NDEL),descb(NDEL),descc(NDEL)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
!=====

 mmat  = SIZE( a_matrix , DIM=1)
 kmat1 = SIZE( a_matrix , DIM=2)
 kmat  = SIZE( b_matrix , DIM=1)
 lmat1 = SIZE( b_matrix , DIM=2)
 lmat  = SIZE( c_matrix , DIM=1)
 nmat  = SIZE( c_matrix , DIM=2)

 if( kmat1 /= kmat ) call die('Dimension error in matmul_ab_scalapack_dp')
 if( mmat  /= lmat ) call die('Dimension error in matmul_ab_scalapack_dp')
 if( lmat1 /= nmat ) call die('Dimension error in matmul_ab_scalapack_dp')


#ifdef HAVE_SCALAPACK
 call select_nprow_npcol(scalapack_block_min,mmat,nmat,nprow,npcol)

 if( nprow /= 1 .OR. npcol /= 1 ) then

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
   write(stdout,'(a,i4,a,i4)') ' Matrix product using SCALAPACK with a grid',nprow,' x ',npcol

   !
   ! Participate to the calculation only if the CPU has been selected
   ! in the grid
   if( cntxt > 0 ) then

     !
     ! Distribute A
     ma = NUMROC(mmat,block_row,iprow,first_row,nprow)
     na = NUMROC(kmat1,block_col,ipcol,first_col,npcol)
     allocate(a_matrix_local(ma,na))
     call DESCINIT(desca,mmat,kmat1,block_row,block_col,first_row,first_col,cntxt,MAX(1,ma),info)
     ! Set up the local copy of the global matrix A
     call create_distributed_copy(a_matrix,desca,a_matrix_local)
     !
     ! Distribute B
     mb = NUMROC(kmat,block_row,iprow,first_row,nprow)
     nb = NUMROC(lmat1,block_col,ipcol,first_col,npcol)
     allocate(b_matrix_local(mb,nb))
     call DESCINIT(descb,kmat,lmat1,block_row,block_col,first_row,first_col,cntxt,MAX(1,mb),info)
     ! Set up the local copy of the global matrix B
     call create_distributed_copy(b_matrix,descb,b_matrix_local)
     !
     ! Prepare C = A * B
     mc = NUMROC(mmat,block_row,iprow,first_row,nprow)
     nc = NUMROC(nmat,block_col,ipcol,first_col,npcol)
     allocate(c_matrix_local(mc,nc))
     call DESCINIT(descc,mmat,nmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mc),info)

     ! Calculate C = A * B
     call PDGEMM('N','N',mmat,nmat,kmat,1.0_dp,a_matrix_local,1,1,desca,    &
                  b_matrix_local,1,1,descb,0.0_dp,c_matrix_local,1,1,descc)

     deallocate(a_matrix_local,b_matrix_local)


   endif

   call gather_distributed_copy(descc,c_matrix_local,c_matrix)

   if( cntxt > 0 ) then
     deallocate(c_matrix_local)
     call BLACS_GRIDEXIT( cntxt )
   endif


 else ! Only one SCALAPACK proc

!   c_matrix(:,:) = MATMUL( a_matrix , b_matrix )
   call DGEMM('N','N',mmat,nmat,kmat,1.0_dp,a_matrix,mmat,b_matrix,kmat,0.0_dp,c_matrix,mmat)

 endif

#else

! c_matrix(:,:) = MATMUL( a_matrix , b_matrix )
 call DGEMM('N','N',mmat,nmat,kmat,1.0_dp,a_matrix,mmat,b_matrix,kmat,0.0_dp,c_matrix,mmat)



#endif


end subroutine matmul_ab_scalapack_dp


!=========================================================================
! Calculate C = A * B for non-distributed COMPLEX matrices
!
!=========================================================================
subroutine matmul_ab_scalapack_cdp(scalapack_block_min,a_matrix,b_matrix,c_matrix)
 implicit none
 integer,intent(in)      :: scalapack_block_min
 complex(dp),intent(in)  :: a_matrix(:,:)
 complex(dp),intent(in)  :: b_matrix(:,:)
 complex(dp),intent(out) :: c_matrix(:,:)
!=====
 integer                   :: mmat,nmat,kmat,lmat
 integer                   :: kmat1,lmat1
 complex(dp),allocatable   :: a_matrix_local(:,:)
 complex(dp),allocatable   :: b_matrix_local(:,:)
 complex(dp),allocatable   :: c_matrix_local(:,:)
 integer :: cntxt
 integer :: ma,na,mb,nb,mc,nc
 integer :: desca(NDEL),descb(NDEL),descc(NDEL)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
 complex(dp),parameter :: ONE  = (1.0_dp,0.0_dp)
 complex(dp),parameter :: ZERO = (0.0_dp,0.0_dp)
!=====

 mmat  = SIZE( a_matrix , DIM=1)
 kmat1 = SIZE( a_matrix , DIM=2)
 kmat  = SIZE( b_matrix , DIM=1)
 lmat1 = SIZE( b_matrix , DIM=2)
 lmat  = SIZE( c_matrix , DIM=1)
 nmat  = SIZE( c_matrix , DIM=2)

 if( kmat1 /= kmat ) call die('Dimension error in matmul_ab_scalapack_cdp')
 if( mmat  /= lmat ) call die('Dimension error in matmul_ab_scalapack_cdp')
 if( lmat1 /= nmat ) call die('Dimension error in matmul_ab_scalapack_cdp')


#ifdef HAVE_SCALAPACK
 call select_nprow_npcol(scalapack_block_min,mmat,nmat,nprow,npcol)

 if( nprow /= 1 .OR. npcol /= 1 ) then

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
   write(stdout,'(a,i4,a,i4)') ' Matrix product using SCALAPACK with a grid',nprow,' x ',npcol

   !
   ! Participate to the calculation only if the CPU has been selected
   ! in the grid
   if( cntxt > 0 ) then

     !
     ! Distribute A
     ma = NUMROC(mmat,block_row,iprow,first_row,nprow)
     na = NUMROC(kmat1,block_col,ipcol,first_col,npcol)
     allocate(a_matrix_local(ma,na))
     call DESCINIT(desca,mmat,kmat1,block_row,block_col,first_row,first_col,cntxt,MAX(1,ma),info)
     ! Set up the local copy of the global matrix A
     call create_distributed_copy(a_matrix,desca,a_matrix_local)
     !
     ! Distribute B
     mb = NUMROC(kmat,block_row,iprow,first_row,nprow)
     nb = NUMROC(lmat1,block_col,ipcol,first_col,npcol)
     allocate(b_matrix_local(mb,nb))
     call DESCINIT(descb,kmat,lmat1,block_row,block_col,first_row,first_col,cntxt,MAX(1,mb),info)
     ! Set up the local copy of the global matrix B
     call create_distributed_copy(b_matrix,descb,b_matrix_local)
     !
     ! Prepare C = A * B
     mc = NUMROC(mmat,block_row,iprow,first_row,nprow)
     nc = NUMROC(nmat,block_col,ipcol,first_col,npcol)
     allocate(c_matrix_local(mc,nc))
     call DESCINIT(descc,mmat,nmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mc),info)

     ! Calculate C = A * B
     call PZGEMM('N','N',mmat,nmat,kmat,ONE,a_matrix_local,1,1,desca,    &
                  b_matrix_local,1,1,descb,ZERO,c_matrix_local,1,1,descc)

     deallocate(a_matrix_local,b_matrix_local)


   endif

   call gather_distributed_copy(descc,c_matrix_local,c_matrix)

   if( cntxt > 0 ) then
     deallocate(c_matrix_local)
     call BLACS_GRIDEXIT( cntxt )
   endif


 else ! Only one SCALAPACK proc

!   c_matrix(:,:) = MATMUL( a_matrix , b_matrix )
   call ZGEMM('N','N',mmat,nmat,kmat,ONE,a_matrix,mmat,b_matrix,kmat,ZERO,c_matrix,mmat)

 endif

#else

! c_matrix(:,:) = MATMUL( a_matrix , b_matrix )
 call ZGEMM('N','N',mmat,nmat,kmat,ONE,a_matrix,mmat,b_matrix,kmat,ZERO,c_matrix,mmat)



#endif


end subroutine matmul_ab_scalapack_cdp


!=========================================================================
! Calculate D = A * B * C for non-distributed matrices
!
!=========================================================================
subroutine matmul_abc_scalapack_dp(scalapack_block_min,a_matrix,b_matrix,c_matrix,d_matrix)
 implicit none
 integer,intent(in)     :: scalapack_block_min
 real(dp),intent(in)    :: a_matrix(:,:)
 real(dp),intent(in)    :: b_matrix(:,:)
 real(dp),intent(in)    :: c_matrix(:,:)
 real(dp),intent(out)   :: d_matrix(:,:)
!=====
 integer                :: mmat,nmat,kmat,lmat
 integer                :: mmat1,nmat1,kmat1,lmat1
 real(dp),allocatable   :: a_matrix_local(:,:)
 real(dp),allocatable   :: b_matrix_local(:,:)
 real(dp),allocatable   :: c_matrix_local(:,:)
 real(dp),allocatable   :: d_matrix_local(:,:)
 real(dp),allocatable   :: m_matrix_local(:,:)
 real(dp),allocatable   :: m_matrix(:,:)
 integer :: cntxt
 integer :: ma,na,mb,nb,mc,nc,md,nd,mm,nm
 integer :: desca(NDEL),descb(NDEL),descc(NDEL),descd(NDEL)
 integer :: descm(NDEL)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
!=====

 mmat  = SIZE( a_matrix , DIM=1)
 kmat1 = SIZE( a_matrix , DIM=2)
 kmat  = SIZE( b_matrix , DIM=1)
 lmat1 = SIZE( b_matrix , DIM=2)
 lmat  = SIZE( c_matrix , DIM=1)
 nmat1 = SIZE( c_matrix , DIM=2)
 mmat1 = SIZE( d_matrix , DIM=1)
 nmat  = SIZE( d_matrix , DIM=2)

 if( mmat1 /= mmat ) call die('Dimension error in matmul_abc_scalapack_dp')
 if( nmat1 /= nmat ) call die('Dimension error in matmul_abc_scalapack_dp')
 if( kmat1 /= kmat ) call die('Dimension error in matmul_abc_scalapack_dp')
 if( lmat1 /= lmat ) call die('Dimension error in matmul_abc_scalapack_dp')


#ifdef HAVE_SCALAPACK
 call select_nprow_npcol(scalapack_block_min,mmat,nmat,nprow,npcol)

 if( nprow /= 1 .OR. npcol /= 1 ) then

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
   write(stdout,'(a,i4,a,i4)') ' Matrix product using SCALAPACK with a grid',nprow,' x ',npcol


   !
   ! Participate to the diagonalization only if the CPU has been selected
   ! in the grid
   if( cntxt > 0 ) then

     !
     ! Distribute A
     ma = NUMROC(mmat,block_row,iprow,first_row,nprow)
     na = NUMROC(kmat,block_col,ipcol,first_col,npcol)
     allocate(a_matrix_local(ma,na))
     call DESCINIT(desca,mmat,kmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,ma),info)
     ! Set up the local copy of the global matrix A
     call create_distributed_copy(a_matrix,desca,a_matrix_local)
     !
     ! Distribute B
     mb = NUMROC(kmat,block_row,iprow,first_row,nprow)
     nb = NUMROC(lmat,block_col,ipcol,first_col,npcol)
     allocate(b_matrix_local(mb,nb))
     call DESCINIT(descb,kmat,lmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mb),info)
     ! Set up the local copy of the global matrix B
     call create_distributed_copy(b_matrix,descb,b_matrix_local)
     !
     ! Prepare M = A * B
     mm = NUMROC(mmat,block_row,iprow,first_row,nprow)
     nm = NUMROC(lmat,block_col,ipcol,first_col,npcol)
     allocate(m_matrix_local(mm,nm))
     call DESCINIT(descm,mmat,lmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mm),info)

     ! Calculate M = A * B
     call PDGEMM('N','N',mmat,lmat,kmat,1.0_dp,a_matrix_local,1,1,desca,    &
                  b_matrix_local,1,1,descb,0.0_dp,m_matrix_local,1,1,descm)

     deallocate(a_matrix_local,b_matrix_local)

     !
     ! Distribute C
     mc = NUMROC(lmat,block_row,iprow,first_row,nprow)
     nc = NUMROC(nmat,block_col,ipcol,first_col,npcol)
     allocate(c_matrix_local(mc,nc))
     call DESCINIT(descc,lmat,nmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mc),info)
     ! Set up the local copy of the global matrix C
     call create_distributed_copy(c_matrix,descc,c_matrix_local)

     !
     ! Prepare D = M * C
     md = NUMROC(mmat,block_row,iprow,first_row,nprow)
     nd = NUMROC(nmat,block_col,ipcol,first_col,npcol)
     allocate(d_matrix_local(md,nd))
     call DESCINIT(descd,mmat,nmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,md),info)

     ! Calculate D = M * C
     call PDGEMM('N','N',mmat,nmat,lmat,1.0_dp,m_matrix_local,1,1,descm,    &
                  c_matrix_local,1,1,descc,0.0_dp,d_matrix_local,1,1,descd)

     deallocate(m_matrix_local,c_matrix_local)

   endif

   call gather_distributed_copy(descd,d_matrix_local,d_matrix)

   if( cntxt > 0 ) then
     deallocate(d_matrix_local)
     call BLACS_GRIDEXIT( cntxt )
   endif


 else ! Only one SCALAPACK proc

   allocate(m_matrix(mmat,lmat))

!   m_matrix(:,:) = MATMUL( a_matrix , b_matrix )
!   d_matrix(:,:) = MATMUL( m_matrix , c_matrix )
   call DGEMM('N','N',mmat,lmat,kmat,1.0_dp,a_matrix,mmat,b_matrix,kmat,0.0_dp,m_matrix,mmat)
   call DGEMM('N','N',mmat,nmat,lmat,1.0_dp,m_matrix,mmat,c_matrix,lmat,0.0_dp,d_matrix,mmat)


   deallocate(m_matrix)

 endif

#else

 allocate(m_matrix(mmat,lmat))

! m_matrix(:,:) = MATMUL( a_matrix , b_matrix )
! d_matrix(:,:) = MATMUL( m_matrix , c_matrix )
 call DGEMM('N','N',mmat,lmat,kmat,1.0_dp,a_matrix,mmat,b_matrix,kmat,0.0_dp,m_matrix,mmat)
 call DGEMM('N','N',mmat,nmat,lmat,1.0_dp,m_matrix,mmat,c_matrix,lmat,0.0_dp,d_matrix,mmat)

 deallocate(m_matrix)


#endif


end subroutine matmul_abc_scalapack_dp


!=========================================================================
! Calculate D = A * B * C for non-distributed COMPLEX matrices
!
!=========================================================================
subroutine matmul_abc_scalapack_cdp(scalapack_block_min,a_matrix,b_matrix,c_matrix,d_matrix)
 implicit none
 integer,intent(in)      :: scalapack_block_min
 complex(dp),intent(in)  :: a_matrix(:,:)
 complex(dp),intent(in)  :: b_matrix(:,:)
 complex(dp),intent(in)  :: c_matrix(:,:)
 complex(dp),intent(out) :: d_matrix(:,:)
!=====
 integer                 :: mmat,nmat,kmat,lmat
 integer                 :: mmat1,nmat1,kmat1,lmat1
 complex(dp),allocatable :: a_matrix_local(:,:)
 complex(dp),allocatable :: b_matrix_local(:,:)
 complex(dp),allocatable :: c_matrix_local(:,:)
 complex(dp),allocatable :: d_matrix_local(:,:)
 complex(dp),allocatable :: m_matrix_local(:,:)
 complex(dp),allocatable :: m_matrix(:,:)
 integer :: cntxt
 integer :: ma,na,mb,nb,mc,nc,md,nd,mm,nm
 integer :: desca(NDEL),descb(NDEL),descc(NDEL),descd(NDEL)
 integer :: descm(NDEL)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
 complex(dp),parameter :: ONE  = (1.0_dp,0.0_dp)
 complex(dp),parameter :: ZERO = (0.0_dp,0.0_dp)
!=====

 mmat  = SIZE( a_matrix , DIM=1)
 kmat1 = SIZE( a_matrix , DIM=2)
 kmat  = SIZE( b_matrix , DIM=1)
 lmat1 = SIZE( b_matrix , DIM=2)
 lmat  = SIZE( c_matrix , DIM=1)
 nmat1 = SIZE( c_matrix , DIM=2)
 mmat1 = SIZE( d_matrix , DIM=1)
 nmat  = SIZE( d_matrix , DIM=2)

 if( mmat1 /= mmat ) call die('Dimension error in matmul_abc_scalapack_cdp')
 if( nmat1 /= nmat ) call die('Dimension error in matmul_abc_scalapack_cdp')
 if( kmat1 /= kmat ) call die('Dimension error in matmul_abc_scalapack_cdp')
 if( lmat1 /= lmat ) call die('Dimension error in matmul_abc_scalapack_cdp')


#ifdef HAVE_SCALAPACK
 call select_nprow_npcol(scalapack_block_min,mmat,nmat,nprow,npcol)

 if( nprow /= 1 .OR. npcol /= 1 ) then

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
   write(stdout,'(a,i4,a,i4)') ' Matrix product using SCALAPACK with a grid',nprow,' x ',npcol


   !
   ! Participate to the diagonalization only if the CPU has been selected
   ! in the grid
   if( cntxt > 0 ) then

     !
     ! Distribute A
     ma = NUMROC(mmat,block_row,iprow,first_row,nprow)
     na = NUMROC(kmat,block_col,ipcol,first_col,npcol)
     allocate(a_matrix_local(ma,na))
     call DESCINIT(desca,mmat,kmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,ma),info)
     ! Set up the local copy of the global matrix A
     call create_distributed_copy(a_matrix,desca,a_matrix_local)
     !
     ! Distribute B
     mb = NUMROC(kmat,block_row,iprow,first_row,nprow)
     nb = NUMROC(lmat,block_col,ipcol,first_col,npcol)
     allocate(b_matrix_local(mb,nb))
     call DESCINIT(descb,kmat,lmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mb),info)
     ! Set up the local copy of the global matrix B
     call create_distributed_copy(b_matrix,descb,b_matrix_local)
     !
     ! Prepare M = A * B
     mm = NUMROC(mmat,block_row,iprow,first_row,nprow)
     nm = NUMROC(lmat,block_col,ipcol,first_col,npcol)
     allocate(m_matrix_local(mm,nm))
     call DESCINIT(descm,mmat,lmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mm),info)

     ! Calculate M = A * B
     call PZGEMM('N','N',mmat,lmat,kmat,ONE,a_matrix_local,1,1,desca,    &
                  b_matrix_local,1,1,descb,ZERO,m_matrix_local,1,1,descm)

     deallocate(a_matrix_local,b_matrix_local)

     !
     ! Distribute C
     mc = NUMROC(lmat,block_row,iprow,first_row,nprow)
     nc = NUMROC(nmat,block_col,ipcol,first_col,npcol)
     allocate(c_matrix_local(mc,nc))
     call DESCINIT(descc,lmat,nmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mc),info)
     ! Set up the local copy of the global matrix C
     call create_distributed_copy(c_matrix,descc,c_matrix_local)

     !
     ! Prepare D = M * C
     md = NUMROC(mmat,block_row,iprow,first_row,nprow)
     nd = NUMROC(nmat,block_col,ipcol,first_col,npcol)
     allocate(d_matrix_local(md,nd))
     call DESCINIT(descd,mmat,nmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,md),info)

     ! Calculate D = M * C
     call PZGEMM('N','N',mmat,nmat,lmat,ONE,m_matrix_local,1,1,descm,    &
                  c_matrix_local,1,1,descc,ZERO,d_matrix_local,1,1,descd)

     deallocate(m_matrix_local,c_matrix_local)

   endif

   call gather_distributed_copy(descd,d_matrix_local,d_matrix)

   if( cntxt > 0 ) then
     deallocate(d_matrix_local)
     call BLACS_GRIDEXIT( cntxt )
   endif


 else ! Only one SCALAPACK proc

   allocate(m_matrix(mmat,lmat))

!   m_matrix(:,:) = MATMUL( a_matrix , b_matrix )
!   d_matrix(:,:) = MATMUL( m_matrix , c_matrix )
   call ZGEMM('N','N',mmat,lmat,kmat,ONE,a_matrix,mmat,b_matrix,kmat,ZERO,m_matrix,mmat)
   call ZGEMM('N','N',mmat,nmat,lmat,ONE,m_matrix,mmat,c_matrix,lmat,ZERO,d_matrix,mmat)


   deallocate(m_matrix)

 endif

#else

 allocate(m_matrix(mmat,lmat))

! m_matrix(:,:) = MATMUL( a_matrix , b_matrix )
! d_matrix(:,:) = MATMUL( m_matrix , c_matrix )
 call ZGEMM('N','N',mmat,lmat,kmat,ONE,a_matrix,mmat,b_matrix,kmat,ZERO,m_matrix,mmat)
 call ZGEMM('N','N',mmat,nmat,lmat,ONE,m_matrix,mmat,c_matrix,lmat,ZERO,d_matrix,mmat)

 deallocate(m_matrix)


#endif


end subroutine matmul_abc_scalapack_cdp


!=========================================================================
! Calculate C = A^T * B * A for non-distributed REAL matrices
!
!=========================================================================
subroutine matmul_transaba_scalapack_dp(scalapack_block_min,a_matrix,b_matrix,c_matrix)
 implicit none
 integer,intent(in)     :: scalapack_block_min
 real(dp),intent(in)    :: a_matrix(:,:)
 real(dp),intent(in)    :: b_matrix(:,:)
 real(dp),intent(out)   :: c_matrix(:,:)
!=====
 integer                :: mmat,kmat
 integer                :: mmat1,kmat1
 integer                :: mmat2,kmat2
 real(dp),allocatable   :: a_matrix_local(:,:)
 real(dp),allocatable   :: b_matrix_local(:,:)
 real(dp),allocatable   :: c_matrix_local(:,:)
 real(dp),allocatable   :: m_matrix_local(:,:)
 real(dp),allocatable   :: m_matrix(:,:)
 integer :: cntxt
 integer :: ma,na,mb,nb,mc,nc,mm,nm
 integer :: desca(NDEL),descb(NDEL),descc(NDEL)
 integer :: descm(NDEL)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
!=====

 kmat  = SIZE( a_matrix , DIM=1)
 mmat  = SIZE( a_matrix , DIM=2)
 kmat1 = SIZE( b_matrix , DIM=1)
 kmat2 = SIZE( b_matrix , DIM=2)
 mmat1 = SIZE( c_matrix , DIM=1)
 mmat2 = SIZE( c_matrix , DIM=2)

 if( mmat1 /= mmat ) then
   write(msg,*) 'mmat1 /= mmat',mmat1,mmat
   call die('Dimension error in matmul_transaba_scalapack_dp'//msg)
 endif
 if( mmat2 /= mmat ) then
   write(msg,*) 'mmat2 /= mmat',mmat2,mmat
   call die('Dimension error in matmul_transaba_scalapack_dp'//msg)
 endif
 if( kmat1 /= kmat ) then
   write(msg,*) 'kmat1 /= kmat',kmat1,kmat
   call die('Dimension error in matmul_transaba_scalapack_dp'//msg)
 endif
 if( kmat2 /= kmat ) then
   write(msg,*) 'kmat2 /= kmat',kmat2,kmat
   call die('Dimension error in matmul_transaba_scalapack_dp'//msg)
 endif


#ifdef HAVE_SCALAPACK
 call select_nprow_npcol(scalapack_block_min,mmat,mmat,nprow,npcol)

 if( nprow /= 1 .OR. npcol /= 1 ) then

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
   write(stdout,'(a,i4,a,i4)') ' Matrix product using SCALAPACK with a grid',nprow,' x ',npcol

   !
   ! Participate to the diagonalization only if the CPU has been selected
   ! in the grid
   if( cntxt > 0 ) then

     !
     ! Distribute A
     ma = NUMROC(kmat,block_row,iprow,first_row,nprow)
     na = NUMROC(mmat,block_col,ipcol,first_col,npcol)
     allocate(a_matrix_local(ma,na))
     call DESCINIT(desca,kmat,mmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,ma),info)
     ! Set up the local copy of the global matrix A
     call create_distributed_copy(a_matrix,desca,a_matrix_local)
     !
     ! Distribute B
     mb = NUMROC(kmat,block_row,iprow,first_row,nprow)
     nb = NUMROC(kmat,block_col,ipcol,first_col,npcol)
     allocate(b_matrix_local(mb,nb))
     call DESCINIT(descb,kmat,kmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mb),info)
     ! Set up the local copy of the global matrix B
     call create_distributed_copy(b_matrix,descb,b_matrix_local)
     !
     ! Prepare M = A^T * B
     mm = NUMROC(mmat,block_row,iprow,first_row,nprow)
     nm = NUMROC(kmat,block_col,ipcol,first_col,npcol)
     allocate(m_matrix_local(mm,nm))
     call DESCINIT(descm,mmat,kmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mm),info)

     ! Calculate M = A^T * B
     call PDGEMM('T','N',mmat,kmat,kmat,1.0_dp,a_matrix_local,1,1,desca,    &
                  b_matrix_local,1,1,descb,0.0_dp,m_matrix_local,1,1,descm)

     deallocate(b_matrix_local)

     !
     ! Prepare C = M * A
     mc = NUMROC(mmat,block_row,iprow,first_row,nprow)
     nc = NUMROC(mmat,block_col,ipcol,first_col,npcol)
     allocate(c_matrix_local(mc,nc))
     call DESCINIT(descc,mmat,mmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mc),info)

     ! Calculate C = M * A
     call PDGEMM('N','N',mmat,mmat,kmat,1.0_dp,m_matrix_local,1,1,descm,    &
                  a_matrix_local,1,1,desca,0.0_dp,c_matrix_local,1,1,descc)

     deallocate(m_matrix_local,a_matrix_local)

   endif

   call gather_distributed_copy(descc,c_matrix_local,c_matrix)

   if( cntxt > 0 ) then
     deallocate(c_matrix_local)
     call BLACS_GRIDEXIT( cntxt )
   endif


 else ! Only one SCALAPACK proc

   allocate(m_matrix(mmat,kmat))

!   m_matrix(:,:) = MATMUL( TRANSPOSE(a_matrix) , b_matrix )
!   c_matrix(:,:) = MATMUL( m_matrix , a_matrix )
   call DGEMM('T','N',mmat,kmat,kmat,1.0_dp,a_matrix,kmat,b_matrix,kmat,0.0_dp,m_matrix,mmat)
   call DGEMM('N','N',mmat,mmat,kmat,1.0_dp,m_matrix,mmat,a_matrix,kmat,0.0_dp,c_matrix,mmat)

   deallocate(m_matrix)

 endif

#else

 allocate(m_matrix(mmat,kmat))

! m_matrix(:,:) = MATMUL( TRANSPOSE(a_matrix) , b_matrix )
! c_matrix(:,:) = MATMUL( m_matrix , a_matrix )
 call DGEMM('T','N',mmat,kmat,kmat,1.0_dp,a_matrix,kmat,b_matrix,kmat,0.0_dp,m_matrix,mmat)
 call DGEMM('N','N',mmat,mmat,kmat,1.0_dp,m_matrix,mmat,a_matrix,kmat,0.0_dp,c_matrix,mmat)

 deallocate(m_matrix)


#endif


end subroutine matmul_transaba_scalapack_dp


!=========================================================================
! Calculate C = A^T * B * A for non-distributed REAL matrices
!
!=========================================================================
subroutine matmul_transaba_scalapack_cdp(scalapack_block_min,a_matrix,b_matrix,c_matrix)
 implicit none
 integer,intent(in)      :: scalapack_block_min
 complex(dp),intent(in)  :: a_matrix(:,:)
 complex(dp),intent(in)  :: b_matrix(:,:)
 complex(dp),intent(out) :: c_matrix(:,:)
!=====
 complex(dp),parameter   :: ONE  = (1.0_dp,0.0_dp)
 complex(dp),parameter   :: ZERO = (0.0_dp,0.0_dp)
 integer                 :: mmat,kmat
 integer                 :: mmat1,kmat1
 integer                 :: mmat2,kmat2
 complex(dp),allocatable :: a_matrix_local(:,:)
 complex(dp),allocatable :: b_matrix_local(:,:)
 complex(dp),allocatable :: c_matrix_local(:,:)
 complex(dp),allocatable :: m_matrix_local(:,:)
 complex(dp),allocatable :: m_matrix(:,:)
 integer :: cntxt
 integer :: ma,na,mb,nb,mc,nc,mm,nm
 integer :: desca(NDEL),descb(NDEL),descc(NDEL)
 integer :: descm(NDEL)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
!=====

 kmat  = SIZE( a_matrix , DIM=1)
 mmat  = SIZE( a_matrix , DIM=2)
 kmat1 = SIZE( b_matrix , DIM=1)
 kmat2 = SIZE( b_matrix , DIM=2)
 mmat1 = SIZE( c_matrix , DIM=1)
 mmat2 = SIZE( c_matrix , DIM=2)

 if( mmat1 /= mmat ) then
   write(msg,*) 'mmat1 /= mmat',mmat1,mmat
   call die('Dimension error in matmul_transaba_scalapack_cdp'//msg)
 endif
 if( mmat2 /= mmat ) then
   write(msg,*) 'mmat2 /= mmat',mmat2,mmat
   call die('Dimension error in matmul_transaba_scalapack_cdp'//msg)
 endif
 if( kmat1 /= kmat ) then
   write(msg,*) 'kmat1 /= kmat',kmat1,kmat
   call die('Dimension error in matmul_transaba_scalapack_cdp'//msg)
 endif
 if( kmat2 /= kmat ) then
   write(msg,*) 'kmat2 /= kmat',kmat2,kmat
   call die('Dimension error in matmul_transaba_scalapack_cdp'//msg)
 endif


#ifdef HAVE_SCALAPACK
 call select_nprow_npcol(scalapack_block_min,mmat,mmat,nprow,npcol)

 if( nprow /= 1 .OR. npcol /= 1 ) then

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
   write(stdout,'(a,i4,a,i4)') ' Matrix product using SCALAPACK with a grid',nprow,' x ',npcol

   !
   ! Participate to the diagonalization only if the CPU has been selected
   ! in the grid
   if( cntxt > 0 ) then

     !
     ! Distribute A
     ma = NUMROC(kmat,block_row,iprow,first_row,nprow)
     na = NUMROC(mmat,block_col,ipcol,first_col,npcol)
     allocate(a_matrix_local(ma,na))
     call DESCINIT(desca,kmat,mmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,ma),info)
     ! Set up the local copy of the global matrix A
     call create_distributed_copy(a_matrix,desca,a_matrix_local)
     !
     ! Distribute B
     mb = NUMROC(kmat,block_row,iprow,first_row,nprow)
     nb = NUMROC(kmat,block_col,ipcol,first_col,npcol)
     allocate(b_matrix_local(mb,nb))
     call DESCINIT(descb,kmat,kmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mb),info)
     ! Set up the local copy of the global matrix B
     call create_distributed_copy(b_matrix,descb,b_matrix_local)
     !
     ! Prepare M = A^T * B
     mm = NUMROC(mmat,block_row,iprow,first_row,nprow)
     nm = NUMROC(kmat,block_col,ipcol,first_col,npcol)
     allocate(m_matrix_local(mm,nm))
     call DESCINIT(descm,mmat,kmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mm),info)

     ! Calculate M = A^T * B
     call PZGEMM('T','N',mmat,kmat,kmat,ONE,a_matrix_local,1,1,desca,    &
                  b_matrix_local,1,1,descb,ZERO,m_matrix_local,1,1,descm)

     deallocate(b_matrix_local)

     !
     ! Prepare C = M * A
     mc = NUMROC(mmat,block_row,iprow,first_row,nprow)
     nc = NUMROC(mmat,block_col,ipcol,first_col,npcol)
     allocate(c_matrix_local(mc,nc))
     call DESCINIT(descc,mmat,mmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mc),info)

     ! Calculate C = M * A
     call PZGEMM('N','N',mmat,mmat,kmat,ONE,m_matrix_local,1,1,descm,    &
                  a_matrix_local,1,1,desca,ZERO,c_matrix_local,1,1,descc)

     deallocate(m_matrix_local,a_matrix_local)

   endif

   call gather_distributed_copy(descc,c_matrix_local,c_matrix)

   if( cntxt > 0 ) then
     deallocate(c_matrix_local)
     call BLACS_GRIDEXIT( cntxt )
   endif


 else ! Only one SCALAPACK proc

   allocate(m_matrix(mmat,kmat))

!   m_matrix(:,:) = MATMUL( TRANSPOSE(a_matrix) , b_matrix )
!   c_matrix(:,:) = MATMUL( m_matrix , a_matrix )
   call ZGEMM('T','N',mmat,kmat,kmat,ONE,a_matrix,kmat,b_matrix,kmat,ZERO,m_matrix,mmat)
   call ZGEMM('N','N',mmat,mmat,kmat,ONE,m_matrix,mmat,a_matrix,kmat,ZERO,c_matrix,mmat)

   deallocate(m_matrix)

 endif

#else

 allocate(m_matrix(mmat,kmat))

! m_matrix(:,:) = MATMUL( TRANSPOSE(a_matrix) , b_matrix )
! c_matrix(:,:) = MATMUL( m_matrix , a_matrix )
 call ZGEMM('T','N',mmat,kmat,kmat,ONE,a_matrix,kmat,b_matrix,kmat,ZERO,m_matrix,mmat)
 call ZGEMM('N','N',mmat,mmat,kmat,ONE,m_matrix,mmat,a_matrix,kmat,ZERO,c_matrix,mmat)

 deallocate(m_matrix)


#endif


end subroutine matmul_transaba_scalapack_cdp


!=========================================================================
! Calculate Trace ( A^T * B ) for non-distributed matrices
!
!=========================================================================
subroutine trace_transab_scalapack(scalapack_block_min,a_matrix,b_matrix,ab_trace)
 use m_tools,only: matrix_trace
 implicit none
 integer,intent(in)     :: scalapack_block_min
 real(dp),intent(in)    :: a_matrix(:,:)
 real(dp),intent(in)    :: b_matrix(:,:)
 real(dp),intent(out)   :: ab_trace
!=====
 integer                :: mmat1,kmat1
 integer                :: mmat2,kmat2
 real(dp),allocatable   :: a_matrix_local(:,:)
 real(dp),allocatable   :: b_matrix_local(:,:)
 real(dp),allocatable   :: m_matrix_local(:,:)
 real(dp),allocatable   :: m_matrix(:,:)
 integer :: cntxt
 integer :: ma,na,mb,nb,mm,nm
 integer :: desca(NDEL),descb(NDEL)
 integer :: descm(NDEL)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
!=====

 kmat1 = SIZE( a_matrix , DIM=1)
 kmat2 = SIZE( a_matrix , DIM=2)
 mmat1 = SIZE( b_matrix , DIM=1)
 mmat2 = SIZE( b_matrix , DIM=2)

 if( mmat1 /= kmat1 ) then
   write(msg,*) 'kmat1 /= mmat1',kmat1,mmat1
   call die('Dimension error in trace_transab_scalapack'//msg)
 endif
 if( mmat2 /= kmat2 ) then
   write(msg,*) 'mmat2 /= kmat2',mmat2,kmat2
   call die('Dimension error in trace_transab_scalapack'//msg)
 endif

#ifdef HAVE_SCALAPACK
 call select_nprow_npcol(scalapack_block_min,kmat1,kmat2,nprow,npcol)

 if( nprow /= 1 .OR. npcol /= 1 ) then

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
   write(stdout,'(a,i4,a,i4)') ' Trace of matrix product using SCALAPACK with a grid',nprow,' x ',npcol

   !
   ! Participate to the calculation only if the CPU has been selected
   ! in the grid
   if( cntxt > 0 ) then

     !
     ! Distribute A
     ma = NUMROC(kmat1,block_row,iprow,first_row,nprow)
     na = NUMROC(kmat2,block_col,ipcol,first_col,npcol)
     allocate(a_matrix_local(ma,na))
     call DESCINIT(desca,kmat1,kmat2,block_row,block_col,first_row,first_col,cntxt,MAX(1,ma),info)
     ! Set up the local copy of the global matrix A
     call create_distributed_copy(a_matrix,desca,a_matrix_local)
     !
     ! Distribute B
     mb = NUMROC(mmat1,block_row,iprow,first_row,nprow)
     nb = NUMROC(mmat2,block_col,ipcol,first_col,npcol)
     allocate(b_matrix_local(mb,nb))
     call DESCINIT(descb,mmat1,mmat2,block_row,block_col,first_row,first_col,cntxt,MAX(1,mb),info)
     ! Set up the local copy of the global matrix B
     call create_distributed_copy(b_matrix,descb,b_matrix_local)
     !
     ! Prepare M = A^T * B
     mm = NUMROC(kmat2,block_row,iprow,first_row,nprow)
     nm = NUMROC(mmat2,block_col,ipcol,first_col,npcol)
     allocate(m_matrix_local(mm,nm))
     call DESCINIT(descm,kmat2,mmat2,block_row,block_col,first_row,first_col,cntxt,MAX(1,mm),info)

     ! Calculate M = A^T * B
     call PDGEMM('T','N',kmat2,mmat2,mmat1,1.0_dp,a_matrix_local,1,1,desca,    &
                  b_matrix_local,1,1,descb,0.0_dp,m_matrix_local,1,1,descm)

     deallocate(a_matrix_local)
     deallocate(b_matrix_local)

     ab_trace = PDLATRA(kmat2,m_matrix_local,1,1,descm) / REAL(nprow*npcol,dp)

     call BLACS_GRIDEXIT( cntxt )

   else
     ab_trace = 0.0_dp
   endif

   call xsum_world(ab_trace)



 else ! Only one SCALAPACK proc

   allocate(m_matrix(kmat2,mmat2))

   m_matrix(:,:) = MATMUL( TRANSPOSE(a_matrix) , b_matrix )
   ab_trace = matrix_trace(m_matrix)

! FIXME the following coding
!   call DGEMM('T','N',mmat,kmat,kmat,1.0_dp,a_matrix,kmat,b_matrix,kmat,0.0_dp,m_matrix,mmat)
!   call DGEMM('N','N',mmat,mmat,kmat,1.0_dp,m_matrix,mmat,a_matrix,kmat,0.0_dp,c_matrix,mmat)

   deallocate(m_matrix)

 endif

#else

 allocate(m_matrix(kmat2,mmat2))

 m_matrix(:,:) = MATMUL( TRANSPOSE(a_matrix) , b_matrix )
 ab_trace = matrix_trace(m_matrix)

! call DGEMM('T','N',mmat,kmat,kmat,1.0_dp,a_matrix,kmat,b_matrix,kmat,0.0_dp,m_matrix,mmat)
! call DGEMM('N','N',mmat,mmat,kmat,1.0_dp,m_matrix,mmat,a_matrix,kmat,0.0_dp,c_matrix,mmat)

 deallocate(m_matrix)


#endif


end subroutine trace_transab_scalapack


!=========================================================================
! Calculate D = A * B * C for distributed matrices
!
!=========================================================================
subroutine matmul_abc_sca(desca,a_matrix_local,descb,b_matrix_local,descc,c_matrix_local,descd,d_matrix_local)
 implicit none
 integer ,intent(in)  :: desca(NDEL),descb(NDEL),descc(NDEL),descd(NDEL)
 real(dp),intent(in)  :: a_matrix_local(:,:)
 real(dp),intent(in)  :: b_matrix_local(:,:)
 real(dp),intent(in)  :: c_matrix_local(:,:)
 real(dp),intent(out) :: d_matrix_local(:,:)
!=====
 integer              :: mmat,nmat,kmat,lmat
 integer              :: mmat1,nmat1,kmat1,lmat1
 real(dp),allocatable :: m_matrix_local(:,:)
 real(dp),allocatable :: m_matrix(:,:)
 integer :: cntxt
 integer :: ma,na,mb,nb,mc,nc,md,nd,mm,nm
 integer :: descm(NDEL)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
!=====

 mmat  = desca(M_)
 kmat1 = desca(N_)
 kmat  = descb(M_)
 lmat1 = descb(N_)
 lmat  = descc(M_)
 nmat1 = descc(N_)
 mmat1 = descd(M_)
 nmat  = descd(N_)

 if( mmat1 /= mmat ) call die('Dimension error in matmul_abc_scalapack')
 if( nmat1 /= nmat ) call die('Dimension error in matmul_abc_scalapack')
 if( kmat1 /= kmat ) call die('Dimension error in matmul_abc_scalapack')
 if( lmat1 /= lmat ) call die('Dimension error in matmul_abc_scalapack')


#ifdef HAVE_SCALAPACK

 cntxt = desca(CTXT_)
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
 write(stdout,'(a,i4,a,i4)') ' Matrix product A * B * C using SCALAPACK with a grid',nprow,' x ',npcol

 !
 ! Participate to the diagonalization only if the CPU has been selected
 ! in the grid
 if( cntxt > 0 ) then

   !
   ! Prepare M = A * B
   mm = NUMROC(mmat,block_row,iprow,first_row,nprow)
   nm = NUMROC(lmat,block_col,ipcol,first_col,npcol)
   allocate(m_matrix_local(mm,nm))
   call DESCINIT(descm,mmat,lmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mm),info)

   ! Calculate M = A * B
   call PDGEMM('N','N',mmat,lmat,kmat,1.0_dp,a_matrix_local,1,1,desca,    &
                b_matrix_local,1,1,descb,0.0_dp,m_matrix_local,1,1,descm)

   ! Calculate D = M * C
   call PDGEMM('N','N',mmat,nmat,lmat,1.0_dp,m_matrix_local,1,1,descm,    &
                c_matrix_local,1,1,descc,0.0_dp,d_matrix_local,1,1,descd)

   deallocate(m_matrix_local)

 endif

#else
  d_matrix_local(:,:) = MATMUL( a_matrix_local(:,:) , MATMUL( b_matrix_local(:,:) , c_matrix_local(:,:) ) )
#endif


end subroutine matmul_abc_sca


!=========================================================================
! Calculate C = A^T * B * A for distributed matrices
!
!=========================================================================
subroutine matmul_transaba_sca(desca,a_matrix_local,descb,b_matrix_local,descc,c_matrix_local)
 implicit none
 integer,intent(in)     :: desca(NDEL),descb(NDEL),descc(NDEL)
 real(dp),intent(in)    :: a_matrix_local(:,:)
 real(dp),intent(in)    :: b_matrix_local(:,:)
 real(dp),intent(out)   :: c_matrix_local(:,:)
!=====
 integer                :: mmat,kmat
 integer                :: mmat1,kmat1
 integer                :: mmat2,kmat2
 real(dp),allocatable   :: m_matrix_local(:,:)
 integer :: cntxt
 integer :: ma,na,mb,nb,mc,nc,mm,nm
 integer :: descm(NDEL)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
!=====


 kmat  = desca(M_)
 mmat  = desca(N_)
 kmat1 = descb(M_)
 kmat2 = descb(N_)
 mmat1 = descc(M_)
 mmat2 = descc(N_)

 if( mmat1 /= mmat ) then
   write(msg,*) 'mmat1 /= mmat',mmat1,mmat
   call die('Dimension error in matmul_transaba_scalapack'//msg)
 endif
 if( mmat2 /= mmat ) then
   write(msg,*) 'mmat2 /= mmat',mmat2,mmat
   call die('Dimension error in matmul_transaba_scalapack'//msg)
 endif
 if( kmat1 /= kmat ) then
   write(msg,*) 'kmat1 /= kmat',kmat1,kmat
   call die('Dimension error in matmul_transaba_scalapack'//msg)
 endif
 if( kmat2 /= kmat ) then
   write(msg,*) 'kmat2 /= kmat',kmat2,kmat
   call die('Dimension error in matmul_transaba_scalapack'//msg)
 endif

#ifdef HAVE_SCALAPACK
 cntxt = desca(CTXT_)
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
 write(stdout,'(a,i4,a,i4)') ' Matrix product A**T * B * A using SCALAPACK with a grid',nprow,' x ',npcol

 !
 ! Participate to the diagonalization only if the CPU has been selected
 ! in the grid
 if( cntxt > 0 ) then

   !
   ! Prepare M = A^T * B
   mm = NUMROC(mmat,block_row,iprow,first_row,nprow)
   nm = NUMROC(kmat,block_col,ipcol,first_col,npcol)
   allocate(m_matrix_local(mm,nm))
   call DESCINIT(descm,mmat,kmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mm),info)

   ! Calculate M = A^T * B
   call PDGEMM('T','N',mmat,kmat,kmat,1.0_dp,a_matrix_local,1,1,desca,    &
                b_matrix_local,1,1,descb,0.0_dp,m_matrix_local,1,1,descm)

   ! Calculate C = M * A
   call PDGEMM('N','N',mmat,mmat,kmat,1.0_dp,m_matrix_local,1,1,descm,    &
                a_matrix_local,1,1,desca,0.0_dp,c_matrix_local,1,1,descc)

   deallocate(m_matrix_local)

 endif

#else
  c_matrix_local(:,:) = MATMUL( TRANSPOSE(a_matrix_local) , MATMUL( b_matrix_local , a_matrix_local ) )
#endif


end subroutine matmul_transaba_sca


!=========================================================================
! Transform a lower or upper triangular matrix into the full representation
!
!=========================================================================
subroutine symmetrize_matrix_sca(uplo,nglobal,desc,matrix,desc_tmp,matrix_tmp)
 implicit none
 character(len=1)       :: uplo
 integer,intent(in)     :: nglobal
 integer,intent(in)     :: desc(NDEL),desc_tmp(NDEL)
 real(dp),intent(inout) :: matrix(:,:)
 real(dp),intent(inout) :: matrix_tmp(:,:)
!=====
 integer                :: mlocal,nlocal
 integer                :: iglobal,jglobal,ilocal,jlocal
!=====

#ifndef HAVE_SCALAPACK

 select case(uplo)
 case('L')
   do jglobal=1,nglobal
     do iglobal=jglobal+1,nglobal
       matrix(jglobal,iglobal) = matrix(iglobal,jglobal)
     enddo
   enddo

 case('U')
   do iglobal=1,nglobal
     do jglobal=iglobal+1,nglobal
       matrix(jglobal,iglobal) = matrix(iglobal,jglobal)
     enddo
   enddo

 case default
   if( uplo /= 'L' .AND. uplo /= 'U' ) then
     call die('symmetrize_matrix_sca: argument uplo should be L or U')
   endif
 end select

#else

 mlocal = SIZE( matrix , DIM=1 )
 nlocal = SIZE( matrix , DIM=2 )

 select case(uplo)
 case('L')
   ! Symmetrize M (lower triangular matrix)
   ! by adding M^T (upper triangular)
   ! be careful about the diagonal terms
   do jlocal=1,nlocal
     jglobal = colindex_local_to_global(desc,jlocal)
     do ilocal=1,mlocal
       iglobal = rowindex_local_to_global(desc,ilocal)

       ! Half the diagonal and erase the upper triangle
       if( iglobal == jglobal ) then
         matrix(ilocal,jlocal) = matrix(ilocal,jlocal) * 0.5_dp
       else if( iglobal < jglobal ) then
         matrix(ilocal,jlocal) = 0.0_dp
       endif

     enddo
   enddo

 case('U')
   ! Symmetrize M (upper triangular matrix)
   ! by adding M^T (lower triangular)
   ! be careful about the diagonal terms
   do jlocal=1,nlocal
     jglobal = colindex_local_to_global(desc,jlocal)
     do ilocal=1,mlocal
       iglobal = rowindex_local_to_global(desc,ilocal)

       ! Half the diagonal and erase the upper triangle
       if( iglobal == jglobal ) then
         matrix(ilocal,jlocal) = matrix(ilocal,jlocal) * 0.5_dp
       else if( iglobal > jglobal ) then
         matrix(ilocal,jlocal) = 0.0_dp
       endif

     enddo
   enddo

 case default
   if( uplo /= 'L' .AND. uplo /= 'U' ) then
     call die('symmetrize_matrix_sca: argument uplo should be L or U')
   endif
 end select

 call PDLACPY('A',nglobal,nglobal,matrix,1,1,desc,matrix_tmp,1,1,desc_tmp)
 call PDGEADD('T',nglobal,nglobal,1.d0,matrix_tmp,1,1,desc,1.d0,matrix,1,1,desc)

#endif


end subroutine symmetrize_matrix_sca


!=========================================================================
subroutine invert_sca(desc,matrix,matrix_inv)
 implicit none

 integer,intent(in)   :: desc(NDEL)
 real(dp),intent(in)  :: matrix(:,:)
 real(dp),intent(out) :: matrix_inv(:,:)
!=====
 real(dp),allocatable :: work(:)
 integer ,allocatable :: iwork(:)
 integer ,allocatable :: ipiv(:)
 integer  :: cntxt
 integer  :: nprow,npcol,iprow,ipcol
 integer  :: mlocal,nlocal
 integer  :: nmat
 integer  :: info
 integer  :: lwork,liwork
!=====

#ifdef HAVE_SCALAPACK
 cntxt = desc(CTXT_)
 call BLACS_GRIDINFO(cntxt,nprow,npcol,iprow,ipcol)

 if( nprow * npcol > 1 ) then

   if( iprow >= 0 ) then

     nmat = desc(M_)
     mlocal = SIZE( matrix , DIM=1 )
     nlocal = SIZE( matrix , DIM=2 )
     matrix_inv(:,:) = matrix(:,:)

     allocate(ipiv(block_row+mlocal))

     call PDGETRF(nmat,nmat,matrix_inv,1,1,desc,ipiv,info)

     if(info/=0) call die('FAILURE in PDGETRF')

     ! Query
     allocate(work(1),iwork(1))
     lwork  = -1
     liwork = -1
     call PDGETRI(nmat,matrix_inv,1,1,desc,ipiv,work, lwork, iwork, liwork, info )
     if(info/=0) call die('FAILURE in PDGETRI')

     lwork  = NINT(work(1))
     liwork = iwork(1)
     deallocate(work,iwork)
     allocate(work(lwork),iwork(liwork))
     call PDGETRI(nmat,matrix_inv,1,1,desc,ipiv,work, lwork, iwork, liwork, info )
     if(info/=0) call die('FAILURE in PDGETRI')

     deallocate(ipiv)
     deallocate(work,iwork)

   endif

 else
   call invert(matrix,matrix_inv)
 endif

#else

 call invert(matrix,matrix_inv)

#endif

end subroutine invert_sca


!=========================================================================
subroutine init_scalapack()
 implicit none

!=====
!=====

#ifdef HAVE_SCALAPACK

 ! Get iproc_sca and nproc_sca
 call BLACS_PINFO( iproc_sca, nproc_sca )

 ! Set nprow, npcol

 ! Squared division of tasks
 nprow_sd = INT(SQRT(REAL(nproc_sca,dp)))
 npcol_sd = nproc_sca / nprow_sd
 do while( npcol_sd <= nproc_sca .AND. nprow_sd >= 1 )
   ! Found a correct distribution
   if( nprow_sd * npcol_sd == nproc_sca ) exit
   npcol_sd = npcol_sd + 1
   nprow_sd = nproc_sca / npcol_sd
 enddo

 if( npcol_sd / nprow_sd > 8 ) then
   call issue_warning('SCALAPACK distribution of processors is much rectangular. ' &
                      // 'This may affect the performance. Try to change the number of processors')
 endif

 call BLACS_GET( -1, 0, cntxt_sd )
 call BLACS_GRIDINIT( cntxt_sd, 'R', nprow_sd, npcol_sd )
 call BLACS_GRIDINFO( cntxt_sd, nprow_sd, npcol_sd, iprow_sd, ipcol_sd )

 ! 3center integrals distribution
 nprow_3center = nprow_sd
 npcol_3center = npcol_sd
 call BLACS_GET( -1, 0, cntxt_3center )
 call BLACS_GRIDINIT( cntxt_3center, 'R', nprow_3center, npcol_3center )
 call BLACS_GRIDINFO( cntxt_3center, nprow_3center, npcol_3center, iprow_3center, ipcol_3center )

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

#else
 ! If SCALAPACK is not present, fill the variable with fake values
 nprow_sd = 1
 npcol_sd = 1
 iprow_sd = 0
 ipcol_sd = 0
 cntxt_3center = 1
 nprow_3center = 1
 npcol_3center = 1
 ipcol_3center = 0
 iprow_3center = 0
#endif

end subroutine init_scalapack


!=========================================================================
subroutine init_scalapack_other(nbf,scalapack_nprow,scalapack_npcol,m_ham,n_ham)
 implicit none

 integer,intent(in)  :: nbf
 integer,intent(in)  :: scalapack_nprow,scalapack_npcol
 integer,intent(out) :: m_ham,n_ham
!=====
 integer :: ier=0
 integer :: color
 integer             :: iproc_auxil
 integer,allocatable :: usermap(:,:)
#ifdef DEBUG
 integer          :: fileunit
 character(len=3) :: ctmp
#endif
!=====

#ifdef HAVE_SCALAPACK
 parallel_ham = scalapack_nprow * scalapack_npcol > 1

 if( parallel_ham ) then

   nprow_ham = scalapack_nprow
   npcol_ham = scalapack_npcol
   if( nprow_ham * npcol_ham > nproc_world ) call die('SCALAPACK manual distribution asks for too many processors')

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

   write(stdout,'(/,a)')           ' ==== SCALAPACK Hamiltonian'
   write(stdout,'(a50,1x,i8)')      'Number of dedicated processors:',nprow_ham * npcol_ham
   write(stdout,'(a50,1x,i8,1x,i8)')   'Grid of dedicated processors:',nprow_ham,npcol_ham

   ! Distribute the remaing procs for auxiliary basis and grid points
   color = MODULO( rank_world , nprow_ham * npcol_ham )
   call MPI_COMM_SPLIT(comm_world,color,rank_world,comm_local,ier);
   call MPI_COMM_SIZE(comm_local,nproc_local,ier)
   call MPI_COMM_RANK(comm_local,rank_local,ier)

   write(stdout,'(a50,1x,i8)')      'Number of local processors:',nproc_local

   call xmax_local(m_ham)
   call xmax_local(n_ham)
   call xmax_local(iprow_ham)
   call xmax_local(ipcol_ham)


   ! Define the transversal communicator
   color = rank_world / ( nprow_ham * npcol_ham )

   call MPI_COMM_SPLIT(comm_world,color,rank_world,comm_trans,ier);
   call MPI_COMM_SIZE(comm_trans,nproc_trans,ier)
   call MPI_COMM_RANK(comm_trans,rank_trans,ier)

   allocate(rank_sca_to_mpi(0:nprow_ham-1,0:npcol_ham-1))
   rank_sca_to_mpi(:,:) = -1
   rank_sca_to_mpi(iprow_ham,ipcol_ham) = rank_trans
   call xmax_trans(rank_sca_to_mpi)


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

#ifdef HAVE_SCALAPACK
 !
 ! Create the SCALAPACK context cntxt_auxil
 ! that precisely matches the MPI_COMMUNICATOR comm_auxil
 !
 call BLACS_GET( -1, 0, cntxt_auxil )

 if( rank_world /= iproc_sca ) then
   call die('init_mpi_other_communicators: coding is valid only if SCALAPACK and MPI order the procs in the same manner')
 endif

 allocate(usermap(nproc_auxil,1))
 do iproc_auxil=0,nproc_auxil-1
   usermap(iproc_auxil+1,1) = iproc_auxil * nproc_ortho
 enddo
 call BLACS_GRIDMAP(cntxt_auxil,usermap,nproc_auxil,nproc_auxil,1)
 deallocate(usermap)
#endif

 call BLACS_GRIDINFO(cntxt_auxil,nprow_auxil,npcol_auxil,iprow_auxil,ipcol_auxil)
 call xmax_ortho(nprow_auxil)
 call xmax_ortho(npcol_auxil)
 call xmax_ortho(iprow_auxil)
 call xmax_ortho(ipcol_auxil)


#ifdef DEBUG
 write(ctmp,'(i3.3)') rank_world
 open(newunit=fileunit,file='DEBUG_mpiinfo_rank'//TRIM(ctmp))
 write(fileunit,*) 'nproc_world:',nproc_world
 write(fileunit,*) 'rank_world:',rank_world
 write(fileunit,*)
 write(fileunit,*) 'nproc_local:',nproc_local
 write(fileunit,*) 'rank_local:',rank_local
 write(fileunit,*)
 write(fileunit,*) 'nproc_trans:',nproc_trans
 write(fileunit,*) 'rank_trans:',rank_trans
 write(fileunit,*)
 write(fileunit,*) 'nproc_auxil:',nproc_auxil
 write(fileunit,*) 'rank_auxil:',rank_auxil
 write(fileunit,*)
 write(fileunit,*) 'nproc_grid:',nproc_grid
 write(fileunit,*) 'rank_grid:',rank_grid
 close(fileunit)
#endif


end subroutine init_scalapack_other


!=========================================================================
function row_block_size(mglobal,iprow,nprow)
 implicit none

 integer,intent(in) :: mglobal,iprow,nprow
 integer            :: row_block_size
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
 integer,intent(out)         :: desc(NDEL),mlocal,nlocal
!=====
 integer :: info
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
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   call die('BUG')
 end select

 write(stdout,'(/,a,i6,a,i6,4x,i6)') ' SCALAPACK info: size of the local matrix for proc #', mlocal,' x ',nlocal,iproc_sca

#else
 desc(:)     = 0
 desc(M_)   = mglobal
 desc(N_)   = nglobal
 desc(LLD_) = mglobal
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
subroutine max_matrix_sca_sd(mmat,nmat,real_matrix)
 implicit none
 integer, intent(in)    :: mmat,nmat
 real(dp),intent(inout) :: real_matrix(mmat,nmat)
!=====
 integer :: idum1(1),idum2(1)
!=====

#ifdef HAVE_SCALAPACK
 call dgamx2d(cntxt_sd,'All',' ',mmat,nmat,real_matrix,mmat,-1,idum1,idum2,-1,-1)
#endif


end subroutine max_matrix_sca_sd


!=========================================================================
function rowindex_global_to_local_distrib(distribution,iglobal) result(ilocal)
 implicit none
 character(len=1),intent(in) :: distribution
 integer,intent(in)          :: iglobal
 integer                     :: ilocal
!=====
!=====
 !
 ! returns the local index if this is proc in charge
 ! else returns 0

#ifdef HAVE_SCALAPACK
 select case(distribution)
 case('S')
   if( iprow_sd == INDXG2P(iglobal,block_row,0,first_row,nprow_sd) ) then
     ilocal = INDXG2L(iglobal,block_row,0,first_row,nprow_sd)
   else
     ilocal = 0
   endif
 case('H')
   if( iprow_ham == INDXG2P(iglobal,block_row,0,first_row,nprow_ham) ) then
     ilocal = INDXG2L(iglobal,block_row,0,first_row,nprow_ham)
   else
     ilocal = 0
   endif
 case('R')
   if( iprow_rd == INDXG2P(iglobal,block_row,0,first_row,nprow_rd) ) then
     ilocal = INDXG2L(iglobal,block_row,0,first_row,nprow_rd)
   else
     ilocal = 0
   endif
 case('C')
   if( iprow_cd == INDXG2P(iglobal,block_row,0,first_row,nprow_cd) ) then
     ilocal = INDXG2L(iglobal,block_row,0,first_row,nprow_cd)
   else
     ilocal = 0
   endif
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   call die('BUG')
 end select
#else
 ilocal = iglobal
#endif

end function rowindex_global_to_local_distrib


!=========================================================================
function colindex_global_to_local_distrib(distribution,jglobal) result(jlocal)
 implicit none
 character(len=1),intent(in) :: distribution
 integer,intent(in)          :: jglobal
 integer                     :: jlocal
!=====
!=====
 !
 ! returns the local index if this is proc in charge
 ! else returns 0

#ifdef HAVE_SCALAPACK
 select case(distribution)
 case('S')
   if( ipcol_sd == INDXG2P(jglobal,block_col,0,first_col,npcol_sd) ) then
     jlocal = INDXG2L(jglobal,block_col,0,first_col,npcol_sd)
   else
     jlocal = 0
   endif
 case('H')
   if( ipcol_ham == INDXG2P(jglobal,block_col,0,first_col,npcol_ham) ) then
     jlocal = INDXG2L(jglobal,block_col,0,first_col,npcol_ham)
   else
     jlocal = 0
   endif
 case('R')
   if( ipcol_rd == INDXG2P(jglobal,block_col,0,first_col,npcol_rd) ) then
     jlocal = INDXG2L(jglobal,block_col,0,first_col,npcol_rd)
   else
     jlocal = 0
   endif
 case('C')
   if( ipcol_cd == INDXG2P(jglobal,block_col,0,first_col,npcol_cd) ) then
     jlocal = INDXG2L(jglobal,block_col,0,first_col,npcol_cd)
   else
     jlocal = 0
   endif
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   call die('BUG')
 end select
#else
 jlocal = jglobal
#endif

end function colindex_global_to_local_distrib


!=========================================================================
function rowindex_local_to_global_distrib(distribution,ilocal) result(iglobal)
 implicit none
 character(len=1),intent(in) :: distribution
 integer,intent(in)          :: ilocal
 integer                     :: iglobal
!=====
!=====

#ifdef HAVE_SCALAPACK
 select case(distribution)
 case('S')
   iglobal = INDXL2G(ilocal,block_row,iprow_sd,first_row,nprow_sd)
 case('H')
   iglobal = INDXL2G(ilocal,block_row,iprow_ham,first_row,nprow_ham)
 case('R')
   iglobal = INDXL2G(ilocal,block_row,iprow_rd,first_row,nprow_rd)
 case('C')
   iglobal = INDXL2G(ilocal,block_row,iprow_cd,first_row,nprow_cd)
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   call die('BUG')
 end select
#else
 iglobal = ilocal
#endif

end function rowindex_local_to_global_distrib


!=========================================================================
function rowindex_local_to_global_procindex(iprow,nprow,ilocal) result(iglobal)
 implicit none
 integer,intent(in)          :: iprow,nprow,ilocal
 integer                     :: iglobal
!=====
!=====

#ifdef HAVE_SCALAPACK
 iglobal = INDXL2G(ilocal,block_row,iprow,first_row,nprow)
#else
 iglobal = ilocal
#endif

end function rowindex_local_to_global_procindex


!=========================================================================
function rowindex_local_to_global_descriptor(desc,ilocal) result(iglobal)
 implicit none
 integer,intent(in)          :: desc(NDEL),ilocal
 integer                     :: iglobal
!=====
#ifdef HAVE_SCALAPACK
 integer          :: iprow,ipcol,nprow,npcol
#endif
!=====

#ifdef HAVE_SCALAPACK
 call BLACS_GRIDINFO(desc(CTXT_),nprow,npcol,iprow,ipcol)
 iglobal = INDXL2G(ilocal,desc(MB_),iprow,desc(RSRC_),nprow)
#else
 iglobal = ilocal
#endif

end function rowindex_local_to_global_descriptor


!=========================================================================
function rowindex_global_to_local_descriptor(desc,iglobal) result(ilocal)
 implicit none
 integer,intent(in)          :: desc(NDEL),iglobal
 integer                     :: ilocal
!=====
#ifdef HAVE_SCALAPACK
 integer          :: iprow,ipcol,nprow,npcol
#endif
!=====

#ifdef HAVE_SCALAPACK
 call BLACS_GRIDINFO(desc(CTXT_),nprow,npcol,iprow,ipcol)
 if( iprow == INDXG2P(iglobal,desc(MB_),0,desc(RSRC_),nprow) ) then
   ilocal = INDXG2L(iglobal,desc(MB_),0,desc(RSRC_),nprow)
 else
   ilocal = 0
 endif
#else
 ilocal = iglobal
#endif

end function rowindex_global_to_local_descriptor


!=========================================================================
function colindex_local_to_global_distrib(distribution,jlocal) result(jglobal)
 implicit none
 character(len=1),intent(in) :: distribution
 integer,intent(in)          :: jlocal
 integer                     :: jglobal
!=====
!=====

#ifdef HAVE_SCALAPACK
 select case(distribution)
 case('S')
   jglobal = INDXL2G(jlocal,block_col,ipcol_sd,first_col,npcol_sd)
 case('H')
   jglobal = INDXL2G(jlocal,block_col,ipcol_ham,first_col,npcol_ham)
 case('R')
   jglobal = INDXL2G(jlocal,block_col,ipcol_rd,first_col,npcol_rd)
 case('C')
   jglobal = INDXL2G(jlocal,block_col,ipcol_cd,first_col,npcol_cd)
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   call die('BUG')
 end select
#else
 jglobal = jlocal
#endif

end function colindex_local_to_global_distrib


!=========================================================================
function colindex_local_to_global_procindex(ipcol,npcol,jlocal) result(jglobal)
 implicit none
 integer,intent(in)          :: ipcol,npcol,jlocal
 integer                     :: jglobal
!=====
!=====

#ifdef HAVE_SCALAPACK
 jglobal = INDXL2G(jlocal,block_col,ipcol,first_col,npcol)
#else
 jglobal = jlocal
#endif

end function colindex_local_to_global_procindex


!=========================================================================
function colindex_local_to_global_descriptor(desc,jlocal) result(jglobal)
 implicit none
 integer,intent(in)          :: desc(NDEL),jlocal
 integer                     :: jglobal
!=====
#ifdef HAVE_SCALAPACK
 integer          :: iprow,ipcol,nprow,npcol
#endif
!=====

#ifdef HAVE_SCALAPACK
 call BLACS_GRIDINFO(desc(CTXT_),nprow,npcol,iprow,ipcol)
 jglobal = INDXL2G(jlocal,desc(NB_),ipcol,desc(CSRC_),npcol)
#else
 jglobal = jlocal
#endif

end function colindex_local_to_global_descriptor


!=========================================================================
function colindex_global_to_local_descriptor(desc,jglobal) result(jlocal)
 implicit none
 integer,intent(in)          :: desc(NDEL),jglobal
 integer                     :: jlocal
!=====
#ifdef HAVE_SCALAPACK
 integer          :: iprow,ipcol,nprow,npcol
#endif
!=====

#ifdef HAVE_SCALAPACK
 call BLACS_GRIDINFO(desc(CTXT_),nprow,npcol,iprow,ipcol)
 if( ipcol == INDXG2P(jglobal,desc(NB_),0,desc(CSRC_),npcol) ) then
   jlocal = INDXG2L(jglobal,desc(NB_),0,desc(CSRC_),npcol)
 else
   jlocal = 0
 endif
#else
 jlocal = jglobal
#endif

end function colindex_global_to_local_descriptor


!=========================================================================
subroutine set_auxil_block_size(block_size_max)
 implicit none
 integer,intent(in) :: block_size_max
!=====
!=====

 if( block_size_max < 1 ) then
   MBLOCK_AUXIL = 1

 else
   MBLOCK_AUXIL = 2**( FLOOR( LOG(REAL(block_size_max,dp)) / LOG( 2.0_dp ) ) )
   MBLOCK_AUXIL = MIN(MBLOCK_AUXIL,block_row)

 endif

 write(stdout,'(/1x,a,i4)') 'SCALAPACK block size for auxiliary basis: ',MBLOCK_AUXIL

end subroutine set_auxil_block_size


!=========================================================================
subroutine finish_scalapack()
 implicit none
!=====

#ifdef HAVE_SCALAPACK
 call BLACS_GRIDEXIT( cntxt_sd )
 call BLACS_GRIDEXIT( cntxt_cd )
 call BLACS_GRIDEXIT( cntxt_rd )
 call BLACS_GRIDEXIT( cntxt_3center )
 if( parallel_ham ) then
   if( cntxt_ham > 0 ) call BLACS_GRIDEXIT( cntxt_ham )
 endif
 if( cntxt_auxil > 0 ) call BLACS_GRIDEXIT( cntxt_auxil )
 call BLACS_EXIT( 0 )
#endif

end subroutine finish_scalapack

#ifdef HAVE_SCALAPACK
!=========================================================================
subroutine diagonalize_davidson_sca(tolerance,desc_ham,ham,neig,eigval,desc_vec,eigvec)
 implicit none

 real(dp),intent(in)  :: tolerance
 real(dp),intent(in)  :: ham(:,:)
 integer,intent(in)   :: desc_ham(NDEL),desc_vec(NDEL)
 integer,intent(in)   :: neig
 real(dp),intent(out) :: eigval(:)
 real(dp),intent(out) :: eigvec(:,:)
!=====
 integer              :: ncycle
 integer              :: mmat,imat
 integer              :: mm,mm_max
 integer              :: ieig,icycle
 real(dp),allocatable :: bb(:,:),atilde(:,:),ab(:,:),qq(:,:)
 real(dp),allocatable :: lambda(:),alphavec(:,:)
 real(dp)             :: residual_norm,norm2_i
 integer              :: desc_bb(NDEL)
 integer              :: mbb,nbb
 integer              :: desc_qq(NDEL)
 integer              :: mqq,nqq
 integer              :: desc_at(NDEL)
 integer              :: mat,nat
 integer              :: info
 integer              :: ilocal,jlocal,iglobal,jglobal
 real(dp),allocatable :: ham_diag(:)
!=====

#if defined(HAVE_SCALAPACK)
 write(stdout,'(/,1x,a,i5)') 'Davidson diago for eigenvector count: ',neig

 eigval(:) = 0.0_dp

 mmat = desc_ham(M_)
 allocate(ham_diag(mmat))

 !
 ! Broadcast the diagonal of H to all procs
 ham_diag(:) = 0.0_dp
 do jlocal=1,SIZE(ham,DIM=2)
   jglobal = colindex_local_to_global(desc_ham,jlocal)
   do ilocal=1,SIZE(ham,DIM=1)
     iglobal = rowindex_local_to_global(desc_ham,ilocal)

     if( iglobal == jglobal ) ham_diag(iglobal) = ham(ilocal,jlocal)
   enddo
 enddo
 call xsum_world(ham_diag)

 ncycle = 30
 mm     = neig
 mm_max = mm * ncycle
 if( mm_max > mmat ) then
   ncycle = mmat / neig
   mm_max = mm * ncycle
 endif

 mbb = NUMROC(mmat  ,block_row,iprow_sd,first_row,nprow_sd)
 nbb = NUMROC(mm_max,block_col,ipcol_sd,first_col,npcol_sd)
 call DESCINIT(desc_bb,mmat,mm_max,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,mbb),info)

 allocate(bb(mbb,nbb))
 allocate(ab(mbb,nbb))

 mqq = NUMROC(mmat,block_row,iprow_sd,first_row,nprow_sd)
 nqq = NUMROC(neig,block_col,ipcol_sd,first_col,npcol_sd)
 call DESCINIT(desc_qq,mmat,neig,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,mqq),info)
 allocate(qq(mqq,nqq))



 !
 ! Initialize with stupid coefficients
! bb(:,1:neig) = 0.01_dp
! forall(ieig=1:neig)
!   bb(ieig,ieig) = 1.0_dp
! end forall
 bb(:,:) = 0.01_dp
 do jlocal=1,nbb
   jglobal = colindex_local_to_global(desc_bb,jlocal)
   do ilocal=1,mbb
     iglobal = rowindex_local_to_global(desc_bb,ilocal)
     bb(ilocal,jlocal) = MIN( EXP( -REAL(iglobal,dp) ) , 0.1_dp )
     if( iglobal == jglobal ) bb(ilocal,jlocal) = 1.0_dp
   enddo
 enddo
 call orthogonalize_sca(desc_bb,1,neig,bb)


! ab(:,1:neig) = MATMUL( ham(:,:) , bb(:,1:neig) )
 call PDGEMM('N','N',mmat,neig,mmat,1.0_dp,ham,1,1,desc_ham,bb,1,1,desc_bb,0.0_dp,ab,1,1,desc_bb)


 do icycle=1,ncycle

   mm = icycle * neig
   allocate(lambda(mm))

   mat = NUMROC(mm,block_row,iprow_sd,first_row,nprow_sd)
   nat = NUMROC(mm,block_col,ipcol_sd,first_col,npcol_sd)
   call DESCINIT(desc_at,mm,mm,block_row,block_col,first_row,first_col,cntxt_sd,MAX(1,mat),info)
   allocate(atilde(mat,nat),alphavec(mat,nat))

   !atilde(1:mm,1:mm) = MATMUL( TRANSPOSE(bb(:,1:mm)) , ab(:,1:mm) )
   call PDGEMM('T','N',mm,mm,mmat,1.0_dp,bb,1,1,desc_bb,ab,1,1,desc_bb,0.0_dp,atilde,1,1,desc_at)


   call diagonalize_sca(mm,desc_at,atilde,lambda,desc_at,alphavec)

   deallocate(atilde)

   !write(stdout,*) 'icycle',icycle,lambda(1:mm)

   ! qq = bb * alphavec
   call PDGEMM('N','N',mmat,neig,mm,1.0_dp,bb,1,1,desc_bb,alphavec,1,1,desc_at,0.0_dp,qq,1,1,desc_qq)
   eigvec(:,:) = qq(:,:)
   eigval(1:neig) = lambda(1:neig)
   ! qq = qq * Lambda
   do ieig=1,neig
     call PDSCAL(mmat,lambda(ieig),qq,1,ieig,desc_qq,1)
   enddo


   ! qq = ab * alphavec - lambda * bb * alphavec
   call PDGEMM('N','N',mmat,neig,mm,1.0_dp,ab,1,1,desc_bb,alphavec,1,1,desc_at,-1.0_dp,qq,1,1,desc_qq)


   deallocate(alphavec)
!   residual_norm = 0.0_dp
!   do ieig=1,neig
!
!     qq(:,ieig) = MATMUL( ab(:,1:mm) ,  alphavec(1:mm,ieig) ) &
!                   - lambda(ieig) * MATMUL( bb(:,1:mm) , alphavec(1:mm,ieig) )
!
!     residual_norm = MAX( residual_norm , NORM2(qq(:,ieig)) )
!   enddo

   residual_norm = 0.0_dp
   do ieig=1,neig
     norm2_i = 0.0_dp
     call PDNRM2(mmat,norm2_i,qq,1,ieig,desc_qq,1)
     residual_norm = MAX( residual_norm , norm2_i )
   enddo
   call xmax_world(residual_norm)


   write(stdout,'(1x,a,i4,1x,i4,1x,es12.4,1x,f18.8)') 'Cycle, Subspace dim, Max residual norm, Electronic energy: ', &
                                                      icycle,mm,residual_norm,lambda(1)

   !
   ! Convergence reached... or not
   if( icycle == ncycle .OR. residual_norm < tolerance ) then
     exit
   endif


   !
   ! New trial vectors
   !
!   forall(imat=1:nmat,ieig=1:neig)
!     bb(imat,mm+ieig) = qq(imat,ieig) / ( lambda(ieig) - ham(imat,imat) )
!   end forall
   do jlocal=1,nqq
     jglobal = colindex_local_to_global(desc_qq,jlocal)
     do ilocal=1,mqq
       iglobal = rowindex_local_to_global(desc_qq,ilocal)
       qq(ilocal,jlocal) = qq(ilocal,jlocal) / ( lambda(jglobal) - ham_diag(iglobal) )
     enddo
   enddo


   call PDGEMR2D(mmat,neig,qq,1,1,desc_qq,bb,1,mm+1,desc_bb,cntxt_sd)

   call orthogonalize_sca(desc_bb,mm+1,mm+neig,bb)

   !ab(:,mm+1:mm+neig) = MATMUL( ham(:,:) , bb(:,mm+1:mm+neig) )
   call PDGEMM('N','N',mmat,neig,mmat,1.0_dp,ham,1,1,desc_ham,bb,1,mm+1,desc_bb,0.0_dp,ab,1,mm+1,desc_bb)



   deallocate(lambda)


 enddo ! icycle

 if( ALLOCATED(lambda) ) deallocate(lambda)
 deallocate(ab,qq,bb,ham_diag)

#else
 call die('diagonalize_davidson_sca: should not be called when not compiled with HAVE_SCALAPACK')
#endif

end subroutine diagonalize_davidson_sca
#endif

#ifdef HAVE_SCALAPACK
!=========================================================================
subroutine orthogonalize_sca(desc_vec,mvec_ortho,nvec_ortho,vec)
 implicit none
!=====
 integer,intent(in)     :: desc_vec(NDEL)
 integer,intent(in)     :: mvec_ortho,nvec_ortho
 real(dp),intent(inout) :: vec(:,:)
!=====
 integer :: ii,ivec,jvec,mglobal
 real(dp) :: dot_prod_ij
 real(dp) :: norm_i
!=====

#if defined(HAVE_SCALAPACK)
 mglobal = desc_vec(M_)

 do ivec=mvec_ortho,nvec_ortho
   ! Orthogonalize to previous vectors
   do jvec=1,ivec-1
     call PDDOT(mglobal,dot_prod_ij,vec,1,ivec,desc_vec,1,vec,1,jvec,desc_vec,1)

!     vec(:,ivec) = vec(:,ivec) - vec(:,jvec) * DOT_PRODUCT( vec(:,ivec) , vec(:,jvec) )
     call PDAXPY(mglobal,-dot_prod_ij,vec,1,jvec,desc_vec,1,vec,1,ivec,desc_vec,1)

   enddo

   ! Normalize
   ! vec(:,ivec) = vec(:,ivec) / NORM2( vec(:,ivec) )
   call PDNRM2(mglobal,norm_i,vec,1,ivec,desc_vec,1)
   call PDSCAL(mglobal,1.0_dp/norm_i,vec,1,ivec,desc_vec,1)

 enddo

#else
 call die('diagonalize_davidson_sca: should not be called when not compiled with HAVE_SCALAPACK')
#endif

end subroutine orthogonalize_sca
#endif

!=========================================================================
subroutine select_nprow_npcol(scalapack_block_min,mmat,nmat,nprow,npcol)
 implicit none

 integer,intent(in)  :: scalapack_block_min,mmat,nmat
 integer,intent(out) :: nprow,npcol
!=====
 integer :: max_dim
!=====

 max_dim = MAX(mmat,nmat)

 ! If mmat ~ nmat (within some tolerance, then use a square distribution!
 if( ABS( REAL(mmat-nmat,dp) ) / REAL(max_dim,dp) < 0.50_dp ) then

   nprow = MIN( FLOOR(SQRT(REAL(nproc_sca,dp))) , max_dim / scalapack_block_min )
   nprow = MAX(nprow,1)
   npcol = MIN( FLOOR(SQRT(REAL(nproc_sca,dp))) , max_dim / scalapack_block_min )
   npcol = MAX(npcol,1)

 else if( nmat < mmat ) then

   npcol = nmat / scalapack_block_min
   npcol = MAX(npcol,1)
   nprow = MIN( nproc_sca / npcol , mmat / scalapack_block_min )
   nprow = MAX(nprow,1)

 else

   nprow = mmat / scalapack_block_min
   nprow = MAX(nprow,1)
   npcol = MIN( nproc_sca / nprow , nmat / scalapack_block_min )
   npcol = MAX(npcol,1)

 endif


 if( nprow * npcol > nproc_sca ) call die('select_nprow_npcol: forbidden SCALAPACK grid')

end subroutine select_nprow_npcol


!=========================================================================
end module m_scalapack
!=========================================================================
