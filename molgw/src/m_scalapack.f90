!=========================================================================
! This file is part of MOLGW.
!
! This file contains a module with some SCALAPACK wrappers
! These wrappers are meant to be independant from the global variables 
! defined in module m_mpi:
! All the required SCALAPACK proc grid information, block size are contained in a descriptor
!
! All the subroutine are meant to be valid with OR without SCALAPACK

!=========================================================================
module m_scalapack
 use m_definitions
 use m_warning,only: die,issue_warning
 use m_tools,only: diagonalize
 use m_mpi
#ifdef HAVE_MPI
 use mpi
#endif


 interface diagonalize_sca
   module procedure diagonalize_inplace_sca
   module procedure diagonalize_outofplace_sca
 end interface


contains


!=========================================================================
! Create a distributed copy of a global matrix owned by each core
!
!=========================================================================
subroutine create_distributed_copy(nglobal,matrix_global,desc,matrix)
 implicit none
 integer,intent(in)   :: nglobal
 integer,intent(in)   :: desc(ndel)
 real(dp),intent(in)  :: matrix_global(nglobal,nglobal)
 real(dp),intent(out) :: matrix(:,:)
!=====
 integer              :: mlocal,nlocal
 integer              :: ilocal,jlocal,iglobal,jglobal
!=====

 mlocal = SIZE( matrix , DIM=1 )
 nlocal = SIZE( matrix , DIM=2 )

 do jlocal=1,nlocal
   jglobal = colindex_local_to_global(desc,jlocal)
   do ilocal=1,mlocal
     iglobal = rowindex_local_to_global(desc,ilocal)
     matrix(ilocal,jlocal) = matrix_global(iglobal,jglobal)
   enddo
 enddo


end subroutine create_distributed_copy


!=========================================================================
! Gather a distributed matrix into a global matrix owned by each core
!
!=========================================================================
subroutine gather_distributed_copy(nglobal,desc,matrix,matrix_global)
 implicit none
 integer,intent(in)   :: nglobal
 integer,intent(in)   :: desc(ndel)
 real(dp),intent(in)  :: matrix(:,:)
 real(dp),intent(out) :: matrix_global(nglobal,nglobal)
!=====
 integer              :: contxt
 integer              :: mlocal,nlocal
 integer              :: ilocal,jlocal,iglobal,jglobal
 integer              :: rank_master
!=====

 contxt = desc(2)

 mlocal = SIZE( matrix , DIM=1 )
 nlocal = SIZE( matrix , DIM=2 )

 matrix_global(:,:) = 0.0_dp
 do jlocal=1,nlocal
   jglobal = colindex_local_to_global(desc,jlocal)
   do ilocal=1,mlocal
     iglobal = rowindex_local_to_global(desc,ilocal)
     matrix_global(iglobal,jglobal) = matrix(ilocal,jlocal)
   enddo
 enddo

 ! Only the master proc (0,0) gets the complete information
 call DGSUM2D(contxt,'A',' ',nglobal,nglobal,matrix_global,nglobal,0,0)


end subroutine gather_distributed_copy


!=========================================================================
! Multiply on the left or on the right a distributed square matrix 
! with a non-distributed diagonal
!=========================================================================
subroutine matmul_diag_sca(side,nglobal,diag,desc,matrix)
 implicit none
 character(len=1),intent(in) :: side
 integer,intent(in)          :: nglobal
 integer,intent(in)          :: desc(ndel)
 real(dp),intent(in)         :: diag(nglobal)
 real(dp),intent(inout)      :: matrix(:,:)
!=====
 integer                     :: mlocal,nlocal
 integer                     :: iglobal,jglobal,ilocal,jlocal
!=====


#ifndef HAVE_SCALAPACK

 select case(side)
 case('L')
   forall(iglobal=1:nglobal)
     matrix(iglobal,:) = matrix(iglobal,:) * diag(iglobal)
   end forall

 case('R')
   forall(jglobal=1:nglobal)
     matrix(:,jglobal) = matrix(:,jglobal) * diag(jglobal)
   end forall

 case default
   if( side /= 'L' .AND. side /= 'R' ) then
     call die('matmul_diag_sca: argument side should be L or R')
   endif
 end select

#else

 mlocal = SIZE( matrix , DIM=1 )
 nlocal = SIZE( matrix , DIM=2 )

 select case(side)
 case('L')
   do ilocal=1,mlocal
     iglobal = rowindex_local_to_global(desc,ilocal)
     matrix(ilocal,:) = matrix(ilocal,:) * diag(iglobal)
   enddo

 case('R')
! Alternative coding
!   do jglobal=1,nglobal
!     call PDSCAL(nglobal,diag(jglobal),matrix,1,jglobal,desc,1)
!   enddo
   do jlocal=1,nlocal
     jglobal = colindex_local_to_global(desc,jlocal)
     matrix(:,jlocal) = matrix(:,jlocal) * diag(jglobal)
   enddo

 case default
   if( side /= 'L' .AND. side /= 'R' ) then
     call die('matmul_diag_sca: argument side should be L or R')
   endif
 end select

#endif


end subroutine matmul_diag_sca


!=========================================================================
! Diagonalize a distributed matrix
!=========================================================================
subroutine diagonalize_inplace_sca(nglobal,desc,matrix,eigval)
 implicit none
 integer,intent(in)     :: desc(ndel),nglobal
 real(dp),intent(inout) :: matrix(:,:)
 real(dp),intent(out)   :: eigval(nglobal)
!=====
 integer              :: desc_eigvec(ndel)
 integer              :: mlocal,nlocal
 integer              :: lwork,info
 real(dp),allocatable :: work(:)
 real(dp),allocatable :: eigvec(:,:)
 integer              :: neigval,neigvec
 integer,allocatable  :: iwork(:)
 integer              :: liwork
#ifdef SELECT_PDSYEVX
 real(dp)             :: ABSTOL
 integer              :: iclustr(2*nprow_sd*npcol_sd)
 real(dp)             :: gap(nprow_sd*npcol_sd)
 integer              :: ifail(nglobal)
 real(dp),external    :: PDLAMCH
#endif
!=====

#ifdef HAVE_SCALAPACK
 desc_eigvec = desc

 mlocal = SIZE( matrix , DIM=1 )
 nlocal = SIZE( matrix , DIM=2 )

 allocate(eigvec(mlocal,nlocal))

 !
 ! First call to get the dimension of the array work
 lwork = -1
 liwork = -1
 allocate(work(1))
 allocate(iwork(1))
#ifdef SELECT_PDSYEVX
 ABSTOL = PDLAMCH(desc(2), 'U')
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
#else
 call PDSYEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)
! call PDSYEVR('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
!              neigval,neigvec,eigval,                                &
!              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
!              info)
#endif


 !
 ! Second call to actually perform the diago
 lwork = NINT(work(1))
 liwork = iwork(1)

 deallocate(work)
 deallocate(iwork)
 allocate(work(lwork))
 allocate(iwork(liwork))
#ifdef SELECT_PDSYEVX
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
#else
 call PDSYEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)
! call PDSYEVR('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
!              neigval,neigvec,eigval,                                &
!              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
!              info)
#endif
 deallocate(work)
 deallocate(iwork)


 matrix(:,:) = eigvec(:,:)

 deallocate(eigvec)

#else

 call diagonalize(nglobal,matrix,eigval)

#endif


end subroutine diagonalize_inplace_sca


!=========================================================================
! Diagonalize a distributed matrix
!=========================================================================
subroutine diagonalize_outofplace_sca(nglobal,desc,matrix,eigval,desc_eigvec,eigvec)
 implicit none
 integer,intent(in)     :: nglobal
 integer,intent(in)     :: desc(ndel)
 integer,intent(in)     :: desc_eigvec(ndel)
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
#ifdef SELECT_PDSYEVX
 ABSTOL = PDLAMCH(desc(2), 'U')
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
#else
 call PDSYEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)
! call PDSYEVR('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
!              neigval,neigvec,eigval,                                &
!              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
!              info)
#endif


 !
 ! Second call to actually perform the diago
 lwork = NINT(work(1))
 liwork = iwork(1)

 deallocate(work)
 deallocate(iwork)
 allocate(work(lwork))
 allocate(iwork(liwork))
#ifdef SELECT_PDSYEVX
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
#else
 call PDSYEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)
! call PDSYEVR('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
!              neigval,neigvec,eigval,                                &
!              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
!              info)
#endif
 deallocate(work)
 deallocate(iwork)


#else

 call diagonalize(nglobal,matrix,eigval,eigvec)

#endif


end subroutine diagonalize_outofplace_sca


!=========================================================================
! Diagonalize a non-distributed matrix
!
!=========================================================================
subroutine diagonalize_scalapack(scalapack_block_min,nmat,matrix_global,eigval)
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
 integer :: descm(ndel),descz(ndel)
 real(dp) :: alpha
 real(dp),allocatable :: matrix(:,:)
 real(dp),allocatable :: work(:)
 integer :: lwork
 integer :: rank_master
 integer :: rank_sca,nprocs_sca
 integer,external :: NUMROC,INDXL2G
!=====

#ifdef HAVE_SCALAPACK
 nprow = MIN(nprow_sd,nmat/scalapack_block_min)
 npcol = MIN(npcol_sd,nmat/scalapack_block_min)
 nprow = MAX(nprow,1)
 npcol = MAX(npcol,1)

 call BLACS_PINFO( rank_sca, nprocs_sca )

 call BLACS_GET( -1, 0, cntxt )
 call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
 write(stdout,'(a,i4,a,i4)') ' Diagonalization using SCALAPACK with a grid',nprow,' x ',npcol

 ! Find the master
 if( iprow == 0 .AND. ipcol == 0 ) then
   rank_master = rank
 else
   rank_master = -1
 endif
 call xmax(rank_master)

 !
 ! Participate to the diagonalization only if the CPU has been selected 
 ! in the grid
 if( cntxt > 0 ) then
   mlocal = NUMROC(nmat,block_row,iprow,first_row,nprow)
   nlocal = NUMROC(nmat,block_col,ipcol,first_col,npcol)

   allocate(matrix(mlocal,nlocal))
    
   call DESCINIT(descm,nmat,nmat,block_row,block_col,first_row,first_col,cntxt,MAX(1,mlocal),info)

   ! Set up the local copy of the matrix_global
   call create_distributed_copy(nmat,matrix_global,descm,matrix)

   call diagonalize_sca(nmat,descm,matrix,eigval)

   call gather_distributed_copy(nmat,descm,matrix,matrix_global)
   deallocate(matrix)

   call BLACS_GRIDEXIT( cntxt )

 endif

 ! Then the master proc (0,0) broadcasts to all the others
 call xbcast(rank_master,matrix_global)
 call xbcast(rank_master,eigval)


#else

 call diagonalize(nmat,matrix_global,eigval)

#endif


end subroutine diagonalize_scalapack


!=========================================================================
! Transform a lower or upper triangular matrix into the full representation
!
!=========================================================================
subroutine symmetrize_matrix_sca(uplo,nglobal,desc,matrix,desc_tmp,matrix_tmp)
 implicit none
 character(len=1)       :: uplo
 integer,intent(in)     :: nglobal
 integer,intent(in)     :: desc(ndel),desc_tmp(ndel)
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
end module m_scalapack
!=========================================================================
