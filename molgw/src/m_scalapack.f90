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
subroutine create_distributed_copy(matrix_global,desc,matrix)
 implicit none
 integer,intent(in)   :: desc(ndel)
 real(dp),intent(in)  :: matrix_global(:,:)
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
subroutine gather_distributed_copy(desc,matrix,matrix_global)
 implicit none
 integer,intent(in)   :: desc(ndel)
 real(dp),intent(in)  :: matrix(:,:)
 real(dp),intent(out) :: matrix_global(:,:)
!=====
 integer              :: contxt
 integer              :: mlocal,nlocal,mglobal,nglobal
 integer              :: ilocal,jlocal,iglobal,jglobal
 integer              :: rank_master
!=====

 contxt = desc(2)

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
 call DGSUM2D(contxt,'A',' ',mglobal,nglobal,matrix_global,nglobal,0,0)


end subroutine gather_distributed_copy


!=========================================================================
! Multiply on the left or on the right a distributed square matrix 
! with a non-distributed diagonal
!=========================================================================
subroutine matmul_diag_sca(side,diag,desc,matrix)
 implicit none
 character(len=1),intent(in) :: side
 integer,intent(in)          :: desc(ndel)
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
 real(dp),allocatable :: matrix(:,:)
 integer :: rank_master
 integer,external :: NUMROC
!=====

#ifdef HAVE_SCALAPACK
 nprow = MIN(nprow_sd,nmat/scalapack_block_min)
 npcol = MIN(npcol_sd,nmat/scalapack_block_min)
 nprow = MAX(nprow,1)
 npcol = MAX(npcol,1)

 if( nprow /= 1 .OR. npcol /= 1 ) then

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
     call create_distributed_copy(matrix_global,descm,matrix)

     call diagonalize_sca(nmat,descm,matrix,eigval)

     call gather_distributed_copy(descm,matrix,matrix_global)
     deallocate(matrix)

     call BLACS_GRIDEXIT( cntxt )

   endif

   ! Then the master proc (0,0) broadcasts to all the others
   call xbcast(rank_master,matrix_global)
   call xbcast(rank_master,eigval)


 else ! Only one SCALAPACK proc

   call diagonalize(nmat,matrix_global,eigval)

 endif

#else

 call diagonalize(nmat,matrix_global,eigval)

#endif


end subroutine diagonalize_scalapack


!=========================================================================
! Calculate D = A * B * C for non-distributed matrices
!
!=========================================================================
subroutine product_abc_scalapack(scalapack_block_min,a_matrix,b_matrix,c_matrix,d_matrix)
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
 integer :: desca(ndel),descb(ndel),descc(ndel),descd(ndel)
 integer :: descm(ndel)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
 integer :: rank_master
 integer,external :: NUMROC
!=====

 mmat  = SIZE( a_matrix , DIM=1)
 kmat1 = SIZE( a_matrix , DIM=2)
 kmat  = SIZE( b_matrix , DIM=1)
 lmat1 = SIZE( b_matrix , DIM=2)
 lmat  = SIZE( c_matrix , DIM=1)
 nmat1 = SIZE( c_matrix , DIM=2)
 mmat1 = SIZE( d_matrix , DIM=1)
 nmat  = SIZE( d_matrix , DIM=2)

 if( mmat1 /= mmat ) call die('Dimension error in product_abc_scalapack')
 if( nmat1 /= nmat ) call die('Dimension error in product_abc_scalapack')
 if( kmat1 /= kmat ) call die('Dimension error in product_abc_scalapack')
 if( lmat1 /= lmat ) call die('Dimension error in product_abc_scalapack')


#ifdef HAVE_SCALAPACK
 nprow = MIN(nprow_sd,mmat/scalapack_block_min)
 npcol = MIN(npcol_sd,mmat/scalapack_block_min)
 nprow = MAX(nprow,1)
 npcol = MAX(npcol,1)

 if( nprow /= 1 .OR. npcol /= 1 ) then

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
   write(stdout,'(a,i4,a,i4)') ' Matrix product using SCALAPACK with a grid',nprow,' x ',npcol
  
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

  
     call gather_distributed_copy(descd,d_matrix_local,d_matrix)
     deallocate(d_matrix_local)
  
     call BLACS_GRIDEXIT( cntxt )
  
   endif
  
   ! Then the master proc (0,0) broadcasts to all the others
   call xbcast(rank_master,d_matrix)


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


end subroutine product_abc_scalapack


!=========================================================================
! Calculate C = A^T * B * A for non-distributed matrices
!
!=========================================================================
subroutine product_transaba_scalapack(scalapack_block_min,a_matrix,b_matrix,c_matrix)
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
 integer :: desca(ndel),descb(ndel),descc(ndel)
 integer :: descm(ndel)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
 integer :: rank_master
 integer,external :: NUMROC
!=====

 kmat  = SIZE( a_matrix , DIM=1)
 mmat  = SIZE( a_matrix , DIM=2)
 kmat1 = SIZE( b_matrix , DIM=1)
 kmat2 = SIZE( b_matrix , DIM=2)
 mmat1 = SIZE( c_matrix , DIM=1)
 mmat2 = SIZE( c_matrix , DIM=2)

 if( mmat1 /= mmat ) call die('Dimension error in product_transaba_scalapack')
 if( mmat2 /= mmat ) call die('Dimension error in product_transaba_scalapack')
 if( kmat1 /= kmat ) call die('Dimension error in product_transaba_scalapack')
 if( kmat2 /= kmat ) call die('Dimension error in product_transaba_scalapack')


#ifdef HAVE_SCALAPACK
 nprow = MIN(nprow_sd,mmat/scalapack_block_min)
 npcol = MIN(npcol_sd,mmat/scalapack_block_min)
 nprow = MAX(nprow,1)
 npcol = MAX(npcol,1)

 if( nprow /= 1 .OR. npcol /= 1 ) then

   call BLACS_GET( -1, 0, cntxt )
   call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
   call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )
   write(stdout,'(a,i4,a,i4)') ' Matrix product using SCALAPACK with a grid',nprow,' x ',npcol
  
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

     call gather_distributed_copy(descc,c_matrix_local,c_matrix)
     deallocate(c_matrix_local)
  

     call BLACS_GRIDEXIT( cntxt )
  
   endif
  
   ! Then the master proc (0,0) broadcasts to all the others
   call xbcast(rank_master,c_matrix)


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


end subroutine product_transaba_scalapack


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
