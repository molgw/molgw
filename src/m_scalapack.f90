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
 use m_tools,only: diagonalize
 use m_mpi
#ifdef HAVE_MPI
 use mpi
#endif

 ! Indexes in the BLACS descriptor
 integer,parameter :: DTYPE_A = 1
 integer,parameter :: CTXT_A  = 2
 integer,parameter :: M_A     = 3
 integer,parameter :: N_A     = 4
 integer,parameter :: MB_A    = 5
 integer,parameter :: NB_A    = 6
 integer,parameter :: RSRC_A  = 7
 integer,parameter :: CSRC_A  = 8
 integer,parameter :: LLD_A   = 9

 !
 ! SCALAPACK variables
 !
 integer,parameter :: NDEL=9
 integer,parameter :: block_col = 32
 integer,parameter :: block_row = 32
 integer,parameter :: first_row = 0
 integer,parameter :: first_col = 0
 
 ! Specific values for the MPI / SCALAPACK transposition
 integer,parameter :: MBLOCK_AUXIL = 1
 integer,parameter :: NBLOCK_AUXIL = 1

 integer,protected :: nproc_sca = 1
 integer,protected :: iproc_sca = 0


 ! SCALAPACK grid: auxil distribution (so to mimic MPI distribution on auxiliary functions)
 integer,protected :: cntxt_auxil
 integer,protected :: nprow_auxil,npcol_auxil,iprow_auxil,ipcol_auxil

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

 interface diagonalize_sca
   module procedure diagonalize_inplace_sca
   module procedure diagonalize_outofplace_sca
 end interface

 interface create_distributed_copy
   module procedure create_distributed_copy_nospin
   module procedure create_distributed_copy_spin
 end interface create_distributed_copy

 interface gather_distributed_copy
   module procedure gather_distributed_copy_nospin
   module procedure gather_distributed_copy_spin
 end interface gather_distributed_copy

#ifdef HAVE_SCALAPACK
 integer,external :: NUMROC,INDXL2G,INDXG2L,INDXG2P,PDLATRA
#endif

contains


!=========================================================================
! Create a distributed copy of a global matrix owned by each core
!
!=========================================================================
subroutine create_distributed_copy_nospin(matrix_global,desc,matrix)
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

 do jlocal=1,nlocal
   jglobal = colindex_local_to_global(desc,jlocal)
   do ilocal=1,mlocal
     iglobal = rowindex_local_to_global(desc,ilocal)
     matrix(ilocal,jlocal) = matrix_global(iglobal,jglobal)
   enddo
 enddo


end subroutine create_distributed_copy_nospin


!=========================================================================
! Create a distributed copy of a global matrix owned by each core
! with spin
!=========================================================================
subroutine create_distributed_copy_spin(matrix_global,desc,matrix)
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

 do idim3=1,ndim3
   do jlocal=1,nlocal
     jglobal = colindex_local_to_global(desc,jlocal)
     do ilocal=1,mlocal
       iglobal = rowindex_local_to_global(desc,ilocal)
       matrix(ilocal,jlocal,idim3) = matrix_global(iglobal,jglobal,idim3)
     enddo
   enddo
 enddo


end subroutine create_distributed_copy_spin


!=========================================================================
! Gather a distributed matrix into a global matrix owned by each core
!
!=========================================================================
subroutine gather_distributed_copy_nospin(desc,matrix,matrix_global)
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

 cntxt = desc(CTXT_A)
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

#endif

end subroutine gather_distributed_copy_nospin


!=========================================================================
! Gather a distributed matrix into a global matrix owned by each core
!
!=========================================================================
subroutine gather_distributed_copy_spin(desc,matrix,matrix_global)
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

 cntxt = desc(CTXT_A)
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

#endif

end subroutine gather_distributed_copy_spin



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
! Diagonalize a distributed matrix
!=========================================================================
subroutine diagonalize_inplace_sca(nglobal,desc,matrix,eigval)
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
#ifdef SELECT_PDSYEVX
 liwork = -1
 allocate(iwork(1))
 allocate(work(3))
 ABSTOL = PDLAMCH(desc(CTXT_A), 'U')
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
#else
 allocate(work(1))
 call PDSYEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)
#endif


 !
 ! Second call to actually perform the diago
 lwork = NINT(work(1))

 deallocate(work)
 allocate(work(lwork))
#ifdef SELECT_PDSYEVX
 deallocate(iwork)
 liwork = iwork(1)
 allocate(iwork(liwork))
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
 deallocate(iwork)
#else
 call PDSYEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)
#endif
 deallocate(work)


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
 real(dp),external    :: PDLAMCH
#endif
!=====

#ifdef HAVE_SCALAPACK

 !
 ! First call to get the dimension of the array work
 lwork = -1
#ifdef SELECT_PDSYEVX
 allocate(work(3))
 liwork = -1
 allocate(iwork(1))
 ABSTOL = PDLAMCH(desc(CTXT_A), 'U')
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
#else
 allocate(work(1))
 call PDSYEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)
#endif


 !
 ! Second call to actually perform the diago
 lwork = NINT(work(1))

 deallocate(work)
 allocate(work(lwork))
#ifdef SELECT_PDSYEVX
 liwork = iwork(1)
 deallocate(iwork)
 allocate(iwork(liwork))
 call PDSYEVX('V','A','L',nglobal,matrix,1,1,desc,0.0_dp,0.0_dp,0,0, &
              ABSTOL,neigval,neigvec,eigval,0.0_dp,                  &
              eigvec,1,1,desc_eigvec,work,lwork,iwork,liwork,        &
              ifail,iclustr,gap,info)
 deallocate(iwork)
#else
 call PDSYEV('V','L',nglobal,matrix,1,1,desc,eigval,eigvec,1,1,desc_eigvec,work,lwork,info)
#endif
 deallocate(work)


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
 integer :: descm(NDEL),descz(NDEL)
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
 integer :: desca(NDEL),descb(NDEL),descc(NDEL),descd(NDEL)
 integer :: descm(NDEL)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
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
 integer :: desca(NDEL),descb(NDEL),descc(NDEL)
 integer :: descm(NDEL)
 integer :: nprow,npcol,iprow,ipcol
 integer :: info
 integer,external :: NUMROC
!=====

 kmat  = SIZE( a_matrix , DIM=1)
 mmat  = SIZE( a_matrix , DIM=2)
 kmat1 = SIZE( b_matrix , DIM=1)
 kmat2 = SIZE( b_matrix , DIM=2)
 mmat1 = SIZE( c_matrix , DIM=1)
 mmat2 = SIZE( c_matrix , DIM=2)

 if( mmat1 /= mmat ) then
   write(msg,*) 'mmat1 /= mmat',mmat1,mmat
   call die('Dimension error in product_transaba_scalapack'//msg)
 endif
 if( mmat2 /= mmat ) then
   write(msg,*) 'mmat2 /= mmat',mmat2,mmat
   call die('Dimension error in product_transaba_scalapack'//msg)
 endif
 if( kmat1 /= kmat ) then
   write(msg,*) 'kmat1 /= kmat',kmat1,kmat
   call die('Dimension error in product_transaba_scalapack'//msg)
 endif
 if( kmat2 /= kmat ) then
   write(msg,*) 'kmat2 /= kmat',kmat2,kmat
   call die('Dimension error in product_transaba_scalapack'//msg)
 endif


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


end subroutine product_transaba_scalapack


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
 nprow = MIN(nprow_sd,kmat1/scalapack_block_min)
 npcol = MIN(npcol_sd,kmat2/scalapack_block_min)
 nprow = MAX(nprow,1)
 npcol = MAX(npcol,1)

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
subroutine product_abc_sca(desca,a_matrix_local,descb,b_matrix_local,descc,c_matrix_local,descd,d_matrix_local)
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
 integer,external :: NUMROC
!=====

 mmat  = desca(M_A)
 kmat1 = desca(N_A)
 kmat  = descb(M_A)
 lmat1 = descb(N_A)
 lmat  = descc(M_A)
 nmat1 = descc(N_A)
 mmat1 = descd(M_A)
 nmat  = descd(N_A)

 if( mmat1 /= mmat ) call die('Dimension error in product_abc_scalapack')
 if( nmat1 /= nmat ) call die('Dimension error in product_abc_scalapack')
 if( kmat1 /= kmat ) call die('Dimension error in product_abc_scalapack')
 if( lmat1 /= lmat ) call die('Dimension error in product_abc_scalapack')


#ifdef HAVE_SCALAPACK

 cntxt = desca(CTXT_A)
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

#endif


end subroutine product_abc_sca


!=========================================================================
! Calculate C = A^T * B * A for distributed matrices
!
!=========================================================================
subroutine product_transaba_sca(desca,a_matrix_local,descb,b_matrix_local,descc,c_matrix_local)
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
 integer,external :: NUMROC
!=====


 kmat  = desca(M_A)
 mmat  = desca(N_A)
 kmat1 = descb(M_A)
 kmat2 = descb(N_A)
 mmat1 = descc(M_A)
 mmat2 = descc(N_A)

 if( mmat1 /= mmat ) then
   write(msg,*) 'mmat1 /= mmat',mmat1,mmat
   call die('Dimension error in product_transaba_scalapack'//msg)
 endif
 if( mmat2 /= mmat ) then
   write(msg,*) 'mmat2 /= mmat',mmat2,mmat
   call die('Dimension error in product_transaba_scalapack'//msg)
 endif
 if( kmat1 /= kmat ) then
   write(msg,*) 'kmat1 /= kmat',kmat1,kmat
   call die('Dimension error in product_transaba_scalapack'//msg)
 endif
 if( kmat2 /= kmat ) then
   write(msg,*) 'kmat2 /= kmat',kmat2,kmat
   call die('Dimension error in product_transaba_scalapack'//msg)
 endif

#ifdef HAVE_SCALAPACK
 cntxt = desca(CTXT_A)
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
 
#endif


end subroutine product_transaba_sca


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
 integer  :: mlocal,nlocal
 integer  :: n
 integer  :: info
 integer  :: lwork,liwork
!=====

#ifdef HAVE_SCALAPACK
 n = desc(M_A)
 mlocal = SIZE( matrix , DIM=1 )
 nlocal = SIZE( matrix , DIM=2 )
 matrix_inv(:,:) = matrix(:,:)

 allocate(ipiv(block_row+mlocal))

 call PDGETRF(n,n,matrix_inv,1,1,desc,ipiv,info)

 if(info/=0) call die('FAILURE in PDGETRF')

 ! Query
 allocate(work(1),iwork(1))
 lwork  = -1
 liwork = -1
 call PDGETRI(n,matrix_inv,1,1,desc,ipiv,work, lwork, iwork, liwork, info )
 if(info/=0) call die('FAILURE in PDGETRI')
 
 lwork  = NINT(work(1))
 liwork = iwork(1)
 deallocate(work,iwork)
 allocate(work(lwork),iwork(liwork))
 call PDGETRI(n,matrix_inv,1,1,desc,ipiv,work, lwork, iwork, liwork, info )
 if(info/=0) call die('FAILURE in PDGETRI')

 deallocate(ipiv)
 deallocate(work,iwork)
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



#ifdef DEBUG
 write(stdout,'(/,a)')           ' ==== SCALAPACK info'
 write(stdout,'(a)')             '   Squared distribution'
 write(stdout,'(a50,1x,i8)')      'Number of proc:',nprow_sd*npcol_sd
 write(stdout,'(a50,1x,i8,1x,i8)') 'Grid of procs:',nprow_sd,npcol_sd
 write(stdout,'(a)')             '       Row distribution'
 write(stdout,'(a50,1x,i8)')      'Number of proc:',nprow_rd*npcol_rd
 write(stdout,'(a50,1x,i8,1x,i8)') 'Grid of procs:',nprow_rd,npcol_rd
 write(stdout,'(a)')             '    Column distribution'
 write(stdout,'(a50,1x,i8)')      'Number of proc:',nprow_cd*npcol_cd
 write(stdout,'(a50,1x,i8,1x,i8)') 'Grid of procs:',nprow_cd,npcol_cd
 write(stdout,'(/)')
#endif
 
#else
 nprow_sd = 1
 npcol_sd = 1
 iprow_sd = 0
 ipcol_sd = 0
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

 call BLACS_GRIDINFO(cntxt_auxil,nprow_auxil,npcol_auxil,iprow_auxil,ipcol_auxil)
 call xmax_ortho(nprow_auxil)
 call xmax_ortho(npcol_auxil)
#endif


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
 integer,intent(out)         :: desc(NDEL),mlocal,nlocal
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
 case default
   write(stdout,*) 'SCALAPACK distribution type does not exist',distribution
   call die('BUG')
 end select

 write(stdout,'(/,a,i6,a,i6,4x,i6)') ' SCALAPACK info: size of the local matrix for proc #', mlocal,' x ',nlocal,iproc_sca

#else
 desc(:)     = 0
 desc(M_A)   = mglobal
 desc(N_A)   = nglobal
 desc(LLD_A) = mglobal
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
 integer,intent(in)          :: desc(NDEL),ilocal
 integer                     :: rowindex_local_to_global_descriptor
!=====
#ifdef HAVE_SCALAPACK
 integer          :: iprow,ipcol,nprow,npcol
#endif
!=====

#ifdef HAVE_SCALAPACK
 call BLACS_GRIDINFO(desc(CTXT_A),nprow,npcol,iprow,ipcol)
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
 integer,intent(in)          :: desc(NDEL),ilocal
 integer                     :: colindex_local_to_global_descriptor
!=====
#ifdef HAVE_SCALAPACK
 integer          :: iprow,ipcol,nprow,npcol
#endif
!=====

#ifdef HAVE_SCALAPACK
 call BLACS_GRIDINFO(desc(CTXT_A),nprow,npcol,iprow,ipcol)
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
 if( parallel_ham ) then
   if( cntxt_ham > 0 ) call BLACS_GRIDEXIT( cntxt_ham )
 endif
 if( cntxt_auxil > 0 ) call BLACS_GRIDEXIT( cntxt_auxil )
 call BLACS_EXIT( 0 )
#endif

end subroutine finish_scalapack


!=========================================================================
end module m_scalapack
!=========================================================================
