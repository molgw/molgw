!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the routines to diagonalize the RPA, TDDFT or BSE matrix
! with the block form
! (  A  B  ) ( X )    ( X )
! ( -A -B  ) ( Y )  = ( Y ) \Omega
!=========================================================================
module m_block_diago
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_mpi
 use m_scalapack
 use m_tools


contains


!=========================================================================
subroutine diago_4blocks_chol(nmat,desc_apb,m_apb,n_apb,amb_matrix,apb_matrix,&
                              bigomega,desc_x,m_x,n_x,xpy_matrix,xmy_matrix)
 implicit none

 integer,intent(in)     :: nmat,m_apb,n_apb,m_x,n_x
 integer,intent(in)     :: desc_apb(NDEL),desc_x(NDEL)
 real(dp),intent(inout) :: amb_matrix(m_apb,n_apb),apb_matrix(m_apb,n_apb)
 real(dp),intent(out)   :: bigomega(nmat)
 real(dp),intent(out)   :: xpy_matrix(m_x,n_x)
 real(dp),intent(out)   :: xmy_matrix(m_x,n_x)
!=====
 integer  :: info
 integer  :: lwork,liwork
 real(dp),allocatable :: work(:)
 integer,allocatable :: iwork(:)
!=====

 call start_clock(timing_diago_h2p)


#ifdef HAVE_SCALAPACK

 write(stdout,'(/,a)') ' Performing the block diago with Cholesky'

 allocate(work(1))
 allocate(iwork(1))
 lwork=-1
 liwork=-1
 call PDBSSOLVER1(nmat,apb_matrix,1,1,desc_apb,amb_matrix,1,1,desc_apb,    &
                  bigomega,xpy_matrix,1,1,desc_x,xmy_matrix,               &
                  work,lwork,iwork,liwork,info)
 if(info/=0) call die('SCALAPACK failed')

 lwork  = NINT(work(1))
 deallocate(work)
 call clean_allocate('Buffer array for SCALAPACK diago',work,lwork)

 liwork = iwork(1)
 deallocate(iwork)
 allocate(iwork(liwork))

 call PDBSSOLVER1(nmat,apb_matrix,1,1,desc_apb,amb_matrix,1,1,desc_apb,    &
                  bigomega,xpy_matrix,1,1,desc_x,xmy_matrix,               &
                  work,lwork,iwork,liwork,info)
 if(info/=0) call die('SCALAPACK failed')

 call clean_deallocate('Buffer array for SCALAPACK diago',work)
 deallocate(iwork)

 !
 ! the subroutine returns X and Y, but we want (X+Y) and (X-Y)
 ! 
 ! Save X in (A+B)
 call PDLACPY('A',nmat,nmat,xpy_matrix,1,1,desc_x,apb_matrix,1,1,desc_apb)
 ! Save Y in (A-B)
 call PDLACPY('A',nmat,nmat,xmy_matrix,1,1,desc_x,amb_matrix,1,1,desc_apb)

 call PDGEADD('N',nmat,nmat,1.0_dp,amb_matrix,1,1,desc_apb, 1.0_dp,xpy_matrix,1,1,desc_x)
 call PDGEADD('N',nmat,nmat,1.0_dp,apb_matrix,1,1,desc_apb,-1.0_dp,xmy_matrix,1,1,desc_x)



#else

 ! Cholevski decomposition of (A+B) = L * L^T
 call DPOTRF('L',nmat,apb_matrix,nmat,info)

 ! Calculate L^T * (A-B) * L
 call DSYGST(3,'L',nmat,amb_matrix,nmat,apb_matrix,nmat,info)

 ! Diagonalize L^T * (A-B) * L
 lwork = -1
 allocate(work(1))
 call DSYEV('V','L',nmat,amb_matrix,nmat,bigomega,work,lwork,info)
 lwork = NINT(work(1))
 deallocate(work)

 allocate(work(lwork))
 call DSYEV('V','L',nmat,amb_matrix,nmat,bigomega,work,lwork,info)
 deallocate(work)

 bigomega(:) = SQRT( bigomega(:) )

 xpy_matrix(:,:) = amb_matrix(:,:)
 xmy_matrix(:,:) = amb_matrix(:,:)

 ! Calculate L * Z
 call DTRMM('L','L','N','N',nmat,nmat,1.0_dp,apb_matrix,nmat,xmy_matrix,nmat)

 ! Calculate L^{-T} * Z
 call DTRSM('L','L','T','N',nmat,nmat,1.0_dp,apb_matrix,nmat,xpy_matrix,nmat)

 !
 ! X-Y = L * Z / Omega^{1/2}
 ! X+Y = L^{-T} * Z * Omega^{1/2}
 block
  integer :: imat
  forall(imat=1:nmat)
    xpy_matrix(:,imat) = xpy_matrix(:,imat) * SQRT( bigomega(imat) )
    xmy_matrix(:,imat) = xmy_matrix(:,imat) / SQRT( bigomega(imat) )
  end forall
 end block


#endif

 call stop_clock(timing_diago_h2p)

end subroutine diago_4blocks_chol


!=========================================================================
subroutine diago_4blocks_rpa_sca(nmat,desc_apb,m_apb,n_apb,amb_diag_rpa,apb_matrix,&
                                 bigomega,desc_x,m_x,n_x,xpy_matrix)
#ifdef HAVE_ELPA
 use elpa1
 use elpa2
 use elpa
#endif
 implicit none

 integer,intent(in)     :: nmat,m_apb,n_apb,m_x,n_x
 integer,intent(in)     :: desc_apb(NDEL),desc_x(NDEL)
 real(dp),intent(in)    :: amb_diag_rpa(nmat)
 real(dp),intent(inout) :: apb_matrix(m_apb,n_apb)
 real(dp),intent(out)   :: bigomega(nmat)
 real(dp),intent(out)   :: xpy_matrix(m_x,n_x)
!=====
 real(dp)             :: amb_diag_sqrt(nmat)
#ifdef HAVE_ELPA
 logical         :: success
 integer         :: comm_row,comm_col
#endif
!=====

 call start_clock(timing_diago_h2p)


 write(stdout,'(/,a)') ' Performing the SCALAPACK block diago when A-B is diagonal'

 ! First symmetrize (A+B) (lower triangular matrix)
 ! Use X+Y as a work buffer
 call symmetrize_matrix_sca('L',nmat,desc_apb,apb_matrix,desc_x,xpy_matrix)

 ! Set (A-B)^1/2
 amb_diag_sqrt(:) = SQRT( amb_diag_rpa(:) )

 !
 ! Prepare (A-B)^1/2 * (A+B) * (A-B)^1/2
 !

 ! Calculate (A+B) * (A-B)^1/2
 call matmul_diag_sca('R',amb_diag_sqrt,desc_apb,apb_matrix)

 ! Calculate (A-B)^1/2 * [ (A+B) * (A-B)^1/2 ]
 call matmul_diag_sca('L',amb_diag_sqrt,desc_apb,apb_matrix)


 ! Diagonalization
#ifdef HAVE_ELPA
 block
   integer :: info
   info = get_elpa_communicators(comm_world,iprow_sd,ipcol_sd,comm_row,comm_col)
   success = elpa_solve_evp_real(nmat,nmat,apb_matrix,m_apb,bigomega,xpy_matrix,m_x,desc_apb(MB_),n_apb, &
                                 comm_row,comm_col,comm_world,method='2stage')
 end block
#else
 call diagonalize_sca_pdsyevr(nmat,desc_apb,apb_matrix,bigomega,desc_x,xpy_matrix)
#endif

 bigomega(:) = SQRT( bigomega(:) )


 !
 ! Prepare  ( X + Y ) = (A-B)^{1/2} * Z * \Omega^{-1/2}
 !

 ! Calculate Z * \Omega^{-1/2}
 call matmul_diag_sca('R',1.0_dp/SQRT(bigomega(:)),desc_x,xpy_matrix)

 ! Calculate (A-B)^{1/2} * [ Z * \Omega^{-1/2} ]
 call matmul_diag_sca('L',amb_diag_sqrt,desc_x,xpy_matrix)


 call stop_clock(timing_diago_h2p)


end subroutine diago_4blocks_rpa_sca


!=========================================================================
subroutine diago_4blocks_davidson(toldav,nexcitation,nmat,amb_diag_rpa, &
                                  desc_apb,m_apb,n_apb,amb_matrix,apb_matrix, &
                                  bigomega,desc_x,m_x,n_x,xpy_matrix,xmy_matrix)
 implicit none

 real(dp),intent(in)    :: toldav
 integer,intent(in)     :: nexcitation
 integer,intent(in)     :: nmat,m_apb,n_apb,m_x,n_x
 integer,intent(in)     :: desc_apb(NDEL),desc_x(NDEL)
 real(dp),intent(in)    :: amb_diag_rpa(nmat)
 real(dp),intent(inout) :: amb_matrix(m_apb,n_apb)
 real(dp),intent(inout) :: apb_matrix(m_apb,n_apb)
 real(dp),intent(out)   :: bigomega(nmat)
 real(dp),intent(out)   :: xpy_matrix(m_x,n_x)
 real(dp),intent(out)   :: xmy_matrix(m_x,n_x)
!=====
 integer              :: descb(NDEL),desce(NDEL)
 integer,parameter    :: SMALL_BLOCK=4
 integer,parameter    :: NCYCLE=20
 integer              :: nbb,nbbc,nbba
 integer              :: ibb,jbb
 integer              :: icycle
 integer              :: imat
 integer              :: iglobal,jglobal,ib,jb
 integer              :: mb,nb
 integer              :: me,ne
 integer              :: info
 real(dp)             :: tolres
 real(dp),allocatable :: ab_local(:,:),bb_local(:,:)
 real(dp),allocatable :: ev_local(:,:)
 real(dp),allocatable :: bb(:,:)
 real(dp),allocatable :: ql(:,:),qr(:,:)
 real(dp),allocatable :: eigvec_left(:,:),eigvec_right(:,:)
 real(dp),allocatable :: bb_apb_bb(:,:),bb_amb_bb(:,:)
 real(dp),allocatable :: apb_bb(:,:),amb_bb(:,:)
 real(dp),allocatable :: apb_tilde(:,:),amb_tilde(:,:),c_tilde(:,:)
 real(dp),allocatable :: amb_sqrt_tilde(:,:),amb_sqrt_inv_tilde(:,:)
 real(dp),allocatable :: bigomega_tmp(:)
 logical,allocatable  :: maskmin(:)
 integer              :: lwork
 real(dp),allocatable :: work(:)
!=====

 call start_clock(timing_diago_h2p)

 write(stdout,'(/,a,i4)') ' Performing the Davidson block diago'

 !
 ! No need to symmetrize (A+B) and (A-B), since we call the specialized DSYMM or PDSYMM

 !
 ! Maximum dimension of the iterative subspace (tilde matrices)
 nbb = MIN( nexcitation * ( 1 + 2 * NCYCLE ) , nmat)

 !
 ! Allocate the small matrices
 allocate(bb(nmat,nbb))
 allocate(amb_bb(nmat,nbb))
 allocate(apb_bb(nmat,nbb))
 allocate(bb_amb_bb(nbb,nbb))
 allocate(bb_apb_bb(nbb,nbb))
 allocate(ql(nmat,nexcitation))
 allocate(qr(nmat,nexcitation))

 ! Current dimension of the iterative subspace
 nbbc = nexcitation 

 !
 ! Initialize with a stupid guess based on the diagonal
 !
 bb(:,1:nbbc)=0.001_dp
 ! Find the nexcitation lowest diagonal elements
 allocate(maskmin(nmat))
 maskmin(:)=.TRUE.
 do ibb=1,nbbc
   jbb = MINLOC( amb_diag_rpa(:) ,DIM=1, MASK=maskmin(:))
   bb(jbb,ibb) = 1.0_dp
   maskmin(jbb) = .FALSE. 
 enddo
 deallocate(maskmin)
 ! Gram-Schmidt orthonormalization of the initial guesses
 do ibb=1,nbbc
   !
   ! Orthogonalize to previous vectors
   do jbb=1,ibb-1
     bb(:,ibb) = bb(:,ibb) - bb(:,jbb) * DOT_PRODUCT( bb(:,ibb) , bb(:,jbb) )
   enddo
   !
   ! Normalize
   bb(:,ibb) = bb(:,ibb) / NORM2( bb(:,ibb) )
 enddo

 ! 
 ! The time-consumming operation:
 ! Calculate (A-B) b  and  (A+B) b
#ifndef HAVE_SCALAPACK
! amb_bb(:,1:nbbc) = MATMUL( amb_matrix(:,:) , bb(:,1:nbbc) )
! apb_bb(:,1:nbbc) = MATMUL( apb_matrix(:,:) , bb(:,1:nbbc) )
 call DSYMM('L','L',nmat,nbbc,1.0_dp,amb_matrix,nmat,bb(:,1:nbbc),nmat,0.0_dp,amb_bb(:,1:nbbc),nmat)
 call DSYMM('L','L',nmat,nbbc,1.0_dp,apb_matrix,nmat,bb(:,1:nbbc),nmat,0.0_dp,apb_bb(:,1:nbbc),nmat)
#else

  !
  ! Distribute bb
  mb = NUMROC(nmat,block_row,iprow_sd,first_row,nprow_sd)
  nb = NUMROC(nbbc,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
  call DESCINIT(descb,nmat,nbbc,block_row,SMALL_BLOCK,first_row,first_col,cntxt_sd,MAX(1,mb),info)

  allocate(bb_local(mb,nb),ab_local(mb,nb))
  do jb=1,nb
    jglobal = INDXL2G(jb,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
    do ib=1,mb
      iglobal = INDXL2G(ib,block_row,iprow_sd,first_row,nprow_sd)
      bb_local(ib,jb) = bb(iglobal,jglobal)
    enddo
  enddo

  !
  ! Calculate (A-B) b
  call PDSYMM('L','L',nmat,nbbc,              &
              1.0_dp,amb_matrix,1,1,desc_apb, &
              bb_local,1,1,descb,             &
              0.0_dp,ab_local,1,1,descb)

  amb_bb(:,1:nbbc) = 0.0_dp
  do jb=1,nb
    jglobal = INDXL2G(jb,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
    do ib=1,mb
      iglobal = INDXL2G(ib,block_row,iprow_sd,first_row,nprow_sd)
      amb_bb(iglobal,jglobal) = ab_local(ib,jb)
    enddo
  enddo
  call xsum_world(amb_bb(:,1:nbbc))

  !
  ! Calculate (A+B) b
  call PDSYMM('L','L',nmat,nbbc,              &
              1.0_dp,apb_matrix,1,1,desc_apb, &
              bb_local,1,1,descb,             &
              0.0_dp,ab_local,1,1,descb)

  apb_bb(:,1:nbbc) = 0.0_dp
  do jb=1,nb
    jglobal = INDXL2G(jb,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
    do ib=1,mb
      iglobal = INDXL2G(ib,block_row,iprow_sd,first_row,nprow_sd)
      apb_bb(iglobal,jglobal) = ab_local(ib,jb)
    enddo
  enddo
  call xsum_world(apb_bb(:,1:nbbc))


  deallocate(bb_local,ab_local)
#endif


 ! Calculate and store   b (A+B) b = b^T [ (A+B) b ]
 call DGEMM('T','N',nbbc,nbbc,nmat,1.0_dp,bb(:,1:nbbc),nmat,apb_bb(:,1:nbbc),nmat,0.0_dp,bb_apb_bb(1:nbbc,1:nbbc),nbbc)
 ! Calculate and store   b (A-B) b = b^T [ (A-B) b ]
 call DGEMM('T','N',nbbc,nbbc,nmat,1.0_dp,bb(:,1:nbbc),nmat,amb_bb(:,1:nbbc),nmat,0.0_dp,bb_amb_bb(1:nbbc,1:nbbc),nbbc)


 do icycle=1,NCYCLE

   allocate(apb_tilde(nbbc,nbbc),amb_tilde(nbbc,nbbc))
   allocate(c_tilde(nbbc,nbbc))
   allocate(amb_sqrt_tilde(nbbc,nbbc),amb_sqrt_inv_tilde(nbbc,nbbc))
   allocate(bigomega_tmp(nbbc))
   allocate(eigvec_left(nbbc,nbbc),eigvec_right(nbbc,nbbc))


   ! Copy the saved matrices into the work matrices
   apb_tilde(:,:) = bb_apb_bb(1:nbbc,1:nbbc)
   amb_tilde(:,:) = bb_amb_bb(1:nbbc,1:nbbc)


   ! Cholevski decomposition of (A-B) = L * L^T
   call DPOTRF('L',nbbc,amb_tilde,nbbc,info)
   if( info /= 0 ) then
     call die('Matrix (A-B) is not positive definite')
   endif

   ! Calculate L^T * (A+B) * L
   call DSYGST(3,'L',nbbc,apb_tilde,nbbc,amb_tilde,nbbc,info)

   ! Diagonalize the matrix L^T * (A+B) * L = Z \Omega Z^T
   lwork = -1
   allocate(work(1))
   call DSYEV('V','L',nbbc,apb_tilde,nbbc,bigomega_tmp,work,lwork,info)
   lwork = NINT(work(1))
   deallocate(work)
  
   allocate(work(lwork))
   call DSYEV('V','L',nbbc,apb_tilde,nbbc,bigomega_tmp,work,lwork,info)
   deallocate(work)

   bigomega_tmp(:) = SQRT( bigomega_tmp(:) )

   eigvec_left(:,:)  = apb_tilde(:,:)
   eigvec_right(:,:) = apb_tilde(:,:)

   ! Calculate L * Z
   call DTRMM('L','L','N','N',nbbc,nbbc,1.0_dp,amb_tilde,nbbc,eigvec_right,nbbc)
  
   ! Calculate L^{-T} * Z
   call DTRSM('L','L','T','N',nbbc,nbbc,1.0_dp,amb_tilde,nbbc,eigvec_left,nbbc)
  
   !
   ! X-Y = L * Z / Omega^{1/2}
   ! X+Y = L^{-T} * Z * Omega^{1/2}
   forall(ibb=1:nexcitation)
     eigvec_left(:,ibb)  = eigvec_left(:,ibb)  * SQRT( bigomega_tmp(ibb) )
     eigvec_right(:,ibb) = eigvec_right(:,ibb) / SQRT( bigomega_tmp(ibb) )
   end forall

   !
   ! Calculate the WL and WR of Stratmann and Scuseria
   do ibb=1,nexcitation
     ql(:,ibb) = MATMUL( apb_bb(:,1:nbbc) ,  eigvec_right(:,ibb) ) &
                   - bigomega_tmp(ibb) * MATMUL( bb(:,1:nbbc) , eigvec_left(:,ibb) )
     qr(:,ibb) = MATMUL( amb_bb(:,1:nbbc) ,  eigvec_left(:,ibb) )  &  
                   - bigomega_tmp(ibb) * MATMUL( bb(:,1:nbbc) , eigvec_right(:,ibb) )
   enddo

   !
   ! Some output at each cycle
   write(stdout,'(/,a,i4)') ' Davidson iteration ',icycle
   write(stdout,'(a,i5)')   ' Iterative subspace dimension ',nbbc
   do ibb=1,nexcitation
     write(stdout,'(a,i4,1x,f13.6,1x,es12.4)') ' Excitation   Energy (eV)  Residual: ',ibb,bigomega_tmp(ibb) * Ha_eV,&
                  MAX( NORM2(ql(:,ibb)) , NORM2(qr(:,ibb)) )
     tolres = MAXVAL( NORM2(ql(:,:),DIM=1) )
   enddo

   !
   ! Stop the iterations here, because one of the stopping criteria has been met
   if( tolres < toldav ) then
     write(stdout,'(a,es12.4,a,es12.4)') ' Davidson diago converged ',tolres,' is lower than ',toldav
     exit
   endif
   if( icycle == NCYCLE ) then
     write(stdout,'(a,1x,i4)')           ' Maximum iteration number reached',NCYCLE
     write(stdout,'(a,es12.4,a,es12.4)') ' Davidson diago not converged ',tolres,' is larger than ',toldav
     call issue_warning('TDDFT or BSE Davidson diago not fully converged')
     exit
   endif
   if( nbbc + 2 * nexcitation > nmat ) then
     write(stdout,'(a,i6)') ' Iterative subspace is larger than the transition space ',nmat
     write(stdout,'(a,es12.4,a,es12.4)') ' Davidson diago not converged ',tolres,' is larger than ',toldav
     call issue_warning('TDDFT or BSE Davidson diago not fully converged')
     exit
   endif

   !
   ! Preconditioning of the residual to get the new trial vectors
   do ibb=1,nexcitation
     do imat=1,nmat
       bb(imat,nbbc+2*ibb-1) = ql(imat,ibb) / ( bigomega_tmp(ibb) - amb_diag_rpa(imat) )
       bb(imat,nbbc+2*ibb  ) = qr(imat,ibb) / ( bigomega_tmp(ibb) - amb_diag_rpa(imat) )
     enddo
   enddo

   !
   ! Orthogonalize to all previous
   do ibb=nbbc+1,nbbc+2*nexcitation
     do jbb=1,ibb-1
       bb(:,ibb) = bb(:,ibb) - bb(:,jbb) * DOT_PRODUCT( bb(:,ibb) , bb(:,jbb) )
     enddo
     ! Normalize
     bb(:,ibb) = bb(:,ibb) / NORM2( bb(:,ibb) )
   enddo

   !
   ! Calculate the new  (A-B) * b   and   (A+B) * b

   !  Dimensions to be added
   nbba = 2 * nexcitation
#ifndef HAVE_SCALAPACK
!   amb_bb(:,nbbc+1:nbbc+nbba) = MATMUL( amb_matrix(:,:) , bb(:,nbbc+1:nbbc+2*nexcitation) )
!   apb_bb(:,nbbc+1:nbbc+nbba) = MATMUL( apb_matrix(:,:) , bb(:,nbbc+1:nbbc+2*nexcitation) )
   call DSYMM('L','L',nmat,nbba,1.0_dp,amb_matrix,nmat,bb(:,nbbc+1:nbbc+nbba),nmat,0.0_dp,amb_bb(:,nbbc+1:nbbc+nbba),nmat)
   call DSYMM('L','L',nmat,nbba,1.0_dp,apb_matrix,nmat,bb(:,nbbc+1:nbbc+nbba),nmat,0.0_dp,apb_bb(:,nbbc+1:nbbc+nbba),nmat)
#else

    !
    ! Distribute bb
    mb = NUMROC(nmat,block_row,iprow_sd,first_row,nprow_sd)
    nb = NUMROC(nbba,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
    call DESCINIT(descb,nmat,nbba,block_row,SMALL_BLOCK,first_row,first_col,cntxt_sd,MAX(1,mb),info)
  
    allocate(bb_local(mb,nb),ab_local(mb,nb))
    do jb=1,nb
      jglobal = INDXL2G(jb,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
      do ib=1,mb
        iglobal = INDXL2G(ib,block_row,iprow_sd,first_row,nprow_sd)
        bb_local(ib,jb) = bb(iglobal,nbbc+jglobal)
      enddo
    enddo
  
    !
    ! Calculate (A-B) b
    call PDSYMM('L','L',nmat,nbba,              &
                1.0_dp,amb_matrix,1,1,desc_apb, &
                bb_local,1,1,descb,             &
                0.0_dp,ab_local,1,1,descb)
  
    amb_bb(:,nbbc+1:nbbc+nbba) = 0.0_dp
    do jb=1,nb
      jglobal = INDXL2G(jb,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
      do ib=1,mb
        iglobal = INDXL2G(ib,block_row,iprow_sd,first_row,nprow_sd)
        amb_bb(iglobal,nbbc+jglobal) = ab_local(ib,jb)
      enddo
    enddo
    call xsum_world(amb_bb(:,nbbc+1:nbbc+nbba))
  
    !
    ! Calculate (A+B) b
    call PDSYMM('L','L',nmat,nbba,              &
                1.0_dp,apb_matrix,1,1,desc_apb, &
                bb_local,1,1,descb,             &
                0.0_dp,ab_local,1,1,descb)
  
    apb_bb(:,nbbc+1:nbbc+nbba) = 0.0_dp
    do jb=1,nb
      jglobal = INDXL2G(jb,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
      do ib=1,mb
        iglobal = INDXL2G(ib,block_row,iprow_sd,first_row,nprow_sd)
        apb_bb(iglobal,nbbc+jglobal) = ab_local(ib,jb)
      enddo
    enddo
    call xsum_world(apb_bb(:,nbbc+1:nbbc+nbba))
  
  
    deallocate(bb_local,ab_local)
#endif

   ! Add the missing b ( A+B) b
   bb_apb_bb(nbbc+1:nbbc+nbba,     1:nbbc+nbba) = MATMUL( TRANSPOSE(bb(:,nbbc+1:nbbc+nbba)) , apb_bb(:,     1:nbbc+nbba) )
   bb_apb_bb(     1:nbbc+nbba,nbbc+1:nbbc+nbba) = MATMUL( TRANSPOSE(bb(:,     1:nbbc+nbba)) , apb_bb(:,nbbc+1:nbbc+nbba) )
   bb_apb_bb(nbbc+1:nbbc+nbba,nbbc+1:nbbc+nbba) = MATMUL( TRANSPOSE(bb(:,nbbc+1:nbbc+nbba)) , apb_bb(:,nbbc+1:nbbc+nbba) )

   ! Add the missing b ( A-B) b
   bb_amb_bb(nbbc+1:nbbc+nbba,     1:nbbc+nbba) = MATMUL( TRANSPOSE(bb(:,nbbc+1:nbbc+nbba)) , amb_bb(:,     1:nbbc+nbba) )
   bb_amb_bb(     1:nbbc+nbba,nbbc+1:nbbc+nbba) = MATMUL( TRANSPOSE(bb(:,     1:nbbc+nbba)) , amb_bb(:,nbbc+1:nbbc+nbba) )
   bb_amb_bb(nbbc+1:nbbc+nbba,nbbc+1:nbbc+nbba) = MATMUL( TRANSPOSE(bb(:,nbbc+1:nbbc+nbba)) , amb_bb(:,nbbc+1:nbbc+nbba) )


   !
   ! Set the new dimension for the tilde matrices
   nbbc = nbbc + nbba


   deallocate(apb_tilde,amb_tilde,bigomega_tmp)
   deallocate(c_tilde)
   deallocate(amb_sqrt_tilde,amb_sqrt_inv_tilde)
   deallocate(eigvec_left,eigvec_right)
 enddo



 ! Calculate X and Y from L and R
 ! Remember
 ! L = | X - Y >
 ! R = | X + Y >
 bigomega(1:nexcitation) = bigomega_tmp(1:nexcitation)

#ifndef HAVE_SCALAPACK

 xpy_matrix(:,1:nexcitation) = MATMUL( bb(:,1:nbbc) , eigvec_right(:,1:nexcitation) )
 xmy_matrix(:,1:nexcitation) = MATMUL( bb(:,1:nbbc) , eigvec_left (:,1:nexcitation) )

#else

 !
 ! 0. Distribute bb
 !
 mb = NUMROC(nmat,block_row,iprow_sd,first_row,nprow_sd)
 nb = NUMROC(nbbc,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
 call DESCINIT(descb,nmat,nbbc,block_row,SMALL_BLOCK,first_row,first_col,cntxt_sd,MAX(1,mb),info)
  
 allocate(bb_local(mb,nb))
 do jb=1,nb
   jglobal = INDXL2G(jb,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
   do ib=1,mb
     iglobal = INDXL2G(ib,block_row,iprow_sd,first_row,nprow_sd)
     bb_local(ib,jb) = bb(iglobal,jglobal)
   enddo
 enddo

 me = NUMROC(nbbc,SMALL_BLOCK,iprow_sd,first_row,nprow_sd)
 ne = NUMROC(nexcitation,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
 call DESCINIT(desce,nbbc,nexcitation,SMALL_BLOCK,SMALL_BLOCK,first_row,first_col,cntxt_sd,MAX(1,me),info)
 allocate(ev_local(me,ne))
 !
 ! 1. Deal with (X-Y)
 !
 do jb=1,ne
   jglobal = INDXL2G(jb,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
   do ib=1,me
     iglobal = INDXL2G(ib,SMALL_BLOCK,iprow_sd,first_row,nprow_sd)
     ev_local(ib,jb) = eigvec_left(iglobal,jglobal)
   enddo
 enddo

 call PDGEMM('N','N',nmat,nexcitation,nbbc,  &
             1.0_dp,bb_local,1,1,descb,      &
             ev_local,1,1,desce,             &
             0.0_dp,xmy_matrix,1,1,desc_x)

 !
 ! 2. Deal with (X+Y)
 !
 do jb=1,ne
   jglobal = INDXL2G(jb,SMALL_BLOCK,ipcol_sd,first_col,npcol_sd)
   do ib=1,me
     iglobal = INDXL2G(ib,SMALL_BLOCK,iprow_sd,first_row,nprow_sd)
     ev_local(ib,jb) = eigvec_right(iglobal,jglobal)
   enddo
 enddo

 call PDGEMM('N','N',nmat,nexcitation,nbbc,  &
             1.0_dp,bb_local,1,1,descb,      &
             ev_local,1,1,desce,             &
             0.0_dp,xpy_matrix,1,1,desc_x)

 deallocate(bb_local,ev_local)


#endif


 deallocate(bb_apb_bb,bb_amb_bb)
 deallocate(apb_tilde,amb_tilde,bigomega_tmp)
 deallocate(c_tilde)
 deallocate(amb_sqrt_tilde,amb_sqrt_inv_tilde)
 deallocate(eigvec_left,eigvec_right)
  

 call stop_clock(timing_diago_h2p)


end subroutine diago_4blocks_davidson


end module m_block_diago
!=========================================================================
