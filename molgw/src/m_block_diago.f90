!=========================================================================
! This file is part of MOLGW.
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


contains


!=========================================================================
subroutine diago_4blocks_sqrt(nmat,amb_matrix,apb_matrix,eigenvalue,bigx,bigy)
 use m_tools 
 implicit none

 integer,intent(in)           :: nmat
 real(prec_td),intent(inout)  :: amb_matrix(nmat,nmat)
 real(prec_td),intent(inout)  :: apb_matrix(nmat,nmat)  ! apb_matrix constains (A+B) in the input, however it used a matrix buffer after
 real(dp),intent(out)         :: eigenvalue(nmat)
 real(prec_td),intent(out)    :: bigx(nmat,nmat),bigy(nmat,nmat)
!=====
 integer                   :: t_ij,t_kl
 real(prec_td),allocatable :: amb_eigval(:),bigomega(:)
!=====

 call start_clock(timing_diago_h2p)

 ! First symmetrize the matrices since only the lower triangle was calculated
 do t_kl=1,nmat
   do t_ij=t_kl+1,nmat
     amb_matrix(t_kl,t_ij) = amb_matrix(t_ij,t_kl)
     apb_matrix(t_kl,t_ij) = apb_matrix(t_ij,t_kl)
   enddo
 enddo



 write(stdout,'(/,a)') ' Performing the block diago with square root of matrices'

 !
 ! Calculate (A-B)^{1/2}
 ! First diagonalize (A-B):
 ! (A-B) R = R D
 ! (A-B) is real symmetric, hence R is orthogonal R^{-1} = tR
 ! (A-B)       = R D tR 
 ! (A-B)^{1/2} = R D^{1/2} tR 
 write(stdout,'(a,i8,a,i8)') ' Diago to get (A - B)^{1/2}                   ',nmat,' x ',nmat
 allocate(amb_eigval(nmat))
 call diagonalize(nmat,amb_matrix,amb_eigval)


 ! bigx contains the (A-B)**1/2
 ! bigy contains the (A-B)**-1/2
 forall(t_kl=1:nmat)
   bigx(:,t_kl) = amb_matrix(:,t_kl)*SQRT(amb_eigval(t_kl))
   bigy(:,t_kl) = amb_matrix(:,t_kl)/SQRT(amb_eigval(t_kl))
 end forall
 deallocate(amb_eigval)

 amb_matrix = TRANSPOSE( amb_matrix )
 bigx(:,:) = MATMUL( bigx(:,:) , amb_matrix(:,:) )
 bigy(:,:) = MATMUL( bigy(:,:) , amb_matrix(:,:) )
 
 ! Use amb_matrix as a temporary matrix here:
 amb_matrix(:,:) = MATMUL( apb_matrix , bigx )
 apb_matrix(:,:)  = MATMUL( bigx, amb_matrix )

! write(stdout,*) 'CC ',matrix_is_symmetric(nmat,apb_matrix)


 write(stdout,'(a,i8,a,i8)') ' Diago (A - B)^{1/2} * (A + B) * (A - B)^{1/2}',nmat,' x ',nmat
 allocate(bigomega(nmat))
 call diagonalize(nmat,apb_matrix,bigomega)

 bigomega(:) = SQRT(bigomega(:))

 forall(t_kl=1:nmat)
   apb_matrix(:,t_kl) = apb_matrix(:,t_kl) / SQRT(bigomega(t_kl))
   eigenvalue(t_kl) = bigomega(t_kl)
 end forall

 ! Save (A-B)**-1/2 in amb_matrix 
 amb_matrix(:,:) = bigy(:,:)

 bigx(:,:) = 0.5_dp * MATMUL( bigx(:,:)   , apb_matrix(:,:) )
 bigy(:,:) = bigx(:,:)

 apb_matrix(:,:) = 0.5_dp * MATMUL( amb_matrix(:,:) , apb_matrix(:,:) )
 forall(t_kl=1:nmat)
   apb_matrix(:,t_kl) = apb_matrix(:,t_kl) * bigomega(t_kl)
 end forall
 deallocate(bigomega)

 ! Finalize Resonant (positive excitations second index from 1 to nmat)
 bigx(:,:) = bigx(:,:) + apb_matrix(:,:)
 bigy(:,:) = bigy(:,:) - apb_matrix(:,:)


 call stop_clock(timing_diago_h2p)


end subroutine diago_4blocks_sqrt


!=========================================================================
subroutine diago_4blocks_chol(nmat,desc_apb,m_apb,n_apb,amb_matrix,apb_matrix,&
                              eigenvalue,desc_x,m_x,n_x,bigx,bigy)
 use m_tools 
 implicit none

 integer,intent(in)     :: nmat,m_apb,n_apb,m_x,n_x
 integer,intent(in)     :: desc_apb(ndel),desc_x(ndel)
 real(dp),intent(inout) :: amb_matrix(m_apb,n_apb),apb_matrix(m_apb,n_apb)
 real(dp),intent(out)   :: eigenvalue(nmat)
 real(dp),intent(out)   :: bigx(m_x,n_x)
 real(dp),intent(out)   :: bigy(m_x,n_x)
!=====
 integer  :: info
 integer  :: lwork,liwork
 real(dp),allocatable :: work(:)
 integer,allocatable :: iwork(:)
!=====

#ifdef HAVE_SCALAPACK
 call start_clock(timing_diago_h2p)

 write(stdout,'(/,a)') ' Performing the block diago with Cholesky'

 allocate(work(1))
 allocate(iwork(1))
 lwork=-1
 liwork=-1
 call pdbssolver1(nmat,apb_matrix,1,1,desc_apb,amb_matrix,1,1,desc_apb,    &
                  eigenvalue,bigx,1,1,desc_x,bigy,                         &
                  work,lwork,iwork,liwork,info)
 if(info/=0) call die('SCALAPACK failed')

 lwork  = NINT(work(1))
 deallocate(work)
 call clean_allocate('Buffer array for SCALAPACK diago',work,lwork)

 liwork = iwork(1)
 deallocate(iwork)
 allocate(iwork(liwork))

 call pdbssolver1(nmat,apb_matrix,1,1,desc_apb,amb_matrix,1,1,desc_apb,    &
                  eigenvalue,bigx,1,1,desc_x,bigy,                         &
                  work,lwork,iwork,liwork,info)
 if(info/=0) call die('SCALAPACK failed')

 call clean_deallocate('Buffer array for SCALAPACK diago',work)
 deallocate(iwork)



 call stop_clock(timing_diago_h2p)

#else
 call die('Cholesky diago cannot run without SCALAPACK')
#endif

end subroutine diago_4blocks_chol


!=========================================================================
subroutine diago_4blocks_rpa_sca(nmat,desc_apb,m_apb,n_apb,amb_diag_rpa,apb_matrix,&
                                 eigenvalue,desc_x,m_x,n_x,bigx)
 use m_tools 
 implicit none

 integer,intent(in)     :: nmat,m_apb,n_apb,m_x,n_x
 integer,intent(in)     :: desc_apb(ndel),desc_x(ndel)
 real(dp),intent(in)    :: amb_diag_rpa(nmat)
 real(dp),intent(inout) :: apb_matrix(m_apb,n_apb)
 real(dp),intent(out)   :: eigenvalue(nmat)
 real(dp),intent(out)   :: bigx(m_x,n_x)
!=====
 integer              :: info
 integer              :: ilocal,jlocal,iglobal,jglobal
 integer              :: lwork
 real(dp),allocatable :: work(:)
 real(dp)             :: amb_diag_sqrt(nmat)
!=====

#ifdef HAVE_SCALAPACK
 call start_clock(timing_diago_h2p)


 write(stdout,'(/,a)') ' Performing the SCALAPACK block diago when A-B is diagonal'

 ! First symmetrize (A+B) 
 ! by adding (A+B)^T
 do jlocal=1,n_apb
   jglobal = colindex_local_to_global('S',jlocal)
   do ilocal=1,m_apb
     iglobal = rowindex_local_to_global('S',ilocal)
     if( iglobal == jglobal ) then
       apb_matrix(ilocal,jlocal) = apb_matrix(ilocal,jlocal) * 0.5_dp
     else if( iglobal < jglobal ) then
       apb_matrix(ilocal,jlocal) = 0.0_dp
     endif
   enddo
 enddo
 bigx(:,:) = apb_matrix(:,:)
 call PDGEADD('T',nmat,nmat,1.d0,bigx,1,1,desc_x,1.d0,apb_matrix,1,1,desc_apb)


 ! Calculate (A-B)^1/2 * (A+B) * (A-B)^1/2

 ! Get (A-B)^1/2
 amb_diag_sqrt(:) = SQRT( amb_diag_rpa(:) )

 ! Calculate (A+B) * (A-B)^1/2
 ! Use bigx as a temporary matrix
 do jlocal=1,n_apb
   jglobal = colindex_local_to_global('S',jlocal)
   bigx(:,jlocal) = apb_matrix(:,jlocal) * amb_diag_sqrt(jglobal)
 enddo

 ! Calculate (A-B)^1/2 * [ (A+B) (A-B)^1/2 ]
 do ilocal=1,m_apb
   iglobal = rowindex_local_to_global('S',ilocal)
   apb_matrix(ilocal,:) = amb_diag_sqrt(iglobal) * bigx(ilocal,:)
 enddo

 ! Diagonalization
 call diagonalize_sca_outofplace(desc_apb,nmat,m_apb,n_apb,apb_matrix,eigenvalue,desc_x,m_x,n_x,bigx)

 eigenvalue(:) = SQRT( eigenvalue(:) )

 ! Normalization
 do jlocal=1,n_apb
   jglobal = colindex_local_to_global('S',jlocal)
   bigx(:,jlocal) = bigx(:,jlocal) / SQRT( eigenvalue(jglobal) )
 enddo

 ! Calculate X + Y = (A-B)^1/2 * Z
 do ilocal=1,m_apb
   iglobal = rowindex_local_to_global('S',ilocal)
   bigx(ilocal,:) = amb_diag_sqrt(iglobal) * bigx(ilocal,:)
 enddo


 call stop_clock(timing_diago_h2p)

#else
 call die('RPA block diago cannot run without SCALAPACK')
#endif

end subroutine diago_4blocks_rpa_sca


!=========================================================================
subroutine diago_4blocks_rpa(nmat,amb_diag_rpa,apb_matrix,eigenvalue,bigx)
 use m_tools 
 implicit none

 integer,intent(in)     :: nmat
 real(dp),intent(in)    :: amb_diag_rpa(nmat)
 real(dp),intent(inout) :: apb_matrix(nmat,nmat)
 real(dp),intent(out)   :: eigenvalue(nmat)
 real(dp),intent(out)   :: bigx(nmat,nmat)
!=====
 integer              :: info
 integer              :: lwork
 real(dp),allocatable :: work(:)
 real(dp)             :: amb_diag_sqrt(nmat)
 integer              :: imat,jmat
!=====

 call start_clock(timing_diago_h2p)


 write(stdout,'(/,a)') ' Performing the block diago when A-B is diagonal'

 ! First symmetrize (A+B)
 do jmat=1,nmat
   do imat=jmat+1,nmat
     apb_matrix(jmat,imat) = apb_matrix(imat,jmat)
   enddo
 enddo



 ! Calculate (A-B)^1/2 * (A+B) * (A-B)^1/2

 ! Get (A-B)^1/2
 amb_diag_sqrt(:) = SQRT( amb_diag_rpa(:) )

 ! Calculate (A+B) * (A-B)^1/2
 ! Use bigx as a temporary matrix
 do jmat=1,nmat
   bigx(:,jmat) = apb_matrix(:,jmat) * amb_diag_sqrt(jmat)
 enddo

 ! Calculate (A-B)^1/2 * [ (A+B) (A-B)^1/2 ]
 do imat=1,nmat
   apb_matrix(imat,:) = amb_diag_sqrt(imat) * bigx(imat,:)
 enddo

 ! Diagonalization
 call diagonalize(nmat,apb_matrix,eigenvalue,bigx)

 eigenvalue(:) = SQRT( eigenvalue(:) )

 ! Normalization
 do jmat=1,nmat
   bigx(:,jmat) = bigx(:,jmat) / SQRT( eigenvalue(jmat) )
 enddo

 ! Calculate X + Y = (A-B)^1/2 * Z
 do imat=1,nmat
   bigx(imat,:) = amb_diag_sqrt(imat) * bigx(imat,:)
 enddo


 call stop_clock(timing_diago_h2p)


end subroutine diago_4blocks_rpa


!=========================================================================
subroutine diago_4blocks_davidson(toldav,nexcitation,nmat,amb_diag_rpa,amb_matrix,apb_matrix,eigenvalue,bigx,bigy)
 use m_tools
 implicit none

 real(dp),intent(in)    :: toldav
 integer,intent(in)     :: nexcitation
 integer,intent(in)     :: nmat
 real(dp),intent(in)    :: amb_diag_rpa(nmat)
 real(dp),intent(inout) :: amb_matrix(nmat,nmat)
 real(dp),intent(inout) :: apb_matrix(nmat,nmat)
 real(dp),intent(out)   :: eigenvalue(nmat)
 real(dp),intent(out)   :: bigx(nmat,nmat)
 real(dp),intent(out)   :: bigy(nmat,nmat)
!=====
 integer,parameter    :: ncycle=20
 integer              :: nbb,nbbc
 integer              :: ibb,jbb
 integer              :: icycle
 integer              :: imat
 real(dp)             :: tolres
 real(dp),allocatable :: bb(:,:)
 real(dp),allocatable :: ql(:,:),qr(:,:)
 real(dp),allocatable :: eigvec_left(:,:),eigvec_right(:,:)
 real(dp),allocatable :: apb_bb(:,:)
 real(dp),allocatable :: amb_bb(:,:)
 real(dp),allocatable :: apb_tilde(:,:),amb_tilde(:,:),c_tilde(:,:)
 real(dp),allocatable :: amb_sqrt_tilde(:,:),amb_sqrt_inv_tilde(:,:)
 real(dp),allocatable :: omega2(:)
 logical,allocatable  :: maskmin(:)
 integer              :: t_ij,t_kl
!=====

 call start_clock(timing_diago_h2p)

 write(stdout,'(/,a,i4)') ' Performing the Davidson block diago for excitation counts ',nexcitation

 ! First symmetrize the matrices since only the lower triangle was calculated
 do t_kl=1,nmat
   do t_ij=t_kl+1,nmat
     amb_matrix(t_kl,t_ij) = amb_matrix(t_ij,t_kl)
     apb_matrix(t_kl,t_ij) = apb_matrix(t_ij,t_kl)
   enddo
 enddo



 ! Maximum dimension of the iterative subspace (tilde matrices)
 nbb = MIN( nexcitation * ( 1 + 2 * ncycle ) , nmat)

 allocate(bb(nmat,nbb))
 allocate(amb_bb(nmat,nbb))
 allocate(apb_bb(nmat,nbb))
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
 do ibb=1,nexcitation
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

 amb_bb(:,1:nbbc) = MATMUL( amb_matrix(:,:) , bb(:,1:nbbc) )
 apb_bb(:,1:nbbc) = MATMUL( apb_matrix(:,:) , bb(:,1:nbbc) )


 do icycle=1,ncycle

   nbbc = nexcitation + 2 * nexcitation * (icycle-1)


   allocate(apb_tilde(nbbc,nbbc),amb_tilde(nbbc,nbbc))
   allocate(c_tilde(nbbc,nbbc))
   allocate(amb_sqrt_tilde(nbbc,nbbc),amb_sqrt_inv_tilde(nbbc,nbbc))
   allocate(omega2(nbbc))
   allocate(eigvec_left(nbbc,nbbc),eigvec_right(nbbc,nbbc))

   apb_tilde(1:nbbc,1:nbbc) = MATMUL( TRANSPOSE(bb(:,1:nbbc)) , apb_bb(:,1:nbbc) )
   amb_tilde(1:nbbc,1:nbbc) = MATMUL( TRANSPOSE(bb(:,1:nbbc)) , amb_bb(:,1:nbbc) )


   ! Calculate the square-root of amb_tilde
   amb_sqrt_tilde(:,:) = amb_tilde(:,:)
   call diagonalize(nbbc,amb_sqrt_tilde,omega2)

   if( ANY( omega2 < 1.0e-10_dp ) ) then
     call die('Matrix (A-B) is not positive definite')
   endif
   amb_sqrt_inv_tilde(:,:) = amb_sqrt_tilde(:,:)
   do ibb=1,nbbc
     amb_sqrt_tilde(:,ibb)     = amb_sqrt_tilde(:,ibb)     * SQRT( SQRT(omega2(ibb)) )
     amb_sqrt_inv_tilde(:,ibb) = amb_sqrt_inv_tilde(:,ibb) / SQRT( SQRT(omega2(ibb)) )
   enddo
   amb_sqrt_tilde(:,:)     = MATMUL( amb_sqrt_tilde     , TRANSPOSE( amb_sqrt_tilde)     )
   amb_sqrt_inv_tilde(:,:) = MATMUL( amb_sqrt_inv_tilde , TRANSPOSE( amb_sqrt_inv_tilde) )


   c_tilde(:,:) = MATMUL( amb_sqrt_tilde , MATMUL( apb_tilde , amb_sqrt_tilde) )

   ! Diagonalize the C matrix
   call diagonalize(nbbc,c_tilde,omega2)


   ! TRANSPOSE the left eigvec
   eigvec_left (:,:) = MATMUL( amb_sqrt_inv_tilde , c_tilde )
   eigvec_right(:,:) = MATMUL( amb_sqrt_tilde , c_tilde )

   do ibb=1,nbbc
     eigvec_left (:,ibb) = eigvec_left (:,ibb) * omega2(ibb)**0.25_dp
     eigvec_right(:,ibb) = eigvec_right(:,ibb) / omega2(ibb)**0.25_dp
   enddo

   !
   ! Calculate the WL and WR of Stratmann and Scuseria
   do ibb=1,nexcitation
     ql(:,ibb) = MATMUL( apb_bb(:,1:nbbc) ,  eigvec_right(:,ibb) ) &
                   - SQRT( omega2(ibb) ) * MATMUL( bb(:,1:nbbc) , eigvec_left(:,ibb) )
     qr(:,ibb) = MATMUL( amb_bb(:,1:nbbc) ,  eigvec_left(:,ibb) )  &  
                   - SQRT( omega2(ibb) ) * MATMUL( bb(:,1:nbbc) , eigvec_right(:,ibb) )
   enddo

   !
   ! Some output
   write(stdout,'(/,a,i4)') ' Davidson iteration ',icycle
   do ibb=1,nexcitation
     write(stdout,'(a,i4,x,f13.6,x,es12.4)') ' Excitation   Energy (eV)  Residual: ',ibb,SQRT(omega2(ibb)) * Ha_eV,&
                  MAX( NORM2(ql(:,ibb)) , NORM2(qr(:,ibb)) )
     tolres = MAXVAL( NORM2(ql(:,:),DIM=1) )
   enddo

   !
   ! Stop the iterations here, because one of the stopping criteria has been met
   if( tolres < toldav ) then
     write(stdout,'(a,es12.4,a,es12.4)') ' Davidson diago converged ',tolres,' is lower than ',toldav
     exit
   endif
   if( icycle == ncycle ) then
     write(stdout,'(a)') ' Maximum iteration number reached',ncycle
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
       bb(imat,nbbc+2*ibb-1) = ql(imat,ibb) / ( SQRT(omega2(ibb)) - amb_diag_rpa(imat) )
       bb(imat,nbbc+2*ibb  ) = qr(imat,ibb) / ( SQRT(omega2(ibb)) - amb_diag_rpa(imat) )
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
   amb_bb(:,nbbc+1:nbbc+2*nexcitation) = MATMUL( amb_matrix(:,:) , bb(:,nbbc+1:nbbc+2*nexcitation) )
   apb_bb(:,nbbc+1:nbbc+2*nexcitation) = MATMUL( apb_matrix(:,:) , bb(:,nbbc+1:nbbc+2*nexcitation) )



   deallocate(apb_tilde,amb_tilde,omega2)
   deallocate(c_tilde)
   deallocate(amb_sqrt_tilde,amb_sqrt_inv_tilde)
   deallocate(eigvec_left,eigvec_right)
 enddo


 ! Calculate X and Y from L and R
 ! Remember
 ! L = | X - Y >
 ! R = | X + Y >
 eigenvalue(1:nexcitation) = SQRT( omega2(1:nexcitation) )
 bigx(:,1:nexcitation) = 0.5_dp *( MATMUL( bb(:,1:nbbc) , eigvec_left(:,1:nexcitation) )  &
                                  +MATMUL( bb(:,1:nbbc) , eigvec_right(:,1:nexcitation) ) )
 bigy(:,1:nexcitation) = 0.5_dp *(-MATMUL( bb(:,1:nbbc) , eigvec_left(:,1:nexcitation) )  &
                                  +MATMUL( bb(:,1:nbbc) , eigvec_right(:,1:nexcitation) ) )


 deallocate(apb_tilde,amb_tilde,omega2)
 deallocate(c_tilde)
 deallocate(amb_sqrt_tilde,amb_sqrt_inv_tilde)
 deallocate(eigvec_left,eigvec_right)
  

 call stop_clock(timing_diago_h2p)


end subroutine diago_4blocks_davidson


end module m_block_diago
!=========================================================================
