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
 ! by adding (A+B)
 do jlocal=1,n_apb
   jglobal = colindex_local_to_global('S',jlocal)
   do ilocal=1,m_apb
     iglobal = rowindex_local_to_global('S',ilocal)
     if( iglobal == jglobal ) then
       apb_matrix(ilocal,jlocal) = apb_matrix(ilocal,jlocal) * 0.5_dp
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


end module m_block_diago
!=========================================================================
