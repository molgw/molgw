!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the linear algebra routines
!
!=========================================================================
module m_linear_algebra
 use m_definitions
 use m_warning,only: die

 interface invert
   module procedure invert_dp
   module procedure invert_inplace_dp
   module procedure invert_cdp
   module procedure invert_inplace_cdp
 end interface

 interface diagonalize_wo_vectors
   module procedure diagonalize_wo_vectors_dp
 end interface

 interface diagonalize
   module procedure diagonalize_cdp
   module procedure diagonalize_dp
   module procedure diagonalize_sp
   module procedure diagonalize_inplace_cdp
   module procedure diagonalize_inplace_dp
   module procedure diagonalize_inplace_sp
 end interface


contains


!=========================================================================
! Override the upper part of a matrix with its lower part
! Final matrix is symmetric
subroutine matrix_lower_to_full(matrix)
 implicit none
 real(dp),intent(inout)  :: matrix(:,:)
!=====
 integer :: imat,jmat,nmat
!=====
 nmat = SIZE(matrix,DIM=1)
 do jmat=1,nmat
   do imat=jmat+1,nmat
     matrix(jmat,imat) = matrix(imat,jmat)
   enddo
 enddo

end subroutine matrix_lower_to_full


!=========================================================================
! Calculate the trace of a real matrix
function matrix_trace(matrix)
 implicit none
 real(dp),intent(in) :: matrix(:,:)
 real(dp)            :: matrix_trace
!=====
 integer :: n1,i1
!=====

 n1 = SIZE( matrix , DIM=1 )
 if( n1 /= SIZE( matrix , DIM=2 ) ) call die('matrix_trace: non square matrix')

 matrix_trace = 0.0_dp
 do i1=1,n1
   matrix_trace = matrix_trace + matrix(i1,i1)
 enddo

end function matrix_trace


!=========================================================================
! Calculate the trace of a complex matrix
function matrix_trace_cmplx(matrix)
 implicit none
 complex(dp),intent(in) :: matrix(:,:)
 complex(dp)            :: matrix_trace_cmplx
!=====
 integer :: n1,i1
!=====

 n1 = SIZE( matrix , DIM=1 )
 if( n1 /= SIZE( matrix , DIM=2 ) ) call die('matrix_trace: non square matrix')

 matrix_trace_cmplx = ( 0.0_dp, 0.0_dp )
 do i1=1,n1
   matrix_trace_cmplx = matrix_trace_cmplx + matrix(i1,i1)
 enddo

end function matrix_trace_cmplx


!=========================================================================
function matrix_is_symmetric(matrix)
 implicit none
 logical             :: matrix_is_symmetric
 real(dp),intent(in) :: matrix(:,:)
!=====
 integer :: imat,jmat
!=====

 matrix_is_symmetric = .TRUE.
 do imat=1,SIZE(matrix,DIM=1)
   do jmat=1,imat-1
     if( ABS( matrix(imat,jmat) - matrix(jmat,imat) ) > 1.0e-5_dp ) then
       matrix_is_symmetric=.FALSE.
       return
     endif
   enddo
 enddo

end function matrix_is_symmetric


!=========================================================================
subroutine invert_dp(matrix,matrix_inv)
 implicit none

 real(dp),intent(in)  :: matrix(:,:)
 real(dp),intent(out) :: matrix_inv(:,:)
!=====
 integer  :: nmat
 real(dp),allocatable :: work(:)
 integer,allocatable  :: ipiv(:)
 integer              :: info
!=====

 nmat = SIZE( matrix(:,:) , DIM=1)
 allocate(work(nmat))
 allocate(ipiv(nmat))

 matrix_inv(:,:) = matrix(:,:)

 call DGETRF(nmat,nmat,matrix_inv,nmat,ipiv,info)
 if(info/=0) call die('FAILURE in DGETRF')

 call DGETRI(nmat,matrix_inv,nmat,ipiv,work,nmat,info)
 if(info/=0) call die('FAILURE in DGETRI')

 deallocate(work,ipiv)

end subroutine invert_dp


!=========================================================================
subroutine invert_inplace_dp(matrix)
 implicit none

 real(dp),intent(inout) :: matrix(:,:)
!=====
 integer              :: nmat
 real(dp),allocatable :: work(:)
 integer,allocatable  :: ipiv(:)
 integer              :: info
!=====

 nmat = SIZE( matrix(:,:) , DIM=1)
 allocate(work(nmat))
 allocate(ipiv(nmat))

 call DGETRF(nmat,nmat,matrix,nmat,ipiv,info)
 if(info/=0) call die('FAILURE in DGETRF')

 call DGETRI(nmat,matrix,nmat,ipiv,work,nmat,info)
 if(info/=0) call die('FAILURE in DGETRI')

 deallocate(work,ipiv)

end subroutine invert_inplace_dp


!=========================================================================
subroutine invert_cdp(matrix,matrix_inv)
 implicit none

 complex(dp),intent(in)  :: matrix(:,:)
 complex(dp),intent(out) :: matrix_inv(:,:)
!=====
 integer                 :: nmat
 complex(dp),allocatable :: work(:)
 integer,allocatable     :: ipiv(:)
 integer                 :: info
!=====

 nmat = SIZE( matrix(:,:) , DIM=1)
 allocate(work(nmat))
 allocate(ipiv(nmat))

 matrix_inv(:,:) = matrix(:,:)

 call ZGETRF(nmat,nmat,matrix_inv,nmat,ipiv,info)
 if(info/=0) call die('FAILURE in ZGETRF')

 call ZGETRI(nmat,matrix_inv,nmat,ipiv,work,nmat,info)
 if(info/=0) call die('FAILURE in ZGETRI')

 deallocate(work,ipiv)

end subroutine invert_cdp


!=========================================================================
subroutine invert_inplace_cdp(matrix)
 implicit none
 complex(dp),intent(in) :: matrix(:,:)
!=====
 integer                 :: nmat
 complex(dp),allocatable :: work(:)
 integer,allocatable     :: ipiv(:)
 integer                 :: info
!=====

 nmat = SIZE(matrix,DIM=1)
 allocate(work(nmat))
 allocate(ipiv(nmat))

 call ZGETRF(nmat,nmat,matrix,nmat,ipiv,info)
 if(info/=0) call die('FAILURE in ZGETRF')

 call ZGETRI(nmat,matrix,nmat,ipiv,work,nmat,info)
 if(info/=0) call die('FAILURE in ZGETRI')


 deallocate(work,ipiv)

end subroutine invert_inplace_cdp


!=========================================================================
subroutine diagonalize_wo_vectors_dp(flavor,matrix,eigval)
 implicit none
 character(len=1),intent(in) :: flavor
 real(dp),intent(inout)      :: matrix(:,:)
 real(dp),intent(out)        :: eigval(:)
!=====
 integer :: nmat
 integer :: info
 real(dp),allocatable :: work(:)
!=====

 nmat = SIZE(matrix,DIM=1)
 allocate(work(3*nmat-1))

 call DSYEV('N','L',nmat,matrix,nmat,eigval,work,3*nmat-1,info)

 deallocate(work)


end subroutine diagonalize_wo_vectors_dp


!=========================================================================
subroutine diagonalize_dp(flavor,matrix,eigval,eigvec)
 implicit none
 character(len=1),intent(in) :: flavor
 real(dp),intent(in)         :: matrix(:,:)
 real(dp),intent(out)        :: eigval(:)
 real(dp),intent(out)        :: eigvec(:,:)
!=====
 integer              :: nmat,lwork
 real(dp),allocatable :: work(:)
 integer              :: info
 integer              :: liwork
 integer,allocatable  :: iwork(:)
 integer,allocatable  :: isuppz(:)
 real(dp),allocatable :: matrix_tmp(:,:)
 real(dp)             :: ABSTOL
 real(dp),external    :: DLAMCH
 integer              :: neig
!=====

 nmat = SIZE(matrix,DIM=1)

 eigvec(:,:) = matrix(:,:)

 lwork = -1
 allocate(work(1))

 select case(flavor)
 case('r','R')
   allocate(iwork(1))
   allocate(matrix_tmp(nmat,nmat))
   allocate(isuppz(2*nmat))
   ABSTOL = DLAMCH('S')
   matrix_tmp(:,:) = matrix(:,:)
   call DSYEVR('V','A','L',nmat,matrix_tmp,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,iwork,liwork,info)
   liwork = iwork(1)
   deallocate(iwork)
 case('d','D')
   allocate(iwork(1))
   call DSYEVD('V','L',nmat,eigvec,nmat,eigval,work,lwork,iwork,liwork,info)
   liwork = iwork(1)
   deallocate(iwork)
 case default
   call DSYEV('V','L',nmat,eigvec,nmat,eigval,work,lwork,info)
 end select

 lwork = NINT(work(1))
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_dp: diago failure 1')

 allocate(work(lwork))

 select case(flavor)
 case('r','R')
   allocate(iwork(liwork))
   call DSYEVR('V','A','L',nmat,matrix_tmp,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,iwork,liwork,info)
   deallocate(iwork,isuppz)
   deallocate(matrix_tmp)
 case('d','D')
   allocate(iwork(liwork))
   call DSYEVD('V','L',nmat,eigvec,nmat,eigval,work,lwork,iwork,liwork,info)
   deallocate(iwork)
 case default
   call DSYEV('V','L',nmat,eigvec,nmat,eigval,work,lwork,info)
 end select

 deallocate(work)

 if( info /= 0 ) call die('diagonalize_dp: diago failure 2')

end subroutine diagonalize_dp


!=========================================================================
subroutine diagonalize_sp(flavor,matrix,eigval,eigvec)
 implicit none
 character(len=1),intent(in) :: flavor
 real(sp),intent(in)         :: matrix(:,:)
 real(sp),intent(out)        :: eigval(:)
 real(sp),intent(out)        :: eigvec(:,:)
!=====
 integer  :: nmat
 real(sp),allocatable :: work(:)
 integer  :: info
!=====

 nmat = SIZE(matrix,DIM=1)
 allocate(work(3*nmat-1))

 eigvec(:,:) = matrix(:,:)

 call SSYEV('V','L',nmat,eigvec,nmat,eigval,work,3*nmat-1,info)

 deallocate(work)

 if( info /= 0 ) call die('diagonalize_sp: diago failure 1')

end subroutine diagonalize_sp


!=========================================================================
subroutine diagonalize_inplace_dp(flavor,matrix,eigval)
 implicit none
 character(len=1),intent(in) :: flavor
 real(dp),intent(inout)      :: matrix(:,:)
 real(dp),intent(out)        :: eigval(:)
!=====
 integer              :: nmat
 real(dp),allocatable :: work(:)
 integer              :: lwork,info
 integer             :: liwork
 integer,allocatable :: iwork(:)
 integer,allocatable :: isuppz(:)
 real(dp),allocatable    :: eigvec(:,:)
 real(dp)                :: ABSTOL
 real(dp),external       :: DLAMCH
 integer                 :: neig
!=====

 nmat = SIZE(matrix,DIM=1)

 lwork = -1
 allocate(work(1))

 select case(flavor)
 case('r','R')
   allocate(iwork(1))
   allocate(eigvec(nmat,nmat))
   ABSTOL=DLAMCH('S')
   allocate(isuppz(2*nmat))
   call DSYEVR('V','A','L',nmat,matrix,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,iwork,liwork,info)
   liwork = iwork(1)
   deallocate(iwork)
 case('d','D')
   allocate(iwork(1))
   call DSYEVD('V','L',nmat,matrix,nmat,eigval,work,lwork,iwork,liwork,info)
   liwork = iwork(1)
   deallocate(iwork)
 case default
   call DSYEV('V','L',nmat,matrix,nmat,eigval,work,lwork,info)
 end select

 lwork = NINT(work(1))
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_inplace_dp: diago failure 1')

 allocate(work(lwork))

 select case(flavor)
 case('r','R')
   allocate(iwork(liwork))
   call DSYEVR('V','A','L',nmat,matrix,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,iwork,liwork,info)
   matrix(:,:) = eigvec(:,:)
   deallocate(eigvec)
   deallocate(iwork,isuppz)
 case('d','D')
   allocate(iwork(liwork))
   call DSYEVD('V','L',nmat,matrix,nmat,eigval,work,lwork,iwork,liwork,info)
   deallocate(iwork)
 case default
   call DSYEV('V','L',nmat,matrix,nmat,eigval,work,lwork,info)
 end select
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_inplace_dp: diago failure 2')

end subroutine diagonalize_inplace_dp


!=========================================================================
subroutine diagonalize_inplace_sp(flavor,matrix,eigval)
 implicit none
 character(len=1),intent(in) :: flavor
 real(sp),intent(inout) :: matrix(:,:)
 real(sp),intent(out)   :: eigval(:)
!=====
 integer              :: nmat
 real(dp),allocatable :: work(:)
 integer              :: lwork,info
!=====

 nmat = SIZE(matrix,DIM=1)

 lwork = -1
 allocate(work(1))
 call SSYEV('V','L',nmat,matrix,nmat,eigval,work,lwork,info)
 lwork = NINT(work(1))
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_inplace_sp: diago failure 1')

 allocate(work(lwork))
 call SSYEV('V','L',nmat,matrix,nmat,eigval,work,lwork,info)
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_inplace_sp: diago failure 2')

end subroutine diagonalize_inplace_sp


!=========================================================================
subroutine diagonalize_cdp(flavor,matrix,eigval,eigvec)
 implicit none
 character(len=1),intent(in) :: flavor
 complex(dp),intent(in)      :: matrix(:,:)
 real(dp),intent(out)        :: eigval(:)
 complex(dp),intent(out)     :: eigvec(:,:)
!=====
 integer                 :: nmat
 complex(dp),allocatable :: work(:)
 real(dp),allocatable    :: rwork(:)
 integer                 :: lwork,lrwork,info
 integer                 :: liwork
 integer,allocatable     :: iwork(:)
 integer,allocatable     :: isuppz(:)
 complex(dp),allocatable :: matrix_tmp(:,:)
 real(dp)                :: ABSTOL
 real(dp),external       :: DLAMCH
 integer                 :: neig
!=====

 nmat = SIZE(matrix,DIM=1)
 eigvec(:,:) = matrix(:,:)
 lwork = -1

 allocate(work(1))
 allocate(rwork(1))

 select case(flavor)
 case('r','R')
   allocate(iwork(1))
   allocate(matrix_tmp(nmat,nmat))
   matrix_tmp(:,:) = matrix(:,:)
   ABSTOL = DLAMCH('S')
   allocate(isuppz(2*nmat))
   call ZHEEVR('V','A','L',nmat,matrix_tmp,nmat,0.d0,0.d0,1,1,ABSTOL,neig, &
               eigval,eigvec,nmat,isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
   liwork = iwork(1)
   deallocate(iwork)
   lrwork = NINT(rwork(1))
 case('d','D')
   allocate(iwork(1))
   call ZHEEVD('V','L',nmat,eigvec,nmat,eigval,work,lwork,rwork,lrwork,iwork,liwork,info)
   liwork = iwork(1)
   deallocate(iwork)
   lrwork = NINT(rwork(1))
 case default
   call ZHEEV('V','L',nmat,eigvec,nmat,eigval,work,lwork,rwork,info)
   lrwork = 3 * nmat - 2
 end select

 lwork = NINT(REAL(work(1),dp))
 deallocate(rwork)
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_cdp: diago failure 1')

 allocate(work(lwork))
 allocate(rwork(lrwork))

 select case(flavor)
 case('r','R')
   allocate(iwork(1))
   call ZHEEVR('V','A','L',nmat,matrix_tmp,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
   deallocate(iwork,isuppz)
   deallocate(matrix_tmp)
 case('d','D')
   allocate(iwork(liwork))
   call ZHEEVD('V','L',nmat,eigvec,nmat,eigval,work,lwork,rwork,lrwork,iwork,liwork,info)
   deallocate(iwork)
 case default
   call ZHEEV('V','L',nmat,eigvec,nmat,eigval,work,lwork,rwork,info)
 end select

 deallocate(work)
 deallocate(rwork)

 if( info /= 0 ) call die('diagonalize_cdp: diago failure 2')

end subroutine diagonalize_cdp


!=========================================================================
subroutine diagonalize_inplace_cdp(flavor,matrix,eigval)
 implicit none
 character(len=1),intent(in) :: flavor
 complex(dp),intent(inout) :: matrix(:,:)
 real(dp),intent(out)      :: eigval(:)
!=====
 integer                 :: nmat
 complex(dp),allocatable :: work(:)
 real(dp),allocatable    :: rwork(:)
 integer                 :: lwork,lrwork,info
 integer                 :: liwork
 integer,allocatable     :: iwork(:)
 integer,allocatable     :: isuppz(:)
 complex(dp),allocatable :: eigvec(:,:)
 real(dp)                :: ABSTOL
 real(dp),external       :: DLAMCH
 integer                 :: neig
!=====

 nmat = SIZE(matrix,DIM=1)

 lwork = -1
 allocate(work(1))
 allocate(rwork(1))

 select case(flavor)
 case('r','R')
   allocate(iwork(1))
   allocate(eigvec(nmat,nmat))
   allocate(isuppz(2*nmat))
   ABSTOL=DLAMCH('S')
   call ZHEEVR('V','A','L',nmat,matrix,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
   liwork = iwork(1)
   deallocate(iwork)
   lrwork = NINT(rwork(1))
 case('d','D')
   allocate(iwork(1))
   call ZHEEVD('V','L',nmat,matrix,nmat,eigval,work,lwork,rwork,lrwork,iwork,liwork,info)
   liwork = iwork(1)
   deallocate(iwork)
   lrwork = NINT(rwork(1))
 case default
   call ZHEEV('V','L',nmat,matrix,nmat,eigval,work,lwork,rwork,info)
   lrwork = 3 * nmat - 2
 end select

 lwork = NINT(REAL(work(1),dp))
 deallocate(rwork)
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_inplace_cdp: diago failure 1')

 allocate(work(lwork))
 allocate(rwork(lrwork))

 select case(flavor)
 case('r','R')
   allocate(iwork(1))
   call ZHEEVR('V','A','L',nmat,matrix,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
   matrix(:,:) = eigvec(:,:)
   deallocate(eigvec)
   deallocate(iwork,isuppz)
 case('d','D')
   allocate(iwork(liwork))
   call ZHEEVD('V','L',nmat,matrix,nmat,eigval,work,lwork,rwork,lrwork,iwork,liwork,info)
   deallocate(iwork)
 case default
   call ZHEEV('V','L',nmat,matrix,nmat,eigval,work,lwork,rwork,info)
 end select

 deallocate(work)
 deallocate(rwork)

 if( info /= 0 ) call die('diagonalize_inplace_cdp: diago failure 2')

end subroutine diagonalize_inplace_cdp


!=========================================================================
subroutine diagonalize_davidson(tolerance,nstep,ham,neig,eigval,eigvec)
 implicit none

 real(dp),intent(in)  :: tolerance
 integer,intent(inout) :: nstep
 real(dp),intent(in)  :: ham(:,:)
 integer,intent(in)   :: neig
 real(dp),intent(out) :: eigval(:)
 real(dp),intent(out) :: eigvec(:,:)
!=====
 integer              :: nmat,imat
 integer              :: mm,mm_max
 integer              :: ieig,icycle
 real(dp),allocatable :: bb(:,:),atilde(:,:),ab(:,:),qq(:,:)
 real(dp),allocatable :: lambda(:),alphavec(:,:)
 real(dp)             :: residual_norm
!=====

 write(stdout,'(/,1x,a,i5)') 'Davidson diago for eigenvector count: ',neig

 nmat = SIZE(ham(:,:),DIM=1)
 eigval(:) = 0.0_dp


 mm     = neig
 mm_max = mm * nstep
 if( mm_max > nmat ) then
   nstep = nmat / neig
   mm_max = mm * nstep
 endif

 allocate(bb(nmat,mm_max))
 allocate(atilde(mm_max,mm_max))
 allocate(ab(nmat,mm_max))

 allocate(qq(nmat,neig))

 !
 ! Initialize with stupid coefficients
 bb(:,1:neig) = 0.01_dp
 forall(ieig=1:neig)
   bb(ieig,ieig) = 1.0_dp
 end forall
 call orthogonalize(bb(:,1:neig))


 ab(:,1:mm) = MATMUL( ham(:,:) , bb(:,1:mm) )


 do icycle=1,nstep

   mm = icycle * neig
   write(stdout,*) 'icycle mm',icycle,mm


   atilde(1:mm,1:mm) = MATMUL( TRANSPOSE(bb(:,1:mm)) , ab(:,1:mm) )


   allocate(lambda(mm),alphavec(mm,mm))
   call diagonalize(' ',atilde(1:mm,1:mm),lambda,alphavec)

   write(stdout,*) 'icycle',icycle,lambda(1:mm)

   residual_norm = 0.0_dp
   do ieig=1,neig

     qq(:,ieig) = MATMUL( ab(:,1:mm) ,  alphavec(1:mm,ieig) ) &
                   - lambda(ieig) * MATMUL( bb(:,1:mm) , alphavec(1:mm,ieig) )

     residual_norm = MAX( residual_norm , NORM2(qq(:,ieig)) )
   enddo

   write(stdout,'(1x,a,1x,i4,1x,es12.4)') 'Max residual norm for cycle: ',icycle,residual_norm

   !
   ! Convergence reached... or not
   if( icycle == nstep .OR. residual_norm < tolerance ) then
     eigval(1:neig) = lambda(1:neig)
     eigvec(:,1:neig) = MATMUL( bb(:,1:mm) , alphavec(1:mm,1:neig) )
     deallocate(lambda,alphavec)
     exit
   endif


   !
   ! New trial vectors
   !
   forall(imat=1:nmat,ieig=1:neig)
     bb(imat,mm+ieig) = qq(imat,ieig) / ( lambda(ieig) - ham(imat,imat) )
   end forall
   call orthogonalize(bb(:,:mm+neig))



   ab(:,mm+1:mm+neig) = MATMUL( ham(:,:) , bb(:,mm+1:mm+neig) )


   deallocate(lambda,alphavec)

 enddo ! icycle

 deallocate(ab,atilde)

end subroutine diagonalize_davidson


!=========================================================================
subroutine orthogonalize(vec)
 implicit none
!=====
 real(dp),intent(inout) :: vec(:,:)
!=====
 integer :: ivec,jvec,nvec
!=====

 nvec = SIZE(vec(:,:),DIM=2)

 do ivec=1,nvec
   ! Orthogonalize to previous vectors
   do jvec=1,ivec-1
     vec(:,ivec) = vec(:,ivec) - vec(:,jvec) * DOT_PRODUCT( vec(:,ivec) , vec(:,jvec) ) / DOT_PRODUCT( vec(:,jvec) , vec(:,jvec) )
   enddo
   ! Normalize
   vec(:,ivec) = vec(:,ivec) / NORM2( vec(:,ivec) )
 enddo

end subroutine orthogonalize


!=========================================================================
subroutine check_unitarity(cmat)
 implicit none
 complex(dp),intent(in) :: cmat(:,:)
!=====
 real(dp),parameter :: tol=1.0e-9_dp
 integer :: nmat
 integer :: imat,jmat
 complex(dp),allocatable :: cmat_tmp(:,:)
!=====

 nmat = SIZE(cmat,DIM=1)
 allocate(cmat_tmp,MOLD=cmat)

 cmat_tmp = MATMUL( cmat , TRANSPOSE(CONJG(cmat)) )
 do imat=1,nmat
   do jmat=1,nmat
     if(imat==jmat) then
      if(ABS(cmat_tmp(imat,jmat)-1.0_dp)>tol) then
        write(stdout,*) imat,jmat,cmat_tmp(imat,jmat)
        call die('MATRIX IS NOT UNITARY/ORTHOGONAL')
      endif
     else
      if(ABS(cmat_tmp(imat,jmat))>tol) then
        write(stdout,*) imat,jmat,cmat_tmp(imat,jmat)
        call die('MATRIX IS NOT UNITARY/ORTHOGONAL')
      endif
     endif
   enddo
 enddo
 cmat_tmp = MATMUL( TRANSPOSE(CONJG(cmat)) , cmat )
 do imat=1,nmat
   do jmat=1,nmat
     if(imat==jmat) then
      if(ABS(cmat_tmp(imat,jmat)-1.0_dp)>tol) then
        write(stdout,*) imat,jmat,cmat_tmp(imat,jmat)
        call die('MATRIX IS NOT UNITARY/ORTHOGONAL')
      endif
     else
      if(ABS(cmat_tmp(imat,jmat))>tol) then
        write(stdout,*) imat,jmat,cmat_tmp(imat,jmat)
        call die('MATRIX IS NOT UNITARY/ORTHOGONAL')
      endif
     endif
   enddo
 enddo

end subroutine check_unitarity


!=========================================================================
subroutine cross_product(u1,u2,u3)
 implicit none
 real(dp),intent(in)  :: u1(3),u2(3)
 real(dp),intent(out) :: u3(3)
!=====
!=====

 u3(1) = u1(2) * u2(3) - u1(3) * u2(2)
 u3(2) = u1(3) * u2(1) - u1(1) * u2(3)
 u3(3) = u1(1) * u2(2) - u1(2) * u2(1)

end subroutine cross_product


!=========================================================================
pure function determinant_3x3_matrix(mat) RESULT(det)
 implicit none

 real(dp) :: det
 real(dp),intent(in) :: mat(3,3)
!=====
!=====

 det =  mat(1,1) * mat(2,2) * mat(3,3)  &
      + mat(1,2) * mat(2,3) * mat(3,1)  &
      + mat(1,3) * mat(2,1) * mat(3,2)  &
      - mat(1,3) * mat(2,2) * mat(3,1)  &
      - mat(1,2) * mat(2,1) * mat(3,3)  &
      - mat(1,1) * mat(2,3) * mat(3,2)

end function determinant_3x3_matrix


!=========================================================================
end module m_linear_algebra


!=========================================================================
