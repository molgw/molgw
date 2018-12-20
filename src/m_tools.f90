!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! some unsorted basic numerical algorithms
!
!=========================================================================
module m_tools
 use m_definitions
 use m_warning,only: die

 integer,save :: idum

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

 interface append_to_list
   module procedure append_to_list_i
   module procedure append_to_list_r
 end interface


contains


!=========================================================================
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
subroutine init_seed(iseed)
 implicit none
 integer,intent(in),optional :: iseed
!=====
 integer :: itmp,jtmp
!=====

 if(PRESENT(iseed)) then
  idum=iseed
 else
  call system_clock(idum,itmp,jtmp)
 endif

 write(stdout,'(a,1x,i12)') 'Random seed set to',idum


end subroutine


!=========================================================================
function random()
 implicit none
 real(dp) :: random

 call random_number(random)

end function random


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
subroutine diagonalize_wo_vectors_dp(matrix,eigval)
 implicit none
 real(dp),intent(inout) :: matrix(:,:)
 real(dp),intent(out) :: eigval(:)
!=====
 integer :: nmat
 integer :: info
 real(dp),allocatable :: work(:)
! integer  :: iwork(5*n),ifail(n)
! real(dp) :: z(1,n)
! real(dp) :: work(8*n)
!=====

 nmat = SIZE(matrix,DIM=1)
 allocate(work(3*nmat-1))

 call DSYEV('N','U',nmat,matrix,nmat,eigval,work,3*nmat-1,info)

 deallocate(work)

! call DSYEVX('N','A','U',nmat,matrix,nmat,0.0_dp,0.0_dp,0,0,&
!                         1.0e-20_dp,nmat,eigval,z,1,work,8*nmat,iwork,&
!                         ifail,info)

end subroutine diagonalize_wo_vectors_dp


!=========================================================================
subroutine diagonalize_dp(matrix,eigval,eigvec)
 implicit none
 real(dp),intent(in)  :: matrix(:,:)
 real(dp),intent(out) :: eigval(:)
 real(dp),intent(out) :: eigvec(:,:)
!=====
 integer :: nmat,lwork
 real(dp),allocatable :: work(:)
 integer :: info
#if !defined(LAPACK_DIAGO_FLAVOR_)
 integer             :: liwork
 integer,allocatable :: iwork(:)
#endif
#if !defined(LAPACK_DIAGO_FLAVOR_) && !defined(LAPACK_DIAGO_FLAVOR_D)
 integer             :: isuppz(2*SIZE(matrix,DIM=1))
 real(dp),allocatable:: matrix_tmp(:,:)
 real(dp)            :: ABSTOL
 real(dp),external   :: DLAMCH
 integer             :: neig
#endif
!=====

 nmat = SIZE(matrix,DIM=1)

 eigvec(:,:) = matrix(:,:)

 lwork = -1
 allocate(work(1))
#if defined(LAPACK_DIAGO_FLAVOR_)
 call DSYEV('V','U',nmat,eigvec,nmat,eigval,work,lwork,info)
#elif defined(LAPACK_DIAGO_FLAVOR_D)
 allocate(iwork(1))
 call DSYEVD('V','U',nmat,eigvec,nmat,eigval,work,lwork,iwork,liwork,info)
 liwork = iwork(1)
 deallocate(iwork)
#else
 allocate(iwork(1))
 allocate(matrix_tmp(nmat,nmat))
 ABSTOL = DLAMCH('S')
 matrix_tmp(:,:) = matrix(:,:)
 call DSYEVR('V','A','U',nmat,matrix_tmp,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,iwork,liwork,info)
 liwork = iwork(1)
 deallocate(iwork)
#endif
 lwork = NINT(work(1))
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_dp: diago failure 1')

 allocate(work(lwork))
#if defined(LAPACK_DIAGO_FLAVOR_)
 call DSYEV('V','U',nmat,eigvec,nmat,eigval,work,lwork,info)
#elif defined(LAPACK_DIAGO_FLAVOR_D)
 allocate(iwork(liwork))
 call DSYEVD('V','U',nmat,eigvec,nmat,eigval,work,lwork,iwork,liwork,info)
 deallocate(iwork)
#else
 allocate(iwork(liwork))
 call DSYEVR('V','A','U',nmat,matrix_tmp,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,iwork,liwork,info)
 deallocate(iwork)
 deallocate(matrix_tmp)
#endif
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_dp: diago failure 2')

end subroutine diagonalize_dp


!=========================================================================
subroutine diagonalize_sp(matrix,eigval,eigvec)
 implicit none
 real(sp),intent(in)  :: matrix(:,:)
 real(sp),intent(out) :: eigval(:)
 real(sp),intent(out) :: eigvec(:,:)
!=====
 integer  :: nmat
 real(sp),allocatable :: work(:)
 integer  :: info
!=====

 nmat = SIZE(matrix,DIM=1)
 allocate(work(3*nmat-1))

 eigvec(:,:) = matrix(:,:)

 call SSYEV('V','U',nmat,eigvec,nmat,eigval,work,3*nmat-1,info)

 deallocate(work)

 if( info /= 0 ) call die('diagonalize_sp: diago failure 1')

end subroutine diagonalize_sp


!=========================================================================
subroutine diagonalize_inplace_dp(matrix,eigval)
 implicit none
 real(dp),intent(inout) :: matrix(:,:)
 real(dp),intent(out)   :: eigval(:)
!=====
 integer              :: nmat
 real(dp),allocatable :: work(:)
 integer              :: lwork,info
#if !defined(LAPACK_DIAGO_FLAVOR_)
 integer             :: liwork
 integer,allocatable :: iwork(:)
#endif
#if !defined(LAPACK_DIAGO_FLAVOR_) && !defined(LAPACK_DIAGO_FLAVOR_D)
 integer                 :: isuppz(2*SIZE(matrix,DIM=1))
 real(dp),allocatable    :: eigvec(:,:)
 real(dp)                :: ABSTOL
 real(dp),external       :: DLAMCH
 integer                 :: neig
#endif
!=====

 nmat = SIZE(matrix,DIM=1)

 lwork = -1
 allocate(work(1))
#if defined(LAPACK_DIAGO_FLAVOR_)
 call DSYEV('V','U',nmat,matrix,nmat,eigval,work,lwork,info)
#elif defined(LAPACK_DIAGO_FLAVOR_D)
 allocate(iwork(1))
 call DSYEVD('V','U',nmat,matrix,nmat,eigval,work,lwork,iwork,liwork,info)
 liwork = iwork(1)
 deallocate(iwork)
#else
 allocate(iwork(1))
 allocate(eigvec(nmat,nmat))
 ABSTOL=DLAMCH('S')
 call DSYEVR('V','A','U',nmat,matrix,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,iwork,liwork,info)
 liwork = iwork(1)
 deallocate(iwork)
#endif
 lwork = NINT(work(1))
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_inplace_dp: diago failure 1')

 allocate(work(lwork))
#if defined(LAPACK_DIAGO_FLAVOR_)
 call DSYEV('V','U',nmat,matrix,nmat,eigval,work,lwork,info)
#elif defined(LAPACK_DIAGO_FLAVOR_D)
 allocate(iwork(liwork))
 call DSYEVD('V','U',nmat,matrix,nmat,eigval,work,lwork,iwork,liwork,info)
 deallocate(iwork)
#else
 allocate(iwork(liwork))
 call DSYEVR('V','A','U',nmat,matrix,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,iwork,liwork,info)
 matrix(:,:) = eigvec(:,:)
 deallocate(eigvec)
 deallocate(iwork)
#endif
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_inplace_dp: diago failure 2')

end subroutine diagonalize_inplace_dp


!=========================================================================
subroutine diagonalize_inplace_sp(matrix,eigval)
 implicit none
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
 call SSYEV('V','U',nmat,matrix,nmat,eigval,work,lwork,info)
 lwork = NINT(work(1))
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_inplace_sp: diago failure 1')

 allocate(work(lwork))
 call SSYEV('V','U',nmat,matrix,nmat,eigval,work,lwork,info)
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_inplace_sp: diago failure 2')

end subroutine diagonalize_inplace_sp


!=========================================================================
subroutine diagonalize_cdp(matrix,eigval,eigvec)
 implicit none
 complex(dp),intent(in)  :: matrix(:,:)
 real(dp),intent(out)    :: eigval(:)
 complex(dp),intent(out) :: eigvec(:,:)
!=====
 integer :: nmat
 complex(dp),allocatable :: work(:)
 real(dp),allocatable    :: rwork(:)
 integer                 :: lwork,lrwork,info
#if !defined(LAPACK_DIAGO_FLAVOR_)
 integer                 :: liwork
 integer,allocatable     :: iwork(:)
#endif
#if !defined(LAPACK_DIAGO_FLAVOR_) && !defined(LAPACK_DIAGO_FLAVOR_D)
 integer                 :: isuppz(2*SIZE(matrix,DIM=1))
 complex(dp),allocatable :: matrix_tmp(:,:)
 real(dp)                :: ABSTOL
 real(dp),external       :: DLAMCH
 integer                 :: neig
#endif
!=====

 nmat = SIZE(matrix,DIM=1)
 eigvec(:,:) = matrix(:,:)
 lwork = -1

 allocate(work(1))
 allocate(rwork(1))
#if defined(LAPACK_DIAGO_FLAVOR_)
 call ZHEEV('V','U',nmat,eigvec,nmat,eigval,work,lwork,rwork,info)
 lrwork = 3 * nmat - 2
#elif defined(LAPACK_DIAGO_FLAVOR_D)
 allocate(iwork(1))
 call ZHEEVD('V','U',nmat,eigvec,nmat,eigval,work,lwork,rwork,lrwork,iwork,liwork,info)
 liwork = iwork(1)
 deallocate(iwork)
 lrwork = NINT(rwork(1))
#else
 allocate(iwork(1))
 allocate(matrix_tmp(nmat,nmat))
 matrix_tmp(:,:) = matrix(:,:)
 ABSTOL = DLAMCH('S')
 call ZHEEVR('V','A','U',nmat,matrix_tmp,nmat,0.d0,0.d0,1,1,ABSTOL,neig, &
             eigval,eigvec,nmat,isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
 liwork = iwork(1)
 deallocate(iwork)
 lrwork = NINT(rwork(1))
#endif
 lwork = NINT(REAL(work(1),dp))
 deallocate(rwork)
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_cdp: diago failure 1')

 allocate(work(lwork))
 allocate(rwork(lrwork))
#if defined(LAPACK_DIAGO_FLAVOR_)
 call ZHEEV('V','U',nmat,eigvec,nmat,eigval,work,lwork,rwork,info)
#elif defined(LAPACK_DIAGO_FLAVOR_D)
 allocate(iwork(liwork))
 call ZHEEVD('V','U',nmat,eigvec,nmat,eigval,work,lwork,rwork,lrwork,iwork,liwork,info)
 deallocate(iwork)
#else
 allocate(iwork(1))
 call ZHEEVR('V','A','U',nmat,matrix_tmp,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
 deallocate(iwork)
 deallocate(matrix_tmp)
#endif

 deallocate(work)
 deallocate(rwork)

 if( info /= 0 ) call die('diagonalize_cdp: diago failure 2')

end subroutine diagonalize_cdp


!=========================================================================
subroutine diagonalize_inplace_cdp(matrix,eigval)
 implicit none
 complex(dp),intent(inout) :: matrix(:,:)
 real(dp),intent(out)      :: eigval(:)
!=====
 integer                 :: nmat
 complex(dp),allocatable :: work(:)
 real(dp),allocatable    :: rwork(:)
 integer                 :: lwork,lrwork,info
#if !defined(LAPACK_DIAGO_FLAVOR_)
 integer                 :: liwork
 integer,allocatable     :: iwork(:)
#endif
#if !defined(LAPACK_DIAGO_FLAVOR_) && !defined(LAPACK_DIAGO_FLAVOR_D)
 integer                 :: isuppz(2*SIZE(matrix,DIM=1))
 complex(dp),allocatable :: eigvec(:,:)
 real(dp)                :: ABSTOL
 real(dp),external       :: DLAMCH
 integer                 :: neig
#endif
!=====

 nmat = SIZE(matrix,DIM=1)

 lwork = -1
 allocate(work(1))
 allocate(rwork(1))
#if defined(LAPACK_DIAGO_FLAVOR_)
 call ZHEEV('V','U',nmat,matrix,nmat,eigval,work,lwork,rwork,info)
 lrwork = 3 * nmat - 2
#elif defined(LAPACK_DIAGO_FLAVOR_D)
 allocate(iwork(1))
 call ZHEEVD('V','U',nmat,matrix,nmat,eigval,work,lwork,rwork,lrwork,iwork,liwork,info)
 liwork = iwork(1)
 deallocate(iwork)
 lrwork = NINT(rwork(1))
#else
 allocate(iwork(1))
 allocate(eigvec(nmat,nmat))
 ABSTOL=DLAMCH('S')
 call ZHEEVR('V','A','U',nmat,matrix,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
 liwork = iwork(1)
 deallocate(iwork)
 lrwork = NINT(rwork(1))
#endif
 lwork = NINT(REAL(work(1),dp))
 deallocate(rwork)
 deallocate(work)

 if( info /= 0 ) call die('diagonalize_inplace_cdp: diago failure 1')

 allocate(work(lwork))
 allocate(rwork(lrwork))
#if defined(LAPACK_DIAGO_FLAVOR_)
 call ZHEEV('V','U',nmat,matrix,nmat,eigval,work,lwork,rwork,info)
#elif defined(LAPACK_DIAGO_FLAVOR_D)
 allocate(iwork(liwork))
 call ZHEEVD('V','U',nmat,matrix,nmat,eigval,work,lwork,rwork,lrwork,iwork,liwork,info)
 deallocate(iwork)
#else
 allocate(iwork(1))
 call ZHEEVR('V','A','U',nmat,matrix,nmat,0.d0,0.d0,1,1,ABSTOL,neig,eigval,eigvec,nmat,isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
 matrix(:,:) = eigvec(:,:)
 deallocate(eigvec)
 deallocate(iwork)
#endif

 deallocate(work)
 deallocate(rwork)

 if( info /= 0 ) call die('diagonalize_inplace_cdp: diago failure 2')

end subroutine diagonalize_inplace_cdp


!=========================================================================
!
! Generalized eigenvalue problem
!
!=========================================================================
subroutine diagonalize_generalized_sym(n,matrix,overlap,eigval,eigvec)
 implicit none
 integer,intent(in) :: n
 real(dp),intent(in) :: matrix(n,n),overlap(n,n)
 real(dp),intent(out) :: eigval(n)
 real(dp),intent(out) :: eigvec(n,n)
!=====
 real(dp) :: work(3*n-1),tmp(n,n)
 integer :: info
!=====

 ! write(stdout,*) 'diagonalize_generalized_sym: Enter'
 eigvec(:,:) = matrix(:,:)
 tmp(:,:) = overlap(:,:)

 ! A*x = lambda * B * x
 call DSYGV(1,'V','U',n,eigvec,n,tmp,n,eigval,work,3*n-1,info)
 if(info/=0) call die('ERROR in the symmetric generalized eigenvalue problem')
! write(stdout,*) 'optimal lwork',REAL(work(1))


 ! write(stdout,*) 'diagonalize_generalized_sym: Exit'
end subroutine diagonalize_generalized_sym


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
   call diagonalize(atilde(1:mm,1:mm),lambda,alphavec)

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
subroutine coeffs_gausslegint(xmin,xmax,x,weights,n)
!
! Compute the coefficients (supports and weights)
! for Gauss-Legendre integration.
! Inspired by a routine due to G. Rybicki.
!
! Input:
! xmin=lower bound of integration
! xmax=upper bound of integration
! n=order of integration
!
! Output:
! x(n)=array of support points
! weights(n)=array of integration weights
!

 implicit none

 integer,intent(in) :: n
 real(dp),intent(in) :: xmin,xmax
 real(dp),intent(out) :: x(n),weights(n)
!=====
 real(dp) :: tol,xl,xmean,z,p1,p2,p3,pp,z1
 integer  :: i,j
!=====

 tol=1.d-13

 xl=(xmax-xmin)*0.50_dp
 xmean=(xmax+xmin)*0.50_dp

 do i=1,(n+1)/2
  z = COS(pi*(i-0.250_dp)/(n+0.50_dp))

  do

    p1=1.0_dp
    p2=0.0_dp

    do j=1,n

     p3=p2
     p2=p1
     p1=((2.0_dp*j - 1.0_dp)*z*p2 - (j-1.0_dp)*p3)/j

    enddo

    pp=n*(p2-z*p1)/(1.0_dp - z**2)
    z1=z
    z=z1-p1/pp

    if(abs(z-z1) < tol) exit

  enddo

  x(i)=xmean-xl*z
  x(n+1-i)=xmean+xl*z
  weights(i)=2.0_dp * xl/((1.0_dp-z**2)*pp**2)
  weights(n+1-i)=weights(i)

 enddo

end subroutine coeffs_gausslegint


!=========================================================================
subroutine gequad(p,w)
 implicit none

!Arguments
 real(dp), intent(out) :: p(89),w(89)
!
!       range [10^(-9),1] and accuracy ~10^(-8);
!
  p(1)=4.96142640560223544d19
  p(2)=1.37454269147978052d19
  p(3)=7.58610013441204679d18
  p(4)=4.42040691347806996d18
  p(5)=2.61986077948367892d18
  p(6)=1.56320138155496681d18
  p(7)=9.35645215863028402d17
  p(8)=5.60962910452691703d17
  p(9)=3.3666225119686761d17
  p(10)=2.0218253197947866d17
  p(11)=1.21477756091902017d17
  p(12)=7.3012982513608503d16
  p(13)=4.38951893556421099d16
  p(14)=2.63949482512262325d16
  p(15)=1.58742054072786174d16
  p(16)=9.54806587737665531d15
  p(17)=5.74353712364571709d15
  p(18)=3.455214877389445d15
  p(19)=2.07871658520326804d15
  p(20)=1.25064667315629928d15
  p(21)=7.52469429541933745d14
  p(22)=4.5274603337253175d14
  p(23)=2.72414006900059548d14
  p(24)=1.63912168349216752d14
  p(25)=9.86275802590865738d13
  p(26)=5.93457701624974985d13
  p(27)=3.5709554322296296d13
  p(28)=2.14872890367310454d13
  p(29)=1.29294719957726902d13
  p(30)=7.78003375426361016d12
  p(31)=4.68148199759876704d12
  p(32)=2.8169955024829868d12
  p(33)=1.69507790481958464d12
  p(34)=1.01998486064607581d12
  p(35)=6.13759486539856459d11
  p(36)=3.69320183828682544d11
  p(37)=2.22232783898905102d11
  p(38)=1.33725247623668682d11
  p(39)=8.0467192739036288d10
  p(40)=4.84199582415144143d10
  p(41)=2.91360091170559564d10
  p(42)=1.75321747475309216d10
  p(43)=1.0549735552210995d10
  p(44)=6.34815321079006586d9
  p(45)=3.81991113733594231d9
  p(46)=2.29857747533101109d9
  p(47)=1.38313653595483694d9
  p(48)=8.32282908580025358d8
  p(49)=5.00814519374587467d8
  p(50)=3.01358090773319025d8
  p(51)=1.81337994217503535d8
  p(52)=1.09117589961086823d8
  p(53)=6.56599771718640323d7
  p(54)=3.95099693638497164d7
  p(55)=2.37745694710665991d7
  p(56)=1.43060135285912813d7
  p(57)=8.60844290313506695d6
  p(58)=5.18000974075383424d6
  p(59)=3.116998193057466d6
  p(60)=1.87560993870024029d6
  p(61)=1.12862197183979562d6
  p(62)=679132.441326077231_dp
  p(63)=408658.421279877969_dp
  p(64)=245904.473450669789_dp
  p(65)=147969.568088321005_dp
  p(66)=89038.612357311147_dp
  p(67)=53577.7362552358895_dp
  p(68)=32239.6513926914668_dp
  p(69)=19399.7580852362791_dp
  p(70)=11673.5323603058634_dp
  p(71)=7024.38438577707758_dp
  p(72)=4226.82479307685999_dp
  p(73)=2543.43254175354295_dp
  p(74)=1530.47486269122675_dp
  p(75)=920.941785160749482_dp
  p(76)=554.163803906291646_dp
  p(77)=333.46029740785694_dp
  p(78)=200.6550575335041_dp
  p(79)=120.741366914147284_dp
  p(80)=72.6544243200329916_dp
  p(81)=43.7187810415471025_dp
  p(82)=26.3071631447061043_dp
  p(83)=15.8299486353816329_dp
  p(84)=9.52493152341244004_dp
  p(85)=5.72200417067776041_dp
  p(86)=3.36242234070940928_dp
  p(87)=1.75371394604499472_dp
  p(88)=0.64705932650658966_dp
  p(89)=0.072765905943708247_dp
!
  w(1)=47.67445484528304247d10
  w(2)=11.37485774750442175d9
  w(3)=78.64340976880190239d8
  w(4)=46.27335788759590498d8
  w(5)=24.7380464827152951d8
  w(6)=13.62904116438987719d8
  w(7)=92.79560029045882433d8
  w(8)=52.15931216254660251d8
  w(9)=31.67018011061666244d8
  w(10)=1.29291036801493046d8
  w(11)=1.00139319988015862d8
  w(12)=7.75892350510188341d7
  w(13)=6.01333567950731271d7
  w(14)=4.66141178654796875d7
  w(15)=3.61398903394911448d7
  w(16)=2.80225846672956389d7
  w(17)=2.1730509180930247d7
  w(18)=1.68524482625876965d7
  w(19)=1.30701489345870338d7
  w(20)=1.01371784832269282d7
  w(21)=7.86264116300379329d6
  w(22)=6.09861667912273717d6
  w(23)=4.73045784039455683d6
  w(24)=3.66928949951594161d6
  w(25)=2.8462050836230259d6
  w(26)=2.20777394798527011d6
  w(27)=1.71256191589205524d6
  w(28)=1.32843556197737076d6
  w(29)=1.0304731275955989d6
  w(30)=799345.206572271448_dp
  w(31)=620059.354143595343_dp
  w(32)=480986.704107449333_dp
  w(33)=373107.167700228515_dp
  w(34)=289424.08337412132_dp
  w(35)=224510.248231581788_dp
  w(36)=174155.825690028966_dp
  w(37)=135095.256919654065_dp
  w(38)=104795.442776800312_dp
  w(39)=81291.4458222430418_dp
  w(40)=63059.0493649328682_dp
  w(41)=48915.9040455329689_dp
  w(42)=37944.8484018048756_dp
  w(43)=29434.4290473253969_dp
  w(44)=22832.7622054490044_dp
  w(45)=17711.743950151233_dp
  w(46)=13739.287867104177_dp
  w(47)=10657.7895710752585_dp
  w(48)=8267.42141053961834_dp
  w(49)=6413.17397520136448_dp
  w(50)=4974.80402838654277_dp
  w(51)=3859.03698188553047_dp
  w(52)=2993.51824493299154_dp
  w(53)=2322.1211966811754_dp
  w(54)=1801.30750964719641_dp
  w(55)=1397.30379659817038_dp
  w(56)=1083.91149143250697_dp
  w(57)=840.807939169209188_dp
  w(58)=652.228524366749422_dp
  w(59)=505.944376983506128_dp
  w(60)=392.469362317941064_dp
  w(61)=304.444930257324312_dp
  w(62)=236.162932842453601_dp
  w(63)=183.195466078603525_dp
  w(64)=142.107732186551471_dp
  w(65)=110.23530215723992_dp
  w(66)=85.5113346705382257_dp
  w(67)=66.3325469806696621_dp
  w(68)=51.4552463353841373_dp
  w(69)=39.9146798429449273_dp
  w(70)=30.9624728409162095_dp
  w(71)=24.018098812215013_dp
  w(72)=18.6312338024296588_dp
  w(73)=14.4525541233150501_dp
  w(74)=11.2110836519105938_dp
  w(75)=8.69662175848497178_dp
  w(76)=6.74611236165731961_dp
  w(77)=5.23307018057529994_dp
  w(78)=4.05937850501539556_dp
  w(79)=3.14892659076635714_dp
  w(80)=2.44267408211071604_dp
  w(81)=1.89482240522855261_dp
  w(82)=1.46984505907050079_dp
  w(83)=1.14019261330527007_dp
  w(84)=0.884791217422925293_dp
  w(85)=0.692686387080616483_dp
  w(86)=0.585244576897023282_dp
  w(87)=0.576182522545327589_dp
  w(88)=0.596688817388997178_dp
  w(89)=0.607879901151108771_dp
!
end subroutine gequad


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
pure function capitalize(str)
 implicit none
 character(*), intent(in) :: str
 character(LEN(str))      :: capitalize
!=====
 character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
 character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
 integer :: ic, ii
!=====

 capitalize = str
 do ii=1,LEN_TRIM(str)
   ic = INDEX(low,str(ii:ii))
   if (ic > 0) capitalize(ii:ii) = cap(ic:ic)
 enddo

end function capitalize


!=========================================================================
function orbital_momentum_number(amc)
 character(len=1),intent(in) :: amc
 integer :: orbital_momentum_number
!=====

 select case(capitalize(amc))
 case('S')
   orbital_momentum_number = 0
 case('P')
   orbital_momentum_number = 1
 case('D')
   orbital_momentum_number = 2
 case('F')
   orbital_momentum_number = 3
 case('G')
   orbital_momentum_number = 4
 case('H')
   orbital_momentum_number = 5
 case('I')
   orbital_momentum_number = 6
 case('K')
   orbital_momentum_number = 7
 case default
   write(stdout,*) amc,capitalize(amc)
   call die('orbital_momentum_number: keyword unknown')
 end select


end function orbital_momentum_number


!=========================================================================
pure function orbital_momentum_name(am)
 integer,intent(in) :: am
 character(len=1) :: orbital_momentum_name
!=====

 select case(am)
 case(0)
   orbital_momentum_name='s'
 case(1)
   orbital_momentum_name='p'
 case(2)
   orbital_momentum_name='d'
 case(3)
   orbital_momentum_name='f'
 case(4)
   orbital_momentum_name='g'
 case(5)
   orbital_momentum_name='h'
 case(6)
   orbital_momentum_name='i'
 case(7)
   orbital_momentum_name='k'
 case(8)
   orbital_momentum_name='l'
 case(9)
   orbital_momentum_name='m'
 case(10)
   orbital_momentum_name='n'
 case(11)
   orbital_momentum_name='o'
 case(12)
   orbital_momentum_name='q'
 case default
   orbital_momentum_name='x'
 end select

end function orbital_momentum_name


!=========================================================================
subroutine append_to_list_i(new_element,list)
 implicit none

 integer,intent(in)                :: new_element
 integer,allocatable,intent(inout) :: list(:)
!=====
 integer :: nsize
 integer,allocatable :: list_old(:)
!=====

 if( ALLOCATED(list) ) then
   nsize = SIZE(list)
 else
   nsize = 0
 endif

 ! Copy old list and free the list
 allocate(list_old(nsize))
 if( nsize > 0 ) then
   list_old(1:nsize) =list(1:nsize)
   deallocate(list)
 endif

 allocate(list(nsize+1))
 if( nsize > 0 ) then
   list(1:nsize) = list_old(1:nsize)
 endif
 list(nsize+1) = new_element


end subroutine append_to_list_i


!=========================================================================
subroutine append_to_list_r(new_element,list)
 implicit none

 real(dp),intent(in)                :: new_element
 real(dp),allocatable,intent(inout) :: list(:)
!=====
 integer :: nsize
 real(dp),allocatable :: list_old(:)
!=====

 if( ALLOCATED(list) ) then
   nsize = SIZE(list)
 else
   nsize = 0
 endif

 ! Copy old list and free the list
 allocate(list_old(nsize))
 if( nsize > 0 ) then
   list_old(1:nsize) =list(1:nsize)
   deallocate(list)
 endif

 allocate(list(nsize+1))
 if( nsize > 0 ) then
   list(1:nsize) = list_old(1:nsize)
 endif
 list(nsize+1) = new_element


end subroutine append_to_list_r


!=========================================================================
function pade(z_in,f_in,z_out) RESULT(f_out)
 implicit none

 complex(dp),intent(in) :: z_in(:),f_in(:)
 complex(dp),intent(in) :: z_out
 complex(dp) :: f_out
!=====
 integer     :: ip,np
 complex(dp),allocatable :: aa(:)
 complex(dp),allocatable :: Az(:), Bz(:)
!=====

 np = SIZE(z_in)
 allocate(aa(np),Az(0:np),Bz(0:np))

 call calculate_pade_a(aa,z_in,f_in)

 Az(0) = (0.0_dp,0.0_dp)
 Az(1) = aa(1)
 Bz(0) = (1.0_dp,0.0_dp)
 Bz(1) = (1.0_dp,0.0_dp)

 do ip=1,np-1
   Az(ip+1) = Az(ip) + ( z_out - z_in(ip) ) * aa(ip+1) * Az(ip-1)
   Bz(ip+1) = Bz(ip) + ( z_out - z_in(ip) ) * aa(ip+1) * Bz(ip-1)
 enddo

 f_out = Az(np) / Bz(np)

 deallocate(aa,Az,Bz)

end function pade


!=========================================================================
subroutine calculate_pade_a(aa,z_in,f_in)
 implicit none

 complex(dp),intent(in)  :: z_in(:),f_in(:)
 complex(dp),intent(out) :: aa(:)
!=====
 integer     :: np,ip,jp
 complex(dp),allocatable :: gtmp(:,:)
!=====

 np = SIZE(z_in)
 allocate(gtmp(np,np))

 gtmp(1,1:np) = f_in(1:np)

 do ip=2,np
   do jp=ip,np
     gtmp(ip,jp) = (gtmp(ip-1,ip-1) - gtmp(ip-1,jp)) / ( (z_in(jp) - z_in(ip-1)) * gtmp(ip-1,jp) )
   enddo
 enddo
 do ip=1,np
   aa(ip) = gtmp(ip,ip)
 enddo

 deallocate(gtmp)

end subroutine calculate_pade_a

!=========================================================================
function get_number_of_lines(filename) result(nlines)
 implicit none
 character(len=*),intent(in)  :: filename
 character(len=100)           :: cur_string
 integer                      :: nlines
!=====
 integer   :: file_unit,io
!=====

 nlines=0
 open(newunit=file_unit, file = filename)
 do
   read(file_unit,'(A)',iostat=io)cur_string
   if ( io /= 0 ) exit
   nlines=nlines+1
 end do
 close(file_unit)

end function get_number_of_lines

!=========================================================================
pure function get_number_of_elements(string) result(num)
 implicit none
 character(len=*),intent(in)  :: string
 integer                      :: num
!=====
 integer   :: ip,pos
!=====

 pos = 1
 num = 0

 do
   ip = VERIFY(string(pos:),' ')  !-- Find next non-blank
   if( ip == 0 ) exit             !-- No word found
   num = num + 1                 !-- Found something
   pos = pos + ip - 1             !-- Move to start of the word
   ip = SCAN(string(pos:),' ')    !-- Find next blank
   if( ip == 0 ) exit             !-- No blank found
   pos = pos + ip - 1             !-- Move to the blank
 end do

end function get_number_of_elements


!=========================================================================
subroutine string_to_integers(string_in,iarray)
 implicit none

 character(len=*),intent(in) :: string_in
 integer,intent(inout)       :: iarray(:)
!=====
 character(LEN(string_in)) :: string
 integer                   :: ilen,inextblank,ii
!=====

 string = string_in

 ilen = LEN(TRIM(string))
 ii = 0
 do while( ilen > 0 )
   string = ADJUSTL(string)
   inextblank = INDEX(string,' ')
   ii = ii + 1
   if( ii > SIZE(iarray) ) exit
   read(string(1:inextblank-1),*) iarray(ii)
   string = string(inextblank+1:)
   ilen = LEN(TRIM(string))
 enddo

end subroutine string_to_integers


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
subroutine string_to_reals(string_in,rarray)
 implicit none

 character(len=*),intent(in) :: string_in
 real(dp),intent(inout)      :: rarray(:)
!=====
 character(LEN(string_in)) :: string
 integer            :: ilen,inextblank,ii
!=====

 string = string_in

 ilen = LEN(TRIM(string))
 ii = 0
 do while( ilen > 0 )
   string = ADJUSTL(string)
   inextblank = INDEX(string,' ')
   ii = ii + 1
   if( ii > SIZE(rarray) ) exit
   read(string(1:inextblank-1),*) rarray(ii)
   string = string(inextblank+1:)
   ilen = LEN(TRIM(string))
 enddo

end subroutine string_to_reals


end module m_tools


!=========================================================================
