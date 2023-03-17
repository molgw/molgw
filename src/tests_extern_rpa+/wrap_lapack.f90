subroutine diagonalize_general_matrix(N,A,WR,VR)

! Diagonalize a non-symmetric square matrix

  implicit none

! Input variables

  integer :: i,j,k
  integer,intent(in)            :: N
  double precision,intent(inout):: A(N,N)
  double precision,intent(out)  :: VR(N,N)
  double precision,intent(out)  :: WR(N)

! Local variables

  integer                       :: lwork,info
  double precision,allocatable  :: work(:),WI(:),VL(:,:)

! Memory allocation

  allocate(work(1),WI(N),VL(N,N))
  
  lwork = -1
  call dgeev('V','V',N,A,N,WR,WI,VL,N,VR,N,work,lwork,info)
  lwork = int(work(1))

  deallocate(work)

  allocate(work(lwork))

  call dgeev('V','V',N,A,N,WR,WI,VL,N,VR,N,work,lwork,info)

  if(info /= 0) then 
    print*,'Problem in diagonalize_general_matrix (dgeev)!!'
  endif

end subroutine diagonalize_general_matrix

subroutine diagonalize_matrix(N,A,e)

! Diagonalize a square matrix

  implicit none

! Input variables

  integer,intent(in)            :: N
  double precision,intent(inout):: A(N,N)
  double precision,intent(out)  :: e(N)

! Local variables

  integer                       :: lwork,info
  double precision,allocatable  :: work(:)

! Memory allocation

  allocate(work(1))
  
  lwork = -1
  call dsyev('V','U',N,A,N,e,work,lwork,info)
  lwork = int(work(1))

  deallocate(work)

  allocate(work(lwork))

  call dsyev('V','U',N,A,N,e,work,lwork,info)
 
  if(info /= 0) then 
    print*,'Problem in diagonalize_matrix (dsyev)!!'
  endif
  
end subroutine diagonalize_matrix

subroutine svd(N,A,U,D,Vt)

  ! Compute A = U.D.Vt
  ! Dimension of A is NxN

  implicit none

  integer, intent(in)             :: N
  double precision,intent(in)     :: A(N,N)
  double precision,intent(out)    :: U(N,N)
  double precision,intent(out)    :: Vt(N,N)
  double precision,intent(out)    :: D(N)
  double precision,allocatable    :: work(:)
  integer                         :: info,lwork

  double precision,allocatable    :: scr(:,:)

  allocate (scr(N,N))

  scr(:,:) = A(:,:)

  ! Find optimal size for temporary arrays

  allocate(work(1))

  lwork = -1
  call dgesvd('A','A',N,N,scr,N,D,U,N,Vt,N,work,lwork,info)
  lwork = int(work(1))

  deallocate(work)

  allocate(work(lwork))

  call dgesvd('A','A',N,N,scr,N,D,U,N,Vt,N,work,lwork,info)

  deallocate(work,scr)

  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif

end

subroutine inverse_matrix(N,A,B)

! Returns the inverse of the square matrix A in B

  implicit none

  integer,intent(in)             :: N
  double precision, intent(in)   :: A(N,N)
  double precision, intent(out)  :: B(N,N)

  integer                        :: info,lwork
  integer, allocatable           :: ipiv(:)
  double precision,allocatable   :: work(:)

  allocate (ipiv(N),work(N*N))
  lwork = size(work)

  B(1:N,1:N) = A(1:N,1:N)

  call dgetrf(N,N,B,N,ipiv,info)

  if (info /= 0) then

    print*,info
    stop 'error in inverse (dgetrf)!!'

  endif

  call dgetri(N,B,N,ipiv,work,lwork,info)

  if (info /= 0) then

    print *,  info
    stop 'error in inverse (dgetri)!!'

  endif

  deallocate(ipiv,work)

end subroutine inverse_matrix

subroutine linear_solve(N,A,b,x,rcond)

! Solve the linear system A.x = b where A is a NxN matrix
! and x and x are vectors of size N

  implicit none

  integer,intent(in)             :: N
  double precision,intent(out)   :: A(N,N),b(N),rcond
  double precision,intent(out)   :: x(N)

  integer                        :: info,lwork
  double precision               :: ferr,berr
  integer,allocatable            :: ipiv(:),iwork(:)
  double precision,allocatable   :: AF(:,:),work(:)

  lwork = 3*N
  allocate(AF(N,N),ipiv(N),work(lwork),iwork(N))

  call dsysvx('N','U',N,1,A,N,AF,N,ipiv,b,N,x,N,rcond,ferr,berr,work,lwork,iwork,info)

! if (info /= 0) then

!   print *,  info
!   stop 'error in linear_solve (dsysvx)!!'

! endif

end subroutine linear_solve

subroutine easy_linear_solve(N,A,b,x)

! Solve the linear system A.x = b where A is a NxN matrix
! and x and x are vectors of size N

  implicit none

  integer,intent(in)             :: N
  double precision,intent(in)    :: A(N,N),b(N)
  double precision,intent(out)   :: x(N)

  integer                        :: info,lwork
  integer,allocatable            :: ipiv(:)
  double precision,allocatable   :: work(:)

  allocate(ipiv(N),work(N*N))
  lwork = size(work)

  x = b

  call dsysv('U',N,1,A,N,ipiv,x,N,work,lwork,info)

  if (info /= 0) then

    print *,  info
    stop 'error in linear_solve (dsysv)!!'

  endif

end subroutine easy_linear_solve

