subroutine linear_response(ispin,dRPA,TDA,eta,nBas,nC,nO,nV,nR,nS,lambda,e,ERI,Ec,Omega,XpY,XmY)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dRPA
  logical,intent(in)            :: TDA
  double precision,intent(in)   :: eta
  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

  ! Local variables

  integer :: i,j,k

  double precision              :: trace_matrix
  double precision,allocatable  :: A(:,:)
  double precision,allocatable  :: B(:,:)
  double precision,allocatable  :: ApB(:,:)
  double precision,allocatable  :: AmB(:,:)
  double precision,allocatable  :: AmBSq(:,:)
  double precision,allocatable  :: AmBIv(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: tmp(:,:)

! Output variables

  double precision,intent(out)  :: Ec
  double precision,intent(out)  :: Omega(nS)
  double precision,intent(out)  :: XpY(nS,nS)
  double precision,intent(out)  :: XmY(nS,nS)

! Memory allocation

  allocate(A(nS,nS),B(nS,nS),ApB(nS,nS),AmB(nS,nS),AmBSq(nS,nS),AmBIv(nS,nS),Z(nS,nS),tmp(nS,nS))

! Build A and B matrices 

  call linear_response_A_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,e,ERI,A)

! Tamm-Dancoff approximation

  if(TDA) then
 
    B(:,:)   = 0d0
    XpY(:,:) = A(:,:)
    call diagonalize_matrix(nS,XpY,Omega)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

  else

    call linear_response_B_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,ERI,B)

    ! Build A + B and A - B matrices 

    ApB = A + B
    AmB = A - B

  ! Diagonalize linear response matrix

    call diagonalize_matrix(nS,AmB,Omega)


!   do ia=1,nS
!     if(Omega(ia) < 0d0) Omega(ia) = 0d0
!   end do

    call ADAt(nS,AmB,1d0*dsqrt(Omega),AmBSq)
    call ADAt(nS,AmB,1d0/dsqrt(Omega),AmBIv)

    call dgemm('N','N',nS,nS,nS,1d0,ApB,size(ApB,1),AmBSq,size(AmBSq,1),0d0,tmp,size(tmp,1))
    call dgemm('N','N',nS,nS,nS,1d0,AmBSq,size(AmBSq,1),tmp,size(tmp,1),0d0,Z,size(Z,1))

   call diagonalize_matrix(nS,Z,Omega)

 
  ! do ia=1,nS
  !   if(Omega(ia) < 0d0) Omega(ia) = 0d0
  ! end do

    Omega = dsqrt(Omega)

    call dgemm('T','N',nS,nS,nS,1d0,Z,size(Z,1),AmBSq,size(AmBSq,1),0d0,XpY,size(XpY,1))
    call DA(nS,1d0/dsqrt(Omega),XpY)

    call dgemm('T','N',nS,nS,nS,1d0,Z,size(Z,1),AmBIv,size(AmBIv,1),0d0,XmY,size(XmY,1))
    call DA(nS,1d0*dsqrt(Omega),XmY)
 
  end if

  ! Compute the RPA correlation energy

    Ec = 0.5d0*(sum(Omega) - trace_matrix(nS,A))

end subroutine linear_response
