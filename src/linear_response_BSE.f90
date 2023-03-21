subroutine linear_response_BSE(ispin,dRPA,TDA,BSE,eta,nBas,nC,nO,nV,nR,nS,lambda,e,ERI,A_BSE,B_BSE,Ec,Omega,XpY,XmY)

! Compute linear response with BSE additional terms

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: dRPA
  logical,intent(in)            :: TDA
  logical,intent(in)            :: BSE
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: A_BSE(nS,nS)
  double precision,intent(in)   :: B_BSE(nS,nS)

! Local variables

  double precision              :: trace_matrix
  double precision,allocatable  :: A(:,:)
  double precision,allocatable  :: B(:,:)
  double precision,allocatable  :: ApB(:,:)
  double precision,allocatable  :: AmB(:,:)
  double precision,allocatable  :: AmBSq(:,:)
  double precision,allocatable  :: AmBIv(:,:)
  double precision,allocatable  :: Z(:,:)

! Output variables

  double precision,intent(out)  :: Ec
  double precision,intent(out)  :: Omega(nS)
  double precision,intent(out)  :: XpY(nS,nS)
  double precision,intent(out)  :: XmY(nS,nS)

! Memory allocation

  allocate(A(nS,nS),B(nS,nS),ApB(nS,nS),AmB(nS,nS),AmBSq(nS,nS),AmBIv(nS,nS),Z(nS,nS))

! Build A and B matrices 

  call linear_response_A_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,e,ERI,A)

  if(BSE) A(:,:) = A(:,:) - A_BSE(:,:)

! Tamm-Dancoff approximation

  if(TDA) then
 
    B(:,:)   = 0d0
    XpY(:,:) = A(:,:)
    call diagonalize_matrix(nS,XpY,Omega)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

  else

    call linear_response_B_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,ERI,B)

    if(BSE) B(:,:) = B(:,:) - B_BSE(:,:)

    ! Build A + B and A - B matrices 

    ApB = A + B
    AmB = A - B

  ! Diagonalize linear response matrix

    call diagonalize_matrix(nS,AmB,Omega)


    call ADAt(nS,AmB,1d0*sqrt(Omega),AmBSq)
    call ADAt(nS,AmB,1d0/sqrt(Omega),AmBIv)

    Z = matmul(AmBSq,matmul(ApB,AmBSq))

   call diagonalize_matrix(nS,Z,Omega)


    Omega = sqrt(Omega)

    XpY = matmul(transpose(Z),AmBSq)
    call DA(nS,1d0/sqrt(Omega),XpY)

    XmY = matmul(transpose(Z),AmBIv)
    call DA(nS,1d0*sqrt(Omega),XmY)
 
  end if

  ! Compute the RPA correlation energy

    Ec = 0.5d0*(sum(Omega) - trace_matrix(nS,A))

end subroutine linear_response_BSE
