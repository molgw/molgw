subroutine ACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,ERI,XpY,XmY,EcAC)

! Compute the correlation energy via the adiabatic connection formula

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: exchange_kernel
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: ia,jb,kc
  double precision              :: delta_spin
  double precision              :: delta_Kx
  double precision,allocatable  :: Ap(:,:)
  double precision,allocatable  :: Bp(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: Y(:,:)
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: EcAC

! Singlet or triplet manifold?

  delta_spin = 0d0
  if(ispin == 1) delta_spin = +1d0
  if(ispin == 2) delta_spin = -1d0

! Exchange kernel

  delta_Kx = 0d0
  if(exchange_kernel) delta_Kx = 1d0
  
! Memory allocation

  allocate(Ap(nS,nS),Bp(nS,nS),X(nS,nS),Y(nS,nS))

! Compute Aiajb = (ia|bj) and Biajb = (ia|jb)

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1

            Ap(ia,jb) = (1d0 + delta_spin)*ERI(i,b,a,j) & 
                      - delta_Kx*ERI(i,b,j,a)

            Bp(ia,jb) = (1d0 + delta_spin)*ERI(i,j,a,b) &
                      - delta_Kx*ERI(i,j,b,a)

        enddo
      enddo
    enddo
  enddo

! Compute Tr(K x P_lambda)

  X(:,:) = 0.5d0*(XpY(:,:) + XmY(:,:))
  Y(:,:) = 0.5d0*(XpY(:,:) - XmY(:,:))

  EcAC = trace_matrix(nS,matmul(X,matmul(Bp,transpose(Y))) + matmul(Y,matmul(Bp,transpose(X)))) &
       + trace_matrix(nS,matmul(X,matmul(Ap,transpose(X))) + matmul(Y,matmul(Ap,transpose(Y)))) &
       - trace_matrix(nS,Ap)

end subroutine ACFDT_correlation_energy

