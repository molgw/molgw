subroutine RPAx(TDA,doACFDT,exchange_kernel,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ENuc,ERHF, & 
                ERI,dipole_int,eHF)

! Perform random phase approximation calculation with exchange (aka TDHF)

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: TDA
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin
  double precision,allocatable  :: Omega(:,:)
  double precision,allocatable  :: XpY(:,:,:)
  double precision,allocatable  :: XmY(:,:,:)

  double precision              :: EcRPAx(nspin)
  double precision              :: EcAC(nspin)

! Hello world

  write(*,*)
  write(*,*)'***********************************************************'
  write(*,*)'|  Random phase approximation calculation with exchange   |'
  write(*,*)'***********************************************************'
  write(*,*)

! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*) ' => RPAx + TDA = CIS '
    write(*,*)
  end if

! Initialization

  EcRPAx(:) = 0d0
  EcAC(:)   = 0d0

! Memory allocation

  allocate(Omega(nS,nspin),XpY(nS,nS,nspin),XmY(nS,nS,nspin))

! Singlet manifold

  if(singlet) then 

    ispin = 1

    call linear_response(ispin,.false.,TDA,eta,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI, &
                         EcRPAx(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call print_excitation('RPAx@HF     ',ispin,nS,Omega(:,ispin))
    call print_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

  endif

! Triplet manifold 

  if(triplet) then 

    ispin = 2

    call linear_response(ispin,.false.,TDA,eta,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI, &
                         EcRPAx(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call print_excitation('RPAx@HF     ',ispin,nS,Omega(:,ispin))
    call print_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

  endif

! if(exchange_kernel) then

    EcRPAx(1) = 0.5d0*EcRPAx(1)
    EcRPAx(2) = 1.5d0*EcRPAx(2)

! end if

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10)') 'Tr@RPAx correlation energy (singlet) =',EcRPAx(1)
  write(*,'(2X,A50,F20.10)') 'Tr@RPAx correlation energy (triplet) =',EcRPAx(2)
  write(*,'(2X,A50,F20.10)') 'Tr@RPAx correlation energy           =',EcRPAx(1) + EcRPAx(2)
  write(*,'(2X,A50,F20.10)') 'Tr@RPAx total energy                 =',ENuc + ERHF + EcRPAx(1) + EcRPAx(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

  if(doACFDT) then

    write(*,*) '-------------------------------------------------------'
    write(*,*) 'Adiabatic connection version of RPAx correlation energy'
    write(*,*) '-------------------------------------------------------'
    write(*,*)

    call ACFDT(exchange_kernel,.false.,.false.,.false.,TDA,.false.,singlet,triplet,eta, &
               nBas,nC,nO,nV,nR,nS,ERI,eHF,eHF,EcAC)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'AC@RPAx correlation energy (singlet) =',EcAC(1)
    write(*,'(2X,A50,F20.10)') 'AC@RPAx correlation energy (triplet) =',EcAC(2)
    write(*,'(2X,A50,F20.10)') 'AC@RPAx correlation energy           =',EcAC(1) + EcAC(2)
    write(*,'(2X,A50,F20.10)') 'AC@RPAx total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine RPAx
