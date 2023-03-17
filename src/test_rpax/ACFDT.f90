subroutine ACFDT(exchange_kernel,doXBS,dRPA,TDA_W,TDA,BSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,eW,e,EcAC)

! Compute the correlation energy via the adiabatic connection fluctuation dissipation theorem

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: doXBS
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: dRPA
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: BSE
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eW(nBas)
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin
  integer                       :: isp_W
  integer                       :: iAC
  double precision              :: lambda
  double precision,allocatable  :: Ec(:,:)

  double precision              :: EcRPA
  double precision,allocatable  :: WA(:,:)
  double precision,allocatable  :: WB(:,:)
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

  double precision,allocatable  :: Omega(:,:)
  double precision,allocatable  :: XpY(:,:,:)
  double precision,allocatable  :: XmY(:,:,:)

! Output variables

  double precision,intent(out)  :: EcAC(nspin)

! Memory allocation

  allocate(Ec(nAC,nspin))
  allocate(WA(nS,nS),WB(nS,nS),OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nBas,nBas,nS))
  allocate(Omega(nS,nspin),XpY(nS,nS,nspin),XmY(nS,nS,nspin))

! Antisymmetrized kernel version

  if(exchange_kernel) then

    write(*,*)
    write(*,*) '*** Exchange kernel version ***'
    write(*,*)

  end if

  EcAC(:) = 0d0
  Ec(:,:) = 0d0

! Compute (singlet) RPA screening 

  isp_W = 1
  EcRPA = 0d0

  call linear_response(isp_W,.true.,TDA_W,eta,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI, &
                       EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

  call BSE_static_kernel_KA(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,WA)
  call BSE_static_kernel_KB(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,WB)

! Singlet manifold

  if(singlet) then

    ispin = 1

    write(*,*) '--------------'
    write(*,*) 'Singlet states'
    write(*,*) '--------------'
    write(*,*) 

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','Ec(lambda)','Tr(K x P_lambda)'
    write(*,*) '-----------------------------------------------------------------------------------'

    do iAC=1,nAC
 
      lambda = rAC(iAC)

      if(doXBS) then

        call linear_response(isp_W,.true.,TDA_W,eta,nBas,nC,nO,nV,nR,nS,lambda,eW,ERI, &
                             EcRPA,OmRPA,XpY_RPA,XmY_RPA)
        call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)
!       call print_excitation('W^lambda:   ',isp_W,nS,OmRPA)

        call BSE_static_kernel_KA(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,OmRPA,rho_RPA,WA)
        call BSE_static_kernel_KB(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,OmRPA,rho_RPA,WB)

      end if

      call linear_response_BSE(ispin,dRPA,TDA,BSE,eta,nBas,nC,nO,nV,nR,nS,lambda,e,ERI,WA,WB, &
                               EcAC(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

      call ACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS, &
                                    ERI,XpY(:,:,ispin),XmY(:,:,ispin),Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    if(exchange_kernel) EcAC(ispin) = 0.5d0*EcAC(ispin)

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

  end if
 
! Triplet manifold

  if(triplet) then

    ispin = 2

    write(*,*) '--------------'
    write(*,*) 'Triplet states'
    write(*,*) '--------------'
    write(*,*) 

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','Ec(lambda)','Tr(K x P_lambda)'
    write(*,*) '-----------------------------------------------------------------------------------'

    do iAC=1,nAC
 
      lambda = rAC(iAC)

      if(doXBS) then

        call linear_response(isp_W,.true.,TDA_W,eta,nBas,nC,nO,nV,nR,nS,lambda,eW,ERI, &
                             EcRPA,OmRPA,XpY_RPA,XmY_RPA)
        call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

        call BSE_static_kernel_KA(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,OmRPA,rho_RPA,WA)
        call BSE_static_kernel_KB(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,OmRPA,rho_RPA,WB)

      end if  

      call linear_response_BSE(ispin,dRPA,TDA,BSE,eta,nBas,nC,nO,nV,nR,nS,lambda,e,ERI,WA,WB, &
                               EcAC(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

      call ACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,ERI,XpY(:,:,ispin),XmY(:,:,ispin),Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    if(exchange_kernel) EcAC(ispin) = 1.5d0*EcAC(ispin)

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine ACFDT
