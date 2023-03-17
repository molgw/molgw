program main
logical::TDA=.false.
logical::doACFDT=.true.
logical::exchange_kernel=.true.
logical::singlet=.true.
logical::triplet=.true.
double precision::eta=0.01
integer::nBas=4
integer::nC=0
integer::nO=1
integer::nV=3
integer::nR ! Dunno this one
integer::nS ! Dunno this one
double precision::ENuc
double precision::ERHF
double precision,allocatable::eHF(:)
double precision,allocatable::ERI(:,:,:,:)
double precision,allocatable::dipole_int(:,:,:)

allocate(eHF(nBas),ERI(nBas,nBas,nBas,nBas),dipole_int(nBas,nBas,ncart))
eHF=0.0e0;ERI=0.0e0;dipole_int=0.0e0;

call RPAx(TDA,doACFDT,exchange_kernel,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

deallocate(eHF,ERI,dipole_int)

end program main
