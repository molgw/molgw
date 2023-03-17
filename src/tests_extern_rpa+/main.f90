program main
logical::TDA=.false.
logical::doACFDT=.true.
! The logical variable exchange_kernel is the key quantity!
! RPA+ (RPAx-II)
! logical::exchange_kernel=.true.
! RPA+ in T2's paper (RPAx-I)
logical::exchange_kernel=.false.
logical::singlet=.true.
logical::triplet=.true.
double precision::eta=0.d0
! nC   = number of core orbitals
! nO   = number of occupied orbitals
! nV   = number of virtual orbitals (see below)
! nR   = number of Rydberg orbitals 
! nBas = number of basis functions (see below)
!      = nO + nV
! nS   = number of single excitation 
!      = nO*nV
integer::i,j,k,l,m,neris
integer::nBas,nC,nO,nV,nR,nS,ncart=3
double precision::ENuc
double precision::ERHF
double precision,allocatable::eHF(:)
double precision,allocatable::ERI(:,:,:,:)
double precision,allocatable::dipole_int(:,:,:)
! Read from MOLGW DUMP file. 
! ./main < DUMP 
read(*,*) nBas
read(*,*) nC
read(*,*) nO
read(*,*) nV
read(*,*) nR
read(*,*) ERHF
ns=nO*nV
allocate(eHF(nBas),ERI(nBas,nBas,nBas,nBas),dipole_int(nBas,nBas,ncart))
eHF=0.0e0;ERI=0.0e0;dipole_int=0.0e0;
do i=1,nBAS
 read(*,*) eHF(i)
enddo
do i=1,nBAS
 do j=1,nBAS
  do k=1,nBAS
   do l=1,nBAS
     read(*,*) ERI(i,j,k,l)
   enddo
  enddo
 enddo
enddo

allocate(eHF(nBas),ERI(nBas,nBas,nBas,nBas),dipole_int(nBas,nBas,ncart))
eHF=0.0e0;ERI=0.0e0;dipole_int=0.0e0;

call RPAx(TDA,doACFDT,exchange_kernel,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

deallocate(eHF,ERI,dipole_int)

end program main
