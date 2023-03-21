subroutine print_transition_vectors(spin_allowed,nBas,nC,nO,nV,nR,nS,dipole_int,Omega,XpY,XmY)

! Print transition vectors for linear response calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: spin_allowed
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision              :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

! Local variables

  integer                       :: ia,jb,j,b
  integer                       :: maxS = 10
  double precision              :: S2
  double precision,parameter    :: thres_vec = 0.1d0
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)
  double precision,allocatable  :: os(:)

! Memory allocation

  maxS = min(nS,maxS)
  allocate(X(nS),Y(nS),os(maxS))

! Compute oscillator strengths

  os(:) = 0d0
  if(spin_allowed) call oscillator_strength(nBas,nC,nO,nV,nR,nS,maxS,dipole_int,Omega,XpY,XmY,os)

! Print details about excitations

  do ia=1,maxS

    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    ! <S**2> values

    if(spin_allowed) then 
      S2 = 0d0
    else
      S2 = 2d0
    end if

    print*,'-------------------------------------------------------------'
    write(*,'(A15,I3,A2,F10.6,A3,A6,F6.4,A11,F6.4)') &
            ' Excitation n. ',ia,': ',Omega(ia)*HaToeV,' eV','  f = ',os(ia),'  <S**2> = ',S2
    print*,'-------------------------------------------------------------'

    jb = 0
    do j=nC+1,nO
      do b=nO+1,nBas-nR
        jb = jb + 1
        if(abs(X(jb)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') j,' -> ',b,' = ',X(jb)/sqrt(2d0)
      end do
    end do
 
    jb = 0
    do j=nC+1,nO
      do b=nO+1,nBas-nR
        jb = jb + 1
        if(abs(Y(jb)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') j,' <- ',b,' = ',Y(jb)/sqrt(2d0)
      end do
    end do
   write(*,*)

  end do

! Thomas-Reiche-Kuhn sum rule

  write(*,'(A30,F10.6)') 'Thomas-Reiche-Kuhn sum rule = ',sum(os(:))
  write(*,*)

end subroutine print_transition_vectors
