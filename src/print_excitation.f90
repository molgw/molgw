subroutine print_excitation(method,ispin,nS,Omega)

! Print excitation energies for a given spin manifold

  implicit none
  include 'parameters.h'

! Input variables

  character*12,intent(in)            :: method
  integer,intent(in)                 :: ispin,nS
  double precision,intent(in)        :: Omega(nS)

! Local variables

  character*14                       :: spin_manifold
  integer,parameter                  :: maxS = 20
  integer                            :: ia

  if(ispin == 1) spin_manifold = 'singlet'
  if(ispin == 2) spin_manifold = 'triplet'
  if(ispin == 3) spin_manifold = 'alpha-beta'
  if(ispin == 4) spin_manifold = 'alpha-alpha' 
  if(ispin == 5) spin_manifold = 'spin-conserved'
  if(ispin == 6) spin_manifold = 'spin-flip'
  if(ispin == 7) spin_manifold = 'beta-beta'

  write(*,*)
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A14,A14,A14,A9)') method,' calculation: ',spin_manifold,' manifold'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'

  do ia=1,min(maxS,nS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') & 
      '|',ia,'|',Omega(ia),'|',Omega(ia)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------'
  write(*,*)

end subroutine print_excitation


