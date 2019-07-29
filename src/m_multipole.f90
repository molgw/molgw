!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the subroutines to calculate the static dipole and static quadrupole tensor
!
!=========================================================================
module m_multipole
 use m_definitions
 use m_basis_set
 use m_inputparam,only: nspin
 use m_hamiltonian_tools,only: setup_density_matrix
 use m_hamiltonian_onebody,only: calculate_dipole_ao,calculate_quadrupole_ao
 use m_atoms




contains


!=========================================================================
! Many optional arguments however not all combinations are possible (be careful!)
subroutine static_dipole(basis,occupation,c_matrix_in,p_matrix_in,dipole_ao_in,dipole_out)
  implicit none

  type(basis_set),intent(in)       :: basis
  real(dp),optional,intent(in)     :: occupation(:,:)
  real(dp),optional,intent(in)     :: c_matrix_in(:,:,:)
  complex(dp),optional,intent(in)  :: p_matrix_in(:,:,:)
  real(dp),optional,intent(in)     :: dipole_ao_in(:,:,:)
  real(dp),optional,intent(out)    :: dipole_out(3)
  !=====
  integer                    :: iatom,idir
  real(dp)                   :: dipole(3)
  real(dp),allocatable       :: dipole_ao(:,:,:)
  real(dp),allocatable       :: p_matrix(:,:,:)
  !=====


  write(stdout,'(/,a)') ' Calculate the static dipole'

  if( .NOT. PRESENT(c_matrix_in) .AND. .NOT. PRESENT(p_matrix_in) ) call die('static_dipole: should provide either C or P matrix')
  if( .NOT. PRESENT(occupation)  .AND. .NOT. PRESENT(p_matrix_in) ) call die('static_dipole: need occupation when P is not here')

  allocate(dipole_ao(basis%nbf,basis%nbf,3))
  if( .NOT. PRESENT(dipole_ao_in) ) then
    !
    ! First precalculate all the needed dipole in the basis set
    !
    call calculate_dipole_ao(basis,dipole_ao)
  else
    dipole_ao(:,:,:) = dipole_ao_in(:,:,:)
  endif

  allocate(p_matrix(basis%nbf,basis%nbf,nspin))
  if( .NOT. PRESENT(p_matrix_in) ) then
    call setup_density_matrix(c_matrix_in,occupation,p_matrix)
  else
    p_matrix(:,:,:) = p_matrix_in(:,:,:)%re
  endif


  ! Calculate electronic part with a minus sign (for electrons)
  do idir=1,3
    dipole(idir) = -SUM( dipole_ao(:,:,idir) * SUM( p_matrix(:,:,:) , DIM=3 ) )
  enddo

  deallocate(dipole_ao,p_matrix)

  ! Add the nuclear part
  do iatom=1,natom
    dipole(:) = dipole(:) + zvalence(iatom) * xatom(:,iatom)
  enddo

  if( .NOT. PRESENT(dipole_out) ) then
    write(stdout,'(1x,a,3(2x,f14.6))') 'Dipole (a.u.):  ',dipole(:)
    write(stdout,'(1x,a,3(2x,f14.6))') 'Dipole (Debye): ',dipole(:) * au_debye
  else
    dipole_out(:) = dipole(:)
  endif


end subroutine static_dipole


!=========================================================================
subroutine static_quadrupole(basis,occupation,c_matrix)
  implicit none

  type(basis_set),intent(in) :: basis
  real(dp),intent(in)        :: occupation(:,:),c_matrix(:,:,:)
  !=====
  integer                    :: iatom,idir,jdir
  real(dp)                   :: trace
  real(dp)                   :: quad(3,3)
  real(dp),allocatable       :: quad_ao(:,:,:,:)
  real(dp)                   :: p_matrix(basis%nbf,basis%nbf,nspin)
  !=====


  write(stdout,'(/,a)') ' Calculate the static quadrupole'


  !
  ! First precalculate all the needed quadrupoles in the basis set
  !
  call calculate_quadrupole_ao(basis,quad_ao)

  call setup_density_matrix(c_matrix,occupation,p_matrix)

  ! Minus sign for electrons
  do jdir=1,3
    do idir=1,3
      quad(idir,jdir) = -SUM( quad_ao(:,:,idir,jdir) * SUM( p_matrix(:,:,:) , DIM=3 ) )
    enddo
  enddo

  deallocate(quad_ao)

  do iatom=1,natom
    do jdir=1,3
      quad(:,jdir) = quad(:,jdir) + zvalence(iatom) * xatom(:,iatom) * xatom(jdir,iatom)
    enddo
  enddo

  write(stdout,'(1x,a,3(2x,f14.6))') 'Quadrupole (a.u.):    ',quad(1,:)
  write(stdout,'(1x,a,3(2x,f14.6))') '                      ',quad(2,:)
  write(stdout,'(1x,a,3(2x,f14.6))') '                      ',quad(3,:)
  write(stdout,*)
  write(stdout,'(1x,a,3(2x,f14.6))') 'Quadrupole (Debye.A): ',quad(1,:) * au_debye * bohr_A
  write(stdout,'(1x,a,3(2x,f14.6))') '                      ',quad(2,:) * au_debye * bohr_A
  write(stdout,'(1x,a,3(2x,f14.6))') '                      ',quad(3,:) * au_debye * bohr_A

  trace = quad(1,1) + quad(2,2) + quad(3,3)
  quad(1,1) = quad(1,1) - trace / 3.0_dp
  quad(2,2) = quad(2,2) - trace / 3.0_dp
  quad(3,3) = quad(3,3) - trace / 3.0_dp

  write(stdout,*)
  write(stdout,'(1x,a,3(2x,f14.6))') 'Traceless quadrupole (a.u.):    ',quad(1,:)
  write(stdout,'(1x,a,3(2x,f14.6))') '                                ',quad(2,:)
  write(stdout,'(1x,a,3(2x,f14.6))') '                                ',quad(3,:)
  write(stdout,*)
  write(stdout,'(1x,a,3(2x,f14.6))') 'Traceless quadrupole (Debye.A): ',quad(1,:) * au_debye * bohr_A
  write(stdout,'(1x,a,3(2x,f14.6))') '                                ',quad(2,:) * au_debye * bohr_A
  write(stdout,'(1x,a,3(2x,f14.6))') '                                ',quad(3,:) * au_debye * bohr_A



end subroutine static_quadrupole


!=========================================================================
subroutine spatial_extension(basis,c_matrix)
  implicit none

  type(basis_set),intent(in) :: basis
  real(dp),intent(in)        :: c_matrix(:,:,:)
  !=====
  integer                    :: nbf,istate,ispin
  real(dp),allocatable       :: dipole_ao(:,:,:)
  real(dp),allocatable       :: quad_ao(:,:,:,:)
  real(dp),allocatable       :: trace_dipole_ao(:,:)
  real(dp),allocatable       :: trace_quad_ao(:,:)
  real(dp)                   :: variance(nspin),mean(nspin)
  !=====

  nbf = SIZE(c_matrix,DIM=1)

  write(stdout,'(/,1x,a)') 'Calculate the wavefunction spatial extension'

  !
  ! First precalculate all the dipole and quadrupole in the AO basis set
  !
  call calculate_dipole_ao(basis,dipole_ao)
  call calculate_quadrupole_ao(basis,quad_ao)
  allocate(trace_dipole_ao(nbf,nbf))
  allocate(trace_quad_ao(nbf,nbf))
  trace_dipole_ao(:,:) =  dipole_ao(:,:,1) + dipole_ao(:,:,2) + dipole_ao(:,:,3)
  trace_quad_ao(:,:)   =  quad_ao(:,:,1,1) + quad_ao(:,:,2,2) + quad_ao(:,:,3,3)
  deallocate(dipole_ao)
  deallocate(quad_ao)

  write(stdout,'(/,1x,70("="))')
  write(stdout,'(1x,a,/)') ' State index          spatial extension (bohr)'
  do istate=1,SIZE(c_matrix,DIM=2)
    do ispin=1,nspin

      variance(ispin) = DOT_PRODUCT( c_matrix(:,istate,ispin) , MATMUL( trace_quad_ao(:,:)   , c_matrix(:,istate,ispin) ) )
      mean(ispin)     = DOT_PRODUCT( c_matrix(:,istate,ispin) , MATMUL( trace_dipole_ao(:,:) , c_matrix(:,istate,ispin) ) )

    enddo

    write(stdout,'(1x,i6,2x,*(2x,f12.3))') istate,SQRT( variance(:) - mean(:)**2 )

  enddo
  write(stdout,'(1x,70("="),/)')


  deallocate(trace_dipole_ao)
  deallocate(trace_quad_ao)


end subroutine spatial_extension


!=========================================================================
end module m_multipole
!=========================================================================
