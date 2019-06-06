!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the subroutines to calculate the static dipole and static quadrupole tensor
!
!=========================================================================
subroutine static_dipole(nstate,basis,occupation,c_matrix)
 use m_definitions
 use m_basis_set
 use m_inputparam,only: nspin
 use m_hamiltonian_tools,only: setup_density_matrix
 use m_hamiltonian_onebody,only: calculate_dipole_basis
 use m_atoms
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),c_matrix(basis%nbf,nstate,nspin)
!=====
 integer                    :: iatom,idir
 real(dp)                   :: dipole(3)
 real(dp),allocatable       :: dipole_basis(:,:,:)
 real(dp)                   :: p_matrix(basis%nbf,basis%nbf,nspin)
!=====


! call start_clock(timing_spectrum)

 write(stdout,'(/,a)') ' Calculate the static dipole'


 !
 ! First precalculate all the needed dipole in the basis set
 !
 call calculate_dipole_basis(basis,dipole_basis)

 call setup_density_matrix(c_matrix,occupation,p_matrix)

 ! Minus sign for electrons
 do idir=1,3
   dipole(idir) = -SUM( dipole_basis(:,:,idir) * SUM( p_matrix(:,:,:) , DIM=3 ) )
 enddo

 deallocate(dipole_basis)

 do iatom=1,natom
   dipole(:) = dipole(:) + zvalence(iatom) * xatom(:,iatom)
 enddo

 write(stdout,'(1x,a,3(2x,f14.6))') 'Dipole (a.u.):  ',dipole(:)
 write(stdout,'(1x,a,3(2x,f14.6))') 'Dipole (Debye): ',dipole(:) * au_debye


end subroutine static_dipole


!=========================================================================
subroutine static_quadrupole(nstate,basis,occupation,c_matrix)
 use m_definitions
 use m_basis_set
 use m_inputparam,only: nspin
 use m_hamiltonian_tools,only: setup_density_matrix
 use m_hamiltonian_onebody,only: calculate_quadrupole_basis
 implicit none

 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin),c_matrix(basis%nbf,nstate,nspin)
!=====
 integer                    :: iatom,idir,jdir
 real(dp)                   :: trace
 real(dp)                   :: quad(3,3)
 real(dp),allocatable       :: quad_basis(:,:,:,:)
 real(dp)                   :: p_matrix(basis%nbf,basis%nbf,nspin)
!=====


! call start_clock(timing_spectrum)

 write(stdout,'(/,a)') ' Calculate the static quadrupole'


 !
 ! First precalculate all the needed quadrupoles in the basis set
 !
 call calculate_quadrupole_basis(basis,quad_basis)

 call setup_density_matrix(c_matrix,occupation,p_matrix)

 ! Minus sign for electrons
 do jdir=1,3
   do idir=1,3
     quad(idir,jdir) = -SUM( quad_basis(:,:,idir,jdir) * SUM( p_matrix(:,:,:) , DIM=3 ) )
   enddo
 enddo

 deallocate(quad_basis)

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
