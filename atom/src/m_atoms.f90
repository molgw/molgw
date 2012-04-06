!=========================================================================
module m_atoms
 use m_definitions

 integer                      :: natom

 real(dp),allocatable         :: zatom(:)
 integer,allocatable          :: basis_element(:)

 real(dp),allocatable         :: x(:,:)

contains

!=========================================================================
subroutine nucleus_nucleus_energy(energy)
 implicit none
 real(dp),intent(out) :: energy
!=====
 integer              :: iatom,jatom
!=====

 energy = 0.0_dp
 do iatom=1,natom
   do jatom=iatom+1,natom
     energy = energy + zatom(iatom) * zatom(jatom) / SQRT( SUM( (x(:,iatom) - x(:,jatom))**2) )
   enddo
 enddo

end subroutine nucleus_nucleus_energy

end module m_atoms
