!=========================================================================
module m_atoms
 use m_definitions

 integer                      :: natom

 real(dp),allocatable         :: zatom(:)
 integer,allocatable          :: basis_element(:)

 real(dp),allocatable         :: x(:,:)


end module m_atoms
