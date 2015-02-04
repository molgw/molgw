!=========================================================================
#include "macros.h"
!=========================================================================
module m_atoms
 use m_definitions
 use m_mpi
 use m_elements

 real(dp),parameter,private     :: tol_geom=1.0e-5_dp

 integer,protected              :: natom
 integer,protected              :: natom_type
 integer,protected              :: nbond

 real(dp),allocatable,protected :: zatom(:)
 integer,allocatable,protected  :: basis_element(:)

 real(dp),allocatable,protected :: x(:,:)

 logical,protected              :: inversion=.TRUE.
 logical,protected              :: linear=.TRUE.
 logical,protected              :: planar=.TRUE.
 real(dp),protected             :: xinversion(3)
 real(dp),protected             :: xnormal(3)


contains


!=========================================================================
subroutine init_atoms(natom_read,zatom_read,x_read)
 use m_tools,only: cross_product
 implicit none
 integer,intent(in)  :: natom_read
 real(dp),intent(in) :: zatom_read(natom_read),x_read(3,natom_read)
!=====
 integer  :: iatom,jatom
 real(dp) :: xtmp(3),x21(3),x31(3),xnormal(3)
 logical  :: found
!=====

 natom = natom_read
 allocate(zatom(natom))
 allocate(basis_element(natom))
 allocate(x(3,natom))

 zatom(:) = zatom_read(:)
 basis_element(:)=NINT(zatom(:))

 x(:,:) = x_read(:,:)

 !
 ! Check for atoms too close
 do iatom=1,natom
   do jatom=iatom+1,natom
     if( NORM2( x(:,iatom)-x(:,jatom) ) < 0.2 ) then
       WRITE_MASTER(*,*) 'Atoms',iatom,jatom
       WRITE_MASTER(*,*) 'are closer than 0.2 bohr'
       stop'stop here'
     endif
   enddo
 enddo

 !
 ! Find the covalent bond is a simple distance criterium
 nbond = 0
 do iatom=1,natom
   do jatom=iatom+1,natom
     if( NORM2( x(:,iatom)-x(:,jatom) ) < 4.0 ) then
       nbond = nbond + 1
     endif
   enddo
 enddo

 !
 ! Does the molecule have inversion symmetry?
 xinversion(:) = 0.0_dp
 do iatom=1,natom
   xinversion(:) = xinversion(:) + x(:,iatom) / REAL(natom,dp)
 enddo

 do iatom=1,natom
   xtmp(:) = 2.0_dp * xinversion(:) - x(:,iatom)
   found = .FALSE.
   do jatom=1,natom
     if( NORM2( xtmp(:) - x(:,jatom) ) < tol_geom ) then
       if( ABS(zatom(iatom)-zatom(jatom)) < tol_geom ) found = .TRUE.
       exit
     endif
   enddo
   inversion = inversion .AND. found
 enddo

 !
 ! Is the molecule linear, planar?
 if( natom > 2 ) then
   x21(:) = x(:,2) - x(:,1)
   do iatom=1,natom
     x31(:) = x(:,iatom) - x(:,1)
     call cross_product(x21,x31,xnormal)
     if( NORM2(xnormal(:)) > tol_geom ) then
       xnormal(:) = xnormal(:) / NORM2(xnormal(:))
       linear=.FALSE.
       exit
     endif
   enddo
   if( .NOT. linear) then
     do iatom=1,natom
       if( ABS(DOT_PRODUCT( x(:,iatom) , xnormal(:) )) > tol_geom ) planar=.FALSE.
     enddo
   else
     planar=.FALSE.
   endif
 else
  ! Molecule is linear
  ! Set planar to FALSE for safety
  linear=.TRUE.
  planar=.FALSE.
 endif





end subroutine init_atoms


!=========================================================================
function atoms_core_states()
 implicit none
 integer :: atoms_core_states
!=====
 integer :: iatom
!=====

 atoms_core_states=0
 do iatom=1,natom
   atoms_core_states = atoms_core_states + element_core(zatom(iatom))
 enddo
 
end function atoms_core_states


!=========================================================================
subroutine get_bondcenter(ibond,xbond)
 implicit none
 integer,intent(in)   :: ibond
 real(dp),intent(out) :: xbond(3)
!=====
 integer :: iatom,jatom,jbond
!=====

 jbond = 0
 do iatom=1,natom
   do jatom=iatom+1,natom
     if( NORM2( x(:,iatom)-x(:,jatom) ) < 4.0 ) then
       jbond = jbond + 1
       if(jbond==ibond) then
         xbond(:) = 0.5_dp * ( x(:,iatom) - x(:,jatom) )
         return
       endif
     endif
   enddo
 enddo


end subroutine get_bondcenter


!=========================================================================
subroutine destroy_atoms()
 implicit none
!=====

 if(allocated(zatom)) deallocate(zatom)
 if(allocated(basis_element)) deallocate(basis_element)
 if(allocated(x)) deallocate(x)

end subroutine destroy_atoms


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
!=========================================================================
