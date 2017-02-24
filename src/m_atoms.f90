!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the structure of the molecule (atomic positions etc)
!
!=========================================================================
module m_atoms
 use m_definitions
 use m_warning,only: die,issue_warning
 use m_elements

 real(dp),parameter,private     :: tol_geom=1.0e-5_dp

 integer,public                 :: natom_basis
 integer,public                 :: natom
 integer,public                 :: nghost
 integer,protected              :: natom_type
 integer,protected              :: nbond

 real(dp),allocatable,public    :: zatom(:)
 integer,allocatable,protected  :: basis_element(:)

 real(dp),allocatable,protected :: x(:,:)
 real(dp),allocatable,protected :: vel(:,:)
 real(dp),allocatable,public    :: force(:,:)

 ! See we keep these arrays in the long-term
 real(dp),allocatable,public    :: force_nuc_nuc(:,:)
 real(dp),allocatable,public    :: force_kin(:,:)
 real(dp),allocatable,public    :: force_nuc(:,:)
 real(dp),allocatable,public    :: force_har(:,:)
 real(dp),allocatable,public    :: force_exx(:,:)
 real(dp),allocatable,public    :: force_exc(:,:)
 real(dp),allocatable,public    :: force_ovp(:,:)
 real(dp),allocatable,public    :: force_hl(:,:)

 logical,protected              :: inversion=.TRUE.
 logical,protected              :: linear=.TRUE.
 logical,protected              :: planar=.TRUE.
 real(dp),protected             :: xcenter(3)
 real(dp),protected             :: xnormal(3)


contains


!=========================================================================
subroutine init_atoms(natom_read,nghost_read,zatom_read,x_read,calculate_forces)
 use m_tools,only: cross_product
 implicit none
 integer,intent(in)  :: natom_read,nghost_read
 real(dp),intent(in) :: zatom_read(natom_read+nghost_read),x_read(3,natom_read+nghost_read)
 logical,intent(in)  :: calculate_forces
!=====
 integer  :: iatom,jatom
 real(dp) :: x21(3),x31(3)
 real(dp) :: bond_length
!=====

 natom  = natom_read
 nghost = nghost_read
 natom_basis = natom + nghost

 allocate(zatom(natom_basis))
 allocate(basis_element(natom_basis))
 allocate(x(3,natom_basis))
 ! For relaxation or dynamics only 
 if( calculate_forces ) then
   allocate(force(3,natom))
   allocate(force_nuc_nuc(3,natom))
   allocate(force_kin(3,natom))
   allocate(force_nuc(3,natom))
   allocate(force_har(3,natom))
   allocate(force_exx(3,natom))
   allocate(force_exc(3,natom))
   allocate(force_ovp(3,natom))
   allocate(force_hl(3,natom))
 endif

 zatom(1:natom)              = zatom_read(1:natom)
 ! Ghost atoms do not have a positive nucleus
 zatom(natom+1:natom+nghost) = 0.0_dp
 ! But ghost atoms have basis functions centered on them.
 basis_element(:)=NINT(zatom_read(:))

 x(:,:) = x_read(:,:)

 !
 ! Check for atoms too close
 do iatom=1,natom
   do jatom=iatom+1,natom
     if( NORM2( x(:,iatom)-x(:,jatom) ) < 0.2 ) then
       write(stdout,*) 'Atoms',iatom,jatom
       write(stdout,*) 'are closer than 0.2 bohr'
       call issue_warning('Some atoms are too close')
     endif
   enddo
 enddo

 !
 ! Find the covalent bonds based a simple distance criterium
 nbond = 0
 do iatom=1,natom
   do jatom=iatom+1,natom
     bond_length =  element_covalent_radius(zatom(iatom)) + element_covalent_radius(zatom(jatom))
     if( NORM2( x(:,iatom)-x(:,jatom) ) <  1.2_dp * bond_length  ) then
       nbond = nbond + 1
     endif
   enddo
 enddo

 !
 ! Does the molecule have inversion symmetry?
 call find_inversion()

 !
 ! Is the molecule linear, planar?
 if( natom > 2 ) then
   x21(:) = x(:,2) - x(:,1)
   do iatom=3,natom
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
   atoms_core_states = atoms_core_states + element_core(zatom(iatom),basis_element(iatom))
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

 if(ALLOCATED(zatom))         deallocate(zatom)
 if(ALLOCATED(basis_element)) deallocate(basis_element)
 if(ALLOCATED(x))             deallocate(x)
 if(ALLOCATED(force))         deallocate(force)
 if(ALLOCATED(force_nuc_nuc)) deallocate(force_nuc_nuc)
 if(ALLOCATED(force_kin))     deallocate(force_kin)
 if(ALLOCATED(force_nuc))     deallocate(force_nuc)
 if(ALLOCATED(force_har))     deallocate(force_har)
 if(ALLOCATED(force_exx))     deallocate(force_exx)
 if(ALLOCATED(force_exc))     deallocate(force_exc)
 if(ALLOCATED(force_ovp))     deallocate(force_ovp)
 if(ALLOCATED(force_hl))      deallocate(force_hl)

end subroutine destroy_atoms


!=========================================================================
subroutine relax_atoms(lbfgs_plan)
 use m_lbfgs,only: lbfgs_state,lbfgs_execute
 implicit none

 type(lbfgs_state),intent(inout) :: lbfgs_plan
!=====

 call lbfgs_execute(lbfgs_plan)

end subroutine relax_atoms


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


!=========================================================================
subroutine nucleus_nucleus_force()
 implicit none
!=====
 integer              :: iatom,jatom
!=====

 force_nuc_nuc(:,:) = 0.0_dp
 do iatom=1,natom
   do jatom=1,natom
     if( iatom == jatom ) cycle
     force_nuc_nuc(:,iatom) = force_nuc_nuc(:,iatom) &
                               + zatom(iatom) * zatom(jatom) / ( SUM( (x(:,iatom) - x(:,jatom))**2) )**1.50_dp &
                                  * ( x(:,iatom) - x(:,jatom) )
   enddo
 enddo

end subroutine nucleus_nucleus_force


!=========================================================================
subroutine find_inversion()
 implicit none
!=====
 integer  :: iatom,jatom
 logical  :: found
 real(dp) :: xtmp(3)
!=====

 xcenter(:) = 0.0_dp
 do iatom=1,natom
   xcenter(:) = xcenter(:) + x(:,iatom) / REAL(natom,dp)
 enddo

 do iatom=1,natom
   xtmp(:) = 2.0_dp * xcenter(:) - x(:,iatom)
   found = .FALSE.
   do jatom=1,natom
     if( NORM2( xtmp(:) - x(:,jatom) ) < tol_geom ) then
       if( ABS(zatom(iatom)-zatom(jatom)) < tol_geom ) found = .TRUE.
       exit
     endif
   enddo
   inversion = inversion .AND. found
 enddo

end subroutine find_inversion


!=========================================================================
function same_element(iatom,jatom)
 implicit none
 integer,intent(in) :: iatom,jatom
 logical :: same_element
!=====

 same_element = ( ABS( zatom(iatom) - zatom(jatom) ) < 1.0e-5 )

end function same_element



end module m_atoms
!=========================================================================
