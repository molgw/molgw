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
 use m_linear_algebra,only: cross_product

 real(dp),parameter,private     :: tol_geom=1.0e-5_dp

 integer,public                 :: natom_basis
 integer,public                 :: natom
 integer,public                 :: nghost
 integer,protected              :: natom_type
 integer,protected              :: nbond
 integer,public                 :: nprojectile

 real(dp),allocatable,public    :: zvalence(:)
 real(dp),allocatable,protected :: zatom(:)
 integer,allocatable,protected  :: zbasis(:)

 real(dp),allocatable,protected :: xatom(:,:)
 real(dp),allocatable,protected :: xbasis(:,:)
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
subroutine init_atoms(zatom_read,x_read,vel_projectile,calculate_forces,excit_name)
 implicit none
 real(dp),intent(in) :: zatom_read(natom+nghost),x_read(3,natom+nghost)
 real(dp),intent(in) :: vel_projectile(3)
 logical,intent(in)  :: calculate_forces
 character(len=12),intent(in) :: excit_name
!=====
 integer  :: iatom,jatom
 real(dp) :: x21(3),x31(3)
 real(dp) :: bond_length
!=====

 ! here x_read contains all coordinates: first (natom-nprojectile) real atoms, then (nghost) ghost atoms and the last but not least (nprojectile) projectile (which is 0 or 1)

 ! xatom and zatom designate the physical nuclei
 allocate(xatom(3,natom))
 allocate(zatom(natom))
 allocate(zvalence(natom))
 ! xbasis and zbasis designate the basis centers and nature
 allocate(zbasis(natom_basis))
 allocate(xbasis(3,natom_basis))

 allocate(vel(3,natom))
 vel(:,:) = 0.0_dp

 if( excit_name == "NUCLEUS" .OR. excit_name == "ANTINUCLEUS" .OR. excit_name == 'ION' ) then
   vel(:,natom) = vel_projectile(:)
 endif
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

 ! List of atoms is organized as follows:
 ! 1. physical atoms   :    nucleus | basis
 ! 2. ghost atoms      :      no    | basis
 ! 3. projectile       :    nucleus |   yes if 'ION' / no if 'NUCLEUS' or 'ANTINUCLEUS'
 !
 ! natom       contains the number of sites having a nucleus: number of physical atoms + number of ionic projectiles (0 or 1)
 ! natom_basis contains the number of sites having basis functions:  number of physical atoms + number of ghost atoms
 !
 if( nprojectile == 0 ) then
   xatom(:,1:natom) = x_read(:,1:natom)
   zatom(1:natom)   = zatom_read(1:natom)
 else
   if( excit_name == 'ION' .OR. excit_name == "ANTIION" ) then
     xatom(:,1:natom-nprojectile) = x_read(:,1:natom-nprojectile)
     zatom(1:natom-nprojectile)   = zatom_read(1:natom-nprojectile)
     xatom(:,natom)                = x_read(:,natom_basis)
     zatom(natom)                  = zatom_read(natom_basis)
   else
     xatom(:,1:natom_basis-nghost) = x_read(:,1:natom_basis-nghost)
     zatom(1:natom_basis-nghost)   = zatom_read(1:natom_basis-nghost)
     xatom(:,natom)                = x_read(:,natom_basis+nprojectile)
     zatom(natom)                  = zatom_read(natom_basis+nprojectile)
   endif
 endif


 if( excit_name == "ANTINUCLEUS" .OR. excit_name == "ANTIION" ) then
   zatom(natom) = -zatom(natom)
 endif

 ! Ghost atoms do not have a positive nucleus
 !zatom(natom+1:natom+nghost) = 0.0_dp
 ! But ghost atoms have basis functions centered on them.
 zbasis(1:natom_basis)   = NINT(zatom_read(1:natom_basis))
 xbasis(:,1:natom_basis) = x_read(:,1:natom_basis)

 !
 ! Check for atoms too close
 do iatom=1,natom
   do jatom=iatom+1,natom
     if( NORM2( xatom(:,iatom)-xatom(:,jatom) ) < 0.2 ) then
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
     bond_length =  element_covalent_radius(NINT(zatom(iatom))) + element_covalent_radius(NINT(zatom(jatom)))
     if( NORM2( xatom(:,iatom)-xatom(:,jatom) ) <  1.2_dp * bond_length  ) then
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
   x21(:) = xatom(:,2) - xatom(:,1)
   do iatom=3,natom-nprojectile
     x31(:) = xatom(:,iatom) - xatom(:,1)
     call cross_product(x21,x31,xnormal)
     if( NORM2(xnormal(:)) > tol_geom ) then
       xnormal(:) = xnormal(:) / NORM2(xnormal(:))
       linear=.FALSE.
       exit
     endif
   enddo
   if( .NOT. linear) then
     do iatom=1,natom-nprojectile
       if( ABS(DOT_PRODUCT( xatom(:,iatom) , xnormal(:) )) > tol_geom ) planar=.FALSE.
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
   atoms_core_states = atoms_core_states + element_core(zvalence(iatom),zatom(iatom))
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
     if( NORM2( xatom(:,iatom)-xatom(:,jatom) ) < 4.0 ) then
       jbond = jbond + 1
       if(jbond==ibond) then
         xbond(:) = 0.5_dp * ( xatom(:,iatom) - xatom(:,jatom) )
         return
       endif
     endif
   enddo
 enddo


end subroutine get_bondcenter


!=========================================================================
subroutine change_position_one_atom(iatom,xposition)
 implicit none
 integer,intent(in)   :: iatom
 real(dp),intent(in)  :: xposition(3)
!=====

 xatom(:,iatom) = xposition(:)

end subroutine change_position_one_atom


!=========================================================================
subroutine destroy_atoms()
 implicit none
!=====

 if(ALLOCATED(zatom))         deallocate(zatom)
 if(ALLOCATED(zvalence))      deallocate(zvalence)
 if(ALLOCATED(zbasis))        deallocate(zbasis)
 if(ALLOCATED(xatom))         deallocate(xatom)
 if(ALLOCATED(xbasis))        deallocate(xbasis)
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
subroutine relax_atoms(lbfgs_plan,etotal)
 use m_lbfgs
 implicit none

 type(lbfgs_state),intent(inout) :: lbfgs_plan
 real(dp),intent(in)             :: etotal
!=====
 integer  :: info,iatom,idir
 real(dp) :: xnew(3,natom)
!=====

 xnew(:,:) = xatom(:,:)

 info = lbfgs_execute(lbfgs_plan,xnew,etotal,-force)

 ! Do not move the atoms by more than 0.20 bohr
 do iatom=1,natom-nprojectile
   do idir=1,3
     if( ABS( xnew(idir,iatom) - xatom(idir,iatom) ) > 0.20_dp ) then
       xnew(idir,iatom) = xatom(idir,iatom) + SIGN( 0.20_dp , xnew(idir,iatom) - xatom(idir,iatom) )
     endif
   enddo
 enddo

 xatom(:,:)       = xnew(:,:)
 ! Ghost atoms never move
 xbasis(:,:natom) = xnew(:,:)

end subroutine relax_atoms


!=========================================================================
subroutine output_positions()
 implicit none
!=====
 integer :: ighost,iatom
!=====

 write(stdout,'(/,1x,a)') '================================'
 write(stdout,*) '      Atom list'
 write(stdout,*) '                       bohr                                        angstrom'
 do iatom=1,natom-nprojectile
   write(stdout,'(1x,a,i3,2x,a2,a,3(1x,f12.6),6x,3(1x,f12.6))') 'atom  ',iatom, &
                                                           element_name(REAL(zatom(iatom),dp)),': ',  &
                                                           xatom(:,iatom),xatom(:,iatom)*bohr_A
 enddo

 if( nghost > 0 ) write(stdout,'(a)') ' == ghost list'
 do ighost=1,nghost
   write(stdout,'(1x,a,i3,2x,a2,a,3(1x,f12.6),6x,3(1x,f12.6))') 'ghost ',iatom, &
                                           element_name(REAL(zbasis(natom-nprojectile+ighost),dp)),': ',  &
                                           xbasis(:,natom-nprojectile+ighost),xbasis(:,natom-nprojectile+ighost)*bohr_A
 enddo
 if( nprojectile > 0 ) then
   write(stdout,'(a)') ' == projectile'
   write(stdout,'(1x,a,i3,2x,a2,a,3(1x,f12.6),6x,3(1x,f12.6))') 'atom  ',natom+nghost, &
                                                           element_name(REAL(zatom(natom),dp)),': ',  &
                                                           xatom(:,natom),xatom(:,natom)*bohr_A
 endif


 write(stdout,'(1x,a,/)') '================================'

end subroutine output_positions


!=========================================================================
subroutine output_projectile_position()
 implicit none
!=====

 write(stdout,*)
 if( nprojectile > 0 ) then
   write(stdout,'(a)') ' === projectile position: ----------bohr---------------    |||   ------------- angstrom----------==='

   if( zatom(natom) > 0 ) then
     write(stdout,'(1x,a,i3,2x,a2,a,3(1x,f12.6),6x,3(1x,f12.6))') 'atom  ',natom+nghost, &
                                                           element_name(REAL(zatom(natom),dp)),': ',  &
                                                           xatom(:,natom),xatom(:,natom)*bohr_A
   else
     write(stdout,'(1x,a,i3,2x,a4,a2,a,3(1x,f12.6),6x,3(1x,f12.6))') 'atom  ',natom+nghost, &
                                                           'anti',element_name(REAL(zatom(natom),dp)),': ',  &
                                                           xatom(:,natom),xatom(:,natom)*bohr_A
   endif

 endif

end subroutine output_projectile_position


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
     energy = energy + zvalence(iatom) * zvalence(jatom) / SQRT( SUM( (xatom(:,iatom) - xatom(:,jatom))**2) )
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
                               + zvalence(iatom) * zvalence(jatom) / ( SUM( (xatom(:,iatom) - xatom(:,jatom))**2) )**1.50_dp &
                                  * ( xatom(:,iatom) - xatom(:,jatom) )
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
 do iatom=1,natom-nprojectile
   xcenter(:) = xcenter(:) + xatom(:,iatom) / REAL(natom,dp)
 enddo

 do iatom=1,natom-nprojectile
   xtmp(:) = 2.0_dp * xcenter(:) - xatom(:,iatom)
   found = .FALSE.
   do jatom=1,natom-nprojectile
     if( NORM2( xtmp(:) - xatom(:,jatom) ) < tol_geom ) then
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
