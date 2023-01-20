!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the structure of the molecule (atomic positions etc)
!
!=========================================================================
#include "molgw.h"
module m_atoms
  use m_definitions
  use m_warning,only: die,issue_warning
  use m_elements
  use m_linear_algebra,only: cross_product

  real(dp),parameter,private     :: tol_geom=1.0e-5_dp

  integer,public                 :: ncenter_basis
  integer,public                 :: ncenter_nuclei
  integer,protected              :: nbond
  integer,public                 :: nprojectile
  integer,private                :: nghost_

  real(dp),allocatable,public    :: zvalence(:)
  real(dp),allocatable,protected :: zatom(:)
  integer,allocatable,protected  :: zbasis(:)

  real(dp),allocatable,protected :: xatom(:,:)
  real(dp),allocatable,protected :: xbasis(:,:)
  real(dp),allocatable,protected :: vel_nuclei(:,:)
  real(dp),allocatable,protected :: vel_basis(:,:)
  real(dp),allocatable,public    :: force(:,:)
  real(dp),public                :: force_projectile(3) = 0.0_dp

  ! See if we keep these arrays in the long-term
  real(dp),allocatable,public    :: force_nuc_nuc(:,:)
  real(dp),allocatable,public    :: force_kin(:,:)
  real(dp),allocatable,public    :: force_nuc(:,:)
  real(dp),allocatable,public    :: force_har(:,:)
  real(dp),allocatable,public    :: force_exx(:,:)
  real(dp),allocatable,public    :: force_exc(:,:)
  real(dp),allocatable,public    :: force_ovp(:,:)
  real(dp),allocatable,public    :: force_hellfeyn(:,:)

  logical,protected              :: inversion=.TRUE.
  logical,protected              :: linear=.TRUE.
  logical,protected              :: planar=.TRUE.
  real(dp),protected             :: xcenter(3)
  real(dp),protected             :: xnormal(3)


contains


!=========================================================================
subroutine init_atoms(natom_in,nghost_in,nucleus_wo_basis,zatom_read,x_read,vel_projectile, &
                      calculate_forces,excit_name,projectile_charge_scaling)
  implicit none

  integer,intent(in)  :: natom_in,nghost_in
  logical,intent(in)  :: nucleus_wo_basis(:)
  real(dp),intent(in) :: zatom_read(:),x_read(:,:)
  real(dp),intent(in) :: vel_projectile(3)
  logical,intent(in)  :: calculate_forces
  character(len=12),intent(in)   :: excit_name
  real(dp),intent(in) :: projectile_charge_scaling
  !=====
  integer  :: natom_read
  integer  :: iatom,jatom,jcenter
  real(dp) :: x21(3),x31(3)
  real(dp) :: bond_length
  !=====

  natom_read = SIZE(zatom_read(:))
  nghost_ = nghost_in
  ! here x_read contains all coordinates: first (natom_in-nprojectile) real atoms, then (nghost) ghost atoms
  ! and the last but not least (nprojectile) projectile (which is 0 or 1)

  ! xatom and zatom designate the physical nuclei that generate a Coulomb potential
  allocate(xatom(3,ncenter_nuclei))
  allocate(zatom(ncenter_nuclei))
  allocate(zvalence(ncenter_nuclei))
  ! xbasis and zbasis designate the basis centers and nature
  allocate(zbasis(ncenter_basis))
  allocate(xbasis(3,ncenter_basis))

  allocate(vel_nuclei(3,ncenter_nuclei))
  allocate(vel_basis(3,ncenter_basis))
  vel_nuclei(:,:) = 0.0_dp
  vel_basis(:,:) = 0.0_dp

  if( nprojectile /= 0 ) then
    vel_nuclei(:,ncenter_nuclei) = vel_projectile(:)
  endif

  if( excit_name == 'ION' .OR. excit_name == 'ANTIION' ) then
    vel_basis(:,ncenter_basis) = vel_projectile(:)
  endif

  allocate(force_nuc_nuc(3,ncenter_nuclei))
  ! For relaxation or dynamics only
  if( calculate_forces .OR. excit_name == "ION" .OR. excit_name == 'ANTIION' ) then
    !if( natom_in /= ncenter_nuclei .OR. natom_in /= ncenter_basis ) then
    !   call die('init_atoms: forces not implemented with ghosts or projectiles')
    !endif
    allocate(force(3,ncenter_nuclei))
    allocate(force_kin(3,ncenter_nuclei))
    allocate(force_nuc(3,ncenter_nuclei))
    allocate(force_har(3,ncenter_nuclei))
    allocate(force_exx(3,ncenter_nuclei))
    allocate(force_exc(3,ncenter_nuclei))
    allocate(force_ovp(3,ncenter_nuclei))
    allocate(force_hellfeyn(3,ncenter_nuclei))
  endif

 ! List of atoms is organized as follows:
 ! 1. physical atoms   :    nucleus | basis
 ! 2. ghost atoms      :      no    | basis
 ! 3. projectile       :    nucleus |   yes if 'ION' / no if 'NUCLEUS' or 'ANTINUCLEUS'
 !
 ! ncenter_nuclei contains the number of sites having a nucleus
 ! ncenter_basis  contains the number of sites having basis functions
 !
 xatom(:,1:ncenter_nuclei-nprojectile) = x_read(:,1:ncenter_nuclei-nprojectile)
 zatom(1:ncenter_nuclei-nprojectile)   = zatom_read(1:ncenter_nuclei-nprojectile)
 if( nprojectile /= 0 ) then
   xatom(:,ncenter_nuclei) = x_read(:,natom_read)
   zatom(ncenter_nuclei)   = zatom_read(natom_read)
 endif

 if( excit_name == "ANTINUCLEUS" .OR. excit_name == "ANTIION" ) then
   zatom(ncenter_nuclei) = -zatom(ncenter_nuclei)
 endif
 !
 ! In case of a projectile excitation, offer the possibility to tweak
 ! the charge of the projectile with an input variable
 ! Remember that the projectile always comes last in the atom list.
 if( excit_name == "NUCLEUS" .OR. excit_name == "ANTINUCLEUS" ) then
   zatom(ncenter_nuclei) = zatom(ncenter_nuclei) * projectile_charge_scaling
 end if

 !! Setup basis identity and centers
 ! Ghost atoms do not have a positive nucleus
 jcenter = 0
 do iatom=1,natom_read
   if( nucleus_wo_basis(iatom) ) cycle
   jcenter = jcenter + 1
   xbasis(:,jcenter) = x_read(:,iatom)
   zbasis(jcenter)   = NINT(zatom_read(iatom))
 enddo

 !
 ! Check for atoms too close
 do iatom=1,ncenter_nuclei
   do jatom=iatom+1,ncenter_nuclei
     if( NORM2( xatom(:,iatom)-xatom(:,jatom) ) < 0.2 ) then
       write(stdout,*) 'Atoms',iatom,jatom
       write(stdout,*) 'are closer than 0.2 bohr'
       call issue_warning('Some atoms are too close')
     endif
   enddo
 enddo

 !
 ! Find the covalent bonds based on simple distance criterium
 nbond = 0
 do iatom=1,ncenter_nuclei
   do jatom=iatom+1,ncenter_nuclei
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
 if( ncenter_nuclei > 2 ) then
   x21(:) = xatom(:,2) - xatom(:,1)
   do iatom=3,ncenter_nuclei
     x31(:) = xatom(:,iatom) - xatom(:,1)
     call cross_product(x21,x31,xnormal)
     if( NORM2(xnormal(:)) > tol_geom ) then
       xnormal(:) = xnormal(:) / NORM2(xnormal(:))
       linear=.FALSE.
       exit
     endif
   enddo
   if( .NOT. linear) then
     do iatom=1,ncenter_nuclei
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
  integer :: icenter
  !=====

  atoms_core_states=0
  do icenter=1,ncenter_nuclei
    atoms_core_states = atoms_core_states + element_core(zvalence(icenter),zatom(icenter))
  enddo

end function atoms_core_states


!=========================================================================
subroutine get_bondcenter(ibond,xbond)
  implicit none

  integer,intent(in)   :: ibond
  real(dp),intent(out) :: xbond(3)
  !=====
  integer :: icenter,jcenter,jbond
  !=====

  jbond = 0
  do icenter=1,ncenter_nuclei
    do jcenter=icenter+1,ncenter_nuclei
      if( NORM2( xatom(:,icenter)-xatom(:,jcenter) ) < 4.0 ) then
        jbond = jbond + 1
        if(jbond==ibond) then
          xbond(:) = 0.5_dp * ( xatom(:,icenter) - xatom(:,jcenter) )
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
  !=====

  xatom(:,iatom) = xposition(:)

end subroutine change_position_one_atom


!=========================================================================
subroutine change_basis_center_one_atom(iatom,xposition)
 implicit none
 integer,intent(in)   :: iatom
 real(dp),intent(in)  :: xposition(3)
!=====

 xbasis(:,iatom) = xposition(:)

end subroutine change_basis_center_one_atom


!=========================================================================
subroutine destroy_atoms()
  implicit none
  !=====

  if(ALLOCATED(zatom))          deallocate(zatom)
  if(ALLOCATED(zvalence))       deallocate(zvalence)
  if(ALLOCATED(zbasis))         deallocate(zbasis)
  if(ALLOCATED(xatom))          deallocate(xatom)
  if(ALLOCATED(xbasis))         deallocate(xbasis)
  if(ALLOCATED(force))          deallocate(force)
  if(ALLOCATED(force_nuc_nuc))  deallocate(force_nuc_nuc)
  if(ALLOCATED(force_kin))      deallocate(force_kin)
  if(ALLOCATED(force_nuc))      deallocate(force_nuc)
  if(ALLOCATED(force_har))      deallocate(force_har)
  if(ALLOCATED(force_exx))      deallocate(force_exx)
  if(ALLOCATED(force_exc))      deallocate(force_exc)
  if(ALLOCATED(force_ovp))      deallocate(force_ovp)
  if(ALLOCATED(force_hellfeyn)) deallocate(force_hellfeyn)

end subroutine destroy_atoms


!=========================================================================
subroutine relax_atoms(lbfgs_plan,etotal)
  use m_lbfgs
  implicit none

  type(lbfgs_state),intent(inout) :: lbfgs_plan
  real(dp),intent(in)             :: etotal
  !=====
  integer  :: info,iatom,idir
  real(dp) :: xnew(3,ncenter_nuclei)
  !=====

  xnew(:,:) = xatom(:,:)

  info = lbfgs_execute(lbfgs_plan,xnew,etotal,-force)

  ! Do not move the atoms by more than 0.20 bohr
  do iatom=1,ncenter_nuclei-nprojectile
    do idir=1,3
      if( ABS( xnew(idir,iatom) - xatom(idir,iatom) ) > 0.20_dp ) then
        xnew(idir,iatom) = xatom(idir,iatom) + SIGN( 0.20_dp , xnew(idir,iatom) - xatom(idir,iatom) )
      endif
    enddo
  enddo

  ! assume here that all atoms are both nuclei centers and basis centers  (no ghost, no projectile)
  xatom(:,:) = xnew(:,:)
  xbasis(:,:) = xnew(:,:)

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
  do iatom=1,ncenter_nuclei-nprojectile
    write(stdout,'(1x,a,i3,2x,a8,a,3(1x,f12.6),6x,3(1x,f12.6))') 'atom  ',iatom, &
                                                            element_name_long(zatom(iatom)),': ',  &
                                                            xatom(:,iatom),xatom(:,iatom)*bohr_A
  enddo

  if( nghost_ > 0 ) write(stdout,'(a)') ' == ghost list'
  do ighost=ncenter_nuclei-nprojectile+1,ncenter_nuclei-nprojectile+nghost_
    write(stdout,'(1x,a,i3,2x,a8,a,3(1x,f12.6),6x,3(1x,f12.6))') 'ghost ',iatom, &
                                            element_name_long(REAL(zbasis(ighost),dp)),': ',  &
                                            xbasis(:,ighost),xbasis(:,ighost)*bohr_A
  enddo
  if( nprojectile > 0 ) then
    write(stdout,'(a)') ' == projectile'
    write(stdout,'(1x,a,i3,2x,a8,a,3(1x,f12.6),6x,3(1x,f12.6))') 'atom  ',ncenter_nuclei, &
                                                            element_name_long(zatom(ncenter_nuclei)),': ',  &
                                                            xatom(:,ncenter_nuclei),xatom(:,ncenter_nuclei)*bohr_A
  endif


  write(stdout,'(1x,a,/)') '================================'

end subroutine output_positions


!=========================================================================
subroutine output_projectile_position()
  implicit none
  !=====
  !=====

  write(stdout,*)
  if( nprojectile > 0 ) then
    write(stdout,'(a)') ' === projectile position: ----------bohr---------------    |||   ------------- angstrom----------==='

    if( zatom(ncenter_nuclei) > 0 ) then
      write(stdout,'(1x,a,i3,2x,a2,a,3(1x,f12.6),6x,3(1x,f12.6))') 'atom  ',ncenter_nuclei, &
                                                            element_name(REAL(zatom(ncenter_nuclei),dp)),': ',  &
                                                            xatom(:,ncenter_nuclei),xatom(:,ncenter_nuclei)*bohr_A
    else
      write(stdout,'(1x,a,i3,2x,a4,a2,a,3(1x,f12.6),6x,3(1x,f12.6))') 'atom  ',ncenter_nuclei, &
                                                            'anti',element_name(REAL(zatom(ncenter_nuclei),dp)),': ',  &
                                                            xatom(:,ncenter_nuclei),xatom(:,ncenter_nuclei)*bohr_A
    endif

  endif

end subroutine output_projectile_position


!=========================================================================
subroutine nucleus_nucleus_energy(energy)
 implicit none

  real(dp),intent(out) :: energy
  !=====
  integer              :: icenter,jcenter
  !=====

  energy = 0.0_dp
  do icenter=1,ncenter_nuclei
    do jcenter=icenter+1,ncenter_nuclei
      energy = energy + zvalence(icenter) * zvalence(jcenter) / SQRT( SUM( (xatom(:,icenter) - xatom(:,jcenter))**2) )
    enddo
  enddo

end subroutine nucleus_nucleus_energy


!=========================================================================
subroutine nucleus_nucleus_force()
 !=====
  implicit none
  integer              :: icenter,jcenter
  !=====

  write(*,*) 'FBFB',ncenter_nuclei

  force_nuc_nuc(:,:) = 0.0_dp
  do icenter=1,ncenter_nuclei
    do jcenter=1,ncenter_nuclei
      if( icenter == jcenter ) cycle
      force_nuc_nuc(:,icenter) = force_nuc_nuc(:,icenter) &
                        + zvalence(icenter) * zvalence(jcenter) / ( SUM( (xatom(:,icenter) - xatom(:,jcenter))**2) )**1.50_dp &
                               * ( xatom(:,icenter) - xatom(:,jcenter) )
    enddo
  enddo
  write(*,*) 'FBFB1',force_nuc_nuc(:,1)
  write(*,*) 'FBFB2',force_nuc_nuc(:,2)

end subroutine nucleus_nucleus_force


!=========================================================================
subroutine find_inversion()
  implicit none
  !=====
  integer  :: icenter,jcenter
  logical  :: found
  real(dp) :: xtmp(3)
  !=====

  xcenter(:) = 0.0_dp
  do icenter=1,ncenter_nuclei-nprojectile
    xcenter(:) = xcenter(:) + xatom(:,icenter) / REAL(ncenter_nuclei-nprojectile,dp)
  enddo

  do icenter=1,ncenter_nuclei-nprojectile
    xtmp(:) = 2.0_dp * xcenter(:) - xatom(:,icenter)
    found = .FALSE.
    do jcenter=1,ncenter_nuclei-nprojectile
      if( NORM2( xtmp(:) - xatom(:,jcenter) ) < tol_geom ) then
        if( ABS(zatom(icenter)-zatom(jcenter)) < tol_geom ) found = .TRUE.
        exit
      endif
    enddo
    inversion = inversion .AND. found
  enddo

end subroutine find_inversion


!=========================================================================
function same_element(icenter,jcenter)
  implicit none
  integer,intent(in) :: icenter,jcenter
  logical :: same_element
  !=====
  !=====

  same_element = ( ABS( zatom(icenter) - zatom(jcenter) ) < 1.0e-5 )

end function same_element



end module m_atoms
!=========================================================================
