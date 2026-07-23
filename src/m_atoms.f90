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
  use m_warning, only: die, issue_warning
  use m_elements
  use m_linear_algebra, only: cross_product, inverse_3x3_matrix, determinant_3x3_matrix

  implicit none

  real(dp), parameter, private     :: tol_geom=1.0e-5_dp

  integer, public                 :: ncenter_basis
  integer, public                 :: ncenter_nuclei
  integer, protected              :: nbond
  integer, public                 :: nprojectile
  integer, private                :: nghost_

  real(dp), allocatable, public    :: zvalence(:)
  real(dp), allocatable, protected :: zatom(:)
  integer, allocatable, protected  :: zbasis(:)

  real(dp), allocatable, protected :: xatom(:, :)
  real(dp), allocatable, protected :: xbasis(:, :)
  real(dp), allocatable, protected :: vel_nuclei(:, :)
  real(dp), allocatable, protected :: vel_basis(:, :)
  real(dp), allocatable, public    :: force(:, :)
  real(dp), public                :: force_projectile(3) = 0.0_dp
  real(dp), public                :: force_projectile_nonconserv(3) = 0.0_dp

  ! See if we keep these arrays in the long-term
  real(dp), allocatable, public    :: force_nuc_nuc(:, :)
  real(dp), allocatable, public    :: force_kin(:, :)
  real(dp), allocatable, public    :: force_nuc(:, :)
  real(dp), allocatable, public    :: force_har(:, :)
  real(dp), allocatable, public    :: force_exx(:, :)
  real(dp), allocatable, public    :: force_exc(:, :)
  real(dp), allocatable, public    :: force_ovp(:, :)
  real(dp), allocatable, public    :: force_hellfeyn(:, :)

  logical, protected              :: inversion=.TRUE.
  logical, protected              :: linear=.TRUE.
  logical, protected              :: planar=.TRUE.
  real(dp), protected             :: xcenter(3)
  real(dp), protected             :: xnormal(3)

  logical, protected              :: pbc_ = .FALSE.
  real(dp), protected             :: aprim(3, 3) = 0.0_dp
  real(dp), protected             :: aprim_inv(3, 3) = 1000.0_dp
  real(dp), protected             :: bprim(3, 3) = 1000.0_dp
  real(dp), protected             :: volume
  real(dp), protected             :: recip_volume
  real(dp), protected             :: minimal_image_distance


contains


!=========================================================================
subroutine init_atoms(natom_in, nghost_in, nucleus_wo_basis, zatom_read, x_read, vel_projectile, &
                      calculate_forces, excit_name, projectile_charge_scaling)

  integer, intent(in)  :: natom_in, nghost_in
  logical, intent(in)  :: nucleus_wo_basis(:)
  real(dp), intent(in) :: zatom_read(:), x_read(:, :)
  real(dp), intent(in) :: vel_projectile(3)
  logical, intent(in)  :: calculate_forces
  character(len=12), intent(in)   :: excit_name
  real(dp), intent(in) :: projectile_charge_scaling
  !=====
  integer  :: natom_read
  integer  :: iatom, jatom, jcenter
  real(dp) :: x21(3), x31(3)
  real(dp) :: bond_length
  !=====

  natom_read = SIZE(zatom_read(:))
  nghost_ = nghost_in
  ! here x_read contains all coordinates: first (natom_in-nprojectile) real atoms, then (nghost) ghost atoms
  ! and the last but not least (nprojectile) projectile (which is 0 or 1)

  ! xatom and zatom designate the physical nuclei that generate a Coulomb potential
  allocate(xatom(3, ncenter_nuclei))
  allocate(zatom(ncenter_nuclei))
  allocate(zvalence(ncenter_nuclei))
  ! xbasis and zbasis designate the basis centers and nature
  allocate(zbasis(ncenter_basis))
  allocate(xbasis(3, ncenter_basis))

  allocate(vel_nuclei(3, ncenter_nuclei))
  allocate(vel_basis(3, ncenter_basis))
  vel_nuclei(:, :) = 0.0_dp
  vel_basis(:, :) = 0.0_dp

  if( nprojectile /= 0 ) then
    vel_nuclei(:, ncenter_nuclei) = vel_projectile(:)
  end if

  if( excit_name == 'ION' .OR. excit_name == 'ANTIION' ) then
    vel_basis(:, ncenter_basis) = vel_projectile(:)
  end if

  ! Allocate force arrays (tiny memory footprint anyway)
  allocate(force(3, ncenter_nuclei))
  allocate(force_nuc_nuc(3, ncenter_nuclei))
  allocate(force_kin(3, ncenter_nuclei))
  allocate(force_nuc(3, ncenter_nuclei))
  allocate(force_har(3, ncenter_nuclei))
  allocate(force_exx(3, ncenter_nuclei))
  allocate(force_exc(3, ncenter_nuclei))
  allocate(force_ovp(3, ncenter_nuclei))
  allocate(force_hellfeyn(3, ncenter_nuclei))

  ! List of atoms is organized as follows:
  ! 1. physical atoms   :    nucleus | basis
  ! 2. ghost atoms      :      no    | basis
  ! 3. projectile       :    nucleus |   yes if 'ION' / no if 'NUCLEUS' or 'ANTINUCLEUS'
  !
  ! ncenter_nuclei contains the number of sites having a nucleus
  ! ncenter_basis  contains the number of sites having basis functions
  !
  xatom(:, 1:ncenter_nuclei-nprojectile) = x_read(:, 1:ncenter_nuclei-nprojectile)
  zatom(1:ncenter_nuclei-nprojectile)   = zatom_read(1:ncenter_nuclei-nprojectile)
  if( nprojectile /= 0 ) then
    xatom(:, ncenter_nuclei) = x_read(:, natom_read)
    zatom(ncenter_nuclei)   = zatom_read(natom_read)
  end if

  if( excit_name == "ANTINUCLEUS" .OR. excit_name == "ANTIION" ) then
    zatom(ncenter_nuclei) = -zatom(ncenter_nuclei)
  end if
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
  do iatom=1, natom_read
    if( nucleus_wo_basis(iatom) ) cycle
    jcenter = jcenter + 1
    xbasis(:, jcenter) = x_read(:, iatom)
    zbasis(jcenter)   = NINT(zatom_read(iatom))
  end do

  !
  ! Check for atoms too close
  do iatom=1, ncenter_nuclei
    do jatom=iatom+1, ncenter_nuclei
      if( NORM2( xatom(:, iatom)-xatom(:, jatom) ) < 0.2 ) then
        write(stdout, *) 'Atoms', iatom, jatom
        write(stdout, *) 'are closer than 0.2 bohr'
        call issue_warning('Some atoms are too close')
      end if
    end do
  end do

  !
  ! Find the covalent bonds based on simple distance criterium
  nbond = 0
  do iatom=1, ncenter_nuclei
    do jatom=iatom+1, ncenter_nuclei
      bond_length =  element_covalent_radius(NINT(zatom(iatom))) + element_covalent_radius(NINT(zatom(jatom)))
      if( NORM2( xatom(:, iatom)-xatom(:, jatom) ) <  1.2_dp * bond_length  ) then
        nbond = nbond + 1
      end if
    end do
  end do

  !
  ! Does the molecule have inversion symmetry?
  call find_inversion()

  !
  ! Is the molecule linear, planar?
  if( ncenter_nuclei > 2 ) then
    x21(:) = xatom(:, 2) - xatom(:, 1)
    do iatom=3, ncenter_nuclei
      x31(:) = xatom(:, iatom) - xatom(:, 1)
      call cross_product(x21, x31, xnormal)
      if( NORM2(xnormal(:)) > tol_geom ) then
        xnormal(:) = xnormal(:) / NORM2(xnormal(:))
        linear=.FALSE.
        exit
      end if
    end do
    if( .NOT. linear) then
      do iatom=1, ncenter_nuclei
        if( ABS(DOT_PRODUCT( xatom(:, iatom) , xnormal(:) )) > tol_geom ) planar=.FALSE.
      end do
    else
      planar=.FALSE.
    end if
  else
    ! Molecule is linear
    ! Set planar to FALSE for safety
    linear=.TRUE.
    planar=.FALSE.
  end if


end subroutine init_atoms


!=========================================================================
function atoms_core_states()

  integer :: atoms_core_states
  !=====
  integer :: icenter
  !=====

  atoms_core_states=0
  do icenter=1, ncenter_nuclei
    atoms_core_states = atoms_core_states + element_core(zvalence(icenter), zatom(icenter))
  end do

end function atoms_core_states


!=========================================================================
function atoms_core_states_gaussian() result(atoms_core_states)

  integer :: atoms_core_states
  !=====
  integer :: icenter
  !=====

  atoms_core_states=0
  do icenter=1, ncenter_nuclei
    atoms_core_states = atoms_core_states + element_core_gaussian(zvalence(icenter), zatom(icenter))
  end do

end function atoms_core_states_gaussian


!=========================================================================
subroutine get_bondcenter(ibond, xbond)

  integer, intent(in)   :: ibond
  real(dp), intent(out) :: xbond(3)
  !=====
  integer :: icenter, jcenter, jbond
  !=====

  jbond = 0
  do icenter=1, ncenter_nuclei
    do jcenter=icenter+1, ncenter_nuclei
      if( NORM2( xatom(:, icenter)-xatom(:, jcenter) ) < 4.0 ) then
        jbond = jbond + 1
        if(jbond==ibond) then
          xbond(:) = 0.5_dp * ( xatom(:, icenter) - xatom(:, jcenter) )
          return
        end if
      end if
    end do
  end do

end subroutine get_bondcenter


!=========================================================================
subroutine change_position_one_atom(iatom, xposition)
  integer, intent(in)   :: iatom
  real(dp), intent(in)  :: xposition(3)
  !=====
  !=====

  xatom(:, iatom) = xposition(:)

end subroutine change_position_one_atom


!=========================================================================
subroutine change_basis_center_one_atom(iatom, xposition)
  integer, intent(in)   :: iatom
  real(dp), intent(in)  :: xposition(3)
  !=====

  xbasis(:, iatom) = xposition(:)

end subroutine change_basis_center_one_atom


!=========================================================================
subroutine destroy_atoms()
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
subroutine relax_atoms(lbfgs_plan, etotal)
  use m_lbfgs

  type(lbfgs_state), intent(inout) :: lbfgs_plan
  real(dp), intent(in)             :: etotal
  !=====
  integer  :: info, iatom, idir
  real(dp) :: xnew(3, ncenter_nuclei)
  !=====

  xnew(:, :) = xatom(:, :)

  info = lbfgs_execute(lbfgs_plan, xnew, etotal, -force)

  ! Do not move the atoms by more than 0.20 bohr
  do iatom=1, ncenter_nuclei-nprojectile
    do idir=1, 3
      if( ABS( xnew(idir, iatom) - xatom(idir, iatom) ) > 0.20_dp ) then
        xnew(idir, iatom) = xatom(idir, iatom) + SIGN( 0.20_dp , xnew(idir, iatom) - xatom(idir, iatom) )
      end if
    end do
  end do

  ! assume here that all atoms are both nuclei centers and basis centers  (no ghost, no projectile)
  xatom(:, :) = xnew(:, :)
  xbasis(:, :) = xnew(:, :)

end subroutine relax_atoms


!=========================================================================
subroutine output_positions()
  !=====
  integer :: ighost, iatom
  !=====

  write(stdout, '(/,1x,a)') '================================'
  write(stdout, *) '      Atom list'
  write(stdout, *) '                       bohr                                        angstrom'
  do iatom=1, ncenter_nuclei-nprojectile
    write(stdout, '(1x,a,i3,2x,a8,a,3(1x,f12.6),6x,3(1x,f12.6))') 'atom  ', iatom, &
                                                            element_name_long(zatom(iatom)), ': ',  &
                                                            xatom(:, iatom), xatom(:, iatom)*bohr_A
  end do

  if( nghost_ > 0 ) write(stdout, '(a)') ' == ghost list'
  do ighost=ncenter_nuclei-nprojectile+1, ncenter_nuclei-nprojectile+nghost_
    write(stdout, '(1x,a,i3,2x,a8,a,3(1x,f12.6),6x,3(1x,f12.6))') 'ghost ', iatom, &
                                            element_name_long(zbasis(ighost)), ': ',  &
                                            xbasis(:, ighost), xbasis(:, ighost)*bohr_A
  end do
  if( nprojectile > 0 ) then
    write(stdout, '(a)') ' == projectile'
    write(stdout, '(1x,a,i3,2x,a8,a,3(1x,f12.6),6x,3(1x,f12.6))') 'atom  ', ncenter_nuclei, &
                                                            element_name_long(zatom(ncenter_nuclei)), ': ',  &
                                                            xatom(:, ncenter_nuclei), xatom(:, ncenter_nuclei)*bohr_A
  end if


  write(stdout, '(1x,a,/)') '================================'

end subroutine output_positions


!=========================================================================
subroutine output_projectile_position()
  !=====
  !=====

  write(stdout, *)
  if( nprojectile > 0 ) then
    write(stdout, '(a)') ' === projectile position: ----------bohr---------------    |||   ------------- angstrom----------==='

    if( zatom(ncenter_nuclei) > 0 ) then
      write(stdout, '(1x,a,i3,2x,a2,a,3(1x,f12.6),6x,3(1x,f12.6))') 'atom  ', ncenter_nuclei, &
                                                            element_name(zatom(ncenter_nuclei)), ': ',  &
                                                            xatom(:, ncenter_nuclei), xatom(:, ncenter_nuclei)*bohr_A
    else
      write(stdout, '(1x,a,i3,2x,a4,a2,a,3(1x,f12.6),6x,3(1x,f12.6))') 'atom  ', ncenter_nuclei, &
                                                            'anti', element_name(zatom(ncenter_nuclei)), ': ',  &
                                                            xatom(:, ncenter_nuclei), xatom(:, ncenter_nuclei)*bohr_A
    end if

  end if

end subroutine output_projectile_position


!=========================================================================
subroutine nucleus_nucleus_energy(energy)

  real(dp), intent(out) :: energy
  !=====
  integer              :: icenter, jcenter
  !=====

  energy = 0.0_dp
  do icenter=1, ncenter_nuclei
    do jcenter=icenter+1, ncenter_nuclei
      energy = energy + zvalence(icenter) * zvalence(jcenter) / SQRT( SUM( (xatom(:, icenter) - xatom(:, jcenter))**2) )
    end do
  end do

end subroutine nucleus_nucleus_energy


!=========================================================================
subroutine nucleus_nucleus_force()
  !=====
  integer              :: icenter, jcenter
  !=====

  force_nuc_nuc(:, :) = 0.0_dp
  do icenter=1, ncenter_nuclei
    do jcenter=1, ncenter_nuclei
      if( icenter == jcenter ) cycle
      force_nuc_nuc(:, icenter) = force_nuc_nuc(:, icenter) &
                        + zvalence(icenter) * zvalence(jcenter) / ( SUM( (xatom(:, icenter) - xatom(:, jcenter))**2) )**1.50_dp &
                               * ( xatom(:, icenter) - xatom(:, jcenter) )
    end do
  end do

end subroutine nucleus_nucleus_force


!=========================================================================
subroutine find_inversion()
  !=====
  integer  :: icenter, jcenter
  logical  :: found
  real(dp) :: xtmp(3)
  !=====

  xcenter(:) = 0.0_dp
  do icenter=1, ncenter_nuclei-nprojectile
    xcenter(:) = xcenter(:) + xatom(:, icenter) / REAL(ncenter_nuclei-nprojectile, dp)
  end do

  do icenter=1, ncenter_nuclei-nprojectile
    xtmp(:) = 2.0_dp * xcenter(:) - xatom(:, icenter)
    found = .FALSE.
    do jcenter=1, ncenter_nuclei-nprojectile
      if( NORM2( xtmp(:) - xatom(:, jcenter) ) < tol_geom ) then
        if( ABS(zatom(icenter)-zatom(jcenter)) < tol_geom ) found = .TRUE.
        exit
      end if
    end do
    inversion = inversion .AND. found
  end do

end subroutine find_inversion


!=========================================================================
function same_element(icenter, jcenter)
  integer, intent(in) :: icenter, jcenter
  logical :: same_element
  !=====
  !=====

  same_element = ( ABS( zatom(icenter) - zatom(jcenter) ) < 1.0e-5 )

end function same_element


!=========================================================================
subroutine setup_periodicity_vectors(length_unit, a1, a2, a3)
  character(len=*), intent(in) :: length_unit
  real(dp), intent(in) :: a1(3), a2(3), a3(3)
  !=====
  real(dp) :: length_factor, x(3)
  integer :: i1, i2, i3
  !=====

  select case(TRIM(length_unit))
  case('A', 'ANGSTROM')
    length_factor = 1.0_dp / bohr_A
  case('BOHR', 'AU','A.U','A.U.')
    length_factor = 1.0_dp
  case default
    call die('units for lengths in input file not understood')
  end select

  if( NORM2(a1) > 1.0e-1_dp ) then
    aprim(:, 1) = a1(:) * length_factor
  end if

  if( NORM2(a2) > 1.0e-1_dp ) then
    aprim(:, 2) = a2(:) * length_factor
  end if

  if( NORM2(a3) > 1.0e-1_dp ) then
    aprim(:, 3) = a3(:) * length_factor
  end if

  pbc_ = NORM2(aprim) > 1.0e-1_dp


  if( pbc_ ) then

    volume = determinant_3x3_matrix(aprim)
    ! Calculate the reciprocal lattice vectors
    if( volume < 0.0_dp ) then
      call die('setup_periodicity_vectors: periodic vectors are in indirect order. Swap two of them to recover direct order.')
    end if
    aprim_inv(:, :) = inverse_3x3_matrix(aprim)
    bprim(:, :) = 2.0_dp * pi * TRANSPOSE(aprim_inv)
    recip_volume = determinant_3x3_matrix(bprim)

    minimal_image_distance = 9999.9_dp
    do i3=-5, 5
      do i2=-5, 5
        do i1=-5, 5
          if( i1 == 0 .AND. i2 == 0 .AND. i3 == 0 ) cycle
          x(:) = MATMUL(aprim, [i1, i2, i3])
          minimal_image_distance = MIN(minimal_image_distance, NORM2(x(:)) )
        end do
      end do
    end do

    write(stdout, '(/,1x,a,/)') 'Periodic boundary conditions are switched on'
    write(stdout, '(1x,a)') 'Primitive cell vectors (bohr):'
    write(stdout, '(1x,3(f12.6,1x))')  aprim(:, 1)
    write(stdout, '(1x,3(f12.6,1x))')  aprim(:, 2)
    write(stdout, '(1x,3(f12.6,1x))')  aprim(:, 3)
    write(stdout, '(1x,a,1x,f12.6)') 'Volume (bohr**3):             ', volume
    write(stdout, '(1x,a,1x,f12.6)') 'Minimal image distance (bohr):', minimal_image_distance
    write(stdout, '(1x,a)') 'Reciprocal lattice vectors (bohr**-1):'
    write(stdout, '(1x,3(f12.6,1x))')  bprim(:, 1)
    write(stdout, '(1x,3(f12.6,1x))')  bprim(:, 2)
    write(stdout, '(1x,3(f12.6,1x))')  bprim(:, 3)
    write(stdout, '(1x,a,1x,f12.6)') 'Volume (bohr**-3):', recip_volume
    write(stdout, '(1x,a)')  '==============================================='

 end if


end subroutine setup_periodicity_vectors


end module m_atoms
!=========================================================================
