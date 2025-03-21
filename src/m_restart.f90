!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the reading / writing methods for the RESTART files
!
!=========================================================================
#include "molgw.h"
module m_restart
  use m_definitions
  use m_timing
  use m_mpi
  use m_inputparam
  use m_atoms
  use m_basis_set
  use m_hamiltonian_tools, only: dump_out_occupation
  use m_hamiltonian_onebody, only: setup_overlap_mixedbasis, setup_overlap
  use m_linear_algebra, only: invert


  !
  ! Restart file types
  integer, parameter ::           NO_RESTART = 0
  integer, parameter ::        SMALL_RESTART = 1
  integer, parameter ::          BIG_RESTART = 2
  integer, parameter ::        BASIS_RESTART = 3



contains


!=========================================================================
subroutine write_restart(restart_type, basis, occupation, c_matrix, energy, hamiltonian_fock)
  implicit none

  integer, intent(in)           :: restart_type
  type(basis_set), intent(in)   :: basis
  real(dp), intent(in)          :: occupation(:, :), energy(:, :)
  real(dp), intent(in)          :: c_matrix(:, :, :)
  real(dp), optional, intent(in) :: hamiltonian_fock(:, :, :)
  !=====
  integer                    :: nstate
  integer, parameter          :: restart_version=201609
  integer                    :: restartfile
  integer                    :: ispin, istate, ibf, nstate_local
  !=====

  !
  ! Only the proc iomaster writes down the RESTART file
  if( .NOT. is_iomaster) return

  call start_clock(timing_restart_file)

  nstate = SIZE(occupation, DIM=1)

  select case(restart_type)
  case(SMALL_RESTART)
    write(stdout, '(/,a)') ' Writing a small RESTART file'
  case(BIG_RESTART)
    write(stdout, '(/,a)') ' Writing a big RESTART file'
  case default
    call die('write_restart: bug')
  end select

  if( restart_type == BIG_RESTART .AND. .NOT. PRESENT(hamiltonian_fock) ) then
    call die('write_restart: an input hamiltonian_fock is needed for a BIG restart')
  endif


  open(newunit=restartfile, file='RESTART', form='unformatted', action='write')

  ! An integer to ensure backward compatibility in the future
  write(restartfile) restart_version
  ! RESTART file type SMALL_RESTART=1 or BIG_RESTART=2
  write(restartfile) restart_type
  ! Atomic structure
  write(restartfile) ncenter_nuclei
  write(restartfile) zatom(1:ncenter_nuclei)
  write(restartfile) xatom(:, 1:ncenter_nuclei)
  ! Calculation type
  write(restartfile) calc_type%scf_name
  ! Basis set
  call write_basis_set(restartfile, basis)
  ! Spin channels
  write(restartfile) nspin
  ! Nstate
  write(restartfile) nstate
  ! Occupations
  write(restartfile) occupation(:, :)
  ! Eigen energies
  write(restartfile) energy(:, :)

  ! Number of states written down in the RESTART file
  if( restart_type == SMALL_RESTART ) then
    ! Identify the highest occupied state in order to
    ! save I-O in SMALL_RESTART
    nstate_local = 0
    do ispin=1, nspin
      do istate=1, nstate
        if( ANY( occupation(istate, :) > completely_empty ) ) nstate_local = istate
      enddo
    enddo
  else
    ! Or write down all the states in BIG_RESTART
    nstate_local = nstate
  endif

  write(restartfile) nstate_local

  ! Wavefunction coefficients C
  do ispin=1, nspin
    do istate=1, nstate_local
      write(restartfile) c_matrix(:, istate, ispin)
    enddo
  enddo

  if(restart_type == BIG_RESTART) then

    do ispin=1, nspin
      do ibf=1, basis%nbf
        write(restartfile) hamiltonian_fock(:, ibf, ispin)
      enddo
    enddo

  endif

  close(restartfile)
  call stop_clock(timing_restart_file)

end subroutine write_restart


!=========================================================================
subroutine read_restart(restart_type, restart_filename, basis, occupation, c_matrix, energy, hamiltonian_fock)
  implicit none

  integer, intent(out)           :: restart_type
  character(len=*), intent(in)   :: restart_filename
  type(basis_set), intent(in)    :: basis
  real(dp), intent(inout)        :: occupation(:, :)
  real(dp), intent(out)          :: c_matrix(:, :, :), energy(:, :)
  real(dp), optional, intent(out) :: hamiltonian_fock(:, :, :)
  !=====
  integer                    :: restartfile
  integer                    :: ispin, istate, ibf, nstate_local
  logical                    :: file_exists, same_scf, same_basis, same_geometry, same_spin
  integer                    :: nstate
  integer                    :: restart_version_read
  integer                    :: restart_type_read
  character(len=100)         :: scf_name_read
  integer                    :: natom_read
  real(dp), allocatable       :: zatom_read(:), x_read(:, :)
  type(basis_set)            :: basis_read
  integer                    :: nspin_read
  integer                    :: nstate_read
  real(dp), allocatable       :: occupation_read(:, :)
  real(dp), allocatable       :: energy_read(:, :)
  real(dp), allocatable       :: c_matrix_read(:, :, :)
  real(dp), allocatable       :: overlapm1(:, :)
  real(dp), allocatable       :: overlap_mixedbasis(:, :)
  !=====

  inquire(file=restart_filename, exist=file_exists)
  if(.NOT. file_exists) then
    write(stdout, '(/,1x,a,1x,a)') TRIM(restart_filename), 'file not found'
    restart_type = NO_RESTART
    return
  endif

  nstate = SIZE(occupation, DIM=1)

  open(newunit=restartfile, file=restart_filename, form='unformatted', status='old', action='read')


  ! An integer to ensure backward compatibility in the future
  read(restartfile) restart_version_read

  if( restart_version_read == 201602 ) then
    write(stdout, '(1x,a,i8)') 'Old RESTART file found. Version: ', restart_version_read
    call issue_warning('RESTART file: Old version. Backward compatibility is not ensured. Skipping the reading')
    restart_type = NO_RESTART
    close(restartfile)
    return
  endif

  if( restart_version_read /= 201609 ) then
    call issue_warning('RESTART file: Version not readable. Skipping the reading')
    restart_type = NO_RESTART
    close(restartfile)
    return
  endif


  ! RESTART file type SMALL_RESTART=1 or BIG_RESTART=2
  read(restartfile) restart_type_read
  if( restart_type_read /= SMALL_RESTART .AND. restart_type_read /= BIG_RESTART ) then
    call issue_warning('RESTART file: Type not readable. Skipping the reading')
    restart_type = NO_RESTART
    close(restartfile)
    return
  endif
  restart_type = restart_type_read

  !
  ! Input keyword ignore_bigrestart enforces a small_restart
  if( ignore_bigrestart_ ) then
    restart_type = SMALL_RESTART
  endif


  ! Atomic structure
  read(restartfile) natom_read
  allocate(zatom_read(natom_read), x_read(3, natom_read))
  read(restartfile) zatom_read(1:natom_read)
  read(restartfile) x_read(:, 1:natom_read)
  if( natom_read /= ncenter_nuclei  &
  .OR. ANY( ABS( zatom_read(1:MIN(natom_read, ncenter_nuclei)) - zatom(1:MIN(natom_read, ncenter_nuclei)) ) > 1.0e-5_dp ) &
  .OR. ANY( ABS(   x_read(:, 1:MIN(natom_read, ncenter_nuclei)) - xatom(:, 1:MIN(natom_read, ncenter_nuclei)) ) > 1.0e-5_dp ) ) then
    same_geometry = .FALSE.
    call issue_warning('RESTART file: Geometry has changed')
  else
    same_geometry = .TRUE.
  endif
  deallocate(zatom_read, x_read)


  ! Calculation type
  read(restartfile) scf_name_read
  same_scf = ( TRIM(scf_name_read) == TRIM(calc_type%scf_name) )
  if( .NOT. same_scf) then
    call issue_warning('RESTART file: SCF type has changed')
    restart_type = SMALL_RESTART
  endif


  ! Basis set
  call read_basis_set(restartfile, basis_read)
  same_basis = compare_basis_set(basis, basis_read)
  if( .NOT. same_basis) then
    call issue_warning('RESTART file: Basis set has changed')
    restart_type = SMALL_RESTART
  endif
  if( basis%gaussian_type /= basis_read%gaussian_type ) then
    write(stdout, *) 'The basis type (cartesian or pure) cannot be changed when restarting from a previous calculation'
    call die('Erase the RESTART file or change the keyword gaussian_type and start the calculation over')
  endif


  ! Spin channels
  read(restartfile) nspin_read
  same_spin = ( nspin == nspin_read )
  if( .NOT. same_spin ) then
    call issue_warning('RESTART file: Number of spin channels has changed')
    restart_type = SMALL_RESTART
  endif


  ! Nstate
  read(restartfile) nstate_read
  if( nstate /= nstate_read ) then
    call issue_warning('RESTART file: Number of states has changed')
    restart_type = SMALL_RESTART
  endif


  ! Occupations
  allocate(occupation_read(nstate_read, nspin_read))
  read(restartfile) occupation_read(:, :)
  if( ANY( ABS( occupation_read(1:MIN(nstate_read, nstate), 1) - occupation(1:MIN(nstate_read, nstate), 1) )  > 1.0e-5_dp ) ) then
    if( temperature > 1.0e-8_dp) then
      occupation(1:MIN(nstate_read, nstate), 1) = occupation_read(1:MIN(nstate_read, nstate), 1)
      occupation(1:MIN(nstate_read, nstate), nspin) = occupation_read(1:MIN(nstate_read, nstate), nspin_read)
      write(stdout, '(1x,a)') "Reading occupations from a RESTART file"
      call dump_out_occupation('=== Occupations ===', occupation)
    else
      call issue_warning('RESTART file: Occupations have changed')
      restart_type = SMALL_RESTART
    endif
  endif
  deallocate(occupation_read)


  ! Eigen energies
  allocate(energy_read(nstate_read, nspin_read))
  read(restartfile) energy_read(:, :)
  energy(:, :) = 1000.0_dp
  energy(1:MIN(nstate, nstate_read), 1) = energy_read(1:MIN(nstate, nstate_read), 1)
  energy(1:MIN(nstate, nstate_read), nspin) = energy_read(1:MIN(nstate, nstate_read), nspin_read)
  deallocate(energy_read)


  ! Number of states written down in the RESTART file
  read(restartfile) nstate_local


  ! Wavefunction coefficients C
  allocate(c_matrix_read(basis_read%nbf, nstate_local, nspin_read))
  do ispin=1, nspin_read
    do istate=1, nstate_local
      read(restartfile) c_matrix_read(:, istate, ispin)
    enddo
  enddo

  if( same_basis ) then
    c_matrix(:, :, :) = 0.0_dp
    do istate=1, MIN(nstate_local, nstate)
      c_matrix(1:MIN(basis_read%nbf, basis%nbf), istate, 1) &
          = c_matrix_read(1:MIN(basis_read%nbf, basis%nbf), istate, 1)
    enddo
    do istate=1, MIN(nstate_local, nstate)
      c_matrix(1:MIN(basis_read%nbf, basis%nbf), istate, nspin) &
          = c_matrix_read(1:MIN(basis_read%nbf, basis%nbf), istate, nspin_read)
    enddo

    ! Fill the rest of the array with identity
    if( nstate_local < nstate ) then
      do istate=nstate_local+1, nstate
        c_matrix(istate, istate, :) = 1.0_dp
      enddo
    endif
  else
    allocate(overlapm1(basis%nbf, basis%nbf))
    allocate(overlap_mixedbasis(basis%nbf, basis_read%nbf))

    ! Calculate the overlap matrix of the final basis set
    call setup_overlap(basis, overlapm1)

    ! Invert the overlap of the final basis set
    call invert(overlapm1)

    ! Get the scalar products between the old and the new basis sets
    ! Be aware: this is a rectangular matrix
    call setup_overlap_mixedbasis(basis, basis_read, overlap_mixedbasis)
    c_matrix(:, 1:nstate_local, 1) = MATMUL(overlapm1(:, :), &
                                         MATMUL(overlap_mixedbasis(:, :) , c_matrix_read(:, 1:nstate_local, 1) ) )
    c_matrix(:, 1:nstate_local, nspin) = MATMUL(overlapm1(:, :), &
                                         MATMUL(overlap_mixedbasis(:, :) , c_matrix_read(:, 1:nstate_local, nspin_read) ) )
    ! Fill the rest of the array with identity
    do istate=nstate_local+1, nstate
      c_matrix(istate, istate, :) = 1.0_dp
    enddo
    deallocate(overlapm1, overlap_mixedbasis)

    restart_type = BASIS_RESTART
    close(restartfile)
    return

  endif


  if( ignore_bigrestart_ .OR. restart_type_read == SMALL_RESTART .OR. .NOT. PRESENT(hamiltonian_fock) &
      .OR. .NOT. same_spin ) then

    close(restartfile)
    return

  else

    if( same_basis ) then

      do ispin=1, nspin_read
        do ibf=1, basis_read%nbf
          read(restartfile) hamiltonian_fock(:, ibf, ispin)
        enddo
      enddo

      if( same_geometry ) then
        restart_type = BIG_RESTART
      endif
      close(restartfile)
      return

    endif


  endif

  ! the code should never reach that point.
  call die('Internal error in the read_restart subroutine')


end subroutine read_restart


!=========================================================================
end module m_restart
!=========================================================================
