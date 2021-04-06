!=========================================================================
! This file is part of MOLGW.
!
! Copyright (C) 2010-2020  Fabien Bruneval
! MOLGW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! MOLGW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with MOLGW.  If not, see <http://www.gnu.org/licenses/>.
!=========================================================================
!
! This is the main of MOLGW
!
! It consists of 3 parts:
!   1. Initialization
!   2. SCF cycles
!   3. Post-processings: GW, MP2, CI, TDDFT etc.
!
!=========================================================================
program molgw
  use m_definitions
  use m_timing
  use m_warning
  use m_memory
  use m_scalapack
  use m_inputparam
  use m_scf
  use m_atoms
  use m_ecp
  use m_gaussian
  use m_basis_set
  use m_lbfgs
  use m_eri
  use m_eri_calculate
  use m_eri_ao_mo
  use m_dft_grid
  use m_spectral_function
  use m_hamiltonian_onebody
  use m_selfenergy_tools
  use m_selfenergy_evaluation
  use m_scf_loop
  use m_tddft_propagator
  use m_tddft_variables
  use m_virtual_orbital_space
  use m_ci
  use m_dm_analysis
  use m_dm_mbpt
  use m_restart
  use m_multipole
  use m_io
  use m_fourier_quadrature
  implicit none

  !=====
  type(basis_set)            :: basis
  type(basis_set)            :: auxil_basis
  type(spectral_function)    :: wpol
  type(lbfgs_state)          :: lbfgs_plan
  type(energy_contributions) :: en_gks,en_mbpt
  integer                 :: restart_type
  integer                 :: nstate
  integer                 :: istep
  logical                 :: is_restart,is_big_restart,is_basis_restart
  logical                 :: restart_tddft_is_correct = .TRUE.
  logical                 :: is_converged
  real(dp)                :: erpa_tmp,egw_tmp
  real(dp),allocatable    :: hamiltonian_tmp(:,:,:)
  real(dp),allocatable    :: hamiltonian_kinetic(:,:)
  real(dp),allocatable    :: hamiltonian_nucleus(:,:)
  real(dp),allocatable    :: hamiltonian_fock(:,:,:)
  real(dp),allocatable    :: s_matrix(:,:)
  real(dp),allocatable    :: x_matrix(:,:)
  real(dp),allocatable    :: c_matrix(:,:,:)
  real(dp),allocatable    :: energy(:,:)
  real(dp),allocatable    :: occupation(:,:)
  real(dp),allocatable    :: exchange_m_vxc(:,:,:)
  !=====

  !
  !
  ! Part 1 / 3 : Initialization
  !
  !

  call init_mpi_world()

  call init_scalapack()
  !
  ! initialize the warning counters
  call init_warning()

  !
  ! start counting time here
  call init_timing()
  call start_clock(timing_total)

  !
  ! Output some welcome messages and compilation options
  call header()

  !
  ! Reading input file: the input parameters are stored in the module m_inputparam
  call read_inputfile_namelist()

  ! Finalize the MPI initialization
  call init_mpi_other_communicators(mpi_nproc_ortho)

  !
  ! Build all the Cartesian to Pure Gaussian transforms
  call setup_cart_to_pure_transforms()


  !
  ! Prepare relaxation with LBFGS
  if( move_nuclei == 'relax' ) then
    call lbfgs_init(lbfgs_plan,3*natom,5,diag_guess=2.0_dp)
  endif


  !
  ! Nucleus motion loop
  !
  do istep=1,nstep

    if( move_nuclei == 'relax' ) then
      write(stdout,'(/,/,1x,a,i5,/)') ' === LBFGS step ',istep
    endif

    call start_clock(timing_prescf)

    !
    ! Nucleus-nucleus repulsion contribution to the energy
    call nucleus_nucleus_energy(en_gks%nuc_nuc)

    !
    ! Build up the basis set
    !
    write(stdout,*) 'Setting up the basis set for wavefunctions'
    call init_basis_set(basis_path,basis_name,ecp_basis_name,gaussian_type,basis)

    !
    ! SCALAPACK distribution that depends on the system specific size, parameters etc.
    call init_scalapack_other(basis%nbf,eri3_nprow,eri3_npcol)

    if( print_rho_grid_ ) call dm_dump(basis)

    !
    ! Calculate overlap matrix S so to obtain "nstate" as soon as possible
    !
    call clean_allocate('Overlap matrix S',s_matrix,basis%nbf,basis%nbf)
    !
    ! Build up the overlap matrix S
    ! S only depends onto the basis set
    call setup_overlap(basis,s_matrix)

    !
    ! Calculate the square root inverse of the overlap matrix S
    ! Eliminate those eigenvalues that are too small in order to stabilize the
    ! calculation
    !
    ! A crucial parameter is defined here: nstate
    call setup_x_matrix(min_overlap,s_matrix,nstate,x_matrix)

    allocate(occupation(nstate,nspin))
    allocate(    energy(nstate,nspin))
    !
    ! Build the first occupation array
    ! as the energy are not known yet, set temperature to zero
    call set_occupation(0.0_dp,electrons,magnetization,energy,occupation)

    !
    !
    ! Precalculate the Coulomb integrals here
    !
    !
    ! ERI are to be stored in the module m_eri
    call prepare_eri(basis)

    !
    ! If an auxiliary basis is given, then set it up now
    if( has_auxil_basis ) then
      write(stdout,'(/,a)') ' Setting up the auxiliary basis set for Coulomb integrals'
      if( TRIM(capitalize(auxil_basis_name(1))) /= 'AUTO' .AND. TRIM(capitalize(auxil_basis_name(1))) /= 'PAUTO' ) then
        call init_basis_set(basis_path,auxil_basis_name,ecp_auxil_basis_name,gaussian_type,auxil_basis)
      else
        call init_auxil_basis_set_auto(auxil_basis_name,basis,gaussian_type,auto_auxil_fsam,auto_auxil_lmaxinc,auxil_basis)
      endif
    endif

    call calculation_parameters_yaml(basis%nbf,auxil_basis%nbf,nstate)

    !
    ! Attempt to evaluate the peak memory
    !
    if( memory_evaluation_ ) call evaluate_memory(basis%nbf,auxil_basis%nbf,nstate,occupation)


    if( .NOT. has_auxil_basis ) then
      !
      ! If no auxiliary basis is given,
      ! then calculate the required 4-center integrals
      call calculate_eri(print_eri_,basis,0.0_dp)
      !
      ! for Range-separated hybrids, calculate the long-range ERI
      if(calc_type%need_exchange_lr) then
        call calculate_eri(print_eri_,basis,rcut)
      endif

    else

      ! 2-center and 3-center integrals
      call calculate_eri_ri(basis,auxil_basis,0.0_dp)


      ! If Range-Separated Hybrid are requested
      ! If is_big_restart, these integrals are NOT needed, I chose code this!
      if(calc_type%need_exchange_lr ) then
        ! 2-center and 3-center integrals
        call calculate_eri_ri(basis,auxil_basis,rcut)
      endif

      call reshuffle_distribution_3center()

    endif
    ! ERI integrals have been computed and stored
    !


    !
    ! Allocate the main arrays
    ! 2D arrays
    call clean_allocate('Kinetic operator T',hamiltonian_kinetic,basis%nbf,basis%nbf)
    call clean_allocate('Nucleus operator V',hamiltonian_nucleus,basis%nbf,basis%nbf)
    call clean_allocate('Fock operator F',hamiltonian_fock,basis%nbf,basis%nbf,nspin)
    call clean_allocate('Wavefunctions C',c_matrix,basis%nbf,nstate,nspin)  ! not distributed right now

    !
    ! Try to read a RESTART file if it exists
    if( read_restart_ ) then
      call read_restart(restart_type,'RESTART',basis,occupation,c_matrix,energy,hamiltonian_fock)
    else
      restart_type = NO_RESTART
    endif
    is_restart       = ( restart_type /= NO_RESTART )
    is_big_restart   = ( restart_type == BIG_RESTART )
    is_basis_restart = ( restart_type == BASIS_RESTART )
    if( is_restart .AND. (.NOT.is_big_restart) .AND. (.NOT.is_basis_restart) ) write(stdout,*) 'Restarting from a RESTART file'
    if( is_big_restart   ) write(stdout,*) 'Restarting from a finalized RESTART file'
    if( is_basis_restart ) write(stdout,*) 'Restarting from a finalized RESTART but with a different basis set'
    ! When a BIG RESTART file is provided, assume it contains converged SCF information
    is_converged     = is_big_restart


    !
    ! Calculate the parts of the hamiltonian that does not change along
    ! with the SCF cycles
    !
    ! Kinetic energy contribution
    call setup_kinetic(basis,hamiltonian_kinetic)

    !
    ! Nucleus-electron interaction
    call setup_nucleus(basis,hamiltonian_nucleus)

    !
    ! Testing the quadrature in Fourier space
    !if( .TRUE. ) then
    !  !                        basis projectile n basis_target
    !  call setup_overlap_fourier(basis,basis,s_matrix)
    !  call setup_kinetic_fourier(basis,basis,hamiltonian_kinetic)
    !  call setup_nucleus_fourier(basis,basis,hamiltonian_nucleus)
    !endif


    if( nelement_ecp > 0 ) then
      call setup_nucleus_ecp(basis,hamiltonian_nucleus)
    endif

    !If RESTART_TDDFT file exists and is correct, skip the SCF loop and start RT-TDDFT simulation
    if( read_tddft_restart_ ) then
      call check_restart_tddft(nstate,occupation,restart_tddft_is_correct)
      ! When restart_tddft_is_correct  is TRUE, then override is_converged
      if( restart_tddft_is_correct ) is_converged = .TRUE.
    end if


    if( restart_tddft_is_correct .AND. read_tddft_restart_ ) exit

    if( is_basis_restart ) then
      !
      ! Setup the initial c_matrix by diagonalizing an approximate Hamiltonian
      call issue_warning('basis restart is not fully implemented: use with care')
      call diagonalize_hamiltonian_scalapack(hamiltonian_fock,x_matrix,energy,c_matrix)
    endif


    !
    ! For self-consistent calculations (QSMP2, QSGW, QSCOHSEX) that depend on empty states,
    ! ignore the restart file if it is not a big one
    if( calc_type%selfenergy_technique == QS ) then
      if( restart_type /= BIG_RESTART ) then
        call issue_warning('RESTART file has been ignored, since it does not contain the required empty states')
        is_restart = .FALSE.
      endif
    endif


    if( .NOT. is_restart) then

      allocate(hamiltonian_tmp(basis%nbf,basis%nbf,1))

      select case(TRIM(init_hamiltonian))
      case('GUESS')
        !
        ! Setup the initial c_matrix by diagonalizing an approximate Hamiltonian
        !
        ! Calculate a very approximate vhxc based on simple gaussians density placed on atoms
        call dft_approximate_vhxc(basis,hamiltonian_tmp(:,:,1))

        hamiltonian_tmp(:,:,1) = hamiltonian_tmp(:,:,1) + hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
      case('CORE')
        hamiltonian_tmp(:,:,1) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)

      case default
        call die('molgw: init_hamiltonian option is not valid')
      end select

      write(stdout,'(/,a)') ' Approximate hamiltonian'
      call diagonalize_hamiltonian_scalapack(hamiltonian_tmp(:,:,1:1),x_matrix,energy(:,1:1),c_matrix(:,:,1:1))

      deallocate(hamiltonian_tmp)

      ! The hamiltonian is still spin-independent:
      c_matrix(:,:,nspin) = c_matrix(:,:,1)

    endif

    call stop_clock(timing_prescf)


    !
    !
    ! Part 2 / 3 : SCF cycles
    !
    !

    !
    ! Big SCF loop is in there
    ! Only do it if the calculation is NOT a big restart
    if( .NOT. is_big_restart .AND. nscf > 0 ) then
      call scf_loop(is_restart,                                     &
                    basis,                                          &
                    x_matrix,s_matrix,                              &
                    hamiltonian_kinetic,hamiltonian_nucleus,        &
                    occupation,energy,                              &
                    hamiltonian_fock,                               &
                    c_matrix,en_gks,is_converged)
    endif

    !
    ! Big RESTART file written if converged
    !
    if( is_converged .AND. print_bigrestart_ ) then
      call write_restart(BIG_RESTART,basis,occupation,c_matrix,energy,hamiltonian_fock)
    else
      if( print_restart_ ) then
        call write_restart(SMALL_RESTART,basis,occupation,c_matrix,energy)
      endif
    endif

    !
    ! If requested, evaluate the forces
    if( move_nuclei == 'relax' ) then
      call calculate_force(basis,nstate,occupation,energy,c_matrix)
      call relax_atoms(lbfgs_plan,en_gks%total)
      call output_positions()

      if( MAXVAL(force(:,:)) < tolforce ) then
        write(stdout,'(1x,a,es16.6,a,es16.6,/)') 'Forces are     converged: ',MAXVAL(force(:,:)) , '   < ',tolforce
        exit
      else
        write(stdout,'(1x,a,es16.6,a,es16.6,/)') 'Forces are not converged: ',MAXVAL(force(:,:)) , '   > ',tolforce
        !
        ! If it is not the last step, then deallocate everything and start over
        if( istep /= nstep ) then
          call deallocate_eri()
          if( has_auxil_basis ) call destroy_eri_3center()
          if( has_auxil_basis .AND. calc_type%need_exchange_lr ) call destroy_eri_3center_lr()
          call clean_deallocate('Overlap matrix S',s_matrix)
          call clean_deallocate('Overlap X * X**H = S**-1',x_matrix)
          call clean_deallocate('Fock operator F',hamiltonian_fock)
          call clean_deallocate('Kinetic operator T',hamiltonian_kinetic)
          call clean_deallocate('Nucleus operator V',hamiltonian_nucleus)
          call clean_deallocate('Wavefunctions C',c_matrix)
          deallocate(energy,occupation)
          call destroy_basis_set(basis)
          if(has_auxil_basis) call destroy_basis_set(auxil_basis)
        endif
      endif
    endif

  enddo ! istep


  if( move_nuclei == 'relax' ) then
    call lbfgs_destroy(lbfgs_plan)
  endif

  ! This overrides the value of is_converged
  if( assume_scf_converged_ ) is_converged = .TRUE.
  if( .NOT. is_converged ) then
    call issue_warning('SCF loop is not converged. The postscf calculations (if any) will be skipped. &
                        Use keyword assume_scf_converged to override this security check')
  endif

  !
  !
  ! Part 3 / 3 : Post-processings
  !
  !
  call start_clock(timing_postscf)

  !
  ! Evaluate spin contamination
  call evaluate_s2_operator(occupation,c_matrix,s_matrix)


  if( print_multipole_ ) then
    !
    ! Evaluate the static dipole
    call static_dipole(basis,occupation,c_matrix)
    !
    ! Evaluate the static quadrupole
    call static_quadrupole(basis,occupation,c_matrix)
  endif

  if( print_wfn_ )  call plot_wfn(basis,c_matrix)
  if( print_wfn_ )  call plot_rho(basis,occupation,c_matrix)
  if( print_cube_ ) call plot_cube_wfn('GKS',basis,occupation,c_matrix)
  if( print_wfn_files_ )  call print_wfn_file('GKS',basis,occupation,c_matrix,en_gks%total,energy)
  if( print_pdos_ ) call mulliken_pdos(basis,s_matrix,c_matrix,occupation,energy)
  if( print_spatial_extension_ ) call spatial_extension(basis,c_matrix)
  if( .FALSE.     ) call plot_rho_list(nstate,basis,occupation,c_matrix)
  if( print_dens_traj_ ) call plot_rho_traj_bunch_contrib(nstate,basis,occupation,c_matrix,0,0.0_dp)
  if( print_dens_traj_points_set_ ) call plot_rho_traj_points_set_contrib(nstate,basis,occupation,c_matrix,0,0.0_dp)
  if( .FALSE. ) call write_cube_from_header('GKS',basis,occupation,c_matrix)
  if( .FALSE. ) call plot_wfn_fourier(basis,c_matrix)



  !
  ! RT-TDDFT Simulation (only if SCF cycles were converged)
  !
  if( calc_type%is_real_time .AND. is_converged ) then
    call calculate_propagation(basis,occupation,c_matrix,restart_tddft_is_correct)
  end if

  !
  ! If RSH calculations were performed, then deallocate the LR integrals which
  ! are not needed anymore
  !

  if( calc_type%need_exchange_lr ) call deallocate_eri_4center_lr()
  if( has_auxil_basis .AND. calc_type%need_exchange_lr ) call destroy_eri_3center_lr()


  !
  ! Calculate or read a correlated density matrix
  !
  if( read_fchk /= 'NO' &
     .OR. TRIM(pt_density_matrix) /= 'NO' &
     .OR. use_correlated_density_matrix_ ) then
    call get_dm_mbpt(basis,occupation,energy,c_matrix,s_matrix,hamiltonian_kinetic,hamiltonian_nucleus,hamiltonian_fock)
  endif


  call clean_deallocate('Overlap matrix S',s_matrix)
  call clean_deallocate('Kinetic operator T',hamiltonian_kinetic)
  call clean_deallocate('Nucleus operator V',hamiltonian_nucleus)
  call clean_deallocate('Overlap X * X**H = S**-1',x_matrix)



  !
  ! Prepare the diagonal of the matrix Sigma_x - Vxc
  ! for the forthcoming GW or PT corrections
  if( calc_type%selfenergy_approx > 0 .AND. calc_type%selfenergy_technique /= QS ) then
    call clean_allocate('Sigx - Vxc',exchange_m_vxc,nstate,nstate,nspin)
    call setup_exchange_m_vxc(basis,occupation,energy,c_matrix,hamiltonian_fock,exchange_m_vxc)
  endif
  call clean_deallocate('Fock operator F',hamiltonian_fock)


  !
  ! CI calculation (only if SCF cycles were converged)
  !
  if( calc_type%is_ci .AND. is_converged ) then
    if( nspin /= 1 ) call die('molgw: CI calculations need spin-restriction. Set nspin to 1')

    !
    ! Set the range of states on which to evaluate the self-energy
    call selfenergy_set_state_range(MIN(nstate,nvirtualg-1),occupation)

    if( is_virtual_fno ) then
      call calculate_virtual_fno(basis,nstate,nsemax,occupation,energy,c_matrix)
    endif
    if(has_auxil_basis) then
      call calculate_eri_3center_eigen(c_matrix)
    else
      call calculate_eri_4center_eigen_uks(c_matrix,1,MIN(nstate,nvirtualg-1))  ! TODO set the nstate_min to a more finely tuned value
    endif

    call prepare_ci(basis,MIN(nstate,nvirtualg-1),ncoreg,c_matrix)

    call full_ci_nelectrons(0,NINT(electrons),ci_spin_multiplicity-1,en_gks%nuc_nuc)

    if(calc_type%is_selfenergy) then
      if( ci_greens_function == 'BOTH' .OR. ci_greens_function == 'HOLES' ) then
        call full_ci_nelectrons( 1,NINT(electrons)-1,1,en_gks%nuc_nuc)
      endif
      if( ci_greens_function == 'BOTH' .OR. ci_greens_function == 'ELECTRONS' ) then
        call full_ci_nelectrons(-1,NINT(electrons)+1,1,en_gks%nuc_nuc)
      endif
      call full_ci_nelectrons_selfenergy()
    endif


    if(has_auxil_basis) then
      call destroy_eri_3center_eigen()
    else
      call destroy_eri_4center_eigen_uks()
    endif

    call destroy_ci()

    if( is_virtual_fno ) then
      call destroy_fno(basis,nstate,energy,c_matrix)
    endif

  endif

  !
  ! final evaluation for MP2 total energy
  !
  if( calc_type%is_mp2 ) then

    if(has_auxil_basis) then
      call mp2_energy_ri(nstate,basis,occupation,energy,c_matrix,en_gks%mp2)
    else
      call mp2_energy(nstate,basis,occupation,c_matrix,energy,en_gks%mp2)
    endif

    write(stdout,'(a,2x,f19.10)') ' MP2 Energy       (Ha):',en_gks%mp2
    write(stdout,*)
    en_gks%total = en_gks%nuc_nuc + en_gks%kinetic + en_gks%nucleus + en_gks%hartree + en_gks%exx + en_gks%mp2

    write(stdout,'(a,2x,f19.10)') ' MP2 Total Energy (Ha):',en_gks%total
    write(stdout,*)

  endif


  !
  ! final evaluation for MP3 total energy
  !
  if( calc_type%is_mp3 ) then
    if(has_auxil_basis) then
      call mp3_energy_ri(nstate,basis,occupation,energy,c_matrix,en_gks%mp3)
    else
      call die('MP3 energy without RI not implemented')
    endif
    write(stdout,'(a,2x,f19.10)') ' MP3 Energy       (Ha):',en_gks%mp3
    write(stdout,*)

    en_gks%total = en_gks%total + en_gks%mp3

    write(stdout,'(a,2x,f19.10)') ' MP3 Total Energy (Ha):',en_gks%total
    write(stdout,*)

  endif


  !
  ! Linear-response time dependent calculations work for BSE and TDDFT
  ! (only if the SCF cycles were converged)
  if( ( calc_type%is_td .OR. calc_type%is_bse ) .AND. is_converged ) then
    call init_spectral_function(nstate,occupation,0,wpol)
    call polarizability(.FALSE.,.FALSE.,basis,nstate,occupation,energy,c_matrix,erpa_tmp,egw_tmp,wpol)
    call destroy_spectral_function(wpol)
  endif

  !
  ! Self-energy calculation: PT2, GW, GWGamma, COHSEX
  ! (only if the SCF cycles were converged)
  if( calc_type%selfenergy_approx > 0 .AND. calc_type%selfenergy_technique /= QS .AND. is_converged ) then
    en_mbpt = en_gks
    call selfenergy_evaluation(basis,auxil_basis,occupation,energy,c_matrix,exchange_m_vxc,en_mbpt)
    call print_energy_yaml('mbpt energy',en_mbpt)
    call clean_deallocate('Sigx - Vxc',exchange_m_vxc)
  endif


  !
  ! Cleanly exiting the code
  !
  call clean_deallocate('Wavefunctions C',c_matrix)
  deallocate(energy,occupation)

  call deallocate_eri()
  if(has_auxil_basis) call destroy_eri_3center()

  call destroy_basis_set(basis)
  if(has_auxil_basis) call destroy_basis_set(auxil_basis)
  call destroy_atoms()

  call destroy_cart_to_pure_transforms()


  call stop_clock(timing_postscf)
  call stop_clock(timing_total)

  call this_is_the_end()


end program molgw


!=========================================================================
