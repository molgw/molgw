!=========================================================================
! This file is part of MOLGW.
!
! Copyright (C) 2010-2016  Fabien Bruneval
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
 use m_tools
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
 use m_hamiltonian
 use m_hamiltonian_sca
 use m_hamiltonian_onebody
 use m_hamiltonian_buffer
 use m_selfenergy_tools
 use m_scf_loop
 use m_tddft_propagator
 use m_tddft_variables
 use m_virtual_orbital_space
 use m_ci
 implicit none

!=====
 type(basis_set)         :: basis
 type(basis_set)         :: auxil_basis
 type(spectral_function) :: wpol
 type(lbfgs_state)       :: lbfgs_plan
 integer                 :: restart_type
 integer                 :: nstate
 integer                 :: istep
 logical                 :: is_restart,is_big_restart,is_basis_restart
 logical                 :: restart_tddft_is_correct = .TRUE.
 real(dp),allocatable    :: hamiltonian_tmp(:,:,:)
 real(dp),allocatable    :: hamiltonian_kinetic(:,:)
 real(dp),allocatable    :: hamiltonian_nucleus(:,:)
 real(dp),allocatable    :: hamiltonian_fock(:,:,:)
 real(dp),allocatable    :: s_matrix(:,:)
 real(dp),allocatable    :: s_matrix_sqrt_inv(:,:)
 real(dp),allocatable    :: c_matrix(:,:,:),c_matrix_tmp(:,:,:)
 real(dp),allocatable    :: energy(:,:)
 real(dp),allocatable    :: occupation(:,:)
 real(dp),allocatable    :: exchange_m_vxc_diag(:,:)
 real(dp),allocatable    :: exchange_m_vxc(:,:,:)
 integer                 :: m_ham,n_ham                  ! distribute a  basis%nbf x basis%nbf   matrix
 integer                 :: m_c,n_c                      ! distribute a  basis%nbf x nstate      matrix
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
 ! Output some welcome message and compilation options
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
   call nucleus_nucleus_energy(en%nuc_nuc)

   !
   ! Build up the basis set
   !
   write(stdout,*) 'Setting up the basis set for wavefunctions'
   call init_basis_set(basis_path,basis_name,ecp_basis_name,gaussian_type,basis)

   !
   ! SCALAPACK distribution that depends on the system specific size, parameters etc.
   call init_scalapack_other(basis%nbf,eri3_nprow,eri3_npcol,scalapack_nprow,scalapack_npcol,m_ham,n_ham)
   if( m_ham /= basis%nbf .OR. n_ham /= basis%nbf ) then
     call issue_warning('SCALAPACK is used to distribute the SCF hamiltonian')
   endif


   !
   !
   ! Precalculate the Coulomb integrals here
   !
   !
   ! ERI are stored "privately" in the module m_eri
   call prepare_eri(basis)

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
     !
     ! If an auxiliary basis is given,
     ! then set it up now and calculate the required ERI: 2- and 3-center integrals
     write(stdout,'(/,a)') ' Setting up the auxiliary basis set for Coulomb integrals'
     call init_basis_set(basis_path,auxil_basis_name,ecp_auxil_basis_name,gaussian_type,auxil_basis)

     ! 2-center integrals
     call calculate_eri_2center_scalapack(auxil_basis,0.0_dp)
     ! 3-center integrals
     call calculate_eri_3center_scalapack(basis,auxil_basis,0.0_dp)


     ! If Range-Separated Hybrid are requested
     ! If is_big_restart, these integrals are NOT needed, I chose code this!
     if(calc_type%need_exchange_lr ) then
       ! 2-center integrals
       call calculate_eri_2center_scalapack(auxil_basis,rcut)
       ! 3-center integrals
       call calculate_eri_3center_scalapack(basis,auxil_basis,rcut)
     endif

   endif
   ! ERI integrals have been computed and stored
   !



   !
   ! Allocate the main arrays
   ! 2D arrays
   call clean_allocate('Overlap matrix S',s_matrix,m_ham,n_ham)
   call clean_allocate('Kinetic operator T',hamiltonian_kinetic,m_ham,n_ham)
   call clean_allocate('Nucleus operator V',hamiltonian_nucleus,m_ham,n_ham)
   call clean_allocate('Fock operator F',hamiltonian_fock,basis%nbf,basis%nbf,nspin) ! Never distributed


   ! Allocate the only complete array buffer basis%nbf x basis%nbf in case of SCALAPACK
   if( parallel_ham .AND. parallel_buffer ) call allocate_parallel_buffer(basis%nbf)

   !
   ! Build up the overlap matrix S
   ! S only depends onto the basis set
   if( parallel_ham ) then
     call setup_overlap_buffer_sca(basis,s_matrix)
   else
     call setup_overlap(basis,s_matrix)
   endif

   !
   ! Calculate the square root inverse of the overlap matrix S
   ! Eliminate those eigenvalues that are too small in order to stabilize the
   ! calculation
   !
   ! Crucial parameters are defined here: nstate
   !                                      m_c and n_c
   !
   if( parallel_ham ) then
     call setup_sqrt_overlap_sca(min_overlap,desc_ham,s_matrix,desc_c,nstate,s_matrix_sqrt_inv)
     m_c    = SIZE( s_matrix_sqrt_inv , DIM=1 )
     n_c    = SIZE( s_matrix_sqrt_inv , DIM=2 )

   else
     call setup_sqrt_overlap(min_overlap,s_matrix,nstate,s_matrix_sqrt_inv)

     m_c = basis%nbf
     n_c = nstate

   endif

   if( m_c /= basis%nbf .OR. n_c /= nstate ) then
     call issue_warning('SCALAPACK is used to distribute the wavefunction coefficients')
   endif


   ! Allocate the nstate arrays: c_matrix, occupation, energy
   ! 2D arrays
   call clean_allocate('Wavefunctions C',c_matrix,basis%nbf,nstate,nspin)  ! not distributed right now
   ! 1D arrays
   allocate(occupation(nstate,nspin))
   allocate(    energy(nstate,nspin))


   !
   ! Build the first occupation array
   ! as the energy are not known yet, set temperature to zero
   call set_occupation(0.0_dp,electrons,magnetization,energy,occupation)

   !
   ! Try to read a RESTART file if it exists
   if( read_restart_ ) then
     call read_restart(restart_type,basis,nstate,occupation,c_matrix,energy,hamiltonian_fock,'RESTART')
   else
     restart_type = NO_RESTART
   endif
   is_restart       = ( restart_type /= NO_RESTART )
   is_big_restart   = ( restart_type == BIG_RESTART )
   is_basis_restart = ( restart_type == BASIS_RESTART )
   if( is_restart .AND. (.NOT.is_big_restart) .AND. (.NOT.is_basis_restart) ) write(stdout,*) 'Restarting from a RESTART file'
   if( is_big_restart   ) write(stdout,*) 'Restarting from a finalized RESTART file'
   if( is_basis_restart ) write(stdout,*) 'Restarting from a finalized RESTART but with a different basis set'

   if( parallel_ham .AND. .NOT. is_big_restart ) then
     call clean_allocate('Wavefunctions C',c_matrix_tmp,m_c,n_c,nspin)
     call create_distributed_copy(c_matrix,desc_c,c_matrix_tmp)
     call clean_deallocate('Wavefunctions C',c_matrix)
     call move_alloc(c_matrix_tmp,c_matrix)
   endif

   !
   ! Calculate the parts of the hamiltonian that does not change along
   ! with the SCF cycles
   !
   ! Kinetic energy contribution
   if( parallel_ham ) then
     call setup_kinetic_buffer_sca(basis,hamiltonian_kinetic)
   else
     call setup_kinetic(basis,hamiltonian_kinetic)
   endif

   !
   ! Nucleus-electron interaction
   if( parallel_ham ) then
     call setup_nucleus_buffer_sca(basis,hamiltonian_nucleus)
     if( nelement_ecp > 0 ) then
       call die('ECP not implemented with SCALAPACK yet. Set scalapack_nprow = scalapack_npcol to 1.')
     endif
   else
     call setup_nucleus(basis,hamiltonian_nucleus)

     if( nelement_ecp > 0 ) then
       call setup_nucleus_ecp(basis,hamiltonian_nucleus)
     endif
   endif

   !If RESTART_TDDFT file exists and is correct, skip the SCF loop and start RT-TDDFT simulation
   if( read_tddft_restart_ ) then
     call check_restart_tddft(nstate,occupation,restart_tddft_is_correct)
   end if


   if( restart_tddft_is_correct .AND. read_tddft_restart_ ) exit

   if( is_basis_restart ) then
     !
     ! Setup the initial c_matrix by diagonalizing an approximate Hamiltonian
     if( parallel_ham ) call die('basis_restart not implemented with distributed hamiltonian')
     call issue_warning('basis restart is not fully implemented: use with care')
     call diagonalize_hamiltonian_scalapack(hamiltonian_fock,s_matrix_sqrt_inv,energy,c_matrix)
   endif


   !
   ! For self-consistent calculations (QSMP2, QSGW, QSCOHSEX) that depend on empty states,
   ! ignore the restart file if it is not a big one
   if( calc_type%selfenergy_technique == QS ) then
     if( restart_type /= EMPTY_STATES_RESTART .AND. restart_type /= BIG_RESTART ) then
       call issue_warning('RESTART file has been ignored, since it does not contain the required empty states')
       is_restart = .FALSE.
     endif
   endif


   if( .NOT. is_restart) then

     allocate(hamiltonian_tmp(m_ham,n_ham,1))

     select case(TRIM(init_hamiltonian))
     case('GUESS')
       !
       ! Setup the initial c_matrix by diagonalizing an approximate Hamiltonian
       !
       ! Calculate a very approximate vhxc based on simple gaussians placed on atoms
       if( parallel_ham ) then
         if( parallel_buffer ) then
           call dft_approximate_vhxc_buffer_sca(basis,m_ham,n_ham,hamiltonian_tmp(:,:,1))
         else
           call dft_approximate_vhxc_sca(basis,m_ham,n_ham,hamiltonian_tmp(:,:,1))
         endif
       else
         call dft_approximate_vhxc(basis,hamiltonian_tmp(:,:,1))
       endif

       hamiltonian_tmp(:,:,1) = hamiltonian_tmp(:,:,1) + hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
     case('CORE')
       hamiltonian_tmp(:,:,1) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
     case default
       call die('molgw: init_hamiltonian option is not valid')
     end select


     write(stdout,'(/,a)') ' Approximate hamiltonian'
     if( parallel_ham ) then
       call diagonalize_hamiltonian_sca(desc_ham,hamiltonian_tmp(:,:,1:1),desc_c,s_matrix_sqrt_inv,energy(:,1:1),c_matrix(:,:,1:1))
     else
       call diagonalize_hamiltonian_scalapack(hamiltonian_tmp(:,:,1:1),s_matrix_sqrt_inv,energy(:,1:1),c_matrix(:,:,1:1))
     endif

     deallocate(hamiltonian_tmp)

     ! The hamiltonian is still spin-independent:
     c_matrix(:,:,nspin) = c_matrix(:,:,1)

   endif


   ! Deallocate the array index_pair when it is not needed
   ! This array can be pretty large indeed.
   ! index_pair is only need for the development version of Luttinger-Ward
   ! or for the SCALAPACK with no buffer part of the code (which is not functional anyway)
   ! or for the 4-center integrals code
   if( .NOT. (calc_type%selfenergy_approx == LW                         &
              .OR. calc_type%selfenergy_approx == LW2                   &
              .OR. calc_type%selfenergy_approx == GSIGMA )              &
       .AND. .NOT. ( parallel_ham .AND. .NOT. parallel_buffer ) &
       .AND. has_auxil_basis ) then
     call deallocate_index_pair()
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
                   nstate,m_ham,n_ham,m_c,n_c,                     &
                   s_matrix_sqrt_inv,s_matrix,                     &
                   hamiltonian_kinetic,hamiltonian_nucleus,        &
                   occupation,energy,                              &
                   hamiltonian_fock,                               &
                   c_matrix)
   else
     if( parallel_ham .AND. parallel_buffer ) call destroy_parallel_buffer()
   endif


   !
   ! If requested, evaluate the forces
   if( move_nuclei == 'relax' ) then
     call calculate_force(basis,nstate,occupation,energy,c_matrix)
     call relax_atoms(lbfgs_plan,en%tot)
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
         call clean_deallocate('Overlap sqrt S^{-1/2}',s_matrix_sqrt_inv)
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




 !
 !
 ! Part 3 / 3 : Post-processings
 !
 !
 call start_clock(timing_postscf)


 ! temporary section for the charge calculation
!  call init_dft_grid(basis,grid_level,dft_xc_needs_gradient,.TRUE.,64)
!  call calc_normalization_r(64,basis,occupation,c_matrix)
!  call destroy_dft_grid()

 if( print_multipole_ ) then
   !
   ! Evaluate the static dipole
   call static_dipole(nstate,basis,occupation,c_matrix)
   !
   ! Evaluate the static quadrupole
   call static_quadrupole(nstate,basis,occupation,c_matrix)
 endif

 if( print_wfn_ )  call plot_wfn(nstate,basis,c_matrix)
 if( print_wfn_ )  call plot_rho(nstate,basis,occupation,c_matrix)
 if( print_cube_ ) call plot_cube_wfn('GKS',nstate,basis,occupation,c_matrix)
 if( print_pdos_ ) call mulliken_pdos(nstate,basis,s_matrix,c_matrix,occupation,energy)
 if( .FALSE.     ) call plot_rho_list(nstate,basis,occupation,c_matrix)
 if( print_dens_traj_ ) call plot_rho_traj_bunch_contrib(nstate,basis,occupation,c_matrix,0,0.0_dp)
 if( .FALSE. ) call read_cube_wfn(nstate,basis,occupation,c_matrix)

 call clean_deallocate('Overlap matrix S',s_matrix)
 call clean_deallocate('Kinetic operator T',hamiltonian_kinetic)
 call clean_deallocate('Nucleus operator V',hamiltonian_nucleus)
 call clean_deallocate('Overlap sqrt S^{-1/2}',s_matrix_sqrt_inv)


 !
 !
 ! Post-processing start here
 !
 !

 !
 ! RT-TDDFT Simulation
 if(calc_type%is_real_time) then
   call calculate_propagation(basis,occupation,c_matrix)
 end if

 !
 ! If RSH calculations were performed, then deallocate the LR integrals which
 ! are not needed anymore
 !

 if( calc_type%need_exchange_lr ) call deallocate_eri_4center_lr()
 if( has_auxil_basis .AND. calc_type%need_exchange_lr ) call destroy_eri_3center_lr()

 ! Performs a distribution strategy change here for the 3-center integrals
 ! ( nprow_3center x npcol_3center ) => ( nprow_auxil x 1 )
 if( has_auxil_basis ) then
   if( calc_type%selfenergy_approx > 0 .OR. calc_type%is_ci .OR. calc_type%is_td &
       .OR. calc_type%is_mp2 .OR. calc_type%is_mp3 .OR. calc_type%is_bse ) then
     call reshuffle_distribution_3center()
   endif
 endif

 !
 !
 ! Post-processing start here
 !
 !

 ! Performs a distribution strategy change here for the 3-center integrals
 ! ( nprow_3center x npcol_3center ) => ( nprow_auxil x 1 )
 if( has_auxil_basis ) then
   if( calc_type%selfenergy_approx > 0 .OR. calc_type%is_ci .OR. calc_type%is_td &
       .OR. calc_type%is_mp2 .OR. calc_type%is_mp3 .OR. calc_type%is_bse ) then
     call reshuffle_distribution_3center()
   endif
 endif

 !
 ! Prepare the diagonal of the matrix Sigma_x - Vxc
 ! for the forthcoming GW or PT corrections
 if( calc_type%selfenergy_approx > 0 .AND. calc_type%selfenergy_technique /= QS ) then

   allocate(exchange_m_vxc_diag(nstate,nspin))

   if( calc_type%selfenergy_static ) then
     !
     ! Calculate the static part of the self-energy at the first order and store it in exchange_m_vxc_diag
     !
     allocate(exchange_m_vxc(nstate,nstate,nspin))

     call setup_exchange_m_vxc(basis,occupation,energy,c_matrix,hamiltonian_fock,exchange_m_vxc_diag,exchange_m_vxc)

     call selfenergy_set_state_range(MIN(nstate,nvirtualg-1),occupation)
     call pt1_selfenergy(nstate,basis,occupation,energy,c_matrix,exchange_m_vxc,exchange_m_vxc_diag)

     deallocate(exchange_m_vxc)

   else
     call setup_exchange_m_vxc(basis,occupation,energy,c_matrix,hamiltonian_fock,exchange_m_vxc_diag)
   endif

 endif
 call clean_deallocate('Fock operator F',hamiltonian_fock)


 !
 ! CI calculation
 !
 if(calc_type%is_ci) then
   if(nspin/=1) call die('molgw: CI calculations need spin-restriction. Set nspin to 1')

   !
   ! Set the range of states on which to evaluate the self-energy
   call selfenergy_set_state_range(MIN(nstate,nvirtualg-1),occupation)

   if( is_virtual_fno ) then
     call virtual_fno(basis,nstate,nsemax,occupation,energy,c_matrix)
   endif
   if(has_auxil_basis) then
     call calculate_eri_3center_eigen(c_matrix)
   else
     call calculate_eri_4center_eigen_uks(c_matrix,1,MIN(nstate,nvirtualg-1))  ! TODO set the nstate_min to a more finely tuned value
   endif

   call prepare_ci(basis,MIN(nstate,nvirtualg-1),ncoreg,c_matrix)

   call full_ci_nelectrons(0,NINT(electrons),ci_spin_multiplicity-1,en%nuc_nuc)

   if(calc_type%is_selfenergy) then
     if( ci_greens_function == 'BOTH' .OR. ci_greens_function == 'HOLES' ) then
       call full_ci_nelectrons( 1,NINT(electrons)-1,1,en%nuc_nuc)
     endif
     if( ci_greens_function == 'BOTH' .OR. ci_greens_function == 'ELECTRONS' ) then
       call full_ci_nelectrons(-1,NINT(electrons)+1,1,en%nuc_nuc)
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
     call mp2_energy_ri(nstate,basis,occupation,energy,c_matrix,en%mp2)
   else
     call mp2_energy(nstate,basis,occupation,c_matrix,energy,en%mp2)
   endif

   write(stdout,'(a,2x,f19.10)') ' MP2 Energy       (Ha):',en%mp2
   write(stdout,*)
   en%tot = en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%mp2

   write(stdout,'(a,2x,f19.10)') ' MP2 Total Energy (Ha):',en%tot
   write(stdout,'(a,2x,f19.10)') ' SE+MP2  Total En (Ha):',en%tot+en%se
   write(stdout,*)

 endif


 !
 ! final evaluation for MP3 total energy
 !
 if( calc_type%is_mp3 ) then
   if(has_auxil_basis) then
     call mp3_energy_ri(nstate,basis,occupation,energy,c_matrix,en%mp3)
   else
     call die('MP3 energy without RI not implemented')
   endif
   write(stdout,'(a,2x,f19.10)') ' MP3 Energy       (Ha):',en%mp3
   write(stdout,*)

   en%tot = en%tot + en%mp3

   write(stdout,'(a,2x,f19.10)') ' MP3 Total Energy (Ha):',en%tot
   write(stdout,'(a,2x,f19.10)') ' SE+MP3  Total En (Ha):',en%tot+en%se
   write(stdout,*)

 endif


 !
 ! Time Dependent calculations
 ! works for DFT, HF, and hybrid
 !
 if(calc_type%is_td .OR. calc_type%is_bse) then
   call init_spectral_function(nstate,occupation,0,wpol)
   call polarizability(.FALSE.,.FALSE.,basis,nstate,occupation,energy,c_matrix,en%rpa,wpol)
   call destroy_spectral_function(wpol)
 endif

 !
 ! Self-energy calculation: PT2, GW, GWGamma, COHSEX
 !
 if( calc_type%selfenergy_approx > 0 .AND. calc_type%selfenergy_technique /= QS ) then
   call selfenergy_evaluation(basis,auxil_basis,nstate,occupation,energy,c_matrix,exchange_m_vxc_diag)
   deallocate(exchange_m_vxc_diag)
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

 call total_memory_statement()

 call stop_clock(timing_postscf)
 call stop_clock(timing_total)
 call output_timing()

 call output_all_warnings()

 write(stdout,'(/,1x,a,/)') 'This is the end'

 call finish_mpi()


end program molgw


!=========================================================================
