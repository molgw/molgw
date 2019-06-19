!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the main SCF loop for Hartree-Fock or generalized Kohn-Sham
!
!=========================================================================
module m_scf_loop
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_linear_algebra
 use m_atoms
 use m_mpi
 use m_scalapack
 use m_basis_set
 use m_inputparam
 use m_scf
 use m_eri_ao_mo
 use m_dft_grid
 use m_spectral_function
 use m_hamiltonian_tools
 use m_hamiltonian_twobodies
 use m_hamiltonian_wrapper
 use m_hamiltonian_cmplx
 use m_selfenergy_tools
 use m_restart



contains


!=========================================================================
subroutine scf_loop(is_restart,&
                    basis,&
                    x_matrix,s_matrix,&
                    hamiltonian_kinetic,hamiltonian_nucleus,&
                    occupation, &
                    energy, &
                    hamiltonian_fock,&
                    c_matrix,en_gks,is_converged)
 implicit none

!=====
 logical,intent(in)                 :: is_restart
 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: x_matrix(:,:)
 real(dp),intent(in)                :: s_matrix(:,:)
 real(dp),intent(in)                :: hamiltonian_kinetic(:,:)
 real(dp),intent(in)                :: hamiltonian_nucleus(:,:)
 real(dp),intent(inout)             :: occupation(:,:)
 real(dp),intent(out)               :: energy(:,:)
 real(dp),allocatable,intent(inout) :: hamiltonian_fock(:,:,:)
 real(dp),allocatable,intent(inout) :: c_matrix(:,:,:)
 type(energy_contributions),intent(inout) :: en_gks
 logical,intent(out)                :: is_converged
!=====
 type(spectral_function) :: wpol
 integer                 :: nstate
 logical                 :: stopfile_found
 integer                 :: file_density_matrix
 integer                 :: ispin,iscf
 real(dp),allocatable    :: p_matrix(:,:,:)
 real(dp),allocatable    :: hamiltonian(:,:,:)
 real(dp),allocatable    :: hamiltonian_hartree(:,:)
 real(dp),allocatable    :: hamiltonian_exx(:,:,:)
 real(dp),allocatable    :: hamiltonian_xc(:,:,:)
 real(dp),allocatable    :: matrix_tmp(:,:,:)
!=====


 call start_clock(timing_scf)

 nstate = SIZE(x_matrix,DIM=2)

 !
 ! Initialize the SCF mixing procedure
 call init_scf(basis%nbf,nstate)

 !
 ! Allocate the main arrays
 call clean_allocate('Total Hamiltonian H',hamiltonian,basis%nbf,basis%nbf,nspin)
 call clean_allocate('Hartree potential Vh',hamiltonian_hartree,basis%nbf,basis%nbf)
 call clean_allocate('Exchange operator Sigx',hamiltonian_exx,basis%nbf,basis%nbf,nspin)
 call clean_allocate('XC operator Vxc',hamiltonian_xc,basis%nbf,basis%nbf,nspin)
 call clean_allocate('Density matrix P',p_matrix,basis%nbf,basis%nbf,nspin)
 hamiltonian_exx(:,:,:) = 0.0_dp


 !
 ! Setup the grids for the quadrature of DFT potential/energy
 if( calc_type%is_dft ) then
   call init_dft_grid(basis,grid_level,dft_xc(1)%needs_gradient,.TRUE.,BATCH_SIZE)
 endif

 !
 ! Setup the density matrix: p_matrix
 call setup_density_matrix(c_matrix,occupation,p_matrix)


 !
 ! Start the big scf loop
 !
 do iscf=1,nscf
   write(stdout,'(/,a)') '-------------------------------------------'
   write(stdout,'(a,1x,i4,/)') ' *** SCF cycle No:',iscf


   en_gks%kin  = SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
   en_gks%nuc  = SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix(:,:,:),DIM=3) )

   !
   ! Setup kinetic and nucleus contributions (that are independent of the
   ! density matrix and therefore of spin channel)
   !
   do ispin=1,nspin
     hamiltonian(:,:,ispin) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
   enddo

   !
   ! Hartree contribution to the Hamiltonian
   !
   call calculate_hartree(basis,p_matrix,hamiltonian_hartree,eh=en_gks%hart)

   ! calc_type%is_core is an inefficient way to get the Kinetic+Nucleus Hamiltonian
   if( calc_type%is_core ) then
     hamiltonian_hartree(:,:) = 0.0_dp
     en_gks%hart = 0.0_dp
   endif
   do ispin=1,nspin
     hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) + hamiltonian_hartree(:,:)
   enddo


   !
   !  XC part of the Hamiltonian
   !
   hamiltonian_xc(:,:,:) = 0.0_dp
   en_gks%exx_hyb = 0.0_dp

   !
   ! DFT XC potential is added here
   ! hamiltonian_xc is used as a temporary matrix
   if( calc_type%is_dft ) then
     call dft_exc_vxc_batch(BATCH_SIZE,basis,occupation,c_matrix,hamiltonian_xc,en_gks%xc)
   endif

   !
   ! LR Exchange contribution to the Hamiltonian
   ! Use hamiltonian_exx as a temporary matrix (no need to save it for later use)
   if(calc_type%need_exchange_lr) then

     call calculate_exchange_lr(basis,p_matrix,hamiltonian_exx,ex=en_gks%exx_hyb,occupation=occupation,c_matrix=c_matrix)
     ! Rescale with alpha_hybrid_lr for range-separated hybrid functionals
     en_gks%exx_hyb = en_gks%exx_hyb * alpha_hybrid_lr
     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:) * alpha_hybrid_lr

   endif

   !
   ! Exchange contribution to the Hamiltonian
   if( calc_type%need_exchange ) then

     call calculate_exchange(basis,p_matrix,hamiltonian_exx,ex=en_gks%exx,occupation=occupation,c_matrix=c_matrix)

     ! Rescale with alpha_hybrid for hybrid functionals
     en_gks%exx_hyb = en_gks%exx_hyb + alpha_hybrid * en_gks%exx
     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:) * alpha_hybrid
   endif


   !
   ! QSGW or COHSEX self energy
   !
   if( ( calc_type%selfenergy_approx == GW .OR. calc_type%selfenergy_approx == COHSEX ) &
        .AND. calc_type%selfenergy_technique == QS  &
        .AND. ( iscf > 5 .OR. is_restart ) ) then


     call init_spectral_function(nstate,occupation,0,wpol)
     call polarizability(.TRUE.,.TRUE.,basis,nstate,occupation,energy,c_matrix,en_gks%rpa,wpol)

     if( ABS(en_gks%rpa) > 1.e-6_dp) then
       en_gks%tot = en_gks%tot + en_gks%rpa
       write(stdout,'(/,a,f19.10)') ' RPA Total energy (Ha): ',en_gks%tot
     endif

     !
     ! Set the range of states on which to evaluate the self-energy
     call selfenergy_set_state_range(nstate,occupation)

     allocate(matrix_tmp(basis%nbf,basis%nbf,nspin))
     call gw_selfenergy_qs(nstate,basis,occupation,energy,c_matrix,s_matrix,wpol,matrix_tmp)

     call dump_out_matrix(.FALSE.,'=== Self-energy ===',basis%nbf,nspin,matrix_tmp)
     call destroy_spectral_function(wpol)

     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + matrix_tmp(:,:,:)
     deallocate(matrix_tmp)

   endif

   !
   ! QSPT2
   !
   if( calc_type%selfenergy_approx == PT2 .AND. calc_type%selfenergy_technique == QS .AND. ( iscf > 5 .OR. is_restart ) ) then

     !
     ! Set the range of states on which to evaluate the self-energy
     call selfenergy_set_state_range(nstate,occupation)

     allocate(matrix_tmp(basis%nbf,basis%nbf,nspin))
     call pt2_selfenergy_qs(nstate,basis,occupation,energy,c_matrix,s_matrix,matrix_tmp,en_gks%mp2)

     write(stdout,'(a,2x,f19.10)') ' MP2 Energy       (Ha):',en_gks%mp2
     write(stdout,*)
     en_gks%tot = en_gks%tot + en_gks%mp2
     write(stdout,'(a,2x,f19.10)') ' MP2 Total Energy (Ha):',en_gks%tot

     call dump_out_matrix(.FALSE.,'=== Self-energy ===',basis%nbf,nspin,matrix_tmp)

     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + matrix_tmp(:,:,:)
     deallocate(matrix_tmp)

   endif


   !
   ! Add the XC part of the hamiltonian to the total hamiltonian
   hamiltonian(:,:,:) = hamiltonian(:,:,:) + hamiltonian_xc(:,:,:)


   ! All the components of the energy have been calculated at this stage
   ! Sum up to get the total energy
   en_gks%tot = en_gks%nuc_nuc + en_gks%kin + en_gks%nuc + en_gks%hart + en_gks%exx_hyb + en_gks%xc

   ! Make sure all the MPI threads have the exact same Hamiltonian
   ! It helps stabilizing the SCF cycles in parallel
   call xsum_world(hamiltonian)
   hamiltonian(:,:,:) = hamiltonian(:,:,:) / REAL(nproc_world,dp)

   !
   ! If requested, the level shifting procedure is triggered:
   ! All the unoccupied states are penalized with an energy =  level_shifting_energy
   if( level_shifting_energy > 1.0e-6_dp ) then
     call level_shifting_up(s_matrix,c_matrix,occupation,level_shifting_energy,hamiltonian)
   endif

   ! DIIS or simple mixing on the hamiltonian
   call hamiltonian_prediction(s_matrix,x_matrix,p_matrix,hamiltonian,en_gks%tot)


   !
   ! Diagonalize the Hamiltonian H
   ! Generalized eigenvalue problem with overlap matrix S
   ! H \varphi = E S \varphi
   ! save the old eigenvalues
   ! This subroutine works with or without scalapack
   call diagonalize_hamiltonian_scalapack(hamiltonian,x_matrix,energy,c_matrix)

   !
   ! When level_shifting is used, the unoccupied state energies have to be brought
   ! back to their original value,
   ! So that the "physical" energies are written down
   if( level_shifting_energy > 1.0e-6_dp ) then
     call level_shifting_down(s_matrix,c_matrix,occupation,level_shifting_energy,energy,hamiltonian)
   endif

   call dump_out_energy('=== Energies ===',nstate,nspin,occupation,energy)

   call output_new_homolumo('gKS',nstate,occupation,energy,1,nstate)


   !
   ! Output the total energy and its components
   write(stdout,*)
   write(stdout,'(a25,1x,f19.10)') 'Nucleus-Nucleus (Ha):',en_gks%nuc_nuc
   write(stdout,'(a25,1x,f19.10)') 'Kinetic Energy  (Ha):',en_gks%kin
   write(stdout,'(a25,1x,f19.10)') 'Nucleus Energy  (Ha):',en_gks%nuc
   write(stdout,'(a25,1x,f19.10)') 'Hartree Energy  (Ha):',en_gks%hart
   if(calc_type%need_exchange) then
     write(stdout,'(a25,1x,f19.10)') 'Exchange Energy (Ha):',en_gks%exx_hyb
   endif
   if( calc_type%is_dft ) then
     write(stdout,'(a25,1x,f19.10)') 'XC Energy       (Ha):',en_gks%xc
   endif
   write(stdout,'(/,a25,1x,f19.10,/)') 'Total Energy    (Ha):',en_gks%tot


   ! If fractional occupancies are allowed, then recalculate the occupations
   if( temperature > 1.0e-8_dp ) then
     call set_occupation(temperature,electrons,magnetization,energy,occupation)
   endif

   !
   ! Setup the new density matrix: p_matrix
   ! Save the old one for the convergence criterium
   call setup_density_matrix(c_matrix,occupation,p_matrix)


   !
   ! p_matrix preconditioning to damp out charge oscillations
   !
   if( kerker_k0 > 1.0e-6_dp .OR. density_matrix_damping > 1.0e-6_dp ) &
      call density_matrix_preconditioning(hamiltonian_kinetic,s_matrix,p_matrix)


   is_converged = check_converged(p_matrix)
   inquire(file='STOP',exist=stopfile_found)

   if( is_converged .OR. stopfile_found ) exit

   !
   ! Write down a "small" RESTART file at each step
   if( print_restart_ ) then
     call write_restart(SMALL_RESTART,basis,occupation,c_matrix,energy)
   endif


 !
 ! end of the big SCF loop
 enddo


 write(stdout,'(/,1x,a)') '=================================================='
 write(stdout,'(1x,a)') 'The SCF loop ends here'
 write(stdout,'(1x,a)') '=================================================='

 !
 ! Cleanly deallocate the integral grid information
 ! and the scf mixing information
 !
 call destroy_scf()
 if( calc_type%is_dft ) call destroy_dft_grid()


 if( print_hartree_ ) then
   call print_hartee_expectation(basis,p_matrix,c_matrix,occupation,hamiltonian_hartree,hamiltonian_exx)
 endif

 if( print_density_matrix_ .AND. is_iomaster ) then
   write(stdout,'(1x,a)') 'Write DENSITY_MATRIX_GKS file'
   open(newunit=file_density_matrix,file='DENSITY_MATRIX_GKS',form='unformatted',action='write')
   do ispin=1,nspin
     write(file_density_matrix) p_matrix(:,:,ispin)
   enddo
   close(file_density_matrix)
 endif


 !
 ! Form the final Fock matrix and store it if needed
 !
 if( is_converged  &
   .AND. ( print_bigrestart_  &
          .OR. TRIM(pt_density_matrix) /= 'NO'   &
          .OR. calc_type%selfenergy_approx > 0  )  ) then
   call get_fock_operator(basis,p_matrix,c_matrix,occupation,en_gks, &
                             hamiltonian,hamiltonian_xc,hamiltonian_exx,hamiltonian_fock)
 endif

 !
 ! Cleanly deallocate the arrays
 !
 call clean_deallocate('Density matrix P',p_matrix)
 call clean_deallocate('Total Hamiltonian H',hamiltonian)
 call clean_deallocate('Hartree potential Vh',hamiltonian_hartree)
 call clean_deallocate('Exchange operator Sigx',hamiltonian_exx)
 call clean_deallocate('XC operator Vxc',hamiltonian_xc)


 write(stdout,'(/,/,a25,1x,f19.10,/)') 'SCF Total Energy (Ha):',en_gks%tot

 if( ABS(en_gks%exx) > 1.0e-6) then
   write(stdout,'(a25,1x,f19.10)')       '      EXX Energy (Ha):',en_gks%exx
   write(stdout,'(a25,1x,f19.10)')       'Total EXX Energy (Ha):',en_gks%nuc_nuc + en_gks%kin + en_gks%nuc + en_gks%hart + en_gks%exx
 endif


 call stop_clock(timing_scf)

end subroutine scf_loop


!=========================================================================
subroutine get_fock_operator(basis,p_matrix,c_matrix,occupation,en, &
                             hamiltonian,hamiltonian_xc,hamiltonian_exx,hamiltonian_fock)
 implicit none

 type(basis_set),intent(in)               :: basis
 type(energy_contributions),intent(inout) :: en
 real(dp),intent(in)        :: p_matrix(:,:,:),c_matrix(:,:,:)
 real(dp),intent(in)        :: occupation(:,:)
 real(dp),intent(in)        :: hamiltonian(:,:,:)
 real(dp),intent(in)        :: hamiltonian_xc(:,:,:)
 real(dp),intent(inout)     :: hamiltonian_exx(:,:,:)
 real(dp),intent(out)       :: hamiltonian_fock(:,:,:)
!=====
!=====

 !
 ! Get the exchange operator if not already calculated
 !
 if( ALL( ABS(hamiltonian_exx(:,:,:)) < 1.0e-6_dp ) ) then
   call calculate_exchange(basis,p_matrix,hamiltonian_exx,ex=en%exx,occupation=occupation,c_matrix=c_matrix)
 endif

 hamiltonian_fock(:,:,:) = hamiltonian(:,:,:) - hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:)

end subroutine get_fock_operator


!=========================================================================
! Print out the Hartree and exchange diagonal expectation values if requested
subroutine print_hartee_expectation(basis,p_matrix,c_matrix,occupation,hamiltonian_hartree,hamiltonian_exx)
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: p_matrix(:,:,:),c_matrix(:,:,:)
 real(dp),intent(in)        :: occupation(:,:)
 real(dp),intent(in)        :: hamiltonian_hartree(:,:)
 real(dp),intent(inout)     :: hamiltonian_exx(:,:,:)
!=====
 integer                 :: restart_type
 integer                 :: nstate,istate,ispin
 real(dp),allocatable    :: c_matrix_restart(:,:,:)
 real(dp),allocatable    :: hartree_ii(:,:),exchange_ii(:,:)
 real(dp),allocatable    :: energy_restart(:,:),occupation_restart(:,:)
!=====

 nstate = SIZE(c_matrix,DIM=2)

 !
 ! Get the exchange operator if not already calculated
 !
 if( ALL( ABS(hamiltonian_exx(:,:,:)) < 1.0e-6_dp ) ) then
   call calculate_exchange(basis,p_matrix,hamiltonian_exx,occupation=occupation,c_matrix=c_matrix)
 endif

 call clean_allocate('RESTART: C',c_matrix_restart,basis%nbf,nstate,nspin)
 allocate(energy_restart(nstate,nspin))
 allocate(occupation_restart(nstate,nspin))

 call read_restart(restart_type,'RESTART_TEST',basis,occupation_restart,c_matrix_restart,energy_restart)

 if( restart_type == NO_RESTART ) then
   c_matrix_restart(:,:,:) = c_matrix(:,:,:)
 else
   write(stdout,'(1x,a,a)') 'RESTART file read: ','RESTART_TEST'
 endif

 do ispin=1,nspin
   do istate=1,nstate
      hartree_ii(istate,ispin)  = DOT_PRODUCT( c_matrix_restart(:,istate,ispin) , MATMUL( hamiltonian_hartree(:,:) , c_matrix_restart(:,istate,ispin) ) )
      exchange_ii(istate,ispin) = DOT_PRODUCT( c_matrix_restart(:,istate,ispin) , MATMUL( hamiltonian_exx(:,:,ispin) , c_matrix_restart(:,istate,ispin) ) )
   enddo
 enddo

 call dump_out_energy('=== Hartree expectation value ===',nstate,nspin,occupation,hartree_ii)
 call dump_out_energy('=== Exchange expectation value ===',nstate,nspin,occupation,exchange_ii)

 deallocate(energy_restart,occupation_restart)
 call clean_deallocate('RESTART: C',c_matrix_restart)

end subroutine print_hartee_expectation


!=========================================================================
end module m_scf_loop
!=========================================================================
