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
 use m_inputparam
 use m_hamiltonian_wrapper

 integer,parameter,private :: BATCH_SIZE = 128

contains


!=========================================================================
subroutine scf_loop(is_restart,&
                    basis,&
                    nstate,m_ham,n_ham,m_c,n_c,&
                    x_matrix,s_matrix,&
                    hamiltonian_kinetic,hamiltonian_nucleus,&
                    occupation, &
                    energy, &
                    hamiltonian_fock,&
                    c_matrix)
 use m_tools
 use m_atoms
 use m_basis_set
 use m_scf
 use m_eri
 use m_eri_calculate
 use m_eri_ao_mo
 use m_dft_grid
 use m_spectral_function
 use m_hamiltonian
 use m_hamiltonian_sca
 use m_hamiltonian_buffer
 use m_selfenergy_tools
 implicit none

!=====
 logical,intent(in)                 :: is_restart
 type(basis_set),intent(in)         :: basis
 integer,intent(in)                 :: nstate,m_ham,n_ham,m_c,n_c
 real(dp),intent(in)                :: x_matrix(m_c,n_c)
 real(dp),intent(in)                :: s_matrix(m_ham,n_ham)
 real(dp),intent(in)                :: hamiltonian_kinetic(m_ham,n_ham)
 real(dp),intent(in)                :: hamiltonian_nucleus(m_ham,n_ham)
 real(dp),intent(inout)             :: occupation(nstate,nspin)
 real(dp),intent(out)               :: energy(nstate,nspin)
 real(dp),allocatable,intent(inout) :: hamiltonian_fock(:,:,:)
 real(dp),allocatable,intent(inout) :: c_matrix(:,:,:)
!=====
 type(spectral_function)    :: wpol
 type(energy_contributions) :: en_dm_corr
 character(len=64)       :: restart_filename
 logical                 :: is_converged,stopfile_found,density_matrix_found
 integer                 :: file_density_matrix
 integer                 :: ispin,iscf,istate
 integer                 :: restart_type
 real(dp)                :: energy_tmp
 real(dp),allocatable    :: p_matrix(:,:,:)
 real(dp),allocatable    :: hamiltonian(:,:,:)
 real(dp),allocatable    :: hamiltonian_hartree(:,:)
 real(dp),allocatable    :: hamiltonian_exx(:,:,:)
 real(dp),allocatable    :: hamiltonian_xc(:,:,:)
 real(dp),allocatable    :: matrix_tmp(:,:,:)
 real(dp),allocatable    :: p_matrix_corr(:,:,:)
 real(dp),allocatable    :: hamiltonian_hartree_corr(:,:)
 real(dp),allocatable    :: hamiltonian_exx_corr(:,:,:)
 real(dp),allocatable    :: hfock_restart(:,:,:)
 real(dp),allocatable    :: c_matrix_restart(:,:,:)
 real(dp),allocatable    :: c_matrix_tmp(:,:,:)
 real(dp),allocatable    :: occupation_tmp(:,:)
 real(dp)                :: hartree_ii(nstate,nspin),exchange_ii(nstate,nspin)
 real(dp)                :: energy_restart(nstate,nspin)
!=====
 real(dp),allocatable    :: energy_exx(:,:)
 real(dp),allocatable    :: c_matrix_exx(:,:,:)
!=====


 call start_clock(timing_scf)

 ! Old Fock operator will be updated
 ! Get rid of it!
 call clean_deallocate('Fock operator F',hamiltonian_fock) ! Never distributed

 !
 ! Initialize the SCF mixing procedure
 call init_scf(m_ham,n_ham,m_c,n_c,basis%nbf,nstate)

 !
 ! Allocate the main arrays
 call clean_allocate('Total Hamiltonian H',hamiltonian,m_ham,n_ham,nspin)
 call clean_allocate('Hartree potential Vh',hamiltonian_hartree,m_ham,n_ham)
 call clean_allocate('Exchange operator Sigx',hamiltonian_exx,m_ham,n_ham,nspin)
 call clean_allocate('XC operator Vxc',hamiltonian_xc,m_ham,n_ham,nspin)
 call clean_allocate('Density matrix P',p_matrix,m_ham,n_ham,nspin)


 if( calc_type%is_dft ) then
   !
   ! Setup the grids for the quadrature of DFT potential/energy
   call init_dft_grid(basis,grid_level,dft_xc_needs_gradient,.TRUE.,BATCH_SIZE)
 endif

 !
 ! Setup the density matrix: p_matrix
 if( parallel_ham ) then
   call setup_density_matrix_sca(c_matrix,occupation,p_matrix)
 else
   call setup_density_matrix(c_matrix,occupation,p_matrix)
 endif


 !
 ! Start the big scf loop
 !
 do iscf=1,nscf
   write(stdout,'(/,a)') '-------------------------------------------'
   write(stdout,'(a,1x,i4,/)') ' *** SCF cycle No:',iscf


   if( cntxt_ham > 0 ) then
     en%kin  = SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
     en%nuc  = SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
   else
     en%kin  = 0.0_dp
     en%nuc  = 0.0_dp
   endif
   call xsum_trans(en%kin)
   call xsum_trans(en%nuc)

   !
   ! Setup kinetic and nucleus contributions (that are independent of the
   ! density matrix and therefore of spin channel)
   !
   hamiltonian(:,:,1) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
   if(nspin==2) hamiltonian(:,:,nspin) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)

   !
   ! Hartree contribution to the Hamiltonian
   !
   call calculate_hartree(basis,p_matrix,hamiltonian_hartree,eh=en%hart)

   ! calc_type%is_core is an inefficient way to get the Kinetic+Nucleus Hamiltonian
   if( calc_type%is_core ) then
     hamiltonian_hartree(:,:) = 0.0_dp
     en%hart = 0.0_dp
   endif
   do ispin=1,nspin
     hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) + hamiltonian_hartree(:,:)
   enddo


   !
   !  XC part of the Hamiltonian
   !
   hamiltonian_xc(:,:,:) = 0.0_dp
   en%exx_hyb = 0.0_dp

   !
   ! DFT XC potential is added here
   ! hamiltonian_xc is used as a temporary matrix
   if( calc_type%is_dft ) then

     if( parallel_ham ) then
       if( parallel_buffer ) then
         call dft_exc_vxc_buffer_sca(BATCH_SIZE,basis,occupation,c_matrix,hamiltonian_xc,en%xc)
       else
         call issue_warning('Exc calculation with SCALAPACK is not coded yet. Just skip it')
         hamiltonian_xc(:,:,:) = 0.0_dp
         en%xc = 0.0_dp
       endif
     else
       call dft_exc_vxc_batch(BATCH_SIZE,basis,occupation,c_matrix,hamiltonian_xc,en%xc)
     endif

   endif

   !
   ! LR Exchange contribution to the Hamiltonian
   ! Use hamiltonian_exx as a temporary matrix (no need to save it for later use)
   if(calc_type%need_exchange_lr) then

     call calculate_exchange_lr(basis,p_matrix,hamiltonian_exx,ex=energy_tmp,occupation=occupation,c_matrix=c_matrix)
     ! Rescale with alpha_hybrid_lr for range-separated hybrid functionals
     en%exx_hyb = en%exx_hyb + alpha_hybrid_lr * energy_tmp
     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:) * alpha_hybrid_lr

   endif

   !
   ! Exchange contribution to the Hamiltonian
   if( calc_type%need_exchange ) then

     call calculate_exchange(basis,p_matrix,hamiltonian_exx,ex=en%exx,occupation=occupation,c_matrix=c_matrix)

     ! Rescale with alpha_hybrid for hybrid functionals
     en%exx_hyb = en%exx_hyb + alpha_hybrid * en%exx
     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:) * alpha_hybrid
   endif


   !
   ! QSGW or COHSEX self energy
   if( ( calc_type%selfenergy_approx == GW .OR. calc_type%selfenergy_approx == COHSEX ) &
        .AND. calc_type%selfenergy_technique == QS  &
        .AND. ( iscf > 5 .OR. is_restart ) ) then

     if( parallel_ham ) call die('QSGW not implemented with parallel_ham')

     call init_spectral_function(nstate,occupation,0,wpol)
     call polarizability(.TRUE.,.TRUE.,basis,nstate,occupation,energy,c_matrix,en%rpa,wpol)

     if( ABS(en%rpa) > 1.e-6_dp) then
       en%tot = en%tot + en%rpa
       write(stdout,'(/,a,f19.10)') ' RPA Total energy (Ha): ',en%tot
     endif

     !
     ! Set the range of states on which to evaluate the self-energy
     call selfenergy_set_state_range(nstate,occupation)

     allocate(matrix_tmp(m_ham,n_ham,nspin))
     call gw_selfenergy_qs(nstate,basis,occupation,energy,c_matrix,s_matrix,wpol,matrix_tmp)

     call dump_out_matrix(.FALSE.,'=== Self-energy ===',basis%nbf,nspin,matrix_tmp)
     call destroy_spectral_function(wpol)

     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + matrix_tmp(:,:,:)
     deallocate(matrix_tmp)

   endif

   !
   ! QSPT2
   if( calc_type%selfenergy_approx == PT2 .AND. calc_type%selfenergy_technique == QS .AND. ( iscf > 5 .OR. is_restart ) ) then

     if( parallel_ham ) call die('QSPT2 not implemented with parallel_ham')

     !
     ! Set the range of states on which to evaluate the self-energy
     call selfenergy_set_state_range(nstate,occupation)

     allocate(matrix_tmp(m_ham,n_ham,nspin))
     call pt2_selfenergy_qs(nstate,basis,occupation,energy,c_matrix,s_matrix,matrix_tmp,en%mp2)

     write(stdout,'(a,2x,f19.10)') ' MP2 Energy       (Ha):',en%mp2
     write(stdout,*)
     en%tot = en%tot + en%mp2
     write(stdout,'(a,2x,f19.10)') ' MP2 Total Energy (Ha):',en%tot

     call dump_out_matrix(.FALSE.,'=== Self-energy ===',basis%nbf,nspin,matrix_tmp)

     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + matrix_tmp(:,:,:)
     deallocate(matrix_tmp)

   endif


   !
   ! Add the XC part of the hamiltonian to the total hamiltonian
   hamiltonian(:,:,:) = hamiltonian(:,:,:) + hamiltonian_xc(:,:,:)


   ! All the components of the energy have been calculated at this stage
   ! Sum up to get the total energy
   en%tot = en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx_hyb + en%xc

   ! Make sure all the MPI threads have the exact same Hamiltonian
   ! It helps stabilizing the SCF cycles in parallel
   if( .NOT. parallel_ham ) then
     call xsum_world(hamiltonian)
     hamiltonian(:,:,:) = hamiltonian(:,:,:) / REAL(nproc_world,dp)
   endif

   !
   ! If requested, the level shifting procedure is triggered:
   ! All the unoccupied states are penalized with an energy =  level_shifting_energy
   if( level_shifting_energy > 1.0e-6_dp ) then
     if( parallel_ham ) call die('level_shifting: not implemented with parallel_ham')
     call level_shifting_up(s_matrix,c_matrix,occupation,level_shifting_energy,hamiltonian)
   endif

   ! DIIS or simple mixing on the hamiltonian
   call hamiltonian_prediction(s_matrix,x_matrix,p_matrix,hamiltonian)


   !
   ! Diagonalize the Hamiltonian H
   ! Generalized eigenvalue problem with overlap matrix S
   ! H \varphi = E S \varphi
   ! save the old eigenvalues
   if( parallel_ham ) then
     call diagonalize_hamiltonian_sca(desc_ham,hamiltonian,desc_c,x_matrix,energy,c_matrix)
   else
     ! This subroutine works with or without scalapack
     call diagonalize_hamiltonian_scalapack(hamiltonian,x_matrix,energy,c_matrix)
   endif

   !
   ! When level_shifting is used, the unoccupied state energies have to be brought
   ! back to their original value,
   ! So that the "physical" energies are written down
   if( level_shifting_energy > 1.0e-6_dp ) then
     if( parallel_ham ) call die('level_shifting: not implemented with parallel_ham')
     call level_shifting_down(s_matrix,c_matrix,occupation,level_shifting_energy,energy,hamiltonian)
   endif

   call dump_out_energy('=== Energies ===',nstate,nspin,occupation,energy)

   call output_new_homolumo('gKS',nstate,occupation,energy,1,nstate)


   !
   ! Output the total energy and its components
   write(stdout,*)
   write(stdout,'(a25,1x,f19.10)') 'Nucleus-Nucleus (Ha):',en%nuc_nuc
   write(stdout,'(a25,1x,f19.10)') 'Kinetic Energy  (Ha):',en%kin
   write(stdout,'(a25,1x,f19.10)') 'Nucleus Energy  (Ha):',en%nuc
   write(stdout,'(a25,1x,f19.10)') 'Hartree Energy  (Ha):',en%hart
   if(calc_type%need_exchange) then
     write(stdout,'(a25,1x,f19.10)') 'Exchange Energy (Ha):',en%exx_hyb
   endif
   if( calc_type%is_dft ) then
     write(stdout,'(a25,1x,f19.10)') 'XC Energy       (Ha):',en%xc
   endif
   write(stdout,'(/,a25,1x,f19.10,/)') 'Total Energy    (Ha):',en%tot


   ! If fractional occupancies are allowed, then recalculate the occupations
   if( temperature > 1.0e-8_dp ) then
     call set_occupation(temperature,electrons,magnetization,energy,occupation)
   endif

   !
   ! Setup the new density matrix: p_matrix
   ! Save the old one for the convergence criterium
   if( parallel_ham ) then
     call setup_density_matrix_sca(c_matrix,occupation,p_matrix)
   else
     call setup_density_matrix(c_matrix,occupation,p_matrix)
   endif


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
   ! Skip writing when parallel_ham since all the information is not available on the master proc.
   if( print_restart_ .AND. .NOT. parallel_ham ) then
     call clean_allocate('Fock operator F',hamiltonian_fock,basis%nbf,basis%nbf,nspin)
     call get_fock_operator(hamiltonian,hamiltonian_xc,hamiltonian_exx,hamiltonian_fock)
     call write_restart(SMALL_RESTART,basis,nstate,occupation,c_matrix,energy,hamiltonian_fock)
     call clean_deallocate('Fock operator F',hamiltonian_fock)
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


 !
 ! Get the exchange operator if not already calculated
 !
 if( ABS(en%exx) < 1.0e-6_dp ) then
   call calculate_exchange(basis,p_matrix,hamiltonian_exx,ex=en%exx,occupation=occupation,c_matrix=c_matrix)
 endif


 !
 ! Print out some expectation values if requested
 !
 if( print_hartree_ ) then

   call clean_allocate('RESTART: C',c_matrix_restart,basis%nbf,nstate,nspin)
   call clean_allocate('RESTART: H',hfock_restart,basis%nbf,basis%nbf,nspin)

   restart_filename='RESTART_TEST'
   call read_restart(restart_type,basis,nstate,occupation,c_matrix_restart,energy_restart,hfock_restart,restart_filename)

   if( restart_type /= NO_RESTART ) then
     write(stdout,'(1x,a,a)') 'RESTART file read: ',restart_filename
     do ispin=1,nspin
       do istate=1,nstate
          hartree_ii(istate,ispin)  = DOT_PRODUCT( c_matrix_restart(:,istate,ispin) , MATMUL( hamiltonian_hartree(:,:) , c_matrix_restart(:,istate,ispin) ) )
          exchange_ii(istate,ispin) = DOT_PRODUCT( c_matrix_restart(:,istate,ispin) , MATMUL( hamiltonian_exx(:,:,ispin) , c_matrix_restart(:,istate,ispin) ) )
       enddo
     enddo

   else
     write(stdout,'(1x,a)') 'no RESTART file read'
     do ispin=1,nspin
       do istate=1,nstate
          hartree_ii(istate,ispin) =  DOT_PRODUCT( c_matrix(:,istate,ispin) , MATMUL( hamiltonian_hartree(:,:) , c_matrix(:,istate,ispin) ) )
          exchange_ii(istate,ispin) =  DOT_PRODUCT( c_matrix(:,istate,ispin) , MATMUL( hamiltonian_exx(:,:,ispin) , c_matrix(:,istate,ispin) ) )
       enddo
     enddo

   endif
   call dump_out_energy('=== Hartree expectation value ===',nstate,nspin,occupation,hartree_ii)
   call dump_out_energy('=== Exchange expectation value ===',nstate,nspin,occupation,exchange_ii)

   call clean_deallocate('RESTART: C',c_matrix_restart)
   call clean_deallocate('RESTART: H',hfock_restart)

 endif

 !
 ! Is there a correlated density matrix to be read or to be calculated
 if( read_fchk /= 'NO' .OR. TRIM(pt_density_matrix) /= 'NO' .OR. use_correlated_density_matrix_ ) then

   call clean_allocate('Correlated density matrix',p_matrix_corr,basis%nbf,basis%nbf,nspin)
   call clean_allocate('Correlated Hartree potential',hamiltonian_hartree_corr,basis%nbf,basis%nbf)
   call clean_allocate('Correlated exchange operator',hamiltonian_exx_corr,basis%nbf,basis%nbf,nspin)
   p_matrix_corr(:,:,:) = 0.0_dp

   !
   ! Three possibilities: read_fchk , pt_density_matrix, DENSITY_MATRIX
   !

   ! Option 1:
   ! Is there a Gaussian formatted checkpoint file to be read?
   if( read_fchk /= 'NO') call read_gaussian_fchk(read_fchk,'gaussian.fchk',basis,p_matrix_corr)

   ! Option 2:
   ! Calculate a MBPT density matrix if requested
   select case(TRIM(pt_density_matrix))
   case('ONE-RING')
     ! This keyword calculates the 1-ring density matrix as it is derived in PT2 theory
     call selfenergy_set_state_range(nstate,occupation)
     call fock_density_matrix(nstate,basis,occupation,energy,c_matrix,hamiltonian_exx,hamiltonian_xc,p_matrix_corr)
     call onering_density_matrix(nstate,basis,occupation,energy,c_matrix,p_matrix_corr)
   case('PT2')
     ! This keyword calculates the PT2 density matrix as it is derived in PT2 theory (differs from MP2 density matrix)
     call selfenergy_set_state_range(nstate,occupation)
     call fock_density_matrix(nstate,basis,occupation,energy,c_matrix,hamiltonian_exx,hamiltonian_xc,p_matrix_corr)
     call pt2_density_matrix(nstate,basis,occupation,energy,c_matrix,p_matrix_corr)
   case('GW','G0W0')
     ! This keyword calculates the GW density matrix as it is derived in the new GW theory
     call init_spectral_function(nstate,occupation,0,wpol)
     call polarizability(.TRUE.,.TRUE.,basis,nstate,occupation,energy,c_matrix,en%rpa,wpol)
     call selfenergy_set_state_range(nstate,occupation)
     call fock_density_matrix(nstate,basis,occupation,energy,c_matrix,hamiltonian_exx,hamiltonian_xc,p_matrix_corr)
     call gw_density_matrix(nstate,basis,occupation,energy,c_matrix,wpol,p_matrix_corr)
     call destroy_spectral_function(wpol)
   case('GW_IMAGINARY','G0W0_IMAGINARY')
     ! This keyword calculates the GW density matrix as it is derived in the new GW theory
     ! using an imaginary axis integral
     call init_spectral_function(nstate,occupation,nomega_imag,wpol)
     call polarizability_grid_scalapack(basis,nstate,occupation,energy,c_matrix,en%rpa,wpol)
     call selfenergy_set_state_range(nstate,occupation)
     call fock_density_matrix(nstate,basis,occupation,energy,c_matrix,hamiltonian_exx,hamiltonian_xc,p_matrix_corr)
     call gw_density_matrix_imag(nstate,basis,occupation,energy,c_matrix,wpol,p_matrix_corr)
     call destroy_spectral_function(wpol)
   end select


   ! Option 3:
   ! If no p_matrix_corr is present yet, then try to read it from a DENSITY_MATRIX file
   if( ALL( ABS(p_matrix_corr(:,:,:)) < 0.01_dp ) ) then
     inquire(file='DENSITY_MATRIX',exist=density_matrix_found)
     if( density_matrix_found) then
       write(stdout,'(/,1x,a)') 'Reading a MOLGW density matrix file: DENSITY_MATRIX'
       open(newunit=file_density_matrix,file='DENSITY_MATRIX',form='unformatted',action='read')
       do ispin=1,nspin
         read(file_density_matrix) p_matrix_corr(:,:,ispin)
       enddo
       close(file_density_matrix)
     else
       call die('m_scf_loop: no correlated density matrix read or calculated though input file suggests you really want one')
     endif

   endif

   if( print_hartree_ .OR. use_correlated_density_matrix_ ) then

     en_dm_corr%nuc_nuc = en%nuc_nuc
     en_dm_corr%kin = SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix_corr(:,:,:),DIM=3) )
     en_dm_corr%nuc = SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix_corr(:,:,:),DIM=3) )

     call calculate_hartree(basis,p_matrix_corr,hamiltonian_hartree_corr,eh=en_dm_corr%hart)

     call calculate_exchange(basis,p_matrix_corr,hamiltonian_exx_corr,ex=en_dm_corr%exx)

     en_dm_corr%tot = en_dm_corr%nuc_nuc + en_dm_corr%kin + en_dm_corr%nuc +  en_dm_corr%hart + en_dm_corr%exx
     write(stdout,'(/,1x,a)') 'Energies from correlated density matrix'
     write(stdout,'(a25,1x,f19.10)')   'Kinetic Energy (Ha):',en_dm_corr%kin
     write(stdout,'(a25,1x,f19.10)')   'Nucleus Energy (Ha):',en_dm_corr%nuc
     write(stdout,'(a25,1x,f19.10)')   'Hartree Energy (Ha):',en_dm_corr%hart
     write(stdout,'(a25,1x,f19.10)')  'Exchange Energy (Ha):',en_dm_corr%exx
     write(stdout,'(a25,1x,f19.10)') 'Total EXX Energy (Ha):',en_dm_corr%tot

     do ispin=1,nspin
       do istate=1,nstate
          hartree_ii(istate,ispin)  =  DOT_PRODUCT( c_matrix(:,istate,ispin) , &
                                                    MATMUL( hamiltonian_hartree_corr(:,:) , c_matrix(:,istate,ispin) ) )
          exchange_ii(istate,ispin) =  DOT_PRODUCT( c_matrix(:,istate,ispin) , &
                                                    MATMUL( hamiltonian_exx_corr(:,:,ispin) , c_matrix(:,istate,ispin) ) )
       enddo
     enddo
     call dump_out_energy('=== Hartree expectation value from correlated density matrix ===',nstate,nspin,occupation,hartree_ii)
     call dump_out_energy('=== Exchange expectation value from correlated density matrix ===',nstate,nspin,occupation,exchange_ii)
   endif

   if( print_multipole_ .OR. print_cube_ ) then
     allocate(c_matrix_tmp,MOLD=p_matrix)
     allocate(occupation_tmp(basis%nbf,nspin))
     call get_c_matrix_from_p_matrix(p_matrix_corr,c_matrix_tmp,occupation_tmp)
     if( print_multipole_ ) then
       call static_dipole(basis%nbf,basis,occupation_tmp,c_matrix_tmp)
       call static_quadrupole(basis%nbf,basis,occupation_tmp,c_matrix_tmp)
     endif
     if( print_cube_ ) then
       call plot_cube_wfn('MBPT',basis%nbf,basis,occupation_tmp,c_matrix_tmp)
     endif
     deallocate(c_matrix_tmp)
     deallocate(occupation_tmp)
   endif

   if( use_correlated_density_matrix_ ) then
     !
     ! Since the density matrix p_matrix is updated,
     ! one needs to recalculate the hartree and the exchange potentials
     ! let us include the old hartree in hamiltonian_xc and the new one in hamiltonian_exchange
     do ispin=1,nspin
       hamiltonian_xc(:,:,ispin)  = hamiltonian_xc(:,:,ispin) + hamiltonian_hartree(:,:)
       hamiltonian_exx(:,:,ispin) = hamiltonian_exx_corr(:,:,ispin) + hamiltonian_hartree_corr(:,:)
     enddo

   endif

   write(stdout,*)
   call clean_deallocate('Correlated density matrix',p_matrix_corr)
   call clean_deallocate('Correlated Hartree potential',hamiltonian_hartree_corr)
   call clean_deallocate('Correlated exchange operator',hamiltonian_exx_corr)

 endif

 !
 ! Form the Fock matrix and store it
 !
 call clean_allocate('Fock operator F',hamiltonian_fock,basis%nbf,basis%nbf,nspin)
 call get_fock_operator(hamiltonian,hamiltonian_xc,hamiltonian_exx,hamiltonian_fock)

 if( print_density_matrix_ .AND. is_iomaster ) then
   write(stdout,'(1x,a)') 'Write DENSITY_MATRIX_GKS file'
   open(newunit=file_density_matrix,file='DENSITY_MATRIX_GKS',form='unformatted',action='write')
   do ispin=1,nspin
     write(file_density_matrix) p_matrix(:,:,ispin)
   enddo
   close(file_density_matrix)
 endif

 !
 ! Cleanly deallocate the arrays
 !
 call clean_deallocate('Density matrix P',p_matrix)
 call clean_deallocate('Total Hamiltonian H',hamiltonian)
 call clean_deallocate('Hartree potential Vh',hamiltonian_hartree)
 call clean_deallocate('Exchange operator Sigx',hamiltonian_exx)
 call clean_deallocate('XC operator Vxc',hamiltonian_xc)

 write(stdout,'(/,/,a25,1x,f19.10,/)') 'SCF Total Energy (Ha):',en%tot
 write(stdout,'(a25,1x,f19.10)')       '      EXX Energy (Ha):',en%exx
 write(stdout,'(a25,1x,f19.10)')       'Total EXX Energy (Ha):',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx


 !
 ! Deallocate the buffer here
 if( parallel_ham .AND. parallel_buffer ) call destroy_parallel_buffer()
 !
 ! At this point, all procs get the complete c_matrix
 !
 call form_c_matrix_global(basis%nbf,nstate,c_matrix)


 !
 ! Single excitation term
 !
 call single_excitations(nstate,basis%nbf,energy,occupation,c_matrix,hamiltonian_fock,en%se)
 if( ABS(en%se) > 1.0e-6_dp )  write(stdout,'(a25,1x,f19.10)') 'Singles correction (Ha):',en%se
 write(stdout,'(a25,1x,f19.10,/)')   'Est. HF Energy (Ha):',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%se

 !
 ! Evaluate spin contamination
 if( .NOT. parallel_ham ) call evaluate_s2_operator(occupation,c_matrix,s_matrix)


 !
 ! Big RESTART file written if converged
 !
 if( is_converged .AND. print_bigrestart_ ) then
   call write_restart(BIG_RESTART,basis,nstate,occupation,c_matrix,energy,hamiltonian_fock)
 else
   if( print_restart_ ) then
     call write_restart(SMALL_RESTART,basis,nstate,occupation,c_matrix,energy,hamiltonian_fock)
   endif
 endif


 call stop_clock(timing_scf)

end subroutine scf_loop


!=========================================================================
subroutine calculate_hamiltonian_hxc(basis,nstate,occupation,c_matrix,p_matrix,hamiltonian_hxc,ehxc)
 use m_scalapack
 use m_basis_set
 use m_hamiltonian
 use m_hamiltonian_sca
 use m_hamiltonian_buffer
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(:,:,:)
 real(dp),intent(in)        :: p_matrix(:,:,:)
 real(dp),intent(out)       :: hamiltonian_hxc(:,:,:)
 real(dp),intent(out)       :: ehxc
!=====
 integer              :: ispin
 real(dp),allocatable :: hamiltonian_tmp(:,:)
 real(dp),allocatable :: hamiltonian_spin_tmp(:,:,:)
 real(dp)             :: ehart,exc,eexx,eexx_hyb
!=====


 allocate(hamiltonian_tmp,MOLD=hamiltonian_hxc(:,:,1))
 allocate(hamiltonian_spin_tmp,MOLD=hamiltonian_hxc(:,:,:))

 !
 ! Hartree contribution to the Hamiltonian
 !
 call calculate_hartree(basis,p_matrix,hamiltonian_tmp,eh=ehart)

 do ispin=1,nspin
   hamiltonian_hxc(:,:,ispin) = hamiltonian_tmp(:,:)
 enddo


 !
 !  XC part of the Hamiltonian
 !

 !
 ! DFT XC potential is added here
 !
 if( calc_type%is_dft ) then
   hamiltonian_spin_tmp(:,:,:) = 0.0_dp

   if( parallel_ham ) then
     if( parallel_buffer ) then
       call dft_exc_vxc_buffer_sca(BATCH_SIZE,basis,occupation,c_matrix,hamiltonian_spin_tmp,exc)
     else
       call issue_warning('Exc calculation with SCALAPACK is not coded yet. Just skip it')
       hamiltonian_spin_tmp(:,:,:) = 0.0_dp
       exc = 0.0_dp
     endif
   else
     call dft_exc_vxc_batch(BATCH_SIZE,basis,occupation,c_matrix,hamiltonian_spin_tmp,exc)
   endif

   hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_tmp(:,:,:)
 endif


 !
 ! LR Exchange contribution to the Hamiltonian
 !
 if(calc_type%need_exchange_lr) then
   hamiltonian_spin_tmp(:,:,:) = 0.0_dp

   call calculate_exchange_lr(basis,p_matrix,hamiltonian_spin_tmp,ex=eexx,occupation=occupation,c_matrix=c_matrix)
   ! Rescale with alpha_hybrid_lr for range-separated hybrid functionals
   eexx_hyb = alpha_hybrid_lr * eexx
   hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_tmp(:,:,:) * alpha_hybrid_lr

 endif


 !
 ! Exchange contribution to the Hamiltonian
 !
 if( calc_type%need_exchange ) then
   hamiltonian_spin_tmp(:,:,:) = 0.0_dp

   call calculate_exchange(basis,p_matrix,hamiltonian_spin_tmp,ex=eexx,occupation=occupation,c_matrix=c_matrix)
   ! Rescale with alpha_hybrid for hybrid functionals
   eexx_hyb = eexx_hyb + alpha_hybrid * eexx
   hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_tmp(:,:,:) * alpha_hybrid

 endif

 ehxc = ehart + eexx_hyb + exc


end subroutine calculate_hamiltonian_hxc


!=========================================================================
subroutine calculate_hamiltonian_hxc_ri_cmplx(basis,                  &
                                              nstate,                 &
                                              nocc,                   &
                                              m_ham,                  &
                                              n_ham,                  &
                                              m_c,                    &
                                              n_c,                    &
                                              occupation,             &
                                              c_matrix_cmplx,         &
                                              p_matrix_cmplx,         &
                                              hamiltonian_hxc_cmplx)
 use m_scalapack
 use m_basis_set
 use m_hamiltonian
 use m_hamiltonian_cmplx
 use m_hamiltonian_sca
 use m_hamiltonian_buffer
 use m_tools,only: matrix_trace_cmplx
 use m_scf
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: m_ham,n_ham
 integer,intent(in)         :: nstate
 integer,intent(in)         :: nocc
 integer,intent(in)         :: m_c,n_c
 real(dp),intent(in)        :: occupation(nstate,nspin)
 complex(dp),intent(in)    :: c_matrix_cmplx(m_c,n_c,nspin)
 complex(dp),intent(in)    :: p_matrix_cmplx(m_ham,n_ham,nspin)
 complex(dp),intent(out)   :: hamiltonian_hxc_cmplx(m_ham,n_ham,nspin)
!=====
 integer         :: ispin
 real(dp)        :: p_matrix(m_ham,n_ham,nspin)
 real(dp)        :: hamiltonian_tmp(m_ham,n_ham,nspin)
!=====

 en%hart    = 0.0_dp
 en%xc      = 0.0_dp
 en%exx     = 0.0_dp
 en%exx_hyb = 0.0_dp

! if ( parallel_ham ) call die('parallel_ham not yet implemented for tddft propagator')

 ! Initialize real arrays

 p_matrix=REAL(p_matrix_cmplx,dp)

 hamiltonian_hxc_cmplx = ( 0.0_dp , 0.0_dp )

 !
 ! Exchange contribution to the Hamiltonian
 !
 if( calc_type%need_exchange ) then
   call setup_exchange_versatile_ri_cmplx(occupation,c_matrix_cmplx,p_matrix_cmplx,hamiltonian_hxc_cmplx,en%exx)

   ! Rescale with alpha_hybrid for hybrid functionals
   en%exx_hyb = alpha_hybrid * en%exx
   hamiltonian_hxc_cmplx(:,:,:) = hamiltonian_hxc_cmplx(:,:,:) * alpha_hybrid
 endif


   !
   ! Hartree contribution to the Hamiltonian
   ! Hartree contribution is real and depends only on real(p_matrix)
   !
   !call calculate_hartree(basis,p_matrix,hamiltonian_tmp(:,:,1),eh=en%hart)
   call start_clock(timing_tddft_hartree)
   call setup_hartree_versatile_ri(p_matrix,hamiltonian_tmp(:,:,1),en%hart)
   call stop_clock(timing_tddft_hartree)

 do ispin=1,nspin
   hamiltonian_hxc_cmplx(:,:,ispin) = hamiltonian_hxc_cmplx(:,:,ispin) + hamiltonian_tmp(:,:,1)
 enddo

 !
 !  XC part of the Hamiltonian
 !

 !
 ! DFT XC potential is added here
 !
 if( calc_type%is_dft ) then
   call dft_exc_vxc_batch(BATCH_SIZE,basis,occupation,c_matrix_cmplx,hamiltonian_tmp,en%xc)

   hamiltonian_hxc_cmplx(:,:,:) = hamiltonian_hxc_cmplx(:,:,:) + hamiltonian_tmp(:,:,:)
 endif

! write(file_time_data,"(6(x,e16.10,2x),'    ')",advance='no') enuc,ekin,ehart, eexx_hyb,exc, enuc+ekin+ehart+eexx_hyb+exc
 !
 ! LR Exchange contribution to the Hamiltonian
 !
 ! if(calc_type%need_exchange_lr) then
 !   hamiltonian_spin_tmp(:,:,:) = 0.0_dp
 !
 !     call setup_exchange_longrange_ri(basis%nbf,nstate,occupation,c_matrix,p_matrix,hamiltonian_spin_tmp,eexx)
 !
 !   ! Rescale with alpha_hybrid_lr for range-separated hybrid functionals
 !   eexx_hyb = alpha_hybrid_lr * eexx
 !   hamiltonian_hxc(:,:,:) = hamiltonian_hxc(:,:,:) + hamiltonian_spin_tmp(:,:,:) * alpha_hybrid_lr
 ! endif


end subroutine  calculate_hamiltonian_hxc_ri_cmplx


!=========================================================================
subroutine get_fock_operator(hamiltonian,hamiltonian_xc,hamiltonian_exx,hamiltonian_fock)
 use m_mpi
 use m_scalapack
 implicit none

 real(dp),intent(in)    :: hamiltonian(:,:,:)
 real(dp),intent(in)    :: hamiltonian_xc(:,:,:)
 real(dp),intent(in)    :: hamiltonian_exx(:,:,:)
 real(dp),intent(out)   :: hamiltonian_fock(:,:,:)
!=====
 real(dp),allocatable   :: hfock_local(:,:,:)
!=====

 if( parallel_ham ) then

   call clean_allocate('Local Fock operator F',hfock_local,SIZE(hamiltonian,DIM=1),SIZE(hamiltonian,DIM=2),nspin)
   hfock_local(:,:,:) = hamiltonian(:,:,:) - hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:)
   call gather_distributed_copy(desc_ham,hfock_local,hamiltonian_fock)
   call clean_deallocate('Local Fock operator F',hfock_local)

 else
   hamiltonian_fock(:,:,:) = hamiltonian(:,:,:) - hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:)

 endif

end subroutine get_fock_operator


!=========================================================================
subroutine form_c_matrix_global(nbf,nstate,c_matrix)
 use m_mpi
 use m_scalapack
 implicit none

 integer,intent(in)                 :: nbf,nstate
 real(dp),allocatable,intent(inout) :: c_matrix(:,:,:)
!=====
 real(dp),allocatable :: c_matrix_local(:,:,:)
!=====

 if( .NOT. parallel_ham ) return    ! Nothing to do

 write(stdout,'(/,1x,a)') 'Form the C matrix on all procs'


 call clean_allocate('Local wfn coeff C',c_matrix_local,SIZE(c_matrix,DIM=1),SIZE(c_matrix,DIM=2),nspin)
! if( cntxt_ham > 0 ) then
 c_matrix_local(:,:,:) = c_matrix(:,:,:)
! endif
 call clean_deallocate('Wavefunctions C',c_matrix)
 call clean_allocate('Wavefunctions C',c_matrix,nbf,nstate,nspin)

 call gather_distributed_copy(desc_c,c_matrix_local,c_matrix)

 call clean_deallocate('Local wfn coeff C',c_matrix_local)


 write(stdout,'(1x,a)') 'C matrix on all procs formed'

end subroutine form_c_matrix_global



!=========================================================================
end module m_scf_loop
!=========================================================================
