!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the main SCF loop for Hartree-Fock or Kohn-Sham
!
!=========================================================================
module m_scf_loop
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_inputparam


contains


!=========================================================================
subroutine scf_loop(is_restart,& 
                    basis,auxil_basis,&
                    nstate,m_ham,n_ham,m_c,n_c,&
                    s_matrix_sqrt_inv,s_matrix,&
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
 use m_eri_calculate_lr
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
 type(basis_set),intent(in)         :: auxil_basis
 integer,intent(in)                 :: nstate,m_ham,n_ham,m_c,n_c
 real(dp),intent(in)                :: s_matrix_sqrt_inv(m_c,n_c)
 real(dp),intent(in)                :: s_matrix(m_ham,n_ham)
 real(dp),intent(in)                :: hamiltonian_kinetic(m_ham,n_ham)
 real(dp),intent(in)                :: hamiltonian_nucleus(m_ham,n_ham)
 real(dp),intent(inout)             :: occupation(nstate,nspin)
 real(dp),intent(out)               :: energy(nstate,nspin)
 real(dp),allocatable,intent(inout) :: hamiltonian_fock(:,:,:)
 real(dp),allocatable,intent(inout) :: c_matrix(:,:,:)
!=====
 type(spectral_function) :: wpol
 logical                 :: is_converged,stopfile_found
 integer                 :: ispin,iscf,istate
 real(dp)                :: energy_tmp
 real(dp),allocatable    :: p_matrix(:,:,:)
 real(dp),allocatable    :: p_matrix_sqrt(:,:,:),p_matrix_occ(:,:)
 real(dp),allocatable    :: hamiltonian(:,:,:)
 real(dp),allocatable    :: hamiltonian_hartree(:,:)
 real(dp),allocatable    :: hamiltonian_exx(:,:,:)
 real(dp),allocatable    :: hamiltonian_xc(:,:,:)
 real(dp),allocatable    :: matrix_tmp(:,:,:)
 real(dp),allocatable    :: p_matrix_old(:,:,:)
 real(dp),allocatable    :: self_energy_old(:,:,:)
 real(dp),allocatable    :: energy_exx(:,:)
 real(dp),allocatable    :: c_matrix_exx(:,:,:)
!=============================


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
 call clean_allocate('Previous density matrix Pold',p_matrix_old,m_ham,n_ham,nspin)
 call clean_allocate('Density matrix sqrt P^{1/2}',p_matrix_sqrt,m_ham,n_ham,nspin)
 allocate(p_matrix_occ(basis%nbf,nspin))

 if( calc_type%is_dft ) then
   !
   ! Setup the grids for the quadrature of DFT potential/energy
   call init_dft_grid(grid_level)
   ! The following is coded but not used... yet!
!   call setup_bf_radius(basis)
 endif

 !
 ! Setup the density matrix: p_matrix
 if( parallel_ham ) then
   call setup_density_matrix_sca(basis%nbf,nstate,m_c,n_c,c_matrix,occupation,m_ham,n_ham,p_matrix)
 else
   call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)
 endif


 !
 ! start the big scf loop
 !
 do iscf=1,nscf
   write(stdout,'(/,a)') '-------------------------------------------'
   write(stdout,'(a,1x,i4,/)') ' *** SCF cycle No:',iscf


   !
   ! Calculate the matrix square-root of the density matrix P
   if( parallel_ham ) then
     call setup_sqrt_density_matrix_sca(basis%nbf,m_ham,n_ham,p_matrix,p_matrix_sqrt,p_matrix_occ)
   else
     call setup_sqrt_density_matrix(basis%nbf,p_matrix,p_matrix_sqrt,p_matrix_occ)
   endif


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
   if( .NOT. is_full_auxil) then
     call setup_hartree(print_matrix_,basis%nbf,p_matrix,hamiltonian_hartree,en%hart)
   else
     if( parallel_ham ) then
       if( parallel_buffer ) then
         call setup_hartree_ri_buffer_sca(print_matrix_,basis%nbf,m_ham,n_ham,p_matrix,hamiltonian_hartree,en%hart)
       else
         call setup_hartree_ri_sca(print_matrix_,basis%nbf,m_ham,n_ham,p_matrix,hamiltonian_hartree,en%hart)
       endif
     else
       call setup_hartree_ri(print_matrix_,basis%nbf,p_matrix,hamiltonian_hartree,en%hart)
     endif
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
         call dft_exc_vxc_buffer_sca(m_ham,n_ham,basis,p_matrix_occ,p_matrix_sqrt,p_matrix,hamiltonian_xc,en%xc)
       else
         call issue_warning('Exc calculation with SCALAPACK is not coded yet. Just skip it')
         hamiltonian_xc(:,:,:) = 0.0_dp
         en%xc = 0.0_dp
       endif
     else
       call dft_exc_vxc(basis,p_matrix_occ,p_matrix_sqrt,p_matrix,hamiltonian_xc,en%xc)
     endif

   endif

   !
   ! LR Exchange contribution to the Hamiltonian
   ! Use hamiltonian_exx as a temporary matrix (no need to save it for later use)
   if(calc_type%need_exchange_lr) then

     if( .NOT. is_full_auxil) then
       call setup_exchange_longrange(print_matrix_,basis%nbf,p_matrix,hamiltonian_exx,energy_tmp)
     else
       call setup_exchange_longrange_ri(print_matrix_,basis%nbf,p_matrix_occ,p_matrix_sqrt,p_matrix,hamiltonian_exx,energy_tmp)
     endif
     ! Rescale with alpha_hybrid_lr for range-separated hybrid functionals
     en%exx_hyb = en%exx_hyb + alpha_hybrid_lr * energy_tmp
     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:) * alpha_hybrid_lr
   endif

   !
   ! Exchange contribution to the Hamiltonian
   if( calc_type%need_exchange ) then

     if( .NOT. is_full_auxil) then
       call setup_exchange(print_matrix_,basis%nbf,p_matrix,hamiltonian_exx,en%exx)
     else
       if( parallel_ham ) then
         if( parallel_buffer ) then
           call setup_exchange_ri_buffer_sca(print_matrix_,basis%nbf,m_ham,n_ham,p_matrix_occ,p_matrix_sqrt,p_matrix,hamiltonian_exx,en%exx)
         else
           call setup_exchange_ri_sca(print_matrix_,basis%nbf,m_ham,n_ham,p_matrix_occ,p_matrix_sqrt,p_matrix,hamiltonian_exx,en%exx)
         endif
       else
         call setup_exchange_ri(print_matrix_,basis%nbf,p_matrix_occ,p_matrix_sqrt,p_matrix,hamiltonian_exx,en%exx)
       endif
     endif
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

     call init_spectral_function(nstate,occupation,wpol)
     call polarizability(basis,auxil_basis,nstate,occupation,energy,c_matrix,en%rpa,wpol)

     if( ABS(en%rpa) > 1.e-6_dp) then
       en%tot = en%tot + en%rpa
       write(stdout,'(/,a,f19.10)') ' RPA Total energy (Ha): ',en%tot
     endif

     !
     ! Set the range of states on which to evaluate the self-energy
     call selfenergy_set_state_range(nstate,occupation)

     allocate(matrix_tmp(m_ham,n_ham,nspin))
     call gw_selfenergy_qs(nstate,basis,occupation,energy,c_matrix,s_matrix,wpol,matrix_tmp)

     if( .NOT. ALLOCATED(self_energy_old) ) then
       allocate(self_energy_old(basis%nbf,basis%nbf,nspin))
       self_energy_old(:,:,:) = 0.0_dp
     endif
     matrix_tmp(:,:,:) = alpha_mixing * matrix_tmp(:,:,:) + (1.0_dp-alpha_mixing) * self_energy_old(:,:,:)
     self_energy_old(:,:,:) = matrix_tmp(:,:,:)
     call dump_out_matrix(print_matrix_,'=== Self-energy ===',basis%nbf,nspin,matrix_tmp)
     call destroy_spectral_function(wpol)

     hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix_tmp(:,:,:)
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

     if( .NOT. ALLOCATED(self_energy_old) ) then
       allocate(self_energy_old(basis%nbf,basis%nbf,nspin))
       self_energy_old(:,:,:) = 0.0_dp
     endif
     matrix_tmp(:,:,:) = alpha_mixing * matrix_tmp(:,:,:) + (1.0_dp-alpha_mixing) * self_energy_old(:,:,:)
     self_energy_old(:,:,:) = matrix_tmp(:,:,:)
     call dump_out_matrix(print_matrix_,'=== Self-energy ===',basis%nbf,nspin,matrix_tmp)
  
     hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix_tmp(:,:,:)
     deallocate(matrix_tmp)

   endif

   !
   ! Add the XC part of the hamiltonian to the total hamiltonian
   hamiltonian(:,:,:) = hamiltonian(:,:,:) + hamiltonian_xc(:,:,:)
   

   !
   ! If requested, the level shifting procedure is triggered: 
   ! All the unoccupied states are penalized with an energy =  level_shifting_energy
   if( level_shifting_energy > 1.0e-6_dp ) then
     if( parallel_ham ) call die('level_shifting: not implemented with parallel_ham')
     call level_shifting(basis%nbf,nstate,s_matrix,c_matrix,occupation,level_shifting_energy,hamiltonian)
   endif


   ! DIIS or simple mixing on the hamiltonian
   call hamiltonian_prediction(s_matrix,s_matrix_sqrt_inv,p_matrix,hamiltonian)

  
   !
   ! Diagonalize the Hamiltonian H
   ! Generalized eigenvalue problem with overlap matrix S
   ! H \varphi = E S \varphi
   ! save the old eigenvalues
   if( parallel_ham ) then
     call diagonalize_hamiltonian_sca(nspin,basis%nbf,nstate,m_ham,n_ham,hamiltonian,s_matrix_sqrt_inv, &
                                      energy,m_c,n_c,c_matrix)
   else
     ! This subroutine works with or without scalapack
     call diagonalize_hamiltonian_scalapack(nspin,basis%nbf,nstate,hamiltonian,s_matrix_sqrt_inv,energy,c_matrix)
   endif

   !
   ! When level_shifting is used, the unoccupied state energies have to be brought
   ! back to their original value,
   ! So that the "physical" energies are written down
   if( level_shifting_energy > 1.0e-6_dp ) then
     do ispin=1,nspin
       do istate=1,nstate
         if( occupation(istate,ispin) < completely_empty ) then
           energy(istate,ispin) = energy(istate,ispin) - level_shifting_energy
         endif
       enddo
     enddo
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
   en%tot = en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx_hyb + en%xc
   write(stdout,'(/,a25,1x,f19.10,/)') 'Total Energy    (Ha):',en%tot


   !
   ! Setup the new density matrix: p_matrix
   ! Save the old one for the convergence criterium
   p_matrix_old(:,:,:) = p_matrix(:,:,:)
   if( parallel_ham ) then
     call setup_density_matrix_sca(basis%nbf,nstate,m_c,n_c,c_matrix,occupation,m_ham,n_ham,p_matrix)
   else
     call setup_density_matrix(basis%nbf,nstate,c_matrix,occupation,p_matrix)
   endif

   is_converged = check_converged(p_matrix_old,p_matrix)
   inquire(file='STOP',exist=stopfile_found)

   if( is_converged .OR. stopfile_found ) exit

   !
   ! Write down a "small" RESTART file at each step
   if( print_restart_ .AND. .NOT. parallel_ham ) then
     call clean_allocate('Fock operator F',hamiltonian_fock,basis%nbf,basis%nbf,nspin)
     call get_fock_operator(hamiltonian,hamiltonian_xc,hamiltonian_exx,hamiltonian_fock)
     call write_restart(SMALL_RESTART,basis,nstate,occupation,c_matrix,energy,hamiltonian_fock)
     call clean_deallocate('Fock operator F',hamiltonian_fock)
   endif

   ! Damping of the density matrix p_matrix
   call simple_mixing_p_matrix(p_matrix_old,p_matrix)
   
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
 if( .NOT. is_full_auxil) then
   if( ABS(en%exx) < 1.0e-6_dp ) call setup_exchange(print_matrix_,basis%nbf,p_matrix,hamiltonian_exx,en%exx)
 else
   if( ABS(en%exx) < 1.0e-6_dp ) then
     if( parallel_ham ) then
       if( parallel_buffer ) then
         call setup_exchange_ri_buffer_sca(print_matrix_,basis%nbf,m_ham,n_ham,p_matrix_occ,p_matrix_sqrt,p_matrix,hamiltonian_exx,en%exx)
       else
         call setup_exchange_ri_sca(print_matrix_,basis%nbf,m_ham,n_ham,p_matrix_occ,p_matrix_sqrt,p_matrix,hamiltonian_exx,en%exx)
       endif
     else
       call setup_exchange_ri(print_matrix_,basis%nbf,p_matrix_occ,p_matrix_sqrt,p_matrix,hamiltonian_exx,en%exx)
     endif
   endif
 endif

 !
 ! Obtain the Fock matrix and store it
 !
 call clean_allocate('Fock operator F',hamiltonian_fock,basis%nbf,basis%nbf,nspin)
 call get_fock_operator(hamiltonian,hamiltonian_xc,hamiltonian_exx,hamiltonian_fock)

 !
 ! Cleanly deallocate the arrays
 !
 call clean_deallocate('Density matrix P',p_matrix)
 call clean_deallocate('Density matrix sqrt P^{1/2}',p_matrix_sqrt)
 call clean_deallocate('Previous density matrix Pold',p_matrix_old)
 call clean_deallocate('Total Hamiltonian H',hamiltonian)
 call clean_deallocate('Hartree potential Vh',hamiltonian_hartree)
 call clean_deallocate('Exchange operator Sigx',hamiltonian_exx)
 call clean_deallocate('XC operator Vxc',hamiltonian_xc)
 if( ALLOCATED(self_energy_old) ) deallocate(self_energy_old)
 if( ALLOCATED(p_matrix_occ) )    deallocate(p_matrix_occ)

 write(stdout,'(/,/,a25,1x,f19.10,/)') 'SCF Total Energy (Ha):',en%tot
 write(stdout,'(a25,1x,f19.10)')       '      EXX Energy (Ha):',en%exx
 write(stdout,'(a25,1x,f19.10)')       'Total EXX Energy (Ha):',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx


 !
 ! At this point, all procs get the complete c_matrix
 !
 call form_c_matrix_global(basis%nbf,nstate,c_matrix)


 if( parallel_buffer ) then  
   !
   ! Spin contamination?
   call evaluate_s2_operator(basis%nbf,nstate,occupation,c_matrix,s_matrix)

   !
   ! Single excitation term
   !
   call single_excitations(nstate,basis%nbf,energy,occupation,c_matrix,hamiltonian_fock,en%se)
   write(stdout,'(a25,1x,f19.10)') 'Singles correction (Ha):',en%se
   write(stdout,'(a25,1x,f19.10,/)')   'Est. HF Energy (Ha):',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%se


   ! A dirty section for the Luttinger-Ward functional
   if(calc_type%selfenergy_approx==LW .OR. calc_type%selfenergy_approx==LW2 .OR. calc_type%selfenergy_approx==GSIGMA) then
     allocate(energy_exx(nstate,nspin))
     allocate(c_matrix_exx(basis%nbf,nstate,nspin))
     call issue_warning('ugly coding here write temp file fort.1000 and fort.1001')

     do ispin=1,nspin
       write(stdout,*) 'Diagonalization H_exx for spin channel',ispin
       call diagonalize_generalized_sym(basis%nbf,&
                                        hamiltonian_fock(:,:,ispin),s_matrix(:,:),&
                                        energy_exx(:,ispin),c_matrix_exx(:,:,ispin))
     enddo
     write(stdout,*) 'FBFB LW sum(      epsilon) + Eii -<vxc> - EH + Ex',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx
     write(stdout,*) 'FBFB LW sum(tilde epsilon) + Eii - EH - Ex       ',SUM( occupation(:,:)*energy_exx(:,:) ) + en%nuc_nuc - en%hart - en%exx
     open(1000,form='unformatted')
     do ispin=1,nspin
       do istate=1,nstate
         write(1000) c_matrix_exx(:,istate,ispin)
       enddo
     enddo
     close(1000)
     open(1001,form='unformatted')
     write(1001) energy_exx(:,:)
     close(1001)
     deallocate(energy_exx,c_matrix_exx)
   endif

 endif


 !
 ! Big RESTART file written if converged
 ! TODO: implement writing with a parallelized hamiltonian
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
 integer                :: rank_master
!=====

 if( parallel_ham ) then
   !
   ! Coding to be moved to low-level subroutines
   if( cntxt_ham > 0 .AND. iprow_ham == 0 .AND. ipcol_ham == 0 ) then
     rank_master = rank_world
   else
     rank_master = -1
   endif
   call xmax_world(rank_master)

   if( cntxt_ham > 0 ) then
     call clean_allocate('Local Fock operator F',hfock_local,SIZE(hamiltonian,DIM=1),SIZE(hamiltonian,DIM=2),nspin)
     hfock_local(:,:,:) = hamiltonian(:,:,:) - hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:)
     call gather_distributed_copy_spin(desc_ham,hfock_local,hamiltonian_fock)
     call clean_deallocate('Local Fock operator F',hfock_local)
   endif
   call xbcast_world(rank_master,hamiltonian_fock)

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
 integer              :: rank_master
!=====

 if( .NOT. parallel_ham ) return    ! Nothing to do

 write(stdout,'(/,1x,a)') 'Form the C matrix on all procs'
 !
 ! Coding to be moved to low-level subroutines
 if( cntxt_ham > 0 .AND. iprow_ham == 0 .AND. ipcol_ham == 0 ) then
   rank_master = rank_world
 else
   rank_master = -1
 endif
 call xmax_world(rank_master)

 if( cntxt_ham > 0 ) then
   call clean_allocate('Local wfn coeff C',c_matrix_local,SIZE(c_matrix,DIM=1),SIZE(c_matrix,DIM=2),nspin)
   c_matrix_local(:,:,:) = c_matrix(:,:,:)
   call clean_deallocate('Wavefunctions C',c_matrix)
   call clean_allocate('Wavefunctions C',c_matrix,nbf,nstate,nspin)

   call gather_distributed_copy_spin(desc_c,c_matrix_local,c_matrix)

   call clean_deallocate('Local wfn coeff C',c_matrix_local)

 else
   call clean_deallocate('Wavefunctions C',c_matrix)
   call clean_allocate('Wavefunctions C',c_matrix,nbf,nstate,nspin)
 endif

 call xbcast_world(rank_master,c_matrix)

 write(stdout,'(1x,a)') 'C matrix on all procs formed'

end subroutine form_c_matrix_global



!=========================================================================
end module m_scf_loop
!=========================================================================