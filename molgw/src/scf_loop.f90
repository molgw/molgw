!=========================================================================
! This file contains
! the main SCF loop for Hartree-Fock or Kohn-Sham
!=========================================================================
subroutine scf_loop(basis,prod_basis,auxil_basis,&
                    nstate,s_matrix_sqrt_inv,&
                    s_matrix,c_matrix,p_matrix,&
                    hamiltonian_kinetic,hamiltonian_nucleus,&
                    hamiltonian_hartree,hamiltonian_exx,hamiltonian_xc,&
                    occupation,energy)
 use m_definitions
 use m_timing
 use m_warning
 use m_inputparam
 use m_tools
 use m_scf
 use m_atoms
! use m_gaussian
 use m_basis_set
 use m_eri
 use m_dft_grid
 use m_spectral_function
 use m_timedependent
#ifdef _OPENMP
 use omp_lib
#endif
 implicit none

!=====
 type(basis_set),intent(in)         :: basis
 type(basis_set),intent(in)         :: prod_basis
 type(basis_set),intent(in)         :: auxil_basis
 integer,intent(in)                 :: nstate
 real(dp),intent(in)                :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(in)                :: s_matrix_sqrt_inv(basis%nbf,nstate)
 real(dp),intent(inout)             :: c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(inout)             :: p_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in)                :: hamiltonian_kinetic(basis%nbf,basis%nbf)
 real(dp),intent(in)                :: hamiltonian_nucleus(basis%nbf,basis%nbf)
 real(dp),intent(inout)             :: hamiltonian_hartree(basis%nbf,basis%nbf)
 real(dp),intent(inout)             :: hamiltonian_exx(basis%nbf,basis%nbf,nspin)
 real(dp),intent(inout)             :: hamiltonian_xc(basis%nbf,basis%nbf,nspin)
 real(dp),intent(inout)             :: occupation(basis%nbf,nspin)
 real(dp),intent(inout)             :: energy(basis%nbf,nspin)
!=====
 type(spectral_function) :: wpol
 integer                 :: ispin,iscf,istate
 integer                 :: fileunit,ncore
 character(len=100)      :: title
 real(dp)                :: energy_tmp
 real(dp),allocatable    :: ehomo(:),elumo(:)
 real(dp),allocatable    :: hamiltonian(:,:,:)
 real(dp),allocatable    :: hamiltonian_vxc(:,:,:)
 real(dp),allocatable    :: matrix_tmp(:,:,:)
 real(dp),allocatable    :: p_matrix_old(:,:,:)
 real(dp),allocatable    :: exchange_m_vxc_diag(:,:)
 real(dp),allocatable    :: self_energy_old(:,:,:)
 logical                 :: is_converged,stopfile_found,file_exists
 real(dp),allocatable    :: energy_exx(:,:)
 real(dp),allocatable    :: c_matrix_exx(:,:,:)
 real(dp),allocatable    :: occupation_tmp(:,:),p_matrix_tmp(:,:,:)
!=============================

 !
 ! Allocate the main arrays
 allocate(hamiltonian(basis%nbf,basis%nbf,nspin))
 allocate(matrix_tmp(basis%nbf,basis%nbf,nspin))
 allocate(p_matrix_old(basis%nbf,basis%nbf,nspin))
 allocate(exchange_m_vxc_diag(basis%nbf,nspin))
 allocate(self_energy_old(basis%nbf,basis%nbf,nspin))
 self_energy_old(:,:,:) = 0.0_dp
 if(calc_type%is_dft) allocate(hamiltonian_vxc(basis%nbf,basis%nbf,nspin))
 allocate(ehomo(nspin))
 allocate(elumo(nspin))


 call start_clock(timing_scf)

 !
 ! start the big scf loop
 !
 do iscf=1,nscf
   write(stdout,'(/,a)') '-------------------------------------------'
   write(stdout,'(a,x,i4,/)') ' *** SCF cycle No:',iscf

   en%kin  = SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
   en%nuc  = SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix(:,:,:),DIM=3) )

   !
   ! Setup kinetic and nucleus contributions (that are independent of the
   ! density matrix and therefore of spin channel)
   !
   hamiltonian(:,:,1) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:) 
   if(nspin==2) hamiltonian(:,:,nspin)    = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:) 

   if( calc_type%read_potential ) then
     call read_potential(print_matrix_,basis%nbf,nspin,p_matrix,matrix_tmp,en%hart)
     hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix_tmp(:,:,:)

   else

     !
     ! Hartree contribution to the Hamiltonian
     !
     if( .NOT. is_full_auxil) then
       call setup_hartree(print_matrix_,basis%nbf,nspin,p_matrix,hamiltonian_hartree,en%hart)
     else
       call setup_hartree_ri(print_matrix_,basis%nbf,nspin,p_matrix,hamiltonian_hartree,en%hart)
     endif
     do ispin=1,nspin
       hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) + hamiltonian_hartree(:,:)
     enddo
   endif

  
   !
   ! Reset XC part of the Hamiltonian
   hamiltonian_xc(:,:,:) = 0.0_dp

   !
   ! Exchange contribution to the Hamiltonian
   if( calc_type%need_exchange ) then

     if( .NOT. is_full_auxil) then
       call setup_exchange(print_matrix_,basis%nbf,p_matrix,hamiltonian_exx,en%exx)
     else
       call setup_exchange_ri(print_matrix_,basis%nbf,p_matrix,hamiltonian_exx,en%exx)
     endif
     ! Rescale with alpha_hybrid for hybrid functionals
     en%exx_hyb = alpha_hybrid * en%exx
     hamiltonian_xc(:,:,:) = hamiltonian_exx(:,:,:) * alpha_hybrid
   endif

   if(calc_type%need_exchange_lr) then
     if( .NOT. is_full_auxil) then
       call setup_exchange_longrange(print_matrix_,basis%nbf,p_matrix,matrix_tmp,energy_tmp)
     else
       call setup_exchange_longrange_ri(print_matrix_,basis%nbf,p_matrix,matrix_tmp,energy_tmp)
     endif
     ! Rescale with alpha_hybrid_lr for range-separated hybrid functionals
     en%exx_hyb = en%exx_hyb + alpha_hybrid_lr * energy_tmp
     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + matrix_tmp(:,:,:) * alpha_hybrid_lr
   endif


   !
   ! DFT XC potential is added here
   if( calc_type%is_dft ) then

     call dft_exc_vxc(basis,p_matrix,ehomo,hamiltonian_vxc,en%xc)

     title='=== DFT XC contribution ==='
     call dump_out_matrix(print_matrix_,title,basis%nbf,nspin,hamiltonian_vxc)

     hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + hamiltonian_vxc(:,:,:)
   endif

   !
   ! QPscGW self energy
   if( calc_type%is_gw .AND. ( calc_type%gwmethod == QS .OR. calc_type%gwmethod == QSCOHSEX) &
       .AND. iscf > 5 ) then

     if(has_auxil_basis) call prepare_eri_3center_eigen(c_matrix)
     call init_spectral_function(basis%nbf,nstate,occupation,wpol)
     call polarizability(basis,prod_basis,auxil_basis,occupation,energy,c_matrix,en%rpa,wpol)

     if( ABS(en%rpa) > 1.e-6_dp) then
       en%tot = en%tot + en%rpa
       write(stdout,'(/,a,f19.10)') ' RPA Total energy [Ha]: ',en%tot
     endif

     exchange_m_vxc_diag(:,:)=0.0_dp
     call gw_selfenergy(calc_type%gwmethod,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix_tmp,en%gw)
     if(has_auxil_basis) call destroy_eri_3center_eigen()

     matrix_tmp(:,:,:) = alpha_mixing * matrix_tmp(:,:,:) + (1.0_dp-alpha_mixing) * self_energy_old(:,:,:)
     self_energy_old(:,:,:) = matrix_tmp(:,:,:)
     title='=== Self-energy ==='
     call dump_out_matrix(print_matrix_,title,basis%nbf,nspin,matrix_tmp)
     call destroy_spectral_function(wpol)

     hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix_tmp(:,:,:)

   endif

   !
   ! QPscMP2
   if( calc_type%is_mp2 .AND. calc_type%gwmethod == QS .AND. iscf > 5 ) then

     exchange_m_vxc_diag(:,:)=0.0_dp

     call mp2_selfenergy(calc_type%gwmethod,basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,matrix_tmp,en%mp2)

     write(stdout,'(a,2x,f19.10)') ' MP2 Energy       [Ha]:',en%mp2
     write(stdout,*) 
     en%tot = en%tot + en%mp2
     write(stdout,'(a,2x,f19.10)') ' MP2 Total Energy [Ha]:',en%tot

     matrix_tmp(:,:,:) = alpha_mixing * matrix_tmp(:,:,:) + (1.0_dp-alpha_mixing) * self_energy_old(:,:,:)
     self_energy_old(:,:,:) = matrix_tmp(:,:,:)
     title='=== Self-energy ==='
     call dump_out_matrix(print_matrix_,title,basis%nbf,nspin,matrix_tmp)
  
     hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix_tmp(:,:,:)

   endif

   !
   ! Add the XC part of the hamiltonian to the total hamiltonian
   hamiltonian(:,:,:) = hamiltonian(:,:,:) + hamiltonian_xc(:,:,:)
   
   title='=== Total Hamiltonian ==='
   call dump_out_matrix(print_matrix_,title,basis%nbf,nspin,hamiltonian)

   !
   ! If requested, the level shifting procedure is triggered: 
   ! All the unoccupied states are penalized with an energy =  level_shifting_energy
   if( level_shifting_energy > 1.0e-6_dp ) then
     call level_shifting(basis%nbf,s_matrix,c_matrix,occupation,level_shifting_energy,hamiltonian)
   endif
  
  
   !
   ! Diagonalize the Hamiltonian H
   ! Generalized eigenvalue problem with overlap matrix S
   ! H \phi = E S \phi
   ! save the old eigenvalues
#if 0
   ! Keep old coding for reference
   ! TODO remove this when sure
   do ispin=1,nspin
     write(stdout,*) 'Diagonalization for spin channel',ispin
     call start_clock(timing_diago_hamiltonian)
     call diagonalize_generalized_sym(basis%nbf,&
                                      hamiltonian(:,:,ispin),s_matrix(:,:),&
                                      energy(:,ispin),c_matrix(:,:,ispin))
     call stop_clock(timing_diago_hamiltonian)
   enddo
#else 
   call diagonalize_hamiltonian(nspin,basis%nbf,nstate,hamiltonian,s_matrix_sqrt_inv,energy,c_matrix)
#endif

   !
   ! When level_shifting is used, the unoccupied state energies have to be brought
   ! back to their original value,
   ! So that the "physical" energies are written down
   if( level_shifting_energy > 1.0e-6_dp ) then
     do ispin=1,nspin
       do istate=1,basis%nbf
         if( occupation(istate,ispin) < completely_empty ) then
           energy(istate,ispin) = energy(istate,ispin) - level_shifting_energy
         endif
       enddo
     enddo
   endif
  
   title='=== Energies ==='
   call dump_out_energy(title,nstate,nspin,occupation(1:nstate,:),energy(1:nstate,:))

   call output_homolumo(basis%nbf,nspin,occupation,energy,ehomo,elumo)


   if(print_matrix_) then
     !
     ! REMEMBER:
     ! \phi_i = \sum_alpha C_{alpha i} \varphi_alpha 
     ! 
     ! hence transpose the c_matrix for a correct output by dump_out_matrix
     do ispin=1,nspin
       matrix_tmp(:,:,ispin) = TRANSPOSE( c_matrix(:,:,ispin) )
     enddo
     title='=== C coefficients ==='
     call dump_out_matrix(print_matrix_,title,basis%nbf,nspin,matrix_tmp)
     matrix_tmp(:,:,1) = MATMUL( c_matrix(:,:,1), MATMUL( s_matrix(:,:), TRANSPOSE(c_matrix(:,:,1)) ) )
     title='=== C S C^T = identity ? ==='
     call dump_out_matrix(print_matrix_,title,basis%nbf,1,matrix_tmp)
     matrix_tmp(:,:,1) = MATMUL( TRANSPOSE(c_matrix(:,:,1)), MATMUL( s_matrix(:,:), c_matrix(:,:,1) ) )
     title='=== C^T S C = identity ? ==='
     call dump_out_matrix(print_matrix_,title,basis%nbf,1,matrix_tmp)
   endif

  
   !
   ! Setup the new density matrix: p_matrix
   ! Save the old one for the convergence criterium
   p_matrix_old(:,:,:) = p_matrix(:,:,:)
   call setup_density_matrix(basis%nbf,nspin,c_matrix,occupation,p_matrix)
   title='=== density matrix P ==='
   call dump_out_matrix(print_matrix_,title,basis%nbf,nspin,p_matrix)
  
   !
   ! Output the total energy and its components
   write(stdout,*)
   write(stdout,'(a25,x,f19.10)') 'Nucleus-Nucleus [Ha]:',en%nuc_nuc
   write(stdout,'(a25,x,f19.10)') 'Kinetic Energy  [Ha]:',en%kin
   write(stdout,'(a25,x,f19.10)') 'Nucleus Energy  [Ha]:',en%nuc
   write(stdout,'(a25,x,f19.10)') 'Hartree Energy  [Ha]:',en%hart
   if(calc_type%need_exchange) then
     write(stdout,'(a25,x,f19.10)') 'Exchange Energy [Ha]:',en%exx_hyb
   endif
   if( calc_type%is_dft ) then
     write(stdout,'(a25,x,f19.10)') 'XC Energy       [Ha]:',en%xc
   endif
   en%tot = en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx_hyb + en%xc
   write(stdout,'(/,a25,x,f19.10,/)') 'Total Energy    [Ha]:',en%tot

   !
   ! Store the history of residuals
   call store_residual(p_matrix_old,p_matrix)
   !
   ! Produce the next density matrix
   call new_p_matrix(p_matrix)

   is_converged = check_converged()
   inquire(file='STOP',exist=stopfile_found)

   if( is_converged .OR. stopfile_found ) exit

   !
   ! Write down a "small" RESTART file at each step
   if( print_restart_ ) then
    call write_restart(SMALL_RESTART,basis,occupation,c_matrix,energy,hamiltonian_hartree,hamiltonian_exx,hamiltonian_xc)
   endif
   
 !
 ! end of the big SCF loop
 enddo


 write(stdout,*)
 write(stdout,*) '=================================================='
 write(stdout,*) 'The SCF loop ends here'
 write(stdout,*) '=================================================='

 call destroy_scf()

 !
 ! Spin contamination?
 call evaluate_s2_operator(nspin,basis%nbf,occupation,c_matrix,s_matrix)

 !
 ! Get the exchange operator if not already calculated
 !
 if( .NOT. is_full_auxil) then
   if( ABS(en%exx) < 1.0e-6_dp ) call setup_exchange(print_matrix_,basis%nbf,p_matrix,hamiltonian_exx,en%exx)
 else
   if( ABS(en%exx) < 1.0e-6_dp ) call setup_exchange_ri(print_matrix_,basis%nbf,p_matrix,hamiltonian_exx,en%exx)
 endif

 write(stdout,'(/,/,a25,x,f19.10,/)') 'SCF Total Energy [Ha]:',en%tot
 write(stdout,'(a25,x,f19.10)')       '      EXX Energy [Ha]:',en%exx
 write(stdout,'(a25,x,f19.10)')       'Total EXX Energy [Ha]:',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx

 !
 ! Single excitation term
 !
 ! Obtain the Fock matrix
 matrix_tmp(:,:,:) = hamiltonian(:,:,:) - hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:)
 ! And pass it to single_excitations
 call single_excitations(basis%nbf,energy,occupation,c_matrix,matrix_tmp)
 write(stdout,'(a25,x,f19.10)') 'Single Excitations [Ha]:',en%se
 write(stdout,'(a25,x,f19.10)')     'Est. HF Energy [Ha]:',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%se

 ! A dirty section for the Luttinger-Ward functional
#if 1
 if(calc_type%gwmethod==LW .OR. calc_type%gwmethod==LW2 .OR. calc_type%gwmethod==GSIGMA) then
   allocate(energy_exx(basis%nbf,nspin))
   allocate(c_matrix_exx(basis%nbf,basis%nbf,nspin))
   call issue_warning('ugly coding here write temp file fort.1000 and fort.1001')
!   call issue_warning('ugly coding here write temp file fort.2000 and fort.2001')
   do ispin=1,nspin
     write(stdout,*) 'Diagonalization H_exx for spin channel',ispin
     call diagonalize_generalized_sym(basis%nbf,&
                                      matrix_tmp(:,:,ispin),s_matrix(:,:),&
                                      energy_exx(:,ispin),c_matrix_exx(:,:,ispin))
   enddo
   write(stdout,*) 'FBFB LW sum(      epsilon) + Eii -<vxc> - EH + Ex',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx
   write(stdout,*) 'FBFB LW sum(tilde epsilon) + Eii - EH - Ex       ',SUM( occupation(:,:)*energy_exx(:,:) ) + en%nuc_nuc - en%hart - en%exx
   open(1000,form='unformatted')
   do ispin=1,nspin
     do istate=1,basis%nbf
       write(1000) c_matrix_exx(:,istate,ispin)
     enddo
   enddo
   close(1000)
!   open(2000,form='unformatted')
!   do ispin=1,nspin
!     do istate=1,basis%nbf
!       write(2000) c_matrix(:,istate,ispin)
!     enddo
!   enddo
!   close(2000)
   open(1001,form='unformatted')
   write(1001) energy_exx(:,:)
   close(1001)
!   open(2001,form='unformatted')
!   write(2001) energy(:,:)
!   close(2001)
!   call issue_warning('hacking the coefficients c_matrix_exx')
!   c_matrix=c_matrix_exx
!   energy=energy_exx
   deallocate(energy_exx,c_matrix_exx)
 endif
#endif


 !
 ! Testing the core/valence splitting
 !
 inquire(file='manual_coresplitting',exist=file_exists)
 if(file_exists) then
   if( alpha_hybrid_lr > 0.001 ) then
     call die('RSH not implemented yet')
   endif
   write(stdout,'(/,a)') ' TESTING CORE-VALENCE SPLITTING'
   open(newunit=fileunit,file='manual_coresplitting',status='old')
   read(fileunit,*) ncore
   close(fileunit)
   write(msg,'(a,i4,2x,i4)') 'core-valence splitting switched on up to state = ',ncore
   call issue_warning(msg)

   allocate(occupation_tmp(basis%nbf,nspin))
   allocate(p_matrix_tmp(basis%nbf,basis%nbf,nspin))
   ! Override the occupation of the core electrons
   occupation_tmp(:,:) = occupation(:,:)
   do istate=1,ncore
     occupation_tmp(istate,:) = 0.0_dp
   enddo
   call setup_density_matrix(basis%nbf,nspin,c_matrix,occupation_tmp,p_matrix_tmp)
   call dft_exc_vxc(basis,p_matrix_tmp,ehomo,hamiltonian_xc,en%xc)

   if( .NOT. is_full_auxil ) then
     call setup_exchange(print_matrix_,basis%nbf,p_matrix_tmp,hamiltonian_exx,en%exx)
   else
     call setup_exchange_ri(print_matrix_,basis%nbf,p_matrix_tmp,hamiltonian_exx,en%exx)
   endif

   deallocate(occupation_tmp,p_matrix_tmp)

   hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + alpha_hybrid * hamiltonian_exx(:,:,:)
 endif


 !
 ! Big RESTART file written if converged
 !
 if( is_converged .AND. print_restart_ ) then
   call write_restart(BIG_RESTART,basis,occupation,c_matrix,energy,hamiltonian_hartree,hamiltonian_exx,hamiltonian_xc)
 endif

 !
 ! Cleanly deallocate the integral grid information
 if( calc_type%is_dft ) call destroy_dft_grid()

 !
 ! Cleanly deallocate the arrays
 !
 deallocate(hamiltonian)
 deallocate(matrix_tmp,p_matrix_old)
 deallocate(exchange_m_vxc_diag)
 deallocate(self_energy_old)
 if( calc_type%is_dft ) deallocate(hamiltonian_vxc)

 call stop_clock(timing_scf)

end subroutine scf_loop


!=========================================================================
