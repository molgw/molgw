#include "macros.h"
!=========================================================================
subroutine scf_loop(basis,prod_basis,auxil_basis,&
                    s_matrix,c_matrix,p_matrix,&
                    hamiltonian_kinetic,hamiltonian_nucleus,&
                    hamiltonian_exx,hamiltonian_xc,&
                    occupation,energy)
 use m_definitions
 use m_mpi
 use m_timing
 use m_warning
 use m_inputparam
 use m_tools
 use m_scf
 use m_atoms
 use m_gaussian
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
 real(dp),intent(in)                :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(inout)             :: c_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(inout)             :: p_matrix(basis%nbf,basis%nbf,nspin)
 real(dp),intent(in)                :: hamiltonian_kinetic(basis%nbf,basis%nbf)
 real(dp),intent(in)                :: hamiltonian_nucleus(basis%nbf,basis%nbf)
 real(dp),intent(out)               :: hamiltonian_exx(basis%nbf,basis%nbf,nspin)
 real(dp),intent(out)               :: hamiltonian_xc(basis%nbf,basis%nbf,nspin)
 real(dp),intent(inout)             :: occupation(basis%nbf,nspin)
 real(dp),intent(inout)             :: energy(basis%nbf,nspin)
!=====
 type(spectral_function) :: wpol
 integer                 :: ispin,iscf
 character(len=100)      :: title
 real(dp)                :: energy_tmp
 real(dp),allocatable    :: ehomo(:),elumo(:)
 real(dp),allocatable    :: hamiltonian(:,:,:)
 real(dp),allocatable    :: hamiltonian_vxc(:,:,:)
 real(dp),allocatable    :: matrix_tmp(:,:,:)
 real(dp),allocatable    :: p_matrix_old(:,:,:)
 real(dp),allocatable    :: exchange_m_vxc_diag(:,:)
 real(dp),allocatable    :: self_energy_old(:,:,:)
 logical                 :: is_converged
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
   WRITE_MASTER(*,'(/,a)') '-------------------------------------------'
   WRITE_MASTER(*,'(a,x,i4,/)') ' *** SCF cycle No:',iscf

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
   else
     !
     ! Hartree contribution to the Hamiltonian
     !
     if( .NOT. is_full_auxil) then
       call setup_hartree(print_matrix_,basis%nbf,nspin,p_matrix,matrix_tmp,en%hart)
     else
       call setup_hartree_ri(print_matrix_,basis%nbf,nspin,p_matrix,matrix_tmp,en%hart)
     endif
   endif

   hamiltonian(:,:,:) = hamiltonian(:,:,:) + matrix_tmp(:,:,:)
  
   !
   ! Reset XC part of the Hamiltonian
   hamiltonian_xc(:,:,:) = 0.0_dp

   !
   ! Exchange contribution to the Hamiltonian
   if( calc_type%need_exchange ) then

     if( .NOT. is_full_auxil) then
       call setup_exchange(print_matrix_,basis%nbf,p_matrix,hamiltonian_exx,en%exx)
     else
       call setup_exchange_ri(print_matrix_,basis%nbf,c_matrix,occupation,p_matrix,hamiltonian_exx,en%exx)
     endif
     ! Rescale with alpha_hybrid for hybrid functionals
     en%exx = alpha_hybrid * en%exx
     hamiltonian_xc(:,:,:) = hamiltonian_exx(:,:,:) * alpha_hybrid

     if(calc_type%need_exchange_lr) then
       if( .NOT. is_full_auxil) then
         call setup_exchange_longrange(print_matrix_,basis%nbf,p_matrix,matrix_tmp,energy_tmp)
       else
         call setup_exchange_longrange_ri(print_matrix_,basis%nbf,c_matrix,occupation,p_matrix,matrix_tmp,energy_tmp)
       endif
       ! Rescale with alpha_hybrid_lr for range-separated hybrid functionals
       en%exx = en%exx + alpha_hybrid_lr * energy_tmp
       hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + matrix_tmp(:,:,:) * alpha_hybrid_lr
     endif

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
     call init_spectral_function(basis%nbf,occupation,wpol)
     call polarizability(basis,prod_basis,auxil_basis,occupation,energy,c_matrix,en%rpa,wpol)

     if( ABS(en%rpa) > 1.e-6_dp) then
       en%tot = en%tot + en%rpa
       WRITE_MASTER(*,'(/,a,f16.10)') ' RPA Total energy [Ha]: ',en%tot
     endif

     exchange_m_vxc_diag(:,:)=0.0_dp
     call gw_selfenergy(calc_type%gwmethod,basis,prod_basis,occupation,energy,exchange_m_vxc_diag,c_matrix,s_matrix,wpol,matrix_tmp)
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

     WRITE_MASTER(*,'(a,2x,f16.10)') ' MP2 Energy       [Ha]:',en%mp2
     WRITE_MASTER(*,*) 
     en%tot = en%tot + en%mp2
     WRITE_MASTER(*,'(a,2x,f16.10)') ' MP2 Total Energy [Ha]:',en%tot

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
   ! Diagonalize the Hamiltonian H
   ! Generalized eigenvalue problem with overlap matrix S
   ! H \phi = E S \phi
   ! save the old eigenvalues
   do ispin=1,nspin
     WRITE_MASTER(*,*) 'Diagonalization for spin channel',ispin
     call start_clock(timing_diago_hamiltonian)
     call diagonalize_generalized_sym(basis%nbf,&
                                      hamiltonian(:,:,ispin),s_matrix(:,:),&
                                      energy(:,ispin),c_matrix(:,:,ispin))
     call stop_clock(timing_diago_hamiltonian)
   enddo
  
   title='=== Energies ==='
   call dump_out_eigenenergy(title,basis%nbf,nspin,occupation,energy)

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
   WRITE_MASTER(*,*)
   WRITE_MASTER(*,'(a25,x,f16.10)') 'Nucleus-Nucleus [Ha]:',en%nuc_nuc
   WRITE_MASTER(*,'(a25,x,f16.10)') 'Kinetic Energy  [Ha]:',en%kin
   WRITE_MASTER(*,'(a25,x,f16.10)') 'Nucleus Energy  [Ha]:',en%nuc
   WRITE_MASTER(*,'(a25,x,f16.10)') 'Hartree Energy  [Ha]:',en%hart
   if(calc_type%need_exchange) then
     WRITE_MASTER(*,'(a25,x,f16.10)') 'Exchange Energy [Ha]:',en%exx
   endif
   if( calc_type%is_dft ) then
     WRITE_MASTER(*,'(a25,x,f16.10)') 'XC Energy       [Ha]:',en%xc
   endif
   en%tot = en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%xc
   WRITE_MASTER(*,'(/,a25,x,f16.10,/)') 'Total Energy    [Ha]:',en%tot

   !
   ! Store the history of residuals
   call store_residual(p_matrix_old,p_matrix)
   !
   ! Produce the next density matrix
   call new_p_matrix(p_matrix)

   is_converged = check_converged()
   if( is_converged ) exit

   !
   ! Write down a "small" RESTART file at each step
   ! TODO: I should check if it is resource consuming
   call write_small_restart(basis%nbf,occupation,c_matrix)
   
 !
 ! end of the big SCF loop
 enddo


 WRITE_MASTER(*,*)
 WRITE_MASTER(*,*) '=================================================='
 WRITE_MASTER(*,*) 'The SCF loop ends here'
 WRITE_MASTER(*,*) '=================================================='

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
   if( ABS(en%exx) < 1.0e-6_dp ) call setup_exchange_ri(print_matrix_,basis%nbf,c_matrix,occupation,p_matrix,hamiltonian_exx,en%exx)
 endif

 WRITE_MASTER(*,'(/,/,a25,x,f16.10,/)') 'SCF Total Energy [Ha]:',en%tot
 WRITE_MASTER(*,'(a25,x,f16.10)')       '      EXX Energy [Ha]:',en%exx
 WRITE_MASTER(*,'(a25,x,f16.10)')       'Total EXX Energy [Ha]:',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx

 !
 ! Single excitation term
 !
 ! Obtain the Fock matrix
 matrix_tmp(:,:,:) = hamiltonian(:,:,:) - hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:)
 ! And pass it to single_excitations
 call single_excitations(basis%nbf,energy,occupation,c_matrix,matrix_tmp)

 WRITE_MASTER(*,'(a25,x,f16.10)') 'Single Excitations [Ha]:',en%se
 WRITE_MASTER(*,'(a25,x,f16.10)')     'Est. HF Energy [Ha]:',en%nuc_nuc + en%kin + en%nuc + en%hart + en%exx + en%se

 !
 ! Big RESTART file written if converged
 !
 if( is_converged ) call write_big_restart(basis%nbf,occupation,c_matrix,energy,hamiltonian_exx,hamiltonian_xc)
 if( calc_type%is_dft ) then
!   ! Output the self-consistent density on the real-space grid
!   call write_density_grid(basis,p_matrix)
   call destroy_dft_grid()
 endif

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
