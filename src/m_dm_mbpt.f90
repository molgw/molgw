!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the reading or the calculation of correlated density matrix
!
!=========================================================================
module m_dm_mbpt
 use m_definitions
 use m_timing
 use m_warning
 use m_memory
 use m_atoms
 use m_basis_set
 use m_inputparam
 use m_spectral_function
 use m_selfenergy_tools
 use m_hamiltonian_tools
 use m_hamiltonian_wrapper
 use m_scf
 use m_multipole


contains


!=========================================================================
subroutine get_dm_mbpt(basis,occupation,energy,c_matrix, &
                       hamiltonian_kinetic,hamiltonian_nucleus,hamiltonian_fock)
 implicit none

 type(basis_set),intent(in)      :: basis
 real(dp),intent(in)             :: occupation(:,:)
 real(dp),intent(in)             :: energy(:,:)
 real(dp),intent(in)             :: c_matrix(:,:,:)
 real(dp),intent(in)             :: hamiltonian_kinetic(:,:)
 real(dp),intent(in)             :: hamiltonian_nucleus(:,:)
 real(dp),intent(inout)          :: hamiltonian_fock(:,:,:)
!=====
 integer                    :: nstate,nocc
 logical                    :: density_matrix_found
 integer                    :: file_density_matrix
 integer                    :: ispin,istate
 type(spectral_function)    :: wpol
 type(energy_contributions) :: en_dm_corr
 real(dp),allocatable       :: h_ii(:,:),exchange_ii(:,:)
 real(dp),allocatable       :: p_matrix_corr(:,:,:)
 real(dp),allocatable       :: hamiltonian_hartree_corr(:,:)
 real(dp),allocatable       :: hamiltonian_exx_corr(:,:,:)
 real(dp),allocatable       :: c_matrix_tmp(:,:,:)
 real(dp),allocatable       :: occupation_tmp(:,:)
!=====

 nstate = SIZE(c_matrix,DIM=2)

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
 if( TRIM(pt_density_matrix) /= 'NO' ) then
   call selfenergy_set_state_range(nstate,occupation)
   call fock_density_matrix(basis,occupation,energy,c_matrix,hamiltonian_fock,p_matrix_corr)

   select case(TRIM(pt_density_matrix))
   case('ONE-RING')
     ! This keyword calculates the 1-ring density matrix as it is derived in PT2 theory
     call onering_density_matrix(nstate,basis,occupation,energy,c_matrix,p_matrix_corr)
   case('PT2')
     ! This keyword calculates the PT2 density matrix as it is derived in PT2 theory (differs from MP2 density matrix)
     call pt2_density_matrix(nstate,basis,occupation,energy,c_matrix,p_matrix_corr)
   case('GW','G0W0')
     ! This keyword calculates the GW density matrix as it is derived in the new GW theory
     call init_spectral_function(nstate,occupation,0,wpol)
     call polarizability(.TRUE.,.TRUE.,basis,nstate,occupation,energy,c_matrix,en_dm_corr%rpa,en_dm_corr%gw,wpol)
     call gw_density_matrix(nstate,basis,occupation,energy,c_matrix,wpol,p_matrix_corr)
     call destroy_spectral_function(wpol)
   case('GW_IMAGINARY','G0W0_IMAGINARY')
     ! This keyword calculates the GW density matrix as it is derived in the new GW theory
     ! using an imaginary axis integral
     call init_spectral_function(nstate,occupation,nomega_imag,wpol)
     call polarizability_grid_scalapack(basis,nstate,occupation,energy,c_matrix,en_dm_corr%rpa,wpol)
     call gw_density_matrix_imag(nstate,basis,occupation,energy,c_matrix,wpol,p_matrix_corr)
     call destroy_spectral_function(wpol)
   case('GW_DYSON','G0W0_DYSON')
     ! This keyword calculates the GW density matrix as it is derived in the new GW theory
     ! using an imaginary axis integral
     call init_spectral_function(nstate,occupation,nomega_imag,wpol)
     call polarizability_grid_scalapack(basis,nstate,occupation,energy,c_matrix,en_dm_corr%rpa,wpol)
     call gw_density_matrix_dyson_imag(nstate,basis,occupation,energy,c_matrix,wpol,p_matrix_corr)
     call destroy_spectral_function(wpol)
   end select
 endif


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
     call die('get_dm_mbpt: no correlated density matrix read or calculated though input file suggests you really want one')
   endif

 endif

 if( print_hartree_ .OR. use_correlated_density_matrix_ ) then

   !
   ! Nucleus-nucleus repulsion contribution to the energy
   call nucleus_nucleus_energy(en_dm_corr%nuc_nuc)
   en_dm_corr%kinetic = SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix_corr(:,:,:),DIM=3) )
   en_dm_corr%nucleus = SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix_corr(:,:,:),DIM=3) )

   call calculate_hartree(basis,p_matrix_corr,hamiltonian_hartree_corr,eh=en_dm_corr%hartree)

   call calculate_exchange(basis,p_matrix_corr,hamiltonian_exx_corr,ex=en_dm_corr%exx)

   en_dm_corr%total = en_dm_corr%nuc_nuc + en_dm_corr%kinetic + en_dm_corr%nucleus +  en_dm_corr%hartree + en_dm_corr%exx
   write(stdout,'(/,1x,a)') 'Energies from correlated density matrix'
   write(stdout,'(a35,1x,f19.10)')   'Kinetic Energy (Ha):',en_dm_corr%kinetic
   write(stdout,'(a35,1x,f19.10)')   'Nucleus Energy (Ha):',en_dm_corr%nucleus
   write(stdout,'(a35,1x,f19.10)')   'Hartree Energy (Ha):',en_dm_corr%hartree
   write(stdout,'(a35,1x,f19.10)')  'Exchange Energy (Ha):',en_dm_corr%exx
   write(stdout,'(a35,1x,f19.10)') 'Total EXX Energy (Ha):',en_dm_corr%total

   if( ABS(en_dm_corr%gw) > 1.0e-8_dp ) then
     write(stdout,'(a35,1x,f19.10)')  'GW correlation Energy (Ha):',en_dm_corr%gw
     en_dm_corr%total = en_dm_corr%total + en_dm_corr%gw
     write(stdout,'(a35,1x,f19.10)')  'Total GM Energy (Ha):',en_dm_corr%total
   endif

   nocc = get_number_occupied_states(occupation)
   allocate(h_ii(nstate,nspin))

   call matrix_ao_to_mo_diag(c_matrix,RESHAPE(hamiltonian_hartree_corr,(/basis%nbf,basis%nbf,1/)),h_ii)
   call dump_out_energy('=== Hartree expectation value from correlated density matrix ===',occupation,h_ii)
   write(stdout,'(1x,a,2(3x,f12.6))') 'Hartree  HOMO expectation (eV):',h_ii(nocc,:) * Ha_eV

   call matrix_ao_to_mo_diag(c_matrix,hamiltonian_exx_corr,h_ii)
   call dump_out_energy('=== Exchange expectation value from correlated density matrix ===',occupation,h_ii)
   write(stdout,'(1x,a,2(3x,f12.6))') 'Exchange HOMO expectation (eV):',h_ii(nocc,:) * Ha_eV
   deallocate(h_ii)

 endif

 if( print_multipole_ .OR. print_cube_ ) then
   allocate(c_matrix_tmp(basis%nbf,basis%nbf,nspin))
   allocate(occupation_tmp(basis%nbf,nspin))
   call get_c_matrix_from_p_matrix(p_matrix_corr,c_matrix_tmp,occupation_tmp)
   if( print_multipole_ ) then
     call static_dipole(basis,occupation_tmp,c_matrix_tmp)
     call static_quadrupole(basis,occupation_tmp,c_matrix_tmp)
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
     hamiltonian_fock(:,:,ispin) = hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:) + hamiltonian_hartree_corr(:,:)  &
                                  + hamiltonian_exx_corr(:,:,ispin)
   enddo

 endif

 write(stdout,*)
 call clean_deallocate('Correlated density matrix',p_matrix_corr)
 call clean_deallocate('Correlated Hartree potential',hamiltonian_hartree_corr)
 call clean_deallocate('Correlated exchange operator',hamiltonian_exx_corr)


end subroutine get_dm_mbpt


!=========================================================================
subroutine fock_density_matrix(basis,occupation,energy,c_matrix,hfock,p_matrix)
 implicit none

 type(basis_set),intent(in)         :: basis
 real(dp),intent(in)                :: occupation(:,:),energy(:,:)
 real(dp),intent(in)                :: c_matrix(:,:,:)
 real(dp),intent(in)                :: hfock(:,:,:)
 real(dp),intent(out)               :: p_matrix(:,:,:)
!=====
 integer  :: nstate
 integer  :: pstate,qstate
 integer  :: istate,jstate
 integer  :: astate,bstate
 integer  :: pqspin
 real(dp),allocatable :: p_matrix_state(:,:,:)
 real(dp),allocatable :: hfock_state(:,:,:)
!=====

 call start_clock(timing_mbpt_dm)
 write(stdout,'(/,1x,a)') 'Calculate the perturbative Fock density matrix'

 nstate = SIZE(occupation,DIM=1)

 call clean_allocate('Density matrix in state basis',p_matrix_state,nstate,nstate,nspin)
 call clean_allocate('Fock matrix in state basis',hfock_state,nstate,nstate,nspin)

 call matrix_ao_to_mo(c_matrix,hfock,hfock_state)

 p_matrix_state(:,:,:) = 0.0_dp
 do pqspin=1,nspin
   do pstate=1,nstate
     p_matrix_state(pstate,pstate,pqspin) = occupation(pstate,pqspin)
   enddo
   do istate=1,nhomo_G
     do astate=nhomo_G+1,nstate
       p_matrix_state(istate,astate,pqspin) = hfock_state(istate,astate,pqspin)  &
                                              / ( energy(istate,pqspin) - energy(astate,pqspin) ) * spin_fact
       p_matrix_state(astate,istate,pqspin) = p_matrix_state(istate,astate,pqspin)
     enddo
   enddo
 enddo

 call matrix_mo_to_ao(c_matrix,p_matrix_state,p_matrix)

 call clean_deallocate('Density matrix in state basis',p_matrix_state)
 call clean_deallocate('Fock matrix in state basis',hfock_state)

 call stop_clock(timing_mbpt_dm)

end subroutine fock_density_matrix



!=========================================================================
end module m_dm_mbpt
!=========================================================================
