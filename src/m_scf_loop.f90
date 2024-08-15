!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval, Mauricio Rodriguez-Mayorga
!
! This module contains
! the main SCF loop for Hartree-Fock or generalized Kohn-Sham
!
!=========================================================================
#include "molgw.h"
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
  use m_selfenergy_tools
  use m_linear_response
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
                    c_matrix,en_gks,scf_has_converged)
  implicit none

  !=====
  logical,intent(in)                 :: is_restart
  type(basis_set),intent(inout)      :: basis
  real(dp),intent(in)                :: x_matrix(:,:)
  real(dp),intent(in)                :: s_matrix(:,:)
  real(dp),intent(in)                :: hamiltonian_kinetic(:,:)
  real(dp),intent(in)                :: hamiltonian_nucleus(:,:)
  real(dp),intent(inout)             :: occupation(:,:)
  real(dp),intent(out)               :: energy(:,:)
  real(dp),allocatable,intent(inout) :: hamiltonian_fock(:,:,:)
  real(dp),allocatable,intent(inout) :: c_matrix(:,:,:)
  type(energy_contributions),intent(inout) :: en_gks
  logical,intent(out)                :: scf_has_converged
  !=====
  type(spectral_function) :: wpol
  integer                 :: nstate
  logical                 :: stopfile_found
  integer                 :: file_density_matrix
  integer                 :: ispin,iscf
  integer                 :: ibf,idir
  real(dp)                :: trace(3)
  real(dp),allocatable    :: tmp_matrix(:,:)
  real(dp),allocatable    :: nabla_ao(:,:,:)
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
    write(stdout,'(/,1x,a)') '-------------------------------------------'
    write(stdout,'(a,1x,i4,/)') ' *** SCF cycle No:',iscf


    en_gks%kinetic  = SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix(:,:,:),DIM=3) )
    en_gks%nucleus  = SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix(:,:,:),DIM=3) )

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
    call calculate_hartree(basis,p_matrix,hamiltonian_hartree,eh=en_gks%hartree)

    ! calc_type%is_core is an inefficient way to get the Kinetic+Nucleus Hamiltonian
    if( calc_type%is_core ) then
      hamiltonian_hartree(:,:) = 0.0_dp
      en_gks%hartree = 0.0_dp
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

      ! Rescale with beta_hybrid for range-separated hybrid functionals
      en_gks%exx_hyb = en_gks%exx_hyb * beta_hybrid
      hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + hamiltonian_exx(:,:,:) * beta_hybrid

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


      call wpol%init(nstate,occupation,0)
      call polarizability(.TRUE.,.TRUE.,basis,occupation,energy,c_matrix,en_gks%rpa,en_gks%gw,wpol)

      if( ABS(en_gks%rpa) > 1.e-6_dp) then
        en_gks%total = en_gks%total + en_gks%rpa
        write(stdout,'(/,a,f19.10)') ' RPA Total energy (Ha): ',en_gks%total
      endif

      !
      ! Set the range of states on which to evaluate the self-energy
      call selfenergy_set_state_range(nstate,occupation,range='all')

      allocate(matrix_tmp(basis%nbf,basis%nbf,nspin))
      call gw_selfenergy_qs(nstate,basis,occupation,energy,c_matrix,s_matrix,wpol,matrix_tmp)

      call dump_out_matrix(.FALSE.,'=== Self-energy ===',matrix_tmp)
      call wpol%destroy()

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
      en_gks%total = en_gks%total + en_gks%mp2
      write(stdout,'(a,2x,f19.10)') ' MP2 Total Energy (Ha):',en_gks%total

      call dump_out_matrix(.FALSE.,'=== Self-energy ===',matrix_tmp)

      hamiltonian_xc(:,:,:) = hamiltonian_xc(:,:,:) + matrix_tmp(:,:,:)
      deallocate(matrix_tmp)

    endif

    !
    ! Add the XC part of the hamiltonian to the total hamiltonian
    hamiltonian(:,:,:) = hamiltonian(:,:,:) + hamiltonian_xc(:,:,:)


    ! All the components of the energy have been calculated at this stage
    ! Sum up to get the total energy
    en_gks%total = en_gks%nuc_nuc + en_gks%kinetic + en_gks%nucleus + en_gks%hartree + en_gks%exx_hyb + en_gks%xc

    ! Make sure all the MPI tasks have the exact same Hamiltonian
    ! It helps stabilizing the SCF cycles in parallel
    call world%sum(hamiltonian)
    hamiltonian(:,:,:) = hamiltonian(:,:,:) / REAL(world%nproc,dp)

    !
    ! If requested, the level shifting procedure is triggered:
    ! All the unoccupied states are penalized with an energy =  level_shifting_energy
    if( level_shifting_energy > 1.0e-6_dp ) then
      call level_shifting_up(s_matrix,c_matrix,occupation,level_shifting_energy,hamiltonian)
    endif

    ! DIIS or simple mixing on the hamiltonian
    call hamiltonian_prediction(s_matrix,x_matrix,p_matrix,hamiltonian,en_gks%total)


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

    call dump_out_energy('=== Energies ===',occupation,energy)

    call output_homolumo('gKS',occupation,energy,1,nstate)


    !
    ! Output the total energy and its components
    write(stdout,*)
    write(stdout,'(a25,1x,f19.10)') 'Nucleus-Nucleus (Ha):',en_gks%nuc_nuc
    write(stdout,'(a25,1x,f19.10)') 'Kinetic Energy  (Ha):',en_gks%kinetic
    write(stdout,'(a25,1x,f19.10)') 'Nucleus Energy  (Ha):',en_gks%nucleus
    write(stdout,'(a25,1x,f19.10)') 'Hartree Energy  (Ha):',en_gks%hartree
    if(calc_type%need_exchange) then
      write(stdout,'(a25,1x,f19.10)') 'Exchange Energy (Ha):',en_gks%exx_hyb
    endif
    if( calc_type%is_dft ) then
      write(stdout,'(a25,1x,f19.10)') 'XC Energy       (Ha):',en_gks%xc
    endif
    write(stdout,'(/,a25,1x,f19.10,/)') 'Total Energy    (Ha):',en_gks%total


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


    scf_has_converged = check_converged(p_matrix)
    inquire(file='STOP',exist=stopfile_found)

    if( scf_has_converged .OR. stopfile_found ) exit

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
    call print_expectations(basis,c_matrix,hamiltonian_kinetic)
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
  ! Form the final Fock matrix and store it only if needed
  !
  if( scf_has_converged  &
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


  write(stdout,'(/,/,a25,1x,f19.10,/)') 'SCF Total Energy (Ha):',en_gks%total

  if( ABS(en_gks%exx) > 1.0e-10) then
    write(stdout,'(a25,1x,f19.10)')       '      EXX Energy (Ha):',en_gks%exx
    en_gks%totalexx = en_gks%nuc_nuc + en_gks%kinetic + en_gks%nucleus + en_gks%hartree + en_gks%exx
    write(stdout,'(a25,1x,f19.10)')       'Total EXX Energy (Ha):',en_gks%totalexx
  endif

  if( print_yaml_ .AND. is_iomaster ) then
    if( scf_has_converged ) then
      write(unit_yaml,'(/,a)') 'scf is converged: True'
    else
      write(unit_yaml,'(/,a)') 'scf is converged: False'
    endif
    call print_energy_yaml('scf energy',en_gks)
    call dump_out_energy_yaml('gks energies',energy)
  endif

  call stop_clock(timing_scf)

end subroutine scf_loop

!=========================================================================
subroutine scf_loop_cmplx(is_restart,&
                    basis,&
                    x_matrix,s_matrix,&
                    hamiltonian_kinetic,hamiltonian_nucleus,&
                    occupation, &
                    energy, &
                    c_matrix,c_matrix_cmplx,en_gks,scf_has_converged)
  implicit none

  !=====
  logical,intent(in)                 :: is_restart
  type(basis_set),intent(inout)      :: basis
  real(dp),intent(in)                :: x_matrix(:,:)
  real(dp),intent(in)                :: s_matrix(:,:)
  real(dp),intent(in)                :: hamiltonian_kinetic(:,:)
  real(dp),intent(in)                :: hamiltonian_nucleus(:,:)
  real(dp),intent(inout)             :: occupation(:,:)
  real(dp),intent(out)               :: energy(:,:)
  real(dp),allocatable,intent(inout) :: c_matrix(:,:,:)
  complex(dp),allocatable,intent(inout) :: c_matrix_cmplx(:,:,:)
  type(energy_contributions),intent(inout) :: en_gks
  logical,intent(out)                :: scf_has_converged
  !=====
  integer                 :: nstate
  logical                 :: stopfile_found
  integer                 :: ispin,iscf,jbf,istate
  real(dp)                :: rms
  real(dp),allocatable    :: s_eigval(:)
  real(dp),allocatable    :: inv_x_matrix(:,:),matrix_tmp(:,:)
  real(dp),allocatable    :: p_matrix_real(:,:)
  complex(dp),allocatable :: hsmall_cmplx(:,:), csmall_cmplx(:,:)
  complex(dp),allocatable :: hamiltonian_cmplx(:,:,:)
  complex(dp),allocatable :: p_matrix_cmplx(:,:,:),p_matrix_cmplx_old(:,:,:)
  complex(dp),allocatable :: ham_hist(:,:,:,:)
  !=====


  call start_clock(timing_scf)

  rms = 1000.0
  nstate = SIZE(x_matrix,DIM=2)
  
  !
  ! Allocate the main arrays
  call clean_allocate('Total Hamiltonian H',hamiltonian_cmplx,basis%nbf,basis%nbf,nspin)
  call clean_allocate('H in orthogonalized basis',hsmall_cmplx,nstate,nstate)
  call clean_allocate('H eigenvectors',csmall_cmplx,nstate,nstate)
  call clean_allocate('Density matrix P',p_matrix_cmplx,basis%nbf,basis%nbf,nspin)
  call clean_allocate('Density matrix P(old)',p_matrix_cmplx_old,basis%nbf,basis%nbf,nspin)
  call clean_allocate('Hamiltonian history',ham_hist,basis%nbf,basis%nbf,nspin,2)
  ham_hist=COMPLEX_ZERO 

  !
  ! Setup the grids for the quadrature of DFT potential/energy
  if( calc_type%is_dft ) then
    call init_dft_grid(basis,grid_level,dft_xc(1)%needs_gradient,.TRUE.,BATCH_SIZE)
  endif

  !
  ! Setup the density matrix: p_matrix
  call setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)


  !
  ! Start the big scf loop
  !
  do iscf=1,nscf
    write(stdout,'(/,1x,a)') '-------------------------------------------'
    write(stdout,'(a,1x,i4,/)') ' *** SCF cycle No:',iscf


    en_gks%kinetic = REAL( SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
    en_gks%nucleus = REAL( SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
    hamiltonian_cmplx=COMPLEX_ZERO

    !--Hamiltonian - Hartree Exchange Correlation---
    call calculate_hamiltonian_hxc_ri_cmplx(basis,                    &
                                            occupation,               &
                                            c_matrix_cmplx,           &
                                            p_matrix_cmplx,           &
                                            hamiltonian_cmplx,en_gks)

    !
    ! Setup kinetic and nucleus contributions (that are independent of the
    ! density matrix and therefore of spin channel)
    !
    do ispin=1,nspin
      hamiltonian_cmplx(:,:,ispin) = hamiltonian_cmplx(:,:,ispin)     &
                      + hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
    enddo

    !! Sum up to get the total energy
    en_gks%total = en_gks%nuc_nuc + en_gks%kinetic + en_gks%nucleus + en_gks%hartree + en_gks%exx_hyb + en_gks%xc

    ! Make sure all the MPI tasks have the exact same Hamiltonian
    ! It helps stabilizing the SCF cycles in parallel
    call world%sum(hamiltonian_cmplx)
    hamiltonian_cmplx(:,:,:) = hamiltonian_cmplx(:,:,:) / REAL(world%nproc,dp)

    !
    ! Simple mixing on the hamiltonian
    ! the newest is 1
    ! the oldest is 2
    if( iscf > 1) then
      ham_hist(:,:,:,2)=ham_hist(:,:,:,1)
    endif
    ham_hist(:,:,:,1)=hamiltonian_cmplx(:,:,:)
    if( iscf > 1) then
      write(stdout,'(/,1x,a,f8.4)') 'Simple mixing with parameter:',alpha_mixing
      hamiltonian_cmplx(:,:,:) = alpha_mixing * ham_hist(:,:,:,1) + (1.0_dp - alpha_mixing) * ham_hist(:,:,:,2)
    endif

    !
    ! H \varphi = S \varphi E
    ! Diagonalize the Hamiltonian S^-1/2 H S^-1/2
    do ispin=1,nspin

      !
      ! H' = X**T * H * X
      ! TODO: slow coding. Use BLAS-level 3 in the future
      hsmall_cmplx(:,:) = MATMUL(TRANSPOSE(x_matrix), MATMUL(hamiltonian_cmplx(:,:,ispin), x_matrix))

      !
      ! H' * C' = C' * E
      call diagonalize(' ',hsmall_cmplx,energy(:,ispin),csmall_cmplx)

      !
      ! C = X * C'
      c_matrix_cmplx(:,:,ispin) = MATMUL(x_matrix, csmall_cmplx)
    enddo

    call dump_out_energy('=== Energies ===',occupation,energy)

    call output_homolumo('gKS',occupation,energy,1,nstate)


    !
    ! Output the total energy and its components
    write(stdout,*)
    write(stdout,'(a25,1x,f19.10)') 'Nucleus-Nucleus (Ha):',en_gks%nuc_nuc
    write(stdout,'(a25,1x,f19.10)') 'Kinetic Energy  (Ha):',en_gks%kinetic
    write(stdout,'(a25,1x,f19.10)') 'Nucleus Energy  (Ha):',en_gks%nucleus
    write(stdout,'(a25,1x,f19.10)') 'Hartree Energy  (Ha):',en_gks%hartree
    if(calc_type%need_exchange) then
      write(stdout,'(a25,1x,f19.10)') 'Exchange Energy (Ha):',en_gks%exx_hyb
    endif
    if( calc_type%is_dft ) then
      write(stdout,'(a25,1x,f19.10)') 'XC Energy       (Ha):',en_gks%xc
    endif
    write(stdout,'(/,a25,1x,f19.10,/)') 'Total Energy    (Ha):',en_gks%total


    !
    ! Setup the new density matrix: p_matrix
    ! Save the old one for the convergence criterium
    call setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)


 ! SCF convergence check
    if( iscf > 1) then
      rms = NORM2( real(p_matrix_cmplx(:,:,:)) - real(p_matrix_cmplx_old(:,:,:)) ) * SQRT( REAL(nspin,dp) ) &
          + NORM2( aimag(p_matrix_cmplx(:,:,:)) - aimag(p_matrix_cmplx_old(:,:,:)) ) * SQRT( REAL(nspin,dp) )
      p_matrix_cmplx_old(:,:,:)=p_matrix_cmplx(:,:,:)
      write(stdout,'(1x,a,es12.5)') 'Convergence criterium on the density matrix: ',rms
    else
      p_matrix_cmplx_old(:,:,:)=p_matrix_cmplx(:,:,:)
    endif

    if( rms < tolscf ) then
      scf_has_converged = .TRUE.
      write(stdout,*) ' ===> convergence has been reached'
      write(stdout,*)
    else
      scf_has_converged = .FALSE.
      write(stdout,*) ' ===> convergence not reached yet'
      write(stdout,*)

      if( iscf == nscf ) then
        if( rms > 1.0e-2_dp ) then
          call issue_warning('SCF convergence is very poor')
        else if( rms > 1.0e-4_dp ) then
          call issue_warning('SCF convergence is poor')
        endif
      endif

    endif

    inquire(file='STOP',exist=stopfile_found)

    ! This typically stops too early. We have to use tolscf<0 !
    if( scf_has_converged .OR. stopfile_found ) exit

    !
    ! end of the big SCF loop
  enddo


  write(stdout,'(/,1x,a)') '=================================================='
  write(stdout,'(1x,a)') 'The SCF loop ends here'
  write(stdout,'(1x,a)') '=================================================='


  !
  ! Store the natural orbital basis representation of the density matrix in c_matrix
  ! and the occupation numbers in occupation(:,1) \in [0:2]
  call clean_allocate('Density matrix P real',p_matrix_real,basis%nbf,basis%nbf)
  if(nspin==2) then
   p_matrix_real(:,:)=real(p_matrix_cmplx(:,:,1)+p_matrix_cmplx(:,:,2))
  else
   p_matrix_real(:,:)=real(p_matrix_cmplx(:,:,1))
  endif
  allocate(s_eigval(basis%nbf),matrix_tmp(basis%nbf,basis%nbf))
  matrix_tmp(:,:) = s_matrix(:,:)
  ! Diagonalization with or without SCALAPACK
  !! S = U*s*U^H
  call diagonalize_scalapack(scf_diago_flavor,scalapack_block_min,matrix_tmp,s_eigval)
  nstate = COUNT( s_eigval(:) > min_overlap )
  call clean_allocate('Overlap INV_X * INV_X**H = S',inv_x_matrix,basis%nbf,nstate)
  inv_x_matrix=0.0_dp
  write(stdout,'(/,a)')       ' Filtering basis functions that induce overcompleteness'
  write(stdout,'(a,es9.2)')   '   Lowest S eigenvalue is           ',MINVAL( s_eigval(:) )
  write(stdout,'(a,es9.2)')   '   Tolerance on overlap eigenvalues ',min_overlap
  write(stdout,'(a,i5,a,i5)') '   Retaining ',nstate,' among ',basis%nbf
  !! INV_X = U*s^(1/2)
  istate = 0
  do jbf=1,basis%nbf
    if( s_eigval(jbf) > min_overlap ) then
      istate = istate + 1
      inv_x_matrix(:,istate) = matrix_tmp(:,jbf) * SQRT( s_eigval(jbf) )
    endif
  enddo
  !  diag[ S^1/2 P S^1/2 ] -> V
  deallocate(matrix_tmp)
  allocate(matrix_tmp(nstate,nstate))
  matrix_tmp=0.0_dp
  matrix_tmp=MATMUL(TRANSPOSE(inv_x_matrix), MATMUL(p_matrix_real, inv_x_matrix))
  call diagonalize(' ',matrix_tmp,occupation(:,1),c_matrix(:,:,1))
  !  C = S^-1/2 V
  c_matrix(:,:,1) = MATMUL(x_matrix, c_matrix(:,:,1))
  deallocate(matrix_tmp,s_eigval)
  if(nspin==2) then
   occupation(:,2)=0.0_dp
   c_matrix(:,:,2)=0.0_dp
  endif

  !
  ! Cleanly deallocate the integral grid information
  !
  if( calc_type%is_dft ) call destroy_dft_grid()


  !
  ! Cleanly deallocate the arrays
  !
  call clean_deallocate('Overlap INV_X * INV_X**H = S',inv_x_matrix)
  call clean_deallocate('Density matrix P real',p_matrix_real)
  call clean_deallocate('Density matrix P',p_matrix_cmplx)
  call clean_deallocate('Density matrix P(old)',p_matrix_cmplx_old)
  call clean_deallocate('Total Hamiltonian H',hamiltonian_cmplx)
  call clean_deallocate('H in orthogonalized basis',hsmall_cmplx)
  call clean_deallocate('H eigenvectors',csmall_cmplx)
  call clean_deallocate('Hamiltonian history',ham_hist)

  write(stdout,'(/,/,a25,1x,f19.10,/)') 'SCF Total Energy (Ha):',en_gks%total

  if( print_yaml_ .AND. is_iomaster ) then
    if( scf_has_converged ) then
      write(unit_yaml,'(/,a)') 'scf is converged: True'
    else
      write(unit_yaml,'(/,a)') 'scf is converged: False'
    endif
    call print_energy_yaml('scf energy',en_gks)
    call dump_out_energy_yaml('gks energies',energy)
  endif

  call stop_clock(timing_scf)

end subroutine scf_loop_cmplx

!=========================================================================
subroutine scf_loop_x2c(basis,&
                    x_matrix,x_matrix_real,&
                    s_matrix,s_matrix_real,&
                    hamiltonian_hcore,&
                    occupation, &
                    energy, &
                    c_matrix_rel,c_matrix,en_gks,scf_has_converged)
  implicit none

  !=====
  type(basis_set),intent(inout)         :: basis
  real(dp),intent(inout)                :: x_matrix_real(:,:)
  real(dp),intent(inout)                :: s_matrix_real(:,:)
  real(dp),intent(inout)                :: occupation(:,:)
  real(dp),intent(out)                  :: energy(:,:)
  real(dp),intent(inout)                :: c_matrix(:,:,:)
  complex(dp),intent(in)                :: x_matrix(:,:)
  complex(dp),intent(in)                :: s_matrix(:,:)
  complex(dp),intent(in)                :: hamiltonian_hcore(:,:)
  complex(dp),allocatable,intent(inout) :: c_matrix_rel(:,:)
  type(energy_contributions),intent(inout) :: en_gks
  logical,intent(out)                   :: scf_has_converged
  !=====
  integer                 :: nstate
  logical                 :: is_x2c,stopfile_found
  integer                 :: iscf,istate,jstate,jbf,nelectrons
  real(dp)                :: rms
  real(dp),allocatable    :: s_eigval(:)
  real(dp),allocatable    :: energy_vec(:)
  real(dp),allocatable    :: inv_x_matrix(:,:),matrix_tmp(:,:)
  real(dp),allocatable    :: p_matrix_real(:,:)
  complex(dp),allocatable :: occ_matrix_rel(:,:)
  complex(dp),allocatable :: hamiltonian_x2c(:,:)
  complex(dp),allocatable :: p_matrix_rel(:,:),p_matrix_rel_old(:,:)  
  complex(dp),allocatable :: hamiltonian_Vhxc(:,:,:)
  complex(dp),allocatable :: hamiltonian_Vhxc2(:,:,:)
  complex(dp),allocatable :: c_matrix_LaorLb(:,:,:)
  complex(dp),allocatable :: p_matrix_LaorLb(:,:,:)
  complex(dp),allocatable :: ham_hist(:,:,:)
  !=====


  call start_clock(timing_scf)

  is_x2c = .true.
  rms = 1000.0
  nstate = 2*basis%nbf ! = 2 basis%nbf
  nelectrons=nint(sum(occupation(:,1))+sum(occupation(:,2)))
  if(nelectrons>basis%nbf) call die("In X2C the number of electrons cannot be > the num. of basis functions")
  occupation=0.0_dp
  occupation(1:nelectrons,1)=1.0_dp
  occupation(1:nelectrons,2)=1.0_dp
  allocate(occ_matrix_rel(nstate,nstate))
  occ_matrix_rel=COMPLEX_ZERO
  do istate=1,nelectrons
    occ_matrix_rel(istate,istate)=1.0_dp
  enddo

  !
  ! Allocate the main arrays
  call clean_allocate('Total Hamiltonian H',hamiltonian_x2c,nstate,nstate)
  call clean_allocate('Hxc operator VHxc',hamiltonian_Vhxc,basis%nbf,basis%nbf,nspin)
  if(calc_type%need_exchange .or. calc_type%need_exchange_lr) then
    call clean_allocate('Hxc operator VHxc2',hamiltonian_Vhxc2,basis%nbf,basis%nbf,nspin)
  endif
  call clean_allocate('Coefs. La or Lb C',c_matrix_LaorLb,basis%nbf,basis%nbf,nspin)
  call clean_allocate('Density matrix P_LaLb',p_matrix_LaorLb,basis%nbf,basis%nbf,nspin)
  call clean_allocate('Density matrix P',p_matrix_rel,nstate,nstate)
  call clean_allocate('Density matrix P(old)',p_matrix_rel_old,nstate,nstate)
  call clean_allocate('Hamiltonian history',ham_hist,nstate,nstate,2)
  call clean_allocate('State energies',energy_vec,nstate)
  ham_hist=COMPLEX_ZERO 

  !
  ! Setup the grids for the quadrature of DFT potential/energy
  if( calc_type%is_dft ) then
    call init_dft_grid(basis,grid_level,dft_xc(1)%needs_gradient,.TRUE.,BATCH_SIZE)
  endif

  !
  ! Setup the density matrices p_matrix_rel and p_matrix_LaorLb
  c_matrix_LaorLb=COMPLEX_ZERO
  do istate=1,basis%nbf
    do jstate=1,nelectrons
      c_matrix_LaorLb(istate,jstate,1)=c_matrix_rel(2*istate-1,jstate) ! L alpha coef. in spin channel 1
      c_matrix_LaorLb(istate,jstate,2)=c_matrix_rel(2*istate  ,jstate) ! L beta  coef. in spin channel 2
    enddo
  enddo
  p_matrix_rel=matmul(c_matrix_rel,matmul(occ_matrix_rel,transpose(conjg(c_matrix_rel))))
  call setup_density_matrix_cmplx(c_matrix_LaorLb,occupation,p_matrix_LaorLb)


  !
  ! Start the big scf loop
  !
  do iscf=1,nscf
    write(stdout,'(/,1x,a)') '-------------------------------------------'
    write(stdout,'(a,1x,i4,/)') ' *** SCF cycle No:',iscf

    hamiltonian_x2c=COMPLEX_ZERO
    ham_hist(:,:,2)=COMPLEX_ZERO 
    en_gks%exx_hyb=0.0_dp
    en_gks%kin_nuc=REAL(SUM(hamiltonian_hcore(:,:)*p_matrix_rel(:,:)),dp)

    !--Hamiltonian - Hartree ---
    hamiltonian_Vhxc=COMPLEX_ZERO
    call calculate_hamiltonian_hartree_x2c(basis,                   &
                                           occupation,              &
                                           p_matrix_LaorLb,         &
                                           hamiltonian_Vhxc,en_gks)
    do istate=1,nstate/2
      do jstate=1,nstate/2
         hamiltonian_x2c(2*istate-1,2*jstate-1)=hamiltonian_Vhxc(istate,jstate,1) ! La spin channel unbar
         hamiltonian_x2c(2*istate  ,2*jstate  )=hamiltonian_Vhxc(istate,jstate,2) ! Lb spin channel   bar
      enddo
    enddo
    hamiltonian_x2c=hamiltonian_x2c+hamiltonian_hcore

    !--Hamiltonian - Exchange Correlation DFT ---
    if( calc_type%is_dft ) then
      hamiltonian_Vhxc=COMPLEX_ZERO
      call calculate_hamiltonian_xc_x2c(basis,                   &
                                        occupation,              &
                                        c_matrix_LaorLb,         &
                                        hamiltonian_Vhxc,en_gks)
      do istate=1,nstate/2
        do jstate=1,nstate/2
           hamiltonian_x2c(2*istate-1,2*jstate-1)=hamiltonian_x2c(2*istate-1,2*jstate-1)+hamiltonian_Vhxc(istate,jstate,1) 
           hamiltonian_x2c(2*istate  ,2*jstate  )=hamiltonian_x2c(2*istate  ,2*jstate  )+hamiltonian_Vhxc(istate,jstate,1) 
        enddo
      enddo
    endif

    !
    ! LR Exchange + Time-rev. LR Exchange ---
    if(calc_type%need_exchange_lr) then
      hamiltonian_Vhxc=COMPLEX_ZERO
      hamiltonian_Vhxc2=COMPLEX_ZERO
      call setup_lr_exchange_ri_x2c_1(occupation,c_matrix_LaorLb,hamiltonian_Vhxc )
      call setup_lr_exchange_ri_x2c_2(occupation,c_matrix_LaorLb,hamiltonian_Vhxc2)
      do istate=1,nstate/2
        do jstate=1,nstate/2
           ham_hist(2*istate-1,2*jstate-1,2)=hamiltonian_Vhxc(istate,jstate,1) ! < La tensor_product ( La |erf(wr)| La ) tensor_product La >
           ham_hist(2*istate  ,2*jstate  ,2)=hamiltonian_Vhxc(istate,jstate,2) ! < Lb tensor_product ( Lb |erf(wr)| Lb ) tensor_product Lb >
           ham_hist(2*istate  ,2*jstate-1,2)=hamiltonian_Vhxc2(istate,jstate,1) ! < Lb tensor_product ( La |erf(wr)| Lb ) tensor_product La >
           ham_hist(2*istate-1,2*jstate  ,2)=hamiltonian_Vhxc2(istate,jstate,2) ! < La tensor_product ( Lb |erf(wr)| La ) tensor_product Lb >
        enddo
      enddo
      en_gks%exx_hyb=0.5_dp*beta_hybrid*REAL(SUM(ham_hist(:,:,2)*p_matrix_rel(:,:)),dp)
      hamiltonian_x2c(:,:)=hamiltonian_x2c(:,:)+beta_hybrid*ham_hist(:,:,2)
      ham_hist(:,:,2)=COMPLEX_ZERO
    endif

    !--Hamiltonian - Exact Exchange + Time-rev. Exchange ---
    if(calc_type%need_exchange) then
      hamiltonian_Vhxc=COMPLEX_ZERO
      hamiltonian_Vhxc2=COMPLEX_ZERO
      call setup_exchange_ri_x2c_1(occupation,c_matrix_LaorLb,hamiltonian_Vhxc )
      call setup_exchange_ri_x2c_2(occupation,c_matrix_LaorLb,hamiltonian_Vhxc2)
      do istate=1,nstate/2
        do jstate=1,nstate/2
           ham_hist(2*istate-1,2*jstate-1,2)=hamiltonian_Vhxc(istate,jstate,1) ! < La tensor_product ( La | La ) tensor_product La >
           ham_hist(2*istate  ,2*jstate  ,2)=hamiltonian_Vhxc(istate,jstate,2) ! < Lb tensor_product ( Lb | Lb ) tensor_product Lb >
           ham_hist(2*istate  ,2*jstate-1,2)=hamiltonian_Vhxc2(istate,jstate,1) ! < Lb tensor_product ( La | Lb ) tensor_product La >
           ham_hist(2*istate-1,2*jstate  ,2)=hamiltonian_Vhxc2(istate,jstate,2) ! < La tensor_product ( Lb | La ) tensor_product Lb >
        enddo
      enddo
      en_gks%exx_hyb=en_gks%exx_hyb+0.5_dp*alpha_hybrid*REAL(SUM(ham_hist(:,:,2)*p_matrix_rel(:,:)),dp)
      hamiltonian_x2c(:,:)=hamiltonian_x2c(:,:)+alpha_hybrid*ham_hist(:,:,2)
    endif

    !! Sum up to get the total energy
    en_gks%total = en_gks%nuc_nuc + en_gks%kin_nuc + en_gks%hartree + en_gks%exx_hyb + en_gks%xc

    ! Make sure all the MPI tasks have the exact same Hamiltonian
    ! It helps stabilizing the SCF cycles in parallel
    call world%sum(hamiltonian_x2c)
    hamiltonian_x2c(:,:) = hamiltonian_x2c(:,:) / REAL(world%nproc,dp)

    !
    ! Simple mixing on the hamiltonian
    ! the newest is 1
    ! the oldest is 2
    if( iscf > 1) then
      ham_hist(:,:,2)=ham_hist(:,:,1)
    endif
    ham_hist(:,:,1)=hamiltonian_x2c(:,:)
    if( iscf > 1) then
      write(stdout,'(/,1x,a,f8.4)') 'Simple mixing with parameter:',alpha_mixing
      hamiltonian_x2c(:,:) = alpha_mixing * ham_hist(:,:,1) + (1.0_dp - alpha_mixing) * ham_hist(:,:,2)
    endif

    !
    ! H' = X**T * H * X
    ! TODO: slow coding. Use BLAS-level 3 in the future
    hamiltonian_x2c = MATMUL(TRANSPOSE(CONJG(x_matrix)), MATMUL(hamiltonian_x2c, x_matrix))

    !
    ! H' * C' = C' * E
    call diagonalize(' ',hamiltonian_x2c,energy_vec,c_matrix_rel)
    do istate=1,nstate/2
      energy(istate,1)=energy_vec(2*istate-1)
      energy(istate,2)=energy_vec(2*istate)
    enddo

    !
    ! C = X * C'
    c_matrix_rel=MATMUL(x_matrix,c_matrix_rel)

    occupation=0.0_dp
    occupation(1:nelectrons/2,1)=1.0_dp
    occupation(1:nelectrons/2,2)=1.0_dp
    call dump_out_energy('=== Energies ===',occupation,energy,is_x2c=is_x2c)

    call output_homolumo('gKS',occupation,energy,1,basis%nbf)
    occupation=0.0_dp
    occupation(1:nelectrons,1)=1.0_dp
    occupation(1:nelectrons,2)=1.0_dp

    !
    ! Output the total energy and its components
    write(stdout,*)
    write(stdout,'(a25,1x,f19.10)') 'Nucleus-Nucleus (Ha):',en_gks%nuc_nuc
    write(stdout,'(a25,1x,f19.10)') 'Kin+Vext Energy (Ha):',en_gks%kin_nuc
    write(stdout,'(a25,1x,f19.10)') 'Hartree Energy  (Ha):',en_gks%hartree
    if(calc_type%need_exchange) then
      write(stdout,'(a25,1x,f19.10)') 'Exchange Energy (Ha):',en_gks%exx_hyb
    endif
    if(calc_type%is_dft) then
      write(stdout,'(a25,1x,f19.10)') 'XC Energy       (Ha):',en_gks%xc
    endif
    write(stdout,'(/,a25,1x,f19.10,/)') 'Total Energy    (Ha):',en_gks%total

    !
    ! Setup the new density matrix: p_matrix_rel
    ! Save the old one for the convergence criterium
    c_matrix_LaorLb=COMPLEX_ZERO
    do istate=1,basis%nbf
      do jstate=1,nelectrons
        c_matrix_LaorLb(istate,jstate,1)=c_matrix_rel(2*istate-1,jstate) ! L alpha coef. in spin channel 1
        c_matrix_LaorLb(istate,jstate,2)=c_matrix_rel(2*istate  ,jstate) ! L beta  coef. in spin channel 2
      enddo
    enddo
    p_matrix_rel=matmul(c_matrix_rel,matmul(occ_matrix_rel,transpose(conjg(c_matrix_rel))))
    call setup_density_matrix_cmplx(c_matrix_LaorLb,occupation,p_matrix_LaorLb)


 ! SCF convergence check
    if( iscf > 1) then
      rms = NORM2( real(p_matrix_rel(:,:))  - real(p_matrix_rel_old(:,:)) )  * SQRT( REAL(nspin,dp) ) &
          + NORM2( aimag(p_matrix_rel(:,:)) - aimag(p_matrix_rel_old(:,:)) ) * SQRT( REAL(nspin,dp) )
      p_matrix_rel_old(:,:)=p_matrix_rel(:,:)
      write(stdout,'(1x,a,es12.5)') 'Convergence criterium on the density matrix: ',rms
    else
      p_matrix_rel_old(:,:)=p_matrix_rel(:,:)
    endif

    if( rms < tolscf ) then
      scf_has_converged = .TRUE.
      write(stdout,*) ' ===> convergence has been reached'
      write(stdout,*)
    else
      scf_has_converged = .FALSE.
      write(stdout,*) ' ===> convergence not reached yet'
      write(stdout,*)

      if( iscf == nscf ) then
        if( rms > 1.0e-2_dp ) then
          call issue_warning('SCF convergence is very poor')
        else if( rms > 1.0e-4_dp ) then
          call issue_warning('SCF convergence is poor')
        endif
      endif

    endif

    inquire(file='STOP',exist=stopfile_found)

    ! This typically stops too early. We have to use tolscf<0 !
    if( scf_has_converged .OR. stopfile_found ) exit

    !
    ! end of the big SCF loop
  enddo

  !
  ! Store the natural orbital basis representation of the density matrix in c_matrix
  ! and the occupation numbers in occupation(:,1) \in [0:2]
  call clean_allocate('Density matrix P real',p_matrix_real,basis%nbf,basis%nbf)
  p_matrix_real=0.0_dp; occupation=0.0_dp; c_matrix=0.0_dp; energy=0.0_dp;
  !do istate=1,nstate/2
  !  do jstate=1,nstate/2
  !     p_matrix_real(istate,jstate)=real(p_matrix_rel(2*istate-1,2*jstate-1)+p_matrix_rel(2*istate,2*jstate))
  !  enddo
  !enddo
  !p_matrix_real=0.5_dp*(p_matrix_real+transpose(p_matrix_real))
  p_matrix_real(:,:)=real(p_matrix_LaorLb(:,:,1)+p_matrix_LaorLb(:,:,2))
  allocate(s_eigval(basis%nbf),matrix_tmp(basis%nbf,basis%nbf))     
  matrix_tmp(:,:) = s_matrix_real(:,:)
  ! Diagonalization with or without SCALAPACK
  !! S = U*s*U^H
  call diagonalize_scalapack(scf_diago_flavor,scalapack_block_min,matrix_tmp,s_eigval)
  nstate = COUNT( s_eigval(:) > min_overlap )
  call clean_allocate('Overlap INV_X * INV_X**H = S',inv_x_matrix,basis%nbf,nstate)
  inv_x_matrix=0.0_dp
  write(stdout,'(/,a)')       ' Filtering basis functions that induce overcompleteness'
  write(stdout,'(a,es9.2)')   '   Lowest S eigenvalue is           ',MINVAL( s_eigval(:) )
  write(stdout,'(a,es9.2)')   '   Tolerance on overlap eigenvalues ',min_overlap
  write(stdout,'(a,i5,a,i5)') '   Retaining ',nstate,' among ',basis%nbf
  !! INV_X = U*s^(1/2)
  istate = 0
  do jbf=1,basis%nbf
    if( s_eigval(jbf) > min_overlap ) then
      istate = istate + 1
      inv_x_matrix(:,istate) = matrix_tmp(:,jbf) * SQRT( s_eigval(jbf) )
    endif
  enddo
  !  diag[ S^1/2 P S^1/2 ] -> V
  deallocate(matrix_tmp)     
  allocate(matrix_tmp(nstate,nstate))     
  matrix_tmp=0.0_dp
  matrix_tmp=MATMUL(TRANSPOSE(inv_x_matrix), MATMUL(p_matrix_real, inv_x_matrix))
  call diagonalize(' ',matrix_tmp,occupation(:,1),c_matrix(:,:,1))
  !  C = S^-1/2 V
  c_matrix(:,:,1) = MATMUL(x_matrix_real, c_matrix(:,:,1))
  deallocate(matrix_tmp,s_eigval)
  if(nspin==2) then
   occupation(:,2)=0.0_dp
   c_matrix(:,:,2)=0.0_dp
  endif

  !
  ! Cleanly deallocate the integral grid information
  !
  if( calc_type%is_dft ) call destroy_dft_grid()


  !
  ! Cleanly deallocate the arrays
  !
  call clean_deallocate('Overlap INV_X * INV_X**H = S real',inv_x_matrix)
  call clean_deallocate('Density matrix P real',p_matrix_real)
  call clean_deallocate('Density matrix P',p_matrix_rel)
  call clean_deallocate('Density matrix P_LaLb',p_matrix_LaorLb)
  call clean_deallocate('Density matrix P(old)',p_matrix_rel_old)
  call clean_deallocate('Coefs. La or Lb C',c_matrix_LaorLb)
  call clean_deallocate('Hxc operator VHxc',hamiltonian_Vhxc)
  call clean_deallocate('Total Hamiltonian H',hamiltonian_x2c)
  call clean_deallocate('Hamiltonian history',ham_hist)
  call clean_deallocate('State energies',energy_vec)
  if(calc_type%need_exchange) then
    call clean_deallocate('Hxc operator VHxc2',hamiltonian_Vhxc2)
  endif

  write(stdout,'(/,/,a25,1x,f19.10,/)') 'SCF Total Energy (Ha):',en_gks%total

  if( print_yaml_ .AND. is_iomaster ) then
    if( scf_has_converged ) then
      write(unit_yaml,'(/,a)') 'scf is converged: True'
    else
      write(unit_yaml,'(/,a)') 'scf is converged: False'
    endif
    call print_energy_yaml('scf energy',en_gks)
    call dump_out_energy_yaml('gks energies',energy)
  endif

  call stop_clock(timing_scf)

end subroutine scf_loop_x2c

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
  if( .NOT. calc_type%need_exchange ) then
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
  integer                 :: nstate,nocc,istate,ispin
  real(dp),allocatable    :: c_matrix_restart(:,:,:)
  real(dp),allocatable    :: h_ii(:,:)
  real(dp),allocatable    :: energy_restart(:,:),occupation_restart(:,:)
  !=====

  nstate = SIZE(c_matrix,DIM=2)
  nocc   = get_number_occupied_states(occupation)

  !
  ! Get the exchange operator if not already calculated
  !
  if( ALL( ABS(hamiltonian_exx(:,:,:)) < 1.0e-6_dp ) ) then
    call calculate_exchange(basis,p_matrix,hamiltonian_exx,occupation=occupation,c_matrix=c_matrix)
  endif

  call clean_allocate('RESTART: C',c_matrix_restart,basis%nbf,nstate,nspin)
  allocate(energy_restart(nstate,nspin))
  allocate(occupation_restart(nstate,nspin))
  allocate(h_ii(nstate,nspin))

  call read_restart(restart_type,'RESTART_TEST',basis,occupation_restart,c_matrix_restart,energy_restart)

  if( restart_type == NO_RESTART ) then
    c_matrix_restart(:,:,:) = c_matrix(:,:,:)
  else
    write(stdout,'(1x,a,a)') 'RESTART file read: ','RESTART_TEST'
  endif

  call matrix_ao_to_mo_diag(c_matrix_restart,hamiltonian_hartree,h_ii)
  call dump_out_energy('=== Hartree expectation value ===',occupation,h_ii)
  call dump_out_energy_yaml('hartree expectation value',h_ii,1,nstate)
  write(stdout,'(1x,a,2(3x,f12.6))') 'Hartree  HOMO expectation (eV):',h_ii(nocc,:) * Ha_eV


  call matrix_ao_to_mo_diag(c_matrix_restart,hamiltonian_exx,h_ii)
  call dump_out_energy('=== Exchange expectation value ===',occupation,h_ii)
  call dump_out_energy_yaml('exchange expectation value',h_ii,1,nstate)
  write(stdout,'(1x,a,2(3x,f12.6))') 'Exchange HOMO expectation (eV):',h_ii(nocc,:) * Ha_eV

  deallocate(h_ii)
  deallocate(energy_restart,occupation_restart)
  call clean_deallocate('RESTART: C',c_matrix_restart)

end subroutine print_hartee_expectation


!=========================================================================
subroutine print_expectations(basis,c_matrix,hkin)
  implicit none

  type(basis_set),intent(inout) :: basis
  real(dp),intent(in)           :: c_matrix(:,:,:)
  real(dp),intent(in)           :: hkin(:,:)
  !=====
  real(dp),allocatable       :: p_matrix(:,:,:)
  real(dp),allocatable       :: hh(:,:),ekin(:,:)
  real(dp)                   :: ehartree_i
  integer                    :: ibf,jbf,nbf,nstate,istate,ispin
  character(len=6)           :: char6
  !=====

  write(stdout,*) 'Print expectations'
  nbf    = SIZE(c_matrix,DIM=1)
  nstate = SIZE(c_matrix,DIM=2)

  allocate(hh(nbf,nbf))
  allocate(p_matrix(nbf,nbf,nspin))
  allocate(ekin(nstate,nspin))

  call matrix_ao_to_mo_diag(c_matrix,hkin,ekin)


  if( print_yaml_ .AND. is_iomaster ) then
    write(unit_yaml,'(/,a)') 'kinetic expectation value:'
    write(unit_yaml,'(4x,a)') 'unit: Ha'
    do ispin=1,nspin
      write(unit_yaml,'(4x,a,i2,a)') 'spin channel',ispin,':'
      do istate=1,nstate
        write(char6,'(i6)') istate
        write(unit_yaml,'(8x,a6,a,1x,es18.8)') ADJUSTL(char6),':',ekin(istate,ispin)
      enddo
    enddo


    write(unit_yaml,'(/,a)') 'hartree self-interaction:'
    write(unit_yaml,'(4x,a)') 'unit: Ha'
  endif

  do ispin=1,nspin
    if( print_yaml_ .AND. is_iomaster ) write(unit_yaml,'(4x,a,i2,a)') 'spin channel',ispin,':'
    do istate=1,nstate

      p_matrix(:,:,:) = 0.0_dp
      call DSYRK('L','N',nbf,1,1.0d0,c_matrix(1,istate,ispin),nbf,0.0d0,p_matrix(1,1,ispin),nbf)
      ! Symmetrize
      do jbf=1,nbf
        do ibf=jbf+1,nbf
          p_matrix(jbf,ibf,ispin) = p_matrix(ibf,jbf,ispin)
        enddo
      enddo

      call calculate_hartree(basis,p_matrix,hh,eh=ehartree_i)

      write(stdout,'(1x,a,i5,es16.6,1x,es16.6)') 'Eh i , T_i (Ha): ',ispin,ehartree_i,ekin(istate,ispin)
      if( print_yaml_ .AND. is_iomaster ) then
        write(char6,'(i6)') istate
        write(unit_yaml,'(8x,a6,a,1x,es18.8)') ADJUSTL(char6),':',ehartree_i
      endif

    enddo
  enddo

  deallocate(hh,p_matrix,ekin)


end subroutine print_expectations


!=========================================================================
end module m_scf_loop
!=========================================================================
