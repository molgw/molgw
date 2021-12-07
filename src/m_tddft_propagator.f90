!=========================================================================
! This file is part of MOLGW.
! Author: Ivan Maliyov
!
! This module contains
! the time propagation of the KS wavefunctions for TDDFT
!
!=========================================================================
module m_tddft_propagator
 use m_definitions
 use m_memory
 use m_warning
 use m_timing
 use m_tddft_variables
 use m_atoms
 use m_string_tools
 use m_multipole
 use m_basis_set
 use m_hamiltonian_tools
 use m_hamiltonian_onebody
 use m_hamiltonian_wrapper
 use m_inputparam
 use m_dft_grid
 use m_linear_algebra
 use m_density_tools
 use m_dft_grid
 use m_scf
 use m_io

 interface propagate_orth
  module procedure propagate_orth_ham_1
  module procedure propagate_orth_ham_2
 end interface propagate_orth

 integer,private                    :: nocc
 real(dp),private                   :: dipole(3)
 real(dp),private                   :: time_read
 real(dp),allocatable,private       :: xatom_start(:,:)
 real(dp),private                   :: excit_field_norm
 !==hamiltonian extrapolation variables==
 real(dp),allocatable,private       :: extrap_coefs(:)
 complex(dp),allocatable,private    :: h_small_hist_cmplx(:,:,:,:)
 complex(dp),allocatable,private    :: c_matrix_orth_hist_cmplx(:,:,:,:)
 integer,private            :: ntau
 !==frozen core==
 real(dp),allocatable       :: energies_start(:,:)
 complex(dp),allocatable    :: a_matrix_orth_start_cmplx(:,:,:)
 !====

 type(energy_contributions),private :: en_tddft

contains


!=========================================================================
subroutine calculate_propagation(basis,occupation,c_matrix,restart_tddft_is_correct)
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: c_matrix(:,:,:)
 real(dp),intent(inout)     :: occupation(:,:)
 logical,intent(in)         :: restart_tddft_is_correct
!=====
 integer                    :: fixed_atom_list(ncenter_nuclei-nprojectile)
 integer                    :: ispin
 integer                    :: istate,nstate_tmp
 integer                    :: nwrite_step
 real(dp)                   :: time_min
 real(dp)                   :: xprojectile(3)
 real(dp),allocatable       :: dipole_ao(:,:,:)
 real(dp),allocatable       :: s_matrix(:,:)
 real(dp),allocatable       :: x_matrix(:,:)
 real(dp),allocatable       :: hamiltonian_kinetic(:,:)
 real(dp),allocatable       :: hamiltonian_nucleus(:,:)
!=====initial values
 integer                    :: nstate
 real(dp),allocatable       :: energy_tddft(:)
 complex(dp),allocatable    :: c_matrix_cmplx(:,:,:)
 complex(dp),allocatable    :: c_matrix_orth_cmplx(:,:,:)
 complex(dp),allocatable    :: h_cmplx(:,:,:)
 complex(dp),allocatable    :: h_small_cmplx(:,:,:)
!=====TDDFT loop variables=============================
 character(len=8)           :: time_key
 integer                    :: iatom
 integer                    :: itau
 integer                    :: iwrite_step
 integer                    :: file_time_data,file_excit_field
 integer                    :: file_dipole_time
 real(dp)                   :: time_cur,time_one_iter
 complex(dp),allocatable    :: p_matrix_cmplx(:,:,:)
 logical                    :: is_identity_ ! keep this varibale
!==cube_diff varibales====================================
 real(dp),allocatable       :: cube_density_start(:,:,:,:)
 integer                    :: nx,ny,nz
 logical                    :: file_exists
!==line density diff variables=============================
 real(dp),allocatable       :: rho_start(:,:)
 real(dp)                   :: point_a(3),point_b(3)
 integer                    :: nr_line_rho
!==qmatrix==
 integer,allocatable        :: istate_cut(:,:)
 integer                    :: file_q_matrix(2)
 integer                    :: iocc
 complex(dp),allocatable    :: c_matrix_orth_start_complete_cmplx(:,:,:)
!====

 call start_clock(timing_tddft_loop)

 write(stdout,'(/,/,1x,a)') '=================================================='
 write(stdout,'(x,a,/)')    'RT-TDDFT simulation'

 ! Here this occupation comes from the set_occupation subroutine with zero temperature
 nstate = SIZE(occupation(:,:),DIM=1)

 if( read_tddft_restart_ .AND. restart_tddft_is_correct ) then
   nocc = get_nocc_from_restart()
 else
   nocc = get_number_occupied_states(occupation)
 end if

 call echo_tddft_variables()

 call clean_allocate('Overlap matrix S for TDDFT',s_matrix,basis%nbf,basis%nbf)
 call clean_allocate('Kinetic operator T for TDDFT',hamiltonian_kinetic,basis%nbf,basis%nbf)
 call clean_allocate('Nucleus operator V for TDDFT',hamiltonian_nucleus,basis%nbf,basis%nbf)

 call setup_overlap(basis,s_matrix)

 ! x_matrix is now allocated with dimension (basis%nbf,nstate))
 call setup_x_matrix(min_overlap,s_matrix,nstate_tmp,x_matrix)
 if( nstate /= nstate_tmp ) then
   call die('Error with nstate in the TDDFT propagator')
 end if

 call nucleus_nucleus_energy(en_tddft%nuc_nuc)

 !
 ! Setup the fixed part of the Hamiltonian: the kinetic energy and the fixed nuclei potential
 !
 call setup_kinetic(basis,hamiltonian_kinetic)

 !
 ! hamiltonian_nucleus contains the contribution from all the fixed atoms (i.e. excluding the projectile)
 ! Remember: the projectile is always the last atom
 do iatom=1,ncenter_nuclei-nprojectile
   fixed_atom_list(iatom) = iatom
 enddo
 call setup_nucleus(basis,hamiltonian_nucleus,fixed_atom_list)

 if( nelement_ecp > 0 ) then
   call setup_nucleus_ecp(basis,hamiltonian_nucleus)
 endif


 if(write_step / time_step - NINT( write_step / time_step ) > 0.0E-10_dp .OR. write_step < time_step ) then
   call die("Tddft error: write_step is not a multiple of time_step or smaller than time_step.")
 end if

 if( calc_type%is_dft ) then
    call init_dft_grid(basis,grid_level,dft_xc(1)%needs_gradient,.TRUE.,BATCH_SIZE)
 endif

 call clean_allocate('Wavefunctions C for TDDFT',c_matrix_cmplx,basis%nbf,nocc,nspin)
 call clean_allocate('Wavefunctions hist. C for TDDFT',c_matrix_orth_cmplx,nstate,nocc,nspin)
 call clean_allocate('Hamiltonian for TDDFT',h_cmplx,basis%nbf,basis%nbf,nspin)
 call clean_allocate('h_small_cmplx for TDDFT',h_small_cmplx,nstate,nstate,nspin)
 call clean_allocate('p_matrix_cmplx for TDDFT',p_matrix_cmplx,basis%nbf,basis%nbf,nspin)


 allocate(xatom_start(3,ncenter_nuclei))

 write(stdout,'(/,1x,a)') "===INITIAL CONDITIONS==="
 ! Getting c_matrix_cmplx(t=0) whether using RESTART_TDDFT file, whether using real c_matrix
 if( read_tddft_restart_ .AND. restart_tddft_is_correct ) then
   ! assign xatom_start, c_matrix_orth_cmplx, time_min with values given in RESTART File
   call read_restart_tddft(nstate,time_read,occupation,c_matrix_orth_cmplx)
   time_min = time_read
   do ispin=1,nspin
     c_matrix_cmplx(:,:,ispin) = MATMUL( x_matrix(:,:) , c_matrix_orth_cmplx(:,:,ispin) )
   end do
   xprojectile = xatom(:,ncenter_nuclei)
   call change_position_one_atom(ncenter_nuclei,xprojectile)
 else
   c_matrix_cmplx(:,:,:) = c_matrix(:,1:nocc,:)
   xatom_start=xatom
   time_min=0.0_dp
 end if

 ! Getting starting value of the Hamiltonian
 ! itau=0 to avoid excitation calculation
 call setup_hamiltonian_cmplx(basis,                   &
                              nstate,                  &
                              0,                       &    ! itau
                              time_min,                &
                              0.0_dp,                  &    ! time_cur
                              occupation,              &
                              c_matrix_cmplx,          &
                              hamiltonian_kinetic,     &
                              hamiltonian_nucleus,     &
                              h_small_cmplx,           &
                              x_matrix,                &
                              dipole_ao,               &
                              h_cmplx,en_tddft)

 ! In case of no restart, find the c_matrix_orth_cmplx by diagonalizing h_small
 if( (.NOT. read_tddft_restart_) .OR. (.NOT. restart_tddft_is_correct)) then
   call clean_allocate('c_matrix_buf for TDDFT',c_matrix_orth_start_complete_cmplx,nstate,nstate,nspin)
   allocate(energy_tddft(nstate))
   do ispin=1, nspin
     call diagonalize(postscf_diago_flavor,h_small_cmplx(:,:,ispin),energy_tddft,c_matrix_orth_start_complete_cmplx(:,:,ispin))
   end do
   ! in order to save the memory, we dont keep inoccupied states (nocc+1:nstate)
   c_matrix_orth_cmplx(1:nstate,1:nocc,1:nspin)=c_matrix_orth_start_complete_cmplx(1:nstate,1:nocc,1:nspin)
   call clean_deallocate('c_matrix_buf for TDDFT',c_matrix_orth_start_complete_cmplx)
   deallocate(energy_tddft)
 end if

 ! Number of iterations
 ntau = NINT( (time_sim-time_min) / time_step )
 nwrite_step = NINT( (time_sim - time_min) / write_step )


 if(excit_type%form==EXCIT_LIGHT) then
   call clean_allocate('Dipole_basis for TDDFT',dipole_ao,basis%nbf,basis%nbf,3)
   call calculate_dipole_ao(basis,dipole_ao)
 end if

 if( print_dens_traj_tddft_ ) then
   call plot_rho_traj_bunch_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,0,0.d0)
 end if


 if( calc_q_matrix_ ) then
   call initialize_q(nstate,nocc,nspin,c_matrix_orth_start_complete_cmplx,h_small_cmplx,istate_cut,file_q_matrix)
 end if

 !==frozen core: energy initialization==
 if( ncore_tddft > 0 ) then
   call clean_allocate('Inial energies for the frozen core',energies_start,nstate,nspin)
   call clean_allocate('a_matrix_orth_start_cmplx for the frozen core',a_matrix_orth_start_cmplx,nstate,nstate,nspin)
   do ispin=1, nspin
     call diagonalize(postscf_diago_flavor,h_small_cmplx(:,:,ispin),energies_start(:,ispin),a_matrix_orth_start_cmplx(:,:,ispin))
   end do
 end if
 !====

!===cube_diff matrix allocation and parameters initialization
 if(print_cube_diff_tddft_) then
   call initialize_cube_diff_cmplx(nx,ny,nz)
   allocate(cube_density_start(nx,ny,nz,nspin))
 end if

 if(print_line_rho_diff_tddft_) then
   call initialize_rho_diff_cmplx(nr_line_rho,point_a,point_b)
   allocate(rho_start(nr_line_rho,nspin))
 end if

 time_min = time_read

 !
 ! Opening files and writing headers in files
 call initialize_files(file_time_data,file_dipole_time,file_excit_field)

 !
 ! Printing initial values of energy and dipole taken from SCF or RESTART_TDDFT
 en_tddft%total = en_tddft%nucleus + en_tddft%kinetic + en_tddft%nuc_nuc &
                 + en_tddft%hartree + en_tddft%exx_hyb + en_tddft%xc + en_tddft%excit

 if( excit_type%form == EXCIT_LIGHT ) then
   call setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)
   call static_dipole(basis,p_matrix_in=p_matrix_cmplx,dipole_ao_in=dipole_ao,dipole_out=dipole)
 endif


 if( print_cube_rho_tddft_ ) call plot_cube_wfn_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,0)

 if( print_cube_diff_tddft_ ) then
   call calc_cube_initial_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,cube_density_start,nx,ny,nz)
   call plot_cube_diff_parallel_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,0,cube_density_start,nx,ny,nz)
 end if

 if(print_line_rho_diff_tddft_) then
   call calc_rho_initial_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,0,time_min,nr_line_rho,point_a,point_b,rho_start)
   call plot_rho_diff_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,0,time_min,nr_line_rho,point_a,point_b,rho_start)
 end if

 if( calc_dens_disc_ )       call calc_density_in_disc_cmplx_dft_grid(basis,occupation,c_matrix_cmplx,0,time_min)
! if( calc_dens_disc_ )       call calc_density_in_disc_cmplx_regular(nstate,nocc,basis,occupation,c_matrix_cmplx,0,time_min)

 if( print_line_rho_tddft_ ) call plot_rho_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,0,time_min)

 call print_tddft_values(time_min,file_time_data,file_dipole_time,file_excit_field,0)
 en_tddft%time = time_min
 write(time_key,'(i8)') 0
 call print_energy_yaml('tddft energy '//TRIM(ADJUSTL(time_key)),en_tddft)

 time_min = time_min + time_step


 !
 ! Extrapolation coefficients and history c_ and h_ matrices (h_small_hist_cmplx)
 call initialize_extrap_coefs(c_matrix_orth_cmplx,h_small_cmplx)

 write(stdout,'(1x,a,/)') "===END OF INITIAL CONDITIONS==="


 !
 ! TDDFT time loop
 !
 time_cur = time_min
 iwrite_step = 1
 itau = 1
 in_tddft_loop = .TRUE.
 do while ( (time_cur - time_sim) < 1.0e-10 )
   if(itau==3) call start_clock(timing_tddft_one_iter)


   !
   ! Use c_matrix_orth_cmplx and h_small_cmplx at (time_cur-time_step) as start values,
   ! than use chosen predictor-corrector sheme to calculate c_matrix_cmplx, c_matrix_orth_cmplx,
   ! h_cmplx and h_small_cmplx and time_cur.
   call predictor_corrector(basis,                  &
                            c_matrix_cmplx,         &
                            c_matrix_orth_cmplx,    &
                            h_cmplx,                &
                            h_small_cmplx,          &
                            x_matrix,               &
                            itau,                   &
                            time_cur,               &
                            occupation,             &
                            hamiltonian_kinetic,    &
                            hamiltonian_nucleus,    &
                            dipole_ao)


   !
   ! debug
   !call check_identity_cmplx(nocc,nocc,MATMUL(MATMUL(TRANSPOSE(CONJG(c_matrix_cmplx(:,:,nspin))),s_matrix(:,:)), c_matrix_cmplx(:,:,nspin) ),is_identity_)
   !if(.NOT. is_identity_) then
   !  write(stdout,*) "C**H*S*C is not identity at itau= ", itau
   !end if

   !
   ! Print tddft values into diferent files: 1) standart output; 2) time_data.dat; 3) dipole_time.dat; 4) excitation_time.dat.
   ! 3) and 4) in case of light excitation
   if( ABS(time_cur / (write_step)- NINT(time_cur / (write_step))) < 1.0e-7 ) then

     if(excit_type%form==EXCIT_LIGHT) then
      call setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)
      call static_dipole(basis,p_matrix_in=p_matrix_cmplx,dipole_ao_in=dipole_ao,dipole_out=dipole)
     end if

     en_tddft%total = en_tddft%nucleus + en_tddft%kinetic + en_tddft%nuc_nuc &
                     + en_tddft%hartree + en_tddft%exx_hyb + en_tddft%xc + en_tddft%excit

     call print_tddft_values(time_cur,file_time_data,file_dipole_time,file_excit_field,itau)
     en_tddft%time = time_cur
     write(time_key,'(i8)') iwrite_step
     call print_energy_yaml('tddft energy '//TRIM(ADJUSTL(time_key)),en_tddft)

     if( print_line_rho_tddft_  )     call plot_rho_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,iwrite_step,time_cur)
     if( print_line_rho_diff_tddft_ ) call plot_rho_diff_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx, &
                                                               iwrite_step,time_cur,nr_line_rho,point_a,point_b,rho_start)
     if( print_cube_rho_tddft_  )     call plot_cube_wfn_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,iwrite_step)
     if( print_cube_diff_tddft_ )     call plot_cube_diff_parallel_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx, &
                                                                         iwrite_step,cube_density_start,nx,ny,nz)
     if( calc_dens_disc_ )            call calc_density_in_disc_cmplx_dft_grid(basis,occupation,c_matrix_cmplx, &
                                                                               iwrite_step,time_cur)
!     if( calc_dens_disc_ )       call calc_density_in_disc_cmplx_regular(nstate,nocc,basis,occupation,c_matrix_cmplx,iwrite_step,time_cur)

     if(calc_q_matrix_) call calculate_q_matrix(occupation,c_matrix_orth_start_complete_cmplx,c_matrix_orth_cmplx, &
                                                istate_cut,file_q_matrix,time_cur)

     iwrite_step = iwrite_step + 1

   end if

   !
   !--TIMING of one iteration--
   if(itau==3) then
     call stop_clock(timing_tddft_one_iter)
     call output_timing_one_iter()
   end if

   ! ---print tdddft restart each n_restart_tddft steps---
   if( print_tddft_restart_ .AND. mod(itau,n_restart_tddft)==0 ) then
     call write_restart_tddft(nstate,time_cur,occupation,c_matrix_orth_cmplx)
   end if

  time_cur = time_min + itau*time_step
  itau = itau + 1
!---
 end do
 in_tddft_loop=.FALSE.
!********end time loop*******************

 if(print_tddft_restart_) then
   !time_cur-time_step to be consistent with the actual last moment of the simulation
   call write_restart_tddft(nstate,time_cur-time_step,occupation,c_matrix_orth_cmplx)
 end if

 if( is_iomaster) then
   close(file_time_data)
   if( excit_type%form==EXCIT_LIGHT) then
     close(file_dipole_time)
     close(file_excit_field)
   end if
 end if

 call clean_deallocate('p_matrix_cmplx for TDDFT',p_matrix_cmplx)
 call clean_deallocate('h_small_hist_cmplx for TDDFT',h_small_hist_cmplx)
 call clean_deallocate('c_matrix_orth_hist_cmplx for TDDFT',c_matrix_orth_hist_cmplx)

 if(ncore_tddft > 0) then
   call clean_deallocate('a_matrix_orth_start_cmplx for the frozen core',a_matrix_orth_start_cmplx)
   call clean_deallocate('Inial energies for frozen core',energies_start)
 end if

 if(ALLOCATED(extrap_coefs)) deallocate(extrap_coefs)
 if(ALLOCATED(cube_density_start)) deallocate(cube_density_start)

 if( calc_type%is_dft ) call destroy_dft_grid()

 call clean_deallocate('Overlap matrix S for TDDFT',s_matrix)
 call clean_deallocate('Overlap sqrt S^{-1/2}',x_matrix)
 call clean_deallocate('Kinetic operator T for TDDFT',hamiltonian_kinetic)
 call clean_deallocate('Nucleus operator V for TDDFT',hamiltonian_nucleus)

 deallocate(xatom_start)

 call clean_deallocate('Dipole_basis for TDDFT',dipole_ao)

 call clean_deallocate('Wavefunctions C for TDDFT',c_matrix_cmplx)
 call clean_deallocate('Wavefunctions hist. C for TDDFT',c_matrix_orth_cmplx)
 call clean_deallocate('Hamiltonian for TDDFT',h_cmplx)
 call clean_deallocate('h_small_cmplx for TDDFT',h_small_cmplx)

 write(stdout,'(/,x,a)') "End of RT-TDDFT simulation"
 write(stdout,'(1x,a,/)') '=================================================='
 call stop_clock(timing_tddft_loop)

end subroutine calculate_propagation


!=========================================================================
subroutine echo_tddft_variables()
 implicit none

 write(stdout,'(/,1x,a)') 'The most important variables of this section:'
 write(stdout,'(2x,a32,2x,es16.8)') 'Simulation time: time_sim',time_sim
 write(stdout,'(2x,a32,2x,es16.8)') 'Time step: time_step',time_step
 write(stdout,'(2x,a32,2x,i8)') 'Number of iterations: ntau',NINT((time_sim)/time_step)
 write(stdout,'(2x,a32,6x,a)')      'Predictor-corrector: pred_corr',pred_corr
 write(stdout,'(2x,a32,6x,a)')      'Propagator: prop_type',prop_type
 write(stdout,'(2x,a32,2x,i8)')     'Number of occupied states: nocc',nocc
 write(stdout,'(2x,a32,2x,i8,/)')     'Hamiltonian history length: n_hist',n_hist

end subroutine echo_tddft_variables


!=========================================================================
subroutine output_timing_one_iter()
 implicit none
 real(dp)           :: time_one_iter
!=====

  time_one_iter = get_timing(timing_tddft_one_iter)
  write(stdout,'(/,1x,a)') '**********************************'
  write(stdout,"(1x,a30,2x,es14.6,1x,a)") "Time of one iteration is", time_one_iter,"s"
  write(stdout,"(1x,a30,2x,3(f12.2,1x,a))") "Estimated calculation time is", time_one_iter*ntau, "s  = ", &
                                            time_one_iter*ntau/60.0_dp, &
                                            "min  = ", time_one_iter*ntau/3600.0_dp, "hrs"
  write(stdout,'(1x,a)') '**********************************'
  flush(stdout)

end subroutine output_timing_one_iter


!=========================================================================
subroutine predictor_corrector(basis,                  &
                               c_matrix_cmplx,         &
                               c_matrix_orth_cmplx,    &
                               h_cmplx,                &
                               h_small_cmplx,          &
                               x_matrix,               &
                               itau,                   &
                               time_cur,               &
                               occupation,             &
                               hamiltonian_kinetic,    &
                               hamiltonian_nucleus,    &
                               dipole_ao)

 implicit none
 type(basis_set),intent(in)      :: basis
 complex(dp),intent(out)         :: c_matrix_cmplx(:,:,:)
 complex(dp),intent(inout)       :: c_matrix_orth_cmplx(:,:,:)
 complex(dp),intent(out)         :: h_cmplx(:,:,:)
 complex(dp),intent(inout)       :: h_small_cmplx(:,:,:)
 real(dp),intent(in)             :: x_matrix(:,:)
 integer,intent(in)              :: itau
 real(dp),intent(in)             :: time_cur
 real(dp),intent(in)             :: occupation(:,:)
 real(dp),intent(in)             :: hamiltonian_kinetic(:,:)
 real(dp),intent(inout)          :: hamiltonian_nucleus(:,:)
 real(dp),allocatable,intent(in) :: dipole_ao(:,:,:)
!=====
 integer      :: nstate,iextr,i_iter,file_iter_norm
!=====

 nstate = SIZE(c_matrix_orth_cmplx,DIM=1)

 write(stdout,'(/,1x,a)') 'PREDICTOR-CORRECTOR BLOCK'

 select case (pred_corr)
 ! ///////////////////////////////////
 case('PC0')
   call propagate_orth(nstate,basis,time_step,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_cmplx,x_matrix,prop_type)
   call setup_hamiltonian_cmplx(basis,                   &
                                nstate,                  &
                                itau,                    &
                                time_cur,                &
                                time_step,               &
                                occupation,              &
                                c_matrix_cmplx,          &
                                hamiltonian_kinetic,     &
                                hamiltonian_nucleus,     &
                                h_small_cmplx,           &
                                x_matrix,                &
                                dipole_ao,               &
                                h_cmplx,en_tddft)


 ! ///////////////////////////////////
 ! Following Cheng and Van Voorhis, Phys Rev B 74 (2006)
 case('PC1')
   !--1--PREDICTOR----| H(2/4),H(6/4)-->H(9/4)
       h_small_cmplx= -3.0_dp/4.0_dp*h_small_hist_cmplx(:,:,:,1)+7.0_dp/4.0_dp*h_small_hist_cmplx(:,:,:,2)
   !--2--PREDICTOR----| C(8/4)---U[H(9/4)]--->C(10/4)
   call propagate_orth(nstate,basis,time_step/2.0_dp,c_matrix_orth_hist_cmplx(:,:,:,1), &
                       c_matrix_cmplx,h_small_cmplx,x_matrix,prop_type)

   !--3--CORRECTOR----| C(10/4)-->H(10/4)
   call setup_hamiltonian_cmplx(basis,                   &
                                nstate,                  &
                                itau,                    &
                                time_cur-time_step/2.0_dp, &
                                time_step,               &
                                occupation,              &
                                c_matrix_cmplx,          &
                                hamiltonian_kinetic,     &
                                hamiltonian_nucleus,     &
                                h_small_cmplx,           &
                                x_matrix,                &
                                dipole_ao,               &
                                h_cmplx,en_tddft)

   !--4--PROPAGATION----| C(8/4)---U[H(10/4)]--->C(12/4)
   call propagate_orth(nstate,basis,time_step,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_cmplx,x_matrix,prop_type)

   !--5--UPDATE----| C(12/4)-->C(8/4); H(6/4)-->H(2/4); H(10/4)-->H(6/4)
   c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
   h_small_hist_cmplx(:,:,:,1)=h_small_hist_cmplx(:,:,:,2)
   h_small_hist_cmplx(:,:,:,2)=h_small_cmplx(:,:,:)


 ! ///////////////////////////////////
 case('PC2B')
   h_small_cmplx=(0.0_dp,0.0_dp)
   !--0--EXTRAPOLATE----
   do iextr=1,n_hist
     h_small_cmplx=h_small_cmplx+extrap_coefs(iextr)*h_small_hist_cmplx(:,:,:,iextr)
   end do
   ! n_hist-2  <---  n_hist; n_hist-3 <-- n_hist-1
   if(n_hist > 2) then
     do iextr=1,n_hist-2
       h_small_hist_cmplx(:,:,:,iextr)=h_small_hist_cmplx(:,:,:,iextr+2)
     end do
   end if

   !--1--PROPAGATE----| C(t)--U[H(1/4dt)]-->C(t+dt/2)
   call propagate_orth(nstate,basis,time_step/2.0_dp,c_matrix_orth_hist_cmplx(:,:,:,1), &
                       c_matrix_cmplx,h_small_cmplx,x_matrix,prop_type)
   !--2--CALCULATE- H(t+dt/4)
   call setup_hamiltonian_cmplx(basis,                   &
                                nstate,                  &
                                itau,                    &
                                time_cur-time_step/2.0_dp, &
                                time_step,               &
                                occupation,              &
                                c_matrix_cmplx,          &
                                hamiltonian_kinetic,     &
                                hamiltonian_nucleus,     &
                                h_small_cmplx,           &
                                x_matrix,                &
                                dipole_ao,               &
                                      h_cmplx,en_tddft)

   if (n_hist > 1) h_small_hist_cmplx(:,:,:,n_hist-1)=h_small_cmplx
   !--3--PROPAGATION----|
   call propagate_orth(nstate,basis,time_step,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_cmplx,x_matrix,prop_type)

   call setup_hamiltonian_cmplx(basis,                   &
                                nstate,                  &
                                itau,                    &
                                time_cur,                &
                                time_step,               &
                                occupation,              &
                                c_matrix_cmplx,          &
                                hamiltonian_kinetic,     &
                                hamiltonian_nucleus,     &
                                h_small_cmplx,           &
                                x_matrix,                &
                                dipole_ao,               &
                                h_cmplx,en_tddft)

   !--5--UPDATE----|
   c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
   h_small_hist_cmplx(:,:,:,n_hist)=h_small_cmplx


 ! ///////////////////////////////////
 ! Iterative propagation with Lagrange interpolation
 ! ---------------|t-dt|-----------------|t|----------|t+dt/2)|----------|t+dt|
 !............|H(n_hist-1)|..........|H(n_hist)|....|H(n_hist+1)|.....|H(n_hist+2)|
 case('PC3','PC4')
   h_cmplx=(0.0_dp,0.0_dp)
   h_small_hist_cmplx(:,:,:,n_hist+2)=(0.0_dp,0.0_dp)
   !--1--EXTRAPOLATION WITH PREVIOUS STEPS----
   do iextr=1,n_hist
     h_small_hist_cmplx(:,:,:,n_hist+2)=h_small_hist_cmplx(:,:,:,n_hist+2)+extrap_coefs(iextr)*h_small_hist_cmplx(:,:,:,iextr)
   end do

   !--2--LOCAL LINEAR INTERPOLATION----| h_cmplx in the 1/2 of the current time interval
   h_small_hist_cmplx(:,:,:,n_hist+1)=0.5_dp*(h_small_hist_cmplx(:,:,:,n_hist)+h_small_hist_cmplx(:,:,:,n_hist+2))

  ! if ( is_iomaster .AND. mod(itau-1,mod_write)==0 ) then
  !   write(name_iter_norm,"(3A,I4.4,A)") "./iter_norm/", TRIM(pred_corr), "_norm_itau_",itau,".dat"
  !   open(newunit=file_iter_norm,file=name_iter_norm)
  ! end if
   do i_iter=1,n_iter
     h_small_cmplx(:,:,:)=h_small_hist_cmplx(:,:,:,n_hist+1)
     !--3--PREDICTOR (propagation of C(0)-->C(1))
     call propagate_orth(nstate,basis,time_step,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx,h_small_cmplx,x_matrix,prop_type)

     !--4--CORRECTOR----| C(1)-->H(1)
     call setup_hamiltonian_cmplx(basis,                   &
                                  nstate,                  &
                                  itau,                    &
                                  time_cur,                &
                                  time_step,               &
                                  occupation,              &
                                  c_matrix_cmplx,          &
                                  hamiltonian_kinetic,     &
                                  hamiltonian_nucleus,     &
                                  h_small_hist_cmplx(:,:,:,n_hist+2),           &
                                  x_matrix,                &
                                  dipole_ao,               &
                                  h_cmplx,en_tddft)

     c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
     !--2B--LOCAL LINEAR INTERPOLATION----| h_cmplx in the 1/2 of the current time interval
     h_small_hist_cmplx(:,:,:,n_hist+1)=0.5_dp*(h_small_hist_cmplx(:,:,:,n_hist)+h_small_hist_cmplx(:,:,:,n_hist+2))

     !**COMPARISON**
   !  if( is_iomaster ) then
   !    if(mod(itau-1,mod_write)==0 ) then
   !      write(file_iter_norm,*) i_iter, NORM2(ABS(h_small_hist_cmplx(:,:,:,n_hist+1)-h_small_cmplx(:,:,:)))
   !    end if
   !  end if

   end do ! i_iter

!    close(file_iter_norm)

   !--5--PROPAGATION----| C(0)---U[H(1/2)]--->C(1)
   call propagate_orth(nstate,basis,time_step,c_matrix_orth_hist_cmplx(:,:,:,1), &
                       c_matrix_cmplx,h_small_hist_cmplx(:,:,:,n_hist+1),x_matrix,prop_type)

   !--6--UPDATE----|C(1)-->C(0)
   c_matrix_orth_cmplx(:,:,:)=c_matrix_orth_hist_cmplx(:,:,:,1)
   do iextr=1,n_hist-1
     h_small_hist_cmplx(:,:,:,iextr)=h_small_hist_cmplx(:,:,:,iextr+1)
   end do
   h_small_hist_cmplx(:,:,:,n_hist)=h_small_hist_cmplx(:,:,:,n_hist+2)


 ! ///////////////////////////////////
 ! Iterative ETRS - enforced time-reversal symmetry
 ! ---------------|t-dt|-----------------|t|--------------------|t+dt|
 !............|H(n_hist-1)|..........|H(n_hist)|.............|H(n_hist+1)|
 case('PC5','PC6')
   h_cmplx=(0.0_dp,0.0_dp)
   h_small_hist_cmplx(:,:,:,n_hist+1)=(0.0_dp,0.0_dp)
   !--1--EXTRAPOLATION WITH PREVIOUS STEPS----
   do iextr=1,n_hist
     h_small_hist_cmplx(:,:,:,n_hist+1)=h_small_hist_cmplx(:,:,:,n_hist+1)+extrap_coefs(iextr)*h_small_hist_cmplx(:,:,:,iextr)
   end do

   do i_iter=1,n_iter
     c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
     h_small_cmplx(:,:,:)=h_small_hist_cmplx(:,:,:,n_hist+1)
     !--3--PREDICTOR (propagation of C(0)-->C(1))
     call propagate_orth(nstate,basis,time_step,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx, &
                         h_small_hist_cmplx(:,:,:,n_hist:n_hist+1),x_matrix,"ETRS")
     !--4--CORRECTOR----| C(1)-->H(1)
     call setup_hamiltonian_cmplx(basis,                   &
                                  nstate,                  &
                                  itau,                    &
                                  time_cur,                &
                                  time_step,               &
                                  occupation,              &
                                  c_matrix_cmplx,          &
                                  hamiltonian_kinetic,     &
                                  hamiltonian_nucleus,     &
                                  h_small_hist_cmplx(:,:,:,n_hist+1) ,    &
                                  x_matrix,                &
                                  dipole_ao,               &
                                  h_cmplx,en_tddft)

      c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)

     !**COMPARISON**
   end do ! i_iter

    !--5--PROPAGATION----| C(0)---U[H(1/2)]--->C(!1)
    call propagate_orth(nstate,basis,time_step,c_matrix_orth_cmplx,c_matrix_cmplx, &
                        h_small_hist_cmplx(:,:,:,n_hist:n_hist+1),x_matrix,"ETRS")

   !--6--UPDATE----|
    c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
   c_matrix_orth_cmplx(:,:,:)=c_matrix_orth_hist_cmplx(:,:,:,1)
   do iextr=1,n_hist
     h_small_hist_cmplx(:,:,:,iextr)=h_small_hist_cmplx(:,:,:,iextr+1)
   end do


 ! ///////////////////////////////////
 ! ETRS proposed by Xavier Andrade
 ! ---------------|t-dt|-----------------|t|--------------------|t+dt|
 !............|---------|...............|H(1)|..................|H(2)|
 case('PC7' )
   h_cmplx=(0.0_dp,0.0_dp)
   h_small_hist_cmplx(:,:,:,2)=(0.0_dp,0.0_dp)

   call propagate_orth(nstate,basis,time_step,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx, &
                       h_small_hist_cmplx(:,:,:,1),x_matrix,"MAG2")

   do i_iter=1,n_iter
     call setup_hamiltonian_cmplx(basis,                   &
                                  nstate,                  &
                                  itau,                    &
                                  time_cur,                &
                                  time_step,               &
                                  occupation,              &
                                  c_matrix_cmplx,          &
                                  hamiltonian_kinetic,     &
                                  hamiltonian_nucleus,     &
                                  h_small_hist_cmplx(:,:,:,2),   &
                                  x_matrix,                &
                                  dipole_ao,               &
                                  h_cmplx,en_tddft)

     if(i_iter/=n_iter) then
       c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
       call propagate_orth(nstate,basis,time_step,c_matrix_orth_hist_cmplx(:,:,:,1),c_matrix_cmplx, &
                           h_small_hist_cmplx(:,:,:,1:2),x_matrix,"ETRS")
     end if
   end do

   !--6--UPDATE----|C(1)-->C(0)
   c_matrix_orth_cmplx(:,:,:)=c_matrix_orth_hist_cmplx(:,:,:,1)
   h_small_hist_cmplx(:,:,:,1)=h_small_hist_cmplx(:,:,:,2)

 end select

 write(stdout,'(/,1x,a)') 'END OF PREDICTOR-CORRECTOR BLOCK'

end subroutine predictor_corrector


!=========================================================================
subroutine initialize_extrap_coefs(c_matrix_orth_cmplx,h_small_cmplx)
 implicit none
 complex(dp),intent(in)    :: c_matrix_orth_cmplx(:,:,:)
 complex(dp),intent(in)    :: h_small_cmplx(:,:,:)
!=====
 integer               :: iextr,ham_hist_dim,nstate
 real(dp)              :: x_pred
 real(dp),allocatable  :: m_nodes(:)
!=====

 nstate = SIZE(c_matrix_orth_cmplx,DIM=1)

 allocate(m_nodes(n_hist),extrap_coefs(n_hist))

 select case (pred_corr)
 case('PC1')
   ham_hist_dim = 2

 case('PC2B')
   ham_hist_dim = n_hist
   do iextr=1,n_hist
     m_nodes(iextr) = (iextr - 1.0_dp) * 0.5_dp
   end do
   x_pred = (n_hist - 1.0_dp) * 0.5_dp + 0.25_dp
   call get_extrap_coefs_lagr(m_nodes,x_pred,extrap_coefs,n_hist)

 case('PC3','PC4')
   ham_hist_dim = n_hist + 2
   do iextr=1,n_hist
     m_nodes(iextr) = iextr - 1.0_dp
   end do
   x_pred = n_hist
   if(pred_corr=='PC3') call get_extrap_coefs_lagr(m_nodes,x_pred,extrap_coefs,n_hist)
   if(pred_corr=='PC4') call get_extrap_coefs_aspc(extrap_coefs,n_hist)

 case('PC5','PC6')
   ham_hist_dim = n_hist + 1
   do iextr=1,n_hist
     m_nodes(iextr) = iextr - 1.0_dp
   end do
   x_pred = n_hist
   if(pred_corr=='PC5') call get_extrap_coefs_lagr(m_nodes,x_pred,extrap_coefs,n_hist)
   if(pred_corr=='PC6') call get_extrap_coefs_aspc(extrap_coefs,n_hist)

 case('PC7' )
   ham_hist_dim = 2

 end select

 if( pred_corr /= 'PC0' ) then
   call clean_allocate('h_small_hist_cmplx for TDDFT',h_small_hist_cmplx,nstate,nstate,nspin,ham_hist_dim)
   call clean_allocate('c_matrix_orth_hist_cmplx for TDDFT',c_matrix_orth_hist_cmplx,nstate,nocc,nspin,1)
   do iextr=1,ham_hist_dim
     h_small_hist_cmplx(:,:,:,iextr)=h_small_cmplx(:,:,:)
   end do
   c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
 end if

 deallocate(m_nodes)

end subroutine initialize_extrap_coefs


!=========================================================================
subroutine print_tddft_values(time_cur,file_time_data,file_dipole_time,file_excit_field,itau)
 implicit none
 integer,intent(in)    :: file_time_data,file_dipole_time,file_excit_field
 real(dp),intent(in)   :: time_cur
 integer,intent(in)    :: itau
!=====

 if( .NOT. is_iomaster ) return

 write(stdout,'(/,1x,a)') &
    '==================================================================================================='
 write(stdout,'(1x,a,i8,a)') &
    '===================== RT-TDDFT values for the iteration  ',itau,' ================================='
 write(stdout,'(a31,1x,f19.10)') 'RT-TDDFT Simulation time (au):', time_cur
 write(stdout,'(a31,1x,f19.10)') 'RT-TDDFT Total Energy    (Ha):', en_tddft%total

 select case(excit_type%form)
 case(EXCIT_PROJECTILE)
   write(file_time_data,"(f10.4,11(2x,es16.8e3))") &
      time_cur, en_tddft%total, xatom(:,ncenter_nuclei), en_tddft%nuc_nuc, en_tddft%nucleus, &
      en_tddft%kinetic, en_tddft%hartree, en_tddft%exx_hyb, en_tddft%xc, &
      en_tddft%excit
   call output_projectile_position()

 case(EXCIT_LIGHT)
   write(file_time_data,"(f10.4,8(2x,es16.8e3))") &
    time_cur, en_tddft%total, en_tddft%nuc_nuc, en_tddft%nucleus, en_tddft%kinetic, &
    en_tddft%hartree, en_tddft%exx_hyb, en_tddft%xc, en_tddft%excit
   write(file_dipole_time,'(4f19.10)') time_cur, dipole(:) * au_debye
   write(file_excit_field,'(2f19.10)') time_cur, excit_field_norm
   write(stdout,'(a31,1x,3f19.10)') 'RT-TDDFT Dipole Moment    (D):', dipole(:) * au_debye
 end select

 write(stdout,'(1x,a,/)') '==================================================================================================='

end subroutine print_tddft_values


!=========================================================================
subroutine initialize_files(file_time_data,file_dipole_time,file_excit_field)
 implicit none
 integer,intent(inout)    :: file_time_data,file_excit_field,file_dipole_time
!=====

 if( .NOT. is_iomaster ) return

 open(newunit=file_time_data,file="time_data.dat")

 if(excit_type%form==EXCIT_LIGHT) then

   open(newunit=file_dipole_time,file="dipole_time.dat")
   open(newunit=file_excit_field,file="excitation_time.dat")

   write(file_excit_field,"(A)") "# time(au)                      E_field_excit_dir(au)"

 end if

!---------------------------------
 select case(excit_type%form)
 case(EXCIT_PROJECTILE)
   write(file_time_data,"(a10,11(a18))") &
           "# time(au)","e_total","x_proj","y_proj","z_proj","enuc_nuc","enuc_wo_proj","ekin","ehart",&
           "eexx_hyb","exc","enuc_proj"

 case(EXCIT_LIGHT)
   write(file_time_data,"(a,a)") &
     " # time(au)     e_total             enuc_nuc             enuc            ekin               ehart   ", &
    "     eexx_hyb            exc             eexcit"

   write(file_dipole_time,"(a)") &
      "# time(au)                      Dipole_x(D)               Dipole_y(D)               Dipole_z(D)"
 end select

end subroutine initialize_files


!=========================================================================
subroutine initialize_q(nstate,nocc,nspin,c_matrix_orth_start_complete_cmplx,h_small_cmplx,istate_cut,file_q_matrix)
 implicit none
 integer,intent(in)                    :: nstate, nocc, nspin
 integer,allocatable,intent(out)       :: istate_cut(:,:)
 complex(dp),allocatable,intent(out)   :: c_matrix_orth_start_complete_cmplx(:,:,:)
 complex(dp),allocatable,intent(in)    :: h_small_cmplx(:,:,:)
 integer,intent(out)                   :: file_q_matrix(2)
!=====
 character(len=50)          :: name_file_q_matrix
 character(len=500)         :: cur_string, header_format
 character(len=500)         :: header
 integer                    :: ispin
 real(dp),allocatable       :: energy_tddft(:)
 logical                    :: file_exists
 integer                    :: file_q_matrix_param,ncut,icut
 integer                    :: num_fields
!=====

 call clean_allocate('c_matrix_orth_start for TDDFT',c_matrix_orth_start_complete_cmplx,nstate,nstate,nspin)
 allocate(energy_tddft(nstate))
 do ispin=1, nspin
   call diagonalize(postscf_diago_flavor,h_small_cmplx(:,:,ispin),energy_tddft,c_matrix_orth_start_complete_cmplx(:,:,ispin))
 end do
 deallocate(energy_tddft)

 inquire(file='manual_q_matrix_param',exist=file_exists)

 if(file_exists) then
   ncut=get_number_of_lines('manual_q_matrix_param')
   allocate(istate_cut(ncut,2))
   open(newunit=file_q_matrix_param,file='manual_q_matrix_param',status='old')
   do icut=1,ncut
     read(file_q_matrix_param,'(A)') cur_string
!     cur_string = ADJUSTL(cur_string)
     num_fields = get_number_of_elements(cur_string)
     if( num_fields == 2 ) then
       read(cur_string,*) istate_cut(icut,1), istate_cut(icut,2)
     else if( num_fields == 1) then
       read(cur_string,*) istate_cut(icut,1)
       istate_cut(icut,2) = nstate
     else
       call die("manual_q_matrix_param must contain one or two fields.")
     end if

   end do
   close(file_q_matrix_param)
 else
   allocate(istate_cut(2,2))
   istate_cut(1,1) = 1; istate_cut(1,2) = 1;
   istate_cut(2,1) = 2; istate_cut(2,2) = nstate;
   call issue_warning('initialize_q: manual_q_matrix_param file was not found')
 endif

! Header for the q_matrix file
 write(header,"(A10,3(A18))") "# time(au)","x_proj","y_proj","z_proj"

 do icut=1,ncut
   write(header,'(a,9x,i4.4,a1,i4.4)') TRIM(header),istate_cut(icut,1),'-',istate_cut(icut,2)
 end do

 if( is_iomaster ) then
   do ispin=1,nspin
     write(name_file_q_matrix,"(a,i1,a)") "q_matrix_", ispin, ".dat"
     open(newunit=file_q_matrix(ispin),file=name_file_q_matrix)

     write(file_q_matrix(ispin),'(a)') TRIM(header)

   end do
 end if

end subroutine initialize_q


!=========================================================================
subroutine calculate_q_matrix(occupation,c_matrix_orth_start_complete_cmplx,c_matrix_orth_cmplx,istate_cut,file_q_matrix,time_cur)
 implicit none
 real(dp),intent(in)      :: occupation(:,:)
 complex(dp),intent(in)   :: c_matrix_orth_start_complete_cmplx(:,:,:)
 complex(dp),intent(in)   :: c_matrix_orth_cmplx(:,:,:)
 integer,intent(in)       :: istate_cut(:,:)
 integer,intent(in)       :: file_q_matrix(:)
 real(dp),intent(in)      :: time_cur
!=====
 integer                  :: istate,jstate,iocc,ispin,icut,ncut,nstate
 complex(dp),allocatable  :: q_matrix_cmplx(:,:)
 real(dp),allocatable     :: q_occ(:)
!=====

 call start_clock(timing_tddft_q_matrix)

 nstate = SIZE(occupation(:,:),DIM=1)
 ncut = SIZE(istate_cut,DIM=1)

 allocate(q_occ(ncut))
 call clean_allocate('q_matrix for TDDFT',q_matrix_cmplx,nstate,nocc)

 q_occ = 0.0_dp

 do ispin=1,nspin

   !q_matrix_cmplx(:,:)= &
   !MATMUL(CONJG(TRANSPOSE(c_matrix_orth_start_complete_cmplx(:,:,ispin)), &
   !                          c_matrix_orth_cmplx(:,:,ispin))

   call ZGEMM('C','N',nstate,nocc,nstate,(1.0d0,0.0d0),c_matrix_orth_start_complete_cmplx(:,:,ispin),nstate, &
                                            c_matrix_orth_cmplx(:,:,ispin),nstate, &
                                         (0.0d0,0.0d0),q_matrix_cmplx,nstate)


   do icut=1,ncut
     do istate=istate_cut(icut,1),istate_cut(icut,2)
       do iocc=1,nocc
         q_occ(icut) = q_occ(icut) + ABS(q_matrix_cmplx(istate,iocc))**2*occupation(iocc,ispin)
       end do
     end do
   end do

   if( is_iomaster) then
     write(file_q_matrix(ispin),"(F10.4,30(2x,es16.8E3))") time_cur, xatom(:,ncenter_nuclei), q_occ(:)
   end if
 end do

 call clean_deallocate('q_matrix for TDDFT',q_matrix_cmplx)

 call stop_clock(timing_tddft_q_matrix)

end subroutine calculate_q_matrix


!=========================================================================
subroutine check_identity_cmplx(n,m,mat_cmplx,ans)
 implicit none
 integer,intent(in)     :: n, m
 complex(dp),intent(in) :: mat_cmplx(n,m)
 logical,intent(inout)  :: ans
!=====
 integer   ::  imat,jmat
 real(dp),parameter :: tol=1.0e-9_dp
!=====

 ans=.TRUE.
 do imat=1,n
   do jmat=1,m
     if(imat==jmat) then
       if(ABS(mat_cmplx(imat,jmat)-1.0_dp)>tol) then
         write(stdout,*) "M(imat,imat)/=1 for: ",imat,jmat,mat_cmplx(imat,jmat)
         ans=.FALSE.
       end if
     else
       if(ABS(mat_cmplx(imat,jmat))>tol) then
         write(stdout,*) "M(imat,jmat)/=0 for: ",imat,jmat,mat_cmplx(imat,jmat)
         ans=.FALSE.
       end if
     end if
   end do
 end do

end subroutine check_identity_cmplx


!=========================================================================
subroutine write_restart_tddft(nstate,time_cur,occupation,c_matrix_orth_cmplx)
 implicit none
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: time_cur
 real(dp),intent(in)        :: occupation(:,:)
 complex(dp),intent(in)     :: c_matrix_orth_cmplx(nstate,nocc,nspin)
!===
 integer                    :: restartfile
 integer                    :: istate,ispin
!=====

 if( .NOT. is_iomaster) return

 call start_clock(timing_restart_tddft_file)

 write(stdout,'(/,a,f19.10)') ' Writing a RESTART_TDDFT file, time_cur= ', time_cur

 open(newunit=restartfile,file='RESTART_TDDFT',form='unformatted',action='write')
 ! nspin
 write(restartfile) nspin
 ! Nstate
 write(restartfile) nstate
 ! Occupations
 write(restartfile) occupation(:,:)
 ! Current Time
 write(restartfile) time_cur
 ! Atomic structure
 write(restartfile) ncenter_nuclei
 write(restartfile) zatom(1:ncenter_nuclei)
 write(restartfile) xatom(:,1:ncenter_nuclei)
 ! Complex wavefunction coefficients C
 do ispin=1,nspin
   do istate=1,nocc
     write(restartfile) c_matrix_orth_cmplx(:,istate,ispin)
   enddo
 enddo

 close(restartfile)

 call stop_clock(timing_restart_tddft_file)

end subroutine write_restart_tddft


!=========================================================================
subroutine check_restart_tddft(nstate,occupation,restart_is_correct)
 implicit none
 logical,intent(out)        :: restart_is_correct
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: occupation(nstate,nspin)
!===
 logical                    :: file_exists
 integer                    :: restartfile
 integer                    :: istate,ispin
 integer                    :: natom_read
 real(dp)                   :: time_cur_read
 real(dp),allocatable       :: occupation_read(:,:)
 real(dp),allocatable       :: zatom_read(:),x_read(:,:)
 integer                    :: nstate_read, nspin_read
!=====

 write(stdout,'(/,a)') ' Checking RESTART_TDDFT file'

 restart_is_correct=.TRUE.

 inquire(file='RESTART_TDDFT',exist=file_exists)
 if(.NOT. file_exists) then
   write(stdout,'(/,a)') ' No RESTART_TDDFT file found'
   restart_is_correct=.FALSE.
   return
 endif

 open(newunit=restartfile,file='RESTART_TDDFT',form='unformatted',status='old',action='read')

 ! Nspin
 read(restartfile) nspin_read
 if(nspin /= nspin_read) then
   call issue_warning('RESTART_TDDFT file: nspin is not the same, restart file will not be used')
   restart_is_correct=.FALSE.
   close(restartfile)
   return
 end if

 ! Nstate
 read(restartfile) nstate_read
 if(nstate /= nstate_read) then
   call issue_warning('RESTART_TDDFT file: nstate is not the same, restart file will not be used')
   restart_is_correct=.FALSE.
   close(restartfile)
   return
 end if

 ! Occupations
 allocate(occupation_read(nstate,nspin))
 read(restartfile) occupation_read(:,:)
 if( (ANY( ABS( occupation_read(:,:) - occupation(:,:) )  > 1.0e-5_dp )) .AND. temperature > 1.0e-8_dp ) then
   call issue_warning('RESTART file: Occupations have changed')
 endif
 deallocate(occupation_read)

 ! Current time
 read(restartfile) time_cur_read

 !Different number of atoms in restart and input files is not provided for tddft restart
 read(restartfile) natom_read
 if( natom_read /= ncenter_nuclei ) then
   call issue_warning('RESTART_TDDFT file: number of atoms has changed.')
   restart_is_correct=.FALSE.
   close(restartfile)
   return
 end if

 allocate(zatom_read(natom_read),x_read(3,natom_read))
 read(restartfile) zatom_read(1:natom_read)
 read(restartfile) x_read(:,1:natom_read)
 if( natom_read /= ncenter_nuclei  &
  .OR. ANY( ABS( zatom_read(1:MIN(natom_read,ncenter_nuclei)) - zatom(1:MIN(natom_read,ncenter_nuclei)) ) > 1.0e-5_dp ) &
  .OR. ANY( ABS(   x_read(:,1:MIN(natom_read,ncenter_nuclei-nprojectile)) &
                   - xatom(:,1:MIN(natom_read,ncenter_nuclei-nprojectile)) ) &
            > 1.0e-5_dp ) ) then
   call issue_warning('RESTART_TDDFT file: Geometry has changed')
 endif
 deallocate(zatom_read,x_read)

 ! Here we do not reead c_matrix_orth_cmplx from the RESTART_TDDFT file

 write(stdout,*) " RESTART_TDDFT file is correct and will be used for the calculation. SCF loop will be omitted."

 close(restartfile)

end subroutine check_restart_tddft


!=========================================================================
subroutine read_restart_tddft(nstate,time_min,occupation,c_matrix_orth_cmplx)
 implicit none
 complex(dp),intent(inout)  :: c_matrix_orth_cmplx(nstate,nocc,nspin)
 real(dp),intent(inout)     :: time_min
 real(dp),intent(inout)     :: occupation(nstate,nspin)
 integer,intent(in)         :: nstate
!===
 logical                    :: file_exists
 integer                    :: restartfile
 integer                    :: istate,ispin
 integer                    :: natom_read
 real(dp),allocatable       :: occupation_read(:,:)
 real(dp),allocatable       :: zatom_read(:),x_read(:,:)
 integer                    :: nstate_read, nspin_read
!=====

 write(stdout,'(/,a)') ' Reading a RESTART_TDDFT file'

 inquire(file='RESTART_TDDFT',exist=file_exists)
 if(.NOT. file_exists) then
   write(stdout,'(/,a)') ' No RESTART_TDDFT file found'
   return
 endif

 open(newunit=restartfile,file='RESTART_TDDFT',form='unformatted',status='old',action='read')

 ! After the subroutine check_restart_tddft, which was called in the molgw.f90, here we are sure that:
 ! nspin==nspin_read
 ! nstate==nstate_read

 ! Nspin
 read(restartfile) nspin_read

 ! Nstate
 read(restartfile) nstate_read

 ! Occupations
 allocate(occupation_read(nstate_read,nspin_read))
 read(restartfile) occupation_read(:,:)
 if( ANY( ABS( occupation_read(:,:) - occupation(:,:) )  > 1.0e-5_dp ) ) then
   if( temperature > 1.0e-8_dp) then
     occupation(:,:) = occupation_read(:,:)
     write(stdout,'(1xa)') "Reading occupations from a RESTART file"
     call dump_out_occupation('=== Occupations ===',occupation)
   else
     call issue_warning('RESTART file: Occupations have changed')
   endif
 endif
 deallocate(occupation_read)

 ! Current Time
 read(restartfile) time_min
 write(stdout,"(1x,a,f7.3)") "time_min= ", time_min

 ! Atomic structure
 read(restartfile) natom_read
 allocate(zatom_read(natom_read),x_read(3,natom_read))
 read(restartfile) zatom_read(1:natom_read)
 read(restartfile) xatom_start(:,1:natom_read)
 deallocate(zatom_read)

 ! Complex wavefunction coefficients C
 do ispin=1,nspin
   do istate=1,nocc
     read(restartfile) c_matrix_orth_cmplx(:,istate,ispin)
   enddo
 enddo

 close(restartfile)

end subroutine read_restart_tddft


!=========================================================================
function get_nocc_from_restart() result(nocc_read)
 implicit none
!=====
 integer                    :: nocc_read
!=====
 logical                    :: file_exists
 integer                    :: restartfile
 integer                    :: nstate_read, nspin_read
 real(dp),allocatable       :: occupation_read(:,:)
!=====

 write(stdout,'(/,a)') ' Reading a RESTART_TDDFT file to get the nocc value'

 inquire(file='RESTART_TDDFT',exist=file_exists)
 if(.NOT. file_exists) then
   write(stdout,'(/,a)') ' No RESTART_TDDFT file found'
   return
 endif

 open(newunit=restartfile,file='RESTART_TDDFT',form='unformatted',status='old',action='read')

 ! Nspin
 read(restartfile) nspin_read

 ! Nstate
 read(restartfile) nstate_read

 ! Occupations
 allocate(occupation_read(nstate_read,nspin_read))
 read(restartfile) occupation_read(:,:)

 nocc_read = get_number_occupied_states(occupation_read)

 deallocate(occupation_read)

 close(restartfile)

end function get_nocc_from_restart


!=========================================================================
subroutine get_extrap_coefs_lagr(m_nodes,x_pred,extrap_coefs,n_hist_cur)
 implicit none
 integer, intent(in)          :: n_hist_cur
 real(dp),intent(in)          :: m_nodes(n_hist_cur)
 real(dp),intent(in)          :: x_pred
 real(dp),intent(inout)       :: extrap_coefs(n_hist_cur)
!=====
 integer      :: inode,jnode
!=====

 extrap_coefs(:)=1.0_dp
 do inode=1,n_hist_cur
   do jnode=1,n_hist_cur
     if(inode==jnode) cycle
     extrap_coefs(inode)=extrap_coefs(inode)*(x_pred-m_nodes(jnode))/(m_nodes(inode)-m_nodes(jnode))
   end do
 end do
end subroutine get_extrap_coefs_lagr


!=========================================================================
subroutine get_extrap_coefs_aspc(extrap_coefs,n_hist_cur)
 implicit none
 integer, intent(in)          :: n_hist_cur!
 real(dp),intent(inout)       :: extrap_coefs(n_hist_cur)
!=====

 if(n_hist_cur==1) then
   extrap_coefs(1)=1.0_dp
 end if


 if(n_hist_cur==2) then
   extrap_coefs(1)=-1.0_dp
   extrap_coefs(2)=2.0_dp
 end if

 if(n_hist_cur==3) then
   extrap_coefs(1)=0.5_dp
   extrap_coefs(2)=-2.0_dp
   extrap_coefs(3)=2.5_dp
 end if

 if(n_hist_cur==4) then
   extrap_coefs(1)=-0.2_dp
   extrap_coefs(2)=1.2_dp
   extrap_coefs(3)=-2.8_dp
   extrap_coefs(4)=2.8_dp
 end if

 if(n_hist_cur==5) then
   extrap_coefs(1)=1.0_dp/14.0_dp
   extrap_coefs(2)=-4.0_dp/7.0_dp
   extrap_coefs(3)=27.0_dp/14.0_dp
   extrap_coefs(4)=-24.0_dp/7.0_dp
   extrap_coefs(5)=3.0_dp
 end if

 if(n_hist_cur==6) then
   extrap_coefs(1)=-1.0_dp/42.0_dp
   extrap_coefs(2)=5.0_dp/21.0_dp
   extrap_coefs(3)=-22.0_dp/21.0_dp
   extrap_coefs(4)=55.0_dp/21.0_dp
   extrap_coefs(5)=-55.0_dp/14.0_dp
   extrap_coefs(6)=22.0_dp/7.0_dp
 end if


end subroutine get_extrap_coefs_aspc


!=========================================================================
subroutine propagate_orth_ham_1(nstate,basis,time_step_cur,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_cmplx,x_matrix,prop_type)
 implicit none
 integer,intent(in)          :: nstate
 type(basis_set),intent(in)  :: basis
 real(dp),intent(in)         :: time_step_cur
 complex(dp),intent(inout)   :: c_matrix_orth_cmplx(nstate,nocc,nspin)
 complex(dp), intent(inout)  :: c_matrix_cmplx(basis%nbf,nocc,nspin)
 ! complex(dp),intent(in)      :: h_small_cmplx(nstate,nstate,nspin)
 ! for ncore>0 we need to modify the h_small
 complex(dp),intent(inout)   :: h_small_cmplx(nstate,nstate,nspin)
 real(dp),intent(in)         :: x_matrix(basis%nbf,nstate)
 character(len=4),intent(in) :: prop_type
!=====
 integer                    :: ispin
 integer                    :: istate,jstate
 complex(dp),allocatable    :: m_tmp_1(:,:)
 complex(dp),allocatable    :: m_tmp_2(:,:)
 complex(dp),allocatable    :: m_tmp_3(:,:)
 real(dp),allocatable       :: m_tmpr1(:,:),m_tmpr2(:,:)
!==variables for the MAG2 propagator
 complex(dp),allocatable    :: a_matrix_orth_cmplx(:,:)
 real(dp)                   :: energy_tddft(nstate)
!==variables for the CN propagator
 complex(dp),allocatable    :: l_matrix_cmplx(:,:) ! Follow the notation of M.A.L.Marques, C.A.Ullrich et al,
 complex(dp),allocatable    :: b_matrix_cmplx(:,:) ! TDDFT Book, Springer (2006), !p205
!==frozen core==
 complex(dp),allocatable    :: h_initial_basis_orth_cmplx(:,:)
 complex(dp),allocatable    :: tmp_conjg_transpose(:,:)
!=====

 call start_clock(timing_tddft_propagation)

 do ispin=1,nspin
   select case(prop_type)
   case('CN')
     allocate(l_matrix_cmplx(nstate,nstate))
     allocate(b_matrix_cmplx(nstate,nstate))
     l_matrix_cmplx(:,:) =  im * time_step_cur / 2.0_dp * h_small_cmplx(:,:,ispin)
     b_matrix_cmplx(:,:) = -l_matrix_cmplx(:,:)
     do istate=1,nstate
       b_matrix_cmplx(istate,istate) = b_matrix_cmplx(istate,istate) + 1.0_dp
       l_matrix_cmplx(istate,istate) = l_matrix_cmplx(istate,istate) + 1.0_dp
     end do
     call invert(l_matrix_cmplx)

     call start_clock(timing_propagate_matmul)
     b_matrix_cmplx(:,:)            = MATMUL( l_matrix_cmplx(:,:),b_matrix_cmplx(:,:))
     c_matrix_orth_cmplx(:,:,ispin) = MATMUL( b_matrix_cmplx(:,:),c_matrix_orth_cmplx(:,:,ispin))

!    c_matrix_orth_cmplx(:,:,ispin) = MATMUL( l_matrix_cmplx(:,:),MATMUL( b_matrix_cmplx(:,:),c_matrix_orth_cmplx(:,:,ispin) ) )
     deallocate(l_matrix_cmplx)
     deallocate(b_matrix_cmplx)
   case('MAG2')
     allocate(a_matrix_orth_cmplx(nstate,nstate))

     !
     ! Frozen core
     if(ncore_tddft > 0) then

       call start_clock(timing_tddft_frozen_core)

       call clean_allocate('h_initial_basis_orth_cmplx for the frozen core',h_initial_basis_orth_cmplx,nstate,nstate)
       call clean_allocate('tmp_conjg_transpose for the frozen core',tmp_conjg_transpose,nstate,nstate)

       ! Express the h_small in the a_matrix_orth_start_cmplx basis
       !h_small_cmplx(:,:,ispin) = MATMUL( CONJG( TRANSPOSE(a_matrix_orth_start_cmplx(:,:,ispin))), MATMUL( h_small_cmplx(:,:,ispin), a_matrix_orth_start_cmplx(:,:,ispin) )  )
       tmp_conjg_transpose = CONJG( TRANSPOSE(a_matrix_orth_start_cmplx(:,:,ispin)))
       call matmul_abc_scalapack(scalapack_block_min,tmp_conjg_transpose, &
                                 h_small_cmplx(:,:,ispin),a_matrix_orth_start_cmplx(:,:,ispin),h_initial_basis_orth_cmplx  )

       ! Modify the h_small in the a_matrix_orth_start_cmplx basis
       do istate=1,nstate
         do jstate=1,nstate
           if( istate > ncore_tddft .AND. jstate > ncore_tddft ) cycle
           h_initial_basis_orth_cmplx(istate,jstate) = 0.d0
         end do
       end do

       ! Return in the _orth_ basis
       !h_small_cmplx(:,:,ispin) = MATMUL( a_matrix_orth_start_cmplx(:,:,ispin), MATMUL( h_small_cmplx(:,:,ispin),CONJG( TRANSPOSE(a_matrix_orth_start_cmplx(:,:,ispin))) )  )
       call matmul_abc_scalapack(scalapack_block_min,a_matrix_orth_start_cmplx(:,:,ispin), &
                                 h_initial_basis_orth_cmplx(:,:),tmp_conjg_transpose,h_small_cmplx(:,:,ispin)  )

       call clean_deallocate('h_initial_basis_orth_cmplx for the frozen core',h_initial_basis_orth_cmplx)
       call clean_deallocate('tmp_conjg_transpose for the frozen core',tmp_conjg_transpose)

       call stop_clock(timing_tddft_frozen_core)

     end if

     !
     ! First part, diagonalize
     ! subroutine is able to take advantage of a real hamiltonian
     call diagonalize_hamiltonian_ortho(h_small_cmplx(:,:,ispin),a_matrix_orth_cmplx,energy_tddft)

     !
     ! Second part, multiply matrices
     call start_clock(timing_propagate_matmul)

     allocate(m_tmp_1(nstate,nstate))
     allocate(m_tmp_2(nstate,nstate))
     ! M1 = A * e^{-idt*e}
     do jstate=1,nstate
       m_tmp_1(:,jstate) = a_matrix_orth_cmplx(:,jstate) * EXP( -im * time_step_cur * energy_tddft(jstate) )
     enddo

     ! M2 = M1 * A**H = (A * e^{-idt*e} ) * A**H
     call ZGEMM('N','C',nstate,nstate,nstate,(1.0d0,0.0d0),m_tmp_1,nstate,             &
                                                           a_matrix_orth_cmplx,nstate, &
                                             (0.0d0,0.0d0),m_tmp_2,nstate)
     deallocate(m_tmp_1)

     allocate(m_tmp_3(nstate,nocc))
     m_tmp_3(:,:) = c_matrix_orth_cmplx(:,:,ispin)

     ! C'^new =  M2 * C'^old = (A * e^{-idt*e} ) * A**H * C'^old
     call ZGEMM('N','N',nstate,nocc,nstate,(1.0d0,0.0d0),m_tmp_2,nstate, &
                                                         m_tmp_3,nstate, &
                                           (0.0d0,0.0d0),c_matrix_orth_cmplx(1,1,ispin),nstate)

     deallocate(m_tmp_3)
     deallocate(m_tmp_2)
     deallocate(a_matrix_orth_cmplx)

     case default
       call die('Invalid choice for the propagation algorithm. Change prop_type or error_prop_types value in the input file')
   end select

   !
   ! C = X * C'
   ! x_matrix is real, let's use this to save time
   allocate(m_tmpr1(nstate,nocc),m_tmpr2(basis%nbf,nocc))
   ! 1. Real part
   m_tmpr1(:,:) = c_matrix_orth_cmplx(:,:,ispin)%re
   call DGEMM('N','N',basis%nbf,nocc,nstate,1.0d0,x_matrix(1,1),basis%nbf, &
                                                  m_tmpr1(1,1),nstate,   &
                                            0.0d0,m_tmpr2(1,1),basis%nbf)
   c_matrix_cmplx(:,:,ispin)%re = m_tmpr2(:,:)
   ! 2. Imaginary part
   m_tmpr1(:,:) = c_matrix_orth_cmplx(:,:,ispin)%im
   call DGEMM('N','N',basis%nbf,nocc,nstate,1.0d0,x_matrix(1,1),basis%nbf, &
                                                  m_tmpr1(1,1),nstate,   &
                                            0.0d0,m_tmpr2(1,1),basis%nbf)
   c_matrix_cmplx(:,:,ispin)%im = m_tmpr2(:,:)

   deallocate(m_tmpr1,m_tmpr2)


   call stop_clock(timing_propagate_matmul)

 end do

 call stop_clock(timing_tddft_propagation)


end subroutine propagate_orth_ham_1


!=========================================================================
subroutine propagate_orth_ham_2(nstate,basis,time_step_cur,c_matrix_orth_cmplx,c_matrix_cmplx, &
                                h_small_hist2_cmplx,x_matrix,prop_type)
 implicit none
 integer,intent(in)          :: nstate
 type(basis_set),intent(in)  :: basis
 real(dp),intent(in)         :: time_step_cur
 complex(dp),intent(inout)   :: c_matrix_orth_cmplx(nstate,nocc,nspin)
 complex(dp), intent(inout)  :: c_matrix_cmplx(basis%nbf,nocc,nspin)
 complex(dp),intent(in)      :: h_small_hist2_cmplx(nstate,nstate,nspin,2)
 real(dp),intent(in)         :: x_matrix(basis%nbf,nstate)
 character(len=4),intent(in) :: prop_type
!=====
 integer             :: ispin,iham
 integer             :: ibf
 complex(dp)         :: a_matrix_orth_cmplx(nstate,nstate,2)
 real(dp)            :: energy_tddft(nstate)
 complex(dp)         :: propagator_eigen(nstate,nstate,2)
!=====

 call start_clock(timing_tddft_propagation)
! a_matrix_cmplx(:,1:nstate) = MATMUL( x_matrix(:,:) , a_matrix_cmplx(:,:) )

 do ispin =1, nspin
   select case (prop_type)
   case('ETRS')
     do iham=1,2
       call diagonalize(postscf_diago_flavor,h_small_hist2_cmplx(:,:,ispin,iham),energy_tddft,a_matrix_orth_cmplx(:,:,iham))
       propagator_eigen(:,:,iham) = ( 0.0_dp , 0.0_dp )
       do ibf=1,nstate
         propagator_eigen(ibf,ibf,iham) = EXP(-im*time_step_cur/2.d0*energy_tddft(ibf))
       end do
     end do
     c_matrix_orth_cmplx(:,:,ispin) = &
         MATMUL(MATMUL(MATMUL(MATMUL( MATMUL( MATMUL( a_matrix_orth_cmplx(:,:,2), propagator_eigen(:,:,2)  ) , &
             CONJG(TRANSPOSE(a_matrix_orth_cmplx(:,:,2)))  ),  &
               a_matrix_orth_cmplx(:,:,1)), propagator_eigen(:,:,1)), &
                  CONJG(TRANSPOSE(a_matrix_orth_cmplx(:,:,1))) ), &
                    c_matrix_orth_cmplx(:,:,ispin) )
   case default
     call die('Invalid choice of the propagation algorithm for the given PC scheme. Change prop_type value in the input file')
   end select
    c_matrix_cmplx(:,:,ispin) = MATMUL( x_matrix(:,:) , c_matrix_orth_cmplx(:,:,ispin) )
 end do

call stop_clock(timing_tddft_propagation)

end subroutine propagate_orth_ham_2


!=========================================================================
subroutine setup_hamiltonian_cmplx(basis,                   &
                                   nstate,                  &
                                   itau,                    &
                                   time_cur,                &
                                   time_step_cur,           &
                                   occupation,              &
                                   c_matrix_cmplx,          &
                                   hamiltonian_kinetic,     &
                                   hamiltonian_nucleus,     &
                                   h_small_cmplx,           &
                                   x_matrix,                &
                                   dipole_ao,               &
                                   hamiltonian_cmplx,       &
                                   en)

 implicit none
!=====
 type(basis_set),intent(in)      :: basis
 integer,intent(in)              :: nstate
 integer,intent(in)              :: itau
 real(dp),intent(in)             :: time_cur
 real(dp),intent(in)             :: time_step_cur
 real(dp),intent(in)             :: occupation(nstate,nspin)
 real(dp),intent(in)             :: hamiltonian_kinetic(basis%nbf,basis%nbf)
 real(dp),intent(in)             :: hamiltonian_nucleus(basis%nbf,basis%nbf)
 real(dp),allocatable,intent(in) :: dipole_ao(:,:,:)
 real(dp),intent(in)             :: x_matrix(basis%nbf,nstate)
 complex(dp),intent(in)          :: c_matrix_cmplx(basis%nbf,nocc,nspin)
 complex(dp),intent(out)         :: hamiltonian_cmplx(basis%nbf,basis%nbf,nspin)
 complex(dp),intent(out)         :: h_small_cmplx(nstate,nstate,nspin)
 type(energy_contributions),intent(inout) :: en
!=====
 logical              :: calc_excit_
 integer              :: ispin, idir
 real(dp)             :: excit_field(3)
 complex(dp)          :: p_matrix_cmplx(basis%nbf,basis%nbf,nspin)
 integer              :: projectile_list(1)
 real(dp),allocatable :: hamiltonian_projectile(:,:)
!=====

 call start_clock(timing_tddft_hamiltonian)

 call setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)

 !--Hamiltonian - Hartree Exchange Correlation---
 call calculate_hamiltonian_hxc_ri_cmplx(basis,                    &
                                         occupation,               &
                                         c_matrix_cmplx,           &
                                         p_matrix_cmplx,           &
                                         hamiltonian_cmplx,en)



 !
 ! Excitation part of the Hamiltonian
 !
 en_tddft%excit = 0.0_dp

 select case(excit_type%form)
 !
 ! Light excitation
 case(EXCIT_LIGHT)
   excit_field = 0.0_dp
   calc_excit_ = .FALSE.
   calc_excit_ = calc_excit_ .OR. ( excit_type%name == 'GAU' )
   calc_excit_ = calc_excit_ .OR. ( excit_type%name == 'HSW'  &
      .AND. ABS(time_cur - excit_type%time0 - excit_omega/2.0_dp)<=excit_omega/2.0_dp )
   calc_excit_ = calc_excit_ .OR. ( excit_type%name == 'STEP' &
      .AND. ABS(time_cur - excit_type%time0 - excit_omega/2.0_dp)<=excit_omega/2.0_dp )
   calc_excit_ = calc_excit_ .OR. ( excit_type%name == 'DEL'  &
      .AND. ABS(time_cur - excit_type%time0)<time_step_cur/2.0_dp )
   if( itau == 0 ) calc_excit_=.FALSE.
   if ( calc_excit_ ) then
     call calculate_excit_field(time_cur,excit_field)
     excit_field_norm = NORM2(excit_field(:))
     do idir=1,3
       do ispin=1, nspin
         hamiltonian_cmplx(:,:,ispin) = hamiltonian_cmplx(:,:,ispin) - dipole_ao(:,:,idir) * excit_field(idir)
         en_tddft%excit = en_tddft%excit + REAL(SUM(dipole_ao(:,:,idir) * excit_field(idir) * p_matrix_cmplx(:,:,ispin)),dp)
       enddo
     end do
   else
     excit_field_norm = 0.0_dp
   end if

 !
 ! Projectile excitation
 case(EXCIT_PROJECTILE)

   !
   ! Move the projectile
   call change_position_one_atom(ncenter_nuclei, &
                                 xatom_start(:,ncenter_nuclei) + vel_nuclei(:,ncenter_nuclei) * ( time_cur - time_read ))

   call nucleus_nucleus_energy(en_tddft%nuc_nuc)

   !
   ! Nucleus-electron interaction due to the projectile only
   projectile_list(1) = ncenter_nuclei
   allocate(hamiltonian_projectile(basis%nbf,basis%nbf))
   call setup_nucleus(basis,hamiltonian_projectile,projectile_list)

   do ispin=1,nspin
     hamiltonian_cmplx(:,:,ispin) = hamiltonian_cmplx(:,:,ispin) + hamiltonian_projectile(:,:)
   enddo
   en_tddft%excit = REAL( SUM( hamiltonian_projectile(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
   deallocate(hamiltonian_projectile)

 end select

 do ispin=1,nspin
   hamiltonian_cmplx(:,:,ispin) = hamiltonian_cmplx(:,:,ispin) + hamiltonian_kinetic(:,:) + hamiltonian_nucleus(:,:)
 enddo

 ! Perform the canonical transform from the original basis to the orthogonal one
 call transform_hamiltonian_ortho(x_matrix,hamiltonian_cmplx,h_small_cmplx)


 ! kinetic and nuclei-electrons energy contributions
 en_tddft%kinetic = REAL( SUM( hamiltonian_kinetic(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)
 en_tddft%nucleus = REAL( SUM( hamiltonian_nucleus(:,:) * SUM(p_matrix_cmplx(:,:,:),DIM=3) ), dp)

 call stop_clock(timing_tddft_hamiltonian)

end subroutine setup_hamiltonian_cmplx


!=========================================================================
subroutine diagonalize_hamiltonian_ortho(h_small_cmplx,a_matrix_orth_cmplx,energy_tddft)
 implicit none

 complex(dp),intent(in)  :: h_small_cmplx(:,:)
 complex(dp),intent(out) :: a_matrix_orth_cmplx(:,:)
 real(dp),intent(out)    :: energy_tddft(:)
!=====
 logical                 :: algo_cmplx
 integer                 :: nstate
 real(dp),allocatable    :: m_tmpr(:,:)
!=====

 call start_clock(timing_propagate_diago)
 algo_cmplx = ANY( ABS(h_small_cmplx(:,:)%im) > 1.0e-12_dp )
 nstate = SIZE(h_small_cmplx(:,:),DIM=1)

 if( algo_cmplx ) then
   !
   ! COMPLEX case
   !
   write(stdout,'(1x,a)') 'Perform a complex matrix diagonalization'
   a_matrix_orth_cmplx(:,:) = h_small_cmplx(:,:)
   call diagonalize_scalapack(postscf_diago_flavor,scalapack_block_min,a_matrix_orth_cmplx,energy_tddft)
 else
   !
   ! REAL case
   !
   write(stdout,'(1x,a)') 'Perform a real matrix diagonalization'
   allocate(m_tmpr(nstate,nstate))

   m_tmpr(:,:) = h_small_cmplx(:,:)%re
   call diagonalize(postscf_diago_flavor,m_tmpr,energy_tddft)
   a_matrix_orth_cmplx(:,:) = m_tmpr(:,:)

   deallocate(m_tmpr)
 endif

 call stop_clock(timing_propagate_diago)

end subroutine diagonalize_hamiltonian_ortho


!=========================================================================
subroutine transform_hamiltonian_ortho(x_matrix,h_cmplx,h_small_cmplx)
 implicit none

 real(dp),intent(in)     :: x_matrix(:,:)
 complex(dp),intent(in)  :: h_cmplx(:,:,:)
 complex(dp),intent(out) :: h_small_cmplx(:,:,:)
!=====
 logical                 :: algo_cmplx
 integer                 :: nbf,nstate,ispin
 real(dp),allocatable    :: m_tmpr1(:,:),m_tmpr2(:,:)
 complex(dp),allocatable :: m_tmp(:,:),x_matrix_cmplx(:,:)
!=====

 call start_clock(timing_tddft_ham_orthobasis)

 ! If any coefficient of H has an imaginary part, use the former (and slower) algo
 algo_cmplx = ANY( ABS(h_cmplx(:,:,:)%im ) > 1.0e-12_dp )

#if defined(HAVE_MKL)
 if( algo_cmplx ) then
   write(stdout,'(1x,a)') 'Transform the hamiltonian into the canonical orthogonal basis: MKL extension, complex '
 else
   write(stdout,'(1x,a)') 'Transform the hamiltonian into the canonical orthogonal basis: MKL extension, real'
 endif
#else
 if( algo_cmplx ) then
   write(stdout,'(1x,a)') 'Transform the hamiltonian into the canonical orthogonal basis: complex '
 else
   write(stdout,'(1x,a)') 'Transform the hamiltonian into the canonical orthogonal basis: real'
 endif
#endif

 nbf    = SIZE(x_matrix,DIM=1)
 nstate = SIZE(x_matrix,DIM=2)

 !   h_small_cmplx(:,:,ispin) = MATMUL( TRANSPOSE(x_matrix(:,:)) , &
 !                   MATMUL( h_cmplx(:,:,ispin) , x_matrix(:,:) ) )

 do ispin=1,nspin
   !call matmul_transaba_scalapack(scalapack_block_min,x_matrix_cmplx,h_cmplx(:,:,ispin),h_small_cmplx(:,:,ispin))

   ! Select between two cases for speed: COMPLEX or REAL
   if( algo_cmplx ) then
     !
     ! COMPLEX CASE
     !

     allocate(x_matrix_cmplx(nbf,nstate))
     allocate(m_tmp(nbf,nstate))
     x_matrix_cmplx(:,:) = x_matrix(:,:)
     call ZHEMM('L','L',nbf,nstate,COMPLEX_ONE,h_cmplx(1,1,ispin),nbf, &
                                               x_matrix_cmplx(1,1),nbf,     &
                                  COMPLEX_ZERO,m_tmp(1,1),nbf)
#if defined(HAVE_MKL)
     call ZGEMMT('L','C','N',nstate,nbf,COMPLEX_ONE,x_matrix_cmplx(1,1),nbf, &
                                                    m_tmp(1,1),nbf,          &
                                        COMPLEX_ZERO,h_small_cmplx(1,1,ispin),nstate)
     call matrix_lower_to_full(h_small_cmplx(:,:,ispin))
#else
     call ZGEMM('C','N',nstate,nstate,nbf,COMPLEX_ONE,x_matrix_cmplx(1,1),nbf, &
                                                      m_tmp(1,1),nbf,  &
                                          COMPLEX_ZERO,h_small_cmplx(1,1,ispin),nstate)
#endif

     deallocate(m_tmp)
     deallocate(x_matrix_cmplx)

   else
     !
     ! REAL CASE
     !

     allocate(m_tmpr1(nbf,nbf))
     allocate(m_tmpr2(nbf,nstate))
     m_tmpr1(:,:) = h_cmplx(:,:,ispin)%re
     call DSYMM('L','L',nbf,nstate,1.0d0,m_tmpr1,nbf, &
                                         x_matrix,nbf,  &
                                   0.0d0,m_tmpr2,nbf)
     deallocate(m_tmpr1)
     allocate(m_tmpr1(nstate,nstate))
#if defined(HAVE_MKL)
     call DGEMMT('L','T','N',nstate,nbf,1.0d0,x_matrix,nbf, &
                                              m_tmpr2,nbf,          &
                                        0.0d0,m_tmpr1,nstate)
     call matrix_lower_to_full(m_tmpr1)
#else
     call DGEMM('T','N',nstate,nstate,nbf,1.0d0,x_matrix,nbf, &
                                                m_tmpr2,nbf,  &
                                          0.0d0,m_tmpr1,nstate)
#endif
     h_small_cmplx(:,:,ispin) = m_tmpr1(:,:)
     deallocate(m_tmpr1)
     deallocate(m_tmpr2)


   endif
 end do ! spin loop

 call stop_clock(timing_tddft_ham_orthobasis)


end subroutine transform_hamiltonian_ortho


!=========================================================================
subroutine calculate_excit_field(time_cur,excit_field)
 implicit none
 real(dp),intent(in)      :: time_cur       ! time in au
 real(dp),intent(inout)   :: excit_field(3) ! electric field in 3 dimensions
!=====
 real(dp)                 :: excit_dir_norm(3)
!=====

 excit_dir_norm(:)=excit_type%dir(:)/NORM2(excit_type%dir(:))

 select case(excit_type%name)
 case('GAU') !Gaussian electic field
   excit_field(:) = excit_type%kappa * EXP( -( time_cur-excit_type%time0 )**2 / 2.0_dp / excit_omega**2 ) * &
                  & excit_dir_norm(:)
 case('HSW') !Hann sine window
   excit_field(:) = excit_type%kappa * SIN( pi / excit_omega * ( time_cur - excit_type%time0  ) )**2 * excit_dir_norm(:)
 case('DEL') ! Delta excitation
   excit_field(:) = excit_type%kappa * excit_dir_norm(:)
 case('STEP') ! Step excitation
   excit_field(:) = excit_type%kappa * excit_dir_norm(:)
 case default
    call die('Invalid choice for the excitation type. Change excit_type value in the input file')
 end select

end subroutine calculate_excit_field

!=========================================================================
end module m_tddft_propagator
!=========================================================================
