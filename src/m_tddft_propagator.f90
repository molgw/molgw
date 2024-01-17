!=========================================================================
! This file is part of MOLGW.
! Author: Ivan Maliyov, Xixi QI
!
! This module contains
! the time propagation of the KS wavefunctions for TDDFT
!
!=========================================================================
#include "molgw.h"
module m_tddft_propagator
  use m_definitions
  use m_memory
  use m_warning
  use m_timing
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
  real(dp),private                   :: time_read !defaut=0.0_dp
  real(dp),allocatable,private       :: xatom_start(:,:)
  real(dp),allocatable,private       :: xbasis_start(:,:)
  real(dp),allocatable,private       :: count_atom_e(:,:), count_atom_e_copy(:,:)
  real(dp),private                   :: excit_field_norm
  logical,private                    :: moving_basis
  !==hamiltonian extrapolation variables==
  real(dp),allocatable,private       :: extrap_coefs(:)
  complex(dp),allocatable,private    :: h_small_hist_cmplx(:,:,:,:)
  complex(dp),allocatable,private    :: c_matrix_orth_hist_cmplx(:,:,:,:)
  complex(dp),allocatable,private    :: h_hist_cmplx(:,:,:,:)
  complex(dp),allocatable,private    :: c_matrix_hist_cmplx(:,:,:,:)
  integer,private                    :: ntau
  !==C(t0) initialization variables==
  integer,allocatable                :: atom_state_occ(:,:)
  real(dp),allocatable               :: m_eigenval(:,:)
  complex(dp),allocatable,private    :: p_matrix_cmplx_hist(:,:,:)
  complex(dp),allocatable            :: m_matrix_small(:,:,:) ! M' = X**H * ( H - i*D ) * X
  complex(dp),allocatable            :: m_eigenvec_small(:,:,:), m_eigenvector(:,:,:)
  !==frozen core==
  real(dp),allocatable               :: energies_start(:,:)
  complex(dp),allocatable            :: a_matrix_orth_start_cmplx(:,:,:)
  !====
  complex(dp),allocatable            :: gos_ao(:,:),gos_mo(:,:)
  complex(dp)                           :: norm

  type(energy_contributions),private :: en_tddft
  type(basis_set),private            :: basis_t,basis_p
  !type(basis_set),private            :: auxil_basis_t,auxil_basis_p

contains


!=========================================================================
subroutine calculate_propagation(basis,auxil_basis,occupation,c_matrix,restart_tddft_is_correct)
  use m_hdf5_tools
  implicit none

  type(basis_set),intent(inout) :: basis
  type(basis_set),intent(inout) :: auxil_basis
  real(dp),intent(in)        :: c_matrix(:,:,:)
  real(dp),intent(inout)     :: occupation(:,:)
  logical,intent(in)         :: restart_tddft_is_correct
  !=====
  integer                    :: fixed_atom_list(ncenter_nuclei-nprojectile)
  integer                    :: ispin
  integer                    :: istate,nstate_tmp
  real(dp)                   :: time_min
  real(dp),allocatable       :: dipole_ao(:,:,:)
  real(dp),allocatable       :: s_matrix(:,:)
  real(dp),allocatable       :: d_matrix(:,:)
  real(dp),allocatable       :: x_matrix(:,:)
  real(dp),allocatable       :: s_matrix_sqrt(:,:)
  real(dp),allocatable       :: hamiltonian_kinetic(:,:)
  real(dp),allocatable       :: hamiltonian_nucleus(:,:)
  !=====initial values
  integer                    :: nstate,info,min_index(1),ind
  real(dp),allocatable       :: energy_tddft(:,:)
  complex(dp),allocatable    :: c_matrix_cmplx(:,:,:), c_matrix_cmplx_scf(:,:,:)
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
  integer                    :: file_mulliken, file_lowdin
  real(dp)                   :: time_cur
  complex(dp),allocatable    :: p_matrix_cmplx(:,:,:)
  logical                    :: is_identity_ ! keep this varibale
  !==cube_diff varibales====================================
  real(dp),allocatable       :: cube_density_start(:,:,:,:)
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
  !==HDF5==
  integer(HID_T)             :: fid, c_mat_group, p_mat_group
  !==DMD==
  character(len=200)         :: snap_name
  complex(dp),allocatable    :: p_matrix_MO_cmplx(:,:,:)
  real(dp),allocatable       :: p_matrix_MO_block(:,:,:)

  call switch_on_rt_tddft_timers()
  call start_clock(timing_tddft_loop)

  write(stdout,'(/,/,1x,a)') '=================================================='
  write(stdout,'(x,a,/)')    'RT-TDDFT simulation'

  nstate = SIZE(occupation(:,:),DIM=1)

  ! Tweak the occupation if the number of electrons has changed from DFT to TDDFT
  ! tddft_charge = -999.0 is the default value that implies tddft_charge=charge
  if( ABS(tddft_charge+999.0_dp) > 1.0e-5_dp .AND. ABS(tddft_charge-charge)> 1.0e-5_dp ) then
    write(stdout,*) 'Set new occupations for TDDFT'
    call set_occupation(0.0_dp,electrons+charge-tddft_charge,magnetization,RESHAPE([0.0_dp],[1,1]),occupation)
  endif

  if( read_tddft_restart_ .AND. restart_tddft_is_correct ) then
    nocc = get_nocc_from_restart()
  else
    nocc = get_number_occupied_states(occupation)
  end if
  !
  !
  ! Initialize main switch moving_basis == .TRUE.  => propagate C
  !                                                => effective hamiltonian = S^{-1} (H - iD)
  !                        moving_basis == .FALSE. => propagate C' = X^{-1} C
  !                                                => effective hamiltonian = H
  moving_basis = excit_type%form == EXCIT_PROJECTILE_W_BASIS .OR. pred_corr(1:2)=='MB'

  write(stdout,*) 'Splitting basis set into TARGET and PROJECTILE basis sets'
  call split_basis_set(basis,basis_t,basis_p)
  call init_libcint(basis_t,basis_p)
  !if( has_auxil_basis ) then
  !  write(stdout,'(/,a)') 'Splitting up the auxiliary basis set'
  !  call split_basis_set(auxil_basis,auxil_basis_t,auxil_basis_p)
  !end if

  call echo_tddft_variables()

  call clean_allocate('Overlap matrix S for TDDFT',s_matrix,basis%nbf,basis%nbf)
  call clean_allocate('Time derivative matrix D for TDDFT',d_matrix,basis%nbf,basis%nbf)
  call clean_allocate('Kinetic operator T for TDDFT',hamiltonian_kinetic,basis%nbf,basis%nbf)
  call clean_allocate('Nucleus operator V for TDDFT',hamiltonian_nucleus,basis%nbf,basis%nbf)

  call clean_allocate('Wavefunctions C for TDDFT',c_matrix_cmplx,basis%nbf,nocc,nspin)
  call clean_allocate('Wavefunctions in ortho base C'' for TDDFT',c_matrix_orth_cmplx,nstate,nocc,nspin)
  call clean_allocate('Hamiltonian for TDDFT',h_cmplx,basis%nbf,basis%nbf,nspin)
  call clean_allocate('h_small_cmplx for TDDFT',h_small_cmplx,nstate,nstate,nspin)
  call clean_allocate('p_matrix_cmplx for TDDFT',p_matrix_cmplx,basis%nbf,basis%nbf,nspin)
  call clean_allocate('Square-Root of Overlap S{1/2}',s_matrix_sqrt,basis%nbf,basis%nbf)

  allocate(xatom_start(3,ncenter_nuclei))
  allocate(xbasis_start(3,ncenter_basis))

  write(stdout,'(/,1x,a)') "===INITIAL CONDITIONS==="

  ! Getting c_matrix_cmplx(t=0) whether using RESTART_TDDFT file, whether using real c_matrix
  if( read_tddft_restart_ .AND. restart_tddft_is_correct ) then
    if( moving_basis ) then
      call read_restart_tddft(nstate,time_read,occupation,c_matrix_cmplx)
    else
      ! assign xatom_start, c_matrix_orth_cmplx, time_min with values given in RESTART File
      call read_restart_tddft(nstate,time_read,occupation,c_matrix_orth_cmplx)
    end if

    time_min = time_read
    call change_position_one_atom(ncenter_nuclei, xatom_start(:,ncenter_nuclei))
    call change_basis_center_one_atom(ncenter_basis, xbasis_start(:,ncenter_basis))
    if( excit_type%form == EXCIT_PROJECTILE_W_BASIS ) call update_basis_eri(basis,auxil_basis)

  else
    c_matrix_cmplx(:,:,:) = c_matrix(:,1:nocc,:)
    !c_matrix_cmplx_scf(:,:,:) = c_matrix(:,1:nocc,:)
    xatom_start=xatom
    xbasis_start=xbasis
    time_min=0.0_dp
  end if

  call setup_overlap(basis,s_matrix)
  call setup_sqrt_overlap(s_matrix,s_matrix_sqrt)
  call setup_x_matrix(min_overlap,s_matrix,nstate_tmp,x_matrix)
  ! x_matrix is now allocated with dimension (basis%nbf,nstate))

  if( excit_type%form == EXCIT_PROJECTILE_W_BASIS ) then
    call setup_d_matrix(basis,d_matrix,.FALSE.)
    call setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)
    en_tddft%id = REAL( SUM( im*d_matrix(:,:) * CONJG(SUM(p_matrix_cmplx(:,:,:),DIM=3)) ), dp)
  else
    d_matrix(:,:) = 0.0_dp
    if( nstate /= nstate_tmp ) then
      call die('Error with nstate in the TDDFT propagator')
    end if
  end if

  if( read_tddft_restart_ .AND. restart_tddft_is_correct ) then
    if( .NOT. moving_basis ) then
      do ispin=1,nspin
        c_matrix_cmplx(:,:,ispin) = MATMUL( x_matrix(:,:) , c_matrix_orth_cmplx(:,:,ispin) )
      end do
    end if
  end if

  call nucleus_nucleus_energy(en_tddft%nuc_nuc)

  ! Setup the fixed part of the Hamiltonian: the kinetic energy and the fixed nuclei potential
  !
  call setup_kinetic(basis,hamiltonian_kinetic)

  ! hamiltonian_nucleus contains the contribution from all the fixed atoms (i.e. excluding the projectile)
  ! Remember: the projectile is always the last atom
  do iatom=1,ncenter_nuclei-nprojectile
    fixed_atom_list(iatom) = iatom
  enddo
  call setup_nucleus(basis,hamiltonian_nucleus,fixed_atom_list)

  if( nelement_ecp > 0 ) then
    call setup_nucleus_ecp(basis,hamiltonian_nucleus)
  endif

  if(write_step / time_step - NINT( write_step / time_step ) > 1.0E-10_dp .OR. write_step < time_step ) then
    call die("Tddft error: write_step is not a multiple of time_step or smaller than time_step.")
  end if

  if( calc_type%is_dft ) then
    call init_dft_grid(basis,tddft_grid_level,dft_xc(1)%needs_gradient,.TRUE.,BATCH_SIZE)
  endif

  ! Getting starting value of the Hamiltonian
  ! itau=0 to avoid excitation calculation
  call setup_hamiltonian_cmplx(basis,                   &
                               nstate,                  &
                               0,                       &    ! itau
                               time_min,                &
                               0.0_dp,                  &    ! time_step
                               occupation,              &
                               c_matrix_cmplx,          &
                               hamiltonian_kinetic,     &
                               hamiltonian_nucleus,     &
                               h_small_cmplx,           &
                               x_matrix,                &
                               dipole_ao,               &
                               h_cmplx,en_tddft)


  !
  ! If not restarting, initialize the coefficients C(t=0)
  if( (.NOT. read_tddft_restart_) .OR. (.NOT. restart_tddft_is_correct)) then
    if( excit_type%form == EXCIT_PROJECTILE_W_BASIS ) then
      select case(capitalize(tddft_wfn_t0))
      case('STATIONARY')
        ! initialize the wavefunctions to be the eigenstates of M = H - i*D + m*v**2*S
        ! which are also that of  U = S**-1 * ( H - i*D )
        call stationary_c_matrix(basis,               &
                                 time_min,            &
                                 s_matrix,            &
                                 x_matrix,            &
                                 d_matrix,            &
                                 occupation ,         &
                                 hamiltonian_kinetic, &
                                 hamiltonian_nucleus, &
                                 dipole_ao,           &
                                 c_matrix_cmplx,      &
                                 c_matrix_orth_cmplx, &
                                 h_cmplx,             &
                                 h_small_cmplx,       &
                                 en_tddft)
      case('SCF')
        write(stdout,'(/,1x,a)') '===== C matrix initialization is skipped ====='
      case default
        call die('calculate_propagation: tddft_wfn_t0 value not recognized')
      end select

    else
      ! In case of no restart, find the c_matrix_orth_cmplx by diagonalizing h_small
      call clean_allocate('c_matrix_buf for TDDFT',c_matrix_orth_start_complete_cmplx,nstate,nstate,nspin)
      allocate(energy_tddft(nstate,nspin))
      do ispin=1, nspin
        call diagonalize(postscf_diago_flavor,h_small_cmplx(:,:,ispin),energy_tddft(:,ispin),&
             c_matrix_orth_start_complete_cmplx(:,:,ispin))
      end do
      ! in order to save the memory, we dont keep inoccupied states (nocc+1:nstate)
      c_matrix_orth_cmplx(1:nstate,1:nocc,1:nspin)=c_matrix_orth_start_complete_cmplx(1:nstate,1:nocc,1:nspin)
      call clean_deallocate('c_matrix_buf for TDDFT',c_matrix_orth_start_complete_cmplx)
      deallocate(energy_tddft)
    end if
  end if

  !!!!! Test for |<psi_scf|e^ivr|psi_t0>|^2 = 1
  !call calculate_gos_ao_mb(basis,gos_ao)
  !allocate(gos_mo(nocc,nocc))
  !gos_mo(:,:) = MATMUL( TRANSPOSE(CONJG( c_matrix_cmplx(:,:,1) )), MATMUL( gos_ao(:,:), c_matrix_cmplx_scf(:,:,1) ) )
  !deallocate(gos_ao)
  !call clean_deallocate('Wavefunctions C for TDDFT',c_matrix_cmplx_scf)
  !norm = matrix_trace_cmplx(MATMUL( CONJG(TRANSPOSE(gos_mo)), gos_mo ))
  !deallocate(gos_mo)
  !write(stdout,*) 'NORM = ', norm

  !!!!! End of test

  ! E_iD = - Tr{P*iD}
  !en_tddft%id = REAL( SUM( im*d_matrix(:,:) * CONJG(SUM(p_matrix_cmplx(:,:,:),DIM=3)) ), dp)

  ! Number of time steps
  ntau = NINT( (time_sim-time_min) / time_step )

  if(excit_type%form==EXCIT_LIGHT) then
    call clean_allocate('Dipole_basis for TDDFT',dipole_ao,basis%nbf,basis%nbf,3)
    call setup_dipole_ao(basis,dipole_ao)
  end if

  if( print_dens_traj_tddft_ ) then
    call plot_rho_traj_bunch_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,0,0.d0)
  end if


  if( calc_q_matrix_ ) then
    call initialize_q(nstate,nocc,nspin,c_matrix_orth_start_complete_cmplx,h_small_cmplx,istate_cut,file_q_matrix)
  end if

  !==frozen core: energy initialization==
  if( ncore_tddft > 0 ) then
    call clean_allocate('Initial energies for the frozen core',energies_start,nstate,nspin)
    call clean_allocate('a_matrix_orth_start_cmplx for the frozen core',a_matrix_orth_start_cmplx,nstate,nstate,nspin)
    do ispin=1, nspin
      call diagonalize(postscf_diago_flavor,h_small_cmplx(:,:,ispin),energies_start(:,ispin),a_matrix_orth_start_cmplx(:,:,ispin))
    end do
  end if
  !====

  !===cube_diff matrix allocation and parameters initialization
  if(print_cube_diff_tddft_) then
    allocate(cube_density_start(cube_nx,cube_ny,cube_nz,nspin))
  end if

  if(print_line_rho_diff_tddft_) then
    call initialize_rho_diff_cmplx(nr_line_rho,point_a,point_b)
    allocate(rho_start(nr_line_rho,nspin))
  end if

  time_min = time_read

  !
  ! Opening files and writing headers in files
  call initialize_files(file_time_data,file_dipole_time,file_excit_field,file_mulliken,file_lowdin)

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
    call plot_cube_diff_cmplx(basis,occupation,c_matrix_cmplx,initialize=.TRUE.)
  end if

  if( print_p_matrix_MO_block_hdf5_ ) then
    allocate(p_matrix_MO_cmplx(nstate,nstate,nspin))
    allocate(p_matrix_MO_block(nocc,nstate-nocc,nspin))

    call setup_density_matrix_MO_cmplx(c_matrix, s_matrix, p_matrix_cmplx, p_matrix_MO_cmplx)
    p_matrix_MO_block(:,:,:) = REAL(p_matrix_MO_cmplx(1:nocc,nocc+1:nstate,:), dp)
  end if

  ! HANDLING HDF5 files here
  if( (print_c_matrix_cmplx_hdf5_ .or. print_p_matrix_MO_block_hdf5_) .and. is_iomaster ) then

    call hdf_open_file(fid, 'rt_tddft.h5', status='NEW')

    call hdf_write_dataset(fid, 'nbf', basis%nbf)
    call hdf_write_dataset(fid, 'nstate', nstate)
    call hdf_write_dataset(fid, 'nocc', nocc)

    call hdf_write_dataset(fid, 'time_step', time_step)
    call hdf_write_dataset(fid, 'occupation', occupation)
    call hdf_write_dataset(fid, 's_matrix', s_matrix)

    ! save the initial complete c_matrix, nstate x nstate
    call hdf_write_dataset(fid, 'c_matrix_complete_0_real', c_matrix)

    if( excit_type%form == EXCIT_LIGHT ) call hdf_write_dataset(fid, 'dipole_ao', dipole_ao)

    if( print_c_matrix_cmplx_hdf5_ ) then
      call hdf_create_group(fid, 'c_matrix')
      call hdf_open_group(fid, 'c_matrix', c_mat_group)
      call dump_matrix_cmplx_hdf5(c_mat_group, c_matrix_cmplx, 0)
    end if

    if( print_p_matrix_MO_block_hdf5_ ) then

      call hdf_create_group(fid, 'p_matrix_MO_block')
      call hdf_open_group(fid, 'p_matrix_MO_block', p_mat_group)
      call hdf_write_dataset(p_mat_group, 'snap_0', p_matrix_MO_block)

    end if

  end if

  if(print_line_rho_diff_tddft_) then
    call calc_rho_initial_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,0,time_min,nr_line_rho,point_a,point_b,rho_start)
    call plot_rho_diff_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,0,time_min,nr_line_rho,point_a,point_b,rho_start)
  end if

  if( calc_dens_disc_ )       call calc_density_in_disc_cmplx_dft_grid(basis,occupation,c_matrix_cmplx,0,time_min)
  ! if( calc_dens_disc_ )       call calc_density_in_disc_cmplx_regular(nstate,nocc,basis,occupation,c_matrix_cmplx,0,time_min)

  if( print_line_rho_tddft_ ) call plot_rho_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,0,time_min)

  if( print_charge_tddft_ ) then
    !call mulliken_pdos_cmplx(basis,s_matrix,c_matrix_cmplx,occupation,file_mulliken,time_min)
    call lowdin_pdos_cmplx(basis,s_matrix_sqrt,c_matrix_cmplx,occupation,file_lowdin,time_min)
  end if

  call print_tddft_values(time_min,file_time_data,file_dipole_time,file_excit_field,0)
  en_tddft%time = time_min
  write(time_key,'(i8)') 0
  call print_energy_yaml('tddft energy '//TRIM(ADJUSTL(time_key)),en_tddft)

  time_min = time_min + time_step


  !
  ! Extrapolation coefficients and history c_ and h_ matrices (h_small_hist_cmplx)
  call initialize_extrap_coefs(c_matrix_orth_cmplx,h_small_cmplx,c_matrix_cmplx,h_cmplx)

  write(stdout,'(/,1x,a,/)') "===END OF INITIAL CONDITIONS==="


  !
  ! TDDFT time loop
  !
  time_cur = time_min
  iwrite_step = 1
  itau = 1

  do while ( (time_cur - time_sim) < 1.0e-10 )
    if ( itau == 3 ) then
      call start_clock(timing_tddft_one_iter)
    end if

    !
    ! Use c_matrix_orth_cmplx and h_small_cmplx at (time_cur-time_step) as start values,
    ! then use chosen predictor-corrector sheme to calculate c_matrix_cmplx, c_matrix_orth_cmplx,
    ! h_cmplx and h_small_cmplx and time_cur.
    call predictor_corrector(basis,                  &
                             auxil_basis,            &
                             c_matrix_cmplx,         &
                             c_matrix_orth_cmplx,    &
                             h_cmplx,                &
                             h_small_cmplx,          &
                             x_matrix,               &
                             s_matrix,               &
                             d_matrix,               &
                             itau,                   &
                             time_cur,               &
                             occupation,             &
                             hamiltonian_kinetic,    &
                             hamiltonian_nucleus,    &
                             dipole_ao)


    !
    ! debug
    !do ispin=1,nspin
    !  is_identity_ = check_identity_cmplx(MATMUL(MATMUL(TRANSPOSE(CONJG( &
    !  c_matrix_cmplx(:,:,ispin))),s_matrix(:,:)), c_matrix_cmplx(:,:,ispin) ))
    !  if( .NOT. is_identity_) then
    !    write(stdout,'(1x,a,i4,a,i7)') 'C**H*S*C is not identity for spin ', &
    !    ispin,' at itau= ', itau
    !  end if
    !enddo

    call setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)
    en_tddft%id = REAL( SUM( im*d_matrix(:,:) * CONJG(SUM(p_matrix_cmplx(:,:,:),DIM=3)) ), dp)

    write(stdout,'(1x,a,2(1x,es15.8))') 'Number of electrons from Tr(PS): ', &
                                        SUM( SUM( p_matrix_cmplx(:,:,:), DIM=3 ) * s_matrix(:,:) )

    !
    ! Print tddft values into diferent files: 1) standart output; 2) time_data.dat; 3) dipole_time.dat; 4) excitation_time.dat;
    ! 5) Mulliken/Lowdin charge file.
    ! 3) and 4) in case of light excitation. 5) in case of charge analysis.

    if ( ABS(time_cur / (write_step)- NINT(time_cur / (write_step))) < 1.0e-7 ) then

      if ( excit_type%form == EXCIT_LIGHT ) then
        call setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)
        call static_dipole(basis,p_matrix_in=p_matrix_cmplx,dipole_ao_in=dipole_ao,dipole_out=dipole)
      end if

      en_tddft%total = en_tddft%nucleus + en_tddft%kinetic + en_tddft%nuc_nuc &
                      + en_tddft%hartree + en_tddft%exx_hyb + en_tddft%xc + en_tddft%excit

      call print_tddft_values(time_cur,file_time_data,file_dipole_time,file_excit_field,itau)
      en_tddft%time = time_cur
      write(time_key,'(i8)') iwrite_step
      call print_energy_yaml('tddft energy '//TRIM(ADJUSTL(time_key)),en_tddft)
      !FBFB moving_basis spurious forces correction
      !call setup_mb_force(basis,s_matrix,c_matrix_cmplx,h_cmplx,occupation,.FALSE.)

      if ( print_line_rho_tddft_  )     call plot_rho_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,iwrite_step,time_cur)
      if ( print_line_rho_diff_tddft_ ) call plot_rho_diff_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx, &
                                                                iwrite_step,time_cur,nr_line_rho,point_a,point_b,rho_start)
      if( print_cube_rho_tddft_  )      call plot_cube_wfn_cmplx(nstate,nocc,basis,occupation,c_matrix_cmplx,iwrite_step)
      if( print_cube_diff_tddft_ )      call plot_cube_diff_cmplx(basis,occupation,c_matrix_cmplx)
      if( calc_dens_disc_ )             call calc_density_in_disc_cmplx_dft_grid(basis,occupation,c_matrix_cmplx, &
                                                                                iwrite_step,time_cur)
      !if ( calc_dens_disc_ )       call calc_density_in_disc_cmplx_regular(nstate,nocc,basis,occupation,c_matrix_cmplx,iwrite_step,time_cur)

      if (calc_q_matrix_) call calculate_q_matrix(occupation,c_matrix_orth_start_complete_cmplx,c_matrix_orth_cmplx, &
                                                 istate_cut,file_q_matrix,time_cur)

      if( print_c_matrix_cmplx_hdf5_ .and. is_iomaster ) call dump_matrix_cmplx_hdf5(c_mat_group, c_matrix_cmplx, iwrite_step)
      if( print_p_matrix_MO_block_hdf5_) then
        call setup_density_matrix_MO_cmplx(c_matrix, s_matrix, p_matrix_cmplx, p_matrix_MO_cmplx)
        p_matrix_MO_block(:,:,:) = REAL(p_matrix_MO_cmplx(1:nocc,nocc+1:nstate,:), dp)
        write(snap_name, '(A,I0)') 'snap_', iwrite_step
        if (is_iomaster) call hdf_write_dataset(p_mat_group, TRIM(snap_name), p_matrix_MO_block)
      end if

      iwrite_step = iwrite_step + 1

    end if

    if ( ABS(time_cur / (calc_charge_step)- NINT(time_cur / (calc_charge_step))) < 1.0e-4 ) then
      if ( print_charge_tddft_ ) then
        if ( excit_type%form == EXCIT_PROJECTILE_W_BASIS ) then
          call setup_sqrt_overlap(s_matrix,s_matrix_sqrt)
        end if
        !call mulliken_pdos_cmplx(basis,s_matrix,c_matrix_cmplx,occupation,file_mulliken,time_cur)
        call lowdin_pdos_cmplx(basis,s_matrix_sqrt,c_matrix_cmplx,occupation,file_lowdin,time_cur)
      end if
    end if

    !
    !--TIMING of one iteration--
    if ( itau == 3 ) then
      call stop_clock(timing_tddft_one_iter)
      call output_timing_one_iter()
    end if

    ! ---print tdddft restart each n_restart_tddft steps---
    if ( print_tddft_restart_ .AND. mod(itau,n_restart_tddft)==0 ) then
      if( moving_basis ) then
        call write_restart_tddft(nstate,time_cur,occupation,c_matrix_cmplx)
      else
        call write_restart_tddft(nstate,time_cur,occupation,c_matrix_orth_cmplx)
      end if
    end if

    time_cur = time_min + itau*time_step
    itau = itau + 1

    !---
  end do

  !********end time loop*******************

  if( (print_c_matrix_cmplx_hdf5_ .or. print_p_matrix_MO_block_hdf5_) .and. is_iomaster ) then

    call hdf_write_dataset(fid, 'nsnap', itau)

    if(print_c_matrix_cmplx_hdf5_) call hdf_close_group(c_mat_group)
    if(print_p_matrix_MO_block_hdf5_) call hdf_close_group(p_mat_group)
    call hdf_close_file(fid)
  end if

  if(print_tddft_restart_) then
    if( moving_basis ) then
      !time_cur-time_step to be consistent with the actual last moment of the simulation
      call write_restart_tddft(nstate,time_cur-time_step,occupation,c_matrix_cmplx)
    else
      call write_restart_tddft(nstate,time_cur-time_step,occupation,c_matrix_orth_cmplx)
    end if
  end if

  if( is_iomaster) then
    close(file_time_data)
    if( excit_type%form==EXCIT_LIGHT) then
      close(file_dipole_time)
      close(file_excit_field)
    end if
    if( print_charge_tddft_ ) then
      !close(file_mulliken)
      close(file_lowdin)
    end if
  end if

  call clean_deallocate('p_matrix_cmplx for TDDFT',p_matrix_cmplx)
  call clean_deallocate('h_small_hist_cmplx for TDDFT',h_small_hist_cmplx)
  call clean_deallocate('c_matrix_orth_hist_cmplx for TDDFT',c_matrix_orth_hist_cmplx)
  call clean_deallocate('h_hist_cmplx for TDDFT',h_hist_cmplx)
  call clean_deallocate('c_matrix_hist_cmplx for TDDFT',c_matrix_hist_cmplx)

  if(ncore_tddft > 0) then
    call clean_deallocate('a_matrix_orth_start_cmplx for the frozen core',a_matrix_orth_start_cmplx)
    call clean_deallocate('Inial energies for frozen core',energies_start)
  end if

  if(ALLOCATED(extrap_coefs)) deallocate(extrap_coefs)
  if(ALLOCATED(cube_density_start)) deallocate(cube_density_start)

  if( calc_type%is_dft ) call destroy_dft_grid()

  call clean_deallocate('Overlap matrix S for TDDFT',s_matrix)
  call clean_deallocate('Time derivative matrix D for TDDFT',d_matrix)
  call clean_deallocate('Transformation matrix X',x_matrix)
  call clean_deallocate('Square-Root of Overlap S{1/2}',s_matrix_sqrt)
  call clean_deallocate('Kinetic operator T for TDDFT',hamiltonian_kinetic)
  call clean_deallocate('Nucleus operator V for TDDFT',hamiltonian_nucleus)

  deallocate(xatom_start)
  deallocate(xbasis_start)

  call clean_deallocate('Dipole_basis for TDDFT',dipole_ao)

  call clean_deallocate('Wavefunctions C for TDDFT',c_matrix_cmplx)
  call clean_deallocate('Wavefunctions in ortho base C'' for TDDFT',c_matrix_orth_cmplx)
  call clean_deallocate('Hamiltonian for TDDFT',h_cmplx)
  call clean_deallocate('h_small_cmplx for TDDFT',h_small_cmplx)

  call destroy_basis_set(basis_t)
  call destroy_basis_set(basis_p)
  !if( has_auxil_basis ) then
  !  call destroy_basis_set(auxil_basis_t)
  !  call destroy_basis_set(auxil_basis_p)
  !end if

  write(stdout,'(/,x,a)') "End of RT-TDDFT simulation"
  write(stdout,'(1x,a,/)') '=================================================='

  call stop_clock(timing_tddft_loop)
  call switch_off_rt_tddft_timers()

end subroutine calculate_propagation


!=========================================================================
subroutine echo_tddft_variables()
  implicit none

  write(stdout,'(/,1x,a)') 'The most important variables of this section:'
  write(stdout,'(2x,a32,2x,es16.8)') 'Simulation time:',time_sim
  write(stdout,'(2x,a32,2x,es16.8)') 'Time step:',time_step
  write(stdout,'(2x,a32,2x,i8)') 'Number of time steps:',NINT((time_sim)/time_step)
  write(stdout,'(2x,a32,6x,l1)') 'Moving basis:',moving_basis
  write(stdout,'(2x,a32,6x,a)')  'Initial wavefunctions:',TRIM(tddft_wfn_t0)
  write(stdout,'(2x,a32,2x,f14.6)') 'Charge:',tddft_charge
  write(stdout,'(2x,a32,6x,a)')      'Predictor-corrector:',TRIM(pred_corr)
  write(stdout,'(2x,a32,6x,a)')      'Propagator:',TRIM(prop_type)
  write(stdout,'(2x,a32,2x,i8)')     'Number of propagated states:',nocc
  write(stdout,'(2x,a32,2x,i8,/)')     'Hamiltonian history length for PC:',n_hist

end subroutine echo_tddft_variables


!=========================================================================
subroutine output_timing_one_iter()
  implicit none
  real(dp)           :: time_one_iter
  !=====
  !=====

  time_one_iter = get_timing(timing_tddft_one_iter)
  write(stdout,'(/,1x,a)') '**********************************'
  write(stdout,"(1x,a32,2x,f14.6)") "Time of one iteration (s): ", time_one_iter
  write(stdout,"(1x,a32,2x,2(f12.2,1x))") "Estimated calculation time (s), (hrs):", time_one_iter*ntau,  &
                                                                                     time_one_iter*ntau/3600.0_dp
  write(stdout,'(1x,a)') '**********************************'
  flush(stdout)

end subroutine output_timing_one_iter


!=========================================================================
subroutine stationary_c_matrix(basis,               &
                               time_min,            &
                               s_matrix,            &
                               x_matrix,            &
                               d_matrix,            &
                               occupation ,         &
                               hamiltonian_kinetic, &
                               hamiltonian_nucleus, &
                               dipole_ao,           &
                               c_matrix_cmplx,      &
                               c_matrix_orth_cmplx, &
                               h_cmplx,             &
                               h_small_cmplx,       &
                               en_tddft)
  implicit none

  type(basis_set),intent(inout)   :: basis
  real(dp),intent(in)             :: time_min
  real(dp),intent(in)             :: x_matrix(:,:)
  real(dp),intent(in)             :: s_matrix(:,:)
  real(dp),intent(in)             :: d_matrix(:,:)
  real(dp),intent(in)             :: occupation(:,:)
  real(dp),intent(inout)          :: hamiltonian_kinetic(:,:)
  real(dp),intent(inout)          :: hamiltonian_nucleus(:,:)
  real(dp),intent(in),allocatable :: dipole_ao(:,:,:)
  complex(dp),intent(inout),allocatable  :: c_matrix_cmplx(:,:,:)
  complex(dp),intent(inout),allocatable  :: c_matrix_orth_cmplx(:,:,:)
  complex(dp),intent(inout)       :: h_cmplx(:,:,:)
  complex(dp),intent(inout)       :: h_small_cmplx(:,:,:)
  type(energy_contributions),intent(inout) :: en_tddft
  !====
  integer,parameter             :: nhist=4
  integer                       :: ihist,icycle,ncycle_max = 50
  integer                       :: ispin,istate,nstate
  complex(dp),allocatable       :: p_matrix_cmplx(:,:,:)
  complex(dp),allocatable       :: m_matrix_small(:,:,:) ! M' = X**H * ( H - i*D ) * X
  complex(dp),allocatable       :: m_eigenvec_small(:,:,:), m_eigenvector(:,:,:)
  complex(dp),allocatable       :: m_tmp(:,:,:)
  real(dp)                      :: rms
  !====

  nstate = SIZE(occupation(:,:),DIM=1)
  allocate( p_matrix_cmplx(basis%nbf,basis%nbf,nspin) )

  call setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)
  allocate( p_matrix_cmplx_hist, SOURCE = p_matrix_cmplx(:,:,:) )
  allocate( h_hist_cmplx(basis%nbf,basis%nbf,nspin,nhist) )
  do ihist=1,nhist
    h_hist_cmplx(:,:,:,ihist) = h_cmplx(:,:,:)
  enddo

  ! Self-consistency loop for C(t0)
  ! M = H - iD + mv**2*S is Hermitian at t0 if the projectile and the target do not overlap
  do icycle = 1, ncycle_max

    write(stdout,'(/,1x,a)')
    write(stdout,*) '=============== Stationary state iteration', icycle, '==============='
    write(stdout,'(/,1x,a)')

    allocate( m_matrix_small, MOLD = h_small_cmplx )
    allocate( m_eigenvec_small, MOLD = h_small_cmplx )
    allocate( m_eigenvector(basis%nbf,nstate,nspin), m_eigenval(nstate,nspin) )
    allocate( m_tmp(basis%nbf,basis%nbf,nspin))

    do ispin=1, nspin
      ! in ortho basis : M' = X**H * ( H-iD+(mv**2+ deltaE)*S ) * X
      m_tmp(:,:,ispin)  = h_cmplx(:,:,ispin) - im*d_matrix(:,:)
      !FBFB Why (1/2) m v**2 ?
      m_tmp(basis_t%nbf + 1:,basis_t%nbf + 1:,ispin)  = m_tmp(basis_t%nbf + 1:,basis_t%nbf + 1:,ispin) &
            + (0.5*SUM(vel_nuclei(:,ncenter_nuclei)**2) + tddft_energy_shift ) * s_matrix(basis_t%nbf + 1:,basis_t%nbf + 1:)

      m_eigenvector(:,:,ispin)  = MATMUL( m_tmp(:,:,ispin), x_matrix(:,:) )
      m_matrix_small(:,:,ispin) = MATMUL( TRANSPOSE(x_matrix(:,:)), m_eigenvector(:,:,ispin) )
      ! diagonalize M'(t0) to get eigenstates C'(t0) for MB propagation

      ! Diagonalization assumes Hermitianity of m_matrix_small
      call diagonalize( postscf_diago_flavor, m_matrix_small(:,:,ispin), &
                        m_eigenval(:,ispin), m_eigenvec_small(:,:,ispin) )
      ! M = X * M'
      m_eigenvector(:,:,ispin) = MATMUL( x_matrix(:,:) , m_eigenvec_small(:,:,ispin) )
    end do

    c_matrix_cmplx(:,1:nocc,:) = m_eigenvector(:,1:nocc,:)
    c_matrix_orth_cmplx(:,1:nocc,:) = m_eigenvec_small(:,1:nocc,:)
    if( nspin > 1 ) then
      write(stdout, '(a30,3(2x,a10))') 'Propagator eigenvalues (eV)', ' ', 'occupation'
      write(stdout, '(a10,4(2x,a10))') 'spin1', 'spin2', 'spin1', 'spin2'
    else
      write(stdout, '(a30,2(2x,a10))') 'Propagator eigenvalues (eV)', 'occupation'
    end if
    do istate=1,MIN(nocc+5,nstate)
      write(stdout, '(4x,4(2x,f12.6))') m_eigenval(istate,:)*Ha_eV, occupation(istate, :)
    end do
    deallocate(m_matrix_small)
    deallocate(m_eigenvec_small, m_eigenvector,m_tmp)
    deallocate(m_eigenval)

    call setup_hamiltonian_cmplx(basis,                    &
                                 nstate,                   &
                                 0,                        &
                                 time_min,                 &
                                 0.0_dp,                   &
                                 occupation,               &
                                 c_matrix_cmplx,           &
                                 hamiltonian_kinetic,      &
                                 hamiltonian_nucleus,      &
                                 h_small_cmplx,            &
                                 x_matrix,                 &
                                 dipole_ao,                &
                                 h_cmplx,en_tddft)

    call setup_density_matrix_cmplx(c_matrix_cmplx,occupation,p_matrix_cmplx)
    en_tddft%id = REAL( SUM( im*d_matrix(:,:) * CONJG(SUM(p_matrix_cmplx(:,:,:),DIM=3)) ), dp)

    if( icycle > 1 ) then
      rms = SQRT( SUM(( p_matrix_cmplx(:,:,:) - p_matrix_cmplx_hist(:,:,:) )**2) ) * SQRT( REAL(nspin,dp) )
      write(stdout,'(1x,a,es14.6)') 'Changes in density matrix: ',rms
      if( rms < tolscf_tddft ) then
        write(stdout,'(1x,a,/)') "=== CONVERGENCE REACHED ==="
        exit
      else
        if( icycle == ncycle_max ) call die("=== TDDFT CONVERGENCE NOT REACHED ===")
      end if
    endif

    !
    ! History mixing to damp the charge oscillations (poor man solution)
    !
    p_matrix_cmplx_hist(:,:,:) = p_matrix_cmplx(:,:,:)
    do ihist=nhist,2,-1
      h_hist_cmplx(:,:,:,ihist) = h_hist_cmplx(:,:,:,ihist-1)
    enddo
    h_hist_cmplx(:,:,:,1) = h_cmplx(:,:,:)
    ! Simple mixing of H over the last nhist iterations
    h_cmplx = SUM( h_hist_cmplx(:,:,:,:) , DIM=4) / REAL(nhist,dp)

  end do
  deallocate(p_matrix_cmplx, p_matrix_cmplx_hist, h_hist_cmplx)

end subroutine stationary_c_matrix


!=========================================================================
subroutine update_basis_eri(basis,auxil_basis)

  implicit none
  type(basis_set),intent(inout)      :: basis
  type(basis_set),intent(inout)      :: auxil_basis
  !=====

  write(stdout,'(/,a)') ' Update moving basis set'
  call moving_basis_set(basis)
  call destroy_libcint(basis)

  if( has_auxil_basis ) then
    write(stdout,'(/,a)') ' Setting up the auxiliary basis set for Coulomb integrals'
    call moving_basis_set(auxil_basis)
    call destroy_libcint(auxil_basis)
    call init_libcint(basis, auxil_basis)
    call init_libcint(auxil_basis)
    !
    ! Setup new eri 2center / 3center
    call calculate_eri_ri(basis,auxil_basis,0.0_dp)
  else
    call init_libcint(basis)
    call deallocate_eri_4center()
    call calculate_eri(print_eri_,basis,0.0_dp)
  endif

end subroutine update_basis_eri


!=========================================================================
subroutine setup_d_matrix(basis,d_matrix,recalc)
  implicit none
  type(basis_set),intent(in)          :: basis
  real(dp),intent(inout)              :: d_matrix(:,:)
  logical,intent(in)                  :: recalc
  !=====
  integer                             :: jbf
  real(dp),allocatable                :: s_matrix_grad(:,:,:)
  !=====

  write(stdout,'(/,a)') ' Setup overlap time derivative matrix D (analytic)'

  allocate(s_matrix_grad(basis%nbf,basis%nbf,3))

  ! This routine returns s_grad = < grad_R phi_alpha | phi_beta >
  if ( recalc ) then
    call recalc_overlap_grad(basis_t,basis_p,s_matrix_grad)
  else
    call setup_overlap_grad(basis,s_matrix_grad)
  end if

  ! We want D:
  !    D = < phi_alpha | d/dt phi_beta >
  !      = < phi_alpha | grad_R phi_beta > . v
  !      = s_grad**T . v

  s_matrix_grad(:,:,1) = TRANSPOSE(s_matrix_grad(:,:,1))
  s_matrix_grad(:,:,2) = TRANSPOSE(s_matrix_grad(:,:,2))
  s_matrix_grad(:,:,3) = TRANSPOSE(s_matrix_grad(:,:,3))

  if ( recalc ) then
    do jbf=basis_t%nbf+1,basis%nbf
      d_matrix(1:basis_t%nbf,jbf) = s_matrix_grad(1:basis_t%nbf,jbf,1) * basis_p%bff(1)%v0(1) &
                                   +s_matrix_grad(1:basis_t%nbf,jbf,2) * basis_p%bff(1)%v0(2) &
                                   +s_matrix_grad(1:basis_t%nbf,jbf,3) * basis_p%bff(1)%v0(3)
    end do
  else
    do jbf=1,basis%nbf
      d_matrix(:,jbf) = s_matrix_grad(:,jbf,1) * basis%bff(jbf)%v0(1) &
                      + s_matrix_grad(:,jbf,2) * basis%bff(jbf)%v0(2) &
                      + s_matrix_grad(:,jbf,3) * basis%bff(jbf)%v0(3)
    end do
  end if

  deallocate(s_matrix_grad)

end subroutine setup_d_matrix


!=========================================================================
subroutine setup_mb_force(basis,s_matrix,c_matrix_cmplx,h_cmplx,occupation,recalc)
  implicit none
  type(basis_set),intent(in)          :: basis
  real(dp),intent(in)                 :: s_matrix(:,:)
  complex(dp),intent(in)              :: c_matrix_cmplx(:,:,:)
  complex(dp),intent(in)              :: h_cmplx(:,:,:)
  real(dp),intent(in)                 :: occupation(:,:)
  !real(dp),intent(inout)              :: force_mb_p(3)
  logical,intent(in)                  :: recalc
  !=====
  integer                             :: ibf,idir,istate,jbf,ispin
  real(dp)                            :: force_mb_p(3)
  real(dp),allocatable                :: s_matrix_hess(:,:,:,:)
  real(dp),allocatable                :: s_matrix_inv(:,:)
  real(dp),allocatable                :: BboldAdagger(:,:,:)
  real(dp),allocatable                :: BAdagger(:,:)
  complex(dp),allocatable             :: DboldA1(:,:,:)
  complex(dp),allocatable             :: DboldA2(:,:,:)
  complex(dp),allocatable             :: DboldA3(:,:,:)
  real(dp),allocatable                :: cc_matrix(:,:,:)   ! cc_matrix corresponds to C^A_\alpha\beta in Kunert-Schmidt EPJ-D (2003)
  complex(dp) :: ctmp(3)
  !=====

  write(stdout,'(/,a)') ' Setup mb force (analytic)'

  allocate(s_matrix_hess(basis%nbf,basis%nbf,3,3))
  allocate(s_matrix_inv,MOLD=s_matrix)
  s_matrix_inv(:,:) = s_matrix(:,:)
  call invert(s_matrix_inv)

  allocate(BboldAdagger(basis%nbf,basis%nbf,3))
  allocate(DboldA1(basis%nbf,basis%nbf,3))
  allocate(DboldA2(basis%nbf,basis%nbf,3))
  allocate(DboldA3(basis%nbf,basis%nbf,3))
  allocate(BAdagger(basis%nbf,basis%nbf))

  call setup_overlap_grad(basis,BboldAdagger)
  do ibf=1,basis%nbf
    ! if A is not a projectile then BboldAdagger is zero
    if( ALL( ABS(basis%bff(ibf)%v0(:)) < 1.0e-6_dp ) ) then
      BboldAdagger(ibf,:,:) = 0.0_dp
    endif
    BAdagger(ibf,:) = basis%bff(ibf)%v0(1) * BboldAdagger(ibf,:,1) &
                    + basis%bff(ibf)%v0(2) * BboldAdagger(ibf,:,2) &
                    + basis%bff(ibf)%v0(3) * BboldAdagger(ibf,:,3)
  enddo

  do idir=1,3
    DboldA1(:,:,idir) = MATMUL( BboldAdagger(:,:,idir), MATMUL( s_matrix_inv, h_cmplx(:,:,1) ) ) &
                        + MATMUL( MATMUL( h_cmplx(:,:,1), s_matrix_inv), TRANSPOSE( BboldAdagger(:,:,idir) ) )
  enddo

  ctmp(:) = 0.0_dp
  ispin=1
  do idir=1,3
    do istate=1,nocc
      do ibf=1,basis%nbf
        do jbf=1,basis%nbf
          ctmp(idir) = ctmp(idir) + DboldA1(ibf,jbf,idir) &
                                       * CONJG(c_matrix_cmplx(ibf,istate,ispin)) * c_matrix_cmplx(jbf,istate,ispin) &
                                       * occupation(istate,ispin)
        enddo
      enddo
    enddo
  enddo
  write(1001,'(3(2x,es16.8))') ctmp(:)%re


  do idir=1,3
    DboldA2(:,:,idir) = im * MATMUL( BAdagger(:,:), MATMUL( s_matrix_inv, TRANSPOSE(BboldAdagger(:,:,idir)) ) ) &
                        - im * MATMUL( BboldAdagger(:,:,idir) , MATMUL( s_matrix_inv, TRANSPOSE(BAdagger(:,:)) )  )
  enddo
  ctmp(:) = 0.0_dp
  ispin=1
  do idir=1,3
    do istate=1,nocc
      do ibf=1,basis%nbf
        do jbf=1,basis%nbf
          ctmp(idir) = ctmp(idir) + DboldA2(ibf,jbf,idir) &
                                       * CONJG(c_matrix_cmplx(ibf,istate,ispin)) * c_matrix_cmplx(jbf,istate,ispin) &
                                       * occupation(istate,ispin)
        enddo
      enddo
    enddo
  enddo
  write(1002,'(3(2x,es16.8))') ctmp(:)%re



  ! This routine returns s_hess = < grad_R_A phi_alpha | grad_R_B phi_beta >
  call setup_overlap_hessian(basis,s_matrix_hess)

  ! Lets calculate < d/dt phi_alpha | grad_R_B phi_beta >
  !                  = < v_R_A . grad_R_A phi_alpha | grad_R_B phi_beta >
  allocate(cc_matrix(basis%nbf,basis%nbf,3))

  do ibf=1,basis%nbf
    cc_matrix(ibf,:,1) = basis%bff(ibf)%v0(1) * s_matrix_hess(ibf,:,1,1) &
                       + basis%bff(ibf)%v0(2) * s_matrix_hess(ibf,:,2,1) &
                       + basis%bff(ibf)%v0(3) * s_matrix_hess(ibf,:,3,1)
    cc_matrix(ibf,:,2) = basis%bff(ibf)%v0(1) * s_matrix_hess(ibf,:,1,2) &
                       + basis%bff(ibf)%v0(2) * s_matrix_hess(ibf,:,2,2) &
                       + basis%bff(ibf)%v0(3) * s_matrix_hess(ibf,:,3,2)
    cc_matrix(ibf,:,3) = basis%bff(ibf)%v0(1) * s_matrix_hess(ibf,:,1,3) &
                       + basis%bff(ibf)%v0(2) * s_matrix_hess(ibf,:,2,3) &
                       + basis%bff(ibf)%v0(3) * s_matrix_hess(ibf,:,3,3)
  enddo

  deallocate(s_matrix_hess)

  do jbf=1,basis%nbf
    ! if beta is not on a projectile, then the gradient is not contributing the for on the projectile
    if( ALL( ABS(basis%bff(jbf)%v0(:)) < 1.0e-6_dp ) ) then
      cc_matrix(:,jbf,:) = 0.0_dp
    endif
  enddo

  ctmp(:) = 0.0_dp
  ispin=1
  do idir=1,3
    do istate=1,nocc
      do ibf=1,basis%nbf
        do jbf=1,basis%nbf
          ctmp(idir) = ctmp(idir) + im * ( cc_matrix(jbf,ibf,idir) - cc_matrix(ibf,jbf,idir) ) &
                                        * CONJG(c_matrix_cmplx(ibf,istate,ispin)) * c_matrix_cmplx(jbf,istate,ispin) &
                                        * occupation(istate,ispin)
        enddo
      enddo
    enddo
  enddo


  deallocate(cc_matrix)

  write(1003,'(3(2x,es16.8))') ctmp(:)%re
  force_mb_p(:) = ctmp(:)%re

end subroutine setup_mb_force


!=========================================================================
subroutine mb_related_updates(basis,                &
                              auxil_basis,need_eri, &
                              time_cur,dt_factor,   &
                              s_matrix,             &
                              d_matrix,             &
                              need_grid)

  implicit none
  type(basis_set),intent(inout)      :: basis
  type(basis_set),intent(inout)      :: auxil_basis
  real(dp),intent(inout)             :: s_matrix(:,:)
  real(dp),intent(inout)             :: d_matrix(:,:)
  real(dp),intent(in)                :: time_cur
  real(dp),intent(in)                :: dt_factor
  logical,intent(in)                 :: need_eri,need_grid
  !=====
  !=====

  ! Update projectile position and its basis center to t+dt/n
  call change_position_one_atom(ncenter_nuclei,xatom_start(:,ncenter_nuclei) &
       + vel_nuclei(:,ncenter_nuclei) * (time_cur - time_read - time_step*dt_factor))
  call change_basis_center_one_atom(ncenter_basis,xbasis_start(:,ncenter_basis) &
       + vel_nuclei(:,ncenter_nuclei) * (time_cur - time_read - time_step*dt_factor))

  call start_clock(timing_update_basis_eri)
  if( need_eri ) then
    ! Update all basis and eri
    call update_basis_eri(basis,auxil_basis)
  else
    ! Update basis only since eri not needed here
    call moving_basis_set(basis)
    call destroy_libcint(basis)
    call init_libcint(basis)
  endif
  call moving_basis_set(basis_p)
  call destroy_libcint(basis_p)
  call destroy_libcint(basis_t)
  call init_libcint(basis_t, basis_p)
  call init_libcint(basis_p)
  call stop_clock(timing_update_basis_eri)

  call start_clock(timing_update_overlaps)
  ! Update S matrix
  !call setup_overlap(basis,s_matrix)
  call recalc_overlap(basis_t,basis_p,s_matrix)

  ! Analytic evaluation of D(t+dt/n)
  call setup_d_matrix(basis,d_matrix,.TRUE.)
  call stop_clock(timing_update_overlaps)

  ! Update DFT grids for H_xc evaluation
  if( need_grid ) then
    if( calc_type%is_dft ) then
      call destroy_dft_grid()
      call init_dft_grid(basis,tddft_grid_level,dft_xc(1)%needs_gradient,.FALSE.,BATCH_SIZE)
    endif
  endif

end subroutine mb_related_updates


!=========================================================================
subroutine predictor_corrector(basis,                  &
                               auxil_basis,            &
                               c_matrix_cmplx,         &
                               c_matrix_orth_cmplx,    &
                               h_cmplx,                &
                               h_small_cmplx,          &
                               x_matrix,               &
                               s_matrix,               &
                               d_matrix,               &
                               itau,                   &
                               time_cur,               &
                               occupation,             &
                               hamiltonian_kinetic,    &
                               hamiltonian_nucleus,    &
                               dipole_ao)
  implicit none

  type(basis_set),intent(inout)      :: basis
  type(basis_set),intent(inout)      :: auxil_basis
  complex(dp),intent(inout)          :: c_matrix_cmplx(:,:,:)
  complex(dp),intent(inout)          :: c_matrix_orth_cmplx(:,:,:)
  complex(dp),intent(inout)          :: h_cmplx(:,:,:)
  complex(dp),intent(inout)          :: h_small_cmplx(:,:,:)
  real(dp),allocatable,intent(inout) :: x_matrix(:,:)
  real(dp),intent(inout)             :: s_matrix(:,:)
  real(dp),intent(inout)             :: d_matrix(:,:)
  integer,intent(in)                 :: itau
  real(dp),intent(in)                :: time_cur
  real(dp),intent(in)                :: occupation(:,:)
  real(dp),intent(inout)             :: hamiltonian_kinetic(:,:)
  real(dp),intent(inout)             :: hamiltonian_nucleus(:,:)
  real(dp),allocatable,intent(in)    :: dipole_ao(:,:,:)
  !=====
  integer              :: nstate,iextr,i_iter,file_iter_norm
  !=====

  nstate = SIZE(c_matrix_orth_cmplx,DIM=1)

  write(stdout,'(/,1x,a)') 'PREDICTOR-CORRECTOR BLOCK'

  select case (pred_corr)
    ! ///////////////////////////////////
  case('MB_PC0')

    if( excit_type%form == EXCIT_PROJECTILE_W_BASIS ) then
      call mb_related_updates(basis,auxil_basis,.TRUE.,&
             time_cur,0.0_dp,s_matrix,d_matrix,.TRUE.)
    endif

    ! Propagate C(t) -> C(t+dt) using M(t) = S(t)^-1 * ( H(t) - i*D(t) )
    call propagate_nonortho(time_step,s_matrix,d_matrix,c_matrix_cmplx,h_cmplx,prop_type)

    ! Evaluate H(t+dt) using C(t+dt)
    call setup_hamiltonian_cmplx(basis,                  &
                                 nstate,                     &
                                 itau,                       &
                                 time_cur,                   &
                                 time_step,                  &
                                 occupation,                 &
                                 c_matrix_cmplx,             &
                                 hamiltonian_kinetic,        &
                                 hamiltonian_nucleus,        &
                                 h_small_cmplx,              &
                                 x_matrix,                   &
                                 dipole_ao,                  &
                                 h_cmplx,en_tddft)


  case('MB_PC1')

    !--1--PREDICTOR----| H(t-3dt/2),H(t-dt/2)-->H(t+dt/4)
    h_cmplx = -3.0_dp/4.0_dp*h_hist_cmplx(:,:,:,1) + 7.0_dp/4.0_dp*h_hist_cmplx(:,:,:,2)


    !--2--PREDICTOR----| C(t)---U[M(t+dt/4)]--->C(t+dt/2)
    if( excit_type%form == EXCIT_PROJECTILE_W_BASIS ) then
      call mb_related_updates(basis,auxil_basis,.FALSE.,&
             time_cur,3.0_dp/4.0_dp,s_matrix,d_matrix,.FALSE.)
    endif


    ! Propagate C(t) -> C(t+dt/2) using M(t+dt/4) = S(t+dt/4)^-1 * ( H(t+td/4) - i*D(t+dt/4) )
    call propagate_nonortho(time_step/2.0_dp,s_matrix,d_matrix,c_matrix_hist_cmplx(:,:,:,1),h_cmplx,prop_type)

    !--3--CORRECTOR----| C(t+dt/2)-->H(t+dt/2)
    if( excit_type%form == EXCIT_PROJECTILE_W_BASIS ) then
      call mb_related_updates(basis,auxil_basis,.TRUE.,&
             time_cur,1.0_dp/2.0_dp,s_matrix,d_matrix,.TRUE.)
    endif

    ! Evaluate H(t+dt/2)
    call setup_hamiltonian_cmplx(basis,                        &
                                 nstate,                       &
                                 itau,                         &
                                 time_cur-time_step/2.0_dp,    &
                                 time_step,                    &
                                 occupation,                   &
                                 c_matrix_hist_cmplx(:,:,:,1), &
                                 hamiltonian_kinetic,          &
                                 hamiltonian_nucleus,          &
                                 h_small_cmplx,                &
                                 x_matrix,                     &
                                 dipole_ao,                    &
                                 h_cmplx,en_tddft)

    !--4--PROPAGATION----| C(t)---U[M(t+dt/2)]--->C(t+dt)
    call propagate_nonortho(time_step,s_matrix,d_matrix,c_matrix_cmplx,h_cmplx,prop_type)

    !--5--UPDATE----| C(t+dt)-->C(t); H(t-dt/2)-->H(t-3dt/2); H(t+dt/2)-->H(t-dt/2)
    c_matrix_hist_cmplx(:,:,:,1) = c_matrix_cmplx(:,:,:)
    h_hist_cmplx(:,:,:,1) = h_hist_cmplx(:,:,:,2)
    h_hist_cmplx(:,:,:,2) = h_cmplx(:,:,:)


    ! ///////////////////////////////////
  case('MB_PC2B')

    h_cmplx = ( 0.0_dp, 0.0_dp )

    !--0--EXTRAPOLATE----> H(t+dt/4)
    do iextr = 1, n_hist
      h_cmplx = h_cmplx + extrap_coefs(iextr) * h_hist_cmplx(:,:,:,iextr)
    end do
    ! Shift for next iteration : n_hist-2  <---  n_hist; n_hist-3 <-- n_hist-1
    if( n_hist > 2 ) then
      do iextr = 1, n_hist-2
        h_hist_cmplx(:,:,:,iextr) = h_hist_cmplx(:,:,:,iextr+2)
      end do
    end if

    !--1--PROPAGATE----| C(t)--U[M(t+dt/4)]-->C(t+dt/2)
    if( excit_type%form == EXCIT_PROJECTILE_W_BASIS ) then
      call mb_related_updates(basis,auxil_basis,.FALSE.,&
             time_cur,3.0_dp/4.0_dp,s_matrix,d_matrix,.FALSE.)
    endif

    ! Propagate C(t) -> C(t+dt/2) using M(t+dt/4) = S(t+dt/4)^-1 * ( H(t+td/4) - i*D(t+dt/4) )
    call propagate_nonortho(time_step/2.0_dp,s_matrix,d_matrix,c_matrix_hist_cmplx(:,:,:,1),h_cmplx,prop_type)

    !--2--EVALUATE----| C(t+dt/2) --> H(t+dt/2)
    if( excit_type%form == EXCIT_PROJECTILE_W_BASIS ) then
      call mb_related_updates(basis,auxil_basis,.TRUE.,&
             time_cur,1.0_dp/2.0_dp,s_matrix,d_matrix,.TRUE.)
    endif

    ! Calculate H(t+dt/2)
    call setup_hamiltonian_cmplx(basis,                        &
                                 nstate,                       &
                                 itau,                         &
                                 time_cur-time_step/2.0_dp,    &
                                 time_step,                    &
                                 occupation,                   &
                                 c_matrix_hist_cmplx(:,:,:,1), &
                                 hamiltonian_kinetic,          &
                                 hamiltonian_nucleus,          &
                                 h_small_cmplx,                &
                                 x_matrix,                     &
                                 dipole_ao,                    &
                                 h_cmplx,en_tddft)

    ! Save in history : n_hist-1
    if (n_hist > 1) h_hist_cmplx(:,:,:,n_hist-1) = h_cmplx

    !--3--PROPAGATION----| C(t)---U[M(t+dt/2)]--->C(t+dt)
    call propagate_nonortho(time_step,s_matrix,d_matrix,c_matrix_cmplx,h_cmplx,prop_type)

    !--4--EVALUATE----| C(t+dt) --> H(t+dt)
    if( excit_type%form == EXCIT_PROJECTILE_W_BASIS ) then
      call mb_related_updates(basis,auxil_basis,.TRUE.,&
             time_cur,0.0_dp,s_matrix,d_matrix,.TRUE.)
    endif

    ! Calculate H(t+dt)
    call setup_hamiltonian_cmplx(basis,                        &
                                 nstate,                       &
                                 itau,                         &
                                 time_cur,                     &
                                 time_step,                    &
                                 occupation,                   &
                                 c_matrix_cmplx,               &
                                 hamiltonian_kinetic,          &
                                 hamiltonian_nucleus,          &
                                 h_small_cmplx,                &
                                 x_matrix,                     &
                                 dipole_ao,                    &
                                 h_cmplx,en_tddft)

    !--5--UPDATE----| C(t+dt) -> C(t); H(t+dt) -> H(t)
    c_matrix_hist_cmplx(:,:,:,1) = c_matrix_cmplx(:,:,:)
    h_hist_cmplx(:,:,:,n_hist) = h_cmplx


    ! ///////////////////////////////////
  case('PC0')
    call propagate_orth(nstate,basis,time_step,c_matrix_orth_cmplx,c_matrix_cmplx,h_small_cmplx,x_matrix,prop_type)
    call setup_hamiltonian_cmplx(basis,               &
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
    call setup_hamiltonian_cmplx(basis,                 &
                                 nstate,                    &
                                 itau,                      &
                                 time_cur-time_step/2.0_dp, &
                                 time_step,                 &
                                 occupation,                &
                                 c_matrix_cmplx,            &
                                 hamiltonian_kinetic,       &
                                 hamiltonian_nucleus,       &
                                 h_small_cmplx,             &
                                 x_matrix,                  &
                                 dipole_ao,                 &
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

    !close(file_iter_norm)

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

  case default
    call die('Invalid choice for the predictor_corrector scheme. Change pred_corr value in the input file')

  end select

  write(stdout,'(/,1x,a)') 'END OF PREDICTOR-CORRECTOR BLOCK'

end subroutine predictor_corrector


!=========================================================================
subroutine initialize_extrap_coefs(c_matrix_orth_cmplx,h_small_cmplx,c_matrix_cmplx,h_cmplx)
  implicit none

  complex(dp),intent(in)    :: c_matrix_orth_cmplx(:,:,:)
  complex(dp),intent(in)    :: h_small_cmplx(:,:,:)
  complex(dp),intent(in)    :: c_matrix_cmplx(:,:,:)
  complex(dp),intent(in)    :: h_cmplx(:,:,:)
  !=====
  integer               :: iextr,ham_hist_dim,nstate,nbf
  real(dp)              :: x_pred
  real(dp),allocatable  :: m_nodes(:)
  !=====

  nstate = SIZE(c_matrix_orth_cmplx,DIM=1)
  nbf = SIZE(h_cmplx,DIM=1)

  allocate(m_nodes(n_hist),extrap_coefs(n_hist))

  select case (pred_corr)
  case('PC0','MB_PC0')
    continue

  case('PC1','MB_PC1')
    ham_hist_dim = 2

  case('PC2B','MB_PC2B')
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

  case default
    call die('Invalid choice for the predictor_corrector scheme. Change pred_corr value in the input file')

  end select

  if( pred_corr /= 'PC0' .OR. pred_corr /= 'MB_PC0' ) then
    call clean_allocate('h_small_hist_cmplx for TDDFT',h_small_hist_cmplx,nstate,nstate,nspin,ham_hist_dim)
    call clean_allocate('c_matrix_orth_hist_cmplx for TDDFT',c_matrix_orth_hist_cmplx,nstate,nocc,nspin,1)
    call clean_allocate('h_hist_cmplx for TDDFT',h_hist_cmplx,nbf,nbf,nspin,ham_hist_dim)
    call clean_allocate('c_matrix_hist_cmplx for TDDFT',c_matrix_hist_cmplx,nbf,nocc,nspin,1)
    do iextr=1,ham_hist_dim
      h_small_hist_cmplx(:,:,:,iextr)=h_small_cmplx(:,:,:)
      h_hist_cmplx(:,:,:,iextr)=h_cmplx(:,:,:)
    end do
    c_matrix_orth_hist_cmplx(:,:,:,1)=c_matrix_orth_cmplx(:,:,:)
    c_matrix_hist_cmplx(:,:,:,1)=c_matrix_cmplx(:,:,:)
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
     '===================== RT-TDDFT values for time step  ',itau,' ================================='
  write(stdout,'(a31,1x,f19.10)') 'RT-TDDFT Simulation time (au):', time_cur
  write(stdout,'(a31,1x,f19.10)') 'RT-TDDFT Total Energy    (Ha):', en_tddft%total

  select case(excit_type%form)
  case(EXCIT_PROJECTILE, EXCIT_PROJECTILE_W_BASIS)
    write(file_time_data,"(f10.4,12(2x,es16.8e3))") &
       time_cur, en_tddft%total, xatom(:,ncenter_nuclei), en_tddft%nuc_nuc, en_tddft%nucleus, &
       en_tddft%kinetic, en_tddft%hartree, en_tddft%exx_hyb, en_tddft%xc, &
       en_tddft%excit, en_tddft%id
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
subroutine initialize_files(file_time_data,file_dipole_time,file_excit_field,file_mulliken,file_lowdin)
  implicit none
  integer,intent(inout)    :: file_time_data,file_excit_field,file_dipole_time,file_mulliken,file_lowdin
  !=====

  if( .NOT. is_iomaster ) return

  open(newunit=file_time_data,file="time_data.dat")

  if(excit_type%form==EXCIT_LIGHT) then

    open(newunit=file_dipole_time,file="dipole_time.dat")
    open(newunit=file_excit_field,file="excitation_time.dat")

    write(file_excit_field,"(A)") "# time(au)                      E_field_excit_dir(au)"

  end if

  if( print_charge_tddft_ ) then
    !open(newunit=file_mulliken, file="mulliken_charge.dat")
    !write(file_mulliken,"(A)") "##### This is the Mulliken charge file #####"
    !write(file_mulliken,"(A)") "# Time (a.u.)        Re{q_A}         Im{q_A}"

    open(newunit=file_lowdin, file="lowdin_charge.dat")
    write(file_lowdin,"(A)") "##### This is the Lowdin charge file #####"
    write(file_lowdin,"(A)") "# Time (a.u.)    x_proj, y_proj, z_proj (bohr)       q_A"
  end if

  !---------------------------------
  select case(excit_type%form)
  case(EXCIT_PROJECTILE, EXCIT_PROJECTILE_W_BASIS)
    write(file_time_data,"(a10,12(a18))") &
            "# time(au)","e_total","x_proj","y_proj","z_proj","enuc_nuc","enuc_wo_proj","ekin","ehart",&
            "eexx_hyb","exc","enuc_proj","e_iD"

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
      !cur_string = ADJUSTL(cur_string)
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
function check_identity_cmplx(matrix_cmplx) RESULT(is_identity)
  implicit none
  complex(dp),intent(in) :: matrix_cmplx(:,:)
  logical                :: is_identity
  !=====
  real(dp),parameter :: tol=1.0e-9_dp
  integer            ::  imat,jmat,mmat,nmat
  !=====

  mmat = SIZE(matrix_cmplx,DIM=1)
  nmat = SIZE(matrix_cmplx,DIM=1)

  is_identity = .TRUE.
  do jmat=1,nmat
    do imat=1,mmat
      if( imat == jmat ) then
        if( ABS(matrix_cmplx(imat,jmat) - 1.0_dp) > tol ) then
          write(stdout,*) "M(imat,imat)/=1 for: ",imat,jmat,matrix_cmplx(imat,jmat)
          is_identity = .FALSE.
        end if
      else
        if( ABS(matrix_cmplx(imat,jmat)) > tol ) then
          write(stdout,*) "M(imat,jmat)/=0 for: ",imat,jmat,matrix_cmplx(imat,jmat)
          is_identity = .FALSE.
        end if
      end if
    end do
  end do

end function check_identity_cmplx


!=========================================================================
subroutine write_restart_tddft(nstate,time_cur,occupation,c_matrix_tddft)
  implicit none
  integer,intent(in)         :: nstate
  real(dp),intent(in)        :: time_cur
  real(dp),intent(in)        :: occupation(:,:)
  complex(dp),intent(in)     :: c_matrix_tddft(:,:,:)
  !=====
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
  write(restartfile) ncenter_basis
  write(restartfile) zatom(1:ncenter_nuclei)
  write(restartfile) xatom(:,1:ncenter_nuclei)
  write(restartfile) xbasis(:,1:ncenter_basis)
  ! Complex wavefunction coefficients C
  do ispin=1,nspin
    do istate=1,nocc
      write(restartfile) c_matrix_tddft(:,istate,ispin)
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
  !=====
  logical                    :: file_exists
  integer                    :: restartfile
  integer                    :: istate,ispin
  integer                    :: natom_read,nbasis_read
  real(dp)                   :: time_cur_read
  real(dp),allocatable       :: occupation_read(:,:)
  real(dp),allocatable       :: zatom_read(:),x_read(:,:),xbasis_read(:,:)
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

  !ncenter_basis
  read(restartfile) nbasis_read
  if( nbasis_read /= ncenter_basis ) then
    call issue_warning('RESTART_TDDFT file: number of basis has changed.')
    restart_is_correct=.FALSE.
    close(restartfile)
    return
  end if

  allocate(zatom_read(natom_read),x_read(3,natom_read),xbasis_read(3,nbasis_read))
  read(restartfile) zatom_read(1:natom_read)
  read(restartfile) x_read(:,1:natom_read)
  read(restartfile) xbasis_read(:,1:nbasis_read)

  if( ANY( ABS( zatom_read(1:natom_read) - zatom(1:natom_read) ) > 1.0e-5_dp ) &
   .OR. ANY( ABS( x_read(:,1:natom_read-nprojectile) &
                 - xatom(:,1:natom_read-nprojectile) ) > 1.0e-5_dp ) ) then
    call issue_warning('RESTART_TDDFT file: Atom geometry has changed')
  endif

  if( ANY( ABS( xbasis_read(:,1:nbasis_read-nprojectile) &
             - xbasis(:,1:nbasis_read-nprojectile) ) > 1.0e-5_dp ) ) then
    call issue_warning('RESTART_TDDFT file: Basis geometry has changed')
  endif

  deallocate(zatom_read,x_read,xbasis_read)

  ! Here we do not read c_matrix_orth_cmplx from the RESTART_TDDFT file

  write(stdout,*) " RESTART_TDDFT file is correct and will be used for the calculation. SCF loop will be omitted."

  close(restartfile)

end subroutine check_restart_tddft


!=========================================================================
subroutine read_restart_tddft(nstate,time_min,occupation,c_matrix_tddft)
  implicit none
  complex(dp),intent(inout)  :: c_matrix_tddft(:,:,:)
  real(dp),intent(inout)     :: time_min
  real(dp),intent(inout)     :: occupation(:,:)
  integer,intent(in)         :: nstate
  !=====
  logical                    :: file_exists
  integer                    :: restartfile
  integer                    :: istate,ispin
  integer                    :: natom_read, nbasis_read
  real(dp),allocatable       :: occupation_read(:,:)
  real(dp),allocatable       :: zatom_read(:)
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
  read(restartfile) nbasis_read
  allocate(zatom_read(natom_read))
  read(restartfile) zatom_read(1:natom_read)
  read(restartfile) xatom_start(:,1:natom_read)
  read(restartfile) xbasis_start(:,1:nbasis_read)
  deallocate(zatom_read)

  ! Complex wavefunction coefficients C
  do ispin=1,nspin
    do istate=1,nocc
      read(restartfile) c_matrix_tddft(:,istate,ispin)
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
subroutine propagate_nonortho(time_step_cur,s_matrix,d_matrix,c_matrix_cmplx,h_cmplx,prop_type)
  implicit none
  real(dp),intent(in)         :: time_step_cur
  complex(dp),intent(inout)   :: c_matrix_cmplx(:,:,:)
  complex(dp),intent(in)      :: h_cmplx(:,:,:)
  real(dp),intent(in)         :: s_matrix(:,:)
  real(dp),intent(in)         :: d_matrix(:,:)
  character(len=4),intent(in) :: prop_type
  !=====
  integer                    :: ispin
  integer                    :: ibf,nbf,nocc
  !==variables for the CN propagator
  complex(dp),allocatable    :: l_matrix_cmplx(:,:) ! Follow the notation of M.A.L.Marques, C.A.Ullrich et al,
  complex(dp),allocatable    :: b_matrix_cmplx(:,:) ! TDDFT Book, Springer (2006), !p205
  complex(dp),allocatable    :: m_matrix_cmplx(:,:) ! M = S**-1 * ( H - i*D )
  complex(dp),allocatable    :: u_matrix_cmplx(:,:),c_matrix_previous_cmplx(:,:)
  real(dp),allocatable       :: s_matrix_inverse(:,:)
  !=====

  nbf = SIZE(s_matrix,DIM=1)
  nocc = SIZE(c_matrix_cmplx,DIM=2)

  call start_clock(timing_tddft_propagation)

  do ispin=1,nspin

    select case(prop_type)

    case('CN')
      !! C(t+dt) = U(t, t+dt) * C(t)
      !! U = [S + i * dt/2 * ( H - i*D )]^-1 * [S - i * dt/2 * ( H - i*D )]

      allocate(l_matrix_cmplx,MOLD=h_cmplx(:,:,1))
      allocate(b_matrix_cmplx,MOLD=h_cmplx(:,:,1))

      l_matrix_cmplx(:,:) =  im * ( h_cmplx(:,:,ispin) - im*d_matrix(:,:) ) * time_step_cur / 2.0_dp
      b_matrix_cmplx(:,:) = -l_matrix_cmplx(:,:)

      l_matrix_cmplx(:,:) = l_matrix_cmplx(:,:) + s_matrix(:,:)
      b_matrix_cmplx(:,:) = b_matrix_cmplx(:,:) + s_matrix(:,:)

      call start_clock(timing_propagate_inverse)
      call invert(l_matrix_cmplx)
      call stop_clock(timing_propagate_inverse)

      allocate(u_matrix_cmplx,MOLD=h_cmplx(:,:,1))
      call start_clock(timing_propagate_matmul)
      !U_matrix(:,:)              = MATMUL( l_matrix_cmplx(:,:),b_matrix_cmplx(:,:))
      call ZGEMM('N','N',nbf,nbf,nbf,(1.0d0,0.0d0),l_matrix_cmplx,nbf, &
                        b_matrix_cmplx,nbf,(0.0d0,0.0d0),u_matrix_cmplx,nbf)
      deallocate(l_matrix_cmplx)
      deallocate(b_matrix_cmplx)

      allocate(c_matrix_previous_cmplx,SOURCE=c_matrix_cmplx(:,:,1))
      !c_matrix_cmplx(:,:,ispin)  = MATMUL( U_matrix(:,:),c_matrix_cmplx(:,:,ispin))
      call ZGEMM('N','N',nbf,nocc,nbf,(1.0d0,0.0d0),u_matrix_cmplx,nbf, &
                        c_matrix_previous_cmplx,nbf,(0.0d0,0.0d0),c_matrix_cmplx(:,:,ispin),nbf)

      deallocate(u_matrix_cmplx,c_matrix_previous_cmplx)

    case('MAG2')
      allocate(b_matrix_cmplx,MOLD=h_cmplx(:,:,1))
      allocate(m_matrix_cmplx,MOLD=h_cmplx(:,:,1))
      allocate(s_matrix_inverse,MOLD=s_matrix)

      call invert(s_matrix,s_matrix_inverse)
      m_matrix_cmplx(:,:) = MATMUL(s_matrix_inverse(:,:),( h_cmplx(:,:,ispin) - im*d_matrix(:,:) ))
      b_matrix_cmplx(:,:) = (-im) * m_matrix_cmplx(:,:) * time_step_cur - 0.5_dp * &
             (time_step_cur**2) * MATMUL( m_matrix_cmplx(:,:), m_matrix_cmplx(:,:) )

      do ibf=1,SIZE(s_matrix,DIM=1)
        b_matrix_cmplx(ibf,ibf) = b_matrix_cmplx(ibf,ibf) + 1.0_dp
      end do
      !call dump_out_matrix(.TRUE.,'===  U*U**H REAL ===',REAL(MATMUL(b_matrix_cmplx(:,:), CONJG(TRANSPOSE(b_matrix_cmplx(:,:))))))
      !call dump_out_matrix(.TRUE.,'===  U*U**H AIMAG ===',AIMAG(MATMUL(b_matrix_cmplx(:,:), CONJG(TRANSPOSE(b_matrix_cmplx(:,:))))))

      call start_clock(timing_propagate_matmul)
      c_matrix_cmplx(:,:,ispin)  = MATMUL( b_matrix_cmplx(:,:), c_matrix_cmplx(:,:,ispin) )

      deallocate(b_matrix_cmplx)
      deallocate(m_matrix_cmplx)
      deallocate(s_matrix_inverse)

    case default
      call die('Invalid choice for the propagation algorithm. Change prop_type or error_prop_types value in the input file')

    end select

    call stop_clock(timing_propagate_matmul)

  end do

  call stop_clock(timing_tddft_propagation)


end subroutine propagate_nonortho


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
      !call dump_out_matrix(.TRUE.,'===  B/L REAL  ===',REAL(b_matrix_cmplx))
      !call dump_out_matrix(.TRUE.,'===  B/L AIMAG ===',AIMAG(b_matrix_cmplx))
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
  !a_matrix_cmplx(:,1:nstate) = MATMUL( x_matrix(:,:) , a_matrix_cmplx(:,:) )

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

  type(basis_set),intent(inout)   :: basis
  integer,intent(in)              :: nstate
  integer,intent(in)              :: itau
  real(dp),intent(in)             :: time_cur
  real(dp),intent(in)             :: time_step_cur
  real(dp),intent(in)             :: occupation(nstate,nspin)
  real(dp),intent(inout)          :: hamiltonian_kinetic(basis%nbf,basis%nbf)
  real(dp),intent(inout)          :: hamiltonian_nucleus(basis%nbf,basis%nbf)
  real(dp),allocatable,intent(in) :: dipole_ao(:,:,:)
  real(dp),intent(in)             :: x_matrix(basis%nbf,nstate)
  complex(dp),intent(in)          :: c_matrix_cmplx(basis%nbf,nocc,nspin)
  complex(dp),intent(out)         :: hamiltonian_cmplx(basis%nbf,basis%nbf,nspin)
  complex(dp),intent(out)         :: h_small_cmplx(nstate,nstate,nspin)
  type(energy_contributions),intent(inout) :: en
  !=====
  logical              :: calc_excit_
  integer              :: ispin, idir, iatom
  real(dp)             :: excit_field(3)
  complex(dp)          :: p_matrix_cmplx(basis%nbf,basis%nbf,nspin)
  integer              :: projectile_list(1),fixed_atom_list(ncenter_nuclei-nprojectile)
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

    !
    ! Projectile excitation with moving basis
  case(EXCIT_PROJECTILE_W_BASIS)

    if ( itau > 0 ) then
      !call setup_kinetic(basis,hamiltonian_kinetic)
      call recalc_kinetic(basis_t,basis_p,hamiltonian_kinetic)
      call nucleus_nucleus_energy(en_tddft%nuc_nuc)
      ! Nucleus-electron interaction due to the fixed target
      call recalc_nucleus(basis_t,basis_p,hamiltonian_nucleus)
      !do iatom=1,ncenter_nuclei-nprojectile
      !  fixed_atom_list(iatom) = iatom
      !enddo
      !call setup_nucleus(basis,hamiltonian_nucleus,fixed_atom_list)
    end if

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

  if( excit_type%form /= EXCIT_PROJECTILE_W_BASIS ) then
    ! Perform the canonical transform from the original basis to the orthogonal one
    call transform_hamiltonian_ortho(x_matrix,hamiltonian_cmplx,h_small_cmplx)
  end if


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
                                         COMPLEX_ZERO,h_small_cmplx(:,:,ispin),nstate)
      call matrix_lower_to_full(h_small_cmplx(:,:,ispin))
#else
      call ZGEMM('C','N',nstate,nstate,nbf,COMPLEX_ONE,x_matrix_cmplx(1,1),nbf, &
                                                       m_tmp(1,1),nbf,  &
                                           COMPLEX_ZERO,h_small_cmplx(:,:,ispin),nstate)
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
    excit_field(:) = excit_type%kappa * EXP( -( time_cur-excit_type%time0 )**2 / 2.0_dp / excit_omega**2 ) &
                     * excit_dir_norm(:)
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
