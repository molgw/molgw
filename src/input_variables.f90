!======================================================================
! This file is part of MOLGW.
!
! The following lines have been generated by the python script input_variables.py
! Do not alter them directly: they will be overriden sooner or later by the script
! To add a new input variable, modify the script directly
! Generated by input_variables.py on 09 January 2024
!======================================================================

  namelist /molgw/   &
    auto_auxil_fsam,       &
    auto_auxil_lmaxinc,       &
    auto_occupation,       &
    auxil_basis,       &
    basis,       &
    comment,       &
    ecp_auxil_basis,       &
    ecp_basis,       &
    ecp_elements,       &
    ecp_quality,       &
    ecp_type,       &
    eri3_genuine,       &
    even_tempered_alpha,       &
    even_tempered_beta,       &
    even_tempered_n_list,       &
    gaussian_type,       &
    incore,       &
    memory_evaluation,       &
    move_nuclei,       &
    nstep,       &
    scf,       &
    tolforce,       &
    eri3_nbatch,       &
    eri3_npcol,       &
    eri3_nprow,       &
    grid_memory,       &
    mpi_nproc_ortho,       &
    scalapack_block_min,       &
    basis_path,       &
    cube_nx,       &
    cube_ny,       &
    cube_nz,       &
    cube_state_min,       &
    cube_state_max,       &
    force_energy_qp,       &
    ignore_bigrestart,       &
    print_bigrestart,       &
    print_cube,       &
    print_transition_density,       &
    print_wfn_files,       &
    print_all_MO_wfn_file,       &
    print_density_matrix,       &
    print_eri,       &
    print_hartree,       &
    print_multipole,       &
    print_pdos,       &
    print_restart,       &
    print_rho_grid,       &
    print_sigma,       &
    print_spatial_extension,       &
    print_w,       &
    print_wfn,       &
    print_yaml,       &
    read_fchk,       &
    read_restart,       &
    calc_dens_disc,       &
    calc_q_matrix,       &
    calc_spectrum,       &
    print_cube_diff_tddft,       &
    print_cube_rho_tddft,       &
    print_c_matrix_cmplx_hdf5,       &
    print_p_matrix_cmplx_hdf5,       &
    print_dens_traj,       &
    print_dens_traj_points_set,       &
    print_dens_traj_tddft,       &
    print_line_rho_diff_tddft,       &
    print_line_rho_tddft,       &
    print_tddft_matrices,       &
    print_charge_tddft,       &
    print_tddft_restart,       &
    read_tddft_restart,       &
    write_step,       &
    calc_charge_step,       &
    analytic_chi,       &
    assume_scf_converged,       &
    acfd_nlambda,       &
    ci_greens_function,       &
    ci_nstate,       &
    ci_nstate_self,       &
    ci_spin_multiplicity,       &
    ci_type,       &
    cphf_cpks_0,       &
    dft_core,       &
    ecp_small_basis,       &
    eta,       &
    frozencore,       &
    gwgwg_skip_vvv,       &
    gwgwg_skip_vv,       &
    gwgwg_static_approximation,       &
    gwgamma_tddft,       &
    mu_origin,       &
    ncoreg,       &
    ncorew,       &
    nexcitation,       &
    nomega_chi_imag,       &
    nomega_chi_real,       &
    nomega_sigma,       &
    nomega_sigma_calc,       &
    nstep_dav,       &
    nstep_gw,       &
    nvel_projectile,       &
    nvirtualg,       &
    nvirtualw,       &
    postscf,       &
    postscf_diago_flavor,       &
    pt3_a_diagrams,       &
    pt_density_matrix,       &
    rcut_mbpt,       &
    scissor,       &
    selfenergy_state_max,       &
    selfenergy_state_min,       &
    selfenergy_state_range,       &
    small_basis,       &
    step_sigma,       &
    step_sigma_calc,       &
    stopping,       &
    stopping_nq,       &
    stopping_dq,       &
    tda,       &
    tddft_grid_quality,       &
    toldav,       &
    triplet,       &
    use_correlated_density_matrix,       &
    virtual_fno,       &
    w_screening,       &
    excit_dir,       &
    excit_kappa,       &
    excit_name,       &
    excit_omega,       &
    excit_time0,       &
    n_hist,       &
    n_iter,       &
    n_restart_tddft,       &
    ncore_tddft,       &
    pred_corr,       &
    prop_type,       &
    projectile_charge_scaling,       &
    r_disc,       &
    tddft_frozencore,       &
    tddft_wfn_t0,       &
    tddft_energy_shift,       &
    tddft_charge,       &
    time_sim,       &
    time_step,       &
    vel_projectile,       &
    tolscf_tddft,       &
    alpha_hybrid,       &
    alpha_mixing,       &
    beta_hybrid,       &
    noft_complex,       &
    noft_dft,       &
    density_matrix_damping,       &
    diis_switch,       &
    gamma_hybrid,       &
    kappa_hybrid,       &
    grid_quality,       &
    init_hamiltonian,       &
    integral_quality,       &
    kerker_k0,       &
    level_shifting_energy,       &
    min_overlap,       &
    mixing_scheme,       &
    tolscf,       &
    npulay_hist,       &
    nscf,       &
    partition_scheme,       &
    scf_diago_flavor,       &
    noft_rsintra,       &
    noft_lowmemERI,       &
    noft_fcidump,       &
    noft_NOTupdateOCC,       &
    noft_NOTupdateORB,       &
    noft_functional,       &
    noft_printdmn,       &
    noft_printswdmn,       &
    noft_printints,       &
    noft_readCOEF,       &
    noft_readFdiag,       &
    noft_readGAMMAS,       &
    noft_readOCC,       &
    noft_sta,       &
    noft_ithresh_lambda,       &
    noft_Lpower,       &
    noft_npairs,       &
    noft_ncoupled,       &
    noft_ndiis,       &
    noft_nscf,       &
    noft_restart,       &
    noft_tolE,       &
    charge,       &
    electric_field_x,       &
    electric_field_y,       &
    electric_field_z,       &
    length_unit,       &
    magnetization,       &
    natom,       &
    nghost,       &
    nspin,       &
    temperature,       &
    xyz_file

!=====

 auto_auxil_fsam=1.5_dp 
 auto_auxil_lmaxinc=1
 auto_occupation='no'
 auxil_basis=''
 basis=''
 comment=''
 ecp_auxil_basis=''
 ecp_basis=''
 ecp_elements=''
 ecp_quality='high'
 ecp_type=''
 eri3_genuine='no'
 even_tempered_alpha=1.0_dp 
 even_tempered_beta=0.5_dp 
 even_tempered_n_list='1'
 gaussian_type='pure'
 incore='yes'
 memory_evaluation='no'
 move_nuclei='no'
 nstep=50
 scf=''
 tolforce=1e-05_dp 
 eri3_nbatch=1
 eri3_npcol=1
 eri3_nprow=1
 grid_memory=400.0_dp 
 mpi_nproc_ortho=1
 scalapack_block_min=100000
 basis_path=''
 cube_nx=30
 cube_ny=30
 cube_nz=30
 cube_state_min=1
 cube_state_max=1
 force_energy_qp='no'
 ignore_bigrestart='no'
 print_bigrestart='yes'
 print_cube='no'
 print_transition_density='no'
 print_wfn_files='no'
 print_all_MO_wfn_file='no'
 print_density_matrix='no'
 print_eri='no'
 print_hartree='no'
 print_multipole='no'
 print_pdos='no'
 print_restart='yes'
 print_rho_grid='no'
 print_sigma='no'
 print_spatial_extension='no'
 print_w='no'
 print_wfn='no'
 print_yaml='yes'
 read_fchk='no'
 read_restart='no'
 calc_dens_disc='no'
 calc_q_matrix='no'
 calc_spectrum='no'
 print_cube_diff_tddft='no'
 print_cube_rho_tddft='no'
 print_c_matrix_cmplx_hdf5='no'
 print_p_matrix_cmplx_hdf5='no'
 print_dens_traj='no'
 print_dens_traj_points_set='no'
 print_dens_traj_tddft='no'
 print_line_rho_diff_tddft='no'
 print_line_rho_tddft='no'
 print_tddft_matrices='no'
 print_charge_tddft='no'
 print_tddft_restart='yes'
 read_tddft_restart='no'
 write_step=1_dp 
 calc_charge_step=1_dp 
 analytic_chi='no'
 assume_scf_converged='no'
 acfd_nlambda=21
 ci_greens_function='holes'
 ci_nstate=1
 ci_nstate_self=1
 ci_spin_multiplicity=1
 ci_type='all'
 cphf_cpks_0='no'
 dft_core=0
 ecp_small_basis=''
 eta=0.001_dp 
 frozencore='no'
 gwgwg_skip_vvv='no'
 gwgwg_skip_vv='no'
 gwgwg_static_approximation='no'
 gwgamma_tddft='no'
 mu_origin=-100.0_dp 
 ncoreg=0
 ncorew=0
 nexcitation=0
 nomega_chi_imag=0
 nomega_chi_real=2
 nomega_sigma=51
 nomega_sigma_calc=1
 nstep_dav=15
 nstep_gw=1
 nvel_projectile=1
 nvirtualg=100000
 nvirtualw=100000
 postscf=''
 postscf_diago_flavor=' '
 pt3_a_diagrams='yes'
 pt_density_matrix='no'
 rcut_mbpt=1.0_dp 
 scissor=0.0_dp 
 selfenergy_state_max=100000
 selfenergy_state_min=1
 selfenergy_state_range=100000
 small_basis=''
 step_sigma=0.01_dp 
 step_sigma_calc=0.03_dp 
 stopping='no'
 stopping_nq=500
 stopping_dq=0.02_dp 
 tda='no'
 tddft_grid_quality='high'
 toldav=0.0001_dp 
 triplet='no'
 use_correlated_density_matrix='no'
 virtual_fno='no'
 w_screening='rpa'
 excit_dir=(/ 1.0_dp , 0.0_dp , 0.0_dp /)
 excit_kappa=2e-05_dp 
 excit_name='no'
 excit_omega=0.2_dp 
 excit_time0=3.0_dp 
 n_hist=2
 n_iter=2
 n_restart_tddft=50
 ncore_tddft=0
 pred_corr='PC2B'
 prop_type='CN'
 projectile_charge_scaling=1.0_dp 
 r_disc=200.0_dp 
 tddft_frozencore='no'
 tddft_wfn_t0='SCF'
 tddft_energy_shift=0.0_dp 
 tddft_charge=-999.0_dp 
 time_sim=10.0_dp 
 time_step=1.0_dp 
 vel_projectile=(/ 0.0_dp , 0.0_dp , 1.0_dp /)
 tolscf_tddft=0.0001_dp 
 alpha_hybrid=0.0_dp 
 alpha_mixing=0.7_dp 
 beta_hybrid=0.0_dp 
 noft_complex='no'
 noft_dft='no'
 density_matrix_damping=0.0_dp 
 diis_switch=0.05_dp 
 gamma_hybrid=1000000.0_dp 
 kappa_hybrid=0.0_dp 
 grid_quality='high'
 init_hamiltonian='guess'
 integral_quality='high'
 kerker_k0=0.0_dp 
 level_shifting_energy=0.0_dp 
 min_overlap=1e-05_dp 
 mixing_scheme='pulay'
 tolscf=1e-07_dp 
 npulay_hist=6
 nscf=50
 partition_scheme='ssf'
 scf_diago_flavor=' '
 noft_rsintra='yes'
 noft_lowmemERI='yes'
 noft_fcidump='no'
 noft_NOTupdateOCC='no'
 noft_NOTupdateORB='no'
 noft_functional='GNOF'
 noft_printdmn='no'
 noft_printswdmn='no'
 noft_printints='no'
 noft_readCOEF='no'
 noft_readFdiag='no'
 noft_readGAMMAS='no'
 noft_readOCC='no'
 noft_sta='no'
 noft_ithresh_lambda=5
 noft_Lpower=0.53_dp 
 noft_npairs=1
 noft_ncoupled=2
 noft_ndiis=5
 noft_nscf=1000
 noft_restart='no'
 noft_tolE=1e-09_dp 
 charge=0.0_dp 
 electric_field_x=0.0_dp 
 electric_field_y=0.0_dp 
 electric_field_z=0.0_dp 
 length_unit='angstrom'
 magnetization=0.0_dp 
 natom=0
 nghost=0
 nspin=1
 temperature=0.0_dp 
 xyz_file=''


!======================================================================
