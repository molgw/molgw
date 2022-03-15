# Input variable list

---

## Mandatory input variables 

[basis](#basis) 
[scf](#scf) 


## Physical system setup input variables 

[charge](#charge) 
[length_unit](#length_unit) 
[magnetization](#magnetization) 
[natom](#natom) 
[nghost](#nghost) 
[nspin](#nspin) 
[temperature](#temperature) 
[xyz_file](#xyz_file) 


## General input variables 

[auto_auxil_fsam](#auto_auxil_fsam) 
[auto_auxil_lmaxinc](#auto_auxil_lmaxinc) 
[auxil_basis](#auxil_basis) 
[basis](#basis) 
[comment](#comment) 
[ecp_auxil_basis](#ecp_auxil_basis) 
[ecp_basis](#ecp_basis) 
[ecp_elements](#ecp_elements) 
[ecp_quality](#ecp_quality) 
[ecp_type](#ecp_type) 
[eri3_genuine](#eri3_genuine) 
[gaussian_type](#gaussian_type) 
[incore](#incore) 
[memory_evaluation](#memory_evaluation) 
[move_nuclei](#move_nuclei) 
[nstep](#nstep) 
[scf](#scf) 
[tolforce](#tolforce) 


## Self-consistency input variables 

[alpha_hybrid](#alpha_hybrid) 
[alpha_mixing](#alpha_mixing) 
[beta_hybrid](#beta_hybrid) 
[noft_complex](#noft_complex) 
[density_matrix_damping](#density_matrix_damping) 
[diis_switch](#diis_switch) 
[noft_fcidump](#noft_fcidump) 
[gamma_hybrid](#gamma_hybrid) 
[grid_quality](#grid_quality) 
[init_hamiltonian](#init_hamiltonian) 
[noft_NOTupdateOCC](#noft_NOTupdateOCC) 
[noft_NOTupdateORB](#noft_NOTupdateORB) 
[noft_functional](#noft_functional) 
[integral_quality](#integral_quality) 
[noft_printdmn](#noft_printdmn) 
[noft_printswdmn](#noft_printswdmn) 
[noft_printints](#noft_printints) 
[noft_readCOEF](#noft_readCOEF) 
[noft_readFdiag](#noft_readFdiag) 
[noft_readGAMMAS](#noft_readGAMMAS) 
[noft_readOCC](#noft_readOCC) 
[noft_sta](#noft_sta) 
[noft_ithresh_lambda](#noft_ithresh_lambda) 
[kerker_k0](#kerker_k0) 
[level_shifting_energy](#level_shifting_energy) 
[noft_lowmemERI](#noft_lowmemERI) 
[noft_Lpower](#noft_Lpower) 
[min_overlap](#min_overlap) 
[mixing_scheme](#mixing_scheme) 
[noft_npairs](#noft_npairs) 
[noft_ncoupled](#noft_ncoupled) 
[npulay_hist](#npulay_hist) 
[nscf](#nscf) 
[noft_ndiis](#noft_ndiis) 
[noft_nscf](#noft_nscf) 
[partition_scheme](#partition_scheme) 
[noft_restart](#noft_restart) 
[scf_diago_flavor](#scf_diago_flavor) 
[noft_tolE](#noft_tolE) 
[tolscf](#tolscf) 


## Correlation and excited states post-treatment input variables 

[assume_scf_converged](#assume_scf_converged) 
[ci_greens_function](#ci_greens_function) 
[ci_nstate](#ci_nstate) 
[ci_nstate_self](#ci_nstate_self) 
[ci_spin_multiplicity](#ci_spin_multiplicity) 
[ci_type](#ci_type) 
[dft_core](#dft_core) 
[ecp_small_basis](#ecp_small_basis) 
[eta](#eta) 
[frozencore](#frozencore) 
[gwgamma_tddft](#gwgamma_tddft) 
[ncoreg](#ncoreg) 
[ncorew](#ncorew) 
[nexcitation](#nexcitation) 
[nomega_chi_imag](#nomega_chi_imag) 
[nomega_sigma](#nomega_sigma) 
[nomega_sigma_calc](#nomega_sigma_calc) 
[nstep_dav](#nstep_dav) 
[nstep_gw](#nstep_gw) 
[nvirtualg](#nvirtualg) 
[nvirtualw](#nvirtualw) 
[postscf](#postscf) 
[postscf_diago_flavor](#postscf_diago_flavor) 
[pt3_a_diagrams](#pt3_a_diagrams) 
[pt_density_matrix](#pt_density_matrix) 
[rcut_mbpt](#rcut_mbpt) 
[scissor](#scissor) 
[selfenergy_state_max](#selfenergy_state_max) 
[selfenergy_state_min](#selfenergy_state_min) 
[selfenergy_state_range](#selfenergy_state_range) 
[small_basis](#small_basis) 
[step_sigma](#step_sigma) 
[step_sigma_calc](#step_sigma_calc) 
[stopping](#stopping) 
[tda](#tda) 
[tddft_grid_quality](#tddft_grid_quality) 
[toldav](#toldav) 
[triplet](#triplet) 
[use_correlated_density_matrix](#use_correlated_density_matrix) 
[virtual_fno](#virtual_fno) 


## IO input variables 

[basis_path](#basis_path) 
[force_energy_qp](#force_energy_qp) 
[ignore_bigrestart](#ignore_bigrestart) 
[print_bigrestart](#print_bigrestart) 
[print_cube](#print_cube) 
[print_wfn_files](#print_wfn_files) 
[print_density_matrix](#print_density_matrix) 
[print_eri](#print_eri) 
[print_hartree](#print_hartree) 
[print_multipole](#print_multipole) 
[print_pdos](#print_pdos) 
[print_restart](#print_restart) 
[print_rho_grid](#print_rho_grid) 
[print_sigma](#print_sigma) 
[print_spatial_extension](#print_spatial_extension) 
[print_w](#print_w) 
[print_wfn](#print_wfn) 
[print_yaml](#print_yaml) 
[read_fchk](#read_fchk) 
[read_restart](#read_restart) 


## Hardware input variables 

[eri3_nbatch](#eri3_nbatch) 
[eri3_npcol](#eri3_npcol) 
[eri3_nprow](#eri3_nprow) 
[grid_memory](#grid_memory) 
[mpi_nproc_ortho](#mpi_nproc_ortho) 
[scalapack_block_min](#scalapack_block_min) 


## Real time TDDFT 

[excit_dir](#excit_dir) 
[excit_kappa](#excit_kappa) 
[excit_name](#excit_name) 
[excit_omega](#excit_omega) 
[excit_time0](#excit_time0) 
[n_hist](#n_hist) 
[n_iter](#n_iter) 
[n_restart_tddft](#n_restart_tddft) 
[ncore_tddft](#ncore_tddft) 
[pred_corr](#pred_corr) 
[prop_type](#prop_type) 
[projectile_charge_scaling](#projectile_charge_scaling) 
[r_disc](#r_disc) 
[tddft_frozencore](#tddft_frozencore) 
[time_sim](#time_sim) 
[time_step](#time_step) 
[vel_projectile](#vel_projectile) 


## IO Real time TDDFT 

[calc_dens_disc](#calc_dens_disc) 
[calc_q_matrix](#calc_q_matrix) 
[calc_spectrum](#calc_spectrum) 
[print_cube_diff_tddft](#print_cube_diff_tddft) 
[print_cube_rho_tddft](#print_cube_rho_tddft) 
[print_dens_traj](#print_dens_traj) 
[print_dens_traj_points_set](#print_dens_traj_points_set) 
[print_dens_traj_tddft](#print_dens_traj_tddft) 
[print_line_rho_diff_tddft](#print_line_rho_diff_tddft) 
[print_line_rho_tddft](#print_line_rho_tddft) 
[print_tddft_matrices](#print_tddft_matrices) 
[print_tddft_restart](#print_tddft_restart) 
[read_tddft_restart](#read_tddft_restart) 
[write_step](#write_step) 


---

## Complete list of input variables 

---
### alpha_hybrid

*Optional* 

**Family:** scf 

**Type:** real 

**Default:** 0.0 

**Description:** 

Only works for Range-Separated hybrid functionals scf='rsh' Sets the amount of range-independent exact-exchange 


---
### alpha_mixing

*Optional* 

**Family:** scf 

**Type:** real 

**Default:** 0.7 

**Description:** 

Sets the amount of output density-matrix for the next iteration. When the SCF cycles have difficulties to converge, one may try to lower this value. 


---
### assume_scf_converged

*Optional* 

**Family:** post 

**Type:** yes/no 

**Default:** no 

**Description:** 

Allows for a post-scf calculation whatever the outcome of the SCF loop. Especially useful when restarting from another SCF calculations with nscf=0 


---
### auto_auxil_fsam

*Optional* 

**Family:** general 

**Type:** real 

**Default:** 1.5 

**Description:** 

Sets the F_SAM parameter in the automatic generation of the auxiliary basis set. The closer to 1.0 the more auxiliary basis functions it will generate. See Yang-Rendell-Frisch, JChemPhys 127, 074102 (2007) for more details. 


---
### auto_auxil_lmaxinc

*Optional* 

**Family:** general 

**Type:** integer 

**Default:** 1 

**Description:** 

Sets the l_MAXINC parameter in the automatic generation of the auxiliary basis set. The larger the more auxiliary basis functions it will generate. See Yang-Rendell-Frisch, JChemPhys 127, 074102 (2007) for more details. 


---
### auxil_basis

*Optional* 

**Family:** general 

**Type:** characters 

**Default:** None 

**Description:** 

Sets the auxiliary basis set. For instance, cc-pVDZ-RI for a Weigend basis set. If present, the auxiliary basis will be used for both the scf cycles and the postscf calculations (TD-DFT, BSE, or GW). 


---
### basis

*Mandatory* 

**Family:** general 

**Type:** characters 

**Default:** None 

**Description:** 

Sets the basis set For Pople sets, use 6-31G for instance or 6-31+G\*. For Dunning sets, use aug-cc-pVTZ for instance. Note that Pople sets are to be used with gaussian_type='cart' One may use ones own basis sets provided that the files are labeled X_mybasisset where X is the element. 


---
### basis_path

*Optional* 

**Family:** io 

**Type:** characters 

**Default:** None 

**Description:** 

Sets the path pointing to the basis functions files. If not specified, then the basis set files will be searched in folder ~molgw/basis/. 


---
### beta_hybrid

*Optional* 

**Family:** scf 

**Type:** real 

**Default:** 0.0 

**Description:** 

Only works for Range-Separated hybrid functionals scf='rsh' Sets the amount of long-range exact-exchange 


---
### calc_dens_disc

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Calculate electronic density in the discs during the real-time dynamics 


---
### calc_q_matrix

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Calculate and print q_matrix which is the projection of a propagated state psi(t) onto the initial state psi(0) in the real-time dynamics 


---
### calc_spectrum

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Calculates absorption spectrum in the real-time dynamics 


---
### charge

*Optional* 

**Family:** system 

**Type:** real 

**Default:** 0.0 

**Description:** 

Sets the total charge of the system. 0 is a neutral system. -2 is a doubly charged anion etc. 


---
### ci_greens_function

**experimental** 

*Optional* 

**Family:** post 

**Type:** characters 

**Default:** holes 

**Description:** 

EXPERIMENTAL. Selects which part of the Green's function is to be calculated: holes, electrons, or both. 


---
### ci_nstate

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 1 

**Description:** 

Selects how many CI states should be calculated in the diagonalization. If ci_nstate is lower than the number of configuration,  a Davidson partial diagonalization is performed, else a full (SCA)LAPACK diagonalization is triggered. 


---
### ci_nstate_self

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 1 

**Description:** 

Selects how many CI states in the N+1 or N-1 electron calculations. If ci_nstate_self is lower than the number of configuration,  a Davidson partial diagonalization is performed, else a full (SCA)LAPACK diagonalization is triggered. 


---
### ci_spin_multiplicity

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 1 

**Description:** 

Spin multiplicity in CI calculations. 


---
### ci_type

*Optional* 

**Family:** post 

**Type:** characters 

**Default:** all 

**Description:** 

Selects which excitations will be included in the CI expansion. Valid choices are 'all', 'CISD', 'CISDT', 'CISDTQ'. 


---
### comment

*Optional* 

**Family:** general 

**Type:** characters 

**Default:** None 

**Description:** 

This is a free expression place. Use it as you wish for commenting, naming, labeling etc. (140 character max just as twitter) 


---
### density_matrix_damping

*Optional* 

**Family:** scf 

**Type:** real 

**Default:** 0.0 

**Description:** 

Adds an additional linear mixing on the density matrix in combination with the Hamiltonian mixing in order to damp out the charge oscillations. Especially useful for metallic systems. 


---
### dft_core

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 0 

**Description:** 

Sets the number of states considered as core in &lt;&Sigma;<sub>x</sub>-<i>v</i><sub>xc</sub>&gt. This options is meant to mimic the pseudopotential approximation. 


---
### diis_switch

*Optional* 

**Family:** scf 

**Type:** real 

**Default:** 0.05 

**Description:** 

When running ADIIS, sets the residue value below which the DIIS method is used to finalize the convergence. 


---
### ecp_auxil_basis

*Optional* 

**Family:** general 

**Type:** characters 

**Default:** None 

**Description:** 

Name of the auxiliary basis set to be used for elements specified in list ecp_elements. 


---
### ecp_basis

*Optional* 

**Family:** general 

**Type:** characters 

**Default:** None 

**Description:** 

Name of the basis set to be used for elements specified in list ecp_elements. 


---
### ecp_elements

*Optional* 

**Family:** general 

**Type:** characters 

**Default:** None 

**Description:** 

Contains the list of elements (separated by spaces) that should be treated with an Effective Core Potential. 


---
### ecp_quality

*Optional* 

**Family:** general 

**Type:** characters 

**Default:** high 

**Description:** 

Sets the number of grid points use to evaluate the Effective Core Potential integrals in real space. Possible values are 'low', 'medium', 'high', 'very high', 'insane'. It could be abbreviated in 'l', 'm', 'h', 'vh', 'i'. 'high' is usually fine. 'insane' is only meant for debugging since it is overdoing a lot. 


---
### ecp_small_basis

**experimental** 

*Optional* 

**Family:** post 

**Type:** characters 

**Default:** None 

**Description:** 

Calls for a smaller basis set used to represent the virtual orbital space with fewer functions. This is the small basis set used for elements with an effective core potential. Only meaningful for GW. 


---
### ecp_type

*Optional* 

**Family:** general 

**Type:** characters 

**Default:** None 

**Description:** 

Name of the Effective Core Potential. For instance, Gold using the cc-pVDZ-PP basis set should have ecp_type='PP', so that MOLGW looks for the file Au_PP in the basis_path folder. 


---
### eri3_genuine

*Optional* 

**Family:** general 

**Type:** yes/no 

**Default:** no 

**Description:** 

If set to 'yes', the 2-center integrals will not be diagonalized and 3-center electron repulsion integrals will remain the genuine 3-center integrals. 


---
### eri3_nbatch

*Optional* 

**Family:** hardware 

**Type:** integer 

**Default:** 1 

**Description:** 

Sets the number of batches when calculating the 3-center integrals. Having a large eri3_nbatch reduces the memory foot print, however it may lower the performance. 


---
### eri3_npcol

*Optional* 

**Family:** hardware 

**Type:** integer 

**Default:** 1 

**Description:** 

Sets number of column processors for the distribution of the 3-center integrals.  eri3_nprow X eri3_npcol must be equal to the number of MPI threads else MOLGW decides on its own. 


---
### eri3_nprow

*Optional* 

**Family:** hardware 

**Type:** integer 

**Default:** 1 

**Description:** 

Sets number of row processors for the distribution of the 3-center integrals.  eri3_nprow X eri3_npcol must be equal to the number of MPI threads else MOLGW decides on its own. 


---
### eta

*Optional* 

**Family:** post 

**Type:** real 

**Default:** 0.001 

**Description:** 

Is a the tiny imaginary part used in the denominator of the Green's function to shift the pole off the axis, so to avoid divergences.This is an energy in Hartree. It should be set to the lowest value possible in theory. However, in practice, a too low value of eta would induce huge and unstable GW corrections. The default value is usually very accurate and there is no need to use a lower value. But for states apart from the band gap, a large value of eta may be beneficial for stability. eta=0.01 is already much more stable. Note that for QSGW increasing eta is most often unavoidable. 


---
### excit_dir

*Optional* 

**Family:** rt_tddft 

**Type:** vector_1d_3 

**Default:** (1.0, 0.0, 0.0) 

**Description:** 

Excitation direction for the real-time dynamics. 


---
### excit_kappa

*Optional* 

**Family:** rt_tddft 

**Type:** real 

**Default:** 2e-05 

**Description:** 

Maximum Gaussian excitation field strength in atomic units. 


---
### excit_name

*Optional* 

**Family:** rt_tddft 

**Type:** characters 

**Default:** no 

**Description:** 

Sets the type of excitation of a system in the real-time dynamics.  'GAU' stands for a linearly polarized uniform Gaussian electric field 


---
### excit_omega

*Optional* 

**Family:** rt_tddft 

**Type:** real 

**Default:** 0.2 

**Description:** 

The excitation pulse width in atomic units for the real-time dynamics. 


---
### excit_time0

*Optional* 

**Family:** rt_tddft 

**Type:** real 

**Default:** 3.0 

**Description:** 

Center of the excitation pulse in atomic units for the real-time dynamics. 


---
### force_energy_qp

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Force the reading of the ENERGY_QP file whatever the postscf choice. 


---
### frozencore

*Optional* 

**Family:** post 

**Type:** yes/no 

**Default:** no 

**Description:** 

Triggers the neglect of core states in GW. H, He, Li, Be have no core states. B-Na have the 1s. Al-Ca have the 1s2s2p. Manual tuning could be achieved with ncoreg, ncorew. 


---
### gamma_hybrid

*Optional* 

**Family:** scf 

**Type:** real 

**Default:** 1000000.0 

**Description:** 

Only works for Range-Separated hybrid functionals scf='rsh' Sets the separation between long-range and short-range. It is input in bohr^-1. 


---
### gaussian_type

*Optional* 

**Family:** general 

**Type:** characters 

**Default:** pure 

**Description:** 

Asks for pure or spherical Gaussian type orbitals with 'pure' or for Cartesian Gaussian orbital with 'cart'. 


---
### grid_memory

*Optional* 

**Family:** hardware 

**Type:** real 

**Default:** 400.0 

**Description:** 

Sets the maximum memory usage in Mb allowed to store the wavefunctions on the quadrature points for XC integrals. 


---
### grid_quality

*Optional* 

**Family:** scf 

**Type:** characters 

**Default:** high 

**Description:** 

Sets the number of grid points use to evaluate the exchange-correlation integrals in real space for the DFT potential and energy. Possible values are 'low', 'medium', 'high', 'very high', 'insane'. It could be abbreviated in 'l', 'm', 'h', 'vh', 'i'. 'high' is usually fine. 'insane' is only meant for debugging since it is overdoing a lot. 


---
### gwgamma_tddft

**experimental** 

*Optional* 

**Family:** post 

**Type:** yes/no 

**Default:** no 

**Description:** 

EXPERIMENTAL. Calculates the vertex using the DFT flavor specified in the ground-state calculation. 


---
### ignore_bigrestart

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Considers a big RESTART as if it was a small RESTART. 


---
### incore

*Optional* 

**Family:** general 

**Type:** yes/no 

**Default:** yes 

**Description:** 

Specify if the 4-center integrals are all calculated at once and stored or if they are calculated on-the-fly. 


---
### init_hamiltonian

*Optional* 

**Family:** scf 

**Type:** characters 

**Default:** guess 

**Description:** 

Selects how to initiate the first hamiltonian for SCF cycles. Today, two options are available: 'guess' for an educated guess based on approximate atomic densities or 'core' for the core hamiltonian. 


---
### integral_quality

*Optional* 

**Family:** scf 

**Type:** characters 

**Default:** high 

**Description:** 

Sets the tolerance value for the screening of the negligible integrals. Possible values are 'low', 'medium', 'high', 'very high', 'insane'. It could be abbreviated in 'l', 'm', 'h', 'vh', 'i'. 'high' is usually fine. 'insane' is only meant for debugging since it is overdoing a lot. 


---
### kerker_k0

**experimental** 

*Optional* 

**Family:** scf 

**Type:** real 

**Default:** 0.0 

**Description:** 

Analog to k0 in Kerker preconditioning for metallic systems. Helps to damp charge oscillations to ensure better SCF convergence. 


---
### length_unit

*Optional* 

**Family:** system 

**Type:** characters 

**Default:** angstrom 

**Description:** 

Chooses the units of the atomic coordinates. Can be 'angstrom' or 'bohr'. Could be abbreviated in 'A' or 'au'. 


---
### level_shifting_energy

*Optional* 

**Family:** scf 

**Type:** real 

**Default:** 0.0 

**Description:** 

Sets the energy shift up of the unoccupied states. Should help the convergence in the case of small HOMO-LUMO gaps. 


---
### magnetization

*Optional* 

**Family:** system 

**Type:** real 

**Default:** 0.0 

**Description:** 

Sets the number of unpaired electrons. In other words, this is the difference between the spin up and spin down occupation. For instance, a spin-doublet calculation is obtained with magnetization=1.0. Only meaningful when nspin=2. 


---
### memory_evaluation

*Optional* 

**Family:** general 

**Type:** yes/no 

**Default:** no 

**Description:** 

Requests a memory evaluation. MOLGW will start normaly, evaluate the size of the arrays, and exit without performing an actual calculation. 


---
### min_overlap

*Optional* 

**Family:** scf 

**Type:** real 

**Default:** 1e-05 

**Description:** 

Sets the minimal eigenvalue of the overlap matrix S. Small eigenvalues imply overcompleteness of the basis set. 


---
### mixing_scheme

*Optional* 

**Family:** scf 

**Type:** characters 

**Default:** pulay 

**Description:** 

Sets the density-matrix update method for SCF cycles. Possible choices are 'pulay' for Pulay DIIS method, 'adiis' for Hu-Yang method, or 'simple' for a simple linear mixing between input and output density-matrices. 


---
### move_nuclei

*Optional* 

**Family:** general 

**Type:** characters 

**Default:** no 

**Description:** 

Tells the code to move or not the position of the nuclei. Available options are 'no' or 'relax'. 


---
### mpi_nproc_ortho

*Optional* 

**Family:** hardware 

**Type:** integer 

**Default:** 1 

**Description:** 

Sets the number of processors left to parallelize on other directions. The main direction (auxiliary basis or DFT grid points) is obtained by <b>mpi_nproc</b> / <b>mpi_nproc_ortho</b>, which must be an integer. 


---
### n_hist

*Optional* 

**Family:** rt_tddft 

**Type:** integer 

**Default:** 2 

**Description:** 

Number of memorised previous hamiltonian values for its extrapolation in the real-time dynamics. n_hist=1 means that H(t_i+1)=H(t_i); n_hist=2 : H(t_i+1)=a*H(t_i)+b*(t_i-1); etc. 


---
### n_iter

*Optional* 

**Family:** rt_tddft 

**Type:** integer 

**Default:** 2 

**Description:** 

Sets the number of iterations for the PC7 in the real-time dynamics 


---
### n_restart_tddft

*Optional* 

**Family:** rt_tddft 

**Type:** integer 

**Default:** 50 

**Description:** 

RESTART_TDDFT file will be written during simulation each n_retart_tddft iteration (provided that print_tddft_restart is yes) 


---
### natom

*Optional* 

**Family:** system 

**Type:** integer 

**Default:** 0 

**Description:** 

Sets the number of atoms in the molecule. This is the number of lines to be read in the following section of the input file if no xyz file is provided. 


---
### ncore_tddft

*Optional* 

**Family:** rt_tddft 

**Type:** integer 

**Default:** 0 

**Description:** 

Sets the number of frozen core states in the real-time dynamics. 


---
### ncoreg

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 0 

**Description:** 

Sets the number of frozen core states in the Green's function G. 


---
### ncorew

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 0 

**Description:** 

Sets the number of frozen core states in the screened Coulomb interaction W, in TD-DFT, and in BSE. 


---
### nexcitation

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 0 

**Description:** 

Sets the number of neutral excitations to be calculated in TD-DFT or BSE.                  0 stands for all the states and triggers the full diagonalization. 


---
### nghost

*Optional* 

**Family:** system 

**Type:** integer 

**Default:** 0 

**Description:** 

Sets the number of ghost atoms in the molecule. Used to place basis function where there is no atom. Useful for Basis Set Superposition Error 


---
### noft_Lpower

*Optional* 

**Family:** scf 

**Type:** real 

**Default:** 0.53 

**Description:** 

Power functional approximation exponent used in NOFT calcs. 


---
### noft_NOTupdateOCC

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** no 

**Description:** 

Do a NOFT optimization but keeping fixed the occ numbers (or GAMMAS) read. 


---
### noft_NOTupdateORB

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** no 

**Description:** 

Do a NOFT optimization but keeping fixed the orbitals read. 


---
### noft_complex

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** no 

**Description:** 

Use complex molecular orb. coeficients in NOFT calcs. (default=no). 


---
### noft_fcidump

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** no 

**Description:** 

Print the FCIDUMP file in NOFT module. 


---
### noft_functional

*Optional* 

**Family:** scf 

**Type:** characters 

**Default:** PNOF7 

**Description:** 

Select the NOFT approx. to use (default=PNOF7). Other options are PNOF5, HF, MULLER, and POWER. 


---
### noft_ithresh_lambda

*Optional* 

**Family:** scf 

**Type:** integer 

**Default:** 5 

**Description:** 

Threshold used to determine [Lambda_pq - Lambda_qp*] hermiticity. 


---
### noft_lowmemERI

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** yes 

**Description:** 

Store the nat. orb. ERI as (all,occ,occ,occ) (default) or as (all,all,all,all) in NOFT module. 


---
### noft_ncoupled

*Optional* 

**Family:** scf 

**Type:** integer 

**Default:** 2 

**Description:** 

Number of coupled orbs. per pair used in NOFT calcs. (default=2 perfect pairing). 


---
### noft_ndiis

*Optional* 

**Family:** scf 

**Type:** integer 

**Default:** 5 

**Description:** 

Number of orb. optimization iterations used in DIIS by NOFT module. 


---
### noft_npairs

*Optional* 

**Family:** scf 

**Type:** integer 

**Default:** 1 

**Description:** 

Number of active electron pairs used in NOFT calcs. (default=1 pair). 


---
### noft_nscf

*Optional* 

**Family:** scf 

**Type:** integer 

**Default:** 1000 

**Description:** 

Maximum number of global iterations used by NOFT module. 


---
### noft_printdmn

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** no 

**Description:** 

Print optimized NOFT 1,2-RDMs (default=0 not to print them). 


---
### noft_printints

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** no 

**Description:** 

Print hCORE and ERImol integrals in the optimized basis (default=0 not to print them). 


---
### noft_printswdmn

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** no 

**Description:** 

Print optimized spin-with NOFT 1,2-RDMs (default=0 not to print them). 


---
### noft_readCOEF

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** no 

**Description:** 

Read NO_COEF file to use those coefficients as initial guess (default=0 not to read them). 


---
### noft_readFdiag

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** no 

**Description:** 

Read F_pp values from F_DIAG file and use them as indep. variables in occ. optimization (default=0 not to read them). 


---
### noft_readGAMMAS

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** no 

**Description:** 

Read Gammas_i from GAMMAS file and use them as indep. variables in occ. optimization (default=0 not to read them). 


---
### noft_readOCC

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** no 

**Description:** 

Read occ. from DM1 file and use them to compute Gammas_i (the indep. variables in occ. optimization). The default=0 not to read them. 


---
### noft_restart

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** no 

**Description:** 

Use binary files to restart NOFT calcs. (default=no). 


---
### noft_sta

*Optional* 

**Family:** scf 

**Type:** yes/no 

**Default:** yes 

**Description:** 

Decide whether to use PNOF7 or PNOF7s (default=1 for PNOF7s). 


---
### noft_tolE

*Optional* 

**Family:** scf 

**Type:** real 

**Default:** 1e-09 

**Description:** 

Threshold used to determine that the energy convergence change in NOFT calcs. is small, hence we have converged. 


---
### nomega_chi_imag

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 0 

**Description:** 

Sets the number of frequencies for the response function used to perform the integral on the imaginary axis 


---
### nomega_sigma

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 51 

**Description:** 

Sets the number of frequencies used to solve the quasiparticle equation in the GW self-energy. 


---
### nomega_sigma_calc

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 1 

**Description:** 

Sets the number of frequencies where the GW self-energy is actually calculated. 


---
### npulay_hist

*Optional* 

**Family:** scf 

**Type:** integer 

**Default:** 6 

**Description:** 

Sets the history record length for Pulay DIIS. 


---
### nscf

*Optional* 

**Family:** scf 

**Type:** integer 

**Default:** 50 

**Description:** 

Sets the maximum number of SCF cycles 


---
### nspin

*Optional* 

**Family:** system 

**Type:** integer 

**Default:** 1 

**Description:** 

Sets the number of spin channels. 1 enforces spin-restricted calculations. 2 means spin-unrestricted. 


---
### nstep

*Optional* 

**Family:** general 

**Type:** integer 

**Default:** 50 

**Description:** 

Sets the number of steps when moving the nuclei. 


---
### nstep_dav

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 15 

**Description:** 

Sets the maximum number of Davidson partial diagonalization steps. Used for TD-DFT, BSE, and full CI. 


---
### nstep_gw

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 1 

**Description:** 

Sets the number of GW iterations for eigenvalue self-consistent GW calculations (GnWn or GnW0). 


---
### nvel_projectile

*Optional* 

**Family:** postscf 

**Type:** integer 

**Default:** 1 

**Description:** 

Number of velocities used in linear-response stopping power. The first velocity is given by vel_projectile. The next ones are multiples of this initial value. 


---
### nvirtualg

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 100000 

**Description:** 

Sets the starting state beyond which states are excluded from the sum in the Green's function G. 


---
### nvirtualw

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 100000 

**Description:** 

Sets the starting state beyond which states are excluded from the sum in the screened Coulomb interaction W, in TD-DFT, and in BSE. 


---
### partition_scheme

*Optional* 

**Family:** scf 

**Type:** characters 

**Default:** ssf 

**Description:** 

Sets the partition scheme for the xc quadrature. Possible choices are 'becke' or 'ssf' (Stratmann-Scuseria-Frisch). 


---
### postscf

*Optional* 

**Family:** post 

**Type:** characters 

**Default:** None 

**Description:** 

Contains the post-processing scheme name. 
TD stands for TD-DFT or TD-HF.
BSE stands for Bethe-Salpeter.
GW stands for perturbative G0W0.
GnW0 stands for GW with eigenvalue self-consistentcy on G.
GnWn stands for GW with eigenvalue self-consistentcy on both G and W.
MP2 stands for guess what.
GWGAMMA (EXPERIMENTAL) stands for vertex corrections. 


---
### postscf_diago_flavor

*Optional* 

**Family:** post 

**Type:** characters 

**Default:**   

**Description:** 

Selects the LAPACK/ScaLAPACK diagonalization routines in the post SCF calculations. Available choices are ' ', 'R', 'D', and 'X'. 


---
### pred_corr

*Optional* 

**Family:** rt_tddft 

**Type:** characters 

**Default:** PC1 

**Description:** 

Sets the predictor-corrector scheme in the real-time dynamics. 


---
### print_bigrestart

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** yes 

**Description:** 

Prints the big RESTART file at the end of the SCF loop. There are two kinds of RESTART files: the small RESTART and the big RESTART. The latter is written only when self-consistency has been reached. It contains all the states and the Hamiltonian and allows one to completely skip the scf loop or to start over with another basis set. 


---
### print_cube

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints some wavefunctions in a 3D volumetric file with cube format 


---
### print_cube_diff_tddft

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints the difference of electronic density with respect to initial density in a 3D volumetric file with cube format for each simulation step in the real-time dynamics 


---
### print_cube_rho_tddft

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints electronic density in a 3D volumetric file with cube format for each simulation step in the real-time dynamics 


---
### print_dens_traj

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints the electronic density along the projectile trajectory for several impact parameters using real wave function 


---
### print_dens_traj_points_set

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints the electronic density between pairs of points given in manual_dens_points_set file. 


---
### print_dens_traj_tddft

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints the electronic density along the projectile trajectory for several impact parameters in the real-time dynamics 


---
### print_density_matrix

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints the density matrix in the DENSITY_MATRIX file 


---
### print_eri

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Dumps the Electron Repulsion Integral on a file. 


---
### print_hartree

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints the Hartree potential and exchange expectation value on eigenstates. 


---
### print_line_rho_diff_tddft

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints electronic density difference along a line, which parameters must be provided in manual_plot_rho_tddft file. 


---
### print_line_rho_tddft

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints electronic density along a line, which parameters must be provided in manual_plot_rho_tddft file. 


---
### print_multipole

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints the electric multipole expansion for the electronic density and the nuclei. 


---
### print_pdos

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints the Mulliken weight of each eigenvector on a given atom or a given series of atoms. 


---
### print_restart

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** yes 

**Description:** 

Prints a small RESTART file at each SCF cycle. There are two kinds of RESTART files: the small RESTART and the big RESTART. The former contains only the information about the occupied wavefunctions. This is a very small file and the writing should not hit too much on performance. 


---
### print_rho_grid

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints the electronic density discretized on the DFT grid into a file 'rho_grid.dat'. The density is calculated from the DENSITY_MATRIX file or from a Gaussian file using 'read_fchk'. 


---
### print_sigma

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints the value of the GW self-energy on the sampling frequencies in files. 


---
### print_spatial_extension

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints the wavefunction extension calculated as &lt;r<sup>2</sup>&gt; - &lt;r&gt;<sup>2</sup> 


---
### print_tddft_matrices

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints some matrices of the real-time dynamics into the file check_matrix.dat. 


---
### print_tddft_restart

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** yes 

**Description:** 

Prints a RESTART_TDDFT file which contains wavefunction coefficients for the last time moment of a simulation. 


---
### print_w

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Dumps the spectral function of the screened Coulomb W. This is necessary for a subsequent BSE run. 


---
### print_wfn

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints some wavefunctions along some selected lines. 


---
### print_wfn_files

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Prints WFN files which can be used as input for post processing. WFN files contain the information about the electronic density and several packages use them to gather chemical information (see programs like: AIMPAC, AIMALL, among others). Setting this variable to yes will produce WFN files for ground-state calculations as well as for GW-corrected densities. Note: currently, only cartesian basis are supported. 


---
### print_yaml

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** yes 

**Description:** 

Creates an output file in YAML format. Easier to read for python post-processing. 


---
### projectile_charge_scaling

*Optional* 

**Family:** rt_tddft 

**Type:** real 

**Default:** 1.0 

**Description:** 

Rescaling of the projectile charge 


---
### prop_type

*Optional* 

**Family:** rt_tddft 

**Type:** characters 

**Default:** MAG2 

**Description:** 

Sets the type of propagation algorithm in the real-time dynamics. 'CN' stands for Crank-Nickolson. 'MAG2' stands for Magnus 2nd order. 


---
### pt3_a_diagrams

*Optional* 

**Family:** post 

**Type:** characters 

**Default:** yes 

**Description:** 

Switch whether to calculate the A diagrams family in PT3. A diagrams are the self-consistent diagrams (PT2 inclusions in the Green's function). Valid choices include: 'yes', 'no', or 'only'. 


---
### pt_density_matrix

*Optional* 

**Family:** post 

**Type:** characters 

**Default:** no 

**Description:** 

Triggers the calculation of a correlated density matrix within MBPT. Valid choices include: 'no', 'PT2', 'ONE-RING', or 'GW'. 


---
### r_disc

*Optional* 

**Family:** rt_tddft 

**Type:** real 

**Default:** 200.0 

**Description:** 

Radius of the disc for density calculations (option calc_dens_disc) for the real-time dynamics. 


---
### rcut_mbpt

**experimental** 

*Optional* 

**Family:** post 

**Type:** real 

**Default:** 1.0 

**Description:** 

EXPERIMENTAL 


---
### read_fchk

*Optional* 

**Family:** io 

**Type:** characters 

**Default:** no 

**Description:** 

Triggers the reading of an external Gaussian formatted checkpoint file (named gaussian.fchk) that contains density matrices. Basis sets have to be precisely the same in MOLGW and in Gaussian, which requires a manual input of the basis set in both codes. Options are 'no' (no reading), 'SCF' (for self-consistent field), 'CC' (for coupled-cluster), or 'MP2' (for MP2). Today, only works for Cartesian Gaussian and for spin restricted calculations. 


---
### read_restart

*Optional* 

**Family:** io 

**Type:** yes/no 

**Default:** no 

**Description:** 

Read the RESTART file and restart from it. 


---
### read_tddft_restart

*Optional* 

**Family:** io_rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Ignore the RESTART_TDDFT file. 


---
### scalapack_block_min

*Optional* 

**Family:** hardware 

**Type:** integer 

**Default:** 100000 

**Description:** 

Sets the minimum block size to distribute a non-distributed matrix with SCALAPACK. If scalapack_block_min=400, then a 900x900 matrix will be distributed on a 2x2 processor grid. If scalapack_block_min=500, then a 900x900 matrix will no be distributed. 


---
### scf

*Mandatory* 

**Family:** general 

**Type:** characters 

**Default:** None 

**Description:** 

Contains the self-consistent scheme name. 
Try LDA, PBE, HSE06, or HF for instance 


---
### scf_diago_flavor

*Optional* 

**Family:** scf 

**Type:** characters 

**Default:**   

**Description:** 

Selects the LAPACK/ScaLAPACK diagonalization routines in the SCF cycles. Available choices are ' ', 'R', 'D', and 'X'. 


---
### scissor

*Optional* 

**Family:** post 

**Type:** real 

**Default:** 0.0 

**Description:** 

Sets a rigid energy shift of the unoccupied states, so to mimick a GW calculation without actually doing it. 


---
### selfenergy_state_max

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 100000 

**Description:** 

Sets the final states for the range of the self-energy evaluation 


---
### selfenergy_state_min

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 1 

**Description:** 

Sets the starting states for the range of the self-energy evaluation 


---
### selfenergy_state_range

*Optional* 

**Family:** post 

**Type:** integer 

**Default:** 100000 

**Description:** 

Sets the range of states around the HOMO level for the self-energy evaluation. For instance, selfenergy_state_range=0 will trigger the calculation of the HOMO only. selfenergy_state_range=1 will trigger the evaluation of the HOMO-1, HOMO, HOMO+1. etc. 


---
### small_basis

**experimental** 

*Optional* 

**Family:** post 

**Type:** characters 

**Default:** None 

**Description:** 

Calls for a smaller basis set used to represent the virtual orbital space with fewer functions. Only meaningful for GW. 


---
### step_sigma

*Optional* 

**Family:** post 

**Type:** real 

**Default:** 0.01 

**Description:** 

Sets the spacing between frequencies in the final GW self-energy output. 


---
### step_sigma_calc

*Optional* 

**Family:** post 

**Type:** real 

**Default:** 0.03 

**Description:** 

Sets the spacing between the frequencies where the GW self-energy is actually calculated. 


---
### stopping

*Optional* 

**Family:** post 

**Type:** characters 

**Default:** no 

**Description:** 

Triggers the calculation of the stopping power within linear-response theory. Only effective when postscf=''td'' or ''bse''. Avaialble values are ''no'', ''spherical'', ''3d''. 


---
### tda

*Optional* 

**Family:** post 

**Type:** yes/no 

**Default:** no 

**Description:** 

Triggers the use of Tamm-Dancoff approximation in TD-DFT or BSE. 


---
### tddft_frozencore

*Optional* 

**Family:** rt_tddft 

**Type:** yes/no 

**Default:** no 

**Description:** 

Do not "propagate" states mentioned in the manual_tddft_frozencore file in the real-time dynamics. 


---
### tddft_grid_quality

*Optional* 

**Family:** post 

**Type:** characters 

**Default:** high 

**Description:** 

Sets the number of grid points use to evaluate the exchange-correlation integrals in real space for the TDDFT kernel. Possible values are 'low', 'medium', 'high', 'very high', 'insane'. It could be abbreviated in 'l', 'm', 'h', 'vh', 'i'. 'high' is usually fine. 'insane' is only meant for debugging since it is overdoing a lot. 


---
### temperature

*Optional* 

**Family:** system 

**Type:** real 

**Default:** 0.0 

**Description:** 

Sets the electronic temperature in the Fermi-Dirac functions. Helps the convergence for some systems. The value is input in Hartree atomic units. 


---
### time_sim

*Optional* 

**Family:** rt_tddft 

**Type:** real 

**Default:** 10.0 

**Description:** 

Duration of a real-time dynamics in atomic units. 


---
### time_step

*Optional* 

**Family:** rt_tddft 

**Type:** real 

**Default:** 1.0 

**Description:** 

Time step for real-time dynamics in atomic units. 


---
### toldav

*Optional* 

**Family:** post 

**Type:** real 

**Default:** 0.0001 

**Description:** 

Sets the tolerance criterium for the maximum norm of the residual in the Davidson diagonalization of TD-DFT, BSE, and full CI. 


---
### tolforce

*Optional* 

**Family:** general 

**Type:** real 

**Default:** 1e-05 

**Description:** 

Sets the target threshold for the maximum force component after nuclei relaxation. 


---
### tolscf

*Optional* 

**Family:** scf 

**Type:** real 

**Default:** 1e-07 

**Description:** 

Sets the residual norm target for the density matrix for the SCF cycles. 


---
### triplet

*Optional* 

**Family:** post 

**Type:** yes/no 

**Default:** no 

**Description:** 

Triggers the calculation of the triplet final state in TD-DFT or BSE. 


---
### use_correlated_density_matrix

*Optional* 

**Family:** post 

**Type:** yes/no 

**Default:** no 

**Description:** 

Chooses to use another density matrix for the Fock hamiltonian to be employed in self-energy calculations.                  Used in conjonction with 'pt_density_matrix' or with 'read_fchk' or read an existing DENSITY_MATRIX file. 


---
### vel_projectile

*Optional* 

**Family:** rt_tddft 

**Type:** vector_1d_3 

**Default:** (0.0, 0.0, 1.0) 

**Description:** 

Projectile initial velocity. Used for real-time tddft and for linear-response stopping power calculations 


---
### virtual_fno

*Optional* 

**Family:** post 

**Type:** yes/no 

**Default:** no 

**Description:** 

Activates the Frozen Natural Orbitals technique to span the virtual orbitals subspace with fewer orbitals. The dimension of the space is set up with the input keyword nvirtualg or nvirtualw. Actually the virtual orbital space is determined by the minimum MIN(nvirtualg,nvirtualw). 


---
### write_step

*Optional* 

**Family:** io_rt_tddft 

**Type:** real 

**Default:** 1 

**Description:** 

Determines the time step for data recording in the real-time dynamics 


---
### xyz_file

*Optional* 

**Family:** system 

**Type:** characters 

**Default:** None 

**Description:** 

Specifies the location of the xyz file that contains the atomic positions. It can be used as an alternate route to set atomic coordinate. 




*Generated by input_variables.py on 15 March 2022* 


