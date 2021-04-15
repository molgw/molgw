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
[density_matrix_damping](#density_matrix_damping) 
[diis_switch](#diis_switch) 
[gamma_hybrid](#gamma_hybrid) 
[grid_quality](#grid_quality) 
[init_hamiltonian](#init_hamiltonian) 
[integral_quality](#integral_quality) 
[kerker_k0](#kerker_k0) 
[level_shifting_energy](#level_shifting_energy) 
[min_overlap](#min_overlap) 
[mixing_scheme](#mixing_scheme) 
[npulay_hist](#npulay_hist) 
[nscf](#nscf) 
[partition_scheme](#partition_scheme) 
[scf_diago_flavor](#scf_diago_flavor) 
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
### auto_auxil_fsam

*Optional* 

Type: real 

Default: 1.5 

Sets the F_SAM parameter in the automatic generation of the auxiliary basis set. The closer to 1.0 the more auxiliary basis functions it will generate. See Yang-Rendell-Frisch, JChemPhys 127, 074102 (2007) for more details. 


---
### auto_auxil_lmaxinc

*Optional* 

Type: integer 

Default: 1 

Sets the l_MAXINC parameter in the automatic generation of the auxiliary basis set. The larger the more auxiliary basis functions it will generate. See Yang-Rendell-Frisch, JChemPhys 127, 074102 (2007) for more details. 


---
### auxil_basis

*Optional* 

Type: characters 

Default: None 

Sets the auxiliary basis set. For instance, cc-pVDZ-RI for a Weigend basis set. If present, the auxiliary basis will be used for both the scf cycles and the postscf calculations (TD-DFT, BSE, or GW). 


---
### basis

*Mandatory* 

Type: characters 

Default: None 

Sets the basis set For Pople sets, use 6-31G for instance or 6-31+G\*. For Dunning sets, use aug-cc-pVTZ for instance. Note that Pople sets are to be used with gaussian_type='cart' One may use ones own basis sets provided that the files are labeled X_mybasisset where X is the element. 


---
### comment

*Optional* 

Type: characters 

Default: None 

This is a free expression place. Use it as you wish for commenting, naming, labeling etc. (140 character max just as twitter) 


---
### ecp_auxil_basis

*Optional* 

Type: characters 

Default: None 

Name of the auxiliary basis set to be used for elements specified in list ecp_elements. 


---
### ecp_basis

*Optional* 

Type: characters 

Default: None 

Name of the basis set to be used for elements specified in list ecp_elements. 


---
### ecp_elements

*Optional* 

Type: characters 

Default: None 

Contains the list of elements (separated by spaces) that should be treated with an Effective Core Potential. 


---
### ecp_quality

*Optional* 

Type: characters 

Default: high 

Sets the number of grid points use to evaluate the Effective Core Potential integrals in real space. Possible values are 'low', 'medium', 'high', 'very high', 'insane'. It could be abbreviated in 'l', 'm', 'h', 'vh', 'i'. 'high' is usually fine. 'insane' is only meant for debugging since it is overdoing a lot. 


---
### ecp_type

*Optional* 

Type: characters 

Default: None 

Name of the Effective Core Potential. For instance, Gold using the cc-pVDZ-PP basis set should have ecp_type='PP', so that MOLGW looks for the file Au_PP in the basis_path folder. 


---
### eri3_genuine

*Optional* 

Type: yes/no 

Default: no 

If set to 'yes', the 2-center integrals will not be diagonalized and 3-center electron repulsion integrals will remain the genuine 3-center integrals. 


---
### gaussian_type

*Optional* 

Type: characters 

Default: pure 

Asks for pure or spherical Gaussian type orbitals with 'pure' or for Cartesian Gaussian orbital with 'cart'. 


---
### incore

*Optional* 

Type: yes/no 

Default: yes 

Specify if the 4-center integrals are all calculated at once and stored or if they are calculated on-the-fly. 


---
### memory_evaluation

*Optional* 

Type: yes/no 

Default: no 

Requests a memory evaluation. MOLGW will start normaly, evaluate the size of the arrays, and exit without performing an actual calculation. 


---
### move_nuclei

*Optional* 

Type: characters 

Default: no 

Tells the code to move or not the position of the nuclei. Available options are 'no' or 'relax'. 


---
### nstep

*Optional* 

Type: integer 

Default: 50 

Sets the number of steps when moving the nuclei. 


---
### scf

*Mandatory* 

Type: characters 

Default: None 

Contains the self-consistent scheme name. 
Try LDA, PBE, HSE06, or HF for instance 


---
### tolforce

*Optional* 

Type: real 

Default: 1e-05 

Sets the target threshold for the maximum force component after nuclei relaxation. 


---
### eri3_nbatch

*Optional* 

Type: integer 

Default: 1 

Sets the number of batches when calculating the 3-center integrals. Having a large eri3_nbatch reduces the memory foot print, however it may lower the performance. 


---
### eri3_npcol

*Optional* 

Type: integer 

Default: 1 

Sets number of column processors for the distribution of the 3-center integrals.  eri3_nprow X eri3_npcol must be equal to the number of MPI threads else MOLGW decides on its own. 


---
### eri3_nprow

*Optional* 

Type: integer 

Default: 1 

Sets number of row processors for the distribution of the 3-center integrals.  eri3_nprow X eri3_npcol must be equal to the number of MPI threads else MOLGW decides on its own. 


---
### grid_memory

*Optional* 

Type: real 

Default: 400.0 

Sets the maximum memory usage in Mb allowed to store the wavefunctions on the quadrature points for XC integrals. 


---
### mpi_nproc_ortho

*Optional* 

Type: integer 

Default: 1 

Sets the number of processors left to parallelize on other directions. The main direction (auxiliary basis or DFT grid points) is obtained by <b>mpi_nproc</b> / <b>mpi_nproc_ortho</b>, which must be an integer. 


---
### scalapack_block_min

*Optional* 

Type: integer 

Default: 100000 

Sets the minimum block size to distribute a non-distributed matrix with SCALAPACK. If scalapack_block_min=400, then a 900x900 matrix will be distributed on a 2x2 processor grid. If scalapack_block_min=500, then a 900x900 matrix will no be distributed. 


---
### basis_path

*Optional* 

Type: characters 

Default: None 

Sets the path pointing to the basis functions files.                  If not specified, then the basis set files will be searched in folder ~molgw/basis/. 


---
### force_energy_qp

*Optional* 

Type: yes/no 

Default: no 

Force the reading of the ENERGY_QP file whatever the postscf choice. 


---
### ignore_bigrestart

*Optional* 

Type: yes/no 

Default: no 

Considers a big RESTART as if it was a small RESTART. 


---
### print_bigrestart

*Optional* 

Type: yes/no 

Default: yes 

Prints the big RESTART file at the end of the SCF loop. There are two kinds of RESTART files: the small RESTART and the big RESTART. The latter is written only when self-consistency has been reached. It contains all the states and the Hamiltonian and allows one to completely skip the scf loop or to start over with another basis set. 


---
### print_cube

*Optional* 

Type: yes/no 

Default: no 

Prints some wavefunctions in a 3D volumetric file with cube format 


---
### print_wfn_files

*Optional* 

Type: yes/no 

Default: no 

Prints WFN files which can be used as input for post processing. WFN files contain the information about the electronic density and several packages use them to gather chemical information (see programs like: AIMPAC, AIMALL, among others). Setting this variable to yes will produce WFN files for ground-state calculations as well as for GW-corrected densities. Note: currently, only cartesian basis are supported. 


---
### print_density_matrix

*Optional* 

Type: yes/no 

Default: no 

Prints the density matrix in the DENSITY_MATRIX file 


---
### print_eri

*Optional* 

Type: yes/no 

Default: no 

Dumps the Electron Repulsion Integral on a file. 


---
### print_hartree

*Optional* 

Type: yes/no 

Default: no 

Prints the Hartree potential and exchange expectation value on eigenstates. 


---
### print_multipole

*Optional* 

Type: yes/no 

Default: no 

Prints the electric multipole expansion for the electronic density and the nuclei. 


---
### print_pdos

*Optional* 

Type: yes/no 

Default: no 

Prints the Mulliken weight of each eigenvector on a given atom or a given series of atoms. 


---
### print_restart

*Optional* 

Type: yes/no 

Default: yes 

Prints a small RESTART file at each SCF cycle. There are two kinds of RESTART files: the small RESTART and the big RESTART. The former contains only the information about the occupied wavefunctions. This is a very small file and the writing should not hit too much on performance. 


---
### print_rho_grid

*Optional* 

Type: yes/no 

Default: no 

Prints the electronic density discretized on the DFT grid into a file 'rho_grid.dat'. The density is calculated from the DENSITY_MATRIX file or from a Gaussian file using 'read_fchk'. 


---
### print_sigma

*Optional* 

Type: yes/no 

Default: no 

Prints the value of the GW self-energy on the sampling frequencies in files. 


---
### print_spatial_extension

*Optional* 

Type: yes/no 

Default: no 

Prints the wavefunction extension calculated as &lt;r<sup>2</sup>&gt; - &lt;r&gt;<sup>2</sup> 


---
### print_w

*Optional* 

Type: yes/no 

Default: no 

Dumps the spectral function of the screened Coulomb W. This is necessary for a subsequent BSE run. 


---
### print_wfn

*Optional* 

Type: yes/no 

Default: no 

Prints some wavefunctions along some selected lines. 


---
### print_yaml

*Optional* 

Type: yes/no 

Default: yes 

Creates an output file in YAML format. Easier to read for python post-processing. 


---
### read_fchk

*Optional* 

Type: characters 

Default: no 

Triggers the reading of an external Gaussian formatted checkpoint file (named gaussian.fchk) that contains density matrices. Basis sets have to be precisely the same in MOLGW and in Gaussian, which requires a manual input of the basis set in both codes. Options are 'no' (no reading), 'SCF' (for self-consistent field), 'CC' (for coupled-cluster), or 'MP2' (for MP2). Today, only works for Cartesian Gaussian and for spin restricted calculations. 


---
### read_restart

*Optional* 

Type: yes/no 

Default: no 

Read the RESTART file and restart from it. 


---
### calc_dens_disc

*Optional* 

Type: yes/no 

Default: no 

Calculate electronic density in the discs during the real-time dynamics 


---
### calc_q_matrix

*Optional* 

Type: yes/no 

Default: no 

Calculate and print q_matrix which is the projection of a propagated state psi(t) onto the initial state psi(0) in the real-time dynamics 


---
### calc_spectrum

*Optional* 

Type: yes/no 

Default: no 

Calculates absorption spectrum in the real-time dynamics 


---
### print_cube_diff_tddft

*Optional* 

Type: yes/no 

Default: no 

Prints the difference of electronic density with respect to initial density in a 3D volumetric file with cube format for each simulation step in the real-time dynamics 


---
### print_cube_rho_tddft

*Optional* 

Type: yes/no 

Default: no 

Prints electronic density in a 3D volumetric file with cube format for each simulation step in the real-time dynamics 


---
### print_dens_traj

*Optional* 

Type: yes/no 

Default: no 

Prints the electronic density along the projectile trajectory for several impact parameters using real wave function 


---
### print_dens_traj_points_set

*Optional* 

Type: yes/no 

Default: no 

Prints the electronic density between pairs of points given in manual_dens_points_set file. 


---
### print_dens_traj_tddft

*Optional* 

Type: yes/no 

Default: no 

Prints the electronic density along the projectile trajectory for several impact parameters in the real-time dynamics 


---
### print_line_rho_diff_tddft

*Optional* 

Type: yes/no 

Default: no 

Prints electronic density difference along a line, which parameters must be provided in manual_plot_rho_tddft file. 


---
### print_line_rho_tddft

*Optional* 

Type: yes/no 

Default: no 

Prints electronic density along a line, which parameters must be provided in manual_plot_rho_tddft file. 


---
### print_tddft_matrices

*Optional* 

Type: yes/no 

Default: no 

Prints some matrices of the real-time dynamics into the file check_matrix.dat. 


---
### print_tddft_restart

*Optional* 

Type: yes/no 

Default: yes 

Prints a RESTART_TDDFT file which contains wavefunction coefficients for the last time moment of a simulation. 


---
### read_tddft_restart

*Optional* 

Type: yes/no 

Default: no 

Ignore the RESTART_TDDFT file. 


---
### write_step

*Optional* 

Type: real 

Default: 1 

Determines the time step for data recording in the real-time dynamics 


---
### assume_scf_converged

*Optional* 

Type: yes/no 

Default: no 

Allows for a post-scf calculation whatever the outcome of the SCF loop. Especially useful when restarting from another SCF calculations with nscf=0 


---
### ci_greens_function

**experimental** 

*Optional* 

Type: characters 

Default: holes 

EXPERIMENTAL. Selects which part of the Green's function is to be calculated: holes, electrons, or both. 


---
### ci_nstate

*Optional* 

Type: integer 

Default: 1 

Selects how many CI states should be calculated in the diagonalization. If ci_nstate is lower than the number of configuration,  a Davidson partial diagonalization is performed, else a full (SCA)LAPACK diagonalization is triggered. 


---
### ci_nstate_self

*Optional* 

Type: integer 

Default: 1 

Selects how many CI states in the N+1 or N-1 electron calculations. If ci_nstate_self is lower than the number of configuration,  a Davidson partial diagonalization is performed, else a full (SCA)LAPACK diagonalization is triggered. 


---
### ci_spin_multiplicity

*Optional* 

Type: integer 

Default: 1 

Spin multiplicity in CI calculations. 


---
### ci_type

*Optional* 

Type: characters 

Default: all 

Selects which excitations will be included in the CI expansion. Valid choices are 'all', 'CISD', 'CISDT', 'CISDTQ'. 


---
### dft_core

*Optional* 

Type: integer 

Default: 0 

Sets the number of states considered as core in &lt;&Sigma;<sub>x</sub>-<i>v</i><sub>xc</sub>&gt. This options is meant to mimic the pseudopotential approximation. 


---
### ecp_small_basis

**experimental** 

*Optional* 

Type: characters 

Default: None 

Calls for a smaller basis set used to represent the virtual orbital space with fewer functions. This is the small basis set used for elements with an effective core potential. Only meaningful for GW. 


---
### eta

*Optional* 

Type: real 

Default: 0.001 

Is a the tiny imaginary part used in the denominator of the Green's function to shift the pole off the axis, so to avoid divergences.This is an energy in Hartree. It should be set to the lowest value possible in theory. However, in practice, a too low value of eta would induce huge and unstable GW corrections. The default value is usually very accurate and there is no need to use a lower value. But for states apart from the band gap, a large value of eta may be beneficial for stability. eta=0.01 is already much more stable. Note that for QSGW increasing eta is most often unavoidable. 


---
### frozencore

*Optional* 

Type: yes/no 

Default: no 

Triggers the neglect of core states in GW. H, He, Li, Be have no core states. B-Na have the 1s. Al-Ca have the 1s2s2p. Manual tuning could be achieved with ncoreg, ncorew. 


---
### gwgamma_tddft

**experimental** 

*Optional* 

Type: yes/no 

Default: no 

EXPERIMENTAL. Calculates the vertex using the DFT flavor specified in the ground-state calculation. 


---
### ncoreg

*Optional* 

Type: integer 

Default: 0 

Sets the number of frozen core states in the Green's function G. 


---
### ncorew

*Optional* 

Type: integer 

Default: 0 

Sets the number of frozen core states in the screened Coulomb interaction W, in TD-DFT, and in BSE. 


---
### nexcitation

*Optional* 

Type: integer 

Default: 0 

Sets the number of neutral excitations to be calculated in TD-DFT or BSE.                  0 stands for all the states and triggers the full diagonalization. 


---
### nomega_chi_imag

*Optional* 

Type: integer 

Default: 0 

Sets the number of frequencies for the response function used to perform the integral on the imaginary axis 


---
### nomega_sigma

*Optional* 

Type: integer 

Default: 51 

Sets the number of frequencies used to solve the quasiparticle equation in the GW self-energy. 


---
### nomega_sigma_calc

*Optional* 

Type: integer 

Default: 1 

Sets the number of frequencies where the GW self-energy is actually calculated. 


---
### nstep_dav

*Optional* 

Type: integer 

Default: 15 

Sets the maximum number of Davidson partial diagonalization steps. Used for TD-DFT, BSE, and full CI. 


---
### nstep_gw

*Optional* 

Type: integer 

Default: 1 

Sets the number of GW iterations for eigenvalue self-consistent GW calculations (GnWn or GnW0). 


---
### nvel_projectile

*Optional* 

Type: integer 

Default: 1 

Number of velocities used in linear-response stopping power. The first velocity is given by vel_projectile. The next ones are multiples of this initial value. 


---
### nvirtualg

*Optional* 

Type: integer 

Default: 100000 

Sets the starting state beyond which states are excluded from the sum in the Green's function G. 


---
### nvirtualw

*Optional* 

Type: integer 

Default: 100000 

Sets the starting state beyond which states are excluded from the sum in the screened Coulomb interaction W, in TD-DFT, and in BSE. 


---
### postscf

*Optional* 

Type: characters 

Default: None 

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

Type: characters 

Default:   

Selects the LAPACK/ScaLAPACK diagonalization routines in the post SCF calculations. Available choices are ' ', 'R', 'D', and 'X'. 


---
### pt3_a_diagrams

*Optional* 

Type: characters 

Default: yes 

Switch whether to calculate the A diagrams family in PT3. A diagrams are the self-consistent diagrams (PT2 inclusions in the Green's function). Valid choices include: 'yes', 'no', or 'only'. 


---
### pt_density_matrix

*Optional* 

Type: characters 

Default: no 

Triggers the calculation of a correlated density matrix within MBPT. Valid choices include: 'no', 'PT2', 'ONE-RING', or 'GW'. 


---
### rcut_mbpt

**experimental** 

*Optional* 

Type: real 

Default: 1.0 

EXPERIMENTAL 


---
### scissor

*Optional* 

Type: real 

Default: 0.0 

Sets a rigid energy shift of the unoccupied states, so to mimick a GW calculation without actually doing it. 


---
### selfenergy_state_max

*Optional* 

Type: integer 

Default: 100000 

Sets the final states for the range of the self-energy evaluation 


---
### selfenergy_state_min

*Optional* 

Type: integer 

Default: 1 

Sets the starting states for the range of the self-energy evaluation 


---
### selfenergy_state_range

*Optional* 

Type: integer 

Default: 100000 

Sets the range of states around the HOMO level for the self-energy evaluation. For instance, selfenergy_state_range=0 will trigger the calculation of the HOMO only. selfenergy_state_range=1 will trigger the evaluation of the HOMO-1, HOMO, HOMO+1. etc. 


---
### small_basis

**experimental** 

*Optional* 

Type: characters 

Default: None 

Calls for a smaller basis set used to represent the virtual orbital space with fewer functions. Only meaningful for GW. 


---
### step_sigma

*Optional* 

Type: real 

Default: 0.01 

Sets the spacing between frequencies in the final GW self-energy output. 


---
### step_sigma_calc

*Optional* 

Type: real 

Default: 0.03 

Sets the spacing between the frequencies where the GW self-energy is actually calculated. 


---
### stopping

*Optional* 

Type: characters 

Default: no 

Triggers the calculation of the stopping power within linear-response theory. Only effective when postscf=''td'' or ''bse''. Avaialble values are ''no'', ''spherical'', ''3d''. 


---
### tda

*Optional* 

Type: yes/no 

Default: no 

Triggers the use of Tamm-Dancoff approximation in TD-DFT or BSE. 


---
### tddft_grid_quality

*Optional* 

Type: characters 

Default: high 

Sets the number of grid points use to evaluate the exchange-correlation integrals in real space for the TDDFT kernel. Possible values are 'low', 'medium', 'high', 'very high', 'insane'. It could be abbreviated in 'l', 'm', 'h', 'vh', 'i'. 'high' is usually fine. 'insane' is only meant for debugging since it is overdoing a lot. 


---
### toldav

*Optional* 

Type: real 

Default: 0.0001 

Sets the tolerance criterium for the maximum norm of the residual in the Davidson diagonalization of TD-DFT, BSE, and full CI. 


---
### triplet

*Optional* 

Type: yes/no 

Default: no 

Triggers the calculation of the triplet final state in TD-DFT or BSE. 


---
### use_correlated_density_matrix

*Optional* 

Type: yes/no 

Default: no 

Chooses to use another density matrix for the Fock hamiltonian to be employed in self-energy calculations.                  Used in conjonction with 'pt_density_matrix' or with 'read_fchk' or read an existing DENSITY_MATRIX file. 


---
### virtual_fno

*Optional* 

Type: yes/no 

Default: no 

Activates the Frozen Natural Orbitals technique to span the virtual orbitals subspace with fewer orbitals. The dimension of the space is set up with the input keyword nvirtualg or nvirtualw. Actually the virtual orbital space is determined by the minimum MIN(nvirtualg,nvirtualw). 


---
### excit_dir

*Optional* 

Type: vector_1d_3 

Default: (1.0, 0.0, 0.0) 

Excitation direction for the real-time dynamics. 


---
### excit_kappa

*Optional* 

Type: real 

Default: 2e-05 

Maximum Gaussian excitation field strength in atomic units. 


---
### excit_name

*Optional* 

Type: characters 

Default: no 

Sets the type of excitation of a system in the real-time dynamics.  'GAU' stands for a linearly polarized uniform Gaussian electric field 


---
### excit_omega

*Optional* 

Type: real 

Default: 0.2 

The excitation pulse width in atomic units for the real-time dynamics. 


---
### excit_time0

*Optional* 

Type: real 

Default: 3.0 

Center of the excitation pulse in atomic units for the real-time dynamics. 


---
### n_hist

*Optional* 

Type: integer 

Default: 2 

Number of memorised previous hamiltonian values for its extrapolation in the real-time dynamics. n_hist=1 means that H(t_i+1)=H(t_i); n_hist=2 : H(t_i+1)=a*H(t_i)+b*(t_i-1); etc. 


---
### n_iter

*Optional* 

Type: integer 

Default: 2 

Sets the number of iterations for the PC7 in the real-time dynamics 


---
### n_restart_tddft

*Optional* 

Type: integer 

Default: 50 

RESTART_TDDFT file will be written during simulation each n_retart_tddft iteration (provided that print_tddft_restart is yes) 


---
### ncore_tddft

*Optional* 

Type: integer 

Default: 0 

Sets the number of frozen core states in the real-time dynamics. 


---
### pred_corr

*Optional* 

Type: characters 

Default: PC1 

Sets the predictor-corrector scheme in the real-time dynamics. 


---
### prop_type

*Optional* 

Type: characters 

Default: MAG2 

Sets the type of propagation algorithm in the real-time dynamics. 'CN' stands for Crank-Nickolson. 'MAG2' stands for Magnus 2nd order. 


---
### projectile_charge_scaling

*Optional* 

Type: real 

Default: 1.0 

Rescaling of the projectile charge 


---
### r_disc

*Optional* 

Type: real 

Default: 200.0 

Radius of the disc for density calculations (option calc_dens_disc) for the real-time dynamics. 


---
### tddft_frozencore

*Optional* 

Type: yes/no 

Default: no 

Do not "propagate" states mentioned in the manual_tddft_frozencore file in the real-time dynamics. 


---
### time_sim

*Optional* 

Type: real 

Default: 10.0 

Duration of a real-time dynamics in atomic units. 


---
### time_step

*Optional* 

Type: real 

Default: 1.0 

Time step for real-time dynamics in atomic units. 


---
### vel_projectile

*Optional* 

Type: vector_1d_3 

Default: (0.0, 0.0, 1.0) 

Projectile initial velocity. Used for real-time tddft and for linear-response stopping power calculations 


---
### alpha_hybrid

*Optional* 

Type: real 

Default: 0.0 

Only works for Range-Separated hybrid functionals scf='rsh' Sets the amount of range-independent exact-exchange 


---
### alpha_mixing

*Optional* 

Type: real 

Default: 0.7 

Sets the amount of output density-matrix for the next iteration. When the SCF cycles have difficulties to converge, one may try to lower this value. 


---
### beta_hybrid

*Optional* 

Type: real 

Default: 0.0 

Only works for Range-Separated hybrid functionals scf='rsh' Sets the amount of long-range exact-exchange 


---
### density_matrix_damping

*Optional* 

Type: real 

Default: 0.0 

Adds an additional linear mixing on the density matrix in combination with the Hamiltonian mixing in order to damp out the charge oscillations. Especially useful for metallic systems. 


---
### diis_switch

*Optional* 

Type: real 

Default: 0.05 

When running ADIIS, sets the residue value below which the DIIS method is used to finalize the convergence. 


---
### gamma_hybrid

*Optional* 

Type: real 

Default: 1000000.0 

Only works for Range-Separated hybrid functionals scf='rsh' Sets the separation between long-range and short-range. It is input in bohr^-1. 


---
### grid_quality

*Optional* 

Type: characters 

Default: high 

Sets the number of grid points use to evaluate the exchange-correlation integrals in real space for the DFT potential and energy. Possible values are 'low', 'medium', 'high', 'very high', 'insane'. It could be abbreviated in 'l', 'm', 'h', 'vh', 'i'. 'high' is usually fine. 'insane' is only meant for debugging since it is overdoing a lot. 


---
### init_hamiltonian

*Optional* 

Type: characters 

Default: guess 

Selects how to initiate the first hamiltonian for SCF cycles. Today, two options are available: 'guess' for an educated guess based on approximate atomic densities or 'core' for the core hamiltonian. 


---
### integral_quality

*Optional* 

Type: characters 

Default: high 

Sets the tolerance value for the screening of the negligible integrals. Possible values are 'low', 'medium', 'high', 'very high', 'insane'. It could be abbreviated in 'l', 'm', 'h', 'vh', 'i'. 'high' is usually fine. 'insane' is only meant for debugging since it is overdoing a lot. 


---
### kerker_k0

**experimental** 

*Optional* 

Type: real 

Default: 0.0 

Analog to k0 in Kerker preconditioning for metallic systems. Helps to damp charge oscillations to ensure better SCF convergence. 


---
### level_shifting_energy

*Optional* 

Type: real 

Default: 0.0 

Sets the energy shift up of the unoccupied states. Should help the convergence in the case of small HOMO-LUMO gaps. 


---
### min_overlap

*Optional* 

Type: real 

Default: 1e-05 

Sets the minimal eigenvalue of the overlap matrix S. Small eigenvalues imply overcompleteness of the basis set. 


---
### mixing_scheme

*Optional* 

Type: characters 

Default: pulay 

Sets the density-matrix update method for SCF cycles. Possible choices are 'pulay' for Pulay DIIS method, 'adiis' for Hu-Yang method, or 'simple' for a simple linear mixing between input and output density-matrices. 


---
### npulay_hist

*Optional* 

Type: integer 

Default: 6 

Sets the history record length for Pulay DIIS. 


---
### nscf

*Optional* 

Type: integer 

Default: 50 

Sets the maximum number of SCF cycles 


---
### partition_scheme

*Optional* 

Type: characters 

Default: ssf 

Sets the partition scheme for the xc quadrature. Possible choices are 'becke' or 'ssf' (Stratmann-Scuseria-Frisch). 


---
### scf_diago_flavor

*Optional* 

Type: characters 

Default:   

Selects the LAPACK/ScaLAPACK diagonalization routines in the SCF cycles. Available choices are ' ', 'R', 'D', and 'X'. 


---
### tolscf

*Optional* 

Type: real 

Default: 1e-07 

Sets the residual norm target for the density matrix for the SCF cycles. 


---
### charge

*Optional* 

Type: real 

Default: 0.0 

Sets the total charge of the system. 0 is a neutral system. -2 is a doubly charged anion etc. 


---
### length_unit

*Optional* 

Type: characters 

Default: angstrom 

Chooses the units of the atomic coordinates. Can be 'angstrom' or 'bohr'. Could be abbreviated in 'A' or 'au'. 


---
### magnetization

*Optional* 

Type: real 

Default: 0.0 

Sets the number of unpaired electrons. In other words, this is the difference between the spin up and spin down occupation. For instance, a spin-doublet calculation is obtained with magnetization=1.0. Only meaningful when nspin=2. 


---
### natom

*Optional* 

Type: integer 

Default: 0 

Sets the number of atoms in the molecule. This is the number of lines to be read in the following section of the input file if no xyz file is provided. 


---
### nghost

*Optional* 

Type: integer 

Default: 0 

Sets the number of ghost atoms in the molecule. Used to place basis function where there is no atom. Useful for Basis Set Superposition Error 


---
### nspin

*Optional* 

Type: integer 

Default: 1 

Sets the number of spin channels. 1 enforces spin-restricted calculations. 2 means spin-unrestricted. 


---
### temperature

*Optional* 

Type: real 

Default: 0.0 

Sets the electronic temperature in the Fermi-Dirac functions. Helps the convergence for some systems. The value is input in Hartree atomic units. 


---
### xyz_file

*Optional* 

Type: characters 

Default: None 

Specifies the location of the xyz file that contains the atomic positions. It can be used as an alternate route to set atomic coordinate. 




*Generated by input_variables.py on 15 April 2021* 


