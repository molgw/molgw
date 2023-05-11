-----------------------------------------
#    MOLGW: Release Notes
-----------------------------------------

----------------------------------------
## What's new in version 3.2
### Overview
- Double-hybrid functionals

### Contributors
- Fabien Bruneval (SRMP, CEA, Université Paris-Saclay, France)
- Mauricio Rodriguez-Mayorga (Universitat d'Alacant, Alicante, Spain)

### Changes affecting the usage

### Changes affecting the compilation
- Makefile has been standardized to help spack packaging
- `HAVE_MPI` automatically implies  `HAVE_SCALAPACK` and vice-versa

### Changes affecting the developers


----------------------------------------
## What's new in version 3.1
### Overview
- Simplified compilation
- Even-tempered basis
- RPA, RPA+, RPAx (RPAx-II), RPAx-I correlation energies
- Double-hybrid functionals (e.g. B2PLYP, PBE0-DH, PBE-QIDH, etc.)
- Inclusion of more RDMFT functionals (e.g. CGA, CA, GU, and GNOF)
- GTH pseudo potentials in CP2K format
- New functionalities in `molgw.py` for running and post-processing series of calculations

### Contributors
- Fabien Bruneval (SRMP, CEA, Université Paris-Saclay, France)
- Mauricio Rodriguez-Mayorga (Universitat d'Alacant, Alicante, Spain)

### Changes affecting the usage
- Even-tempered basis sets can be set with input variables: `even_tempered_alpha`, `even_tempered_beta`, `even_tempered_n_list`.
- RPA, RPA+, RPAx (RPAx-II) are triggered with `postscf` values: `RPA`, `RPA+`, `RPAx-II`, `RPAx-I`
- Double-hybrids (e.g. B2PLYP, PBE0-DH, PBE-QIDH, etc.) require the postscf='MP2' keyword and the amount of EXX and MP2 correlation
- GNOF has become the default one in RDMFT calculations.
- GTH pseudo potentials are assumed when the ECP root name in `ecp_type` contains "GTH" or "gth".

### Changes affecting the compilation
- Compilation is by default with LIBXC and LIBCINT. Use preprocessor variable `-DNO_LIBXC` or `-DNO_LIBCINT` if you want to do otherwise.
Do not use `-DHAVE_LIBCINT` any longer, else a series of warnings may be issued.
- A global Makefile has been created in the root.

### Changes affecting the developers
- Use `#if !defined(NO_LIBXC)` instead of `#if defined(HAVE_LIBXC)`


-----------------------------------------
## What's new in version 3.0
### Overview
- LIBCINT library as an alternative to LIBINT
- Natural Orbital Functional Theory (NOFT)
- GW calculations with W from TDDFT
- pseudopotential in numerical format PSP6 or PSP8 can be used 
- LIBXC functionals can be called directly if one knows their unique LIBXC index
- expansion of the python utilities in molgw.py
- new basis functions (Dunning 7Z)
- various bug fixes, typos

### Contributors
- Fabien Bruneval (SRMP, CEA, Université Paris-Saclay, France)
- Mauricio Rodriguez-Mayorga (Vrije Universiteit, Amsterdam, Netherlands)
- Nike Dattani (HPQC Labs, Waterloo, Canada)
- Zohreh Hashemi (University of Bayreuth, Germany)

### Changes affecting the usage
- Natural orbital functional theory (a.k.a. reduced density matrix functional theory) approximations are now available for singlet-states.
NOFT calculations are triggered with `postscf='NOFT'`.
The corresponding input variables start with `noft_`. Available functionals include: Muller, power, PNOF5, and PNOF7.
- Dynamical self-energy correlation contribution and Spectral weight Z are now reported for omega=E_qp instead of omega=E_gKS.
- `postscf='GWTDDFT'` triggers the calculation of W with HF or TDDFT kernel included. 
- MOLGW can use "solid-state" norm-conserving pseudopotentials in the PSP6 or PSP8 format (from [pseudo-dojo.org](http://www.pseudo-dojo.org/) for instance).
This is slow but functional. It is intended for tests.
- Python script to extract the basis set from a Gaussian formatted checkpoint file (.fchk): `create_basis_from_gaussian_fchk.py`.
- All the LIBXC functionals can be called directly with syntax: `scf='LIBXC:101+130'` for PBE for instance. 
- Possibility to tune the cube file output with new input variables `cube_nx`, `cube_ny`, `cube_nz`, `cube_state_min`, `cube_state_max`.
- Possibility to output the transition density in BSE/TDDFT.
- Possibility to read Gaussian cube files with molgw.py

### Changes affecting the compilation
- The C library [LIBCINT](https://github.com/sunqm/libcint) can replace LIBINT with noticeable advantages.
LIBCINT is easy to compile on any architecture (it's coded in plain C) and its compilation time is low.
According to our tests, LIBCINT appears as twice faster than LIBINT at runtime.
To use LIBCINT, `my_machine.arch` should contain `LIBCINT=-lcint` and preprocessing option `-DHAVE_LIBCINT`

### Changes affecting the developers
- LIBCINT library offers many integral types that were not available with LIBINT. New opportunities are open.


-----------------------------------------
## What's new in version 2.F
### Overview
- MOLGW is now compatible with LIBXC 5
- MOLGW automatically detects the LIBINT configuration. Easier compilation
- Possibility to add point charges in the structure

### Contributors
- Fabien Bruneval (CEA SRMP, France)

### Changes affecting the usage
- Fractional point charges (without basis functions) can be specified in the structure using the syntax:
 0.100   0.000  0.000  0.000   none  none   #   q   x y z  basis  auxiliary_basis

### Changes affecting the compilation
- MOLGW can be linked against LIBXC 5
- MOLGW detects LIBINT configuration to know wheter the one-body integrals and the gradients are available. Preprocessor instructions such as `-DHAVE_LIBNIT_ONEDOBY` are not needed anymore.


-----------------------------------------
## What's new in version 2.E
### Overview
- MOLGW proposes automatically an extrapolated GW energy to the  Complete Basis Set limit when using Dunning basis sets
- GW with analytic continuation is now robust for the HOMO-LUMO gap region. Tested for C60 in aug-cc-pV5Z (>7500 basis functions)
- small bug fixes, speed-ups, memory reductions

### Contributors
- Fabien Bruneval (CEA SRMP, France)
- Xixi Qi (CEA SRMP, France)
- Mauricio Rodriguez-Mayorga (CEA SRMP, France)

### Changes affecting the usage
- Automatic suggestion of an extrapolation to the Complete Basis Set limit for GW energies
- GW analytic continuation technique is fully functional. Use postscf='G0W0_pade'. It is much faster than the analytic formula but it is mostly reliable close to the HOMO-LUMO gap.
- Reduced memory consumption in the Pulay history (SCALAPACK distribution of the large arrays)
- Improved OPENMP

### Changes affecting the compilation
- Assuming now that all Fortran compilers have Fortran 2008 capabilities. Preprocessor key FORTRAN2008 has been removed.

### Changes affecting the developers
- Introduce high-level mpi tools for instance, world%sum() for reduction, world%nproc, world%rank for information


-----------------------------------------
## What's new in version 2.D
### Overview
- Compatibility with gcc/gfortran 10
- Basis files location can be set from an environment variable MOLGW_BASIS_PATH
- Printing of standard WFN files

### Contributors
- Fabien Bruneval (CEA SRMP, France)
- Mauricio Rodriguez-Mayorga (CEA SRMP, France)

### Changes affecting the usage
- Keyword `print_wfn_file` triggers the output of a standard WFN file that can be read with external visualization softwares.
- Environment variable `MOLGW_BASIS_PATH` sets the path to the basis files. It is still be overridden by the input keyword `basis_path`.
- New default value for `postscf_diago_flavor=' '`. Though faster, the former default value was not stable enough for large systems.

### Changes affecting the compilation
- GCC 10 is very picky on the routine calls without an interfaces. Many existing calls to BLAS/LAPACK/SCALAPACK had to be fixed.
- Makefile, my_machine.arch use more standard `FCFLAGS` and `CXXFLAGS` variables instead of `FCOPTS` and `CXXOPTS`
- Fortran long lines have been chopped into pieces so to comply with the 132 character limit of Fortran.
Compiler options such as `-ffree-line-length-none` are not needed any more.

### Changes affecting the developers
- Please respect the 132-character limit of Fortran.


-----------------------------------------
## What's new in version 2.C
### Overview
- Real-time TDDFT is made available
- speed-up in the Hartree, Exchange and AO to MO transform
- calculation of the generalized oscillator strengths (q-dependent) and linear-response stopping power
- use of LIBXC through the C interface
- compatibility with the latest LIBINT versions restored
- creation of a YAML output file that gathers many information that can be easily post-processed via python
- bug fixes

### Contributors
- Fabien Bruneval (CEA SRMP, France)
- Ivan Maliyov (CEA SRMP, France)
- Young-Moo Byun (U. Illinois@Chicago, USA)

### Changes affecting the usage
- Post-processing is not performed if the SCF cycles are not converged within `tolscf` (save user CPU time when a job went wrong)
- Keywords `scalapack_nprow` and `scalapack_npcol` have been eliminated
- Keyword `stopping` triggers the linear-response stopping power calculation
- Value `postscf='real_time'` triggers real-time TDDFT (RT-TDDFT)
- New default value for `postscf_diago_flavor='D'`
- OPENMP parallelisation effort is pursued
- Bug fix of the output of the COHSEX energies (bug discovered by Arjan Berger)

### Changes affecting the compilation
- LIBXC is now linked through the C interface. Therefore, LIBXC compilation does not need to be consistent with MOLGW compilation.
The latest version of LIBXC can be used. Preprocessing flag `-DLIBXC4` is no longer needed.
- Preprocessing option `-DHAVE_MKL` allows for the use of MKL extensions and in particular of `DGEMMT`.
- Use of modern Fortran2008 syntax, such as c%re to obtain the real part of a complex number. Code imcompatible with older compilers (gfortran > 9.0 necessary)

### Changes affecting the developers
- The list of all the input variables is now stored in a YAML file ~molgw/src/input_variables.yaml that is processed with the python script ~molgw/utils/input_variables.py


-----------------------------------------
## What's new in version 2.B
### Overview
- automatic generation of an auxiliary basis following the "Auto" and "PAuto" recipes of Gaussian
- output the Galitskii-Migdal correlation energy = 1/2 Tr[ Sigmac * G ]
- non HF starting points for perturbation theory density matrices
- speed-up in RT-TDDFT
- bug fixes

### Contributors
- Fabien Bruneval (CEA SRMP, France)
- Young-Moo Byun (U. Illinois@Chicago, USA)
- Ivan Maliyov (CEA SRMP, France)

### Changes affecting the usage
- Keyword `auxil_basis='auto'` or `auxil_basis='pauto'` triggers automatic generation of an auxiliary basis


-----------------------------------------
## What's new in version 2.A
### Overview
- GW approximation to the density matrix (only for spin-restricted calculations)
- Third-order perturbation theory (PT3) self-energy (only for spin-restricted calculations)
- OPENMP parallelization
- better graphical solution to the QP equation
- complex wavefunctions and Hamiltonian can be calculated (for real-time TDDFT)
- possibility to use different diagonalization routines
- reduced memory foot print in the 3-center integral set up
- possibility to read formatted checkpoint files from Gaussian (.fchk) and use the read density matrix
- speed-up
- bug fixes

### Contributors
- Fabien Bruneval (CEA SRMP, France)
- Young-Moo Byun (U. Illinois@Chicago, USA)
- Ivan Maliyov (CEA SRMP, France)

### Changes affecting the results
- The graphical solution of the QP equation now selects the fixed point energy that maximises the spectral weight of the peak.
This may change some QP energies away from the HOMO-LUMO gap. As the consequence, the BSE excitation energies are slightly affected.

### Changes affecting the usage
- Keyword `pt_density_matrix` triggers the GW or PT2 density matrix calculation
- Keyword `postscf='PT3'` triggers the PT3 self-energy calculation [see Cederbaum's papers](https://doi.org/10.1016/0167-7977%2884%2990002-9)
- Keyword `read_fchk` allows the user to read a Gaussian formatted file (.fchk)
and resuse the density matrix `read_fchk='SCF'`, `read_fchk='MP2'`, and `read_fchk='CC'`. Only works for Cartesian Gaussian and make sure to use the exact same basis set in the two codes. Gaussian uses modification of the Dunning basis set for instance. A few are available in MOLGW with `basis='cc-pVQZ_gaussian'`.
- OPENMP parallelization is now available for the no-RI part with a special care about the reduction of the NUMA effect.
- Hybrid parallelization (OPENMP/MPI) is now available for most the RI code. The coding heavily relies on the availability of threaded BLAS routines.
When using Intel MKL, please export the corresponding environment variables:
`export OMP_NUM_THREADS=4` and `export MKL_NUM_THREADS=4`.
MOLGW has been shown to run efficiently on Intel KNL architecture
- Memory foot print can be reduced when the memory peak was in the calculation of the 3-center integrals.
Use keyword `eri3_nbatch`. Setting `eri3_nbatch` to a value larger than 1 will decrease the memory consumption in the 3-center integrals calculation and hopefully will not hit too much on the performance.
- Possibility to change the diagonalization subroutines with input variables `scf_diago_flavor` and `postscf_diago_flavor`.
Based on our experience, flavor 'R' pointing to (P)DSYEVR is faster but can induce SCF cycle issues. Flavor ' ' pointing to (P)DYSEV is the safest possibility.

### Changes affecting the compilation
- Fortran compiler needs to be capable of using the polymorphic `class(*)` declarations (Fortran2003).
- `-DDEBUG` produces an outfile for each MPI thread.


-----------------------------------------
## What's new in version 1.H
### Overview
Bug fixes

### Changes affecting the compilation
- MOLGW is now compatible with LIBINT versions 2.4.x
- Possibility to compile MOLGW with LIBINT having the one-body integrals, but not the gradients.
Use the preprocessor flags -DHAVE_LIBINT_ONEBODY and/or -DHAVE_LIBINT_GRADIENTS


-----------------------------------------
## What's new in version 1.G
### Overview
Bug fixes, cleaning, and speed-up.

### Changes affecting the compilation
- Still not possible to link with the versions 4 of LIBXC, due to missing functions on their side

### Changes affecting the usage
- The default value for **read_restart** is now set to 'no'
- Speed-up in the semi-local DFT Vxc calculations thanks to the use of batches


-----------------------------------------
## What's new in version 1.F
### Overview
A few bugs specific to the recent versions of the Intel compilers have been fixed.


-----------------------------------------
## What's new in version 1.E
### Overview
Bug fix with respect to previous version for high angular momenta (L>5)
Considerable speed-up in the diagonalization of the RPA matrix thanks to the use of PDSYEVR instead of PDSYEV


-----------------------------------------
## What's new in version 1.D
### Overview
Simple bug fix with respect to previous version


-----------------------------------------
## What's new in version 1.C
### Overview
This release makes better use of the latest version of LIBINT ( >= v2.2.0).
Together with some internal refactoring of the code and small bug fixes.

### Changes affecting the compilation
- Linking with a recent version of LIBINT (>=v2.2.0) is compulsory
- Compilation flag **-DHAVE_LIBINT_ONEBODY** activates the use of LIBINT one-body integrals and integral gradients
- Use of libtool is no longer necessary

### Changes affecting the usage
- Input variable **no_4center** has been removed.
From now on, the calculations use or do not use resolution-of-identity from the beginning to end.
