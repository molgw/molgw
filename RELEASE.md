-----------------------------------------
#    MOLGW: Release Notes
-----------------------------------------

-----------------------------------------
## What's new in version 2.C.beta
### Overview
- speed-up in the Hartree, Exchange and AO to MO transform
- bug fixes

### Contributors
- Fabien Bruneval (CEA SRMP, France)
- Ivan Maliyov (CEA SRMP, France)

### Changes affecting the usage
- Post-processing are not performed if the SCF cycles are not converged within `tolscf` (save user CPU time when a job went wrong)
- Keywords `scalapack_nprow` and `scalapack_npcol` have been eliminated


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
