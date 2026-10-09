-----------------------------------------
#                 MOLGW
-----------------------------------------


Many-body perturbation theory for atoms, molecules, and clusters


## Getting started

This is a minimalistic README file.
Many more details can be found on the website [molgw.org](http://www.molgw.org/).
A tutorial is available there.

See also:

- `INSTALL.md` for compilation tips
- `RELEASE.md` for the release notes
- `docs/input_variables.md` for the list of all input variables


## Features

MOLGW implements the following schemes:

- Hartree-Fock
- LDA (PW, VWN)
- GGA (PBE, PW91, BLYP)
- potential-only meta-GGA (BJ, RPP)
- hybrid functionals (PBE0, B3LYP)
- double-hybrid functionals (PBE0-DH, PBE-QIDH, B2PLYP, and RSX-QIDH)
- screened hybrid functionals (HSE03, HSE06, CAM-B3LYP, and LC-BLYP)
- any user-developed range-separated hybrid based on wPBEH
- GW@HF or GW@DFT
- GW+SOSEX and GW+2SOSEX (including a positive semi-definite version)
- GW/PT2 density-matrix
- QSGW
- QSMP2
- MP2@HF (MP2 correlation energy used in double-hybrid DFT functionals)
- PT2@HF or PT2@DFT
- PT3@HF or PT3@DFT
- CI for a few electrons
- Linear-response TDDFT or TDHF
- Bethe-Salpeter equation
- real-time TDDFT
- X2C relativistic HF/DFT calculations
- HF/DFT Bogoliubov calculations
- HF/DFT+NOFT (MULLER, CGA, CA, GU, POWER, PNOF5, PNOF7, and GNOF)
- Molecules subject to external finite electric fields
- Forces and geometry optimization
- Interface with CC4S (export and import of Coulomb integrals and HF eigenvalues)

The Python3 package `molgw` is available for automation in `~molgw/utils/molgw/`.
Install it with `pip install -e utils/molgw/`

## Installation

MOLGW needs Fortran 2008 and C++ compilers.
MOLGW is routinely tested with `gfortran`/`g++` and `flang`/`clang` (version 19). Intel `ifx` is also supported.
MOLGW can run in parallel using OPENMP and MPI parallelization.

All the machine dependent variables should be set in file `~molgw/my_machine.arch`
Examples for this file can be found in the folder `~molgw/config/`.
Then
`cd ~molgw`
`make`

- BLAS and LAPACK linear algebra libraries are required.
- [LIBCINT](https://github.com/sunqm/libcint/releases) is required for Gaussian integrals ([LIBINT](https://github.com/evaleev/libint/releases) is still supported but deprecated)
- [LIBXC](https://www.tddft.org/programs/libxc/download/) is required for DFT calculations (else only HF is available)

To run on multi-node computers

- MPI and SCALAPACK are both required

Optional libraries

- [ELPA](https://elpa.mpcdf.mpg.de/) for faster diagonalizations on CPU and GPU
- [FFTW3](https://www.fftw.org/) for fast Fourier transforms
- [HDF5](https://www.hdfgroup.org/) for HDF5 output of some matrices (see `~molgw/HDF5_INSTALL.md`)
- [QCINT](https://github.com/sunqm/qcint) can replace LIBCINT for faster integrals


## Basis sets

Many standard Gaussian basis sets are shipped with MOLGW.

More basis sets can be obtained from [Basis Set Exchange](https://www.basissetexchange.org)
The file can be generated from a NWChem file using the script
`~molgw/utils/basisset_nwchem2molgw.py aug-cc-pVDZ.nwchem`
or from a CP2K file using the script `~molgw/utils/basisset_cp2k2molgw.py`

You may even create your own.


## Usage

`/path/to/molgw/molgw helium.in > helium.out`

Many example input files can be found in `~molgw/tests/inputs/`


## Testing

Run the test suite from the root folder:

- `make test` (sequential)
- `make test NC=4` (4 OpenMP threads per test)
- `make test NP=4` (4 MPI processes per test)

Or run a selection of tests from `~molgw/tests/`:

- `python3 run_testsuite.py --list` to list all the tests
- `python3 run_testsuite.py --input be_qsgw.in` to run a single test
- `python3 run_testsuite.py --keyword GW BSE` to run tests by keyword

The tests and their reference values are defined in `~molgw/tests/inputs/testsuite.yaml`.


## Known issues

- The QSGW SCF loop can be quite unstable for large basis sets; use a large eta.
- The TDDFT GGA kernel can induce very large numerical values that hinder numerical stability and break comparisons with other codes.


## Bug reporting

Please use the [issues](https://github.com/molgw/molgw/issues) section on MOLGW github.


## Information for developers

MOLGW is entirely written in Fortran 2008.
Fortran C bindings are used to call LIBCINT and LIBXC.
The source files can be found in `~molgw/src/`.

### Coding Rules

Fortran intent (in/out/inout) declarations are required for all subroutine arguments.
Single-character variable names are discouraged.

The careful developer should try

- to follow the overall layout and the conventions of the code (double space indent, separation of the list of variables arguments/local, loop counters naming, etc.)
- to protect the data contained in a module with private or protected attribute as much as possible.
- to avoid cascading object access, such as `a%b%c` (Create methods instead)
- to hide the MPI statements with a generic wrapper in module `src/m_mpi.f90`.
- to hide the SCALAPACK statements with a generic wrapper in module `src/m_scalapack.f90` (not yet fully implemented).

### Automatically generated files

A few Fortran source files are generated by Python scripts:

- `src/basis_path.f90`
- `src/git_sha.f90`

are generated by `src/prepare_sourcecode.py` (run at each `make` operation),
`src/noft/gitver.F90` is generated by `src/noft/gitversion.sh` (run at each `make` operation),
and

- `src/input_variable_declaration.f90`
- `src/input_variables.f90`
- `src/echo_input_variables.f90`
- `src/echo_input_variables_yaml.f90`
- `src/input_variables_logical.f90`

are generated by `utils/input_variables.py` from the YAML file `src/input_variables.yaml`.
Do not attempt to edit the Fortran files directly; edit the YAML file instead.

To add a new input variable, append a new variable description in the YAML file `src/input_variables.yaml`.
Then execute the Python script `utils/input_variables.py` (or `make input` in `src/`).
This will generate automatically the Fortran source file `src/input_variables.f90`
and the HTML and markdown documentation files `docs/input_variables.html` and `docs/input_variables.md`.

### Adding a new source file

It requires the manual editing of `src/Makefile` (sorry).
Please check the module dependencies carefully and add the new file to the correct level in the Makefile.
The code should compile properly in parallel with `make -j`.


## Contributors

- Fabien Bruneval
- Ivan Maliyov
- Mauricio Rodriguez-Mayorga
- Xixi Qi
- Young-Moo Byun
- Meiyue Shao
- Marc Sarraute
- Guillaume Piris
- Caroline McKeon
