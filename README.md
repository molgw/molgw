-----------------------------------------
#                 MOLGW
-----------------------------------------


Many-body perturbation theory for atoms, molecules, and clusters


## Getting started

This is a minimalistic README file.
Many more details can be found on the web site [molgw.org](http://www.molgw.org/).
A tutorial section exists there.


## Features

MOLGW implements the following schemes:
- Hartree-Fock
- LDA (PW, VWN)
- GGA (PBE, PW91, BLYP)
- potential-only meta-GGA (BJ, RPP)
- hybrid functionals (PBE0, B3LYP)
- double-hybrid functionals (PBE0-DH, PBE-QIDH, B2PLYP, and RSX-QIDH)
- screened hybrid functionals (HSE03, HSE06, CAM-B3LYP, and LC-BLYP)
- any user-developped range-separated hybrid based on wPBEH
- GW@HF or GW@DFT
- GW/PT2 density-matrix
- QSGW
- QSMP2
- MP2@HF (MP2 correlation energy used in double-hybrid DFT functionals)
- PT2@HF or PT2@DFT
- PT3@HF or PT3@DFT
- CI for few electrons 
- Linear-response TDDFT or TDHF
- Bethe-Salpeter equation
- real-time TDDFT
- HF/DFT+NOFT (MULLER, CGA, CA, GU, POWER, PNOF5, PNOF7, and GNOF)
- Molecules subject to external finite electric fields

The Python3 module `molgw.py` is available for automation.

## Installation

MOLGW needs Fortran 2003 (and a few Fortran2008 features) and C++ compilers.
MOLGW is being tested with `gfortran`, `g++` (version 11.x.x) and `ifort` (version 21).
MOLGW can run in parallel using OPENMP and MPI parallelization.

All the machine dependent variables should be set in file `~molgw/src/my_machine.arch`
Examples for this file can be found in the folder `~molgw/config/`.
Then
`cd ~molgw/src`
`make`

- BLAS and LAPACK linear algebra libraries are required.
- [LIBINT](https://github.com/evaleev/libint/releases) or [LIBCINT](https://github.com/sunqm/libcint/releases) is required for Gaussian integrals
- [LIBXC](https://www.tddft.org/programs/libxc/download/) is required for DFT calculations (else only HF is available)

To run on multi-node computers
- MPI and SCALAPACK are both required


## Basis sets

Many standard Gaussian basis sets are shipped with MOLGW.

More basis sets can be obtained from [Basis Set Exchange](https://bse.pnl.gov/bse/portal)
The file can be generated from a NWChem file using the script
`~molgw/utils/basisset_nwchem2molgw.py aug-cc-pVDZ.nwchem`

You may even create your own.


## Usage

`/path/to/molgw/molgw helium.in > helium.out`

Many example input files can be found in `~molgw/tests/inputs/`


## Known issues

- QSGW scf loop might be quite unstable for large basis sets, use a large eta
- TDDFT GGA kernel can induce very large numerical values that hinders the numerical stability and breaks some comparison with other codes.


## Bug reporting

Please use the [issues](https://github.com/bruneval/molgw/issues) section on MOLGW github.


## Information for developers

Besides the wrapper calls to the LIBINT library, MOLGW is entirely written in Fortran2003/2008.
Fortran C bindings are used to call LIBXC and LIBCINT.
The source files can be found in ~molgw/src/.

### Coding Rules 

The Fortran intent in/out/inout is compulsory for the arguments of a subroutine.
One character variable names are discouraged.

The careful developer should try
- to follow the overall layout and the conventions of the code (double space indent, separation of the list of variables arguments/local, loop counters naming, etc.)
- to protect the data contained in a module with private or protected attribute as much as possible.
- to avoid cascading object access, such as a%b%c (Create methods instead)
- to hide the MPI statements with a generic wrapper in subroutine src/m_mpi.f90.
- to hide the SCALAPACK statements with a generic wrapper in subroutine src/m_scalapack.f90 (not implemented as of today).

### Automatically generated files

A few fortran source files are generated by python scripts:
- src/basis_path.f90
- src/revision.f90 
are generated by src/prepare_sourcecode.py (that is run at each "make" operation)
and
- src/input_variables.f90
is generated by utils/input_variables.py from a YAML file src/input_variables.yaml .
Do not attempt to edit the fortran files. You should rather edit the yaml file.

To add a new input variable, append a new variable description in the YAML file src/input_variables.yaml.
Then execute the python script utils/input_variables.py.
This will generate automatically the Fortran source file src/input_variables.f90
and the HTML and markdown documentation files docs/input_variables.html docs/input_variables.md.

### Adding a new source file

It requires the manual editing of the src/Makefile (sorry).
Please check carefully the module dependence so to compile and add it to the right "level" of the Makefile.
The code should compile properly in parallel with "make -j".


## Contributors

- Fabien Bruneval
- Ivan Maliyov
- Mauricio Rodriguez-Mayorga 
- Xixi Qi
- Young-Moo Byun
- Meiyue Shao
