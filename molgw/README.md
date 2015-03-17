=========================================
#                 MOLGW
=========================================

Many-body perturbation theory and (TD)-DFT calculations for atoms and small molecules

##Features
- Hartree-Fock
- LDA (PW, VWN)
- GGA (PBE, PW91, BLYP)
- potential-only meta-GGA (BJ, RPP)
- hybrid functionals (PBE0, B3LYP)
- screened hybrid functionals (HSE03, HSE06)
- HF+GW
- DFT+GW
- Hybrid+GW
- QPscGW
- HF+MP2
- DFT+MP2
- Hybrid+MP2
- QPscMP2
- CI for 2 electrons 
- TD-HF
- TD-DFT
- BSE


##Installation

MOLGW needs Fortran 2003, c and c++ compilers.
The machine dependent variables should be set in file `~molgw/src/my_machine.arch`
Examples for this file are given in the folder `~molgw/config/`.
Then
`cd ~molgw/src`
`make`

- BLAS and LAPACK linear algebra libraries are required.
- libint is required:
http://sourceforge.net/projects/libint/files/libint-for-beginners/libint-2.0.0-prealpha.tgz/download
- libxc is required: (version > 2.0.1) for DFT calculations
http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.0.1.tar.gz


##Basis sets
Basis sets can be obtained from https://bse.pnl.gov/bse/portal
The file can be generated from a NWChem file using the script
`~molgw/util/basis_nwchem2molgw.py B_aug-cc-pVDZ.nwchem`


##Usage
The basis file needs to be located in the working directory.

`./molgw < helium.in > helium.out`

Example input files can be found in `~molgw/tests/`


##To be done
- pseudopotentials (e.g. ECP or Goedecker flavor)
- parallelization!

##Known issues
- QPscGW scf loop is quite unstable for large basis sets, use a low alpha (<= 0.50), use a large eta
- TD-DFT GGA kernel can induce very large numerical values which limits the numerical stability and breaks some comparison with other codes.
Especially when compiling with gfortran/gcc. ifort/icc behaves much better.

##Author

Fabien Bruneval

Service de Recherches de Métallurgie Physique
CEA Saclay, F-91191 Gif-sur-Yvette, France
