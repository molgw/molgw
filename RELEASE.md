-----------------------------------------
#    MOLGW: Release Notes
-----------------------------------------

## What's new in version 2.A
### Overview


### Contributors
Ivan Maliyov (CEA SRMP, France)
Young-Moo Byun (U. Illinois@Chicago, USA)
Fabien Bruneval (CEA SRMP, France)

### Changes affecting the results
- The graphical solution of the QP equation now selects the fixed point energy that maximises the spectral weight of the peak.
This may change some QP energies away from the HOMO-LUMO gap. As the consequence, the BSE excitation energies are slightly affected.

### Changes affecting the usage
- OPENMP parallelization is nwo available for the no-RI part with a special care about the reduction of the NUMA effect
- Hybrid parallelization is now available for most the RI part
Please export the corresponding environnement variables.
For instance
```
> export OMP_NUM_THREADS=4
> export MKL_NUM_THREADS=4
```
MOLGW has been shown to run efficiently on Intel KNL architecture

### Changes affecting the compilation
- Possibility to change the diagonalization subroutine with -DLAPACK_DIAGO_FLAVOR_D or -DLAPACK_DIAGO_FLAVOR_R.
In our experience, LAPACK_DIAGO_FLAVOR_R is robust and quick.


-----------------------------------------
## What's new in version 1.H
### Overview
Bug fixes

### Changes affecting the compilation
- MOLGW is now compatible with LIBINT versions 2.4.x
- Possibility to compile MOLGW with LIBINT having the one-body integrals, but not the gradients.
Use the preprocessor flags -DHAVE_LIBINT_ONEBODY and/or -DHAVE_LIBINT_GRADIENTS


## What's new in version 1.G
### Overview
Bug fixes, cleaning, and speed-up.

### Changes affecting the compilation
- Still not possible to link with the versions 4 of LIBXC, due to missing functions on their side

### Changes affecting the usage
- The default value for **read_restart** is now set to 'no'
- Speed-up in the semi-local DFT Vxc calculations thanks to the use of batches


## What's new in version 1.F
### Overview
A few bugs specific to the recent versions of the Intel compilers have been fixed.


## What's new in version 1.E
### Overview
Bug fix with respect to previous version for high angular momenta (L>5)
Considerable speed-up in the diagonalization of the RPA matrix thanks to the use of PDSYEVR instead of PDSYEV


## What's new in version 1.D
### Overview
Simple bug fix with respect to previous version


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
