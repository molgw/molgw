=========================================
#    MOLGW: Release Notes
=========================================


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


