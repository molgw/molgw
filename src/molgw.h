/* !=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the preprocessing options
!
!========================================================================= */

#if defined(HAVE_LIBCINT)
#if !defined(NO_LIBINT)
#define NO_LIBINT
#endif
#endif

#if !defined(NO_LIBXC)
#if !defined(HAVE_LIBXC)
#define HAVE_LIBXC
#endif
#endif



