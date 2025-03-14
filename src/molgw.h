#if 0
/*
!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! the preprocessing options
!
!========================================================================= 
*/
#endif

#if !defined(MOLGW_VERSION)
#define MOLGW_VERSION "3.3"
#endif

#if !defined(NO_LIBCINT)
#define HAVE_LIBCINT
#endif

#if defined(HAVE_LIBCINT)
#if !defined(NO_LIBINT)
#define NO_LIBINT
#endif
#endif

#if 0
/*
! Enforce that HAVE_MPI and HAVE_SCALAPACK are both defined when a single one is defined
*/
#endif
#if defined(HAVE_MPI)
#if !defined(HAVE_SCALAPACK)
#define HAVE_SCALAPACK
#endif
#endif

#if defined(HAVE_SCALAPACK)
#if !defined(HAVE_MPI)
#define HAVE_MPI
#endif
#endif

