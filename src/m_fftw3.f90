!=========================================================================
! This file is part of MOLGW.
! Authors: Fabien Bruneval
!
! This module contains
!
!=========================================================================
#include "molgw.h"
#if !defined(NO_LIBINT)
#include<libint2/libint2_params.h>
#end if


module m_fftw3
  use, intrinsic :: ISO_C_BINDING

#if defined(HAVE_FFTW3)
#include"fftw3.f03"
#end if


end module m_fftw3


!=========================================================================
