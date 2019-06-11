!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to interface with the C library LIBXC
! * explicit interfaces with C binding
!
!=========================================================================
module m_libxc_tools
 use m_definitions

#if defined(HAVE_LIBXC)

#include <xc_funcs.h>

 integer(C_INT),parameter :: XC_FAMILY_UNKNOWN  = -1
 integer(C_INT),parameter :: XC_FAMILY_LDA      =  1
 integer(C_INT),parameter :: XC_FAMILY_GGA      =  2
 integer(C_INT),parameter :: XC_FAMILY_MGGA     =  4
 integer(C_INT),parameter :: XC_FAMILY_LCA      =  8
 integer(C_INT),parameter :: XC_FAMILY_OEP      = 16
 integer(C_INT),parameter :: XC_FAMILY_HYB_GGA  = 32
 integer(C_INT),parameter :: XC_FAMILY_HYB_MGGA = 64

 integer(C_INT),parameter :: XC_FLAGS_HAVE_FXC  =  4

!=========================================================================
! interfaces to libxc in C with iso_c_binding
!=========================================================================

 interface

   function xc_func_type_malloc() BIND(C)
     import :: C_PTR
     type(C_PTR) :: xc_func_type_malloc
   end function xc_func_type_malloc

   subroutine xc_version(major,minor,micro) BIND(C)
     import :: C_INT
     integer(C_INT),intent(out) :: major,minor,micro
   end subroutine xc_version

   function xc_func_init(p,functional,nspin) BIND(C)
     import :: C_PTR,C_INT
     type(C_PTR),intent(inout) :: p
     integer(C_INT),value      :: functional,nspin
     integer(C_INT)            :: xc_func_init
   end function xc_func_init

   subroutine xc_func_end(p) BIND(C)
     import :: C_PTR
     type(C_PTR),intent(inout) :: p
   end subroutine xc_func_end

   function xc_func_get_info(p) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     type(C_PTR),intent(in)    :: p
     type(C_PTR)               :: xc_func_get_info
   end function xc_func_get_info

   function xc_func_info_get_name(info) BIND(C)
     import :: C_PTR,C_CHAR
     type(C_PTR),intent(in)    :: info
     character(KIND=C_CHAR)    :: xc_func_info_get_name
   end function xc_func_info_get_name

   function xc_func_info_get_family(info) BIND(C)
     import :: C_PTR,C_INT
     type(C_PTR),intent(in)    :: info
     integer(C_INT)            :: xc_func_info_get_family
   end function xc_func_info_get_family

   function xc_func_info_get_flags(info) BIND(C)
     import :: C_PTR,C_INT
     type(C_PTR),intent(in)    :: info
     integer(C_INT)            :: xc_func_info_get_flags
   end function xc_func_info_get_flags

   subroutine xc_func_set_ext_params(p,ext_params) BIND(C)
     import :: C_PTR,C_DOUBLE
     type(C_PTR),intent(in)    :: p
     real(C_DOUBLE),intent(in) :: ext_params(*)
   end subroutine xc_func_set_ext_params

   subroutine xc_lda_exc_vxc(p,np,rho,exc,vrho) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     type(C_PTR),intent(in)     :: p
     integer(C_INT),value       :: np
     real(C_DOUBLE),intent(in)  :: rho(*)
     real(C_DOUBLE),intent(out) :: exc(*)
     real(C_DOUBLE),intent(out) :: vrho(*)
   end subroutine xc_lda_exc_vxc

   subroutine xc_lda_fxc(p,np,rho,v2rho2) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     type(C_PTR),intent(in)     :: p
     integer(C_INT),value       :: np
     real(C_DOUBLE),intent(in)  :: rho(*)
     real(C_DOUBLE),intent(out) :: v2rho2(*)
   end subroutine xc_lda_fxc

   subroutine xc_gga_exc_vxc(p,np,rho,sigma,exc,vrho,vsigma) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     type(C_PTR),intent(in)     :: p
     integer(C_INT),value       :: np
     real(C_DOUBLE),intent(in)  :: rho(*)
     real(C_DOUBLE),intent(in)  :: sigma(*)
     real(C_DOUBLE),intent(out) :: exc(*)
     real(C_DOUBLE),intent(out) :: vrho(*)
     real(C_DOUBLE),intent(out) :: vsigma(*)
   end subroutine xc_gga_exc_vxc

   subroutine xc_gga_vxc(p,np,rho,sigma,vrho,vsigma) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     type(C_PTR),intent(in)     :: p
     integer(C_INT),value       :: np
     real(C_DOUBLE),intent(in)  :: rho(*)
     real(C_DOUBLE),intent(in)  :: sigma(*)
     real(C_DOUBLE),intent(out) :: vrho(*)
     real(C_DOUBLE),intent(out) :: vsigma(*)
   end subroutine xc_gga_vxc

   subroutine xc_gga_fxc(p,np,rho,sigma,v2rho2,v2rhosigma,v2sigma2) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     type(C_PTR),intent(in)     :: p
     integer(C_INT),value       :: np
     real(C_DOUBLE),intent(in)  :: rho(*)
     real(C_DOUBLE),intent(in)  :: sigma(*)
     real(C_DOUBLE),intent(out) :: v2rho2(*)
     real(C_DOUBLE),intent(out) :: v2rhosigma(*)
     real(C_DOUBLE),intent(out) :: v2sigma2(*)
   end subroutine xc_gga_fxc

end interface

#endif

end module m_libxc_tools


!=========================================================================
