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
     implicit none
     type(C_PTR) :: xc_func_type_malloc
   end function xc_func_type_malloc

   function xc_func_alloc() BIND(C)
     import :: C_PTR
     implicit none
     type(C_PTR) :: xc_func_alloc
   end function xc_func_alloc

   function get_family_id(func) BIND(C)
     import :: C_PTR,C_INT
     implicit none
     type(C_PTR)    :: func
     integer(C_INT) :: get_family_id
   end function get_family_id  !TODO FBFB remove

   function get_nspin(func) bind(c)
     import :: c_ptr,c_int
     implicit none
     type(c_ptr)            :: func
     integer(c_int)         :: get_nspin
   end function get_nspin  !todo fbfb remove

   function set_nspin(func,nspin) bind(c)
     import :: c_ptr,c_int
     implicit none
     type(c_ptr)               :: func
     integer(C_INT),value      :: nspin
     integer(c_int)            :: set_nspin
   end function set_nspin  !todo fbfb remove

   subroutine xc_version(major,minor,micro) BIND(C)
     import :: C_INT
     implicit none
     integer(C_INT)             :: major,minor,micro
   end subroutine xc_version

   function xc_family_from_id(id,family,n) BIND(C)
     import :: C_INT,C_PTR
     implicit none
     integer(C_INT),value       :: id
     type(C_PTR),value          :: family,n
     integer(C_INT)             :: xc_family_from_id
   end function xc_family_from_id

   function xc_func_init(func,functional,nspin) BIND(C)
     import :: C_PTR,C_INT
     type(C_PTR),intent(inout) :: func
     integer(C_INT),value      :: functional,nspin
     integer(C_INT)            :: xc_func_init
   end function xc_func_init

   subroutine xc_func_end(func) BIND(C)
     import :: C_PTR
     implicit none
     type(C_PTR)               :: func
   end subroutine xc_func_end

   function xc_func_get_info(func) BIND(C)
     import :: C_PTR
     implicit none
     type(C_PTR)               :: func
     type(C_PTR)               :: xc_func_get_info
   end function xc_func_get_info

   function xc_func_info_get_name(info) BIND(C)
     import :: C_PTR,C_CHAR
     implicit none
     type(C_PTR)               :: info
     character(KIND=C_CHAR)    :: xc_func_info_get_name
   end function xc_func_info_get_name

   function xc_func_info_get_family(info) BIND(C)
     import :: C_PTR,C_INT
     implicit none
     type(C_PTR)               :: info
     integer(C_INT)            :: xc_func_info_get_family
   end function xc_func_info_get_family

   function xc_func_info_get_flags(info) BIND(C)
     import :: C_PTR,C_INT
     implicit none
     type(C_PTR)               :: info
     integer(C_INT)            :: xc_func_info_get_flags
   end function xc_func_info_get_flags

   subroutine xc_func_set_ext_params(func,ext_params) BIND(C)
     import :: C_PTR,C_DOUBLE
     implicit none
     type(C_PTR)               :: func
     real(C_DOUBLE)            :: ext_params(*)
   end subroutine xc_func_set_ext_params

   subroutine xc_lda_exc(func,np,rho,exc) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     implicit none
     type(C_PTR)                :: func
     integer(C_INT),value       :: np
     real(C_DOUBLE)             :: rho(*)
     real(C_DOUBLE)             :: exc(*)
   end subroutine xc_lda_exc

   subroutine xc_lda_exc_vxc(func,np,rho,exc,vrho) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     implicit none
     type(C_PTR)                :: func
     integer(C_INT),value       :: np
     real(C_DOUBLE)             :: rho(*)
     real(C_DOUBLE)             :: exc(*)
     real(C_DOUBLE)             :: vrho(*)
   end subroutine xc_lda_exc_vxc

   subroutine xc_lda_fxc(func,np,rho,v2rho2) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     implicit none
     type(C_PTR)                :: func
     integer(C_INT),value       :: np
     real(C_DOUBLE)             :: rho(*)
     real(C_DOUBLE)             :: v2rho2(*)
   end subroutine xc_lda_fxc

   subroutine xc_gga_exc(func,np,rho,sigma,exc) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     implicit none
     type(C_PTR)                :: func
     integer(C_INT),value       :: np
     real(C_DOUBLE)             :: rho(*)
     real(C_DOUBLE)             :: sigma(*)
     real(C_DOUBLE)             :: exc(*)
   end subroutine xc_gga_exc

   subroutine xc_gga_exc_vxc(func,np,rho,sigma,exc,vrho,vsigma) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     implicit none
     type(C_PTR)                :: func
     integer(C_INT),value       :: np
     real(C_DOUBLE)             :: rho(*)
     real(C_DOUBLE)             :: sigma(*)
     real(C_DOUBLE)             :: exc(*)
     real(C_DOUBLE)             :: vrho(*)
     real(C_DOUBLE)             :: vsigma(*)
   end subroutine xc_gga_exc_vxc

   subroutine xc_gga_vxc(func,np,rho,sigma,vrho,vsigma) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     implicit none
     type(C_PTR)                :: func
     integer(C_INT),value       :: np
     real(C_DOUBLE)             :: rho(*)
     real(C_DOUBLE)             :: sigma(*)
     real(C_DOUBLE)             :: vrho(*)
     real(C_DOUBLE)             :: vsigma(*)
   end subroutine xc_gga_vxc

   subroutine xc_gga_fxc(func,np,rho,sigma,v2rho2,v2rhosigma,v2sigma2) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     implicit none
     type(C_PTR)                :: func
     integer(C_INT),value       :: np
     real(C_DOUBLE)             :: rho(*)
     real(C_DOUBLE)             :: sigma(*)
     real(C_DOUBLE)             :: v2rho2(*)
     real(C_DOUBLE)             :: v2rhosigma(*)
     real(C_DOUBLE)             :: v2sigma2(*)
   end subroutine xc_gga_fxc

end interface

#endif

end module m_libxc_tools


!=========================================================================
