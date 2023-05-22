!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to interface with the C library LIBXC
! * explicit interfaces with C binding
!
!=========================================================================
#include "molgw.h"
module m_libxc_tools
 use m_definitions
 use m_warning


#if !defined(NO_LIBXC)
#include <xc_funcs.h>
#include <xc_version.h>
#endif

 integer(C_INT),parameter :: XC_FAMILY_UNKNOWN  = -1
 integer(C_INT),parameter :: XC_FAMILY_LDA      =  1
 integer(C_INT),parameter :: XC_FAMILY_GGA      =  2
 integer(C_INT),parameter :: XC_FAMILY_MGGA     =  4
 integer(C_INT),parameter :: XC_FAMILY_LCA      =  8
 integer(C_INT),parameter :: XC_FAMILY_OEP      = 16
 integer(C_INT),parameter :: XC_FAMILY_HYB_GGA  = 32
 integer(C_INT),parameter :: XC_FAMILY_HYB_MGGA = 64

 integer(C_INT),parameter :: XC_FLAGS_HAVE_FXC  =  4

 ! Object that contains the information of a LIBXC functional
 type dft_xc_info
   logical                    :: needs_gradient         ! TRUE if any of the components is not LDA
   integer                    :: nxc = 0                ! number of components of the Vxc
   integer(C_INT)             :: id                     ! ID number in the LIBXC list
   real(dp)                   :: coeff                  ! coefficient that scales this Vxc component
   real(C_DOUBLE)             :: gamma                  ! range-separation parameter (if needed)
   integer(C_INT)             :: nspin                  ! number of spin channels (can be different from the global nspin)
   integer                    :: family                 ! LDA, GGA, HYB_GGA
   type(C_PTR),pointer        :: func => NULL()         ! pointer to the LIBXC functional C structure
 end type

!=========================================================================
! interfaces to libxc in C with iso_c_binding
!=========================================================================

 interface

   function xc_func_alloc() BIND(C)
     import :: C_PTR
     implicit none
     type(C_PTR) :: xc_func_alloc
   end function xc_func_alloc

   subroutine xc_func_free(func) BIND(C)
     import :: C_PTR
     implicit none
     type(C_PTR) :: func
   end subroutine xc_func_free

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

   subroutine xc_func_set_ext_params_name(func,name,param) BIND(C)
     import :: C_PTR,C_DOUBLE,C_CHAR
     implicit none
     type(C_PTR),intent(inout) :: func
     character(kind=C_CHAR)    :: name
     real(C_DOUBLE),value      :: param
   end subroutine xc_func_set_ext_params_name

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

   subroutine xc_hyb_cam_coef(func,omega,alpha,beta) BIND(C)
     import :: C_PTR,C_INT,C_DOUBLE
     implicit none
     type(C_PTR)                :: func
     real(C_DOUBLE)             :: omega,alpha,beta
   end subroutine xc_hyb_cam_coef

   function xc_functional_get_number(name) BIND(C)
     import :: C_INT,C_CHAR
     integer(C_INT) :: xc_functional_get_number
     character(KIND=C_CHAR),intent(in) :: name
   end function

   function xc_functional_get_name(number) BIND(C)
     import :: C_INT,C_CHAR
     character(KIND=C_CHAR) :: xc_functional_get_name
     integer(C_INT),intent(in) :: number
   end function

end interface


contains


!=========================================================================
subroutine init_libxc_info(dft_xc)
 implicit none

 type(dft_xc_info),allocatable,intent(inout) :: dft_xc(:)
!=====
 integer :: ixc
 type(C_PTR)          :: cptr_tmp
!=====

 if( .NOT. ALLOCATED(dft_xc) ) then
   call die('init_libxc_info: dft_xc_info type should have been allocated already at this stage')
 endif

 dft_xc(:)%nxc = COUNT( dft_xc(:)%id /= 0 )


#if !defined(NO_LIBXC)

 !
 ! Initialize the DFT objects for LIBXC
 !

 do ixc=1,dft_xc(1)%nxc

   cptr_tmp = xc_func_alloc()
   call c_f_pointer(cptr_tmp,dft_xc(ixc)%func)

   if( xc_func_init(dft_xc(ixc)%func,dft_xc(ixc)%id,dft_xc(ixc)%nspin) /= 0 ) then
     write(stdout,'(1x,a,i6)') 'Libxc failure when initializing functional: ',dft_xc(ixc)%id
     call die('init_libxc_info: error in LIBXC xc_func_init')
   endif

   !
   ! Tune the range for range separated hybrids
   if( dft_xc(ixc)%id == XC_GGA_X_HJS_PBE .or. dft_xc(ixc)%id == XC_GGA_X_HJS_B88 ) then 
     call xc_func_set_ext_params_name(dft_xc(ixc)%func,C_CHAR_'_omega'//C_NULL_CHAR,dft_xc(ixc)%gamma)
     write(stdout,'(1x,a,f12.4)') 'Tuning range-separation in LIBXC with value: ',dft_xc(ixc)%gamma
   endif
   if( dft_xc(ixc)%id == XC_GGA_X_WPBEH ) then
     call xc_func_set_ext_params_name(dft_xc(ixc)%func,C_CHAR_'_omega'//C_NULL_CHAR,dft_xc(ixc)%gamma)
     write(stdout,'(1x,a,f12.4)') 'Tuning range-separation in LIBXC with value: ',dft_xc(ixc)%gamma
   endif
   if( dft_xc(ixc)%id == XC_GGA_X_ITYH_PBE ) then
     call xc_func_set_ext_params_name(dft_xc(ixc)%func,C_CHAR_'_omega'//C_NULL_CHAR,dft_xc(ixc)%gamma)
     write(stdout,'(1x,a,f12.4)') 'Tuning range-separation in LIBXC with value: ',dft_xc(ixc)%gamma
   endif
   if( dft_xc(ixc)%id == XC_GGA_X_ITYH ) then
     call xc_func_set_ext_params_name(dft_xc(ixc)%func,C_CHAR_'_omega'//C_NULL_CHAR,dft_xc(ixc)%gamma)
     write(stdout,'(1x,a,f12.4)') 'Tuning range-separation in LIBXC with value: ',dft_xc(ixc)%gamma
   endif

 enddo


 dft_xc(:)%needs_gradient = .FALSE.
 do ixc=1,dft_xc(1)%nxc
   dft_xc(ixc)%family = xc_family_from_id(dft_xc(ixc)%id,C_NULL_PTR,C_NULL_PTR)
   dft_xc(:)%needs_gradient = dft_xc(:)%needs_gradient .OR. dft_xc(ixc)%family == XC_FAMILY_GGA &
                                                       .OR. dft_xc(ixc)%family == XC_FAMILY_HYB_GGA
 enddo


! dft_xc(:)%needs_gradient = ANY( ( dft_xc(:)%family == XC_FAMILY_GGA     ) .AND. ( ABS(dft_xc(:)%coeff) > 1.0e-6_dp ) ) &
!                    .OR. ANY( ( dft_xc(:)%family == XC_FAMILY_HYB_GGA ) .AND. ( ABS(dft_xc(:)%coeff) > 1.0e-6_dp ) )

#endif


end subroutine init_libxc_info


!=========================================================================
subroutine copy_libxc_info(dft_xc_in,dft_xc_out)
 implicit none

 type(dft_xc_info),intent(in)                :: dft_xc_in(:)
 type(dft_xc_info),allocatable,intent(inout) :: dft_xc_out(:)
!=====
 integer :: nxc,ixc
!=====

 nxc = SIZE(dft_xc_in)
 allocate(dft_xc_out(nxc))

 do ixc=1,nxc
   dft_xc_out(ixc)%needs_gradient =  dft_xc_in(ixc)%needs_gradient
   dft_xc_out(ixc)%nxc            =  dft_xc_in(ixc)%nxc
   dft_xc_out(ixc)%id             =  dft_xc_in(ixc)%id
   dft_xc_out(ixc)%coeff          =  dft_xc_in(ixc)%coeff
   dft_xc_out(ixc)%gamma          =  dft_xc_in(ixc)%gamma
   dft_xc_out(ixc)%nspin          =  dft_xc_in(ixc)%nspin
   dft_xc_out(ixc)%family         =  dft_xc_in(ixc)%family
   dft_xc_out(ixc)%func           => dft_xc_in(ixc)%func
 enddo

end subroutine copy_libxc_info


!=========================================================================
subroutine destroy_libxc_info(dft_xc)
 implicit none

 type(dft_xc_info),allocatable,intent(inout) :: dft_xc(:)
!=====
 integer :: ixc
!=====

 do ixc=1,dft_xc(1)%nxc
#if !defined(NO_LIBXC)
   call xc_func_end(dft_xc(ixc)%func)
   call xc_func_free(dft_xc(ixc)%func)
#endif
 enddo
 deallocate(dft_xc)

end subroutine destroy_libxc_info


!=========================================================================
end module m_libxc_tools


!=========================================================================
