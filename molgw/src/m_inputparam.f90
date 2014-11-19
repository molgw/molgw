!=========================================================================
#include "macros.h"
!=========================================================================
module m_inputparam
 use m_definitions

 integer,protected            :: nspin,nscf
 real(dp),protected           :: alpha_mixing
 integer,protected            :: print_volume
 character(len=100),protected :: basis_name
 integer,protected            :: gaussian_type
 integer,protected            :: nangular_grid,nradial_grid
 real(dp),protected           :: electrons
 real(dp),protected           :: magnetization

 logical,protected            :: print_matrix
 logical,protected            :: print_basisset
 logical,protected            :: print_eri
 logical,protected            :: print_densitymatrix
 logical,protected            :: print_wfn
 logical,protected            :: print_specfunc

end module m_inputparam

