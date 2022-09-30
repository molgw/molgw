!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! basic tools to work with LIBCINT library
!=========================================================================
#include "molgw.h"
module m_libcint_tools
  use m_definitions
  use m_cart_to_pure
  use m_basis_set
  use m_atoms

  integer,private,parameter :: LMAX_LIBCINT = 7
  ! LIBCINT normalization is difficult to guess for pure gaussians
  ! So we obtain it from calls to the overlap rountine
  real(dp),private :: libcint_pure_norm(LMAX_LIBCINT)
!                  [ 1.0_dp ,                           & ! s
!                    1.0_dp ,                           & ! p
!                    1.0_dp / 1.092548430592079070_dp,  & ! d
!                    1.0_dp / 2.890611442640554055_dp,  & ! f
!                    1.0_dp / 8.671834327921662164_dp,  & ! g
!                    1.0_dp / 28.76122070983994253_dp,  & ! h
!                    1.0_dp / 103.8 ]   ! i

  ! LIBCINT_atm array meaning
  integer,private,parameter :: LIBCINT_CHARGE_OF  = 1
  integer,private,parameter :: LIBCINT_PTR_COORD  = 2
  integer,private,parameter :: LIBCINT_NUC_MOD_OF = 3
  integer,private,parameter :: LIBCINT_PTR_ZETA   = 4
  integer,private,parameter :: LIBCINT_ATM_SLOTS  = 6

  !
  ! LIBCINT_bas array meaning
  integer,private,parameter :: LIBCINT_ATOM_OF    = 1
  integer,private,parameter :: LIBCINT_ANG_OF     = 2
  integer,private,parameter :: LIBCINT_NPRIM_OF   = 3
  integer,private,parameter :: LIBCINT_NCTR_OF    = 4
  integer,private,parameter :: LIBCINT_KAPPA_OF   = 5
  integer,private,parameter :: LIBCINT_PTR_EXP    = 6
  integer,private,parameter :: LIBCINT_PTR_COEFF  = 7
  integer,private,parameter :: LIBCINT_BAS_SLOTS  = 8

  integer,private,parameter :: LIBCINT_PTR_COMMON_ORIG = 1
  integer,private,parameter :: LIBCINT_PTR_RINV_ORIG   = 4
  integer,private,parameter :: LIBCINT_PTR_RANGE_OMEGA = 8
  integer,private,parameter :: LIBCINT_PTR_ENV_START   = 20

  integer,private,parameter :: LIBCINT_env_size=100000

  logical,protected :: libcint_has_range_separation

  integer,external  :: cint2e_cart
  integer,external  :: cint2c2e_sph
  integer,external  :: cint2c2e_cart
  integer,external  :: cint3c2e_cart
  integer,external  :: cint3c2e_sph
  integer,external  :: cint3c1e_cart

  interface

    integer(C_INT) function cint1e_ovlp_cart(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_ovlp_cart

    integer(C_INT) function cint1e_ovlp_sph(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_ovlp_sph

    integer(C_INT) function cint1e_r2_origj_cart(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_r2_origj_cart

    integer(C_INT) function cint1e_r4_origj_cart(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_r4_origj_cart

    integer(C_INT) function cint1e_kin_cart(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_kin_cart

    integer(C_INT) function cint1e_nuc_cart(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_nuc_cart

    integer(C_INT) function cint1e_rinv_cart(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_rinv_cart

    integer(C_INT) function cint1e_ipovlp_cart(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_ipovlp_cart

    integer(C_INT) function cint1e_r_cart(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_r_cart

    integer(C_INT) function cint1e_rr_cart(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_rr_cart

    integer(C_INT) function cint1e_cg_irxp_cart(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_cg_irxp_cart

    integer(C_INT) function cint1e_giao_irjxp_cart(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_giao_irjxp_cart

    !FIXME this interface does not work and I don't know why
    !integer(C_INT) function cint2e_cart(array_cart, shls, atm, natm, bas, nbas, env, opt) bind(C)
    !  import :: C_INT,C_DOUBLE,C_PTR
    !  integer(C_INT),value  :: natm,nbas
    !  real(C_DOUBLE),intent(in) :: env(*)
    !  integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
    !  real(C_DOUBLE),intent(out) :: array_cart(*)
    !  type(C_PTR) :: opt(*)
    !end function cint2e_cart

  end interface


  interface transform_libcint_to_molgw
    module procedure transform_libcint_to_molgw_2d
    module procedure transform_libcint_to_molgw_3d
  end interface


contains


!=========================================================================
! Check with a known integral whether LIBCINT considers LIBCINT_PTR_RANGE_OMEGA
!
subroutine check_capability_libcint(lmax)
  implicit none

  integer,intent(inout) :: lmax
  !=====
  real(C_DOUBLE) :: fake_env(50),integral(1)
  integer(C_INT) :: fake_atm(LIBCINT_ATM_SLOTS,1)
  integer(C_INT) :: fake_bas(LIBCINT_BAS_SLOTS,1)
  integer(C_INT) :: shls(2)
  integer        :: info,off
  real(dp),parameter :: reference_value = 3.8052728203379578E-002_dp
  integer        :: il
  real(C_DOUBLE),allocatable :: ovlp(:)
  !=====

  ! LIBCINT goes up to lmax=7
  if( lmax > 0 ) then
    lmax = MIN(lmax,LMAX_LIBCINT)
  else
    lmax = LMAX_LIBCINT
  endif

  !
  ! Check if LIBCINT has ERF range-separation
  fake_env(LIBCINT_PTR_RANGE_OMEGA+1) = 0.11_dp  ! omega = 0.11 bohr**-1
  off = LIBCINT_PTR_ENV_START

  fake_atm(LIBCINT_CHARGE_OF,1) = 1
  fake_atm(LIBCINT_PTR_COORD,1) = off
  fake_env(off+1:off+3) = 0.0_dp
  off = off + 3
  fake_bas(LIBCINT_ATOM_OF  ,1)  = 0 ! C convention starts with 0
  fake_bas(LIBCINT_ANG_OF   ,1)  = 0
  fake_bas(LIBCINT_NPRIM_OF ,1)  = 1
  fake_bas(LIBCINT_NCTR_OF  ,1)  = 1
  fake_bas(LIBCINT_PTR_EXP  ,1)  = off ! note the 0-based index
  fake_env(off+1) = 2.0_dp  ! alpha = 2.0
  off = off + 1
  fake_bas(LIBCINT_PTR_COEFF,1) = off
  fake_env(off+1) = 1.0    ! coefficient = 1.0


  shls(:) = 0
#if defined(HAVE_LIBCINT)
  info = cint2c2e_cart(integral, shls, fake_atm, 1_C_INT, fake_bas, 1_C_INT, fake_env, 0_C_LONG)
#endif

  libcint_has_range_separation = ABS( integral(1) - reference_value ) < 1.0e-12_dp


  !
  ! Find normalization for LIBCINT spherical
  do il=0,lmax
    fake_env(LIBCINT_PTR_RANGE_OMEGA+1) = 0.11_dp  ! omega = 0.11 bohr**-1
    off = LIBCINT_PTR_ENV_START

    fake_atm(LIBCINT_CHARGE_OF,1) = 1
    fake_atm(LIBCINT_PTR_COORD,1) = off
    fake_env(off+1:off+3) = 0.0_dp
    off = off + 3
    fake_bas(LIBCINT_ATOM_OF  ,1)  = 0 ! C convention starts with 0
    fake_bas(LIBCINT_ANG_OF   ,1)  = 0
    fake_bas(LIBCINT_NPRIM_OF ,1)  = 1
    fake_bas(LIBCINT_NCTR_OF  ,1)  = 1
    fake_bas(LIBCINT_PTR_EXP  ,1)  = off ! note the 0-based index
    fake_env(off+1) = 1.0_dp  ! alpha = 1.0
    off = off + 1
    fake_bas(LIBCINT_PTR_COEFF,1) = off
    select case(il)
    case(0,1)
      fake_env(off+1) = SQRT(4.0_dp * pi) / SQRT( 2.0_dp * il + 1 )  &
                       * ( 2.0_dp / pi )**0.75_dp * 2.0_dp**il
    case default
      fake_env(off+1) = ( 2.0_dp / pi )**0.75_dp * 2.0_dp**il
    end select

    allocate(ovlp((2*il+1)**2))
    fake_bas(LIBCINT_ANG_OF   ,1)  = il
#if defined(HAVE_LIBCINT)
    info = cint1e_ovlp_sph(ovlp, shls, fake_atm, 1_C_INT, fake_bas, 1_C_INT, fake_env)
#endif
    libcint_pure_norm(il) = 1.0_dp / SQRT(ovlp(1))
    write(stdout,*) il,ovlp(1),libcint_pure_norm(il)

    deallocate(ovlp)
  enddo


end subroutine check_capability_libcint


!=========================================================================
subroutine set_rinv_origin_libcint(x0,env_local)
  implicit none
  real(dp),intent(in)          :: x0(3)
  real(C_DOUBLE),intent(inout) :: env_local(:)
  !=====
  !=====

  env_local(LIBCINT_PTR_RINV_ORIG+1:LIBCINT_PTR_RINV_ORIG+3) = x0(:)

end subroutine set_rinv_origin_libcint


!=========================================================================
subroutine set_erf_screening_length_libcint(basis,rcut)
  implicit none
  type(basis_set),intent(inout) :: basis
  real(dp),intent(in)           :: rcut
  !=====
  !=====

  ! LIBCINT convention: a positive omega = 1/rcut means long-range only
  if( rcut > 1.0e-12_dp ) then
    ! Ensure that LIBCINT can calculated "erf" range-separated Coulomb interaction
    if( .NOT. libcint_has_range_separation ) then
       call die('set_erf_screening_length_libcint: LIBCINT compilation does not support range separation')
    endif
    basis%LIBCINT_env(LIBCINT_PTR_RANGE_OMEGA+1) = 1.0_C_DOUBLE / rcut
  else
    basis%LIBCINT_env(LIBCINT_PTR_RANGE_OMEGA+1) = 0.0_C_DOUBLE
  endif

end subroutine set_erf_screening_length_libcint


!=========================================================================
subroutine init_libcint(basis1,basis2)
  implicit none
  type(basis_set),intent(inout)       :: basis1
  type(basis_set),optional,intent(in) :: basis2
  !=====
  integer :: off,icenter_basis,ishell
  !=====

  if( PRESENT(basis2) ) then
    write(stdout,'(/,1x,a)') 'Initialize LIBCINT internal data for 2 basis sets (basis and auxiliary basis)'
  else
    write(stdout,'(/,1x,a)') 'Initialize LIBCINT internal data for 1 basis set (basis or auxiliary basis)'
  endif

  basis1%LIBCINT_natm = ncenter_basis
  basis1%LIBCINT_nbas = basis1%nshell
  if( PRESENT(basis2) ) then
    basis1%LIBCINT_nbas = basis1%LIBCINT_nbas + basis2%nshell
  endif

  allocate(basis1%LIBCINT_atm(LIBCINT_ATM_SLOTS,basis1%LIBCINT_natm))
  allocate(basis1%LIBCINT_bas(LIBCINT_BAS_SLOTS,basis1%LIBCINT_nbas))
  allocate(basis1%LIBCINT_env(LIBCINT_env_size))
  ! Zero-ing the env(:) array is compulsory
  ! Fortan may put garbage (ifort) and C does not like it
  basis1%LIBCINT_env(:) = 0.0_C_DOUBLE

  ! real space origin
  basis1%LIBCINT_env(LIBCINT_PTR_COMMON_ORIG+1) = 0.0_C_DOUBLE
  basis1%LIBCINT_env(LIBCINT_PTR_COMMON_ORIG+2) = 0.0_C_DOUBLE
  basis1%LIBCINT_env(LIBCINT_PTR_COMMON_ORIG+3) = 0.0_C_DOUBLE
  off = LIBCINT_PTR_ENV_START

  do icenter_basis=1,ncenter_basis
    basis1%LIBCINT_atm(LIBCINT_CHARGE_OF,icenter_basis) = zbasis(icenter_basis)
    basis1%LIBCINT_atm(LIBCINT_PTR_COORD,icenter_basis) = off ! note the 0-based index
    basis1%LIBCINT_env(off+1:off+3) = xbasis(:,icenter_basis)
    off = off + 3
  enddo

  do ishell=1,basis1%nshell
    basis1%LIBCINT_bas(LIBCINT_ATOM_OF  ,ishell)  = basis1%shell(ishell)%icenter - 1 ! C convention starts with 0
    basis1%LIBCINT_bas(LIBCINT_ANG_OF   ,ishell)  = basis1%shell(ishell)%am
    basis1%LIBCINT_bas(LIBCINT_NPRIM_OF ,ishell)  = basis1%shell(ishell)%ng
    basis1%LIBCINT_bas(LIBCINT_NCTR_OF  ,ishell)  = 1   ! so far take shells one at a time (even if they share exponents)
    basis1%LIBCINT_bas(LIBCINT_PTR_EXP  ,ishell)  = off ! note the 0-based index
    basis1%LIBCINT_env(off+1:off+basis1%shell(ishell)%ng) = basis1%shell(ishell)%alpha(:)
    off = off + basis1%shell(ishell)%ng
    basis1%LIBCINT_bas(LIBCINT_PTR_COEFF,ishell) = off
    select case(basis1%shell(ishell)%am)
    case(0,1)
      basis1%LIBCINT_env(off+1:off+basis1%shell(ishell)%ng) = basis1%shell(ishell)%coeff(:) &
                * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * basis1%shell(ishell)%am + 1 )
    case default
      basis1%LIBCINT_env(off+1:off+basis1%shell(ishell)%ng) = basis1%shell(ishell)%coeff(:)
    end select
    off = off + basis1%shell(ishell)%ng
  enddo

  basis1%LIBCINT_offset = basis1%nshell

  if( PRESENT(basis2) ) then
    do ishell=1,basis2%nshell
      basis1%LIBCINT_bas(LIBCINT_ATOM_OF  ,basis1%LIBCINT_offset + ishell)  = basis2%shell(ishell)%icenter - 1 ! C convention starts with 0
      basis1%LIBCINT_bas(LIBCINT_ANG_OF   ,basis1%LIBCINT_offset + ishell)  = basis2%shell(ishell)%am
      basis1%LIBCINT_bas(LIBCINT_NPRIM_OF ,basis1%LIBCINT_offset + ishell)  = basis2%shell(ishell)%ng
      basis1%LIBCINT_bas(LIBCINT_NCTR_OF  ,basis1%LIBCINT_offset + ishell)  = 1   ! so far take shells one at a time (even if they share exponents)
      basis1%LIBCINT_bas(LIBCINT_PTR_EXP  ,basis1%LIBCINT_offset + ishell)  = off ! note the 0-based index
      basis1%LIBCINT_env(off+1:off+basis2%shell(ishell)%ng) = basis2%shell(ishell)%alpha(:)
      off = off + basis2%shell(ishell)%ng
      basis1%LIBCINT_bas(LIBCINT_PTR_COEFF,basis1%LIBCINT_offset + ishell) = off
      select case(basis2%shell(ishell)%am)
      case(0,1)
        basis1%LIBCINT_env(off+1:off+basis2%shell(ishell)%ng) = basis2%shell(ishell)%coeff(:) &
                  * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * basis2%shell(ishell)%am + 1 )
      case default
        basis1%LIBCINT_env(off+1:off+basis2%shell(ishell)%ng) = basis2%shell(ishell)%coeff(:)
      end select
      off = off + basis2%shell(ishell)%ng
    enddo
  endif
  write(stdout,'(1x,a,i8)') 'LIBCINT environment maximum index: ',off

end subroutine init_libcint


!=========================================================================
subroutine destroy_libcint(basis)
  implicit none
  type(basis_set),intent(inout) :: basis
  !=====
  !=====

  basis%LIBCINT_natm = 0
  basis%LIBCINT_nbas = 0
  basis%LIBCINT_offset = 0
  if( ALLOCATED(basis%LIBCINT_atm) ) deallocate(basis%LIBCINT_atm)
  if( ALLOCATED(basis%LIBCINT_bas) ) deallocate(basis%LIBCINT_bas)
  if( ALLOCATED(basis%LIBCINT_env) ) deallocate(basis%LIBCINT_env)

end subroutine destroy_libcint


!=========================================================================
subroutine libcint_3center(amA,contrdepthA,A,alphaA,cA, &
                           amC,contrdepthC,C,alphaC,cC, &
                           amD,contrdepthD,D,alphaD,cD, &
                           rcut,eriACD)
  implicit none

  integer(C_INT),intent(in)    :: amA,contrdepthA
  real(C_DOUBLE),intent(in)    :: A(:)
  real(C_DOUBLE),intent(in)    :: alphaA(:)
  real(C_DOUBLE),intent(in)    :: cA(:)
  integer(C_INT),intent(in)    :: amC,contrdepthC
  real(C_DOUBLE),intent(in)    :: C(:)
  real(C_DOUBLE),intent(in)    :: alphaC(:)
  real(C_DOUBLE),intent(in)    :: cC(:)
  integer(C_INT),intent(in)    :: amD,contrdepthD
  real(C_DOUBLE),intent(in)    :: D(:)
  real(C_DOUBLE),intent(in)    :: alphaD(:)
  real(C_DOUBLE),intent(in)    :: cD(:)
  real(C_DOUBLE),intent(in)    :: rcut
  real(C_DOUBLE),intent(inout) :: eriACD(:)
  !=====
  real(C_DOUBLE) :: tmp_env(1000)
  integer(C_INT) :: tmp_atm(LIBCINT_ATM_SLOTS,3)
  integer(C_INT) :: tmp_bas(LIBCINT_BAS_SLOTS,3)
  integer(C_INT) :: shls(3)
  integer        :: info,off
  !=====

  tmp_env(:) = 0.0_dp

  if( rcut < 1.0e-12_dp ) then
    tmp_env(LIBCINT_PTR_RANGE_OMEGA+1) = 0.0_C_DOUBLE
  else
    tmp_env(LIBCINT_PTR_RANGE_OMEGA+1) = 1.0_C_DOUBLE / rcut
  endif
  off = LIBCINT_PTR_ENV_START

  tmp_atm(LIBCINT_CHARGE_OF,1) = 1
  tmp_atm(LIBCINT_PTR_COORD,1) = off
  tmp_env(off+1:off+3) = D(1:3)
  off = off + 3
  tmp_atm(LIBCINT_CHARGE_OF,2) = 1
  tmp_atm(LIBCINT_PTR_COORD,2) = off
  tmp_env(off+1:off+3) = C(1:3)
  off = off + 3
  tmp_atm(LIBCINT_CHARGE_OF,3) = 1
  tmp_atm(LIBCINT_PTR_COORD,3) = off
  tmp_env(off+1:off+3) = A(1:3)
  off = off + 3

  !
  ! D
  tmp_bas(LIBCINT_ATOM_OF  ,1)  = 0 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,1)  = amD
  tmp_bas(LIBCINT_NPRIM_OF ,1)  = contrdepthD
  tmp_bas(LIBCINT_NCTR_OF  ,1)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,1)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthD) = alphaD(1:contrdepthD)
  off = off + contrdepthD
  tmp_bas(LIBCINT_PTR_COEFF,1) = off
  select case(amD)
  case(0,1)
    tmp_env(off+1:off+contrdepthD) = cD(1:contrdepthD) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amD + 1 )
  case default
    tmp_env(off+1:off+contrdepthD) = cD(1:contrdepthD)
  end select
  off = off + contrdepthD
  !
  ! C
  tmp_bas(LIBCINT_ATOM_OF  ,2)  = 1 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,2)  = amC
  tmp_bas(LIBCINT_NPRIM_OF ,2)  = contrdepthC
  tmp_bas(LIBCINT_NCTR_OF  ,2)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,2)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthC) = alphaC(1:contrdepthC)
  off = off + contrdepthC
  tmp_bas(LIBCINT_PTR_COEFF,2) = off
  select case(amC)
  case(0,1)
    tmp_env(off+1:off+contrdepthC) = cC(1:contrdepthC) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amC + 1 )
  case default
    tmp_env(off+1:off+contrdepthC) = cC(1:contrdepthC)
  end select
  off = off + contrdepthC
  !
  ! A
  tmp_bas(LIBCINT_ATOM_OF  ,3)  = 2 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,3)  = amA
  tmp_bas(LIBCINT_NPRIM_OF ,3)  = contrdepthA
  tmp_bas(LIBCINT_NCTR_OF  ,3)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,3)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthA) = alphaA(1:contrdepthA)
  off = off + contrdepthA
  tmp_bas(LIBCINT_PTR_COEFF,3) = off
  select case(amA)
  case(0,1)
    tmp_env(off+1:off+contrdepthA) = cA(1:contrdepthA) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amA + 1 )
  case default
    tmp_env(off+1:off+contrdepthA) = cA(1:contrdepthA)
  end select
  off = off + contrdepthA


  shls(1) = 0
  shls(2) = 1
  shls(3) = 2

#if defined(HAVE_LIBCINT)
  info = cint3c2e_cart(eriACD, shls, tmp_atm, 3_C_INT, tmp_bas, 3_C_INT, tmp_env, 0_C_LONG)
#endif


end subroutine libcint_3center


!=========================================================================
subroutine libcint_overlap(amA,contrdepthA,A,alphaA,cA, &
                           amC,contrdepthC,C,alphaC,cC, &
                           ovlpAC)
  implicit none

  integer(C_INT),intent(in)    :: amA,contrdepthA
  real(C_DOUBLE),intent(in)    :: A(:)
  real(C_DOUBLE),intent(in)    :: alphaA(:)
  real(C_DOUBLE),intent(in)    :: cA(:)
  integer(C_INT),intent(in)    :: amC,contrdepthC
  real(C_DOUBLE),intent(in)    :: C(:)
  real(C_DOUBLE),intent(in)    :: alphaC(:)
  real(C_DOUBLE),intent(in)    :: cC(:)
  real(C_DOUBLE),intent(inout)    :: ovlpAC(:)
  !=====
  real(C_DOUBLE) :: tmp_env(1000)
  integer(C_INT) :: tmp_atm(LIBCINT_ATM_SLOTS,2)
  integer(C_INT) :: tmp_bas(LIBCINT_BAS_SLOTS,2)
  integer(C_INT) :: shls(2)
  integer        :: info,off
  !=====

  tmp_env(:) = 0.0_dp

  off = LIBCINT_PTR_ENV_START

  tmp_atm(LIBCINT_CHARGE_OF,1) = 1
  tmp_atm(LIBCINT_PTR_COORD,1) = off
  tmp_env(off+1:off+3) = C(1:3)
  off = off + 3
  tmp_atm(LIBCINT_CHARGE_OF,2) = 1
  tmp_atm(LIBCINT_PTR_COORD,2) = off
  tmp_env(off+1:off+3) = A(1:3)
  off = off + 3

  !
  ! C
  tmp_bas(LIBCINT_ATOM_OF  ,1)  = 0 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,1)  = amC
  tmp_bas(LIBCINT_NPRIM_OF ,1)  = contrdepthC
  tmp_bas(LIBCINT_NCTR_OF  ,1)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,1)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthC) = alphaC(1:contrdepthC)
  off = off + contrdepthC
  tmp_bas(LIBCINT_PTR_COEFF,1) = off
  select case(amC)
  case(0,1)
    tmp_env(off+1:off+contrdepthC) = cC(1:contrdepthC) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amC + 1 )
  case default
    tmp_env(off+1:off+contrdepthC) = cC(1:contrdepthC)
  end select
  off = off + contrdepthC
  !
  ! A
  tmp_bas(LIBCINT_ATOM_OF  ,2)  = 1 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,2)  = amA
  tmp_bas(LIBCINT_NPRIM_OF ,2)  = contrdepthA
  tmp_bas(LIBCINT_NCTR_OF  ,2)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,2)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthA) = alphaA(1:contrdepthA)
  off = off + contrdepthA
  tmp_bas(LIBCINT_PTR_COEFF,2) = off
  select case(amA)
  case(0,1)
    tmp_env(off+1:off+contrdepthA) = cA(1:contrdepthA) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amA + 1 )
  case default
    tmp_env(off+1:off+contrdepthA) = cA(1:contrdepthA)
  end select
  off = off + contrdepthA


  shls(1) = 0
  shls(2) = 1

#if defined(HAVE_LIBCINT)
  info = cint1e_ovlp_cart(ovlpAC, shls, tmp_atm, 2_C_INT, tmp_bas, 2_C_INT, tmp_env)
#endif


end subroutine libcint_overlap


!=========================================================================
subroutine libcint_kinetic(amA,contrdepthA,A,alphaA,cA, &
                           amC,contrdepthC,C,alphaC,cC, &
                           kinAC)
  implicit none

  integer(C_INT),intent(in)    :: amA,contrdepthA
  real(C_DOUBLE),intent(in)    :: A(:)
  real(C_DOUBLE),intent(in)    :: alphaA(:)
  real(C_DOUBLE),intent(in)    :: cA(:)
  integer(C_INT),intent(in)    :: amC,contrdepthC
  real(C_DOUBLE),intent(in)    :: C(:)
  real(C_DOUBLE),intent(in)    :: alphaC(:)
  real(C_DOUBLE),intent(in)    :: cC(:)
  real(C_DOUBLE),intent(inout) :: kinAC(:)
  !=====
  real(C_DOUBLE) :: tmp_env(1000)
  integer(C_INT) :: tmp_atm(LIBCINT_ATM_SLOTS,2)
  integer(C_INT) :: tmp_bas(LIBCINT_BAS_SLOTS,2)
  integer(C_INT) :: shls(2)
  integer        :: info,off
  !=====

  tmp_env(:) = 0.0_dp

  off = LIBCINT_PTR_ENV_START

  tmp_atm(LIBCINT_CHARGE_OF,1) = 1
  tmp_atm(LIBCINT_PTR_COORD,1) = off
  tmp_env(off+1:off+3) = C(1:3)
  off = off + 3
  tmp_atm(LIBCINT_CHARGE_OF,2) = 1
  tmp_atm(LIBCINT_PTR_COORD,2) = off
  tmp_env(off+1:off+3) = A(1:3)
  off = off + 3

  !
  ! C
  tmp_bas(LIBCINT_ATOM_OF  ,1)  = 0 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,1)  = amC
  tmp_bas(LIBCINT_NPRIM_OF ,1)  = contrdepthC
  tmp_bas(LIBCINT_NCTR_OF  ,1)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,1)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthC) = alphaC(1:contrdepthC)
  off = off + contrdepthC
  tmp_bas(LIBCINT_PTR_COEFF,1) = off
  select case(amC)
  case(0,1)
    tmp_env(off+1:off+contrdepthC) = cC(1:contrdepthC) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amC + 1 )
  case default
    tmp_env(off+1:off+contrdepthC) = cC(1:contrdepthC)
  end select
  off = off + contrdepthC
  !
  ! A
  tmp_bas(LIBCINT_ATOM_OF  ,2)  = 1 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,2)  = amA
  tmp_bas(LIBCINT_NPRIM_OF ,2)  = contrdepthA
  tmp_bas(LIBCINT_NCTR_OF  ,2)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,2)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthA) = alphaA(1:contrdepthA)
  off = off + contrdepthA
  tmp_bas(LIBCINT_PTR_COEFF,2) = off
  select case(amA)
  case(0,1)
    tmp_env(off+1:off+contrdepthA) = cA(1:contrdepthA) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amA + 1 )
  case default
    tmp_env(off+1:off+contrdepthA) = cA(1:contrdepthA)
  end select
  off = off + contrdepthA


  shls(1) = 0
  shls(2) = 1

#if defined(HAVE_LIBCINT)
  info = cint1e_kin_cart(kinAC, shls, tmp_atm, 2_C_INT, tmp_bas, 2_C_INT, tmp_env)
#endif


end subroutine libcint_kinetic


!=========================================================================
subroutine libcint_elecpot(amA,contrdepthA,A,alphaA,cA, &
                           amC,contrdepthC,C,alphaC,cC, &
                           D,elecpotAC)
  implicit none

  integer(C_INT),intent(in)    :: amA,contrdepthA
  real(C_DOUBLE),intent(in)    :: A(:)
  real(C_DOUBLE),intent(in)    :: alphaA(:)
  real(C_DOUBLE),intent(in)    :: cA(:)
  integer(C_INT),intent(in)    :: amC,contrdepthC
  real(C_DOUBLE),intent(in)    :: C(:)
  real(C_DOUBLE),intent(in)    :: alphaC(:)
  real(C_DOUBLE),intent(in)    :: cC(:)
  real(C_DOUBLE),intent(in)    :: D(:)
  real(C_DOUBLE),intent(inout) :: elecpotAC(:)
  !=====
  real(C_DOUBLE) :: tmp_env(1000)
  integer(C_INT) :: tmp_atm(LIBCINT_ATM_SLOTS,2)
  integer(C_INT) :: tmp_bas(LIBCINT_BAS_SLOTS,2)
  integer(C_INT) :: shls(2)
  integer        :: info,off
  !=====

  tmp_env(:) = 0.0_dp
  tmp_env(LIBCINT_PTR_RINV_ORIG+1:LIBCINT_PTR_RINV_ORIG+3) = D(:)

  off = LIBCINT_PTR_ENV_START

  tmp_atm(LIBCINT_CHARGE_OF,1) = 1
  tmp_atm(LIBCINT_PTR_COORD,1) = off
  tmp_env(off+1:off+3) = C(1:3)
  off = off + 3
  tmp_atm(LIBCINT_CHARGE_OF,2) = 1
  tmp_atm(LIBCINT_PTR_COORD,2) = off
  tmp_env(off+1:off+3) = A(1:3)
  off = off + 3

  !
  ! C
  tmp_bas(LIBCINT_ATOM_OF  ,1)  = 0 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,1)  = amC
  tmp_bas(LIBCINT_NPRIM_OF ,1)  = contrdepthC
  tmp_bas(LIBCINT_NCTR_OF  ,1)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,1)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthC) = alphaC(1:contrdepthC)
  off = off + contrdepthC
  tmp_bas(LIBCINT_PTR_COEFF,1) = off
  select case(amC)
  case(0,1)
    tmp_env(off+1:off+contrdepthC) = cC(1:contrdepthC) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amC + 1 )
  case default
    tmp_env(off+1:off+contrdepthC) = cC(1:contrdepthC)
  end select
  off = off + contrdepthC
  !
  ! A
  tmp_bas(LIBCINT_ATOM_OF  ,2)  = 1 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,2)  = amA
  tmp_bas(LIBCINT_NPRIM_OF ,2)  = contrdepthA
  tmp_bas(LIBCINT_NCTR_OF  ,2)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,2)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthA) = alphaA(1:contrdepthA)
  off = off + contrdepthA
  tmp_bas(LIBCINT_PTR_COEFF,2) = off
  select case(amA)
  case(0,1)
    tmp_env(off+1:off+contrdepthA) = cA(1:contrdepthA) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amA + 1 )
  case default
    tmp_env(off+1:off+contrdepthA) = cA(1:contrdepthA)
  end select
  off = off + contrdepthA


  shls(1) = 0
  shls(2) = 1

#if defined(HAVE_LIBCINT)
  info = cint1e_rinv_cart(elecpotAC, shls, tmp_atm, 2_C_INT, tmp_bas, 2_C_INT, tmp_env)
#endif


end subroutine libcint_elecpot


!=========================================================================
subroutine libcint_overlap_3center(amA,contrdepthA,A,alphaA,cA, &
                           amC,contrdepthC,C,alphaC,cC, &
                           amD,contrdepthD,D,alphaD,cD, &
                           ovlpACD)
  implicit none

  integer(C_INT),intent(in)    :: amA,contrdepthA
  real(C_DOUBLE),intent(in)    :: A(:)
  real(C_DOUBLE),intent(in)    :: alphaA(:)
  real(C_DOUBLE),intent(in)    :: cA(:)
  integer(C_INT),intent(in)    :: amC,contrdepthC
  real(C_DOUBLE),intent(in)    :: C(:)
  real(C_DOUBLE),intent(in)    :: alphaC(:)
  real(C_DOUBLE),intent(in)    :: cC(:)
  integer(C_INT),intent(in)    :: amD,contrdepthD
  real(C_DOUBLE),intent(in)    :: D(:)
  real(C_DOUBLE),intent(in)    :: alphaD(:)
  real(C_DOUBLE),intent(in)    :: cD(:)
  real(C_DOUBLE),intent(inout) :: ovlpACD(:)
  !=====
  real(C_DOUBLE) :: tmp_env(1000)
  integer(C_INT) :: tmp_atm(LIBCINT_ATM_SLOTS,3)
  integer(C_INT) :: tmp_bas(LIBCINT_BAS_SLOTS,3)
  integer(C_INT) :: shls(3)
  integer        :: info,off
  !=====

  tmp_env(:) = 0.0_dp

  off = LIBCINT_PTR_ENV_START

  tmp_atm(LIBCINT_CHARGE_OF,1) = 1
  tmp_atm(LIBCINT_PTR_COORD,1) = off
  tmp_env(off+1:off+3) = D(1:3)
  off = off + 3
  tmp_atm(LIBCINT_CHARGE_OF,2) = 1
  tmp_atm(LIBCINT_PTR_COORD,2) = off
  tmp_env(off+1:off+3) = C(1:3)
  off = off + 3
  tmp_atm(LIBCINT_CHARGE_OF,3) = 1
  tmp_atm(LIBCINT_PTR_COORD,3) = off
  tmp_env(off+1:off+3) = A(1:3)
  off = off + 3

  !
  ! D
  tmp_bas(LIBCINT_ATOM_OF  ,1)  = 0 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,1)  = amD
  tmp_bas(LIBCINT_NPRIM_OF ,1)  = contrdepthD
  tmp_bas(LIBCINT_NCTR_OF  ,1)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,1)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthD) = alphaD(1:contrdepthD)
  off = off + contrdepthD
  tmp_bas(LIBCINT_PTR_COEFF,1) = off
  select case(amD)
  case(0,1)
    tmp_env(off+1:off+contrdepthD) = cD(1:contrdepthD) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amD + 1 )
  case default
    tmp_env(off+1:off+contrdepthD) = cD(1:contrdepthD)
  end select
  off = off + contrdepthD
  !
  ! C
  tmp_bas(LIBCINT_ATOM_OF  ,2)  = 1 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,2)  = amC
  tmp_bas(LIBCINT_NPRIM_OF ,2)  = contrdepthC
  tmp_bas(LIBCINT_NCTR_OF  ,2)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,2)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthC) = alphaC(1:contrdepthC)
  off = off + contrdepthC
  tmp_bas(LIBCINT_PTR_COEFF,2) = off
  select case(amC)
  case(0,1)
    tmp_env(off+1:off+contrdepthC) = cC(1:contrdepthC) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amC + 1 )
  case default
    tmp_env(off+1:off+contrdepthC) = cC(1:contrdepthC)
  end select
  off = off + contrdepthC
  !
  ! A
  tmp_bas(LIBCINT_ATOM_OF  ,3)  = 2 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,3)  = amA
  tmp_bas(LIBCINT_NPRIM_OF ,3)  = contrdepthA
  tmp_bas(LIBCINT_NCTR_OF  ,3)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,3)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthA) = alphaA(1:contrdepthA)
  off = off + contrdepthA
  tmp_bas(LIBCINT_PTR_COEFF,3) = off
  select case(amA)
  case(0,1)
    tmp_env(off+1:off+contrdepthA) = cA(1:contrdepthA) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amA + 1 )
  case default
    tmp_env(off+1:off+contrdepthA) = cA(1:contrdepthA)
  end select
  off = off + contrdepthA


  shls(1) = 0
  shls(2) = 1
  shls(3) = 2

#if defined(HAVE_LIBCINT)
  info = cint3c1e_cart(ovlpACD, shls, tmp_atm, 3_C_INT, tmp_bas, 3_C_INT, tmp_env, 0_C_LONG)
#endif


end subroutine libcint_overlap_3center


!=========================================================================
! Special routine for gth projectors
! p_i^lm (r) = N_i^l Y_lm( (r-RC)^) (r-RC)**(l+2i−2) exp(-alphaC * (r-RC)**2)
! Strategy
! if iproj=1, it reduces to an overlap
! if iproj=2, it reduces to an int1e_r2_origj
! if iproj=3, it reduces to an int1e_r4_origj
subroutine libcint_gth_projector(amA,contrdepthA,A,alphaA,cA, &
                                 amC,contrdepthC,C,alphaC,cC, &
                                 iproj,ovlpAC)
  implicit none

  integer(C_INT),intent(in)    :: amA,contrdepthA
  real(C_DOUBLE),intent(in)    :: A(:)
  real(C_DOUBLE),intent(in)    :: alphaA(:)
  real(C_DOUBLE),intent(in)    :: cA(:)
  integer(C_INT),intent(in)    :: amC,contrdepthC
  real(C_DOUBLE),intent(in)    :: C(:)
  real(C_DOUBLE),intent(in)    :: alphaC(:)
  real(C_DOUBLE),intent(in)    :: cC(:)
  integer(C_INT),intent(in)    :: iproj
  real(C_DOUBLE),intent(inout)    :: ovlpAC(:)
  !=====
  real(C_DOUBLE) :: tmp_env(1000)
  integer(C_INT) :: tmp_atm(LIBCINT_ATM_SLOTS,2)
  integer(C_INT) :: tmp_bas(LIBCINT_BAS_SLOTS,2)
  integer(C_INT) :: shls(2)
  integer        :: info,off
  !=====

  tmp_env(:) = 0.0_dp

  off = LIBCINT_PTR_ENV_START

  tmp_atm(LIBCINT_CHARGE_OF,1) = 1
  tmp_atm(LIBCINT_PTR_COORD,1) = off
  tmp_env(off+1:off+3) = C(1:3)
  off = off + 3
  tmp_atm(LIBCINT_CHARGE_OF,2) = 1
  tmp_atm(LIBCINT_PTR_COORD,2) = off
  tmp_env(off+1:off+3) = A(1:3)
  off = off + 3

  !
  ! C
  tmp_bas(LIBCINT_ATOM_OF  ,1)  = 0 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,1)  = amC
  tmp_bas(LIBCINT_NPRIM_OF ,1)  = contrdepthC
  tmp_bas(LIBCINT_NCTR_OF  ,1)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,1)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthC) = alphaC(1:contrdepthC)
  off = off + contrdepthC
  tmp_bas(LIBCINT_PTR_COEFF,1) = off
  tmp_env(off+1:off+contrdepthC) = cC(1:contrdepthC)
  off = off + contrdepthC
  !
  ! A
  tmp_bas(LIBCINT_ATOM_OF  ,2)  = 1 ! C convention starts with 0
  tmp_bas(LIBCINT_ANG_OF   ,2)  = amA
  tmp_bas(LIBCINT_NPRIM_OF ,2)  = contrdepthA
  tmp_bas(LIBCINT_NCTR_OF  ,2)  = 1
  tmp_bas(LIBCINT_PTR_EXP  ,2)  = off ! note the 0-based index
  tmp_env(off+1:off+contrdepthA) = alphaA(1:contrdepthA)
  off = off + contrdepthA
  tmp_bas(LIBCINT_PTR_COEFF,2) = off
  select case(amA)
  case(0,1)
    tmp_env(off+1:off+contrdepthA) = cA(1:contrdepthA) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * amA + 1 )
  case default
    tmp_env(off+1:off+contrdepthA) = cA(1:contrdepthA)
  end select
  off = off + contrdepthA


  shls(1) = 0
  shls(2) = 1

#if defined(HAVE_LIBCINT)
  select case(iproj)
  case(1)
    info = cint1e_ovlp_cart(ovlpAC, shls, tmp_atm, 2_C_INT, tmp_bas, 2_C_INT, tmp_env)
  case(2)
    info = cint1e_r2_origj_cart(ovlpAC, shls, tmp_atm, 2_C_INT, tmp_bas, 2_C_INT, tmp_env)
  case(3)
    info = cint1e_r4_origj_cart(ovlpAC, shls, tmp_atm, 2_C_INT, tmp_bas, 2_C_INT, tmp_env)
  case default
  end select
#endif


end subroutine libcint_gth_projector


!=========================================================================
subroutine transform_libcint_to_molgw_2d(gaussian_type,am1,am2,array_in,matrix_out)
  implicit none
  character(len=4),intent(in)      :: gaussian_type
  integer,intent(in)               :: am1,am2
  real(C_DOUBLE),intent(in)        :: array_in(:)
  real(dp),allocatable,intent(out) :: matrix_out(:,:)
  !=====
  integer :: n1,n2,n1c,n2c
  integer :: ii,i1,i2
  integer :: gt_tag
  !=====

  gt_tag = get_gaussian_type_tag(gaussian_type)
  n1c = number_basis_function_am('CART',am1)
  n2c = number_basis_function_am('CART',am2)
  n1  = number_basis_function_am(gaussian_type,am1)
  n2  = number_basis_function_am(gaussian_type,am2)

  if( .NOT. ALLOCATED(matrix_out) ) allocate(matrix_out(n1,n2))

  if( gt_tag == CARTG ) then
    ii = 0
    do i1=1,n1
      do i2=1,n2
        ii = ii + 1
        matrix_out(i1,i2) = array_in(ii) * cart_to_pure_norm(am1,CARTG)%matrix(i1,i1) &
                                         * cart_to_pure_norm(am2,CARTG)%matrix(i2,i2)
      enddo
    enddo
  else
    ii = 0
    do i1=1,n1
      do i2=1,n2
        ii = ii + 1
        matrix_out(i1,i2) = array_in(ii) * libcint_pure_norm(am1) * libcint_pure_norm(am2)
      enddo
    enddo
  endif


end subroutine transform_libcint_to_molgw_2d


!=========================================================================
subroutine transform_libcint_to_molgw_3d(gaussian_type_left,am1,gaussian_type_right,am2,am3,array_in,matrix_out)
  implicit none
  character(len=4),intent(in)      :: gaussian_type_left,gaussian_type_right
  integer,intent(in)               :: am1,am2,am3
  real(C_DOUBLE),intent(in)        :: array_in(:)
  real(dp),allocatable,intent(out) :: matrix_out(:,:,:)
  !=====
  integer :: n1,n2,n3,n1c,n2c,n3c
  integer :: i1,i2,i3,ii
  integer :: gt_tagl,gt_tagr
  real(dp),allocatable :: matrix_tmp1(:,:)
  real(dp),allocatable :: matrix_tmp2(:,:)
  !=====

  gt_tagl = get_gaussian_type_tag(gaussian_type_left)
  gt_tagr = get_gaussian_type_tag(gaussian_type_right)
  n1c = number_basis_function_am('CART',am1)
  n2c = number_basis_function_am('CART',am2)
  n3c = number_basis_function_am('CART',am3)
  n1  = number_basis_function_am(gaussian_type_left,am1)
  n2  = number_basis_function_am(gaussian_type_right,am2)
  n3  = number_basis_function_am(gaussian_type_right,am3)

  if( .NOT. ALLOCATED(matrix_out) ) allocate(matrix_out(n1,n2,n3))

  if( gt_tagl /= gt_tagr ) then
    call die('transform_libcint_to_molgw_3d: mixed pure/cart integrals not coded')
  endif

  ! When CART -> CART, just normalize (no unitary transform needed)
  if( gt_tagl == CARTG .AND. gt_tagr == CARTG ) then
    ii = 0
    do i1=1,n1
      do i2=1,n2
        do i3=1,n3
          ii = ii + 1
          matrix_out(i1,i2,i3) = array_in(ii) * cart_to_pure_norm(am1,CARTG)%matrix(i1,i1) &
                                              * cart_to_pure_norm(am2,CARTG)%matrix(i2,i2) &
                                              * cart_to_pure_norm(am3,CARTG)%matrix(i3,i3)
        enddo
      enddo
    enddo
  else
    ! PURE PURE
    ii = 0
    do i1=1,n1
      do i2=1,n2
        do i3=1,n3
          ii = ii + 1
          matrix_out(i1,i2,i3) = array_in(ii) * libcint_pure_norm(am1) &
                                              * libcint_pure_norm(am2) &
                                              * libcint_pure_norm(am3)
        enddo
      enddo
    enddo
  endif

end subroutine transform_libcint_to_molgw_3d


!=========================================================================
end module m_libcint_tools
!=========================================================================
