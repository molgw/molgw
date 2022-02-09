module m_libcint_tools
  use m_definitions
  use m_cart_to_pure
  use m_basis_set
  use m_atoms

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

  integer,private,parameter :: LIBCINT_env_size=10000

  integer(C_INT),protected :: LIBCINT_natm,LIBCINT_nbas
  integer(C_INT),protected :: LIBCINT_AUXIL_BASIS_START
  integer(C_INT),protected,allocatable :: atm(:,:)
  integer(C_INT),protected,allocatable :: bas(:,:)
  real(C_DOUBLE),protected,allocatable :: env(:)

  logical,protected :: libcint_has_range_separation
  integer,external  :: cint2e_cart
  integer,external  :: cint2c2e_cart
  integer,external  :: cint3c2e_cart

  interface

    integer(C_INT) function cint1e_ovlp_cart(array_cart, shls, atm, natm, bas, nbas, env) bind(C)
      import :: C_INT,C_DOUBLE
      integer(C_INT),value  :: natm,nbas
      real(C_DOUBLE),intent(in) :: env(*)
      integer(C_INT),intent(in) :: shls(*),atm(*),bas(*)
      real(C_DOUBLE),intent(out) :: array_cart(*)
    end function cint1e_ovlp_cart

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



contains


!=========================================================================
! Check with a known integral whether LIBCINT considers LIBCINT_PTR_RANGE_OMEGA
!
subroutine check_capability_libcint()
  implicit none
  !=====
  real(C_DOUBLE) :: fake_env(50),integral(1)
  integer(C_INT) :: fake_atm(LIBCINT_ATM_SLOTS,1)
  integer(C_INT) :: fake_bas(LIBCINT_BAS_SLOTS,1)
  integer(C_INT) :: shls(2)
  integer        :: info,off
  real(dp),parameter :: reference_value = 3.8052728203379578E-002_dp
  !=====

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
  info = cint2c2e_cart(integral, shls, fake_atm, 1_C_INT, fake_bas, 1_C_INT, fake_env, 0_C_LONG)

  libcint_has_range_separation = ABS( integral(1) - reference_value ) < 1.0e-12_dp 


end subroutine check_capability_libcint


!=========================================================================
subroutine set_rinv_origin_libcint(x0)
  implicit none
  real(dp),intent(in) :: x0(3)
  !=====
  !=====

  env(LIBCINT_PTR_RINV_ORIG+1:LIBCINT_PTR_RINV_ORIG+3) = x0(:)

end subroutine set_rinv_origin_libcint


!=========================================================================
subroutine set_erf_screening_length_libcint(rcut)
  implicit none
  real(dp),intent(in) :: rcut
  !=====
  !=====

  ! LIBCINT convention: a positive omega = 1/rcut means long-range only
  if( rcut > 1.0e-12_dp ) then
    ! Ensure that LIBCINT can calculated "erf" range-separated Coulomb interaction
    if( .NOT. libcint_has_range_separation ) then
       call die('')
    endif
    env(LIBCINT_PTR_RANGE_OMEGA+1) = 1.0_dp / rcut
  else
    env(LIBCINT_PTR_RANGE_OMEGA+1) = 0.0_dp 
  endif

end subroutine set_erf_screening_length_libcint


!=========================================================================
subroutine init_libcint(basis,auxil_basis)
  implicit none
  type(basis_set),intent(in)          :: basis
  type(basis_set),optional,intent(in) :: auxil_basis
  !=====
  integer :: off,icenter_basis,ishell
  !=====

  if( PRESENT(auxil_basis) ) then
    write(stdout,'(/,1x,a)') 'Initialize LIBCINT internal data: basis and auxiliary basis'
  else
    write(stdout,'(/,1x,a)') 'Initialize LIBCINT internal data: basis'
  endif

  LIBCINT_natm = ncenter_basis
  LIBCINT_nbas = basis%nshell
  if( PRESENT(auxil_basis) ) then
    LIBCINT_nbas = LIBCINT_nbas + auxil_basis%nshell
  endif

  allocate(atm(LIBCINT_ATM_SLOTS,LIBCINT_natm))
  allocate(bas(LIBCINT_BAS_SLOTS,LIBCINT_nbas))
  allocate(env(LIBCINT_env_size))

  ! real space origin
  env(LIBCINT_PTR_COMMON_ORIG+1) = 0.0_C_DOUBLE
  env(LIBCINT_PTR_COMMON_ORIG+2) = 0.0_C_DOUBLE
  env(LIBCINT_PTR_COMMON_ORIG+3) = 0.0_C_DOUBLE
  off = LIBCINT_PTR_ENV_START

  do icenter_basis=1,ncenter_basis
    atm(LIBCINT_CHARGE_OF,icenter_basis) = zbasis(icenter_basis)
    atm(LIBCINT_PTR_COORD,icenter_basis) = off ! note the 0-based index
    env(off+1:off+3) = xbasis(:,icenter_basis)
    off = off + 3
  enddo

  do ishell=1,basis%nshell
    bas(LIBCINT_ATOM_OF  ,ishell)  = basis%shell(ishell)%icenter - 1 ! C convention starts with 0
    bas(LIBCINT_ANG_OF   ,ishell)  = basis%shell(ishell)%am
    bas(LIBCINT_NPRIM_OF ,ishell)  = basis%shell(ishell)%ng
    bas(LIBCINT_NCTR_OF  ,ishell)  = 1   ! so far take shells one at a time (even if they share exponents)
    bas(LIBCINT_PTR_EXP  ,ishell)  = off ! note the 0-based index
    env(off+1:off+basis%shell(ishell)%ng) = basis%shell(ishell)%alpha(:)
    off = off + basis%shell(ishell)%ng
    bas(LIBCINT_PTR_COEFF,ishell) = off
    select case(basis%shell(ishell)%am)
    case(0,1)
      env(off+1:off+basis%shell(ishell)%ng) = basis%shell(ishell)%coeff(:) &
                * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * basis%shell(ishell)%am + 1 )
    case default
      env(off+1:off+basis%shell(ishell)%ng) = basis%shell(ishell)%coeff(:)
    end select
    off = off + basis%shell(ishell)%ng
  enddo

  LIBCINT_AUXIL_BASIS_START = basis%nshell

  if( PRESENT(auxil_basis) ) then
    do ishell=1,auxil_basis%nshell
      bas(LIBCINT_ATOM_OF  ,LIBCINT_AUXIL_BASIS_START + ishell)  = auxil_basis%shell(ishell)%icenter - 1 ! C convention starts with 0
      bas(LIBCINT_ANG_OF   ,LIBCINT_AUXIL_BASIS_START + ishell)  = auxil_basis%shell(ishell)%am
      bas(LIBCINT_NPRIM_OF ,LIBCINT_AUXIL_BASIS_START + ishell)  = auxil_basis%shell(ishell)%ng
      bas(LIBCINT_NCTR_OF  ,LIBCINT_AUXIL_BASIS_START + ishell)  = 1   ! so far take shells one at a time (even if they share exponents)
      bas(LIBCINT_PTR_EXP  ,LIBCINT_AUXIL_BASIS_START + ishell)  = off ! note the 0-based index
      env(off+1:off+auxil_basis%shell(ishell)%ng) = auxil_basis%shell(ishell)%alpha(:)
      off = off + auxil_basis%shell(ishell)%ng
      bas(LIBCINT_PTR_COEFF,LIBCINT_AUXIL_BASIS_START + ishell) = off
      select case(auxil_basis%shell(ishell)%am)
      case(0,1)
        env(off+1:off+auxil_basis%shell(ishell)%ng) = auxil_basis%shell(ishell)%coeff(:) &
                  * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * auxil_basis%shell(ishell)%am + 1 )
      case default
        env(off+1:off+auxil_basis%shell(ishell)%ng) = auxil_basis%shell(ishell)%coeff(:)
      end select
      off = off + auxil_basis%shell(ishell)%ng
    enddo
  endif


end subroutine init_libcint


!=========================================================================
subroutine destroy_libcint()
  implicit none

  LIBCINT_natm = 0
  LIBCINT_nbas = 0
  LIBCINT_AUXIL_BASIS_START = 0
  if( ALLOCATED(atm) ) deallocate(atm)
  if( ALLOCATED(bas) ) deallocate(bas)
  if( ALLOCATED(env) ) deallocate(env)

end subroutine destroy_libcint


!=========================================================================
end module m_libcint_tools
!=========================================================================
