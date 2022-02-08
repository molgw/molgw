module m_libcint_tools
  use m_definitions
  use m_cart_to_pure
  use m_basis_set
  use m_atoms


  integer,private,parameter :: LIBCINT_CHARGE_OF  = 1
  integer,private,parameter :: LIBCINT_PTR_COORD  = 2
  integer,private,parameter :: LIBCINT_NUC_MOD_OF = 3
  integer,private,parameter :: LIBCINT_PTR_ZETA   = 4
  integer,private,parameter :: LIBCINT_ATM_SLOTS  = 6

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
  integer,private,parameter :: LIBCINT_PTR_ENV_START = 20

  integer,private,parameter :: LIBCINT_env_size=10000

  integer(C_INT),protected :: LIBCINT_natm,LIBCINT_nbas
  integer(C_INT),protected :: LIBCINT_AUXIL_BASIS_START
  integer(C_INT),protected,allocatable :: atm(:,:)
  integer(C_INT),protected,allocatable :: bas(:,:)
  real(C_DOUBLE),protected,allocatable :: env(:)

  integer,external :: cint2e_cart
  integer,external :: cint2c2e_cart,cint3c2e_cart

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
    env(off+1:off+3) =  xbasis(:,icenter_basis)
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
