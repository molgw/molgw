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
  integer(C_INT),protected,allocatable :: atm(:,:)
  integer(C_INT),protected,allocatable :: bas(:,:)
  real(C_DOUBLE),protected,allocatable :: env(:)


  integer,external :: cint1e_ovlp_cart


 interface transform_libcint_to_molgw
   module procedure transform_libcint_to_molgw_2d
   !module procedure transform_libcint_to_molgw_3d
   !module procedure transform_libcint_to_molgw_4d
 end interface


contains


!=========================================================================
subroutine init_libcint(basis)
  implicit none
  type(basis_set),intent(in) :: basis
  !=====
  integer :: off,icenter_basis,ishell
  !=====

  write(stdout,'(/,1x,a)') 'Initialize LIBCINT internal data'

  LIBCINT_natm = ncenter_basis
  LIBCINT_nbas = basis%nshell
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
    !env(off+1:off+basis%shell(ishell)%ng) = basis%shell(ishell)%coeff(:) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * basis%shell(ishell)%am + 1 )
    select case(basis%shell(ishell)%am)
    case(0,1)
      env(off+1:off+basis%shell(ishell)%ng) = basis%shell(ishell)%coeff(:) * SQRT(4.0_dp * pi) / SQRT( 2.0_dp * basis%shell(ishell)%am + 1 )
    case default
      env(off+1:off+basis%shell(ishell)%ng) = basis%shell(ishell)%coeff(:)
    end select
    !CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+1))
    !env(off + 1) = .7 * CINTgto_norm(bas(ANG_OF,ishell), env(bas(PTR_EXP,ishell)+1))
    !env(off + 2) = .6 * CINTgto_norm(bas(ANG_OF,ishell), env(bas(PTR_EXP,ishell)+2))
    !env(off + 3) = .5 * CINTgto_norm(bas(ANG_OF,ishell), env(bas(PTR_EXP,ishell)+3))
    !env(off + 4) = .4 * CINTgto_norm(bas(ANG_OF,ishell), env(bas(PTR_EXP,ishell)+1))
    !env(off + 5) = .3 * CINTgto_norm(bas(ANG_OF,ishell), env(bas(PTR_EXP,ishell)+2))
    !env(off + 6) = .2 * CINTgto_norm(bas(ANG_OF,ishell), env(bas(PTR_EXP,ishell)+3))
    off = off + basis%shell(ishell)%ng
  enddo



end subroutine init_libcint


#if FALSE
!=========================================================================
! DEMO
subroutine setup_overlap_libcint(basis,s_matrix)
 implicit none
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: s_matrix(basis%nbf,basis%nbf)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable :: array_cart(:)
 integer(C_INT)             :: amA,contrdepthA
 real(C_DOUBLE)             :: A(3)
 real(C_DOUBLE),allocatable :: alphaA(:)
 real(C_DOUBLE),allocatable :: cA(:)
 integer(C_INT)             :: amB,contrdepthB
 real(C_DOUBLE)             :: B(3)
 real(C_DOUBLE),allocatable :: alphaB(:)
 real(C_DOUBLE),allocatable :: cB(:)
!=====
 integer :: i_cart,j_cart,ij
 integer :: ibf_cart,jbf_cart
 integer :: info
 integer :: shls(2)
!=====
 call start_clock(timing_overlap)
#if defined(HAVE_LIBCINT)
 write(stdout,'(/,a)') ' Setup overlap matrix S (LIBCINT)'
#else
 call die('HAVE_LIBCINT is false')
#endif

 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=jshell,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart(ni_cart*nj_cart))

#if defined(HAVE_LIBCINT)
     shls(1) = ishell-1  ! C convention starts with 0
     shls(2) = jshell-1  ! C convention starts with 0
     info = cint1e_ovlp_cart(array_cart, shls, atm, LIBCINT_natm, bas, LIBCINT_nbas, env)


     !call transform_libint_to_molgw(basis%gaussian_type,li,lj,RESHAPE(TRANSPOSE(array_cart),[ni_cart*nj_cart]),matrix)
     call transform_libcint_to_molgw_2d(basis%gaussian_type,li,lj,array_cart,matrix)
#endif

     deallocate(alphaA,cA)


     s_matrix(ibf1:ibf2,jbf1:jbf2) = matrix(:,:)
     s_matrix(jbf1:jbf2,ibf1:ibf2) = TRANSPOSE(matrix(:,:))

     deallocate(array_cart,matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 title='=== Overlap matrix S (LIBCINT) ==='
 call dump_out_matrix(.FALSE.,title,s_matrix)


 call stop_clock(timing_overlap)


end subroutine setup_overlap_libcint
#endif


!=========================================================================
subroutine transform_libcint_to_molgw_2d(gaussian_type,am1,am2,array_in,matrix_out)
 implicit none
 character(len=4),intent(in)      :: gaussian_type
 integer,intent(in)               :: am1,am2
 real(C_DOUBLE),intent(in)        :: array_in(:)
 real(dp),allocatable,intent(out) :: matrix_out(:,:)
!=====
 integer :: n1,n2,n1c,n2c
 integer :: gt_tag
 real(dp),allocatable :: matrix_tmp(:,:)
!=====

 gt_tag = get_gaussian_type_tag(gaussian_type)
 n1c = number_basis_function_am('CART',am1)
 n2c = number_basis_function_am('CART',am2)
 n1  = number_basis_function_am(gaussian_type,am1)
 n2  = number_basis_function_am(gaussian_type,am2)

 if( .NOT. ALLOCATED(matrix_out) ) allocate(matrix_out(n1,n2))

 allocate(matrix_tmp(n1c,n2))
 ! Transform the right index
 matrix_tmp(:,:) = MATMUL( RESHAPE(array_in(:),[n1c,n2c]) , cart_to_pure_norm(am2,gt_tag)%matrix(:,:) )
 ! Transform the left index
 matrix_out(:,:) = MATMUL( TRANSPOSE(cart_to_pure_norm(am1,gt_tag)%matrix(1:n1c,1:n1)), matrix_tmp(:,:) )

 deallocate(matrix_tmp)

end subroutine transform_libcint_to_molgw_2d


end module m_libcint_tools
