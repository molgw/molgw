!=========================================================================
! Test implementation of quadrature in Fourier space to get the integrals for
! - overlap()
! - kinetic
! - electron-nucleus()
!=========================================================================

module m_fourier_quadrature
  use m_definitions
  use m_warning
  use m_timing
  use m_atoms
  use m_basis_set
  use m_inputparam
  use m_hamiltonian_onebody



contains

!=========================================================================
subroutine setup_overlap_fourier(basis_p,basis_t,reference)
  implicit none
  type(basis_set),intent(in) :: basis_p    ! basis set for projectile
  type(basis_set),intent(in) :: basis_t    ! basis set for target
  real(dp),intent(in)        :: reference(basis_p%nbf,basis_t%nbf)
  !======
  real(dp)                :: velocity(3)
  complex(dp),allocatable :: s_matrix_v(:,:)
  !=====


  ! TODO the velocity should be read from the basis_p%bff(1)%velocity(:)
  velocity(:) = 0.0_dp

  !
  !  setup_gos_ao evaluates  < \phi_a | e^{i q.r} | \phi_b >
  !
  call setup_gos_ao(basis_p,basis_t,-velocity,s_matrix_v)


  if( basis_t%nbf < 6 .OR. basis_p%nbf < 6 ) return

  write(stdout,*) ' === Home-made analytic evaluation for velocity',velocity(:)
  write(stdout,'(*(1x,f12.8))') s_matrix_v(1,1:6)%re
  write(stdout,'(*(1x,f12.8))') s_matrix_v(2,1:6)%re
  write(stdout,'(*(1x,f12.8))') s_matrix_v(3,1:6)%re
  write(stdout,'(*(1x,f12.8))') s_matrix_v(4,1:6)%re
  write(stdout,'(*(1x,f12.8))') s_matrix_v(5,1:6)%re
  write(stdout,'(*(1x,f12.8))') s_matrix_v(6,1:6)%re
  write(stdout,*)

  write(stdout,*) ' === Reference'
  write(stdout,'(*(1x,f12.8))') reference(1,1:6)
  write(stdout,'(*(1x,f12.8))') reference(2,1:6)
  write(stdout,'(*(1x,f12.8))') reference(3,1:6)
  write(stdout,'(*(1x,f12.8))') reference(4,1:6)
  write(stdout,'(*(1x,f12.8))') reference(5,1:6)
  write(stdout,'(*(1x,f12.8))') reference(6,1:6)
  write(stdout,*)

  deallocate(s_matrix_v)

end subroutine setup_overlap_fourier


!=========================================================================
subroutine setup_nucleus_fourier(basis_p,basis_t,reference)
  implicit none

  type(basis_set),intent(in) :: basis_p    ! basis set for projectile
  type(basis_set),intent(in) :: basis_t    ! basis set for target
  real(dp),intent(in)        :: reference(basis_p%nbf,basis_t%nbf)
  !======
  real(dp) :: velocity(3)
  integer,parameter :: n1=230
  integer,parameter :: nqradial=60
  integer,parameter :: nq=n1*nqradial
  integer :: ix1,iqradial,iq,iatom
  real(dp) :: xtmp,weight
  real(dp) :: qlist(3,nq),wq(nq),qvec(3),qpvvec(3)
  real(dp) :: x1(n1),y1(n1),z1(n1),w1(n1)
  real(dp) :: xa(nqradial),wxa(nqradial)
  real(dp),parameter :: alpha= 6.0_dp
  complex(dp) :: structure_factor
  complex(dp),allocatable :: s_matrix_mqmv(:,:)
  complex(dp) :: enucl(basis_p%nbf,basis_t%nbf)
  !=====

  call start_clock(timing_tmp1)

  ! TODO the velocity should be read from the basis_t%bff(1)%velocity(:)
  velocity(:) = 0.0_dp
  !velocity(3) = 0.5_dp

  write(stdout,*) 'Fourier-space evaluation of nucleus-electron'
  write(stdout,*) 'velocity:',velocity(:)

  do iqradial=1,nqradial
    xtmp = ( iqradial - 0.5_dp ) / REAL(nqradial,dp)
    xa(iqradial)   = -alpha * log( 1.0_dp - xtmp**3)
    wxa(iqradial)  = 3.0_dp * alpha * xtmp**2 / ( 1.0_dp - xtmp**3 ) / REAL(nqradial,dp)
  enddo
  call ld0230(x1,y1,z1,w1,iq)
  iq = 0
  do ix1=1,n1
    do iqradial=1,nqradial
      iq = iq + 1
      wq(iq) = wxa(iqradial) * w1(ix1) * 4.0_dp * pi * xa(iqradial)**2
      qlist(1,iq) = xa(iqradial) * x1(ix1)
      qlist(2,iq) = xa(iqradial) * y1(ix1)
      qlist(3,iq) = xa(iqradial) * z1(ix1)
    enddo
  enddo


  enucl(:,:) = 0.0_dp
  do iq=1,nq
    qvec(:) = qlist(:,iq)
    qpvvec(:) = qvec(:) + velocity(:)
    weight = wq(iq)

    structure_factor = (0.0_dp,0.0_dp)
    do iatom=1,natom
      structure_factor = structure_factor &
                - zvalence(iatom) * EXP( im * DOT_PRODUCT(qvec(:),xatom(:,iatom)) )
    enddo

    call setup_gos_ao(basis_p,basis_t,-qpvvec,s_matrix_mqmv)

    enucl(:,:) = enucl(:,:) + weight * 4.0_dp * pi / NORM2(qvec)**2 &
                                   * s_matrix_mqmv(:,:) * structure_factor

    !write(1000,*) NORM2(qvec), ABS( weight * 4.0_dp * pi / NORM2(qvec)**2 * s_matrix_mqmv(1,1) * structure_factor )

    deallocate(s_matrix_mqmv)

  !enddo
  !write(1000,*)
  enddo
  enucl(:,:) = enucl(:,:) / (2.0_dp * pi)**3

  call stop_clock(timing_tmp1)

  if( basis_t%nbf < 6 .OR. basis_p%nbf < 6 ) return

  write(stdout,*) ' === Numerical evaluation with points:',n1,nqradial
  write(stdout,*) ' Real part'
  write(stdout,'(*(1x,f12.8))') enucl(1,1:6)%re
  write(stdout,'(*(1x,f12.8))') enucl(2,1:6)%re
  write(stdout,'(*(1x,f12.8))') enucl(3,1:6)%re
  write(stdout,'(*(1x,f12.8))') enucl(4,1:6)%re
  write(stdout,'(*(1x,f12.8))') enucl(5,1:6)%re
  write(stdout,'(*(1x,f12.8))') enucl(6,1:6)%re
  write(stdout,*) ' Imaginary part'
  write(stdout,'(*(1x,f12.8))') enucl(1,1:6)%im
  write(stdout,'(*(1x,f12.8))') enucl(2,1:6)%im
  write(stdout,'(*(1x,f12.8))') enucl(3,1:6)%im
  write(stdout,'(*(1x,f12.8))') enucl(4,1:6)%im
  write(stdout,'(*(1x,f12.8))') enucl(5,1:6)%im
  write(stdout,'(*(1x,f12.8))') enucl(6,1:6)%im
  write(stdout,*)

  write(stdout,*) ' === Reference'
  write(stdout,'(*(1x,f12.8))') reference(1,1:6)
  write(stdout,'(*(1x,f12.8))') reference(2,1:6)
  write(stdout,'(*(1x,f12.8))') reference(3,1:6)
  write(stdout,'(*(1x,f12.8))') reference(4,1:6)
  write(stdout,'(*(1x,f12.8))') reference(5,1:6)
  write(stdout,'(*(1x,f12.8))') reference(6,1:6)
  write(stdout,*)

end subroutine setup_nucleus_fourier


!=========================================================================
subroutine setup_kinetic_fourier(basis_p,basis_t,reference)
  implicit none

  type(basis_set),intent(in) :: basis_p    ! basis set for projectile
  type(basis_set),intent(in) :: basis_t    ! basis set for target
  real(dp),intent(in)        :: reference(basis_p%nbf,basis_t%nbf)
  !======
  integer,parameter :: n1=230
  integer,parameter :: nqradial=60
  integer,parameter :: nq=n1*nqradial
  real(dp) :: velocity(3)
  integer :: ix1,iqradial,iq,ishell,ibf,jbf
  real(dp) :: xtmp,weight
  real(dp) :: qlist(3,nq),wq(nq),qvec(3),qmvvec(3)
  real(dp) :: x1(n1),y1(n1),z1(n1),w1(n1)
  real(dp) :: xa(nqradial),wxa(nqradial)
  real(dp),parameter :: alpha= 5.0_dp
  integer                    :: gt,li,ni_cart,ibf1,ibf1_cart,ibf2,i_cart
  complex(dp)                :: basis_function_q(basis_t%nbf)
  complex(dp)                :: basis_function_qmv(basis_p%nbf)
  complex(dp)                :: ekin(basis_p%nbf,basis_t%nbf)
  complex(dp),allocatable    :: basis_function_cart(:)
  !=====

  gt = get_gaussian_type_tag(basis_p%gaussian_type)

  ! TODO the velocity should be read from the basis_t%bff(1)%velocity(:)
  velocity(:) = 0.0_dp
  !velocity(3) = 0.5_dp

  write(stdout,*) 'Fourier-space evaluation of kinetic energy'
  write(stdout,*) 'velocity:',velocity(:)

  ! Setup the 3D grid in Fourier space
  do iqradial=1,nqradial
    xtmp = ( iqradial - 0.5_dp ) / REAL(nqradial,dp)
    xa(iqradial)   = -alpha * log( 1.0_dp - xtmp**3)
    wxa(iqradial)  = 3.0_dp * alpha * xtmp**2 / ( 1.0_dp - xtmp**3 ) / REAL(nqradial,dp)
  enddo
  call ld0230(x1,y1,z1,w1,iq)
  iq = 0
  do iqradial=1,nqradial
    do ix1=1,n1
      iq = iq + 1
      wq(iq) = wxa(iqradial) * w1(ix1) * 4.0_dp * pi * xa(iqradial)**2
      qlist(1,iq) = xa(iqradial) * x1(ix1)
      qlist(2,iq) = xa(iqradial) * y1(ix1)
      qlist(3,iq) = xa(iqradial) * z1(ix1)
    enddo
  enddo


  ekin(:,:) = 0.0_dp
  do iq=1,nq
    qvec(:) = qlist(:,iq)
    qmvvec(:) = qvec(:) - velocity(:)
    weight = wq(iq)

    !
    ! Fourier transform of the right-hand basis function: \phi_beta(q)
    do ishell=1,basis_t%nshell
      li        = basis_t%shell(ishell)%am
      ni_cart   = number_basis_function_am('CART',li)
      ibf1      = basis_t%shell(ishell)%istart
      ibf1_cart = basis_t%shell(ishell)%istart_cart
      ibf2      = basis_t%shell(ishell)%iend

      allocate(basis_function_cart(ni_cart))

      do i_cart=1,ni_cart
        basis_function_cart(i_cart) = basis_function_fourier(basis_t%bfc(ibf1_cart+i_cart-1),qvec)
      enddo
      basis_function_q(ibf1:ibf2) = MATMUL( TRANSPOSE(cart_to_pure(li,gt)%matrix(:,:)), basis_function_cart(:) )
      deallocate(basis_function_cart)

    enddo

    !
    ! Fourier transform of the left-hand basis function: \phi_alpha(q-v)
    do ishell=1,basis_p%nshell
      li        = basis_p%shell(ishell)%am
      ni_cart   = number_basis_function_am('CART',li)
      ibf1      = basis_p%shell(ishell)%istart
      ibf1_cart = basis_p%shell(ishell)%istart_cart
      ibf2      = basis_p%shell(ishell)%iend

      allocate(basis_function_cart(ni_cart))

      do i_cart=1,ni_cart
        basis_function_cart(i_cart) = basis_function_fourier(basis_p%bfc(ibf1_cart+i_cart-1),qmvvec)
      enddo
      basis_function_qmv(ibf1:ibf2) = MATMUL( TRANSPOSE(cart_to_pure(li,gt)%matrix(:,:)), basis_function_cart(:) )
      deallocate(basis_function_cart)

    enddo

    !ekin(ibf,jbf) = ekin(ibf,jbf) + weight * NORM2(qvec)**2 &
    !                                * basis_function_qmv(ibf) * CONJG( basis_function_q(jbf) )

    call  ZGERC(basis_p%nbf,basis_t%nbf,CMPLX(weight*NORM2(qvec)**2,0.0_dp), &
                basis_function_qmv,1,basis_function_q,1,ekin,basis_p%nbf)


  enddo

  ekin(:,:) = 0.5_dp * CONJG( ekin(:,:) ) * (2.0_dp * pi)**3

  if( basis_t%nbf < 6 .OR. basis_p%nbf < 6 ) return

  write(stdout,*) ' === Numerical evaluation with points:',n1,nqradial
  write(stdout,*) ' Real part'
  write(stdout,'(*(1x,f12.8))') ekin(1,1:6)%re
  write(stdout,'(*(1x,f12.8))') ekin(2,1:6)%re
  write(stdout,'(*(1x,f12.8))') ekin(3,1:6)%re
  write(stdout,'(*(1x,f12.8))') ekin(4,1:6)%re
  write(stdout,'(*(1x,f12.8))') ekin(5,1:6)%re
  write(stdout,'(*(1x,f12.8))') ekin(6,1:6)%re
  write(stdout,*)
  write(stdout,*) ' Imaginary part'
  write(stdout,'(*(1x,f12.8))') ekin(1,1:6)%im
  write(stdout,'(*(1x,f12.8))') ekin(2,1:6)%im
  write(stdout,'(*(1x,f12.8))') ekin(3,1:6)%im
  write(stdout,'(*(1x,f12.8))') ekin(4,1:6)%im
  write(stdout,'(*(1x,f12.8))') ekin(5,1:6)%im
  write(stdout,'(*(1x,f12.8))') ekin(6,1:6)%im
  write(stdout,*)

  write(stdout,*) ' === Reference'
  write(stdout,'(*(1x,f12.8))') reference(1,1:6)
  write(stdout,'(*(1x,f12.8))') reference(2,1:6)
  write(stdout,'(*(1x,f12.8))') reference(3,1:6)
  write(stdout,'(*(1x,f12.8))') reference(4,1:6)
  write(stdout,'(*(1x,f12.8))') reference(5,1:6)
  write(stdout,'(*(1x,f12.8))') reference(6,1:6)
  write(stdout,*)

end subroutine setup_kinetic_fourier


!=========================================================================
subroutine setup_gos_ao(basis_p,basis_t,qvec,gos_ao)
  implicit none
  type(basis_set),intent(in)          :: basis_p
  type(basis_set),intent(in)          :: basis_t
  real(dp),intent(in)                 :: qvec(3)
  complex(dp),allocatable,intent(out) :: gos_ao(:,:)
  !=====
  integer                 :: gt
  integer                 :: ishell,jshell
  integer                 :: ibf1,ibf2,jbf1,jbf2,ibf1_cart,jbf1_cart
  integer                 :: li,lj,ni_cart,nj_cart,i_cart,j_cart
  complex(dp),allocatable :: gos_cart(:,:)
  !=====

  gt = get_gaussian_type_tag(basis_p%gaussian_type)

  allocate(gos_ao(basis_p%nbf,basis_t%nbf))

  !$OMP PARALLEL PRIVATE(li,lj,ni_cart,nj_cart,ibf1,ibf1_cart,ibf2,jbf1,jbf1_cart,jbf2,gos_cart)
  !$OMP DO
  do jshell=1,basis_t%nshell
    lj        = basis_t%shell(jshell)%am
    nj_cart   = number_basis_function_am('CART',lj)
    jbf1      = basis_t%shell(jshell)%istart
    jbf1_cart = basis_t%shell(jshell)%istart_cart
    jbf2      = basis_t%shell(jshell)%iend

    do ishell=1,basis_p%nshell
      li        = basis_p%shell(ishell)%am
      ni_cart   = number_basis_function_am('CART',li)
      ibf1      = basis_p%shell(ishell)%istart
      ibf1_cart = basis_p%shell(ishell)%istart_cart
      ibf2      = basis_p%shell(ishell)%iend


      allocate(gos_cart(ni_cart,nj_cart))

      do i_cart=1,ni_cart
        do j_cart=1,nj_cart
          call basis_function_gos(basis_p%bfc(ibf1_cart+i_cart-1),basis_t%bfc(jbf1_cart+j_cart-1),qvec,gos_cart(i_cart,j_cart))
        enddo
      enddo

      gos_ao(ibf1:ibf2,jbf1:jbf2) = MATMUL( TRANSPOSE( cart_to_pure(li,gt)%matrix(:,:) ) , &
              MATMUL(  gos_cart(:,:) , cart_to_pure(lj,gt)%matrix(:,:) ) )

      deallocate(gos_cart)

    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine setup_gos_ao


!=========================================================================
subroutine evaluate_gos_cart(ga,gb,qvec,gos_cart)
  implicit none

  type(gaussian),intent(in) :: ga,gb
  real(dp),intent(in)         :: qvec(3)
  complex(dp),intent(out)     :: gos_cart
  !=====
  real(dp)    :: apb,q2,fa,fb,fab
  complex(dp) :: factor,ff,gg
  complex(dp) :: intx,inty,intz
  !=====

  ! this routine is buggy ! Be careful

  apb = ga%alpha + gb%alpha

  q2 = SUM(qvec(:)**2)


  factor = ( pi / apb )**1.5_dp * EXP( (-ga%alpha * gb%alpha * SUM((ga%x0(:)-gb%x0(:))**2) &
                  + im * DOT_PRODUCT( qvec , ga%alpha * ga%x0(:) + gb%alpha * gb%x0(:) ) &
                  - 0.25_dp * q2 ) / apb )

  fa  = 0.5_dp / ga%alpha
  fb  = 0.5_dp / gb%alpha
  fab = 2.0_dp * ga%alpha * gb%alpha / apb

  ! X direction
  ff = fab * ( gb%x0(1) - ga%x0(1) ) + im * ga%alpha / apb * qvec(1)
  gg = fab * ( ga%x0(1) - gb%x0(1) ) + im * gb%alpha / apb * qvec(1)
  intx = auxiliary(ga%nx,gb%nx,ff,gg,fa,fb,fab)
  ! Y direction
  ff = fab * ( gb%x0(2) - ga%x0(2) ) + im * ga%alpha / apb * qvec(2)
  gg = fab * ( ga%x0(2) - gb%x0(2) ) + im * gb%alpha / apb * qvec(2)
  inty = auxiliary(ga%ny,gb%ny,ff,gg,fa,fb,fab)
  ! Z direction
  ff = fab * ( gb%x0(3) - ga%x0(3) ) + im * ga%alpha / apb * qvec(3)
  gg = fab * ( ga%x0(3) - gb%x0(3) ) + im * gb%alpha / apb * qvec(3)
  intz = auxiliary(ga%nz,gb%nz,ff,gg,fa,fb,fab)

  gos_cart = factor * intx * inty * intz * ga%norm_factor * gb%norm_factor

contains

function auxiliary(ia,ib,ff,gg,fa,fb,fab)
  integer,intent(in) :: ia,ib
  real(dp),intent(in) :: fa,fb,fab
  complex(dp),intent(out) :: ff,gg
  complex(dp) :: auxiliary
  !=====
  complex(dp),parameter :: one = (1.0_dp,0.0_dp)
  !=====

  select case(ia*10+ib)
  case(00)
    auxiliary = one
  case(01)
    auxiliary = gg * fb
  case(10)
    auxiliary = ff * fa
  case(11)
    auxiliary = ff * gg * fa * fb + fa * fb * fab
  case default
    call die('not implemented')
  end select

end function auxiliary

end subroutine evaluate_gos_cart

end module m_fourier_quadrature


!=========================================================================
