module m_solid_harmonics
  use m_definitions
  implicit none


contains

  pure function eval_l0(r, alpha, c) result(shell)
    real(dp), intent(in) :: r(3), alpha(:), c(:)
    real(dp) :: x, y, z
    real(dp) :: shell(1)
    real(dp) :: r2, expo
    integer :: i, n

    n = SIZE(alpha)
    x = r(1)
    y = r(2)
    z = r(3)
    r2 = x*x + y*y + z*z
    expo = 0.0_dp
    do i=1, n
      expo = expo + exp(-alpha(i) * r2) * c(i)
    enddo
    shell(1) = expo

  end function eval_l0

  pure function eval_l1(r, alpha, c) result(shell)
    real(dp), intent(in) :: r(3), alpha(:), c(:)
    real(dp) :: x, y, z
    real(dp) :: shell(3)
    real(dp) :: r2, expo
    integer :: i, n

    n = SIZE(alpha)
    x = r(1)
    y = r(2)
    z = r(3)
    r2 = x*x + y*y + z*z
    expo = 0.0_dp
    do i=1, n
      expo = expo + exp(-alpha(i) * r2) * c(i)
    enddo

    shell(1) = y * expo
    shell(2) = z * expo
    shell(3) = x * expo

  end function eval_l1

  pure function eval_l2(r, alpha, c) result(shell)
    real(dp), intent(in) :: r(3), alpha(:), c(:)
    real(dp) :: x, y, z
    real(dp) :: shell(5)
    real(dp) :: r2, expo, N2
    real(dp) :: ang0,ang1,ang2c,ang2s
    integer :: i, n

    n = SIZE(alpha)
    x = r(1)
    y = r(2)
    z = r(3)
    r2 = x*x + y*y + z*z
    expo = 0.0_dp
    do i=1, n
      expo = expo + exp(-alpha(i) * r2) * c(i)
    enddo

    ang0  = SQRT(1.0 /12.0)  ! sqrt(5.0_dp /(16*pi))
    ang1  = SQRT(1.0 / 1.0)  ! sqrt(15.0_dp/(4*pi))
    ang2c = SQRT(1.0 / 4.0)  ! sqrt(15.0_dp/(16*pi))
    ang2s = SQRT(1.0 / 1.0)  ! sqrt(15.0_dp/(4*pi))

    shell(1) = ang0  * (3*z*z - r2) * expo
    shell(2) = ang1  * x*z         * expo
    shell(3) = ang1  * y*z         * expo
    shell(4) = ang2c * (x*x - y*y) * expo
    shell(5) = ang2s * x*y         * expo
  end function eval_l2

  pure function eval_l3(r, alpha, c) result(shell)
    real(dp), intent(in) :: r(3), alpha(:), c(:)
    real(dp) :: x, y, z
    real(dp) :: shell(7)
    real(dp) :: r2, expo, N3
    real(dp) :: a0,a1,a2,a2s,a3
    integer  :: i, n

    n = SIZE(alpha)
    x = r(1)
    y = r(2)
    z = r(3)
    r2 = x*x + y*y + z*z
    expo = 0.0_dp
    do i=1, n
      expo = expo + exp(-alpha(i) * r2) * c(i)
    enddo

    a0  = sqrt(1.0_dp/60.)
    a1  = sqrt(1.0_dp/40.)
    a2  = sqrt(1.0_dp/4.)
    a2s = sqrt(1.0_dp/1)
    a3  = sqrt(1.0_dp/24)

    shell(1) = a0  * (5*z*z*z - 3*z*r2)      * expo
    shell(2) = a1  * x*(5*z*z - r2)          * expo
    shell(3) = a1  * y*(5*z*z - r2)          * expo
    shell(4) = a2  * z*(x*x - y*y)           * expo
    shell(5) = a2s * z*x*y                   * expo
    shell(6) = a3  * (x*x*x - 3*x*y*y)       * expo
    shell(7) = a3  * (3*x*x*y - y*y*y)       * expo
  end function eval_l3

  pure function eval_l4(r, alpha, c) result(shell)
    real(dp), intent(in) :: r(3), alpha(:), c(:)
    real(dp) :: x, y, z
    real(8) :: shell(9)
    real(8) :: r2, expo
    real(8) :: a0,a1,a1s,a2,a2s,a3,a3s,a4
    integer  :: i, n
    
    n = SIZE(alpha)
    x = r(1)
    y = r(2)
    z = r(3)
    r2 = x*x + y*y + z*z
    expo = 0.0_dp
    do i=1, n
      expo = expo + exp(-alpha(i) * r2) * c(i)
    enddo
    
    a0 =  sqrt(1.0_dp/6720)
    a1 =  sqrt(1.0_dp/168)
    a2 =  sqrt(1.0_dp/336)
    a3 =  sqrt(1.0_dp/24) 
    a4 =  sqrt(1.0_dp/192)
    
    
    shell(1) = a0 * (35*z**4 - 30*z*z*r2 + 3*r2*r2) * expo
    shell(2) = a1 * z*(7*z*z - 3*r2)*x * expo
    shell(3) = a1 * z*(7*z*z - 3*r2)*y * expo
    shell(4) = a2 * (x*x - y*y)*(7*z*z - r2) * expo
    shell(5) = a2 * (2*x*y) *(7*z*z - r2) * expo
    shell(6) = a3 * z*(x*x*x - 3*x*y*y) * expo
    shell(7) = a3 * z*(3*x*x*y - y*y*y) * expo
    shell(8) = a4 * (x**4 - 6*x*x*y*y + y**4) * expo
    shell(9) = a4 * (4*x*x*x*y - 4*x*y*y*y) * expo
  end function eval_l4
  
  
  pure function eval_l5(r, alpha, c) result(shell)
    real(dp), intent(in) :: r(3), alpha(:), c(:)
    real(dp) :: x, y, z
    real(8) :: shell(11)
    real(8) :: r2, expo
    real(8) :: a0,a1,a1s,a2,a2s,a3,a3s,a4,a4s,a5,a5s
    integer  :: i, n
    
    n = SIZE(alpha)
    x = r(1)
    y = r(2)
    z = r(3)
    r2 = x*x + y*y + z*z
    expo = 0.0_dp
    do i=1, n
      expo = expo + exp(-alpha(i) * r2) * c(i)
    enddo
    
    a0 = sqrt(1.0_dp/60480)
    a1 = sqrt(1.0_dp/36288)
    a2 = sqrt(1.0_dp/13392)
    a3 = sqrt(1.0_dp/23040)
    a4 = sqrt(1.0_dp/192)
    a5 = sqrt(1.0_dp/1920)
    
    shell(1) = a0 * z*(63*z**4 - 70*z*z*r2 + 15*r2*r2) * expo
    shell(2) = a1 * x*(63*z**4 - 42*z*z*r2 + 3*r2*r2) * expo
    shell(3) = a1 * y*(63*z**4 - 42*z*z*r2 + 3*r2*r2) * expo
    shell(4) = a2 * z*(x*x - y*y)*(21*z*z - 3*r2) * expo
    shell(5) = a2 * z*(2*x*y) *(21*z*z - 3*r2) * expo
    shell(6) = a3 * (x*x*x - 3*x*y*y)*(21*z*z - r2) * expo
    shell(7) = a3 * (3*x*x*y - y*y*y)*(21*z*z - r2) * expo
    shell(8) = a4 * z*(x**4 - 6*x*x*y*y + y**4) * expo
    shell(9) = a4 * z*(4*x*x*x*y - 4*x*y*y*y) * expo
    shell(10) = a5 * (x**5 - 10*x*x*x*y*y + 5*x*y*y*y*y) * expo
    shell(11) = a5 * (5*x*x*x*x*y - 10*x*x*y*y*y + y**5) * expo
  end function eval_l5

end module m_solid_harmonics

