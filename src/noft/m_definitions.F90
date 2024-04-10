!!****m* DoNOF/m_definitions
!! NAME
!! m_definitions
!!
!! FUNCTION
!! Define size vars e.g. real(dp).
!!
!! COPYRIGHT
!!
!! NOTES
!!
!! SOURCE

module m_definitions

 implicit none
 integer, parameter :: dp = kind(1.d0)
 real(dp), parameter :: zero = 0.0d0
 real(dp), parameter :: half = 0.5d0
 real(dp), parameter :: one  = 1.0d0
 real(dp), parameter :: two  = 2.0d0
 real(dp), parameter :: four = 4.0d0
 real(dp), parameter :: eight= 8.0d0
 real(dp), parameter :: ten  = 1.0d1
 real(dp), parameter :: twelve = 1.2d1
 real(dp), parameter :: tol1 = 1.0d-1
 real(dp), parameter :: tol3 = 1.0d-3
 real(dp), parameter :: tol5 = 1.0d-5
 real(dp), parameter :: tol6 = 1.0d-6
 real(dp), parameter :: tol8 = 1.0d-8
 real(dp), parameter :: tol9 = 1.0d-9
 real(dp), parameter :: tol10 = 1.0d-10
 real(dp), parameter :: tol16 = 1.0d-16
 real(dp), parameter :: tol20 = 1.0d-20
 real(dp), parameter :: thousand = 1.0d3
 real(dp), parameter :: pi = dacos(-one)
 real(dp), parameter :: Ha_eV = 27.21138505
 complex(dp),parameter  :: im    = (0.0_dp,1.0_dp)
 complex(dp),parameter  :: complex_one  = (1.0_dp,0.0_dp)
 complex(dp),parameter  :: complex_zero = (0.0_dp,0.0_dp)

end module m_definitions
!!***
