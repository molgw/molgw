!!****m* DoNOF/m_vars
!! NAME
!! m_vars
!!
!! FUNCTION
!! Define size vars e.g. real(dp).
!!
!! COPYRIGHT
!!
!! NOTES
!!
!! SOURCE

module m_vars

 implicit none
 integer, parameter :: dp = kind(1.d0)
 real(dp), parameter :: zero = 0.0d0
 real(dp), parameter :: half = 0.5d0
 real(dp), parameter :: one  = 1.0d0
 real(dp), parameter :: two  = 2.0d0
 real(dp), parameter :: four = 4.0d0
 real(dp), parameter :: ten  = 1.0d1
 real(dp), parameter :: tol1 = 1.0d-1
 real(dp), parameter :: tol3 = 1.0d-3
 real(dp), parameter :: tol5 = 1.0d-5
 real(dp), parameter :: tol6 = 1.0d-6
 real(dp), parameter :: tol8 = 1.0d-8
 real(dp), parameter :: tol9 = 1.0d-9
 real(dp), parameter :: tol16 = 1.0d-16
 real(dp), parameter :: tol20 = 1.0d-20
 real(dp), parameter :: thousand = 1.0d3
 real(dp), parameter :: pi   = dacos(-one)

end module m_vars
!!***
