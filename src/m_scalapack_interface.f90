!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This file contains
! explicit interfaces to ensure correct SCALAPACK calls
!
!=========================================================================
module m_scalapack_interface
  use,intrinsic :: ISO_FORTRAN_ENV, only: REAL32, REAL64

  !interface pxsyevd
  !  module procedure PSSYEVD
  !  module procedure PDSYEVD
  !end interface

  interface
    subroutine pssyevd( JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, INFO )
      import REAL32
      character(len=1), intent(in) :: JOBZ, UPLO
      integer, intent(in) :: N, IA, JA, IZ, JZ, LWORK, LIWORK
      integer, intent(out) :: INFO
      integer, intent(in) :: DESCA(*), DESCZ(*), IWORK(*)
      real(kind=REAL32), intent(inout) :: A(*), W(*), Z(*), WORK(*)
    end subroutine pssyevd

    subroutine pdsyevd( JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, INFO )
      import REAL64
      character(len=1), intent(in) :: JOBZ, UPLO
      integer, intent(in) :: N, IA, JA, IZ, JZ, LWORK, LIWORK
      integer, intent(out) :: INFO
      integer, intent(in) :: DESCA(*), DESCZ(*), IWORK(*)
      real(kind=REAL64), intent(inout) :: A(*), W(*), Z(*), WORK(*)
    end subroutine pdsyevd

    subroutine pdgemm( TRANSA, TRANSB, M, N, K, &
               ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB, &
               BETA, C, IC, JC, DESCC )
      import REAL64
      character(len=1), intent(in) :: TRANSA, TRANSB
      integer, intent(in) :: M, N, K, IA, JA, IB, JB, IC, JC
      integer, intent(in) :: DESCA(*), DESCB(*), DESCC(*)
      real(kind=REAL64), intent(in) :: ALPHA, BETA
      real(kind=REAL64), intent(in) :: A(*), B(*)
      real(kind=REAL64), intent(out) :: C(*)
    end subroutine pdgemm
  end interface

end module m_scalapack_interface
