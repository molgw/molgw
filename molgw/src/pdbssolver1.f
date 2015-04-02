      SUBROUTINE PDBSSOLVER1( N, M, IM, JM, DESCM, K, IK, JK, DESCK,
     $                        LAMBDA, X, IX, JX, DESCX, WORK, LWORK,
     $                        IWORK, LIWORK, INFO )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N, IM, JM, IK, JK, IX, JX, LWORK, LIWORK, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            DESCM( * ), DESCK( * ), DESCX( * ), IWORK( * )
      DOUBLE PRECISION   M( * ), K( * ), LAMBDA( * ), X( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDBSSOLVER1() computes all eigenvalues and (both right and
*  left) eigenvectors of 2n-by-2n real matrix
*
*     H = [ A,  B;
*          -B, -A ],
*
*  where both A and B are n-by-n symmetric.
*
*  On entry, the information of H is provided in the lower triangular
*  parts of two n-by-n matrices M = A+B and K = A-B.
*
*  The matrices M and K are required to be positive definite.
*
*  The structure of H leads to the following properties.
*
*     1) H is diagonalizable with n pairs of real eigenvalues
*        (lambda_i, -lambda_i).
*
*     2) The eigenvectors of H has the block structure
*
*           X = [ X_1, X_2;     Y = [ X_1, -X_2;
*                 X_2, X_1 ],        -X_2,  X_1 ],
*
*        and satisfy that
*
*           X_1**T * X_2 = X_2**T * X_1,
*           X_1**T * X_1 - X_2**T * X_2 = I,
*           Y**T * X = I,
*           H * X = X * diag(lambda, -lambda),
*           Y**T * H = diag(lambda, -lambda) * Y**T.
*
*  On exit, only the positive eigenvalues and the corresponding right
*  eigenvectors are returned.  The eigenvalues are sorted in ascending
*  order.  The eigenvectors are normalized (i.e., X = [ X_1; X_2 ] with
*  X_1**T * X_1 - X_2**T * X_2 = I).
*
*  M and K are destroyed on exit.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The number of rows and columns of A and B.
*          The order of the distributed submatrices sub( M ), sub( X ),
*          and sub( Y ) are 2*N.
*          N >= 0.
*
*  M       (local input and output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension
*          (LLD_M, LOCc(JM+N-1)).
*          On entry, the symmetric positive definite matrix M. Only the
*          lower triangular part of M is used to define the elements of
*          the symmetric matrix.
*          On exit, all entries of M are destroyed.
*
*  IM      (global input) INTEGER
*          The row index in the global array M indicating the first
*          row of sub( M ).
*
*  JM      (global input) INTEGER
*          The column index in the global array M indicating the
*          first column of sub( M ).
*
*  DESCM   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix M.
*
*  K       (local input and output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension
*          (LLD_K, LOCc(JK+N-1)).
*          On entry, the symmetric positive definite matrix K. Only the
*          lower triangular part of K is used to define the elements of
*          the symmetric matrix.
*          On exit, all entries of K are destroyed.
*
*  IK      (global input) INTEGER
*          The row index in the global array K indicating the first
*          row of sub( K ).
*
*  JK      (global input) INTEGER
*          The column index in the global array K indicating the
*          first column of sub( K ).
*
*  DESCK   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix K.
*
*  LAMBDA  (global output) DOUBLE PRECISION array, dimension (N)
*          On normal exit LAMBDA contains the positive eigenvalues of H
*          in ascending order.
*
*  X       (local output) DOUBLE PRECISION array,
*          global dimension (2N, N),
*          local dimension ( LLD_X, LOCc(JX+N-1) )
*          On normal exit X contains the right eigenvectors of H.
*
*  IX      (global input) INTEGER
*          X's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*          In this version, only IX = JX = 1 is supported.
*
*  JX      (global input) INTEGER
*          X's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*          In this version, only IX = JX = 1 is supported.
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*          DESCX( CTXT_ ) must equal DESCA( CTXT_ )
*
*  WORK    (local workspace/output) DOUBLE PRECISION array,
*          dimension (LWORK)
*          On output, WORK( 1 ) returns the minimal amount of workspace
*          needed to guarantee completion.
*          If the input parameters are incorrect, WORK( 1 ) may also be
*          incorrect.
*
*  LWORK   (local input) INTEGER
*          The length of the workspace array WORK.
*          If LWORK = -1, the LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          size for the WORK/IWORK array. The required workspace is
*          returned as the first element of WORK/IWORK and no error
*          message is issued by PXERBLA.
*
*  IWORK   (local workspace/output) INTEGER array,
*          dimension (LIWORK)
*          On output, IWORK( 1 ) returns the minimal amount of workspace
*          needed to guarantee completion.
*          If the input parameters are incorrect, IWORK( 1 ) may also be
*          incorrect.
*
*  LIWORK   (local input) INTEGER
*          The length of the workspace array IWORK.
*          If LIWORK = -1, the LIWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          size for the WORK/IWORK array. The required workspace is
*          returned as the first element of WORK/IWORK and no error
*          message is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  The eigensolver did not converge.
*
*  Alignment requirements
*  ======================
*
*  This subroutine requires M and K to be distributed identically in the
*  sense that DESCM( : ) = DESCK( : ).
*  In addition, square blocks ( i.e., MB = NB ) are required.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            ICTXT, NPROCS, NPROW, NPCOL, MYROW, MYCOL, NB,
     $                   TWON, I, J, LWKOPT, LLWORK, LOCALMAT, MROWS,
     $                   MCOLS, LLDM, INDPHI, INDPSI, INDV, INDWORK,
     $                   ITMP, DIMV, NZ, LIWKOPT
      DOUBLE PRECISION   DTMP
      DOUBLE PRECISION   T_CHOL, T_FORMW, T_DIAG, T_VEC1, T_VEC2, T_PREP
*     ..
*     .. Local Arrays ..
      INTEGER            DESCPHI( DLEN_ ), DESCPSI( DLEN_ ),
     $                   DESCV( DLEN_ )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, SQRT
*     ..
*     .. External Functions ..
      EXTERNAL           NUMROC, MPI_WTIME
      INTEGER            NUMROC
      DOUBLE PRECISION   MPI_WTIME
*     ..
*     .. External Subroutines ..
      EXTERNAL           PDAXPY, PDCOPY, PDSCAL, PDLACPY, PDPOTRF,
     $                   PDSYEVR, PDSYGST, PXERBLA, BLACS_GRIDINFO,
     $                   CHK1MAT, PCHK2MAT
*     ..
*     .. Executable Statements ..
*
      T_PREP = MPI_WTIME()
      INFO = 0
      TWON = 2*N
      ICTXT = DESCM( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW*NPCOL
      IF ( NPROW .EQ. -1 ) THEN
         INFO = -( 200+CTXT_ )
      END IF
*
*     Test the input arguments.
*
      IF ( INFO .EQ. 0 .AND. N .LT. 0 ) THEN
         INFO = -1
      END IF
      IF ( INFO .EQ. 0 )
     $   CALL CHK1MAT( N, 3, N, 3, IM, JM, DESCM, 5, INFO )
      IF ( INFO .EQ. 0 .AND. DESCM( MB_ ) .NE. DESCM( NB_ ) )
     $   INFO = -( 500+MB_ )
      IF ( INFO .EQ. 0 )
     $   CALL CHK1MAT( N, 3, N, 3, IK, JK, DESCK, 9, INFO )
      IF ( INFO .EQ. 0 .AND. DESCK( MB_ ) .NE. DESCK( NB_ ) )
     $   INFO = -( 900+MB_ )
      IF ( INFO .EQ. 0)
     $   CALL PCHK2MAT( N, 3, N, 3, IM, JM, DESCM, 5, N, 3, N, 3,
     $        IK, JK, DESCK, 9, 0, ITMP, ITMP, INFO )
      IF ( INFO .EQ. 0 )
     $   CALL CHK1MAT( TWON, 3, N, 3, IX, JX, DESCX, 14, INFO )
      IF ( INFO .EQ. 0 .AND. DESCX( MB_ ) .NE. DESCX( NB_ ) )
     $   INFO = -( 1000+MB_ )
*
*     Compute required workspace.
*
      IF ( INFO .EQ. 0 ) THEN
*
*        Set up local indices for the workspace.
*
         LQUERY = LWORK .EQ. -1
         NB = DESCM( NB_ )
         LLDM = DESCM( LLD_ )
         MROWS = NUMROC( TWON, NB, MYROW, 0, NPROW )
         MCOLS = NUMROC( TWON, NB, MYCOL, 0, NPCOL )
         LOCALMAT = MAX( NB, LLDM )*MCOLS
         DESCPHI( 1:DLEN_ ) = DESCM( 1:DLEN_ )
         DESCPSI( 1:DLEN_ ) = DESCM( 1:DLEN_ )
         DESCV( 1:DLEN_ ) = DESCM( 1:DLEN_ )
         INDPHI = 1
         INDPSI = INDPHI + LOCALMAT
         INDWORK = INDPSI + LOCALMAT
         INDV = INDPHI
         LLWORK = LWORK - INDWORK + 1
*
*        Estimate the workspace required by external subroutines.
*
!         CALL PDSYEV( 'V', 'L', N, DTMP, IK, JK, DESCK, DTMP, DTMP, 1,
!     $        1, DESCV, WORK, -1, ITMP )
         CALL PDSYEVR( 'V', 'A', 'L', N, DTMP, IK, JK, DESCK, ZERO,
     $        ZERO, 1, N, DIMV, NZ, DTMP, DTMP, 1, 1, DESCV, WORK, -1,
     $        IWORK, -1, ITMP )
         LWKOPT = INT( WORK( 1 ) )
         LIWKOPT = IWORK( 1 )
*
         LWKOPT = INDWORK - 1 + MAX( LWKOPT, LOCALMAT )
         IF ( .NOT. LQUERY .AND. LWORK .LT. LWKOPT )
     $      INFO = -16
         IF ( INFO .EQ. 0 .AND. .NOT. LQUERY .AND. LIWORK .LT. LIWKOPT )
     $      INFO = -18
      END IF
*
      IF ( INFO .NE. 0 ) THEN
         CALL PXERBLA( ICTXT, 'PDBSSOLVER1', -INFO )
         RETURN
      END IF
      WORK( 1 ) = DBLE( LWKOPT )
      IF ( LQUERY )
     $   RETURN
*
*     Quick return if possible.
*
      IF ( N .EQ. 0 )
     $   RETURN
      T_PREP = MPI_WTIME() - T_PREP
      IF ( MYROW+MYCOL .EQ. 0 )
     $   WRITE( *, * ) 't_prep = ', T_PREP, ';'
*
*     Compute the Cholesky factorization M = L * L**T.
*     In case of failure, a general solver is needed.
*
      T_CHOL = MPI_WTIME()
      CALL PDPOTRF( 'L', N, M, IM, JM, DESCM, ITMP )
      T_CHOL = MPI_WTIME() - T_CHOL
      IF ( MYROW+MYCOL .EQ. 0 )
     $   WRITE( *, * ) 't_chol = ', T_CHOL, ';'
!      CALL PDLAPRNT( N, N, M, IM, JM, DESCM, 0, 0, 'L', 6,
!     $     WORK( INDWORK ) )
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) "L=tril(L);"
      IF ( ITMP .NE. 0 ) THEN
         INFO = -2
         RETURN
      END IF
*
*     Explicitly formulate W = L**T * K * L.
*     Only the lower triangular part of W is filled by the correct
*     values.
*
      T_FORMW = MPI_WTIME()
      CALL PDSYGST( 3, 'L', N, K, IK, JK, DESCK, M, IM, JM, DESCM, DTMP,
     $     ITMP )
      T_FORMW = MPI_WTIME() - T_FORMW
      IF ( MYROW+MYCOL .EQ. 0 )
     $   WRITE( *, * ) 't_formw = ', T_FORMW, ';'
!      CALL PDLAPRNT( N, N, K, IK, JK, DESCK, 0, 0, 'W',
!     $     6, WORK( INDWORK ) )
*
*     Diagonalization: V**T * (L**T * K * L) * V = diag(lambda).
*
      T_DIAG = MPI_WTIME()
!      CALL PDSYEV( 'V', 'L', N, K, IK, JK, DESCK, LAMBDA, WORK( INDV ),
!     $     1, 1, DESCV, WORK( INDWORK ), LLWORK, ITMP )
      CALL PDSYEVR( 'V', 'A', 'L', N, K, IK, JK, DESCK, ZERO, ZERO,
     $     1, N, DIMV, NZ, LAMBDA, WORK( INDV ), 1, 1, DESCV,
     $     WORK( INDWORK ), LLWORK, IWORK, LIWORK, ITMP )
      T_DIAG = MPI_WTIME() - T_DIAG
      IF ( MYROW+MYCOL .EQ. 0 )
     $   WRITE( *, * ) 't_diag = ', T_DIAG, ';'
      IF ( ITMP .NE. 0 ) THEN
         INFO = ITMP
         WRITE( *, * ), '% PDSYEVR fails with INFO =', INFO
         RETURN
      END IF
      DO I = 1, N
         LAMBDA( I ) = DSQRT( LAMBDA( I ) )
      END DO
!      IF ( MYROW+MYCOL .EQ. 0 ) THEN
!         DO I = 1, N
!            WRITE( *, * ) 'lambda0(', I, ', 1) =', LAMBDA( I ), ';'
!         END DO
!      END IF
!      CALL PDLAPRNT( N, N, WORK( INDV ), 1, 1, DESCV, 0, 0, 'V',
!     $     6, WORK( INDWORK ) )
*
*     Recover the eigenvectors:
*
*        X = [ ( Psi * Lambda**{1/2} + Phi * Lambda**{-1/2} ) / 2
*              ( Psi * Lambda**{1/2} - Phi * Lambda**{-1/2} ) / 2 ]
*
*     where
*
*        Phi = L * V, Psi = L**{-T} * V.
*
      T_VEC1 = MPI_WTIME()
      CALL PDLACPY( 'A', N, N, WORK( INDV ), 1, 1, DESCV,
     $     WORK( INDPSI ), 1, 1, DESCPSI )
      CALL PDTRMM( 'L', 'L', 'N', 'N', N, N, ONE, M, IM, JM, DESCM,
     $     WORK( INDPHI ), 1, 1, DESCPHI )
      CALL PDTRSM( 'L', 'L', 'T', 'N', N, N, ONE, M, IM, JM, DESCM,
     $     WORK( INDPSI ), 1, 1, DESCPSI )
*
*     Scale Psi and Phi.
*
      DO I = 1, N
         CALL PDSCAL( N, HALF*DSQRT( LAMBDA( I ) ),
     $        WORK( INDPSI ), 1, I, DESCPSI, 1 )
      END DO
      DO I = 1, N
         CALL PDSCAL( N, HALF/DSQRT( LAMBDA( I ) ),
     $        WORK( INDPHI ), 1, I, DESCPHI, 1 )
      END DO
      T_VEC1 = MPI_WTIME() - T_VEC1
      IF ( MYROW+MYCOL .EQ. 0 )
     $   WRITE( *, * ) 't_vec1 = ', T_VEC1, ';'
!      CALL PDLAPRNT( N, N, WORK( INDPHI ), 1, 1, DESCPHI, 0, 0, 'Phi',
!     $     6, WORK( INDWORK ) )
!      CALL PDLAPRNT( N, N, WORK( INDPSI ), 1, 1, DESCPSI, 0, 0, 'Psi',
!     $     6, WORK( INDWORK ) )
*
*     Construct X.
*
      T_VEC2 = MPI_WTIME()
      CALL PDGEADD( 'N', N, N, ONE, WORK( INDPSI ), 1, 1, DESCPSI,
     $     ZERO, X, IX, JX, DESCX )
      CALL PDGEADD( 'N', N, N, ONE, WORK( INDPSI ), 1, 1, DESCPSI,
     $     ZERO, X, IX+N, JX, DESCX )
      CALL PDGEADD( 'N', N, N, ONE, WORK( INDPHI ), 1, 1, DESCPHI,
     $     ONE, X, IX, JX, DESCX )
      CALL PDGEADD( 'N', N, N, -ONE, WORK( INDPHI ), 1, 1, DESCPHI,
     $     ONE, X, IX+N, JX, DESCX )
      T_VEC2 = MPI_WTIME() - T_VEC2
      IF ( MYROW+MYCOL .EQ. 0 )
     $   WRITE( *, * ) 't_vec2 = ', T_VEC2, ';'
!      CALL PDLAPRNT( TWON, N, X, IX, JX, DESCX, 0, 0, 'X', 6,
!     $     WORK( INDWORK ) )
*
      WORK( 1 ) = DBLE( LWKOPT )
*
      RETURN
*
*     End of PDBSSOLVER1().
*
      END
