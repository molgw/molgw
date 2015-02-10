#ifdef HAVE_SCALAPACK
      SUBROUTINE PDBSSOLVER1_SVD( N, M, IM, JM, DESCM, K, IK, JK, DESCK,
     $                            LAMBDA, X, IX, JX, DESCX, Y, IY, JY,
     $                            DESCY, WORK, LWORK, INFO )
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            N, IM, JM, IK, JK, IX, JX, IY, JY, LWORK, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            DESCM( * ), DESCK( * ), DESCX( * ), DESCY( * )
      DOUBLE PRECISION   M( * ), K( * ), LAMBDA( * ), X( * ), Y( * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDBSSOLVER1_SVD() computes all eigenvalues and (both right and
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
*  On exit, the outputs X, Y, and lambda satisfy that
*
*     Y**T * X = I,
*     H * X = X * diag(lambda),
*     Y**T * H = diag(lambda) * Y**T.
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
*  LAMBDA  (global output) DOUBLE PRECISION array, dimension (2*N)
*
*  X       (local output) DOUBLE PRECISION array,
*          global dimension (2N, 2N),
*          local dimension ( LLD_X, LOCc(JX+2*N-1) )
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
*  Y       (local output) DOUBLE PRECISION array,
*          global dimension (2N, 2N),
*          local dimension ( LLD_Y, LOCc(JY+2*N-1) )
*          On normal exit Y contains the left eigenvectors of H.
*
*  IY      (global input) INTEGER
*          Y's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*          In this version, only IY = JY = 1 is supported.
*
*  JY      (global input) INTEGER
*          Y's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*          In this version, only IY = JY = 1 is supported.
*
*  DESCY   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Y.
*          DESCY( CTXT_ ) must equal DESCA( CTXT_ )
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
*          size for the WORK array. The required workspace is returned
*          as the first element of WORK and no error message is issued
*          by PXERBLA.
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
*  This subroutine requires (M,K) and (X,Y), respectively, to be
*  distributed identically in the sense that DESCM( : ) = DESCK( : ),
*  DESCX( : )= DESCY( : ), and DESCM( MB_ ) = DESCX( MB_ ).
*  In addition, square blocks ( i.e., MB_M = NB_M ) are required.
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
     $                   TWON, I, J, LWKOPT, LLWORK, LOCALMAT, WROWS,
     $                   WCOLS, LLD, INDS, INDU, INDVT, INDW, INDWORK,
     $                   ITMP
      DOUBLE PRECISION   DTMP
      DOUBLE PRECISION   T_CHOL, T_FORMW, T_DIAG, T_VEC, T_SYM, T_PREP,
     $                   T
*     ..
*     .. Local Arrays ..
      INTEGER            DESCU( DLEN_ ), DESCVT( DLEN_ ), DESCW( DLEN_ )
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
      EXTERNAL           PDAXPY, PDCOPY, PDSWAP, PDTRMM, PDTRSM, PDSCAL,
     $                   PDLACPY, PDLASET, PDGEADD, PDPOTRF, PDGESVD,
     $                   PXERBLA, BLACS_GRIDINFO, DESCSET, CHK1MAT,
     $                   PCHK2MAT
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
     $   CALL CHK1MAT( TWON, 3, TWON, 3, IX, JX, DESCX, 14, INFO )
      IF ( INFO .EQ. 0 .AND. DESCX( MB_ ) .NE. DESCX( NB_ ) )
     $   INFO = -( 1000+MB_ )
      IF ( INFO .EQ. 0 )
     $   CALL CHK1MAT( TWON, 3, TWON, 3, IY, JY, DESCY, 18, INFO )
      IF ( INFO .EQ. 0 .AND. DESCY( MB_ ) .NE. DESCY( NB_ ) )
     $   INFO = -( 1400+MB_ )
      IF ( INFO .EQ. 0 )
     $   CALL PCHK2MAT( TWON, 3, TWON, 3, IX, JX, DESCX, 14, TWON, 3,
     $        TWON, 3, IY, JY, DESCY, 18, 0, ITMP, ITMP, INFO )

*
*     Compute required workspace.
*
      IF ( INFO .EQ. 0 ) THEN
*
*        Set up local indices for the workspace.
*
         LQUERY = LWORK .EQ. -1
         NB = DESCM( NB_ )
         WROWS = NUMROC( N, NB, MYROW, 0, NPROW )
         WCOLS = NUMROC( N, NB, MYCOL, 0, NPCOL )
         LLD = MAX( WROWS, 1 )
         LOCALMAT = MAX( NB, LLD )*WCOLS
         CALL DESCSET( DESCU, N, N, NB, NB, 0, 0, ICTXT, LLD )
         CALL DESCSET( DESCVT, N, N, NB, NB, 0, 0, ICTXT, LLD )
         CALL DESCSET( DESCW, N, N, NB, NB, 0, 0, ICTXT, LLD )
         INDW = 1
         INDU = INDW + LOCALMAT
         INDVT = INDU + LOCALMAT
         INDS = INDVT + LOCALMAT
         INDWORK = INDS + N
         LLWORK = LWORK - INDWORK + 1
*
*        Estimate the workspace required by external subroutines.
*
         CALL PDGESVD( 'V', 'V', N, N, DTMP, 1, 1, DESCW, DTMP, DTMP, 1,
     $        1, DESCU, DTMP, 1, 1, DESCVT, WORK, -1, ITMP )
         LWKOPT = INT( WORK( 1 ) )
*
         LWKOPT = INDWORK - 1 + MAX( LWKOPT, LOCALMAT )
         IF ( .NOT. LQUERY .AND. LWORK .LT. LWKOPT )
     $      INFO = -20
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
*     Compute the Cholesky factorizations
*        M = L1 * L1**T, K = L2 * L2**T
*     In case of failure, a general solver is needed.
*
      T = MPI_WTIME()
      CALL PDPOTRF( 'L', N, M, IM, JM, DESCM, ITMP )
      T_CHOL = MPI_WTIME() - T
      IF ( ITMP .NE. 0 ) THEN
         INFO = -2
         IF ( MYROW+MYCOL .EQ. 0 )
     $      WRITE( *, * ) 't_chol = ', T_CHOL, ';'
         RETURN
      END IF
      T = MPI_WTIME()
      CALL PDPOTRF( 'L', N, K, IK, JK, DESCK, ITMP )
      T_CHOL = T_CHOL + MPI_WTIME() - T
      IF ( MYROW+MYCOL .EQ. 0 )
     $   WRITE( *, * ) 't_chol = ', T_CHOL, ';'
      IF ( ITMP .NE. 0 ) THEN
         INFO = -6
      END IF
!      CALL PDLAPRNT( N, N, M, IM, JM, DESCM, 0, 0, 'L1', 6,
!     $     WORK( INDWORK ) )
!      CALL PDLAPRNT( N, N, K, IK, JK, DESCK, 0, 0, 'L2', 6,
!     $     WORK( INDWORK ) )
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) "L1=tril(L1); L2=tril(L2);"
*
*     Explicitly formulate W = L2**T * L1.
*
      T_FORMW = MPI_WTIME()
      CALL PDLASET( 'U', N, N, ZERO, ZERO, WORK( INDW ), 1, 1, DESCW )
      CALL PDLACPY( 'L', N, N, M, IM, JM, DESCM, WORK( INDW ), 1, 1,
     $     DESCW )
      CALL PDTRMM( 'L', 'L', 'T', 'N', N, N, ONE, K, IK, JK, DESCK,
     $     WORK( INDW ), 1, 1, DESCW )
      T_FORMW = MPI_WTIME() - T_FORMW
      IF ( MYROW+MYCOL .EQ. 0 )
     $   WRITE( *, * ) 't_formw = ', T_FORMW, ';'
!      CALL PDLAPRNT( N, N, WORK( INDW ), 1, 1, DESCW, 0, 0, 'W',
!     $     6, WORK( INDWORK ) )
*
*     Compute the SVD: W = U * S * V**T.
*
      T_DIAG = MPI_WTIME()
      CALL PDGESVD( 'V', 'V', N, N, WORK( INDW ), 1, 1, DESCW,
     $     WORK( INDS ), WORK( INDU ), 1, 1, DESCU, WORK( INDVT ), 1, 1,
     $     DESCVT, WORK( INDWORK ), LLWORK, ITMP )
      T_DIAG = MPI_WTIME() - T_DIAG
      IF ( ITMP .NE. 0 ) THEN
!        It has been observed that ITMP may be equal to N+1, indicating
!        inconsistency among processors.
         IF ( MYROW+MYCOL .EQ. 0 )
     $      WRITE( *, * ) 't_diag = ', T_DIAG, ';'
         INFO = ITMP
         WRITE( *, * ), '% PDGESVD fails with INFO =', INFO
!         RETURN
      END IF
!      CALL PDLAPRNT( N, N, WORK( INDU ), 1, 1, DESCU, 0, 0, 'U',
!     $     6, WORK( INDWORK ) )
!      DO I = 1, N
!         WRITE( *, * ) 'S(', I, ')= ', WORK( INDS+I-1 ), ';'
!      END DO
!      CALL PDLAPRNT( N, N, WORK( INDVT ), 1, 1, DESCVT, 0, 0, 'VT',
!     $     6, WORK( INDWORK ) )
*
*     Copy S to lambda in ascending order.
*     Adjust the order in U and V accordingly.
*
      T = MPI_WTIME()
      DO I = 1, N
         LAMBDA( N+1-I ) = WORK( INDS+I-1 )
      END DO
      DO I = 1, N/2
         CALL PDSWAP( N, WORK( INDU ), 1, I, DESCU, 1,
     $        WORK( INDU ), 1, N+1-I, DESCU, 1 )
         CALL PDSWAP( N, WORK( INDVT ), I, 1, DESCVT, N,
     $        WORK( INDVT ), N+1-I, 1, DESCVT, N )
      END DO
      T_DIAG = T_DIAG + MPI_WTIME() - T
      IF ( MYROW+MYCOL .EQ. 0 )
     $   WRITE( *, * ) 't_diag = ', T_DIAG, ';'
*
*     Recover the eigenvectors X and Y:
*
*        X = [ L1**(-T) * V + L2**(-T) * U
*              L1**(-T) * V - L2**(-T) * U ],
*
*        Y = [ L1 * V + L2 * U
*              L1 * V - L2 * U ].
*
      T = MPI_WTIME()
      CALL PDLACPY( 'A', N, N, WORK( INDVT ), 1, 1, DESCVT,
     $     WORK( INDW ), 1, 1, DESCW )
      CALL PDTRSM( 'R', 'L', 'N', 'N', N, N, HALF, M, IM, JM, DESCM,
     $     WORK( INDW ), 1, 1, DESCW )
      CALL PDGEADD( 'T', N, N, ONE, WORK( INDW ), 1, 1, DESCW, ZERO, X,
     $     IX, JX, DESCX )
      CALL PDGEADD( 'T', N, N, ONE, WORK( INDW ), 1, 1, DESCW, ZERO, X,
     $     IX+N, JX, DESCX )
      CALL PDLACPY( 'A', N, N, WORK( INDU ), 1, 1, DESCU,
     $     WORK( INDW ), 1, 1, DESCW )
      CALL PDTRSM( 'L', 'L', 'T', 'N', N, N, HALF, K, IK, JK, DESCK,
     $     WORK( INDW ), 1, 1, DESCW )
      CALL PDGEADD( 'N', N, N, ONE, WORK( INDW ), 1, 1, DESCW, ONE, X,
     $     IX, JX, DESCX )
      CALL PDGEADD( 'N', N, N, -ONE, WORK( INDW ), 1, 1, DESCW, ONE, X,
     $     IX+N, JX, DESCX )
*
      CALL PDLACPY( 'A', N, N, WORK( INDVT ), 1, 1, DESCVT,
     $     WORK( INDW ), 1, 1, DESCW )
      CALL PDTRMM( 'R', 'L', 'T', 'N', N, N, HALF, M, IM, JM, DESCM,
     $     WORK( INDW ), 1, 1, DESCW )
      CALL PDGEADD( 'T', N, N, ONE, WORK( INDW ), 1, 1, DESCW, ZERO, Y,
     $     IY, JY, DESCY )
      CALL PDGEADD( 'T', N, N, ONE, WORK( INDW ), 1, 1, DESCW, ZERO, Y,
     $     IY+N, JY, DESCY )
      CALL PDLACPY( 'A', N, N, WORK( INDU ), 1, 1, DESCU,
     $     WORK( INDW ), 1, 1, DESCW )
      CALL PDTRMM( 'L', 'L', 'N', 'N', N, N, HALF, K, IK, JK, DESCK,
     $     WORK( INDW ), 1, 1, DESCW )
      CALL PDGEADD( 'N', N, N, ONE, WORK( INDW ), 1, 1, DESCW, ONE, Y,
     $     IY, JY, DESCY )
      CALL PDGEADD( 'N', N, N, -ONE, WORK( INDW ), 1, 1, DESCW, ONE, Y,
     $     IY+N, JY, DESCY )
      T_VEC = T_VEC + MPI_WTIME() - T
      IF ( MYROW+MYCOL .EQ. 0 )
     $   WRITE( *, * ) 't_vec = ', T_VEC, ';'
*
*     Recover negative eigenvalues and the eigenvectors.
*     If this part is removed, the size of X and Y are 2N x N.
*
      T_SYM = MPI_WTIME()
      DO I = 1, N
         LAMBDA( N+I ) = -LAMBDA( I )
      END DO
      CALL PDGEADD( 'N', N, N, ONE, X, IX, JX, DESCX, ZERO,
     $     X, IX+N, JX+N, DESCX )
      CALL PDGEADD( 'N', N, N, ONE, X, IX+N, JX, DESCX, ZERO,
     $     X, IX, JX+N, DESCX )
      CALL PDGEADD( 'N', N, N, ONE, Y, IY, JY, DESCY, ZERO,
     $     Y, IY+N, JY+N, DESCY )
      CALL PDGEADD( 'N', N, N, ONE, Y, IY+N, JY, DESCY, ZERO,
     $     Y, IY, JY+N, DESCY )
      T_SYM = MPI_WTIME() - T_SYM
!      IF ( MYROW+MYCOL .EQ. 0 )
!     $   WRITE( *, * ) 't_sym = ', T_SYM, ';'
*
      WORK( 1 ) = DBLE( LWKOPT )
*
      RETURN
*
*     End of PDBSSOLVER1_SVD().
*
      END
#endif
