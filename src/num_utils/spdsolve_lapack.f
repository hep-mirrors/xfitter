C----------------------------------------------------------------------
C
C> @file spdsolve_lapack.f
C> @brief Cholesky-based solvers for the (symmetric pos-def)
C>        nuisance-parameter system A*x = C.
C>
C> The nuisance matrix A = P + Gamma V^-1 Gamma^T is symmetric and, for
C> P positive, positive definite. DPOTRF+DPOTRS (Cholesky)
C> should be quicker than DGESV. If DPOTRF fails (possible when some
C> SysPriorScale is zero), we fall back to DGESV with a
C> warning; with OpenMP.useLapackSolver=false the CERNLIB routines are
C> used as before.
C>
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C> @brief Solve A*x=C for SPD A: Cholesky with LU fallback
C>
C> @param N      Order of matrix A
C> @param A      On entry: NxN symmetric matrix (both triangles set)
C>               On exit: overwritten by its factorization
C> @param IDIM   Leading dimension of A
C> @param IFAIL  Output: 0 = success, -1 = error
C> @param B      On entry: right-hand side; on exit: solution x
C----------------------------------------------------------------------
      SUBROUTINE SPDSOLVE_LAPACK(N, A, IDIM, IFAIL, B)

      IMPLICIT NONE

      INTEGER N, IDIM, IFAIL
      DOUBLE PRECISION A(IDIM, N)
      DOUBLE PRECISION B(N)

      INTEGER INFO
      INTEGER IPIV(N)

      EXTERNAL DPOTRF, DPOTRS, DGESV

      IFAIL = 0

      IF (N .LT. 1 .OR. N .GT. IDIM) THEN
         CALL hf_errlog(14072420, 'E: SPDSOLVE_LAPACK PARAMETER ERROR')
         IFAIL = -1
         RETURN
      ENDIF

      CALL DPOTRF('L', N, A, IDIM, INFO)

      IF (INFO .EQ. 0) THEN
         CALL DPOTRS('L', N, 1, A, IDIM, B, N, INFO)
         IF (INFO .NE. 0) THEN
            CALL hf_errlog(14072421,
     $           'E: SPDSOLVE_LAPACK: DPOTRS failed')
            IFAIL = -1
         ENDIF
      ELSEIF (INFO .GT. 0) THEN
C Not positive definite : general solve instead.
C DPOTRF('L') only touched the lower triangle; restore it by symmetry.
         CALL hf_errlog(14072422,
     $        'W: SPDSOLVE_LAPACK: nuisance matrix not positive'//
     $        ' definite, falling back to DGESV')
         CALL RESTORE_LOWER_FROM_UPPER(N, A, IDIM)
         CALL DGESV(N, 1, A, IDIM, IPIV, B, N, INFO)
         IF (INFO .NE. 0) THEN
            CALL hf_errlog(14072423,
     $           'E: SPDSOLVE_LAPACK: DGESV failed, singular matrix')
            IFAIL = -1
         ENDIF
      ELSE
         CALL hf_errlog(14072424,
     $        'E: SPDSOLVE_LAPACK: Illegal argument to DPOTRF')
         IFAIL = -1
      ENDIF

      RETURN
      END


C----------------------------------------------------------------------
C> @brief Copy the (intact) strict upper triangle onto the lower one
C----------------------------------------------------------------------
      SUBROUTINE RESTORE_LOWER_FROM_UPPER(N, A, IDIM)
      IMPLICIT NONE
      INTEGER N, IDIM, I, J
      DOUBLE PRECISION A(IDIM, N)
      DO J = 1, N
         DO I = J+1, N
            A(I,J) = A(J,I)
         ENDDO
      ENDDO
      RETURN
      END


C----------------------------------------------------------------------
C> @brief Wrapper: SPD LAPACK solve or CERNLIB DEQN
C>
C> Same calling sequence as DEQN_AUTO/DEQN so it can replace them at the
C> nuisance-solve sites. Dispatches on OpenMP.useLapackSolver via the
C> lapack_settings common (set by init_lapack_solver()).
C>
C> @param N      Order of matrix A
C> @param A      Matrix (does not survive)
C> @param IDIM   Leading dimension
C> @param R      Workspace (CERNLIB DEQN only)
C> @param IFAIL  Error flag
C> @param K      Number of RHS (must be 1 for the LAPACK path here)
C> @param B      RHS on input, solution on output
C----------------------------------------------------------------------
      SUBROUTINE SPDSOLVE_AUTO(N, A, IDIM, R, IFAIL, K, B)

      IMPLICIT NONE

      INTEGER N, IDIM, IFAIL, K
      DOUBLE PRECISION A(IDIM, N)
      DOUBLE PRECISION B(IDIM, K)
      REAL R(N)

      LOGICAL use_lapack
      COMMON /lapack_settings/ use_lapack

      IF (use_lapack .AND. K .EQ. 1) THEN
         CALL SPDSOLVE_LAPACK(N, A, IDIM, IFAIL, B)
      ELSE
         CALL DEQN(N, A, IDIM, R, IFAIL, K, B)
      ENDIF

      RETURN
      END


C----------------------------------------------------------------------
C> @brief Wrapper for CERNLIB DEQINV: solve A*x=C and return A^-1 in A
C>
C> Used once per fit (iflag=3) where both the solution and the full
C> inverse (for nuisance errors and correlations) are required.
C> LAPACK path: invert via DINV_LAPACK,
C> then x = A^-1 * C by DSYMV.
C----------------------------------------------------------------------
      SUBROUTINE DEQINV_AUTO(N, A, IDIM, R, IFAIL, K, B)

      IMPLICIT NONE

      INTEGER N, IDIM, IFAIL, K
      DOUBLE PRECISION A(IDIM, N)
      DOUBLE PRECISION B(IDIM, K)
      REAL R(N)

      DOUBLE PRECISION BSAV(N)
      INTEGER I

      LOGICAL use_lapack
      COMMON /lapack_settings/ use_lapack

      EXTERNAL DSYMV

      IF (use_lapack .AND. K .EQ. 1) THEN
         DO I = 1, N
            BSAV(I) = B(I,1)
         ENDDO
         CALL DINV_LAPACK(N, A, IDIM, IFAIL)
         IF (IFAIL .EQ. 0) THEN
            CALL DSYMV('L', N, 1.0D0, A, IDIM, BSAV, 1, 0.0D0, B, 1)
         ENDIF
      ELSE
         CALL DEQINV(N, A, IDIM, R, IFAIL, K, B)
      ENDIF

      RETURN
      END
