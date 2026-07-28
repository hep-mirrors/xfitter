C----------------------------------------------------------------------
C
C> @file dinv_lapack.f
C> @brief LAPACK-based matrix inversion to replace CERNLIB DINV similar
C> to DEQN. 
C>
C> This provides a drop-in replacement for DINV that uses LAPACK's
C> Cholesky factorization (DPOTRF/DPOTRI) instead of CERNLIB's scalar
C> LU code (DFACT/DFINV). This is written for pos-def symmetric
C> convariance matrices in particular. This
C> uses level-3 BLAS, so it can benefit from
C> multi-threading (OpenMP.blasThreads in parameters.yaml).
C>
C> If the matrix is not positive definite (e.g. a bad 
C> correlation matrix), it falls back to the symmetric indefinite
C> factorization DSYTRF/DSYTRI with a warning.
C>
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C> @brief Invert a symmetric matrix in place using LAPACK
C>
C> @param N      Order of matrix A
C> @param A      On entry: NxN symmetric matrix (full storage)
C>               On exit: its inverse (full storage, both triangles)
C> @param IDIM   Leading dimension of A
C> @param IFAIL  Output: 0 = success, -1 = error
C>
C> @note Drop-in replacement for CERNLIB DINV (which returns the full
C>       inverse); the R workspace argument of DINV is not needed.
C----------------------------------------------------------------------
      SUBROUTINE DINV_LAPACK(N, A, IDIM, IFAIL)

      IMPLICIT NONE

C Arguments
      INTEGER N, IDIM, IFAIL
      DOUBLE PRECISION A(IDIM, N)

C Local variables
      INTEGER INFO, LWORK, I, J
      INTEGER IPIV(N)
      DOUBLE PRECISION WQUERY(1)
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)

C External LAPACK routines
      EXTERNAL DPOTRF, DPOTRI, DSYTRF, DSYTRI

C Initialize
      IFAIL = 0

C Parameter validation
      IF (N .LT. 1 .OR. N .GT. IDIM) THEN
         CALL hf_errlog(14072410, 'E: DINV_LAPACK PARAMETER ERROR')
         IFAIL = -1
         RETURN
      ENDIF

C Try Cholesky first: covariance matrices should be positive definite
      CALL DPOTRF('U', N, A, IDIM, INFO)

      IF (INFO .EQ. 0) THEN
         CALL DPOTRI('U', N, A, IDIM, INFO)
         IF (INFO .NE. 0) THEN
            CALL hf_errlog(14072411,
     $           'E: DINV_LAPACK: DPOTRI failed, singular matrix')
            IFAIL = -1
            RETURN
         ENDIF
      ELSEIF (INFO .GT. 0) THEN
C Not positive definite: warn and retry with symmetric indefinite
C factorization. A was overwritten by the partial factor, so the
C caller's matrix must be restored by symmetry: DPOTRF('U') only
C touches the upper triangle, the lower one is still intact.
         CALL hf_errlog(14072412,
     $        'W: DINV_LAPACK: covariance matrix not positive'//
     $        ' definite, using symmetric indefinite inversion')
         DO J = 1, N
            DO I = 1, J-1
               A(I,J) = A(J,I)
            ENDDO
         ENDDO
C Workspace query for DSYTRF
         CALL DSYTRF('U', N, A, IDIM, IPIV, WQUERY, -1, INFO)
         LWORK = MAX(INT(WQUERY(1)), N)
         ALLOCATE(WORK(LWORK))
         CALL DSYTRF('U', N, A, IDIM, IPIV, WORK, LWORK, INFO)
         IF (INFO .EQ. 0) THEN
            CALL DSYTRI('U', N, A, IDIM, IPIV, WORK, INFO)
         ENDIF
         DEALLOCATE(WORK)
         IF (INFO .NE. 0) THEN
            CALL hf_errlog(14072413,
     $           'E: DINV_LAPACK: DSYTRF/DSYTRI failed,'//
     $           ' singular matrix')
            IFAIL = -1
            RETURN
         ENDIF
      ELSE
         CALL hf_errlog(14072414,
     $        'E: DINV_LAPACK: Illegal argument to DPOTRF')
         IFAIL = -1
         RETURN
      ENDIF

C Both routines return only the upper triangle; fill the lower one
C (callers read the full inverse, as CERNLIB DINV provided)
      DO J = 1, N
         DO I = 1, J-1
            A(J,I) = A(I,J)
         ENDDO
      ENDDO

      RETURN
      END


C----------------------------------------------------------------------
C> @brief Wrapper to switch between DINV and DINV_LAPACK
C>
C> Which routine is used is controlled by the OpenMP.useLapackSolver
C> flag in parameters.yaml, read once at startup by init_lapack_solver()
C> (shared with DEQN_AUTO via the lapack_settings common block).
C>
C> @param N      Order of matrix A
C> @param A      Matrix on input, inverse on output
C> @param IDIM   Leading dimension
C> @param R      Workspace (used by CERNLIB DINV only)
C> @param IFAIL  Error flag
C----------------------------------------------------------------------
      SUBROUTINE DINV_AUTO(N, A, IDIM, R, IFAIL)

      IMPLICIT NONE

      INTEGER N, IDIM, IFAIL
      DOUBLE PRECISION A(IDIM, N)
      REAL R(N)

      LOGICAL use_lapack
      COMMON /lapack_settings/ use_lapack

      IF (use_lapack) THEN
         CALL DINV_LAPACK(N, A, IDIM, IFAIL)
      ELSE
         CALL DINV(N, A, IDIM, R, IFAIL)
      ENDIF

      RETURN
      END
