C----------------------------------------------------------------------
C
C> @file deqn_lapack.f
C> @brief LAPACK-based linear solver to replace CERNLIB DEQN
C>
C> This provides a drop-in replacement for DEQN that uses LAPACK's DGESV
C> instead of CERNLIB's custom code.
C> DGESV internally uses level-3 BLAS (DGEMM, DTRSM) which can benefit
C> from multi-threading, but in my testing, leads to slowdowns. 
C>
C> For now, it operates in single-threaded mode.
C> The relevant setting is OPENBLAS_NUM_THREADS.
C>
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C> @brief Solve A*X=B using LAPACK DGESV
C>
C> @param N      Order of matrix A (number of equations)
C> @param A      On entry: NxN coefficient matrix (stored column-major)
C>               On exit: L and U factors from the factorization A = P*L*U
C> @param IDIM   Leading dimension of A
C> @param R      Workspace (unused, kept for DEQN compatibility)
C> @param IFAIL  Output: 0 = success, -1 = error
C> @param K      Number of right-hand sides
C> @param B      On entry: NxK right-hand side matrix B
C>               On exit: NxK solution matrix X
C>
C> @note This is a drop-in replacement for CERNLIB DEQN
C----------------------------------------------------------------------
      SUBROUTINE DEQN_LAPACK(N, A, IDIM, R, IFAIL, K, B)
      
      IMPLICIT NONE
      
C Arguments
      INTEGER N, IDIM, IFAIL, K
      DOUBLE PRECISION A(IDIM, N)
      DOUBLE PRECISION B(IDIM, K)
      REAL R(N)  ! Unused, kept for DEQN interface compatibility
      
C Local variables
      INTEGER IPIV(N)   ! Pivot indices for LAPACK
      INTEGER INFO
      
C External LAPACK routine
      EXTERNAL DGESV
      
C Initialize
      IFAIL = 0
      
C Parameter validation
      IF (N .LT. 1 .OR. N .GT. IDIM .OR. K .LT. 1) THEN
         CALL hf_errlog(14072406, 'E: DEQN_LAPACK PARAMETER ERROR')
         IFAIL = -1
         RETURN
      ENDIF
      
C Call LAPACK DGESV to solve A*X = B
      CALL DGESV(N, K, A, IDIM, IPIV, B, IDIM, INFO)
      
      IF (INFO .NE. 0) THEN
         IF (INFO .LT. 0) THEN
            CALL hf_errlog(14072407, 
     $           'E: DEQN_LAPACK: Illegal argument to DGESV')
         ELSE
            CALL hf_errlog(14072408,
     $           'E: DEQN_LAPACK: Singular matrix detected')
         ENDIF
         IFAIL = -1
      ENDIF
      
      RETURN
      END

  
C----------------------------------------------------------------------
C> @brief Wrapper to switch between DEQN and DEQN_LAPACK
C>
C> Which solver is used is controlled by the OpenMP.useLapackSolver
C> flag in parameters.yaml, read once at startup by init_lapack_solver().
C>
C> @param N      Order of matrix A
C> @param A      Coefficient matrix
C> @param IDIM   Leading dimension
C> @param R      Workspace
C> @param IFAIL  Error flag
C> @param K      Number of RHS
C> @param B      RHS on input, solution on output
C----------------------------------------------------------------------
      SUBROUTINE DEQN_AUTO(N, A, IDIM, R, IFAIL, K, B)

      IMPLICIT NONE

      INTEGER N, IDIM, IFAIL, K
      DOUBLE PRECISION A(IDIM, N)
      DOUBLE PRECISION B(IDIM, K)
      REAL R(N)

      LOGICAL use_lapack
      COMMON /lapack_settings/ use_lapack

      IF (use_lapack) THEN
         CALL DEQN_LAPACK(N, A, IDIM, R, IFAIL, K, B)
      ELSE
         CALL DEQN(N, A, IDIM, R, IFAIL, K, B)
      ENDIF

      RETURN
      END


C----------------------------------------------------------------------
C> @brief Initialize LAPACK solver settings from parameters.yaml
C>
C> Reads OpenMP.useLapackSolver via GetParamI (C shim getparami_). The
C> value is stashed into gParametersI["UseLapackSolver"] at startup by
C> apply_openmp_settings() in xfitter_pars.cc.
C----------------------------------------------------------------------
      SUBROUTINE init_lapack_solver()

      IMPLICIT NONE

      LOGICAL use_lapack
      COMMON /lapack_settings/ use_lapack

      INTEGER GetParamI
      EXTERNAL GetParamI

      LOGICAL lfirst
      DATA lfirst /.true./
      SAVE lfirst

      IF (.NOT. lfirst) RETURN
      lfirst = .false.

      use_lapack = (GetParamI('UseLapackSolver') .NE. 0)

      IF (use_lapack) THEN
         PRINT '(A)', ' Matrix solver: LAPACK DGESV'
      ELSE
         PRINT '(A)', ' Matrix solver: CERNLIB DEQN (scalar code)'
      ENDIF

      RETURN
      END
