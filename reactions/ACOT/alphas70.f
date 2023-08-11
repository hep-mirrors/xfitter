C ================================================================
      function alQCD(Q,iset)
C========================================================================
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      r2=q**2  !*** SCALE IS Q^2
      alfs= hf_get_alphas(r2)
      alQCD=alfs

      RETURN
      END
