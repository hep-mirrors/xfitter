C ================================================================
      function alQCD(Q,iset)
C========================================================================
C     *** THIS IS A DUMMY FUNCTION **********
C     This links into a dummy routine for the Cteq4 alphas values
C========================================================================
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

C      alQCD = ascteq4(Q,iset)
      r2=q**2  !*** SCALE IS Q^2
C     call getalf(alfs,r2)  !**** NOT THE RIGHT FUNCTION

      alfs= hf_get_alphas(r2)

      alQCD=alfs
cv      print*,'alphas',alqcd
      RETURN
      END
