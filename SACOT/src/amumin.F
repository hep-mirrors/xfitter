      FUNCTION AMUMIN()
C                                                   -=-=- amumin
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Returns the minimum mu allowed by perturbative theory.
      logical AlphaSet, AmHatSet
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
      COMMON /QCDPAR_CT/ FSWITCH, AL, NFL, NORDER, ALPHASET

cpn11      IF (.NOT.ALPHASET) CALL LAMCWZ
      if (.not.alphaset) stop 'AMUMIN: Lambdas are not set'

      AMUMIN = AMN
      RETURN
      END
C
C************************************************************************
C
