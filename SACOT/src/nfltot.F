      FUNCTION NFLTOT()
C                                                   -=-=- nfltot
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C          Returns the total number of flavors.
      COMMON / QCDPAR_CT / FSWITCH, AL, NF, NORDER, ALPHASET
      LOGICAL ALPHASET
      NFLTOT=NF
      RETURN
      END
C
C***********************************************************
C
