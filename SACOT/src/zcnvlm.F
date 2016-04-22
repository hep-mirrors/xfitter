      FUNCTION ZCNVLM (VLAM)
C                                                   -=-=- zcnvlm
C  Auxiliary function for CNVL1, which solves ZCNVLAM=0.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      COMMON / LAMCNV / AMU, ULAM, NFL, IRD, JRD
 
      ZCNVLM= ALPQCD (IRD,NFL,AMU/ULAM,I) - ALPQCD (JRD,NFL,AMU/VLAM,I)
 
      END
C
C**************************************************************
C
