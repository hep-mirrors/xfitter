      FUNCTION ZBRLAM (WLLN)
C                                                   -=-=- zbrlam
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / TRNCOM / VMULM, JRDR, N, N1
 
      WMULM = EXP (WLLN)
      TEM1 = 1./ ALPQCD(JRDR, N1, WMULM, I)
      TEM2 = 1./ ALPQCD(JRDR, N,  VMULM, I)
 
      ZBRLAM = TEM1 - TEM2
 
      END
 
 
C************************************************************
 
