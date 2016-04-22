      FUNCTION RTALF (EFLLN)
C                                                   -=-=- rtalf
C  Auxiliary function for ALFSET, which solves equation RTALF=0.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (PI = 3.141592653589)
      COMMON / RTALFC / ALFST, JORD, NEFF
 
      EFMULM = EXP (EFLLN)
      TEM1 = PI / ALFST
      TEM2 = 1. / ALPQCD (JORD, NEFF, EFMULM, I)
 
      RTALF = TEM1 - TEM2
C************************************************************ 
      END

C                                                          =-=-= Qcdpar
