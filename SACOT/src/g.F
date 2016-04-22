      FUNCTION G (AMU)
C                                                   -=-=- g
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C               Returns effective coupling.
      PARAMETER (PI = 3.1415927)
      G=2.*PI*SQRT(ALPI(AMU))
      RETURN
C                      ******************
      END
C
