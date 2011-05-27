      LOGICAL FUNCTION INTRAC(DUMMY)
C==============================================================================
C     A bypass version of the real INTRAC(DUMMY) of CERN LIB
C     -- to set INTRAC = .FALSE.
C==============================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTRAC = .FALSE.
      RETURN
      END
