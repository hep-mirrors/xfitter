      SUBROUTINE LAMCWZ
C                                                   -=-=- lamcwz
C                       Set /CWZPRM/ from /QCDPAR/
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /QCDPAR_CT/ FSWITCH, AL, NF, NORDER, ALPHASET
      LOGICAL ALPHASET
      CALL SETL1 (NF, AL)
      END
C
C***********************
C
