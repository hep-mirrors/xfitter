      SUBROUTINE SETLAM (NEF, WLAM, IRDR)
C                                                   -=-=- setlam
C     The values of LAMBDA=WLAM with NEF effective flavors is given. The
C     coupling is assumed to be given by the IRDR formula.
C     First lambda is converted to a value that gives approximately the
C     same value of alpha_s when the formula with the current value of
C     NORDER is used.  Then SETL1 is called to update the rest of the
C     internal arrays.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / QCDPAR_CT / FSWITCH, AL, NF, NORDER, ALPHASET
      LOGICAL ALPHASET
 
      IF ((NEF .LT. 0) .OR. (NEF .GT. NF)) THEN
         WRITE(NOUT,*)'NEF out of range in SETLAM, NEF, NF=', NEF,NF
         STOP
      ENDIF
C                         Adjust Lamda value if ORDER parameters don't match.
C                         NORDER is not changed.
      VLAM = WLAM
      IF (IRDR .NE. NORDER) CALL CNVL1 (IRDR, NORDER, NEF, VLAM)
      CALL SETL1 (NEF, VLAM)
      END
C
C************************************************************
C
