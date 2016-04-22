      FUNCTION CH(I)
C                                                   -=-=- ch
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                     Returns charge of parton I.
C                        I=0 is gluon.
C                        I>0 is quark.
C                        I<0 is antiquark.
C                     See  BLOCK DATA for code.
C
      COMMON /QCDPAR_CT/ FSWITCH, AL, NF, NORDER, ALPHASET
      COMMON /COMQCH/ VALCH(9)
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      LOGICAL ALPHASET
      CH=0.
      IF (IABS(I).GT.NF) THEN
          WRITE (NOUT, *) 'I IS OUT OF RANGE IN FUNCTION CH'
      ELSEIF (I.NE.0)  THEN
          CH = VALCH(IABS(I))
          IF (I.LT.0)  CH = -CH
      ENDIF
      RETURN
      END
C
C**********************************************************************
C
C                                                          =-=-= Parqcd
