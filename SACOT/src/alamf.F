      FUNCTION ALAMF(N)
C                                                   -=-=- alamf
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                 Returns the value of LambdaCWZ in the energy range
C                 with N "light" quarks.
 
      LOGICAL ALPHASET, AmHatSet
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
      COMMON /QCDPAR_CT/ FSWITCH, AL, NF, NORDER, ALPHASET
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
cpn11      IF (.NOT.ALPHASET) CALL LAMCWZ
      if (.not.alphaset) stop 'ALAMF: Lambdas are not set'

      IF ((N.LT.0) .OR. (N.GT.9)) THEN
         WRITE (NOUT, *) ' N IS OUT OF RANGE IN ALAMF'
         ALAMF=0.
      ELSE
         ALAMF = ALAM(MAX(N, NF-NHQ))
      ENDIF
      RETURN
      END
C
C***************************************************************
C
