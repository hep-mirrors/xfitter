      FUNCTION ALPIOR (AMU, NL)
C                                                   -=-=- alpior
C               Returns effective g**2/(4pi**2) = alpha/pi.
C               Use formula with NL loops for beta function.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      logical alphaset, AmHatSet
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
      COMMON /QCDPAR_CT/ FSWITCH, AL, NF, NORDER, ALPHASET
 
C                      Use the following as subroutine argument of type
C                      set by IMPLICIT statement:
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
 
      DATA IW1, IW2 / 2*0 /
 
cpn11      IF (.NOT.ALPHASET) CALL LAMCWZ
      if (.not.alphaset) stop 'ALPIOR: Lambdas are not set'

      NEFF = NFL(AMU)
      ALM  = ALAM(NEFF)
      ALPIOR = ALPQCD (NL, NEFF, AMU/ALM, IRT)
 
      IF (IRT .EQ. 1) THEN
         CALL QWARN (IW1, NWRT, 'AMU < ALAM in ALPIOR', 'MU', AMU,
     >              ALM, BIG, 1)
         WRITE (NWRT, '(A,I4,F15.3)') 'NEFF, LAMDA = ', NEFF, ALM
      ELSEIF (IRT .EQ. 2) THEN
         CALL QWARN (IW2,NWRT,'ALPIOR > 3; Be aware!','ALPIOR',ALPIOR,
     >             D0, D1, 0)
         WRITE (NWRT, '(A,I4,2F15.3)') 'NF, LAM, MU= ', NEFF, ALM, AMU
      ENDIF
      END
C
C*****************************************************************
C
