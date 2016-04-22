      FUNCTION ALEPI (AMU, NEF)
C                                                   -=-=- alepi
C                   Returns ALPHA/PI using the Effective Lamda appropriate for
C                             NEF flavors without regard to the value of AMU.
C                   Appropriate for Renormalization Schemes with fixed NEF.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      LOGICAL AlphaSet, AmHatSet
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
      COMMON /QCDPAR_CT/ FSWITCH, AL, NF, NORDER, AlphaSet
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
 
      DATA IW1, IW2 / 2*0 /
C
cpn11      IF (.NOT.ALPHASET) CALL LAMCWZ
      if (.not.alphaset) stop 'ALEPI: Lambdas are not set'
 
      ALM = ALAM(NEF)
      ALEPI = ALPQCD (NORDER, NEF, AMU/ALM, IRT)
 
      IF     (IRT .EQ. 1) THEN
         CALL QWARN (IW1, NWRT, 'AMU < ALAM in ALEPI', 'MU', AMU,
     >             ALM, BIG, 1)
         WRITE (NWRT, '(A,I4,F15.3)') 'NEFF, LAMDA = ', NEF, ALM
      ELSEIF (IRT .EQ. 2) THEN
         CALL QWARN (IW2, NWRT, 'ALPI > 3; Be aware!', 'ALEPI', ALEPI,
     >             D0, D1, 0)
         WRITE (NWRT, '(A,I4,2F15.3)') 'NF, LAM, MU= ', NEF, ALM, AMU
      ENDIF
 
      RETURN
      END
C
C                      **************************************
C
