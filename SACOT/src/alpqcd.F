      FUNCTION ALPQCD (IRDR, NF, RML, IRT)
C                                                   -=-=- alpqcd
 
C                                 Returns the QCD alpha/pi for RML = MU / LAMDA
C                                 using the standard perturbative formula for
C                                 NF flavors and to IRDR th order in 1/LOG(RML)
 
C                                 Return Code:  IRT
C                                                0:   O.K.
C                                                1:   Mu < Lamda; returns 99.
C                                                2:   Alpha > 10 ; be careful!
C                                                3:   IRDR out of range
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
      PARAMETER (CG = 3.0, TR = 0.5, CF = 4.0/3.0)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT
 
      DATA IW1, Iw2 / 0, 0 /
 
      IRT = 0
 
      IF (IRDR .LT. 1 .OR. IRDR .GT. 3) THEN
        WRITE(NOUT, *)
     >  'Order parameter out of range in ALPQCD; IRDR = ', IRDR
        IRT = 3
        STOP
      ENDIF
 
      B0 = (11.* CG  - 2.* NF) / 3.
      B1 = (34.* CG**2 - 10.* CG * NF - 6.* CF * NF) / 3.
      B2 = 2857./2. -5033./18.*NF + 325./54.*NF**2
      RM2 = RML ** 2

C           AlpQcd is used mainly as a mathematical function, for
C            inversion, as well as evaluation of the physical Alpi.
C           Warning should be deferred to the calling program.
 
      IF (RM2 .LE. 1.) THEN
         IRT = 1
C         CALL QWARN(IW1, NWRT,
C     >    'RM2 (=MU/LAMDA) < 1. not allowed in ALPQCD; Alpha->99.',
C     >    'RM2', RM2, D1, BIG, 1)
         ALPQCD = 99
         RETURN
      ENDIF
 
      ALN = LOG (RM2)
      AL = 4./ B0 / ALN
      IF (IRDR .EQ. 2) THEN
        AL = AL * (1.- B1 * LOG(ALN) / ALN / B0**2)

      ELSEIF (IRDR .EQ. 3) THEN 
        AL = AL * (1.- B1 * LOG(ALN) / ALN / B0**2 
     >    + (B1**2*(LOG(ALN)**2 -LOG(ALN) -1.)+B0*B2)/ ALN**2 / B0**4)
      
      ELSEIF (AL .GE. 3.) THEN
         IRT = 2
C         CALL QWARN(IW2, NWRT, 'ALPQCD > 3. in ALPQCD', 'ALPQCD', AL,
C     >              D0, D1, 1)
      ENDIF
 
      ALPQCD = AL
 
      RETURN
C                       *********************
      END
C
