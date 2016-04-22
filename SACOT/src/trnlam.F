      SUBROUTINE TRNLAM (IRDR, NF, IACT, IRT)
C                                                   -=-=- trnlam
 
C     This routine transforms LAMDA (N) to LAMDA (N+IACT) where IACT = 1/-1
C     The transformation is obtained by requiring the coupling constant to
C                be continuous at the scale Mu = Mass of the (N+1)th quark.
 
C                                         IRT is an return code.
C                                            (0 for OK)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      logical AmHatSet
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
      COMMON / TRNCOM / VMULM, JRDR, N, N1
 
      EXTERNAL ZBRLAM
 
      DATA ALM0, BLM0, RERR / 0.01, 10.0, 0.0001 /
      DATA IR1, SML / 0, 1.E-5 /
 
      IRT = 0
 
      N = NF
      JRDR = IRDR
      JACT = IACT
      VLAM = ALAM(N)
 
      IF (JACT .GT. 0) THEN
         N1 = N + 1
         THMS = AMHAT(N1)
         ALM = LOG (THMS/VLAM)
         BLM = BLM0
      ELSE
         N1 = N -1
         THMS = AMHAT(N)
         ALM = ALM0
         THMS = MAX (THMS, SML)
         BLM = LOG (THMS/VLAM)
      ENDIF
C                          Fix up for light quark:
      IF (VLAM .GE. 0.7 * THMS) THEN
         IF (JACT . EQ. 1) THEN
            AMHAT(N1) = 0
         ELSE
            AMHAT(N) = 0
         ENDIF
         IRT = 4
         ALAM(N1) = VLAM
         RETURN
      ENDIF
 
C             QZBRNT is the root-finding function to solve ALPHA(N) = ALPHA(N1)
C             Since 1/Alpha is roughly linear in Log(Mu/Lamda), we use the
C             former in ZBRLAM and the latter as the function variable.
      IF (ALM .GE. BLM) THEN
         WRITE (NOUT, *) 'TRNLAM has ALM >= BLM: ', ALM, BLM
         WRITE (NOUT, *) 'I do not know how to continue'
         STOP
         ENDIF
      VMULM = THMS/VLAM
      ERR = RERR * LOG (VMULM)
      WLLN = QZBRNT (ZBRLAM, ALM, BLM, ERR, IR1)
      ALAM(N1) = THMS / EXP (WLLN)
 
      IF (IR1 .NE. 0) THEN
         WRITE (NOUT, *) 'QZBRNT failed to find VLAM in TRNLAM; ',
     >        'NF, VLAM =', NF, VLAM
         WRITE (NOUT, *) 'I do not know how to continue'
        STOP
      ENDIF
      RETURN
      END
C                             *************************
 
