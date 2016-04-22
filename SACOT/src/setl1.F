      SUBROUTINE SETL1  (NEF, VLAM)
C                                                   -=-=- setl1
C     Given LAMDA = VLAM for NEF flavors:
C                    (i) fills the array  ALAM (0:NF) with effective LAMDA's;
C                    (ii) fills the array AMHAT (NF) with threshold masses;
C                    (iii) count the # of "heavy quarks" (QMS > EFFLAM);
C                    (iv) fix the parameter AMN defined as MAX (ALAM),
C                         times safety factor;
C                    (v) set AL in / QCDPAR / equal to ALAM (NF);
C       Uses formula with NORDER (1 or 2) -- see /QCDPAR/

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      LOGICAL ALPHASET, AmHatSet
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
      COMMON / QCDPAR_CT / FSWITCH, AL, NF, NORDER, ALPHASET
      COMMON / IOUNIT / NIN, NOUT, NWRT

C                    Use NLO prescriptions for setting up mass thresholds for simplicity.
C                    Nrder is a local variable.  Norder remains unchanged.
      Nrder = Min(Norder, 2)

      IF (NEF .LT. 0 .OR. NEF .GT. NF) THEN
        WRITE(NOUT,*)'NEF out of range in SETL1, NEF, NF =',NEF,NF
        STOP
      ENDIF

      call SetAmhatf            !pn10 set switching points given on-shell masses
      ALAM(NEF) = VLAM
      DO 10 N = NEF, 1, -1
        CALL TRNLAM(Nrder, N, -1, IR1)
 10   CONTINUE
      DO 20 N = NEF, NF-1
        CALL TRNLAM(Nrder, N, 1, IR1)
 20   CONTINUE
cpn10 Distinguish between light and heavy flavors more precisely 
c     than in SetAmhatf based on the Lambda_QCD values                        
      DO 30, N = NF, 1, -1
        IF ((ALAM(N) .GE. 0.7 * AMHAT(N))
     >    .OR. (ALAM(N-1) .GE. 0.7 * AMHAT(N)))THEN
          NHQ = NF - N
          GOTO 40
        ENDIF
 30   CONTINUE
      NHQ = NF
 40   CONTINUE
      DO 50, N = NF-NHQ, 1, -1
        AMHAT(N) = 0
        ALAM(N-1) = ALAM(N)
 50   CONTINUE
C========================               Find minimum mu
      AMN = ALAM(NF)
      DO 60, N = 0, NF-1
        IF (ALAM(N) .GT. AMN)  AMN = ALAM(N)
 60   CONTINUE
      AMN = AMN * 1.0001
      AL = ALAM(NF)

      ALPHASET=.TRUE.
      RETURN
C**************************************************************
      END

C                                                          =-=-= Setaux
