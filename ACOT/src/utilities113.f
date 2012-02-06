C========================================================================
      function alphas(Q,iset)
C========================================================================
C     *** THIS IS A DUMMY FUNCTION **********
C     This links into a dummy routine for the Cteq4 alphas values
C========================================================================
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

C      alphas = ascteq4(Q,iset)
      r2=q*q
      call getalf(alfs,r2)
      alphas=alfs


      RETURN
      END


C**************************************************************
C**************************************************************
C**************************************************************
      FUNCTION TRNGLE (X,Y,Z, IRT)
C   "Triangle Function" of the three sides. 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DATA RDOFF / 1E-11 /

      IRT = 0

      AMX = MAX (X*X, Y*Y, Z*Z) 

      TMP = X*X + Y*Y + Z*Z - 2.* (X*Y + Y*Z + Z*X)
      ATMP= ABS(TMP) 
      
      IF     (ATMP .LT. AMX * RDOFF) THEN
        TMP = ATMP
        IRT = 1
      ELSEIF (TMP .LT. 0.) THEN
        PRINT '(A, 4(1PE12.3))', 'X,Y,Z, TMP =', X, Y, Z, TMP
        CALL HF_ERRLOG(103,'F: Negative argument in'//
     +                     'TRNGLE function. Check for errors!')
      ENDIF

      TRNGLE = SQRT (TMP)
 
      RETURN
C			****************************
      END
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C      PROGRAM MAIN
       SUBROUTINE TESTDILOG()
C-----------------------------------------------------------------------------
C      SIMPLE PROGRAM TO TEST DiLog
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       DOUBLE COMPLEX X1,TEST,DiLog
       EXTERNAL DiLog

        DATA X1 /(-1.0,0.0)/

1      CONTINUE
       WRITE(6,*)          X1
       WRITE(6,*) ' ENTER: X1'
       READ (5,*)          X1


       TEST =DiLog(X1)
       WRITE(6,*) TEST
       GOTO 1

       END
C ***************************************************************************
C ***************************************************************************
C     This subroutine calculates the DiLog function for any complex
C     value of the argument with double precision. 15 significant
C     figures are guaranteed for any ABS(X) < 1.E6
C ***************************************************************************
C ***************************************************************************
      DOUBLE COMPLEX FUNCTION DiLog(X)
      DOUBLE COMPLEX X,XX
      DOUBLE PRECISION BER(10),B,P,PI
      DATA PI /3.14159265358979323846/
      DATA P  /1.6449340668482264365/
      DATA BER /0.166666666666666666667,  -0.033333333333333333333,
     >          0.0238095238095238095238, -0.033333333333333333333,
     >          0.075757575757575757576, -0.253113553113553113553,
     >          1.16666666666666666667,   -7.0921568627450980392,
     >          54.971177944862155388,   -529.12424242424242424/
C ***************************************************************************
             NFLAG=0
             XX=X
             LAND=1
             CALL COUNTRY(XX,NFLAG)

         IF (NFLAG.EQ.0) THEN
             XX=1.-X
             LAND=2
             CALL COUNTRY(XX,NFLAG)
         END IF

         IF (NFLAG.EQ.0) THEN
             XX=1./X
             LAND=3
             CALL COUNTRY(XX,NFLAG)
         END IF

         IF (NFLAG.EQ.0) THEN
             XX=1./(1.-X)
             LAND=4
             CALL COUNTRY(XX,NFLAG)
         END IF

      CALL SERIES(XX,DiLog)

      IF (LAND.EQ.1) RETURN

      IF (LAND.EQ.2) THEN
           DiLog=-DiLog+P
           B=ABS(1.-X)
           IF(B.GT.0.D0) DiLog=DiLog-LOG(X)*LOG(1.-X)
      END IF

      IF (LAND.EQ.3) THEN
          DiLog=-DiLog-P-0.5*LOG(-X)*LOG(-X)
      END IF

      IF (LAND.EQ.4) THEN
          DiLog=DiLog+2.*P-LOG(X)*LOG(1.-X)+0.5*(LOG(X-1.))**2
      END IF

      RETURN
      END

C ***************************************************************************
      SUBROUTINE COUNTRY(X,NFLAG)
      DOUBLE COMPLEX X,Z
      DOUBLE PRECISION R
      Z=1.-X
      IF(ABS(Z).EQ.0.) RETURN
      Z=-LOG(Z)
      R=ABS(Z)
      IF (R.LE.DLOG(2.D0)) NFLAG=1
      RETURN
      END
C ***************************************************************************
      SUBROUTINE SERIES(X,SP)
      DOUBLE COMPLEX SP,TP,X,Z,ZZ
      DOUBLE PRECISION BER(10),FAC,A,P,PI
      DATA PI /3.14159265358979323846/
      DATA P  /1.6449340668482264365/
      DATA BER /0.166666666666666666667,  -0.033333333333333333333,
     >          0.0238095238095238095238, -0.033333333333333333333,
     >          0.075757575757575757576, -0.253113553113553113553,
     >          1.16666666666666666667,   -7.0921568627450980392,
     >          54.971177944862155388,   -529.12424242424242424/

      FAC=0.5
      Z=-LOG(1.-X)
      SP=Z*(1.-0.25*Z)
      ZZ=Z*Z*Z
      DO 1 I=1,10
      A=FAC*BER(I)/DFLOAT(2*I+1)
      TP=SP+A*ZZ
      SP=TP
      ZZ=ZZ*Z*Z
  1   FAC=FAC/DFLOAT(2*(I+1)*(2*I+1))
  2   FORMAT(2E27.18,/,2E27.18)
      RETURN
      END
C**************************************************************
C**************************************************************

C************
C                                 LIBDINT.FOR
C     DEC 02 89
C                        ----------------------------
C                        ****************************
C           MODIFIED VERSION USING NEW ADAPTIVE END-POINT PROCEDURE 

C WKT    AUG 08 91      --  Common blocks re-aligned to avoid warning
C WKT    FEB 22 89      --  based on earlier versions of ADPINT by WKT & JCC.
C                        ----------------------------
C List of Subprograms:

C     FUNCTION   ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZSPL (F, I, IER)
C     SUBROUTINE ADZCAL (F,I)
C     SUBROUTINE SGLINT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOTALZ
C     FUNCTION   INTUSZ ()
C
C     COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), NUMINT,
C    > ICTA, ICTB, FA, FB, IB
C                        ****************************

      FUNCTION ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
 
C     Adaptive integration routine which allows the integrand to be 
C     indeterminant at the lower and/or the upper ends of integration. 

C     Can self-adjust to any integrable singularity at the ends and compute 
C     the closest approximant, hence achieve the required accuracy efficiently
C     (provided the switch(s) IACTA (IACTB) are set to 2).
 
C     Input switches for end-treatment:
C        IACTA = 0 :   Use closed lower-end algorithm 
C                1 :   Open lower-end -- use open quadratic approximant
C                2 :   Open lower-end -- use adaptive singular approximant

C        IACTB = 0, 1, 2   (same as above, for the upper end)
 
C                Integral of F(X) from A to B, with error
C                less than ABS(AERR) + ABS(RERR*INTEGRAL)
C                Best estimate of error returned in ERREST.
C                Error code is IER:               0 :  o.k.
C                1 :  maximum calls to function reached before the 
C                     error criteria are met;
C                2 :  IACTA out of range, set to 1;
C                3 :  IACTB out of range, set to 1.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      DATA SMLL / 1E-20 /
     
      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL) THEN
        CALL HF_ERRLOG(110,'F: ADZINT - Both Aerr and Rerr are zero!')
C       STOP 'Both Aerr and Rerr are zero in ADZINT!'
      ENDIF
  
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call', 
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB
 
      NUMINT = 3
      DX = (B-A)/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A 
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1))
             ELSE 
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2.)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I))
             ELSE
               IB = I
               FB = F(B-DX/2.)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I))
          ENDIF
          CALL ADZCAL(F,I)
   10     CONTINUE
       CALL TOTALZ
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
          IF (ERS .GT. TARGET)  THEN
              NUMOLD = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*NUMOLD .GT. TARGET) CALL ADZSPL(F,I,IER)
   40         CONTINUE
              IF (IER.EQ.0 .AND. NUMINT.NE.NUMOLD)  GOTO 30
              ENDIF
      ADZINT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END


      SUBROUTINE ADZSPL (F, I, IER)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      DATA TINY / 1.D-20 /
     
      IF (NUMINT .GE. MAXINT)  THEN
          IER = 1
          RETURN
          ENDIF
      NUMINT = NUMINT + 1
C                                                         New interval NUMINT
      IF (I .EQ. IB) IB = NUMINT
      U(NUMINT) = (U(I) + V(I)) / 2.
      V(NUMINT) = V(I)
 
      FU(NUMINT) = FW(I)
      FV(NUMINT) = FV(I)
C                                                             New interval I
       V(I) =  U(NUMINT)
      FV(I) = FU(NUMINT)
C                                                    Save old Result and Error
      OLDRES = RESULT(I)
      OLDERR = ERR(I)
     
      CALL ADZCAL (F, I)
      CALL ADZCAL (F, NUMINT)
C                                                               Update result
      DELRES = RESULT(I) + RESULT(NUMINT) - OLDRES
      RES = RES + DELRES
C                                  Good error estimate based on Simpson formula
      GODERR = ABS(DELRES) 
C                                                             Update new global 
      ERS = ERS + GODERR - OLDERR
C                                  Improve local error estimates proportionally
      SUMERR = ERR(I) + ERR(NUMINT)
      IF (SUMERR .GT. TINY) THEN
         FAC = GODERR / SUMERR 
      ELSE
         FAC = 1.
      ENDIF
      
      ERR(I)      = ERR(I) * FAC
      ERR(NUMINT) = ERR(NUMINT) * FAC
 
      RETURN
C                        ****************************
      END
 

      SUBROUTINE ADZCAL (F,I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB
 
      SAVE / ADZWRK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.
     
      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL SGLINT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGLINT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                       Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF
 
      RESULT(I) = TEM         
      ERR   (I) = ABS (ER)
 
      RETURN
C                        ****************************
      END

      SUBROUTINE SGLINT (IACT, F1, F2, F3, DX, FINT, ESTER)

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
C                                                         Use quadratic formula
      TEM = DX * (4.*F1 + 3.*F2 + 2.*F3) / 9.
C                 Error est based on Diff between quadratic and linear integrals
      ER  = DX * (4.*F1 - 6.*F2 + 2.*F3) / 9.

C                          Invoke adaptive singular parametrization if IACT = 2
C                      Algorithm is based on the formula F(x) = AA + BB * x **CC
C                 where AA, BB & CC are determined from F(Dx/4), F(Dx/2) & F(Dx)

      IF (IACT .EQ. 2) THEN
          T1 = F2 - F1
          T2 = F3 - F2
          IF (T1*T2 .LE. 0.) GOTO 7
          T3  = T2 - T1
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -D1)  GOTO 7
          BB  = T1**2 / T3
          AA  = (F1*F3 - F2**2) / T3
C                                          Estimated integral based on A+Bx**C
          TMP = DX * (AA + BB* 4.**CC / (CC + 1.))
C                                       Error estimate based on the difference
          ER = TEM - TMP
C                                              Use the improved integral value
          TEM= TMP 
      ENDIF

    7 FINT = TEM
      ESTER= ER
      RETURN
C                        ****************************
      END
     

      FUNCTION INTUSZ ()
C                    Return number of intervals used in last call to ADZINT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      INTUSZ = NUMINT
      RETURN
C                        ****************************
      END
C

      SUBROUTINE TOTALZ
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END
