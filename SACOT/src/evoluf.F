      SUBROUTINE EVOLUF(FQ2,FQ1,Q2,Q1,N)
C                                                   -=-=- evoluf
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                        ----------------------------
C?????????????????????????TO CHECK ????????????????????????????
C                        ----------------------------
C            Returns moments of parton distribution function at Q2 from Q1.
C            by integrating
C                D(F(Q,N,I))/D(LN(Q**2))=-GAMA(I,J)*F(Q,N,J)
C            where N is rank of the moments.
C            GAMA is aN NF X NF  matrix.
C            F(I=1) corresponds to gluon.
C            F(I=2 to NF+1) corresponds to NF quarks.
C
      logical AlphaSet, AmHatSet      
      COMMON /QCDPAR_CT/ FSWITCH, AL, NF, NORDER, ALPHASET
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      DIMENSION FQ1(11),FQ2(11),U(11,11),UINV(11,11),
     >          GAMAD(11,11),TEMP(11,11)
      DATA U,UINV,GAMAD/363*0./
C
cpn11      IF (.NOT.ALPHASET) CALL LAMCWZ
      if (.not.alphaset) stop 'EVOLUF: Lambdas are not set'

      IF((Q1.LE.AMN).OR.(Q2.LE.AMN)) GOTO 300
      IF(N.LT.2) GOTO 600
C
      NMMIN=NF-NFL(Q2)
      NMMAX=NF-NFL(Q1)
C
      T2=2.* LOG(Q1)
      B1=FLOAT(33-2*(NF-NMMAX-1))/4.
      SUM=0.
      DO 10 I=2,N
         SUM=SUM+1./FLOAT(I)
 10   CONTINUE
      GAGG=.75-9./FLOAT((N+1)*(N+2))-9./FLOAT(N*(N-1))
     >        +FLOAT(NF-NMMAX-1)/2.+SUM*9.
      GAFG=3./FLOAT(N)-6./FLOAT((N+1)*(N+2))
      GAGF=1./FLOAT(N+1)+2./FLOAT(N*(N-1))
      GAFF=1.-2./FLOAT(N*(N+1))+4.*SUM
      FQ2(1)=FQ1(1)
      DO 20 I=2,NF+1
         DO 30 J=3,NF+1
            UINV(J,I)=1./FLOAT(NF)
  30     CONTINUE
         FQ2(I)=FQ1(I)
         UINV(I,I)=UINV(I,I)-1.
         U(I,I)=-1.
         U(2,I)=1.
         U(I,1)=1.
         U(I,2)=1.
 20   CONTINUE
C
      DO 100 K=NMMAX+1,NMMIN+1,-1
         T1=T2
         B1=B1-.5
         GAGG=GAGG+.5
         IF(K.EQ.(NMMIN+1)) THEN
            T2=2.* LOG(Q2)
         ELSE
            T2=2.* LOG(AMHAT(NF-K+2))
         ENDIF
         A=SQRT((GAGG-GAFF)**2+4.*NF*GAGF*GAFG)
         U(1,1)=.5*(GAGG-GAFF+A)
         U(1,2)=.5*(GAGG-GAFF-A)
         DO 40 L=2,NF+1
            UINV(1,L)=-U(1,2)/A/FLOAT(NF)
            UINV(2,L)=U(1,1)/A/FLOAT(NF)
  40     CONTINUE
         UINV(1,1)=1./A
         UINV(1,2)=-1./A
         D=2* LOG(ALAM(NF+1-K))
         C=(T2-D)/(T1-D)
         DO 50 M=3,NF+1
            GAMAD(M,M)=C**(-GAFF/B1)
  50     CONTINUE
         GAMAD(1,1)=C**(-(U(1,1)+GAFF)/B1)
         GAMAD(2,2)=C**(-(U(1,2)+GAFF)/B1)
         CALL MTMUL(NF+1,NF+1,NF+1,U,GAMAD,TEMP)
         CALL MTMUL(NF+1,NF+1,NF+1,TEMP,UINV,TEMP)
         CALL MTMUL(NF+1,NF+1,1,TEMP,FQ2,FQ2)
 100  CONTINUE
      RETURN
C
300   WRITE (NOUT,400)
400   FORMAT('Q1 OR/AND Q2 IS TOO SMALL IN EVOLUF')
600   DO 500 I=1,11
500      FQ2(I)=0.
      RETURN
      END
C
C                             *************************
