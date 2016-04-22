      FUNCTION SUD(Q1,Q2,GARRAY,N)
C                                                   -=-=- sud
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Returns integral from Q1 to Q2 of
C                 (dmu/mu) * gam(g(mu)) * ln(Q2/mu)
C              where N terms in anomalous dim. gam(g) are used:
C                 gam(g) = sum(i=1 to N) of GAMMA(i)*((g**2/(4*pi**2))**i).
C                 GAMMA(i) = sum(j=1 t0 i) of gammnf(i,j)*(nlfl**(j-1))
C                 gammnf(i,j) = GARRAY(j+i*(i-1)/2)
C                        where nlfl is the number of "light" flavors.
C
       DIMENSION GAMMA(5),GARRAY(10)
       LOGICAL ALPHASET, AmHatSet
       COMMON /QCDPAR_CT/ FSWITCH, AL, NF, NORDER, ALPHASET
       COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
       COMMON / IOUNIT / NIN, NOUT, NWRT
       GAMMNF(K,L)=GARRAY(L+K*(K-1)/2)
C

cpn11      IF (.NOT.ALPHASET) CALL LAMCWZ
      if (.not.alphaset) stop 'SUD: Lambdas are not set'
 
      SUD=0.
      IF(N.LE.0.) RETURN
      IF((Q1.LE.AMN).OR.(Q2.LE.AMN)) THEN
         WRITE(NOUT, *) 'Q1 OR/AND Q2 IS TOO SMALL IN SUD'
         RETURN
      ENDIF
C
      NMMIN=NF-NFL(Q2)
      NMMAX=NF-NFL(Q1)
      T10=2.* LOG(Q2)
      B1=FLOAT(33-2*(NF-NMMIN+1))/12.
      B2B1SQ=FLOAT(153-19*(NF-NMMIN+1))/(24.*B1*B1)
      GAMMA(2)=GAMMNF(2,1)+GAMMNF(2,2)*(NF-NMMIN+1)
      DO 200 J=NMMIN+1,NMMAX+1
         T20=T10
         B2B1SQ=(B2B1SQ*B1*B1+19./24.)/(B1+1./6.)**2
         B1=B1+1./6.
         IF (J.EQ.(NMMAX+1)) THEN
            T10=2.* LOG(Q1)
         ELSE
            T10=2.* LOG(AMHAT(NF-J+1))
         ENDIF
         ALLAM=2.* LOG(ALAM(NF+1-J))
         TQ2=2.* LOG(Q2)-ALLAM
         T1=T10-ALLAM
         T2=T20-ALLAM
         ALNT2= LOG(T2)
         ALNT1= LOG(T1)
         SUD=SUD+GARRAY(1)*0.25/B1*(TQ2*(ALNT2-ALNT1
     1             +B2B1SQ*((ALNT2+1.)/T2-(ALNT1+1.)/T1))
     2             +T1-T2+B2B1SQ*(ALNT2**2-ALNT1**2)/2.)
C
         IF (N.GE.2) THEN
            GAMMA(2)=GAMMA(2)-GAMMNF(2,2)
            SUD=SUD+0.25*GAMMA(2)/(B1*B1)*(ALNT1-ALNT2
     >                                +TQ2*(1./T1-1./T2))
         ENDIF
200   CONTINUE
      RETURN
      END
C
C**********************
C
