      FUNCTION ANOM(Q1,Q2,GARRAY,N)
C                                                   -=-=- anom
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Returns integral from Q1 to Q2 of
C                 (dmu/mu) * gam(g(mu))
C              where N terms in anomalous dim. gam(g) are used:
C                 gam(g) = sum(i=1 to N) of gamma(i)*((g**2/(4*pi**2))**i).
C                 GAMMA(i) = sum(j=1 to i) of gammnf(i,j)*(nlfl**(j-1))
C                   except GAMMA(1)=gammnf(1,1)+gammnf(1,2)*nlfl
C                 gammnf(i,j) = GARRAY(1+j+i*(i-1)/2)
C                   except gammnf(1,1) = GARRAY(1) & gammnf(1,2)= GARRAY(2)
C                       where nlfl is the number of "light" flavor.
      logical AlphaSet, AmHatSet
      DIMENSION GAMMA(5),GARRAY(10)
      COMMON /QCDPAR_CT/ FSWITCH, AL, NF, NORDER, ALPHASET
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
      COMMON /IOUNIT/ NIN, NOUT, NWRT
C
cpn11      IF (.NOT.ALPHASET) CALL LAMCWZ
      if (.not.alphaset) stop 'ANOM: Lambdas are not set'

      ANOM=0.
      IF (N.LE.0.) RETURN
      IF ((Q1.LE.AMN).OR.(Q2.LE.AMN)) THEN
         WRITE (NOUT, *) 'Q1 OR/AND Q2 IS TOO SMALL IN ANOM'
         RETURN
      ENDIF
C
      NMMIN=NF-NFL(Q2)
      NMMAX=NF-NFL(Q1)
C
  90  T1= LOG(Q2)
      B1=FLOAT(33-2*(NF-NMMIN+1))/12.
      GAMMA(1)=GARRAY(1)+GARRAY(2)*(NF-NMMIN+1)
      DO 200 J=NMMIN+1,NMMAX+1
         T2=T1
         B1=B1+1./6.
         GAMMA(1)=GAMMA(1)-GARRAY(2)
         IF (J.EQ.(NMMAX+1)) THEN
            T1= LOG(Q1)
         ELSE
            T1= LOG(AMHAT(NF-J+1))
         ENDIF
         ANOM=ANOM+0.5*GAMMA(1)* LOG((T2- LOG(ALAM(NF+1-J)))/
     >        (T1- LOG(ALAM(NF+1-J))))/B1
 200  CONTINUE
      RETURN
      END
C
C********************
C
