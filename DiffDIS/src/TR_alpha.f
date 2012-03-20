c=======================================================
      DOUBLE PRECISION FUNCTION ALPHAzo(T)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION mc2,mb2
      COMMON/TRSFPAR/alambda,flavor,qsct,qsdt,iord
c-- qsdt = 4*mc^2, qsct = 4*mb^2 -- really!
c-- flavor not used
c-- alambda = 4-flavour Lambda_QCD
c-- iord = 0/1  LO/NLO

      parameter (pi = 3.1415926535897932384626433832795D0)
      DATA TOL/.0005/
      ITH=0
      TT=T
      mc2=qsdt/4.
      mb2=qsct/4.
      AL=ALAMBDA
      AL2=AL*AL
      FLAV=4.
      QS=AL2*dEXP(T)

      if(qs.lt.0.5d0) then   !!  running stops below 0.5
          qs=0.5d0
          t=dlog(qs/al2)
          tt=t
      endif

      IF(QS.gt.mb2) then
        ITH=1
        T=dLOG(mb2/AL2)
      ELSE IF(QS.lt.mc2) THEN
        ITH=1
        T=dLOG(mc2/AL2)
        GO TO 311
      END IF
      
   11 CONTINUE
      B0=11-2.*FLAV/3.
      IF(IORD.le.0) then
        ALPHAzo=4.*PI/B0/T
        RETURN
      END IF
      
      X1=4.*PI/B0
      B1=102.-38.*FLAV/3.
      X2=B1/B0**2
      AS=X1/T*(1.-X2*dLOG(T)/T)

    5 CONTINUE
      F=-T+X1/AS-X2*dLOG(X1/AS+X2)
      FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
      AS2=AS-F/FP
      DEL=ABS(F/FP/AS)
      IF((DEL-TOL).le.0) then
        ALPHAzo=AS2
        IF(ITH.EQ.0) RETURN
        GO TO (13,14,15) ITH
      END IF

    4 CONTINUE
      AS=AS2
      GO TO 5
      
   13 ALFQC4=ALPHAzo
      FLAV=5.
      ITH=2
      GO TO 11
   14 ALFQC5=ALPHAzo
      ITH=3
      T=TT
      GO TO 11
   15 ALFQS5=ALPHAzo
      ALFINV=1./ALFQS5+1./ALFQC4-1./ALFQC5
      ALPHAzo=1./ALFINV
      RETURN

  311 CONTINUE
      B0=11-2.*FLAV/3.
      IF(IORD.le.0) then
        ALPHAzo=4.*PI/B0/T
        RETURN
      END IF

      X1=4.*PI/B0
      B1=102.-38.*FLAV/3.
      X2=B1/B0**2
      AS=X1/T*(1.-X2*dLOG(T)/T)

   35 CONTINUE
      F=-T+X1/AS-X2*dLOG(X1/AS+X2)
      FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
      AS2=AS-F/FP
      DEL=ABS(F/FP/AS)
      IF(DEL-TOL)33,33,34
   33 CONTINUE
      ALPHAzo=AS2
      IF(ITH.EQ.0) RETURN
      GO TO (313,314,315) ITH
   34 CONTINUE
      AS=AS2
      GO TO 35
  313 ALFQC4=ALPHAzo
      FLAV=3.
      ITH=2
      GO TO 311
  314 ALFQC3=ALPHAzo
      ITH=3
      T=TT
      GO TO 311
  315 ALFQS3=ALPHAzo
      ALFINV=1./ALFQS3+1./ALFQC4-1./ALFQC3
      ALPHAzo=1./ALFINV
      RETURN
      END

