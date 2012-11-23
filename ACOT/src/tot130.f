C2345678901234567890123456789012345678901234567890123456789012345678901234567890
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C        PROGRAM test
        Subroutine testtot_OLD()
C-----------------------------------------------------------------------------
C      Program to COMPUTE K FACTORS
C      
C      08/30/05 FIO. Fix arguments of HADNLOQ1230: f1m<==>f2m
C      02/28/05 FIO. Include ACOT-CHI via "Iflg"   (not yet used)
C      02/28/05 FIO. Fix IHadron call; eliminate common block PDFTx
C      
C      02/06/01 FIO. modify to include F123
C      08/16/02 FIO. modify to include more schemes
C      05/07/2007  Include Z-Z, and G-Z terms
C      23 Jan 2008: Include Short-Cut-2: 
C      
C      
c          If(Ischeme.eq.0)   !*** NLO MASSLESS MS-BAR
c          If(Ischeme.eq.1)   !*** FULL ACOT
c          If(Ischeme.eq.2)   !*** FFS
c          If(Ischeme.eq.3)   !*** SHORT-CUT TOT (S-ACOT)
c          If(Ischeme.eq.4)   !*** Test ACOT (NO NLO Q)
c          If(Ischeme.eq.5)   !*** LO
c          If(Ischeme.eq.6)   !*** Massless LO
c          If(Ischeme.eq.7)   !*** Short-cut2: ACOT w/ Massless NLO-Q
c          If(Ischeme.eq.8)   !*** S-ACOT(Chi)  (DISABLED: 24 JAN 2011: AVOID POTENTIAL CONFUSION)
c          If(Ischeme.eq.9)   !*** S-ACOT(Chi) Voica's version
c          
c          If(Ischeme.eq.51)   !*** lo
c          If(Ischeme.eq.52)   !*** gluon sub
c          If(Ischeme.eq.53)   !*** nlo-gluon
c          If(Ischeme.eq.54)   !*** quark sub
c          If(Ischeme.eq.55)   !*** nlo-quark
c          If(Ischeme.eq.56)   !*** Massless NLO gluon   
c          If(Ischeme.eq.57)   !*** Massless NLO quark   
c          If(Ischeme.eq.58)   !*** UNUSED               
c          If(Ischeme.eq.59)   !*** UNUSED
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       CHARACTER HEADER*78
       CHARACTER TITSCH(0:10)*20

       Dimension XMARRAY(6),Xout123( 3)

       PARAMETER (PI=3.14159265359)
       Common  / ActInt /  AERR, RERR, iActL, iActU
C-----------------------------------------------------------------------------
C--- FOR NEUTRAL CURRENT
C      DATA SBIG,XBJ,Q,F1M,F2M,GLQ,GRQ,GLLEP,GRLEP,XMU,ISET,HMASS 
C    >   /  98596.,0.1, 10.,1.6,1.6,0.5,0.5,0.5,0.5,10, 1,0.938/
C      DATA   XLEPOL /2.0/
C      DATA IPARTIN, IPARTOUT,SCALE  /  4, 4,  -1  /
C-----------------------------------------------------------------------------
C--- FOR CHARGED CURRENT
      DATA SBIG,XBJ,Q,F1M,F2M,GLQ,GRQ,GLLEP,GRLEP,XMU,ISET,HMASS 
     >   /  98596.,0.1, 10.,0.5,1.6,1.0,0.0,1.0,0.0,10.0, 1,0.938/
      DATA   XLEPOL /1.0/
      DATA IPARTIN, IPARTOUT,SCALE,Ihad  /  4, 5,  -1 , 1 /
C-----------------------------------------------------------------------------

       DATA XMARRAY  /  0.1, 0.1, 0.2, 1.6, 5.0, 170.0/
       DATA ICHANNEL / 1 /  !*** ANTI-QUARK ICHANNEL
C ********* ICHANNEL is mainly for debugging to turn off Anti-quark Channel
C ********* There is no physics process where you want to run this way
  
      DATA TITSCH /"MS-BAR M=0","Full ACOT","FFS",
     >   "S-ACOT","ACOT (NO q)","LO","LO M=0",
     >   " "," "," "," "/
C----------------------------------------------------------------------
C SETUP INTEGRATION PARAMETERS
C----------------------------------------------------------------------
      AERR=0.0    
      RERR=1.E-4  !*** Setting this too small will yield unstable results *** 
      iactl=2      
      iactu=2      
      WRITE(6,202)  AERR, RERR, iActL, iActU
202   FO RMAT(' AERR= ',1PG14.7,' RERR= ',1PG14.7,
     >       ' IACTL= ',I3,' IACTU= ',I3)
      WRITE(6,*) ' ENTER:  AERR, RERR, iActL, iActU'
      READ (5,*)           AERR, RERR, iActL, iActU
C----------------------------------------------------------------------
C SETUP PARTON NUMBERS 
C----------------------------------------------------------------------
      WRITE(6,203) IPARTIN, IPARTOUT
203   FORMAT(' IPARTIN= ',I3,' IPARTOUT= ',I3)
      WRITE(6,*) ' ENTER: IPARTIN, IPARTOUT '
      READ (5,*)          IPARTIN, IPARTOUT

       F1M=XMARRAY(ABS(IPARTIN))
       F2M=XMARRAY(ABS(IPARTOUT))

       IF(IPARTIN.EQ.IPARTOUT) THEN
         GLQ=0.5
         GRQ=0.5
         GLLEP=0.5
         GRLEP=0.5
         XLEPOL=2.0
       ENDIF

       WRITE(6,*)    GLQ,GRQ,GLLEP,GRLEP,SBIG,XLEPOL
       WRITE(6,*) '> GLQ,GRQ,GLLEP,GRLEP,SBIG,XLEPOL'
       READ( 5,*)    GLQ,GRQ,GLLEP,GRLEP,SBIG,XLEPOL

       E=(SBIG-HMASS**2)/(2.0*HMASS) 
C----------------------------------------------------------------------
C SET THE LOOP
C----------------------------------------------------------------------
1     WRITE(6,102) XBJ,Q
      WRITE(6,101) E,F1M,F2M,ISET,HMASS,SCALE,ICHANNEL,Ihad
102   FORMAT('   XBJ= ',1PG14.7,'     Q= ',1PG14.7)
101   FORMAT('     E= ',1PG14.7,'   F1M= ',1PG14.7,
     >       '   F2M= ',1PG14.7,'  ISET= ',I3, /,
     >       ' HMASS= ',1PG14.7,' SCALE= ',1PG14.7,/,
     >       ' ANTI-QUARK SWITCH = ',I3,' Ihad = ',I3)
C23456789012345678901234567890123456789012345678901234567890123456789012
      WRITE(6,*) 'ENTER: XBJ,Q,E,F1M,F2M,ISET,HMASS,SCALE,ICHANNEL,Ihad'
      READ (5,*)         XBJ,Q,E,F1M,F2M,ISET,HMASS,SCALE,ICHANNEL,Ihad
C----------------------------------------------------------------------
C SET THE SCALES
C----------------------------------------------------------------------
      SBIG = 2.0*E*HMASS + HMASS*HMASS  !*** 11/11/97 RESET SBIG FOR SIG CALC
      Q2=Q*Q
      Y=Q2/(2.0*HMASS*E*XBJ)

C     Q=SQRT(2.0*HMASS*E*XBJ*Y)
      Wsq = HMASS**2 + Q**2 * (1./XBJ - 1.)
      W=SQRT(Wsq)
      Wmin = F1M + F2M + HMASS

      RHO=1.0
      IF(HMASS.GT.0.0) THEN 
        XNU=Q**2/(2.0*HMASS*XBJ)
        RHO=SQRT(1.0+(Q/XNU)**2)  
      ENDIF
      RHO2=RHO*RHO

C----------------------------------------------------
C   ***  COMPUTE MU SCALE
C----------------------------------------------------
       IF(SCALE.GT.0.0) THEN
           XMU = SCALE
       ELSEIF(SCALE.LT.0.0) THEN
           XMU = Q * Abs(scale)
       ENDIF


      ZERO=0.0
C----------------------------------------------------
c          If(Ischeme.eq.0)   !*** NLO MASSLESS MS-BAR
c          If(Ischeme.eq.1)   !*** FULL ACOT
c          If(Ischeme.eq.2)   !*** FFS
c          If(Ischeme.eq.3)   !*** SHORT-CUT TOT (S-ACOT)
c          If(Ischeme.eq.4)   !*** Test ACOT (NO NLO Q)
c          If(Ischeme.eq.5)   !*** LO
c          If(Ischeme.eq.6)   !*** Massless LO
c          If(Ischeme.eq.7)   !*** Short-cut2: ACOT w/ Massless NLO-Q
c          If(Ischeme.eq.8)   !*** S-ACOT(Chi)
c          If(Ischeme.eq.9)   !*** S-ACOT(Chi) Voica's version
C----------------------------------------------------
         write(6,*) 
     >   '  ============================================ '
         write(6,*) 
     >   '  SCHEME          F1           F2           F3 '
C *********** include 2 possible couplings:
      GLQ1=GLQ
      GLQ2=GLQ
      GRQ1=GRQ
      GRQ2=GRQ
      Do Ischeme=0,9,1
         call TOTF(XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,IPARTIN,IPARTOUT,
     >   XMU,ISET,Ihad,HMASS,Ischeme,Iflg,Ichannel,Xout123) 
C----------------------------------------------------
         write(6,1101) titsch(Ischeme),(Xout123(i),i=1,3,1)
1101     Format(A12,3x,3(1x,1pg12.5))
      ENDDO


      GOTO 1
      END
C----------------------------------------------------------------------



C2345678901234567890123456789012345678901234567890123456789012345678901234567890
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
C-----------------------------------------------------------------------------
      SUBROUTINE TOTF(XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   IPARTIN,IPARTOUT,
     >   XMU,ISET,Ihad,HMASS,Ischeme,Iflg,Ichannel,Xout123) 
C      
C      
C      
C      03/01/08 FIO  Move S-ACOT(CHI) INTO TOT-F
C      02/06/01 FIO. modify to includ F123
C      08/16/02 FIO. modify to include more schemes
C      05/07/2007  Include Z-Z, and G-Z terms
C-----------------------------------------------------------------------------
C----------------------------------------------------
c          If(Ischeme.eq.0)   !*** NLO MASSLESS MS-BAR
c          If(Ischeme.eq.1)   !*** FULL ACOT
c          If(Ischeme.eq.2)   !*** FFS
c          If(Ischeme.eq.3)   !*** SHORT-CUT TOT (S-ACOT)
c          If(Ischeme.eq.4)   !*** Test ACOT (NO NLO Q)
c          If(Ischeme.eq.5)   !*** LO
c          If(Ischeme.eq.6)   !*** Massless LO
c          If(Ischeme.eq.7)   !*** Short-cut2: ACOT w/ Massless NLO-Q
c          If(Ischeme.eq.8)   !*** S-ACOT(Chi)
c          If(Ischeme.eq.9)   !*** S-ACOT(Chi) Voica's version
C----------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       CHARACTER HEADER*78
       logical isnan

       Dimension  
     >            XMARRAY(6), Xout123(3)
     >           ,XLO123( 3),XLO123a( 3),XLO123b( 3)
     >           ,XLO1230(3),XLO123a0(3),XLO123b0(3)
     >           ,XLO123chk( 3),XLO123achk( 3),XLO123bchk( 3)
     >           ,XLO1230chk(3),XLO123a0chk(3),XLO123b0chk(3)
     >           ,XSUB123( 3),XSUB123a( 3),XSUB123b( 3)
     >           ,XSUB1230(3),XSUB123a0(3),XSUB123b0(3)
     >           ,XSUB123Q(3),XSUB123Qa(3),XSUB123Qb(3)
     >           ,XNLO123(3),XTOT123(3),XDIFF123(3),XDIFF123q(3)
     >           ,XNLO123Qz(3),XNLO123Qza(3),XNLO123Qzb(3)
     >           ,XG123(3),XGR123(3) 
     >           ,XHq123( 3),XHq123a( 3),XHq123b( 3)
     >           ,XHq1230(3),XHq123a0(3),XHq123b0(3)
     >           ,XHqSUBQ123(3),XHqSUBQ123a(3),XHqSUBQ123b(3)
     >           ,XGR123Q(3),XTOT1230(3),XTOT123S(3),XRATQG(3)
     >           ,XTOT123X(3),XTOT123Y(3),Xtot123Schi(3)
       DIMENSION F123(3)
       DATA F123 / +1.d0, +1.d0, +1.d0 /  !*** SIGN OF ANTI-QUARK TERMS
C  ********* SIGN OF ANTI-QUARK TERMS IS AUTOMATIC WITH GL<=>GR
       PARAMETER (PI=3.14159265359)
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C  IF WE WANT LO CALC: DO THIS AND RETURN QUICKLY:   23 JAN 2008 fio
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       IF((Ischeme.EQ.5).OR.(Ischeme.EQ.6)) then  !*** THIS IS LO

       imass=0
       IF(Ischeme.EQ.5) imass=1 !***  MASSIVE LO
       
       CALL  HADLO123(imass,Ischeme,Ihad,+IPARTIN ,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XLO123a)

       CALL  HADLO123(imass,Ischeme,Ihad,-IPARTOUT,XBJ,Q,F2M,F1M,
     >   GRQ1,GLQ1,GRQ2,GLQ2,XMU,ISET,HMASS,XLO123b)

        DO I=1,3,1
         XLO123(   I) = XLO123a(   I) + XLO123b(   I) *F123(I)  !*** Leading Order
         Xout123(I) =  XLO123(I)     !***  Leading Order
        ENDDO

        RETURN
        endif !*** End of special LO Section
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      aqk=0.0                   !*********  Anti-quark Channel
      IF(ICHANNEL.NE.0) aqk=1.0
C ********* ICHANNEL is mainly for debugging to turn off Anti-quark Channel
C ********* There is no physics process where you want to run this way

C============================================================
C   ***  If S-ACOT(Chi) setup Xbj shift to XbjChi
C  Note: only need to shift Lo0, Sub0,NLOq0; NLO has full m1 m2 dependence
c
c  Note: We need to fix F2 as this has a factor of x: 
c        XbjChi will give this an extra factor of fix
C============================================================
         XbjChi=Xbj  !*** default
      if(Ischeme.eq.8) then 
         WRITE(6,*)  ' DISABLED: 24 JAN 2011: AVOID POTENTIAL CONFUSION'
         call HF_errlog(106,'F: TOTF - Ischeme=8 DISABLED: 24 JAN 2011')
c        STOP
         Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
         ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))
         fix=Eta*(1.0d0+(f1m+f2m)**2/Q**2)/xbj
         XbjChi=fix * xbj   !*** For S-ACOT(CHI)
         XbjChi=Min(1.0d0,XbjChi)
         fix=XbjChi/xbj
      endif

C============================================================
C   ***  FIRST LOOK AT THE MASSLESS PARTS
C============================================================
C %%% THIS IS THE MASSLESS MS-BAR FORMULA
       imass=0
C----------------------------------------------------
C   ***  COMPUTE STRUCTURE FUNCTIONS: PICK UP T-CHANNEL GRAPHS
C----------------------------------------------------
       CALL  HADLO123(imass,Ischeme,Ihad,IPARTIN,XBJchi,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XLO123a0)
       CALL  HADLO123chk(imass,Ischeme,Ihad,IPARTIN,XBJchi,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XLO123a0chk)
       if(Ischeme.eq.8) XLO123a0(2)= XLO123a0(2)/fix
       CALL HADSUB123(imass,Ischeme,Ihad,IPARTIN,XBJchi,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XSUB123a0)   
       if(Ischeme.eq.8) XSUB123a0(2)= XSUB123a0(2)/fix
       CALL HADNLOQ1230(Ischeme,Ihad,XBJchi,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,HMASS,+IPARTIN ,ISET,XNLO123Qza)
       if(Ischeme.eq.8) XNLO123Qza(2)= XNLO123Qza(2)/fix
C----------------------------------------------------
C   ***  GO BACK AND PICK UP U-CHANNEL GRAPHS
C----------------------------------------------------
C  switch +IPARTIN -> -IPARTOUT and F2M <=> F1M and GLQ <=> GRQ

       CALL  HADLO123(imass,Ischeme,Ihad,-IPARTOUT,XBJchi,Q,F2M,F1M,
     >   GRQ1,GLQ1,GRQ2,GLQ2,XMU,ISET,HMASS,XLO123b0)
       CALL  HADLO123chk(imass,Ischeme,Ihad,-IPARTOUT,XBJchi,Q,F2M,F1M,
     >   GRQ1,GLQ1,GRQ2,GLQ2,XMU,ISET,HMASS,XLO123b0chk)
       CALL  HADLO123(imass,Ischeme,Ihad,-IPARTOUT,XBJchi,Q,F2M,F1M,
     >   GRQ1,GLQ1,GRQ2,GLQ2,XMU,ISET,HMASS,XLO123b0)
       if(Ischeme.eq.8) XLO123b0(2)= XLO123b0(2)/fix
       CALL HADSUB123(imass,Ischeme,Ihad,-IPARTOUT,XBJchi,Q,F2M,F1M,
     >   GRQ1,GLQ1,GRQ2,GLQ2,XMU,ISET,HMASS,XSUB123b0)  
       if(Ischeme.eq.8) XSUB123b0(2)= XSUB123b0(2)/fix
       CALL HADNLOQ1230(Ischeme,Ihad,XBJchi,Q,F2M,F1M,
     >   GRQ1,GLQ1,GRQ2,GLQ2,XMU,HMASS,-IPARTOUT,ISET,XNLO123Qzb)  
       if(Ischeme.eq.8) XNLO123Qzb(2)= XNLO123Qzb(2)/fix

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C============================================================
C   ***  NEXT LOOK AT THE NON-MASSLESS PARTS
C============================================================
C %%% THIS IS THE MASSIVE ACOT FORMULA
       imass=1
C----------------------------------------------------
C   ***  COMPUTE STRUCTURE FUNCTIONS: PICK UP T-CHANNEL GRAPHS
C----------------------------------------------------
       CALL HADLO123(imass,Ischeme,Ihad,IPARTIN,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XLO123a)
       CALL HADLO123chk(imass,Ischeme,Ihad,IPARTIN,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XLO123aCHK)
       CALL HADSUB123(imass,Ischeme,Ihad,IPARTIN,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XSUB123a)  
       CALL HADNLOQ123K(Ischeme,Ihad, XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,HMASS,IPARTIN,ISET,XHq123a)
       CALL HADSUBQ123(Ischeme,Ihad,IPARTIN,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XHqSUBQ123a)
C----------------------------------------------------
C   ***  GO BACK AND PICK UP U-CHANNEL GRAPHS
C----------------------------------------------------
C  switch +IPARTIN -> -IPARTOUT and F2M <=> F1M and GLQ <=> GRQ

       CALL  HADLO123(imass,Ischeme,Ihad,-IPARTOUT,XBJ,Q,F2M,F1M,
     >   GRQ1,GLQ1,GRQ2,GLQ2,XMU,ISET,HMASS,XLO123b)
       CALL  HADLO123chk(imass,Ischeme,Ihad,-IPARTOUT,XBJ,Q,F2M,F1M,
     >   GRQ1,GLQ1,GRQ2,GLQ2,XMU,ISET,HMASS,XLO123bCHK)
       CALL HADSUB123(imass,Ischeme,Ihad,-IPARTOUT,XBJ,Q,F2M,F1M,
     >   GRQ1,GLQ1,GRQ2,GLQ2,XMU,ISET,HMASS,XSUB123b) 
       CALL HADNLOQ123K(Ischeme,Ihad, XBJ,Q,F2M,F1M,
     >   GRQ1,GLQ1,GRQ2,GLQ2,XMU,HMASS,-IPARTout,ISET,XHq123b)
       CALL HADSUBQ123(Ischeme,Ihad,-IPARTout,XBJ,Q,F2M,F1M,
     >   GRQ1,GLQ1,GRQ2,GLQ2,XMU,ISET,HMASS,XHqSUBQ123b)
C============================================================
C   ***  NOW LOOK AT THE OTHER PARTS
C============================================================
C----------------------------------------------------
C   ***  COMPUTE STRUCTURE FUNCTIONS: These have both T and and U channel
C----------------------------------------------------
C %%% THIS IS THE MASSIVE ACOT FORMULA
       CALL HADNLO123(Ischeme,Ihad,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XNLO123)
C----------------------------------------------------
C   ***  COMPUTE MASSLESS GLUON CONTRIBUTION FOR COMPARISON
C----------------------------------------------------
C %%% THIS IS THE MASSLESS MS-BAR FORMULA
       CALL HADNLOG123z(Ischeme,Ihad,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XG123)

       IF((XMU.LE.F1M).OR.(XMU.LE.F2M)) THEN  !****  Patch: 3 Feb. 2009
          XG123(1)=0.0D0
          XG123(2)=0.0D0
          XG123(3)=0.0D0
       ENDIF


C----------------------------------------------------
C   *** PROTECT NAN FROM NLO-Q PROGRAM: HADNLOQ123K
C   If (m/q)>1e4, we could see NaN's comming from HADNLOQ123K
c   use massless result: XHq123a - XHqSUBQ123a = XNLO123Qza
C----------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
       if(((f1m/q).le.1.e-4).or.((f2m/q).le.1.e-4)) then
         DO I=1,3,1
!            if(IsNaN(XHq123a(i))) XHq123a(i)=XNLO123Qza(i)+XHqSUBQ123a(i)
!            if(IsNaN(XHq123b(i))) XHq123b(i)=XNLO123Qzb(i)+XHqSUBQ123b(i)
           XHq123a(i)=XNLO123Qza(i)+XHqSUBQ123a(i)
           XHq123b(i)=XNLO123Qzb(i)+XHqSUBQ123b(i)
         ENDDO
        endif
C----------------------------------------------------
C   ***  ADD TOGETHER
C   ***  INCLUDE SIGN "F123" FOR ANTI-QUARK TERMS
C----------------------------------------------------

        DO I=1,3,1
         XLO1230(  I) = XLO123a0(  I) + XLO123b0(  I) *AQK*F123(I)  !*** MASSLESS MS-BAR
         XSUB1230( I) = XSUB123a0( I) + XSUB123b0( I) *AQK*F123(I)  !*** MASSLESS MS-BAR

         XLO123(   I) = XLO123a(   I) + XLO123b(   I) *AQK*F123(I)  !*** FULL ACOT
         XSUB123(  I) = XSUB123a(  I) + XSUB123b(  I) *AQK*F123(I)  !*** FULL ACOT

         XNLO123Qz(I) = XNLO123Qza(I) + XNLO123Qzb(I) *AQK*F123(I)  !*** MASSLESS MS-BAR

         XHq123(    I)= XHq123a(    I)+ XHq123b(    I)*AQK*F123(I)  !*** FULL ACOT
         XHqSUBQ123(I)= XHqSUBQ123a(I)+ XHqSUBQ123b(I)*AQK*F123(I)  !*** FULL ACOT

        ENDDO
C23456789012345678901234567890123456789012345678901234567890123456789012


C----------------------------------------------------
C   ***  COMPUTE TOTAL STRUCTURE FUNCTIONS
C----------------------------------------------------
       DO N=1,3,1

          XTOT123(N) = XLO123(N) + XNLO123(N) - XSUB123(N) +   !*** FULL ACOT
     >                        1.0* (XHq123(N) - XHqSUBQ123(N) ) 

          XTOT123S(N)= XLO1230(N)+ XNLO123(N) - XSUB1230(N) +  !*** SHORT-CUT TOT (S-ACOT)
     >                             XNLO123Qz(N)   

          XTOT1230(N)= XLO1230(N)+ XG123(N)    +               !*** NLO MASSLESS MS-BAR
     >                             XNLO123Qz(N)            

          XTOT123X(N) = XLO123(N) + XNLO123(N) - XSUB123(N) +  !*** Test ACOT (-NLO Q)
     >                        0.0* (XHq123(N) - XHqSUBQ123(N) ) 

         XTOT123Y(N) = XLO123(N) + XNLO123(N) - XSUB123(N) +  !*** Short-cut2: ACOT w/ Massless NLO-Q
     >                             XNLO123Qz(N)   

         XTOT123Schi(N)=XLO1230(N)+(XNLO123(N)-XSUB1230(N))+XNLO123Qz(N)    !*** S-ACOT(Chi)


C----------------------------------------------------
C   ***  Set output
C----------------------------------------------------
          Xout123(N) = 0.0d0
          If(Ischeme.eq.0) Xout123(N) =  XTOT1230(N)   !*** NLO MASSLESS MS-BAR
          If(Ischeme.eq.1) Xout123(N) =  XTOT123( N)   !*** FULL ACOT
          If(Ischeme.eq.2) Xout123(N) =  XNLO123(N)    !*** FFS
          If(Ischeme.eq.3) Xout123(N) =  XTOT123S(N)   !*** SHORT-CUT TOT (S-ACOT)
          If(Ischeme.eq.4) Xout123(N) =  XTOT123X(N)   !*** Test ACOT (NO NLO Q)
          If(Ischeme.eq.5) Xout123(N) =  XLO123(N)     !*** LO
          If(Ischeme.eq.6) Xout123(N) =  XLO1230(N)    !*** Massless LO
          If(Ischeme.eq.7) Xout123(N) =  XTOT123Y(N)   !*** Short-cut2: ACOT w/ Massless NLO-Q
          If(Ischeme.eq.8) Xout123(N) =  XTOT123Schi(N)   !*** S-ACOT(CHI)
          If(Ischeme.eq.9) Xout123(N) =  XTOT123Schi(N)   !*** S-ACOT(CHI) Voica's version

          If(Ischeme.eq.51) Xout123(N) =  XLO123(N)   !***  lo       
          If(Ischeme.eq.52) Xout123(N) =  XSUB123(N)   !***  gluon sub
          If(Ischeme.eq.53) Xout123(N) =  XNLO123(N)   !***  nlo-gluon
          If(Ischeme.eq.54) Xout123(N) =  XHqSUBQ123(N)   !***  quark sub
          If(Ischeme.eq.55) Xout123(N) =  XHq123(N)   !***  nlo-quark
          If(Ischeme.eq.55) Xout123(N) =  XG123(N)   !*** Massless NLO gluon
          If(Ischeme.eq.56) Xout123(N) =  XNLO123Qz(N)   !***  Massless NLO quark
          If(Ischeme.eq.57) Xout123(N) =  XG123(N)   !*** UNUSED
          If(Ischeme.eq.58) Xout123(N) =  XG123(N)   !*** UNUSED
          If(Ischeme.eq.59) Xout123(N) =  XG123(N)   !*** UNUSED

       ENDDO
C----------------------------------------------------------------------
C
       return
C----------------------------------------------------------------------
       END
C----------------------------------------------------------------------


C========================================================================
      function  IsNaN(x)
C========================================================================
C     *** Home-made function to detect NAN: FIO  1 April 2010
C========================================================================
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      logical tmp,isnan

      tmp=x.ne.x

      IsNaN=tmp

      RETURN
      END



