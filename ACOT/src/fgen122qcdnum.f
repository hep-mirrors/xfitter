C =========================================================================
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C =========================================================================
       SUBROUTINE Fgen123LK(idata,icharge,Mode, XBJ, Q, XMU, F123L)
C-----------------------------------------------------------------------------
C      This is a front-end for Fgen123L
C      On first call of idata point number, it computes and stores the K-factor
C      25 April 2011: 
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      PARAMETER (ndata = 1000,n123L=4,nMode=3,nTotal=ndata*n123L*nMode)  !*** Number of data points for K-factor array
      Dimension XKFACTOR(ndata,n123L,nMode)  !*** ndata points,  (F123L)=(1,2,3,4), Mode= (F,Fc,Fb)=(1,2,3)
      data XKFACTOR /nTotal*0.d0/  
      save XKFACTOR

      Dimension F123L(4),F123LLO(4)
      Common /Ischeme/ Isch, Iset, Iflg, Ihad  !*** pass info out to Fnc123 and Fcc123
      common /fred/ xmc,xmb,HMASS

      Character*80 Message ! Error message text

C-----------------------------------------------------------------------------
c    if idata>ndata, increase k-factor table
C-----------------------------------------------------------------------------
      if(idata.gt.ndata)  then  !**** over-ride and use full calculation:
         write(6,*) ' Error: idata =',idata,' > ',ndata
         write(6,*) ' Increase ndata '
         write(Message,*) 
     +   'F: Fgen123LK - idata =',idata,' > ',ndata,' Increase ndata!'
         call HF_errlog(101,Message)
c        stop
      endif

C-----------------------------------------------------------------------------
c    if idata=0 skip k-factor table and use full calculation
C-----------------------------------------------------------------------------
      if(idata.eq.0)  then  !**** over-ride and use full calculation:
         call Fgen123L(icharge,Mode,xbj,q,xmu,F123L)
         return
      endif

C-----------------------------------------------------------------------------
c     FIRST TIME THROUGH: FILL K-FACTOR      
C-----------------------------------------------------------------------------
      if(XKFACTOR(idata,1,mode).eq.0.d0)  then  !***  FIRST TIME THROUGH: FILL K-FACTOR  ===
         call Fgen123L(icharge,Mode,xbj,q,xmu,F123L)
         IschORIG=Isch
         Isch=5  !*** Massive LO Calculation
         Call Fgen123L(icharge,Mode,XBJ,Q,XMU,F123Llo)
         Isch=IschORIG  !*** Reset Ischeme
C     Generate K-Factor
         do i=1,4
            if(F123Llo(i).eq.0.0d0) then
            XKFACTOR(Idata,i,MODE)=1.0d0  !**** Default if denom is zero
            else
            XKFACTOR(Idata,i,MODE)=F123l(i)/F123llo(i)
            endif
         enddo
c
C-----------------------------------------------------------------------------
         else  !***  NOT FIRST TIME THROUGH: USE K-FACTOR ======================
            IschORIG=Isch
            Isch=5  !*** Massive LO Calculation
            Call Fgen123L(icharge,Mode, XBJ, Q, XMU, F123Llo)
            Isch=IschORIG  !*** Reset Ischeme
C     Use K-Factor
            do i=1,4
               F123L(i)=F123Llo(i) * XKFACTOR(Idata,i,mode)
            enddo
         endif
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
      return
      end

CC =========================================================================
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C =========================================================================
       SUBROUTINE Fgen123L(icharge,Mode, XBJ, Q, XMU, F123L)
C-----------------------------------------------------------------------------
C      This is a front-end for Fgen123. "Fgen123L" simply adds on the "L" piece
C      Program to compute both CC and NC F123
C      25 April 2011: Call Fgen123L and compute "L" term
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension F123L(4)
      Dimension F123( 3)
      Common /Ischeme/ Isch, Iset, Iflg, Ihad  !*** pass info out to Fnc123 and Fcc123
      common /fred/ xmc,xmb,HMASS

      call Fgen123(icharge,Mode,xbj,q,xmu,F123)

c     copy arrays
      F123L(1)=F123(1)
      F123L(2)=F123(2)
      F123L(3)=F123(3)
C----------------------------------------------------------------------
C COMPUTE  FL 
C----------------------------------------------------------------------
      rho=Sqrt(1.0d0+(2.0d0*hmass*xbj/Q)**2)  !*** Get Hmass from /fred/ common block 
      FL=rho**2*F123(2)- 2.0d0*xbj*F123(1)
      F123L(4)=FL

      return
      end

C =========================================================================
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C =========================================================================
       SUBROUTINE Fgen123(icharge,Mode, XBJ, Q, XMU, F123)
C-----------------------------------------------------------------------------
C      Program to compute both CC and NC F123
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C      
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension F123(3)
      Common /Ischeme/ Isch, Iset, Iflg, Ihad  !*** pass info out to Fnc123 and Fcc123

      Character*80 Message ! Error message text

      if(icharge.eq.0) then !*** Neutral Current  (only photon at this point)
          call Fnc123(icharge,Mode,xbj,q,xmu,F123)

      elseif(icharge.eq. 4) then !*** Neutral Current BOTH GAMMA & Z
       Call Fnc123(icharge,Mode, XBJ, Q,XMU, F123)

      elseif(icharge.eq.+1) then !*** Charged Current (W+) 
       Call Fcc123(icharge,Mode, XBJ, Q,XMU, F123)
      
      elseif(icharge.eq.-1) then !*** Charged Current (W-)
       Call Fcc123(icharge,Mode, XBJ, Q,XMU, F123)

      else
c        write(6,*) ' error: icharge =',icharge,' not implemented'
         write(Message,*)
     +   'F: Fgen123 - icharge =',icharge,' is not implemented'
         call HF_errlog(102,Message)
c        stop
      endif


      return
      end

C-----------------------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C =========================================================================
       SUBROUTINE Fnc123(icharge,Mode, XBJ, Q,XMU, F123)
C-----------------------------------------------------------------------------
C      Program to COMPUTE K FACTORS
C      
C      
C      
C      06/25/99 FIO. 
C      02/01/08 FIO Update couplings; fac of 2 in gz+zg 
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension 
     >    XTOT123(3),XCHARM(3),XBOTTOM(3)
     >   ,XXTOT123(3),XXCHARM(3),XXBOTTOM(3)
     >   ,xmarray(6), charge(6), charge3(6)
     >   ,F123(3)
      Dimension  
     >   T3F(6),
     >   qVECTORg(6), qAXIALg(6), qRightg(6), qLeftg(6),
     >   qVECTORz(6), qAXIALz(6), qRightz(6), qLeftz(6),
     >   XTOT123gz(3), XTOT123zz(3),XXTOT123gz(3), XXTOT123zz(3),
     >   term(3), facgg(3),facgz(3),faczz(3)
      Parameter(Iset4F4=14)
      PARAMETER (PI=3.14159265359)
      Common /Ischeme/ Isch, Iset, Iflg, Ihad
      Common  / ActInt /  AERR, RERR, iActL, iActU

C-----------------------------------------------------------------------------
C     DATA SINW2, XMW, XMZ   / 0.23D0,   80.4D0,  91.2D0 /
! pull values from herafitter
cv      DATA SINW2, XMW, XMZ   / 0.2303D0, 80.0D0,  91.188D0 /   !*** MATCH QCDNUM
      common /fredew/ SINW2, XMW, XMZ

      DATA idebugZ,ifirstZ,iLRdebug,xDebugTmp  /1,0,+1, 1.0d0/
C-----------------------------------------------------------------------------
C--- FOR CHARGED CURRENT
C       DATA GLQ,GRQ,GLLEP,GRLEP,HMASS 
C      >   /  1.0,0.0,1.0,0.0,0.938/
C       DATA   XLEPOL /1.0/
C       DATA IPARTIN, IPARTOUT,SCALE  /  3, 4,  -1  /
C-----------------------------------------------------------------------------
C--- FOR NEUTRAL CURRENT  !*** This is for photon only. Modify for Z
       DATA GLQ,GRQ,GLLEP,GRLEP,HMASS 
     >   /  0.5,0.5,0.5,0.5,0.938/
C-----------------------------------------------------------------------------
      common /fred/ xmc,xmb,Hmass  !*** PULL VALUES FROM QCDNUM  fio 14 FEB. 2011
C                          U      D      S      C      B      T
C      DATA XMARRAY  /    0.1,   0.1,   0.2,   1.6,   5.0, 175.0/  !****  quarks:   U,D,S,C,B,T  CTEQ4 ===
       DATA XMARRAY  /    0.1,   0.1,   0.2,   1.3,   4.5, 175.0/  !****  quarks:   U,D,S,C,B,T  CTEQ6 ===
       DATA CHARGE3  /   +2.0,  -1.0,  -1.0,  +2.0,  -1.0,  +2.0/  !****  quarks:   U,D,S,C,B,T
       DATA T3F      /   +0.5,  -0.5,  -0.5,  +0.5,  -0.5,  +0.5/  !****  quarks:   U,D,S,C,B,T
       DATA ICHANNEL  /  1 /
C----------------------------------------------------------------------
C PULL MC AND MB FROM QCDNUM USING 
C      common /fred/ xmc,xmb
C----------------------------------------------------------------------
      XMARRAY(4)=XMC !*** PULL VALUES FROM QCDNUM  fio 14 FEB. 2011
      XMARRAY(5)=XMB !*** PULL VALUES FROM QCDNUM  fio 14 FEB. 2011
C----------------------------------------------------------------------
C WEINBERG ANGLE:
C----------------------------------------------------------------------
       THETAW=ASIN(SQRT(SINW2))
       hmass=0.0d0 !********* PATCH FOR TESTING
C----------------------------------------------------------------------
C INITIALIZATION OF COUPLINGS
C----------------------------------------------------------------------
       iConv=-1   !*** CAUTION: Use ACOT Convention: Not Standard Convention
                  !*** This controls the sign of the F3 term 

       DO I=1,6,1   !*** Loop over quarks: TOTF picks up anti-quarks
C------- Photon couplings
       CHARGE(I)=CHARGE3(I)/3

       qVECTORg(I) = CHARGE(I)
       qAXIALg(I)  = 0.0d0
       qRightg(i)  = (qVECTORg(I) + iConv * qAXIALg(I) )/2.0d0
       qLeftg(i)   = (qVECTORg(I) - iConv * qAXIALg(I) )/2.0d0  !*** CAUTION:Acot Conv. Non-standard

C------- Z couplings
       ii=Abs(I)  !**** For TF3(ii) 
c      qVECTORz(I)=(T3F(II) -2.0D0 * CHARGE(I)* SINW2)/SIN(2.0D0*THETAW) !*** den. is in kappa FIO: 1 feb 08
c      qAXIALz(I) =(T3F(II)                          )/SIN(2.0D0*THETAW) !*** den. is in kappa FIO: 1 feb 08
C
       qVECTORz(I)=(T3F(II) -2.0D0 * CHARGE(I)* SINW2)
       qAXIALz(I) =(T3F(II)                          )
       qRightz(i) = (qVECTORz(I) + iConv * qAXIALz(I) )/2.0d0
       qLeftz(i)  = (qVECTORz(I) - iConv * qAXIALz(I) )/2.0d0  !*** CAUTION:Acot Conv. Non-standard
      ENDDO
C----------------------------------------------------------------------
C INITIALIZATION OF LEPTON COUPLINGS
C----------------------------------------------------------------------
C------- Z-ELECTRON COUPLINGS:
       eleVec=   -0.5d0 +  2.0d0 * SINW2
       eleAxial= -0.5d0

C------- Combinations that multiply F123
       facgg(2) = 1.0d0
       facgz(2) = - eleVec 
       faczz(2) = + (eleVec**2 + eleAxial**2)
 
       facgg(1) = facgg(2) 
       facgz(1) = facgz(2) 
       faczz(1) = faczz(2) 

       facgg(3) = 0.0d0
       facgz(3) = -  eleAxial
       faczz(3) = + 2.0d0 * (eleVec * eleAxial)


C----------------------------------------------------------------------
C INITIALIZATION OF PROPAGATOR FACTORS
C REFERENCE: H1 COLLAB: Eur. Phys. J. C 30, P.1 (2003) 
C----------------------------------------------------------------------
C  xkappa = 1/[ 4 Sin[w]^2 Cos[w]^2 ]
c
      XKAPPAold=1.0D0/(4.0D0*XMW**2/XMZ**2 * (1.0D0- (XMW**2/XMZ**2)))
      XKAPPA   =1.0D0/(4.0D0* Sin(thetaw)**2 *  Cos(thetaw)**2 )
      ZPROP = XKAPPA * (Q**2/(Q**2 + XMZ**2))

C----------------------------------------------------------------------
C SET THE LOOP
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C INITIALIZATION
C----------------------------------------------------------------------
          DO I=1,3,1
             XXTOT123( I)  = 0.0
             XXTOT123gz( I)  = 0.0
             XXTOT123zz( I)  = 0.0
             XXCHARM(  I)  = 0.0
             XXbottom( I)  = 0.0
          ENDDO
C----------------------------------------------------
C   *************************************************
C----------------------------------------------------
C   ***  LOOP OVER PARTON FLAVORS: 
C----------------------------------------------------
      DO IPARTIN=1,5,1
C Note: we loop over quarks; TOTF picks up both Q and Q-bar

      IPARTOUT=IPARTIN         !*** NEUTRAL CURRENT 
      F1M=XMARRAY(ABS(IPARTIN))
      F2M=XMARRAY(ABS(IPARTOUT))
C----------------------------------------------------
C   ***  COMPUTE STRUCTURE FUNCTIONS
C----------------------------------------------------
C      Ischeme= 0  !*** NLO Massless MS-Bar
C      Ischeme= 1  !*** Full ACOT Scheme 
C      Ischeme= 2  !*** FFS 
C      Ischeme= 3  !*** Simplified ACOT Scheme 
C      Ischeme= 4  !*** Test Full ACOT Scheme (no NLO Q)
C      Ischeme= 5  !*** LO
C      Ischeme= 6  !*** Massless LO
c      Ischeme= 7  !*** Short-cut2: ACOT w/ Massless NLO-Q
c      Ischeme= 8  !*** S-ACOT(Chi)
C----------------------------------------------------
      Ischeme= Isch  !***  Read in from common block
      
C----------------------------------------------------
C   ***  DO PHOTON-PHOTON TERMS
C----------------------------------------------------
      GLQ1= qLeftg(IPARTIN)
      GLQ2= qLeftg(IPARTIN)
      GRQ1= qRightg(IPARTIN)
      GRQ2= qRightg(IPARTIN)
      call TOTF(XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,IPARTIN,IPARTOUT,
     >   XMU,ISET,Ihad,HMASS,Ischeme,Iflg,Ichannel,XTOT123) 

C----------------------------------------------------
C   ***  NOW DO PHOTON-Z INTEREFERENCE TERMS
C----------------------------------------------------
      GLQ1= qLeftg(IPARTIN)
      GLQ2= qLeftz(IPARTIN)
      GRQ1= qRightg(IPARTIN)
      GRQ2= qRightz(IPARTIN)
      call TOTF(XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,IPARTIN,IPARTOUT,
     >   XMU,ISET,Ihad,HMASS,Ischeme,Iflg,Ichannel,XTOT123gz) 

C----------------------------------------------------
C   ***  NOW DO Z-Z TERMS
C----------------------------------------------------
      GLQ1= qLeftz(IPARTIN)
      GLQ2= qLeftz(IPARTIN)
      GRQ1= qRightz(IPARTIN)
      GRQ2= qRightz(IPARTIN)
      
cv
cv      hmass = 0.938
cv      ischeme=1
cv      xmu=q
      call TOTF(XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,IPARTIN,IPARTOUT,
     >   XMU,ISET,Ihad,HMASS,Ischeme,Iflg,Ichannel,XTOT123zz) 

cv      print*,'???????',XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,IPARTIN,IPARTOUT,
cv     >   XMU,ISET,Ihad,HMASS,Ischeme,Iflg,Ichannel,(XTOT123zz(i),i=1,3)
cv      stop
C----------------------------------------------------
C   ***  ADD PARTON CONTRIBUTION TO TOTAL STRUCTURE FUNCTIONS
c   ***  Note: quark coupling already included in GLQ GRQ: FIO 5/7/07
C----------------------------------------------------
      if(ifirstZ.eq.0) then
         ifirstZ=1
         write(6,*) ' GZ and ZZ are for testing '
C         write(6,*) ' enter: 0 to turn off  GZ and ZZ terms '
C         read(5,*) idebugZ  !*** PATCH: TURN ON BY DEFAULT FOR NOW
         if(idebugZ.ne.0) then
            idebugZ=1           !*** error checking
            write(6,*) '  GZ and ZZ terms are ON ',idebugZ
         else
            write(6,*) '  GZ and ZZ terms are OFF ',idebugZ
         endif
      endif
      xz=0.0d0  !**** TURN OF ZZ AND GZ FOR ICHARGE.NE.4
      if(icharge.eq.4) xz=1.0d0 * idebugZ  !*** idebugZ is debug over-ride

C----------------------------------------------------
           DO I=1,3,1
                 tmpGG= facgg(i) * XTOT123(I)                   !*** GG Term 
                 tmpGZ= facgz(i) * zprop * XTOT123gz(I)         !*** GZ Term 
                 tmpZZ= faczz(i) * zprop*zprop * XTOT123zz(I)   !*** ZZ Term 
                 term(i)  =tmpGG +  2.0d0* xz* tmpGZ + xz* tmpZZ !*** pick up both GZ and ZG FIO: 1 Feb 2008
           ENDDO

C----------------------------------------------------
C   ***  ADD PARTON CONTRIBUTION TO TOTAL STRUCTURE FUNCTIONS
C----------------------------------------------------
           DO I=1,3,1
             XXTOT123(I)  = XXTOT123(I) +   term(i)    
           ENDDO
C----------------------------------------------------

C----------------------------------------------------
C   ***  PICK OUT TOT F-CHARM AND SAVE THIS
C----------------------------------------------------
       IF(ABS(IPARTIN).EQ.4) THEN    
           DO I=1,3,1
             XXCHARM(I)  = XXCHARM(I) +   term(i) 
           ENDDO
        ENDIF
C----------------------------------------------------
C   ***  PICK OUT TOT F-BOTTOM AND SAVE THIS
C----------------------------------------------------
       IF(ABS(IPARTIN).EQ.5) THEN    
           DO I=1,3,1
             XXBOTTOM(I)  = XXBOTTOM(I) +  term(i) 
           ENDDO
        ENDIF

888   ENDDO
C----------------------------------------------------
C   ***  END LOOP OVER PARTON FLAVORS: 
C----------------------------------------------------
 
C----------------------------------------------------
C   ***  SET RETURN VALUES: 
C   ***  NOTE: BY RESTRICTING LOOP=[1,5] ON QUARKS, 
C   ***     WE TAKE CARE OF FACTOR OF 2 NORM IN PREVIOUS VERSION 
C----------------------------------------------------
      Do i=1,3,1
      If     (Mode.eq.1) Then
         F123(i)= XXTOT123(i)    !*** Match standard normalization
      Elseif (Mode.eq.2) Then
         F123(i)= XXCHARM(i)     !*** Match standard normalization
      Elseif (Mode.eq.3) Then
         F123(i)= XXBOTTOM(i)    !*** Match standard normalization
      Endif

      enddo

      RETURN
      END 

C----------------------------------------------------------------------
