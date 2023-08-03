
C =========================================================================
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C =========================================================================
      SUBROUTINE Fgen123Lxcb(icharge, XBJ, Q, XMU, F123Lxcb,
     > polar)
C-----------------------------------------------------------------------------
C      Program to compute both CC and NC F123
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C      
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension F123Lxcb(3,4), F123L(4),F123Lc(4),F123Lb(4)
      Common /Ischeme/ Isch, Iset, Iflg, Ihad !*** pass info out to Fnc123 and Fcc123

      Character*80 Message ! Error message text

C-----------------------------------------------------------------------------
      if(icharge.eq.0) then !*** Neutral Current  (only photon for icharge=0)
c         call Fnc123Lxcb( icharge,xbj,q,xmu,F123Lxcb, polar)
          call Fnc123Lxcb2(icharge,xbj,q,xmu,F123Lxcb, polar) !*** 

C-----------------------------------------------------------------------------
      elseif(icharge.eq. 4.or.icharge.eq.5) then !*** Neutral Current BOTH GAMMA & Z
c       Call Fnc123Lxcb( icharge, XBJ, Q,XMU, F123Lxcb, polar)
        Call Fnc123Lxcb2(icharge, XBJ, Q,XMU, F123Lxcb, polar)
 
C-----------------------------------------------------------------------------
      elseif((icharge.eq.+1).or.(icharge.eq.-1)) then !*** Charged Current (W+ or W-) 
       Call Fcc123L(icharge,1, XBJ, Q,XMU, F123L)
c      Call Fcc123L(icharge,2, XBJ, Q,XMU, F123Lc)  !*** not yet implemented
c      Call Fcc123L(icharge,3, XBJ, Q,XMU, F123Lb)  !*** not yet implemented
     
      Do j=1,4,1  !*** 3='xcb', 4='123L'
         F123Lxcb(1,j)=F123L( j)    
         F123Lxcb(2,j)=0.0d0    !*** not yet implemented
         F123Lxcb(3,j)=0.0d0    !*** not yet implemented
      enddo
C-----------------------------------------------------------------------------
      else
c        write(6,*) ' error: icharge =',icharge,' not implemented'
         write(Message,*)
     +   'F: Fgen123 - icharge =',icharge,' is not implemented'
         call HF_errlog(102,Message)
c        stop
      endif


      return
      end

C =========================================================================
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C =========================================================================
      SUBROUTINE Fgen123LxcbQCDNUM(index,icharge, XBJ,Q,XMU, F123Lxcb,
     > polarity)
C-----------------------------------------------------------------------------
C      07 Nov. 2012: Use QCDNUM to compute the denominator
C      
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension F123Lxcb(3,4), F123L(4),F123Lc(4),F123Lb(4)
      Common /Ischeme/ Isch, Iset, Iflg, Ihad  !*** pass info out to Fnc123 and Fcc123

      Character*80 Message ! Error message text
      character*(10) XSecType
c      character*(5) XSecType
c      character*(*)  dummy

C-----------------------------------------------------------------------------
c    ZERO ARRAY

      do i=1,3
         do j=1,4
            F123Lxcb(i,j)=0.0d0 
         enddo
      enddo

C-----------------------------------------------------------------------------
c     CALL QCDNUM ROUTINE:
c    
c     icharge_in: 0 NC: photon exchange only
c     icharge_in: 4 NC: e+ gamma+gammaZ+Z 
c     icharge_in: 5 NC: e- gamma+gammaZ+Z 
c     icharge_in:-1 CC e-
c     icharge_in:+1 CC e+

C     THIS IS NC-DIS
      if(icharge.eq.4) then
         charge=+1
      elseif(icharge.eq.5) then
         charge=-1
      else     
         write(6,*) ' ERROR: ACOT ICHARGE = ',ICHARGE
         charge=0
         write(6,*) ' SET CHARGE = O',CHARGE
c         stop
      endif

      q2=q*q
      npts=1
      XSecType='NCDIS'
      
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
      local_hfscheme=0  !*** PATCH
      call UseZmvnsScheme(F2, FL, xF3, F2gamma, FLgamma,
     $     q2, xbj, npts, polarity, charge, XSecType,local_hfscheme)
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
    
C-----------------------------------------------------------------------------
c     COPY QCDNUM RESULTS TO ARRAY
c     !*** 3='xcb', 4='123L'

      F123Lxcb(1,1)=(F2-FL)/(2.0D0*XBJ)  !*** Patch for F1
      F123Lxcb(1,2)=F2
      F123Lxcb(1,3)=xF3/xbj  !*** Convention for xF3
      F123Lxcb(1,4)=FL

c      write(6,*) F2, FL, xF3, F2gamma, FLgamma
c      write(6,*) (F123Lxcb(1,j),j=1,4)
C-----------------------------------------------------------------------------

      return
      end
C-----------------------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C =========================================================================
       SUBROUTINE Fnc123Lxcb2(icharge, XBJ, Q,XMU, F123Lxcb, polar)
C-----------------------------------------------------------------------------
C      Program to COMPUTE K FACTORS
C      
C      this adds N2LO and N3LO to Fnc123Lxcb
C      
C      27 OCT 2012: FIO
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension 
     >    XTOT123(3),XCHARM(3),XBOTTOM(3)
     >   ,XXTOT123(3),XXCHARM(3),XXBOTTOM(3)
     >   ,xmarray(6), charge(6), charge3(6)
     >   ,F123(3),  F123Lxcb(3,4), F123Lxcb2(3,4), ratio(3,4)
      double precision fij(0:6,1:6)  !*** this is a dummy; we don't use detailed info

      Dimension  
     >   T3F(6),
     >   qVECTORg(6), qAXIALg(6), qRightg(6), qLeftg(6),
     >   qVECTORz(6), qAXIALz(6), qRightz(6), qLeftz(6),
     >   XTOT123gz(3), XTOT123zz(3),XXTOT123gz(3), XXTOT123zz(3),
     >   term(3), facgg(3),facgz(3),faczz(3)
      Parameter(Iset4F4=14)
      PARAMETER (PI=3.14159265359)

c      Common /Iacot/  nord  !*** pass nord to ACOT module from Subroutine SetHFSCHEME
      Common /Ischeme/ Isch, Iset, Iflg, Ihad
      Common  / ActInt /  AERR, RERR, iActL, iActU
      logical ifirst
      data ifirst /.true./
      save ifirst
#include "steering.inc"

C-----------------------------------------------------------------------------
C     DATA SINW2, XMW, XMZ   / 0.23D0,   80.4D0,  91.2D0 /
!C-----------------------------------------------------------------------------

C-----------------------------------------------------------------------------
      if((icharge.eq.0).or.(icharge.eq.4).or.(icharge.eq.5)) then !*** Neutral Current  (only photon for icharge=0)
          call Fnc123Lxcb(icharge,xbj,q,xmu,F123Lxcb2, polar) !*** 
      endif


      do i=1,3
         do j=1,4
            F123Lxcb(i,j)=F123Lxcb2(i,j)
         enddo
      enddo

c     ==============================================================
C     NEED BETTER LOGIC HERE TO SKIP FOR LOW ORDER CALC
c     if(kord.le.1) return 
      if(isch.eq.5)         return !*** If doing massive LO, return without N3LO


c     ==============================================================
c     ==============================================================
      nord=3                    !*********** PATCH NEED TO LINK DYNAMIC: FIO: 31 MAR 2022
      nord=2                    !*********** PATCH NEED TO LINK DYNAMIC: FIO: 31 MAR 2022
      nord=1                    !*********** PATCH NEED TO LINK DYNAMIC: FIO: 31 MAR 2022
      nord=1                    !*********** PATCH NEED TO LINK DYNAMIC: FIO: 31 MAR 2022
c     ==============================================================
c     ==============================================================
      

      if(ifirst) then 
c     *** UPDATED 19 APRIL 2016: FIO: 
C     *** GET "NORD" FROM COMMON BLOCK PASSED FROM STEERING.TXT
c         nord=1
         write(6,*) ' First time in ACOT module ' 
         write(6,*) ' set NORD =2,3 FOR N2LO OR N3LO  ' 
         write(6,*) ' set NORD to any other value to skip  N2LO OR N3LO' 
         write(6,*) ' NORD =',nord
c         read( 5,*)  nord
C         if(nord.ne.1) then 
            open(63,file='output/KfactorsACOT3.txt')
            write(63,*) ' OUTPUT NLO AND N3LO K-FACTORS: NORD = ',nord
            write(63,*) 
     >   '  icharge, XBJ, Q,XMU, polar, ',
     >   '  ratios-F123L NxLO/NLO, NLO F123L, NxLO F123L '
C          endif
         ifirst=.false.
      endif

c     SKIP IF NOT N2LO OR N3LO
      if((nord.ne.2).and.(nord.ne.3)) return
         


c     icharge_in: 0 NC: photon exchange only
c     icharge_in: 4 NC: e+ gamma+gammaZ+Z 
c     icharge_in: 5 NC: e- gamma+gammaZ+Z 
c 
c     integer iord ! perturbative order: including terms up to O(alpha_s^iord)
c     integer iboson ! chose exchange boson: 0: full NC ew, 1:gamma gamma
c     integer isf ! choose structure function: 0: FL; 1,2,3: F_1,2,3
c
      if(icharge.eq.0) then 
         iboson=1
      elseif(icharge.eq.4) then 
         iboson=0 
      elseif(icharge.eq.5) then 
         iboson=0
      else
         write(6,*) ' error: icharge = ',icharge
         stop
      endif

          xmuf2=xmu**2
          xmur2=xmu**2
          xnCHI=2 !*** Use ACOT-Chi type scaling
          q2=q**2
          x=xbj

          hmass=0.938d0      !*** Get Hmass from /fred/ common block 
          rho=Sqrt(1.0d0+(2.0d0*hmass*x/q)**2)  !*** used for F1<=>F2,FL conversion

        mord = 1 ! O(alpha_s) !-----------------------------------------------
        call  zmCHI(xnCHI,x,Q2,xmuf2,xmur2,mord,iboson,0,
     >   fLres1,fLc1,fLb1,fLt1,fij)
        call  zmCHI(xnCHI,x,Q2,xmuf2,xmur2,mord,iboson,2,
     >   f2res1,f2c1,f2b1,f2t1,fij)

         f1res1=(rho**2*f2res1-fLres1)/(2.0d0*x)
         f1c1  =(rho**2*f2c1  -fLc1  )/(2.0d0*x)
         f1b1  =(rho**2*f2b1  -fLb1  )/(2.0d0*x)

C       nord = 3 ! O(alpha_s^3) !-----------------------------------------------
c       nord is set above:
        call  zmCHI(xnCHI,x,Q2,xmuf2,xmur2,nord,iboson,0,
     >   fLres3,fLc3,fLb3,fLt3,fij)
        call  zmCHI(xnCHI,x,Q2,xmuf2,xmur2,nord,iboson,2,
     >   f2res3,f2c3,f2b3,f2t3,fij)

         f1res3=(rho**2*f2res3-fLres3)/(2.0d0*x)
         f1c3  =(rho**2*f2c3  -fLc3  )/(2.0d0*x)
         f1b3  =(rho**2*f2b3  -fLb3  )/(2.0d0*x)


C     Get extra contribution for N3LO
        dfLres31=fLres3-fLres1
        dfLc31=fLc3-fLc1
        dfLb31=fLb3-fLb1

        df1res31=f1res3-f1res1
        df1c31=f1c3-f1c1
        df1b31=f1b3-f1b1

        df2res31=f2res3-f2res1
        df2c31=f2c3-f2c1
        df2b31=f2b3-f2b1


C      F123Lxcb(i,j): i=Tot,C,B,  j=1,2,3,L

        F123Lxcb(1,1)= F123Lxcb2(1,1) + df1res31
        F123Lxcb(1,2)= F123Lxcb2(1,2) + df2res31
        F123Lxcb(1,3)= F123Lxcb2(1,3) + 0.0d0   !*** not implemented
        F123Lxcb(1,4)= F123Lxcb2(1,4) + dfLres31

        F123Lxcb(2,1)= F123Lxcb2(2,1) + df1c31
        F123Lxcb(2,2)= F123Lxcb2(2,2) + df2c31
        F123Lxcb(2,3)= F123Lxcb2(2,3) + 0.0d0  !*** not implemented
        F123Lxcb(2,4)= F123Lxcb2(2,4) + dfLc31

        F123Lxcb(3,1)= F123Lxcb2(3,1) + df1b31
        F123Lxcb(3,2)= F123Lxcb2(3,2) + df2b31
        F123Lxcb(3,3)= F123Lxcb2(3,3) + 0.0d0  !*** not implemented
        F123Lxcb(3,4)= F123Lxcb2(3,4) + dfLb31


        do i=1,3
           do j=1,4
              if(F123Lxcb2(i,j).ne.0.0d0) then
                 ratio(i,j)= F123Lxcb(i,j)/ F123Lxcb2(i,j)
              else
                 ratio(i,j)=0.0d0
              endif
           enddo
        enddo

        
        write(6,181) (( ratio(i,j),j=1,4),i=1,3)
 181    format('  N3LO/NLO =  '14x,12(f7.2,1x))  !*** need 28x for alignment

         write(63,*) icharge, XBJ, Q,XMU, polar,
     >   (( ratio(i,j),j=1,4),i=1,3), 
     >   (( F123Lxcb2(i,j),j=1,4),i=1,3),
     >   (( F123Lxcb( i,j),j=1,4),i=1,3)


      RETURN
      END 

C----------------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C =========================================================================
       SUBROUTINE Fnc123Lxcb(icharge, XBJ, Q,XMU, F123Lxcb, polar)
C-----------------------------------------------------------------------------
C      Program to COMPUTE K FACTORS
C      
C      
C      
C      06/25/99 FIO. 
C      02/01/08 FIO Update couplings; fac of 2 in gz+zg 
c      04/13/12 FIO: Compute Tot,c,b, in single pass
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension 
     >    XTOT123(3),XCHARM(3),XBOTTOM(3)
     >   ,XXTOT123(3),XXCHARM(3),XXBOTTOM(3)
     >   ,xmarray(6), charge(6), charge3(6)
     >   ,F123(3),  F123Lxcb(3,4)
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
      Data small / 1.0d-16 /

C-----------------------------------------------------------------------------
C     DATA SINW2, XMW, XMZ   / 0.23D0,   80.4D0,  91.2D0 /
! pull values from herafitter
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
       DATA GLQ,GRQ,GLLEP,GRLEP   !**** ,HMASS 
     >   /  0.5,0.5,0.5,0.5/   !*** get HMASS from common block         FIO  13 April 2012
C-----------------------------------------------------------------------------
C                          U      D      S      C      B      T
       DATA XMARRAY  / 0.1d0, 0.1d0, 0.2d0, 1.3d0, 4.5d0, 175.0d0/  !****  quarks:   U,D,S,C,B,T  CTEQ6 ===
       DATA CHARGE3  /   +2.0,  -1.0,  -1.0,  +2.0,  -1.0,  +2.0/  !****  quarks:   U,D,S,C,B,T
       DATA T3F      /   +0.5,  -0.5,  -0.5,  +0.5,  -0.5,  +0.5/  !****  quarks:   U,D,S,C,B,T
       DATA ICHANNEL  /  1 /
C----------------------------------------------------------------------
C PULL MC AND MB FROM QCDNUM USING 
      common /fred/ xmc,xmb,Hmass  !*** PULL VALUES FROM QCDNUM  fio 14 FEB. 2011
C----------------------------------------------------------------------
       polarity = polar

       XMARRAY(4)=XMC           !*** PULL VALUES FROM QCDNUM  fio 14 FEB. 2011
       XMARRAY(5)=XMB           !*** PULL VALUES FROM QCDNUM  fio 14 FEB. 2011
C----------------------------------------------------------------------
C WEINBERG ANGLE:
C----------------------------------------------------------------------
       THETAW=ASIN(SQRT(SINW2))
!       hmass=0.0d0 !********* PATCH FOR TESTING
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
       if (icharge.eq.4) then
          facgz(2) = - eleVec-polarity*eleAxial
          faczz(2) = eleVec**2 + eleAxial**2 +2.0d0 *
     >        polarity*eleAxial*eleVec 
c23456
ch          print*,'faczz(2), facgz(2)!',faczz(2), facgz(2)
       elseif (icharge.eq.5) then
          facgz(2) = - eleVec+polarity*eleAxial 
          faczz(2) = eleVec**2 + eleAxial**2-2.0d0 *
     >        polarity*eleAxial*eleVec 
c23456    
       endif

       facgg(1) = facgg(2) 
       facgz(1) = facgz(2) 
       faczz(1) = faczz(2) 

       facgg(3) = 0.0d0
       if (icharge.eq.4) then
          facgz(3) = + eleAxial + polarity*eleVec
          faczz(3) = - 2.0d0 * (eleVec * eleAxial) -
     $         polarity*(eleVec**2+eleAxial**2)
ch          print*,'faczz(3), facgz(3)!',faczz(3), facgz(3)
       elseif (icharge.eq.5) then
          facgz(3) = - eleAxial + Polarity*eleVec
          faczz(3) = + 2.0d0 * (eleVec * eleAxial) -
     $         polarity*(eleVec**2+eleAxial**2)
       endif

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
             term( I)  = 0.0
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
            call TOTF(XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,IPARTIN,IPARTOUT,
     >   XMU,ISET,Ihad,HMASS,Ischeme,Iflg,Ichannel,XTOT123zz) 
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
      if(icharge.eq.4.or.icharge.eq.5) xz=1.0d0 * idebugZ  !*** idebugZ is debug over-ride

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
      Do j=1,3,1  !*** 3='xcb', 4='123L'
         F123Lxcb(1,j)= XXTOT123(j)    
         F123Lxcb(2,j)= XXCHARM( j)    
         F123Lxcb(3,j)= XXBOTTOM(j)    
      enddo

C----------------------------------------------------------------------
C COMPUTE  FL 
C----------------------------------------------------------------------
      rho=Sqrt(1.0d0+(2.0d0*hmass*xbj/Q)**2)  !*** Get Hmass from /fred/ common block 
      Do i=1,3,1  !*** 3='xcb', 4='123L'
         F123Lxcb(i,4)=rho**2*F123Lxcb(i,2)- 2.0d0*xbj*F123Lxcb(i,1)
      enddo

C----------------------------------------------------
C   ***  Protect small values:
C----------------------------------------------------
      Do j=1,4,1  !*** 3='xcb', 4='123L'
         Do i=1,3,1             !*** 3='xcb', 4='123L'
            if(Abs(F123Lxcb(i,j)).lt.small) F123Lxcb(i,j)=0.0d0
         enddo
      enddo

C
      RETURN
      END 

C----------------------------------------------------------------------
