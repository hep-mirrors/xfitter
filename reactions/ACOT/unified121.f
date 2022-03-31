C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE HADLO123(IMASS,Isch,Ihad,IPARTIN,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XLO123)
      Implicit Double Precision (A-H, O-Z)
      Dimension XLO123(3), XLOLR(-1:1)

       CALL HADLO(IMASS,Isch,Ihad,IPARTIN,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XLOLR)
       CALL HEL2TEN(XBJ,Q,HMASS,XLOLR,XLO123)

      RETURN
      END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE  HADSUB123(IMASS,Isch,Ihad,IPARTIN,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XSUB123)
      Implicit Double Precision (A-H, O-Z)
      Dimension XSUB123( 3), XSUBLR( -1:1)

       CALL HADSUB(IMASS,Isch,Ihad,IPARTIN,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XSUBLR)
       CALL HEL2TEN(XBJ,Q,HMASS,XSUBLR,XSUB123)

      RETURN
      END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE HADNLO123(Isch,Ihad,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XNLO123)
      Implicit Double Precision (A-H, O-Z)
      Dimension XNLO123(3), XNLOLR(-1:1)

       CALL HADNLO(Isch,Ihad,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XNLOLR)
       CALL HEL2TEN(XBJ,Q,HMASS,XNLOLR,XNLO123)

      RETURN
      END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      Double Precision Function TEN2SIG_OLD(XSIGN,ICHARGED,
     >   XBJ,Q,GLLEP,GRLEP,XLEPOL,
     >   HMASS,SBIG,F123)
C-----------------------------------------------------------------------------
C     Compute dSig/dx/dy 
C     Arguments modified: FIO 8/15/02
C      05/07/2007  Include Z-Z, and G-Z terms
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       Dimension F123(3)
       PARAMETER (PI=3.14159265359)

         TEN2SIG= 0.0

C       ICHARGED     =  0 FOR NC;  NON-ZERO FOR CC
C       XSIGN= +1.0  !*** THIS IS FOR      LEPTON SCATTERING.
C       XSIGN= -1.0  !*** THIS IS FOR ANTI-LEPTON SCATTERING.
C  **** ANTI-LEPTON SCATTERING IS OBTAINED BY GLLEP <=> GRLEP _OR_ BY FLIPPING XSIGN
C
C  ***** Z TERMS NOT YET IMPLEMENT HERE !!!!!!!!!!!!!
C

         E1=(SBIG-HMASS**2)/(2.0*HMASS)
         XNU=Q**2/(2.0*HMASS*XBJ)
         Y=XNU/E1
         E2=E1*(1.0-Y)

         IF(Y.GE.1.0) RETURN

         ALPHA=1.0/137.
         EEM2=ALPHA*4.0*PI
         G1PHOTON= EEM2/Q**2

         VEV = 246.0
         WMASS=80.0
         G= (2.0*WMASS)/VEV
         GBW= G/(2.0*SQRT(2.0))
         G1W= GBW**2/(Q**2+WMASS**2)

 
      IF(ICHARGED.LE.0.0) THEN
C        ******** FOR NEUTRAL CURRENT ONLY *********
C  ***** Z TERMS NOT YET IMPLEMENT HERE !!!!!!!!!!!!!
         XFAC= 2.0 * HMASS * E1/(Pi) * G1PHOTON**2 /XLEPOL
       ELSE
C        ******** FOR CHARGED CURRENT ONLY *********
         XFAC= 2.0 * HMASS * E1/(Pi) * G1W**2 /XLEPOL 
C        XFAC= XFAC /(1+Q**2/WMASS**2)**2  !*** remove: FIO 18 mar 2002
       ENDIF

         GPLUS  = GLLEP**2 + GRLEP**2
         GMINUS = GLLEP**2 - GRLEP**2

         TERM1= XBJ * F123(1) * Y*Y 
         TERM2= F123(2) *( (1.-Y)-(HMASS*XBJ*Y)/(2.0*E1) )
         TERM3= XBJ * F123(3) * Y * (1.-Y/2.) 
         SUM= GPLUS * (TERM1+TERM2) + XSIGN * GMINUS * TERM3

         XSEC= XFAC * SUM
         TEN2SIG= XSEC

        RETURN
        END


C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE HEL2TEN(XBJ,Q,HMASS,XHLR,XH123)
C-----------------------------------------------------------------------------
C      Computes F123 HADRON Helicity Amps from CAOT paper
C      
C      
C      
C      
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension XH123(3), XHLR(-1:1)
      PARAMETER (PI=3.14159265359)
      data isave /0/
      save isave

cv==========================
cv set hevay mass thresholds
C      IF(isave.eq.0) THEN
C         isave=1
C         write(6,*) " PATCH FOR TESTING: HMASS=0 "
C      ENDIF
C         hmass=0.0d0
cv==========================


      XPLUS=XHLR( 1)
      XZERO=XHLR( 0)
      XMINU=XHLR(-1)

      RHO=1.0
      IF(HMASS.GT.0.0) THEN 
        XNU=Q**2/(2.0*HMASS*XBJ)
        RHO=SQRT(1.0+(Q/XNU)**2)  
      ENDIF

      XH123(1)=(1./2.)*(XPLUS+XMINU)
      XH123(2)=(XBJ/RHO**2)*(XPLUS+XMINU+2.0*XZERO)
      XH123(3)=(1./RHO)*(-XPLUS+XMINU)

      RETURN
      END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE TEN2HEL_OLD(XBJ,Q,HMASS,XHLR,XH123)
C-----------------------------------------------------------------------------
C      Computes F123 HADRON Helicity Amps from CAOT paper
C      
C      
C      
C      
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension XH123(3), XHLR(-1:1)
      PARAMETER (PI=3.14159265359)


      F1=XH123(1) 
      F2=XH123(2) 
      F3=XH123(3) 

      RHO=1.0
      IF(HMASS.GT.0.0) THEN 
        XNU=Q**2/(2.0*HMASS*XBJ)
        RHO=SQRT(1.0+(Q/XNU)**2)  
      ENDIF

      XPLUS= F1-(RHO/2.0)*F3
      XMINU= F1+(RHO/2.0)*F3
      XZERO=-F1+(RHO**2/(2.0*XBJ))*F2

      XHLR( 1)=XPLUS
      XHLR( 0)=XZERO
      XHLR(-1)=XMINU

      RETURN
      END


C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      FUNCTION FLMAKE_OLD(XBJ,Q,HMASS,XH123)
C-----------------------------------------------------------------------------
C      Computes F123 HADRON Helicity Amps from CAOT paper
C      
C      
C      
C      
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension XH123(3), XHLR(-1:1)
      PARAMETER (PI=3.14159265359)

      CALL TEN2HEL_OLD(XBJ,Q,HMASS,XHLR,XH123)
      FLMAKE= XHLR( 0)

      RETURN
      END


CC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      FUNCTION FMAKER_OLD(XBJ,Q,HMASS,F123)
C-----------------------------------------------------------------------------
C      Computes R-RATIO given F123 HADRON Helicity Amps  
C      
C      
C      
C      
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension F123(3) 
      PARAMETER (PI=3.14159265359)

      FMAKER=0.0

      RHO=1.0
      IF(HMASS.GT.0.0) THEN 
        XNU=Q**2/(2.0*HMASS*XBJ)
        RHO=SQRT(1.0+(Q/XNU)**2)  
      ENDIF
      RHO2=RHO*RHO


      RATIO=0.0
      IF( F123(1).NE.0.0 ) 
     >     RATIO=  F123(2) / ( 2.0*XBJ* F123(1) ) * RHO2-1.0

      FMAKER=RATIO

      RETURN
      END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE HADLO (IMASS,Isch,Ihad,IPARTON,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XHADLO)
C-----------------------------------------------------------------------------
C      Computes LO Hadron Helicities
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C      
C      
C      
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
        Dimension OmGLO(-1:1), XHADLO(-1:1)
        PARAMETER (PI=3.14159265359)
        DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

        Q2=Q**2

        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))
 
        SHATTH= F2M**2
        XITH=ETA*(Q2-F1M**2+SHATTH +
     >           DELTA(-Q2,F1M**2,SHATTH))/(2.0*Q2)

   
         if(isch.eq.9) xith=Min(eta*(1.0d0+(f1m+f2m)**2/Q2),1.0d0) !*** For S-ACOT(CHI)
         ZERO=0.0
         IF(IMASS.EQ.0) THEN
            XITH=XBJ   !*** MASSLESS KINEMATICS 
            if(isch.eq.9) xith=Min(eta*(1.0d0+(f1m+f2m)**2/Q2),1.0d0) !*** For S-ACOT(CHI)
             CALL  PARLO (Q,ZERO,ZERO,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)
         ELSE
             CALL  PARLO (Q,F1M ,F2M ,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)
         ENDIF

        IHADRON= Ihad              !*** 2/28/05 FIO Fixed
        PdfTMP  = PDF(iSet, IHADRON, IPARTON, XITH, XMU, iRet)

        XHADLO( 1) = OmGLO( 1) * PDFTMP 
        XHADLO( 0) = OmGLO( 0) * PDFTMP 
        XHADLO(-1) = OmGLO(-1) * PDFTMP 

        RETURN
        END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE HADLO123chk(IMASS,Isch,Ihad,IPARTON,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XHADLO123)
C-----------------------------------------------------------------------------
C      Computes LO Hadron Helicities
C      
C      01/24/2011  Covert to 123 basis
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C      
C      
C      
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
        Dimension OmGLO123(3), XHADLO123(3)
        PARAMETER (PI=3.14159265359)
        COMMON /PARCOMMON/ XBJout,Qout,HMASSout
        DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

        XBJout=XBJ          !*** Fill common block for (123)->(LR0) Conversion
        Qout=Q
        HMASSout=HMASS

        Q2=Q**2

        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))
 
        SHATTH= F2M**2
        XITH=ETA*(Q2-F1M**2+SHATTH +
     >           DELTA(-Q2,F1M**2,SHATTH))/(2.0*Q2)

   
         if(isch.eq.9) xith=Min(eta*(1.0d0+(f1m+f2m)**2/Q2),1.0d0) !*** For S-ACOT(CHI)
         ZERO=0.0
         IF(IMASS.EQ.0) THEN
            XITH=XBJ   !*** MASSLESS KINEMATICS 
            if(isch.eq.9) xith=Min(eta*(1.0d0+(f1m+f2m)**2/Q2),1.0d0) !*** For S-ACOT(CHI)
             CALL  PARLO123 (Q,ZERO,ZERO,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO123)
         ELSE
             CALL  PARLO123 (Q,F1M ,F2M ,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO123)
         ENDIF

        IHADRON= Ihad              !*** 2/28/05 FIO Fixed
        PdfTMP  = PDF(iSet, IHADRON, IPARTON, XITH, XMU, iRet)

        XHADLO123( 1) = OmGLO123( 1) * PDFTMP 
        XHADLO123( 2) = OmGLO123( 2) * PDFTMP 
        XHADLO123( 3) = OmGLO123( 3) * PDFTMP 

        RETURN
        END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE PARLO123(Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO123)
C-----------------------------------------------------------------------------
C      Computes LO Parton 123 Amps from AOT paper
C      
C      01/24/2011  Covert to 123 basis
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension OmGLO123(3), OmGLO(-1:1)
      PARAMETER (PI=3.14159265359)
      COMMON /PARCOMMON/ XBJ,Qin,HMASS 
             
      CALL  PARLO (Q,F1M ,F2M ,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)
      CALL HEL2TEN(XBJ,Qin,HMASS,OmGLO,OmGLO123)
 
C-----------------------------------------------------------------------------
      Return
      End
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE PARLO (Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)
C-----------------------------------------------------------------------------
C      Computes LO Parton Helicity Amps from AOT paper
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C      
C      
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension OmGLO(-1:1), OHelLO(-1:1,-1:1)
      PARAMETER (PI=3.14159265359)
      Data iLeft, iLong, iRight / -1, 0, 1 /
     >     igL2, igRgL, igR2 / -1, 0, 1 /

      DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))
C-----------------------------------------------------------------------------
      F1m2    = F1m**2
      F2m2    = F2m**2
      Q2      = Q**2
      DEL   = DELTA(-Q2,F1M2,F2M2)
C-----------------------------------------------------------------------------
      OHelLO(iRight, igR2)  =   (Q2+F1m2+F2m2 + DEL)/DEL
      OHelLO(iRight, igRgL) =   -2*F1m*F2m/DEL
      OHelLO(iRight, igL2)  =   (Q2+F1m2+F2m2 - DEL)/DEL

      OHelLO(iLeft, igR2)  =  OHelLO(iRight, igL2) 
      OHelLO(iLeft, igRgL) =  OHelLO(iRight, igRgL)
      OHelLO(iLeft, igL2)  =  OHelLO(iRight, igR2) 

      OHelLO(iLong, igR2)  = ((F1m2+F2m2) + ((F1m2-F2m2)**2/Q2))/DEL
      OHelLO(iLong, igRgL) = +2*F1m*F2m/DEL
      OHelLO(iLong, igL2)  = OHelLO(iLong, igR2)
C-----------------------------------------------------------------------------
C      WR=     GRQ**2 *OHelLO(iRight, igR2)  + 
C     >    2.0*GRQ*GLQ*OHelLO(iRight, igRgL) + 
C     >        GLQ**2 *OHelLO(iRight, igL2)  
C
C      WZ=     GRQ**2 *OHelLO(iLong,  igR2)  + 
C     >    2.0*GRQ*GLQ*OHelLO(iLong,  igRgL) + 
C     >        GLQ**2 *OHelLO(iLong,  igL2)  
C
C      WL=     GRQ**2 *OHelLO(iLeft,  igR2)  + 
C     >    2.0*GRQ*GLQ*OHelLO(iLeft,  igRgL) + 
C     >        GLQ**2 *OHelLO(iLeft,  igL2)  
C-----------------------------------------------------------------------------
c    MODIFY TO ALLOW FOR Z-GAMMA TERMS   5/7/07
C-----------------------------------------------------------------------------
      WR=     GRQ1*GRQ2         *OHelLO(iRight, igR2)  + 
     >    (GRQ1*GLQ2+GRQ2*GLQ1) *OHelLO(iRight, igRgL) + 
     >        GLQ1*GLQ2         *OHelLO(iRight, igL2)  

      WZ=    GRQ1*GRQ2          *OHelLO(iLong,  igR2)  + 
     >   (GRQ1*GLQ2+GRQ2*GLQ1)  *OHelLO(iLong,  igRgL) + 
     >       GLQ1*GLQ2          *OHelLO(iLong,  igL2)  

      WL=    GRQ1*GRQ2          *OHelLO(iLeft,  igR2)  + 
     >   (GRQ1*GLQ2+GRQ2*GLQ1)  *OHelLO(iLeft,  igRgL) + 
     >       GLQ1*GLQ2          *OHelLO(iLeft,  igL2)  

      OmGLO( 1)=WR  
      OmGLO( 0)=WZ  
      OmGLO(-1)=WL  

C-----------------------------------------------------------------------------
      Return
      End
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE HADSUB (IMASS,Isch,IhadIn,IPARTON,XBJ,Q,F1Min,F2Min,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMUin,ISETin,HMASS,XHADSUB)
C-----------------------------------------------------------------------------
C      Computes SUB Hadron Helicities:
C           Gluon initiated only; can not do quark with parallel structure: FIO 2/17/99
C      
C      25 June 98: Moved XHADSUBQ into separate module
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
        Dimension OmGLO(-1:1), XHADSUB(-1:1) 
        EXTERNAL XINTSUB
        PARAMETER (PI=3.14159265359)
        Common  / ActInt /  AERR, RERR, iActL, iActU
        Common  / CINTSUB /  XMU, XiTH, iSet, IPARTONout, Ihad
        DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

        XHADSUB( 1) = 0.0
        XHADSUB( 0) = 0.0
        XHADSUB(-1) = 0.0

C        write(6,*) 'xmuin, f1min, XMUin.LE.F1Min',
C     >   xmuin, f1min, XMUin.LE.F1Min 
        IF(XMUin.LE.F1Min) RETURN

        F1M    =F1Min
        F2M    =F2Min
        IPARTONout=IPARTON
        iSet =ISETin
        XMU  =XMUin
        Ihad =IhadIn

        Q2=Q**2

        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))
 
        SHATTH= F2M**2
        XITH=ETA*(Q2-F1M**2+SHATTH +
     >           DELTA(-Q2,F1M**2,SHATTH))/(2.0*Q2)

       ALPIX = AlQCD(XMU,iset)/pi
       DivFct = ALPIX * Log((XMU/F1m)**2) / 2.0  
C  *** changed from 4.0 to 2.0 and move 2.0 into XINTSUB  FIO 12/7/95

         if(isch.eq.9) xith=Min(eta*(1.0d0+(f1m+f2m)**2/Q2),1.0d0) !*** For S-ACOT(CHI)

         ZERO=0.0
         IF(IMASS.EQ.0) THEN
             XITH=XBJ   !*** MASSLESS KINEMATICS 
             if(isch.eq.9) xith=Min(eta*(1.0d0+(f1m+f2m)**2/Q2),1.0d0) !*** For S-ACOT(CHI)
             CALL  PARLO (Q,ZERO,ZERO,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)
         ELSE
             CALL  PARLO (Q,F1M ,F2M ,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)
         ENDIF

       XIMAX=1.0
       IF(XITH.LT.XIMAX) THEN
          Convol = AdzInt (XINTSUB, XITH, XIMAX, 
     >                 AERR, RERR, ErrEst, iErr, iActL, iActU)
       ELSE
          Convol  = 0.0
       ENDIF

        PDFSUB = DIVFCT * CONVOL


        XHADSUB( 1) = OmGLO( 1) * PDFSUB 
        XHADSUB( 0) = OmGLO( 0) * PDFSUB 
        XHADSUB(-1) = OmGLO(-1) * PDFSUB 

        RETURN
        END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Function XINTSUB(X)
C-----------------------------------------------------------------------------
C The splitting function for G->Q
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
        Common  / CINTSUB /  XMU, XiTH, iSet, IPARTON , Ihad

      iGluon=0
      GluPdf = Pdf(iSet,iHad,iGluon,XiTH/X, XMU, iRet)
      SpltFn = ( X**2 + (1.0-X)**2 )/2.0
C  ***  moved 2.0 into XINTSUB  FIO 12/7/95
      XINTSUB = SpltFn * GluPdf / X

      Return
      End
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE HADNLOG123z(Isch,IhadIn,XBJin,Qin,F1Min,F2Min,
     >   GLQin1,GRQin1,GLQin2,GRQin2,XMUin,ISETin,HMASS,FHAD123G2z)
C-----------------------------------------------------------------------------
C      Computes Hadron Helicities
C      Calls XINTNLOG123z to do convolution
C      
C      14 DECEMBER 1998: IMPLEMENT THE MASSLESS CASE ALA FURMANSKI & PETRONZIO
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       DIMENSION FHAD123G2z(3)
       PARAMETER (PI=3.141592653589793)
       Common  / ActInt /  AERR, RERR, iActL, iActU
       COMMON /CXINTNLOG123z/ 
     >   XBJ,XIMIN,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   Q,XMU,ISET,IRET,INDEX,Ihad
       EXTERNAL XINTNLOG123z
       DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

        F1M = F1Min
        F2M = F2Min
        GLQ1 = GLQin1
        GRQ1 = GRQin1
        GLQ2 = GLQin2
        GRQ2 = GRQin2
        XBJ = XBJin
        Q   = Qin
        XMU = XMUin
        ISET= ISETin
        Ihad=IhadIn


        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))

        XIMIN= xbj   !*** debugging patch (Remove later)
        XIMIN= ETA   !*** Massless Parton Case:
        XIMIN= ETA*((F1M+F2M)**2+Q**2)/Q**2  !*** Full Case:
        SHATTH= (F1M+F2M)**2
        Q2=Q*Q
        XIMIN=ETA*(Q2-0.0**2+SHATTH +
     >           DELTA(-Q2,0.0d0**2,SHATTH))/(2.0*Q2)

        XIMIN= xbj   !*** use massless kinematics
        XIMAX= 1.0

        DO INDEX=1,3,1
        Conv  = 0.0
        IF(XIMIN.LT.XIMAX)  
     >    Conv = AdzInt (XINTNLOG123z, XIMIN, XIMAX, 
     >                   AERR, RERR, ErrEst, iErr, iActL, iActU)
          FHAD123G2z(INDEX) = CONV
        ENDDO

        RETURN
        END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
        Function XINTNLOG123z(XI)
C-----------------------------------------------------------------------------
C      COMPUTES CONVOLUTION OF PARTON HELICITIES WITH PDF'S
C      USED BY HADNLOG1230 PROGRAM
C      
C      14 DECEMBER 1998: IMPLEMENT THE MASSLESS CASE ALA FURMANSKI & PETRONZIO
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       DOUBLE COMPLEX  DiLog 
       DOUBLE COMPLEX CARG,CSPENCE
       PARAMETER (PI=3.141592653589793)
       DATA SMALL /1.E-14/
       COMMON /CXINTNLOG123z/ 
     >   XBJ,XIMIN,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   Q,XMU,ISET,IRET,INDEX,Ihad

        XINTNLOG123z=0.0

        IHADRON=Ihad             !*** 2/28/05 FIO Fixed
        IGLUON=0

        Z = XI 
        Z0 = XIMIN 
        ZHAT = Z0/Z
 
        IF(Z.GE.1.0-SMALL) Z=(1.0-SMALL)

        ALPHAS= AlQCD(XMU,iset)  

C-----------------------------------------------------------------------------
C --- DELTA FUNCTION TERMS: 

C   there are none for the gluon case  

C-----------------------------------------------------------------------------
C --- REGULAR TERMS: 

         FREG2  = 2 * 1/2. * 
     >          ((z**2 + (1-z)**2)*LOG((1-z)/z) - 1 + 8*z*(1-z))

         FREG21 = + 2 * 1/2. * (4 * z * (1-z))


         FREG1= (FREG2 - FREG21)/2.0    
         FREG3= 0.0

C-----------------------------------------------------------------------------
C --- PLUS FUNCTION TERMS: 

C   there are none for the gluon case  

C-----------------------------------------------------------------------------
C --- SUM TERMS:   (Standard F2 is defined with an overall XBJ factor)

        igluon=0
         Pdfz  = PDF(iSet, IHADRON, igluon, ZHAT , XMU, iRet)

          FTEMP1= FREG1*Pdfz/Z
          FTEMP2= FREG2*Pdfz/Z * XBJ
          FTEMP3= FREG3*Pdfz/Z

        IF    (INDEX.EQ. 1)  THEN
             FTEMP=  FTEMP1 * ( GLQ1*GLQ2 + GRQ1*GRQ2 )
           ELSEIF(INDEX.EQ. 2)  THEN
             FTEMP=  FTEMP2 * ( GLQ1*GLQ2 + GRQ1*GRQ2 )
           ELSEIF(INDEX.EQ. 3)  THEN
             FTEMP=  FTEMP3 * ( GLQ1*GLQ2 - GRQ1*GRQ2 )
           ELSE 
             WRITE(6,*) ' BAD INDEX ',INDEX
             CALL HF_ERRLOG(107,'F: XINTNLOG123z - Bad index')
             STOP
        ENDIF

        FTEMP= FTEMP * ALPHAS/(2.0*PI) 
        FTEMP= FTEMP *2.0   !*** TO MATCH OUR NORMALIZATION

        XINTNLOG123z= FTEMP

        RETURN
        END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE HADNLO(Isch,Ihad,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XHLR)
C-----------------------------------------------------------------------------
C      Computes F(LR0) HADRON Helicity Amps from CAOT paper
C      
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C      
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension XH123(3), XHLR(-1:1)
      PARAMETER (PI=3.14159265359)

      XHLR( 1)=FHADHEL( 1,Isch,Ihad,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS)
      XHLR( 0)=FHADHEL( 0,Isch,Ihad,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS)
      XHLR(-1)=FHADHEL(-1,Isch,Ihad,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS)

      RETURN
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Function FHADHEL (IHEL,Isch,IhadIn,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS)
C-----------------------------------------------------------------------------
C      Computes Hadron Helicities
C      Calls XINTNLO to do convolution
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C      
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
       PARAMETER (PI=3.14159265359)
       Common  / ActInt /  AERR, RERR, iActL, iActU
       COMMON /CXINTNLO/ ETA,Qout,F1Mout,F2Mout,
     >   GLQout1,GRQout1,GLQout2,GRQout2,
     >   XMUout,ISETout,IRETout,IHELout,Ihad
       EXTERNAL XINTNLO

        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))

        IHELout=IHEL
        Qout   =Q
        F1Mout =F1M
        F2Mout =F2M
        GLQout1  =GLQ1
        GRQout1  =GRQ1
        GLQout2  =GLQ2
        GRQout2  =GRQ2
        XMUout =XMU
        ISETout=ISET
        Ihad   =Ihadin

        XIMIN= ETA*((F1M+F2M)**2+Q**2)/Q**2
        XIMAX= 1.0

        Conv  = 0.0
        IF(XIMIN.LT.XIMAX)  
     >  Conv = AdzInt (XINTNLO, XIMIN, XIMAX, 
     >                   AERR, RERR, ErrEst, iErr, iActL, iActU)

        FHADHEL = CONV

        RETURN
        END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
        Function XINTNLO(XI)
C-----------------------------------------------------------------------------
C      COMPUTES CONVOLUTION OF PARTON HELICITIES WITH PDF'S
C      USED BY HADHEL PROGRAM
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C      
C      
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       Dimension OmGNLO(-1:1)
       PARAMETER (PI=3.14159265359)
       COMMON /CXINTNLO/ 
     >   ETA,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,IRET,IHEL,Ihad
       COMMON /CXINTNLOTMP/ JSET

        JSET=ISET  !*** THIS IS A PATCH TO PASS ALQCD THE PDF SET: FIO 22SEP99

        IHADRON=1
        IGLUON=0

CCC ******* NOTE, WE USE ETA HERE INSTEAD OF XIMIN BECAUSE THE SLOW-RESCALING
CCC ******* MASS CORRECTION IS BUILT INTO THE MATRIX ELEMENTS AND NEED NOT BE
CCC ******* PUT IN BY HAND

        XHAT = ETA/XI
        CALL PARNLO (XHAT,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,XMU,OmGNLO)
        FTEMP = OmGNLO(IHEL)
        PdfTMP  = PDF(iSet, Ihad, IGLUON, XI, XMU, iRet)
        TEMP = FTEMP * PDFTMP / XI
        XINTNLO= TEMP

        RETURN
        END


C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE PARNLO (X,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,XMU,OmGNLO)
C-----------------------------------------------------------------------------
C      Computes Parton Helicity Amps from CAOT paper
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C      
C      
C      
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension OmGNLO(-1:1)
      PARAMETER (PI=3.14159265359)
C !*** THIS IS A PATCH TO PASS ALQCD THE PDF SET: FIO 22SEP99
      COMMON /CXINTNLOTMP/ ISET

      DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

      OmGNLO( 1)=0.0
      OmGNLO( 0)=0.0
      OmGNLO(-1)=0.0

      ALPHAS= AlQCD(xMu,iset)  

      Q2= Q**2
      F1M2 = F1M ** 2
      F2M2 = F2M ** 2

CCCC   THIS CHECK IS THE SAME AS THE Q CHECK BELOW
      XMAX=Q2/((F1M+F2M)**2+Q2)
      IF((X.GE.1.0).OR.(X.GE.XMAX)) RETURN

CCCC     HMASS=0  !*** TEMP PATCH: THIS IS PARTON LEVEL
CCCC     SMIN=((F1M+F2M+HMASS)**2-HMASS**2)
      SMIN=(F1M+F2M)**2
      QMIN=SQRT( SMIN*X/(1.0-X) )
      IF(Q.LE.QMIN) RETURN

C     Cm-Energy for the hard process for this SHAT=S
      S = Q2 * (1./X - 1.)
      IF(S.LE.SMIN) RETURN
      RS= SQRT(S)
 
      DEL = DELTA(S, F1M2, F2M2)
      TLOG = + Log(4*F1M2*s/(s+F1M2-F2M2+Del)**2) 
      ULOG = + Log(4*F2M2*s/(s-F1M2+F2M2+Del)**2) 

      XLAM =  ( Q2 + S)          /(2.0*RS)
      BET =  DEL                /(2.0*RS)
      E1  =  ( F1M2 - F2M2 + S) /(2.0*RS)
      E2  =  (-F1M2 + F2M2 + S) /(2.0*RS)
      EQ  =  (-Q2 + S)          /(2.0*RS)

      GSPLUS = -((1./2 + e1*(-1. + e1/XLAM)/XLAM)*TLOG) - 
     >   (1./2 + e2*(-1. + e2/XLAM)/XLAM)*ULOG 
     >   - 2*bet*eq**2/(XLAM**2*rs)

      GXPLUS = -(F1M*F2M*(4*bet*rs + (TLOG + ULOG)*(-F1M2 - F2M2 + s)))/
     >   (4*XLAM**2*s)

      GXPLUS = (2.0)* GXPLUS        !*** ERROR IN ORIGINAL FIO (4/19/93)

      GAPLUS = bet*(-F1M2 + F2M2)/(XLAM**2*rs) - 
     >    TLOG*(1./2 + e1*(-1. + e1/XLAM)/XLAM + 
     >    F1M2*(-F1M2 + F2M2)/(2*XLAM**2*s)) + 
     >    ULOG*(1./2 + e2*(-1. + e2/XLAM)/XLAM - 
     >    F2M2*(-F1M2 + F2M2)/(2*XLAM**2*s))

      GSZERO = eq*(-F1M2 + F2M2)*(TLOG - ULOG)/(XLAM**2*rs) + 
     >   bet*((-F1M2 + F2M2)**2 +
     >   Q2*(-F1M2 - F2M2 + 2*Q2))/(XLAM**2*Q2*rs) + 
     >   (TLOG + ULOG)*(-(eq*(F1M2 + F2M2 - (-F1M2 + F2M2)**2/Q2))/
     >   (2*XLAM**2*rs) + F1M2*F2M2/(XLAM**2*s) + 
     >   (F1M2 + F2M2)*(-(-F1M2 + F2M2)**2 + Q2**2 - 2*XLAM**2*s)/
     >   (4*XLAM**2*Q2*s))

      GXZERO = bet*F1M*F2M/(XLAM**2*rs) + F1M*F2M*(TLOG + ULOG)*
     >      (1/Q2 + (-F1M2 - F2M2 + s)/(2*XLAM**2*s))/2

      GXZERO = (2.0)* GXZERO        !*** ERROR IN ORIGINAL FIO (4/19/93)

C-----------------------------------------------------------------------------
C                                                               The Amplitudes
C                                                               --------------
C                                                    The RR helicity amplitude
       ORRGR2 = GSPLUS + GAPLUS
       ORRGLR = GXPLUS 
       ORRGL2 = GSPLUS - GAPLUS
C                                             The LONG-LONG helicity amplitude
       OZZGR2 = GSZERO 
       OZZGLR = GXZERO 
       OZZGL2 = GSZERO 
C                                                    The LL helicity amplitude
       OLLGR2 = GSPLUS - GAPLUS
       OLLGLR = GXPLUS 
       OLLGL2 = GSPLUS + GAPLUS
C-----------------------------------------------------------------------------
C                                                       Assemble the integrand
C                                                       ----------------------
c     XFAC=ALPHAS/PI/2.0
c     OmGNLO( 1)=XFAC*(ORRGR2*GRQ**2 +2.0*ORRGLR*GLQ*GRQ +ORRGL2*GLQ**2)
c     OmGNLO( 0)=XFAC*(OZZGR2*GRQ**2 +2.0*OZZGLR*GLQ*GRQ +OZZGL2*GLQ**2)
c     OmGNLO(-1)=XFAC*(OLLGR2*GRQ**2 +2.0*OLLGLR*GLQ*GRQ +OLLGL2*GLQ**2)
C-----------------------------------------------------------------------------
C                                                       Assemble the integrand
C                                                       ----------------------
      XFAC=ALPHAS/PI/2.0
      OmGNLO( 1)=XFAC*(ORRGR2*GRQ1*GRQ2+ 
     >   ORRGLR*(GLQ1*GRQ2+GLQ2*GRQ1) +ORRGL2*GLQ1*GLQ2)
      OmGNLO( 0)=XFAC*(OZZGR2*GRQ1*GRQ2+ 
     >   OZZGLR*(GLQ1*GRQ2+GLQ2*GRQ1) +OZZGL2*GLQ1*GLQ2)
      OmGNLO(-1)=XFAC*(OLLGR2*GRQ1*GRQ2+ 
     >   OLLGLR*(GLQ1*GRQ2+GLQ2*GRQ1) +OLLGL2*GLQ1*GLQ2)
C
      Return
      End
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE HADNLOQ123(Isch,IhadIn,XBJin,Qin,F1Min,F2Min,
     >  GLQin1,GRQin1,GLQin2,GRQin2,XMUin,HMASS,IPARTin,ISETin,FHAD123Q)
C-----------------------------------------------------------------------------
C      Computes Hadron Helicities
C      Calls XINTNLO to do convolution
C      
C      FOR CHARGED CURRENT V-A ONLY, M1 IS PURELY A REGULATOR
C      ALA GOTTSCHALK, VAN DER BIJ ..., KRAMER ...
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       DIMENSION FHAD123Q(3)
       PARAMETER (PI=3.141592653589793)
       Common  / ActInt /  AERR, RERR, iActL, iActU
       COMMON /CXINTNLOQ123/ 
     >   XBJ,XIMIN,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   Q,XMU,ISET,IRET,INDEX,IPART,Ihad
       EXTERNAL XINTNLOQ123
       DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

        F1M = F1Min
        F2M = F2Min
        GLQ1 = GLQin1
        GRQ1 = GRQin1
        GLQ2 = GLQin2
        GRQ2 = GRQin2
        XBJ = XBJin
        Q   = Qin
        XMU = XMUin
        IPART= IPARTin
        ISET= ISETin
        Ihad= IhadIn


        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))

        SHATTH= F2M**2
        Q2=Q*Q

        XIMIN= xbj   !*** debugging patch (Remove later)
        XIMIN= ETA*((F1M+F2M)**2+Q**2)/Q**2
        XIMIN=ETA*(Q2-F1M**2+SHATTH +
     >           DELTA(-Q2,F1M**2,SHATTH))/(2.0*Q2)
        XIMAX= 1.0

        DO INDEX=1,3,1
        Conv  = 0.0
        IF(XIMIN.LT.XIMAX)  
     >    Conv = AdzInt (XINTNLOQ123, XIMIN, XIMAX, 
     >                   AERR, RERR, ErrEst, iErr, iActL, iActU)
          FHAD123Q(INDEX) = CONV
        ENDDO

        RETURN
        END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
        Function XINTNLOQ123(XI)
C-----------------------------------------------------------------------------
C      COMPUTES CONVOLUTION OF PARTON HELICITIES WITH PDF'S
C      USED BY HADHEL PROGRAM
C      
C      FIO 1998
C      ALA GOTTSCHALK, VAN DER BIJ ..., KRAMER ...
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
      DOUBLE COMPLEX  DiLog 
      DOUBLE COMPLEX CARG,CSPENCE
       Dimension OmGNLO(-1:1)
       PARAMETER (PI=3.141592653589793)
       DATA SMALL /1.E-14/
       COMMON /CXINTNLOQ123/ 
     >   XBJ,XIMIN,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   Q,XMU,ISET,IRET,INDEX,IPART,Ihad

        XINTNLOQ123=0.0

        IHADRON=Ihad              !*** 2/28/05 FIO Fixed
        IGLUON=0

        XHAT = XIMIN/Xi
        XHAT0= XIMIN 
 
        Z = XI  

        IF(Z.GE.1.0-SMALL) Z=(1.0-SMALL)
C       IF (XHAT.GE.1.0-SMALL) XHAT=(1.0-SMALL)
C       IF (XHAT.LE.    SMALL) XHAT=(    SMALL)

        ALPHAS= AlQCD(XMU,iset)  

        Q2=Q*Q
        X0=Q2/(Q2+F2M**2)
        IF(X0.GE.1.0-SMALL) X0=(1.0-SMALL)

        XLOG =  Log( Q2/( F1M**2 * X0 * Z * (1-X0*Z)) )
        XLOG1=  Log( Q2/( F1M**2 * X0     * (1-X0  )) )

C-----------------------------------------------------------------------------
C --- DELTA FUNCTION TERMS: 
  
        CARG= DCMPLX( -X0/(1.-X0) )
        CSPENCE= DILOG(CARG)  
        SPENCE= DBLE(CSPENCE)

         FDEL2= (3./2. + 2.* LOG(1.-XIMIN)) * (XLOG1 - 2.) 
     >  + 1. - PI*PI/3. - 2.* SPENCE  +2.*LOG(1.-X0)

         FDEL21=  (1.-X0)/X0 * LOG(1.-X0)
         FDEL23=  (1.-X0)/X0 * LOG(1.-X0)

         FDEL1= (FDEL2 - FDEL21)/2.
         FDEL3= (FDEL2 - FDEL23) 
 
         XMEASURE= (1.-XIMIN)   !*** INTEGRATION MEASURE

C-----------------------------------------------------------------------------
C --- REGULAR TERMS: 

         FREG2= 2. + Z - (1.-X0)*Z/( 2.*(1.-X0*Z)**2 )
     >         -(1.-X0)*(Z+2)/(1.-X0*Z) + 1./(2.*(1.-X0*Z))
         FREG21= (1.+X0)*Z - 2.*(1.-X0)/(1.-X0*Z) 
     >         + (1.-X0)**2 * Z**2 /( 1.-X0*Z)
         FREG23= 1.+Z-2.*(1.-X0)/(1.-X0*Z) 

         FREG1= (FREG2 - FREG21)/2.
         FREG3= (FREG2 - FREG23) 

C-----------------------------------------------------------------------------
C --- PLUS FUNCTION TERMS: 

         FPLSZ= (1.+ Z**2)/Z  * (XLOG  -2.) 
         FPLS1= (1.+1.**2)/1. * (XLOG1 -2.) 

         Pdfz  = PDF(iSet, IHADRON, IPART, XHAT , XMU, iRet)
         Pdfz0 = PDF(iSet, IHADRON, IPART, XHAT0, XMU, iRet)

         FPLS2 = (FPLSZ*Pdfz - FPLS1*Pdfz0)/(1.-Z)

         FPLS1= (FPLS2 - 0.0)/2.
         FPLS3= (FPLS2 - 0.0) 
C-----------------------------------------------------------------------------
C --- SUM TERMS:   (F2 is defined with an overall XBJ factor)

          FTEMP1= FDEL1/XMEASURE*Pdfz0 + FREG1*Pdfz/Z + FPLS1
          FTEMP2=(FDEL2/XMEASURE*Pdfz0 + FREG2*Pdfz/Z + FPLS2) * XBJ
          FTEMP3= FDEL3/XMEASURE*Pdfz0 + FREG3*Pdfz/Z + FPLS3

        IF    (INDEX.EQ. 1)  THEN
          FTEMP=  FTEMP1 * ( GLQ1*GLQ2 + GRQ1*GRQ2 )
        ELSEIF(INDEX.EQ. 2)  THEN
          FTEMP=  FTEMP2 * ( GLQ1*GLQ2 + GRQ1*GRQ2 )
        ELSEIF(INDEX.EQ. 3)  THEN
          FTEMP=  FTEMP3 * ( GLQ1*GLQ2 - GRQ1*GRQ2 )
        ELSE 
          WRITE(6,*) ' BAD INDEX ',INDEX
          CALL HF_ERRLOG(108,'F: XINTNLOQ123 - Bad index')
          STOP
        ENDIF

        FTEMP= FTEMP * ALPHAS/(2.0*PI) * (4./3.) 
        FTEMP= FTEMP *2.0   !*** TO MATCH OUR NORMALIZATION

        XINTNLOQ123= FTEMP

        RETURN
        END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE HADNLOQ1230(Isch,IhadIn,XBJin,Qin,F1Min,F2Min,
     >  GLQin1,GRQin1,GLQin2,GRQin2,XMUin,HMASS,IPARTin,ISETin,FHAD123Q)
C-----------------------------------------------------------------------------
C      Computes Hadron Helicities
C      Calls XINTNLOQ1230 to do convolution
C      
C      25 JUNE 1998: IMPLEMENT THE MASSLESS CASE ALA FURMANSKI & PETRONZIO
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       DIMENSION FHAD123Q(3)
       PARAMETER (PI=3.141592653589793)
       Common  / ActInt /  AERR, RERR, iActL, iActU
       COMMON /CXINTNLOQ1230/ 
     >   XBJ,XIMIN,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   Q,XMU,ISET,IRET,INDEX,IPART,Ihad
       EXTERNAL XINTNLOQ1230
       DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

        F1M = F1Min
        F2M = F2Min
        GLQ1 = GLQin1
        GRQ1 = GRQin1
        GLQ2 = GLQin2
        GRQ2 = GRQin2
        XBJ = XBJin
        Q   = Qin
        XMU = XMUin
        IPART= IPARTin
        ISET= ISETin
        Ihad =IhadIn


        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))

        XIMIN= xbj   !*** debugging patch (Remove later)
        XIMIN= ETA   !*** Massless Parton Case:
        XIMIN= ETA*((F1M+F2M)**2+Q**2)/Q**2
        SHATTH= F2M**2
        Q2=Q*Q
        XIMIN=ETA*(Q2-F1M**2+SHATTH +
     >           DELTA(-Q2,F1M**2,SHATTH))/(2.0*Q2)

        XIMIN= xbj   !*** use massless kinematics
        if(isch.eq.9) XIMIN= Min(ETA*((F1M+F2M)**2+Q**2)/Q**2,1.0d0)
         XIMAX= 1.0

        DO INDEX=1,3,1
        Conv  = 0.0
        IF(XIMIN.LT.XIMAX)  
     >    Conv = AdzInt (XINTNLOQ1230, XIMIN, XIMAX, 
     >                   AERR, RERR, ErrEst, iErr, iActL, iActU)
          FHAD123Q(INDEX) = CONV
        ENDDO

        RETURN
        END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
        Function XINTNLOQ1230(XI)
C-----------------------------------------------------------------------------
C      COMPUTES CONVOLUTION OF PARTON HELICITIES WITH PDF'S
C      USED BY HADNLOQ1230 PROGRAM
C      
C      25 JUNE 1998: IMPLEMENT THE MASSLESS CASE ALA FURMANSKI & PETRONZIO
C      
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       DOUBLE COMPLEX  DiLog 
       DOUBLE COMPLEX CARG,CSPENCE
       PARAMETER (PI=3.141592653589793)
       DATA SMALL /1.E-14/
       COMMON /CXINTNLOQ1230/ 
     >   XBJ,XIMIN,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   Q,XMU,ISET,IRET,INDEX,IPART,Ihad

        XINTNLOQ1230=0.0

        IHADRON=Ihad           !*** 2/28/05 FIO Fixed
        IGLUON=0

        XHAT = XIMIN/Xi
        XHAT0= XIMIN
 
        Z = XI  
        z0= XIMIN 

        IF(Z.GE.1.0-SMALL) Z=(1.0-SMALL)

        ALPHAS= AlQCD(XMU,iset)  

C-----------------------------------------------------------------------------
C --- DELTA FUNCTION TERMS: 
  
        CARG= DCMPLX( 1.-z0 )
        CSPENCE= DILOG(CARG)  
        SPENCE= DBLE(CSPENCE)

         FDEL2=(-1.0)*(
     >    PI**2/6. + 9./2.*Z0 + 5./4.*Z0**2 + 3.*LOG(1.-Z0)
     >   -(1./2.)* Z0 * (2.+Z0)*LOG(1.-Z0) - (LOG(1.-Z0))**2 - SPENCE
     >    )

         FDEL2=(-1.0)*(
     >   PI*PI/3. + 7.*z0/2. + z0**2 + 3.*Log(1.-z0) - z0*Log(1.-z0) - 
     >   z0**2*Log(1. - z0)/2. - ( Log(1. - z0) )**2 + z0*Log(z0) + 
     >   z0**2*Log(z0)/2. - 2.*SPENCE
     >   )

         FDEL1= (FDEL2 - 0.0)/2.
         FDEL3= (FDEL2 - 0.0) 
 
         XMEASURE= (1.-z0)   !*** INTEGRATION MEASURE

C-----------------------------------------------------------------------------
C --- REGULAR TERMS: 

         FREG2  = 0.0
         FREG21 = (2.0*Z)
         FREG23 = (1.0+Z)

         FREG1= (FREG2 - FREG21)/2.
         FREG3= (FREG2 - FREG23) 

C-----------------------------------------------------------------------------
C --- PLUS FUNCTION TERMS: 

         FPLS2X= ((1.+ Z**2)/(1.-Z) * ( LOG((1.-Z)/Z) -3./4. ) + 
     >   (9.+5.*Z)/4. )

         Pdfz  = PDF(iSet, IHADRON, IPART, XHAT , XMU, iRet)
         Pdfz0 = PDF(iSet, IHADRON, IPART, XHAT0, XMU, iRet)

         FPLS2  = FPLS2X*(Pdfz/Z  - Pdfz0/1.)

         FPLS21 = 0.0
         FPLS23 = 0.0

         FPLS1= (FPLS2 - FPLS21)/2.
         FPLS3= (FPLS2 - FPLS23) 
C-----------------------------------------------------------------------------
C --- SUM TERMS:   (Standard F2 is defined with an overall XBJ factor)

          FTEMP1= FDEL1/XMEASURE*Pdfz0 + FREG1*Pdfz/Z + FPLS1
          FTEMP2=(FDEL2/XMEASURE*Pdfz0 + FREG2*Pdfz/Z + FPLS2) * XBJ
          FTEMP3= FDEL3/XMEASURE*Pdfz0 + FREG3*Pdfz/Z + FPLS3

        IF    (INDEX.EQ. 1)  THEN
          FTEMP=  FTEMP1 * ( GLQ1*GLQ2 + GRQ1*GRQ2 )
        ELSEIF(INDEX.EQ. 2)  THEN
          FTEMP=  FTEMP2 * ( GLQ1*GLQ2 + GRQ1*GRQ2 )
        ELSEIF(INDEX.EQ. 3)  THEN
          FTEMP=  FTEMP3 * ( GLQ1*GLQ2 - GRQ1*GRQ2 )
        ELSE 
          WRITE(6,*) ' BAD INDEX ',INDEX
          CALL HF_ERRLOG(109,'F: XINTNLOQ1230 - Bad index')
          STOP
        ENDIF

        FTEMP= FTEMP * ALPHAS/(2.0*PI) * (4./3.) 
        FTEMP= FTEMP *2.0   !*** TO MATCH OUR NORMALIZATION

        XINTNLOQ1230= FTEMP

        RETURN
        END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE HADSUBQ123(Isch,IhadIn,IPARTONin,XBJin,Q,F1Min,F2Min,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMUin,ISETin,HMASS,XHADSUBQ123)
C-----------------------------------------------------------------------------
C      Computes SUB Hadron Helicities
C      
C      
C      
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C      
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
        Dimension OmGLO(-1:1),OmGLO123(3),XHADSUBQ123(3)
        EXTERNAL XINTSUBQ123
        PARAMETER (PI=3.14159265359)
        Common  / ActInt /  AERR, RERR, iActL, iActU
        Common /CINTSUBQ/ F1M,XBJ, XMU, XiTH, iSet, IPARTON,Ihad
        DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

        XHADSUBQ123(1) = 0.0
        XHADSUBQ123(2) = 0.0
        XHADSUBQ123(3) = 0.0

        IF(XMUin.LE.F1Min) RETURN

        F1M    =F1Min
        F2M    =F2Min
        XBJ    =XBJin
        XMU    =XMUin
        iSet   =ISETin
        IPARTON=IPARTONin
        Ihad   =Ihadin

        Q2=Q**2

        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))
 
        SHATTH= F2M**2
        XIMIN= xbj   !*** debugging patch (Remove later)
        XIMIN= ETA*((F1M+F2M)**2+Q**2)/Q**2
        XIMIN=ETA*(Q2-F1M**2+SHATTH +
     >           DELTA(-Q2,F1M**2,SHATTH))/(2.0*Q2)
        XIMAX= 1.0

        xith=XIMIN

       ALPHAS = AlQCD(XMU,iset)
       IF(ALPHAS.GT.0.5) ALPIX=0.5
       ComFct = ALPHAS/(2.0*PI) 

       XIMAX=1.0
       XIMIN=XITH
 
C----------------------------------------------------------------------
        Conv = 0.0
        IF(XITH.LT.XIMAX) 
     >    Conv = AdzInt (XINTSUBQ123, XIMIN, XIMAX, 
     >                   AERR, RERR, ErrEst, iErr, iActL, iActU)

        Conv = Conv * ComFct
C----------------------------------------------------------------------

        CALL  PARLO (Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)

        XFAC= 1./(1.+(F2M/Q)**2)

        OmGLO123(1)= ( OmGLO(+1) + OmGLO(-1) )/2.
        OmGLO123(2)= ( OmGLO(+1) + OmGLO(-1) + 2.0*OmGLO(0) )*XFAC
        OmGLO123(3)=  -OmGLO(+1) + OmGLO(-1)
C----------------------------------------------------------------------

      RHO=1.0
      IF(HMASS.GT.0.0) THEN 
        XNU=Q**2/(2.0*HMASS*XBJ)
        RHO=SQRT(1.0+(Q/XNU)**2)  
      ENDIF


        XHADSUBQ123(1) = OmGLO123(1) * Conv 
        XHADSUBQ123(2) = OmGLO123(2) * Conv * XBJ    !*** Def of F2
        XHADSUBQ123(3) = OmGLO123(3) * Conv 


        RETURN
        END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Function XINTSUBQ123(Xin)
C-----------------------------------------------------------------------------
C The splitting function for Q->Q
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      DATA SMALL /1.E-14/
        Common /CINTSUBQ/ 
     >   F1Min,XBJin, XMUIN, XiTHIN, iSetIN, IPARTONin,IhadIn

      F1M=F1Min
      XBJ=XBJin
      XMU=XMUin
      XiTH=XiTHin
      iSet=iSetin
      IPARTON=IPARTONin
      Ihad=IhadIN

      X=Xin
      IF(X.GE.1.0-SMALL) X=1.0-SMALL

       Z =X
      Z0=XiTH

C-----------------------------------------------------------------------------
C --- PDF's

      iHadr=Ihad          !*** 2/28/05 FIO Fixed
      QrkPdf = Pdf(iSet,iHadr,IPARTON,XiTH/X, XMU, iRet)
      QrkPdf0= Pdf(iSet,iHadr,IPARTON,XiTH  , XMU, iRet)

C-----------------------------------------------------------------------------
C --- DELTA FUNCTION TERMS: 
  

        XLOG1= LOG( (F1M/XMU)**2 )
C ***MODIFIED     \/
         delSUB= (+1.)*(  2. - 3./2. *XLOG1  - 2.*XLOG1*Log(1. - z0) 
     >          - 2.*Log(1. - z0) 
     >          - 2.*Log(1. - z0)**2
     >          )

         XMEASURE= (1.-XiTH)   !*** INTEGRATION MEASURE

C-----------------------------------------------------------------------------
C --- PLUS FUNCTION TERMS: 

        regSUB=
     > (
     >  (-XLOG1 -1.)*( (  (1.+ z**2)*QrkPdf/X 
     >                   -(1.+1.**2)*QrkPdf0/1.0) * 1./(1.-z) )
     >  +      (-2.)*( (  (1.+ z**2)*QrkPdf/X 
     >                  - (1.+1.**2)*QrkPdf0/1.0) * 
     >                                       Log(1.-z)/(1.-z) )  !*** modified
     > )


      TOTAL = regSUB + delSUB/XMEASURE* QrkPdf0 /1.0

      TOTAL = (4./3.) * TOTAL   !*** COLOR FACTOR

      XINTSUBQ123 = TOTAL

      Return
      End

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C2345678901234567890123456789012345678901234567890123456789012345678901234567890
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       SUBROUTINE KRETZER1_OLD()
C      PROGRAM KRETZER1
C-----------------------------------------------------------------------------
C      Program to COMPUTE DIS HQ PRODUCTION
C      
C      
C      
C      
C      11/05/99 FIO. MODIFY FRONT END 
C      01/01/99 FIO. CORE CODE BY STEFAN KRETZER
C      S. Kretzer, I. Schienbein. Phys.Rev.D58:094035,1998  hep-ph/9805233 
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s
      dimension acotlo(3), rterm(3), vterm(3)
      DIMENSION FHAD123Q(3)
      DIMENSION XMARRAY(6)
      PARAMETER (PI=3.141592653589793)
C-----------------------------------------------------------------------------
C--- ONLY COMMON BLOCK USED
      Common  / ActInt /  AERR, RERR, iActL, iActU
C-----------------------------------------------------------------------------
C--- FOR NEUTRAL CURRENT
C      DATA SBIG,XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
c     >GLLEP,GRLEP,XMU,ISET,HMASS 
C    >   /  98596.,0.1, 10.,1.6,1.6,0.5,0.50.5,0.5,,0.5,0.5,10,1,0.938/
C      DATA   XLEPOL /2.0/
C      DATA IPARTIN, IPARTOUT,SCALE  /  4, 4,  -1  /
C-----------------------------------------------------------------------------
C--- FOR CHARGED CURRENT
      DATA SBIG,XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   GLLEP,GRLEP,XMU,ISET,HMASS 
     >   / 98596.,0.1, 10.,0.5,1.6,1.0,0.0,1.0,0.0,1.0,0.0,10.0,1,0.938/
      DATA   XLEPOL /1.0/
      DATA IPARTIN, IPARTOUT,SCALE  /  3, 4,  -1  /
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
       DATA XMARRAY  /  0.1, 0.1, 0.2, 1.6, 5.0, 170.0/  
       DATA  AERR, RERR, iActL, iActU    /  0.0, 0.001,1,1/  
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C SETUP INTEGRATION PARAMETERS
C----------------------------------------------------------------------
C     RERR=1.E-8   !*** Setting this too small will yield unstable results ***
      WRITE(6,202)  AERR, RERR, iActL, iActU
202   FORMAT(' AERR= ',1PG14.7,' RERR= ',1PG14.7,
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
         GLQ1=0.5
         GRQ1=0.5
         GLQ2=0.5
         GRQ2=0.5
         GLLEP=0.5
         GRLEP=0.5
         XLEPOL=2.0
       ENDIF

       WRITE(6,*)    GLQ1,GRQ1,GLQ2,GRQ2,GLLEP,GRLEP,SBIG,XLEPOL
       WRITE(6,*) '> GLQ1,GRQ1,GLQ2,GRQ2,GLLEP,GRLEP,SBIG,XLEPOL'
       READ( 5,*)    GLQ1,GRQ1,GLQ2,GRQ2,GLLEP,GRLEP,SBIG,XLEPOL

       E=(SBIG-HMASS**2)/(2.0*HMASS) 

C----------------------------------------------------------------------
C SET THE LOOP
C----------------------------------------------------------------------
1     WRITE(6,102) XBJ,Q, E,F1M,F2M,ISET,HMASS,SCALE
102   FORMAT('   XBJ= ',1PG14.7,'     Q= ',1PG14.7,/,
     >       '     E= ',1PG14.7,'   F1M= ',1PG14.7,
     >       '   F2M= ',1PG14.7,'  ISET= ',I3, /,
     >       ' HMASS= ',1PG14.7,' SCALE= ',1PG14.7)
      WRITE(6,*) ' ENTER:  XBJ,Q, E,F1M,F2M,ISET,HMASS,SCALE  '
      READ (5,*)           XBJ,Q, E,F1M,F2M,ISET,HMASS,SCALE
C----------------------------------------------------------------------
C SET THE SCALES
C----------------------------------------------------------------------
* ... kinematics:
      xb = XBJ
      Q2 = Q*Q  
      m1 = F1M
      m2 = F2M
      m1s = m1 * m1
      m2s = m2 * m2
C----------------------------------------------------
C   ***  COMPUTE MU SCALE
C----------------------------------------------------
       IF(SCALE.GT.0.0) THEN
           XMU = SCALE
       ELSEIF(SCALE.LT.0.0) THEN
           XMU = Q * Abs(scale)
       ENDIF
      Q2f = XMU*XMU     
C      if ( Q2f .le. m2s ) Q2f = m2s
      Q2r = Q2f


C----------------------------------------------------
C   ***  call function
C----------------------------------------------------

      Call HADNLOQ123K(Isch,Ihad, XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,XMU
     >                ,HMASS,IPARTIN,ISET,FHAD123Q)


C----------------------------------------------------
C   ***  PRINT OUT
C----------------------------------------------------

       WRITE(6 ,103) E,XBJ,Y,HMASS,ISET,
     >             SCALE,Q,Q*Q,W,XMU,RHO2
     >        ,(   FHAD123Q(N),N=1,3,1),   STOT


C----------------------------------------------------
103    FORMAT(/,
     >   '---------------------------------------------------',/,
     >   ' E    =',1PG14.7,' XBJ  =',1PG14.7,' Y   =',1PG14.7,/,
     >   ' HMASS=',1PG14.7,' JSET =',i4,10x ,' SCALE=',1PG14.7,/,
     >   ' Q    =',1PG14.7,' Q2   =',1PG14.7,' W    =',1PG14.7,/,
     >   ' FMU  =',1PG14.7,' RHO2 =',1PG14.7,/,
     >   '---------------------------------------------------',/,
     >   ' FHAD123Q   (123): ',4(1PG14.7,1X),/,
     >   '---------------------------------------------------'
     >   )
C----------------------------------------------------------------------
      GOTO 1
C----------------------------------------------------------------------

      end


C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE HADNLOQ123K(Isch,Ihad,XBJin,Qin,F1Min,F2Min,
     >  GLQin1,GRQin1,GLQin2,GRQin2,XMUin,HMASS,IPARTin,ISETin,FHAD123Q)
C-----------------------------------------------------------------------------
C      Computes Hadron Helicities
C      Calls XINTNLO to do convolution
C      
C      HEAVY QUARK INITIATED CONTRIBUTIONS TO DEEP INELASTIC STRUCTURE FUNCTIONS.
C      S. Kretzer, I. Schienbein. Phys.Rev.D58:094035,1998  hep-ph/9805233 
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s
      Double Precision nlo,nloq,nlog
      dimension acotlo(3), rterm(3), vterm(3)
      DIMENSION FHAD123Q(3)
      PARAMETER (PI=3.141592653589793)
C----------------------------------------------------------------------
      external subgconv, subqconv, gintegrand, qintegrand
      external PDF
C----------------------------------------------------------------------
      Common  / CFRED /  iSet, IPARTONin, IPARTONout,Ihadout
      Common  / ActInt /  AERR, RERR, iActL, iActU
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s
      common /str/ str 
      common /counter/ i
      common /delta1/ delta1
      common /comsud/ sud
      common /subq/ subq
      common /acotlo/ acotlo
      common /vecax/ Sp, Sm, Rqp, Rqm
C----------------------------------------------------------------------
      DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     INITIALIZATION
C----------------------------------------------------------------------
      DO I=1,3,1
         FHAD123Q(I)=0.0D0
      ENDDO
C----------------------------------------------------------------------
C     TEST FOR Q<F1M,F2M:   FIO   14 NOV 2012
C----------------------------------------------------------------------
      IF((Qin.LT.F1Min).OR.(Qin.LT.F2Min)) RETURN
C----------------------------------------------------------------------

* ... colour factor
      cf = 4.d0 / 3.d0

      iset=isetIN
      ipart=ipartIN
      IPARTONin=ipart
      IPARTONout=99  !*** Not used
      ihadout=ihad

      xb=XBJin
      Q2= Qin*Qin
      m1= F1Min
      m2= F2Min
      m1s=m1*m1
      m2s=m2*m2
     
      Q2f = XMUin*XMUin  
      Q2r = Q2f
C----------------------------------------------------------------------
C     TEST
C     PUT TARGET MASS CORRECTION HERE : 5/7/07
C----------------------------------------------------------------------
        xbj=xb
c        hmass=0.938d0   !*** Don't override 
        q=sqrt(q2)
        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        ETA   = 2.0/(1.0 + Sqrt(Discr))
C        write(6,*) ' eta = ',eta
        xb=xb*eta
C----------------------------------------------------------------------
C----------------------------------------------------------------------


      xi = xb * (1.d0+m2s/Q2)
      chi = xb/(2.d0*Q2) * ( (Q2-m1s+m2s) + delta(m1s,m2s,-Q2) )
      chit= xb/(2.d0*Q2) * ( (Q2-m1s+m2s) - delta(m1s,m2s,-Q2) )
      axb = ( 1.d0 + (m1+m2)**2/Q2 ) * xb
C----------------------------------------------------------------------
C     TEST
C----------------------------------------------------------------------
      IF(CHI.GE.1.0D0) RETURN


C----------------------------------------------------------------------
C*... vector and axial couplings:  General CURRENT PROCESS:
C      v1 = GLQin + GRQin
C      v2 = GLQin + GRQin
C      a1 = GLQin - GRQin
C      a2 = GLQin - GRQin
C----------------------------------------------------------------------
C----------------------------------------------------------------------
*... vector and axial couplings:  General CURRENT PROCESS:
      v1 = GLQin1 + GRQin1
      v2 = GLQin2 + GRQin2
      a1 = GLQin1 - GRQin1
      a2 = GLQin2 - GRQin2
C----------------------------------------------------------------------

      call get_couplings(a1,a2,v1,v2,Sp,Sm,Rqp,Rqm) 

* ... LO amplitudes including kinematic mass effects:
      acotlo(1) = ( Sp*(Q2+m1s+m2s) - 2.d0*Sm*m1*m2 ) 
     >          / 2.d0 / delta(m1s,m2s,-Q2)
      acotlo(2) = Sp * xb * delta(m1s,m2s,-Q2)/Q2 
      acotlo(3) = 2.d0 * Rqp

C      write(6,*) '(acotlo(i),i=1,3,1),xi,chi,chit,axb',
C     >   (acotlo(i),i=1,3,1),xi,chi,chit,axb


* ... real soft and virtual contributions from NLO quark graphs
      call svnlo(rterm,vterm,Sp,Sm,Rqp,Rqm)

C      str = Ctq3Pd(iset, ipart, chi, dsqrt(Q2f),irt)
       str = PDF(iset, Ihad, ipart, chi, dsqrt(Q2f),irt)
     
C       write(6,*) str,chi, ipart,iset, dsqrt(Q2f)

C----------------------------------------------------
* ... for structure function F_i 
      do 10 i = 1, 3

      sud = fq(i)

* ... NLO contributions: gluon AND quark initiated
C      nlog = dinteg(gintegrand,axb,1.d0,1.d-3)
       nlog = 0.0
C **** NEED PROTECT THE UPPER ENDPOINT: iActU= 1 OR 2, NOT 0
C      nloq = dinteg(qintegrand,chi,1.d0,1.d-4)
       nloq= 0.0d0
       if(chi.le.1.d0) then 
         nloq = AdzInt(qintegrand,chi,1.d0,
     >              AERR, RERR, ErrEst, iErr, iActL, iActU)
       else
         Write(6,*) " error: integration limits >1 in HADNLOQ123K"
         call HF_errlog(105,'W: HADNLOQ123K - integration limits > 1')
       endif

* ... sum of quark/gluon convolutions and soft/virtual terms
      nlo = nlog + acotlo(i) * ( nloq
     >    + cf * (rterm(i)+vterm(i)) * str
* ... including the Sudakov log from the 1/(1-xi)+ distibution  
     >    + cf * sud * str * dlog(1.d0-chi) )

      nloq = nlo - nlog

comment: F_i's normalized to +acotlo(i)*s(chi)+O(alpha_s) for all i
      alphas= alQCD(Sqrt(Q2r),iset)
      nloq  = nloq  * alphas/2.d0/PI

      FHAD123Q(i)=   nloq

C----------------------------------------------------
C  *** ADJUST NORMALIZATION CONVENTIONS TO MATCH OURS: 
      IF(I.EQ.2)  FHAD123Q(i) =  FHAD123Q(i)/(1.0+M2S/Q2) 
C----------------------------------------------------

10    continue  

C----------------------------------------------------

C      write(6,*)  'FHAD123Q', (FHAD123Q(i),i=1,3,1)

      return
      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function subgconv(x)
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s
      external PDF
      Common  / CFRED /  iSet, IPARTONin, IPARTONout,Ihad
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s
C      gl =  Ctq3Pd(1,    0, x, dsqrt(Q2f),irt)
       gl  = PDF(iSet, Ihad, 0, x, dsqrt(Q2f),irt)

      subgconv = 1.d0/x * gl * pqg(xi/x)
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function subqconv(x)
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s
      external PDF
      Common  / CFRED /  iSet, IPARTONin, IPARTONout,Ihad
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s
      common /str/ str
C      sb =  Ctq3Pd(1,           3, chi/x, dsqrt(Q2f), irt)
       sb  = PDF(iSet, Ihad,IPARTONin, chi/x ,dsqrt(Q2f), irt)

      subqconv = 1.d0/x * ( (1.d0+x**2)/(1.d0-x) 
     >         * ( dlog(Q2f/m1s) - 1.d0-2.d0*dlog(1.d0-x) )
     >         * ( sb - x*str ) )

      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function gintegrand(xhat)
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s
      external PDF
      Common  / CFRED /  iSet, IPARTONin, IPARTONout,Ihad
      common /counter/ i
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s

C      gl =  xhat * Ctq3Pd (1,    0, xhat, dsqrt(Q2f), irt)
       gl  =  xhat * PDF(iSet, Ihad, 0, xhat, dsqrt(Q2f), irt)

      if ( i .ne. 2 ) gl = gl / xhat  
      gintegrand = 1.d0/xhat * fsub(xb/xhat,Q2,i,m1,m2) * gl

      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function fsub(x,Q2,i,m1,m2)
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s

      if ( i .eq. 1) fsub = f1sub(x,Q2,m1,m2)
      if ( i .eq. 2) fsub = f2sub(x,Q2,m1,m2)
      if ( i .eq. 3) fsub = f3sub(x,Q2,m1,m2)

      fsub = 2.d0 * fsub

      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function f1sub(z,Q2,m1,m2)
      implicit double precision (a-z)
      integer k
      dimension fp(3)
      common /vecax/ Sp, Sm, Rqp, Rqm

      mi = m1
      mo = m2
c      qp = 2.d0
c      qm = 0.d0
      qp = Sp
      qm = Sm

      v = sqrt(1.d0-(mi+mo)**2/q2*z/(1.d0-z))
      vb = sqrt(1.d0-(mi-mo)**2/q2*z/(1.d0-z))

      fp(1)=-qp/4.d0*v*vb*((1.d0-2.d0*z)**2+qm/qp*mi*mo/q2
     >     * 4.d0*z*(1.d0-z))
      do 20 k = 2, 3
      Ll = dlog( (1.d0+(mi**2-mo**2)/q2*z/(1.d0-z)+v*vb)
     >  /     (1.d0+(mi**2-mo**2)/q2*z/(1.d0-z)-v*vb) )
      fp(k)=qp/4.d0*(.5d0-z*(1.d0-z)+(mi**2-mo**2)/q2*z
     >     * (1.d0-2.d0*z)
     > +(mi**2-mo**2)**2/q2**2*z**2+qm/qp*mi*mo/q2*2.d0*z
     > *(1.d0-z-z*(mi**2+mo**2)/q2))*Ll
      mi = m2
      mo = m1
  20  continue

      f1sub = fp(1) + fp(2) + fp(3)

      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function f2sub(z,Q2,m1,m2)
      implicit double precision (a-z)
      integer  k
      dimension fp(3)
      common /vecax/ Sp, Sm, Rqp, Rqm

      mi = m1
      mo = m2
c      qp = 2.d0
c      qm = 0.d0
      qp = Sp
      qm = Sm

      v = sqrt(1.d0-(mi+mo)**2/q2*z/(1.d0-z))
      vb = sqrt(1.d0-(mi-mo)**2/q2*z/(1.d0-z))

      fp(1)=qp*z*(v*vb*(-.5d0+4.d0*z*(1.d0-z)-((mi**2+mo**2)/q2-
     > (mi**2-mo**2)**2/q2**2)*z*(1.d0-z)))
      do 20 k = 2, 3
      Ll = dlog((1.d0+(mi**2-mo**2)/q2*z/(1.d0-z)+v*vb)
     >  /     (1.d0+(mi**2-mo**2)/q2*z/(1.d0-z)-v*vb))
       fp(k)=qp*z*.25d0*(1.d0-2.d0*z*(1.d0-z)+mi**2/q2
     > *(1.d0+8.d0*z-18.d0*z**2)+mo**2/q2*(1.d0-4.d0*z+6.d0*z**2)-qm/qp
     > * 2.d0*mi*mo/q2
     > -(mi**4+mo**4)/q2**2*2.d0*z*(1.d0-3.d0*z)+mi**2*mo**2/q2**2
     > *4.d0*z
     > *(1.d0-5.d0*z)+(mi**6-mi**4*mo**2-mi**2*mo**4+mo**6)/q2**3
     > *2.d0*z**2)*Ll
      mi = m2
      mo = m1
  20  continue

      f2sub = fp(1) + fp(2) + fp(3)

      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function f3sub(z,Q2,m1,m2)
      implicit double precision (a-z)
      integer k
      dimension fp(3)
      common /vecax/ Sp, Sm, Rqp, Rqm

      mi = m1
      mo = m2
      qp = Sp
      qm = Sm

      v = sqrt(1.d0-(mi+mo)**2/q2*z/(1.d0-z))
      vb = sqrt(1.d0-(mi-mo)**2/q2*z/(1.d0-z))

      fp(1)=v*vb*(mi**2-mo**2)/q2*2.d0*z*(1.d0-z)
      do 20 k = 2, 3
      Ll = dlog((1.d0+(mi**2-mo**2)/q2*z/(1.d0-z)+v*vb)
     >  /     (1.d0+(mi**2-mo**2)/q2*z/(1.d0-z)-v*vb))
       fp(k)=-(.5d0-z*(1.d0-z)+(mi**2-mo**2)/q2*z*(1.d0-2.d0*z)
     > -(mi**4-mo**4)/q2**2*z**2)*Ll
      mi = m2
      mo = m1
  20  continue

      f3sub = -( fp(1) + fp(2) - fp(3) )
      f3sub = Rqp * f3sub

      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c ... g -> q  splitting function at LO:
      function pqg(x)
      implicit double precision (a-z)
      pqg = .5d0 * ( x**2+(1.d0-x)**2 )
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function qintegrand(xhat)
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s
      external PDF
      Common  / CFRED /  iSet, IPARTONin, IPARTONout,Ihad
      common /counter/ i
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s
      common /str/ str
      common /comsud/ sud
      common /subq/ subq
      common /delta1/ delta1

c ... partonic Mandelstam s=(p_quark+q)**2 
      sh = shat(xhat)  

c.......................................................................

C      sb =  Ctq3Pd(1,            3, chi/xhat, dsqrt(Q2f), irt)
       sb  = PDF(iSet, Ihad, IPARTONin, chi/xhat, dsqrt(Q2f), irt)

      qintegrand = 4.d0/3.d0 /xhat * 1.d0/(1.d0-xhat) 
     >           * (  sb*hq1(xhat,sh) - xhat*str*sud  ) 

      if ( i .eq. 2 ) qintegrand = qintegrand
     >                           + 4.d0/3.d0/xhat 
     >                           * sb * hq21(xhat,sh) 

      if ( i .eq. 3 ) qintegrand = qintegrand
     >                           + 4.d0/3.d0/xhat 
     >                           * sb * hq31(xhat,sh)

* .... alternatively subtract convolution part here:
c      qintegrand = qintegrand - subqconv(xhat)*4.d0/3.d0

      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function hq1(xhat,s)
      implicit double precision (a-z)
      dimension acotlo(3)
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s
      common /powers/ Q4,Q6,Q8,m14,m16,m18,m24,m26,m28
      common /acotlo/ acotlo
      common /vecax/ Sp, Sm, Rqp, Rqm
      common /delta1/ delta1

      Q4 = Q2 * Q2
      Q6 = Q4 * Q2
      Q8 = Q6 * Q2
      m14  = m1s * m1s
      m16  = m14 * m1s
      m18  = m16 * m1s
      m24  = m2s * m2s
      m26  = m24 * m2s
      m28  = m26 * m2s
      delta2 = delta(m1s,s,-Q2)
      delta12 = delta1 * delta1
      delta22 = delta2 * delta2
      spp = Q2 + m2s + m1s
      smm = Q2 - m2s - m1s
      spm = Q2 + m2s - m1s
      smp = Q2 - m2s + m1s
      diff = m1s+s+Q2-delta2
      sum  = m1s+s+Q2+delta2
      LN = dlog(diff/sum) 
    
      hq1 = (1.d0-xhat)/8.d0 * (s-m2s)/s *
     > ((8.d0*(-(delta12*(m2s/(s-m2s)**2+s/(s-m2s)**2+
     > (LN*s*spp)/(delta2*(s-m2s)**2))*
     > (-2.d0*m1*m2*Sm+Sp*spp))+
     > 2.d0*m1*m2*Sm*((s*((s-m2s)+2.d0*spm))/(s-m2s)+
     > (delta12+(-3.d0*m1s+2.d0*m2s+Q2)*(s-m2s)+2.d0*m2s*spm)/(s-m2s)+
     > (LN*s*(3.d0*delta12+(s-m2s)**2-4.d0*m1s*smp+(s-m2s)*
     > (m1s+3.d0*spm)))/(delta2*(s-m2s))+
     > ((m2s+(s-m2s))*((s-m2s)+spp))/(2.d0*s))+
     > Sp*(4.d0*m14+2.d0*m1s*(s-m2s)-spm*(m2s+spm)-
     > ((delta12+2.d0*m2s*spm)*spp)/(s-m2s)-
     > (s*spm*((s-m2s)+2.d0*spp))/(s-m2s)+
     > (((s-m2s)+spp)*(delta12-4.d0*m1s*(s-m2s)+
     > (s-m2s)**2-2.d0*m2s*spp))/(4.d0*s)+
     > (LN*s*(-(s-m2s)**3-4.d0*(s-m2s)**2*spm+(s-m2s)*
     > (4.d0*m1s*m2s-7.d0*spm*spp)+
     > 2.d0*spp*(-delta12-2.d0*spm*spp)))/
     > (2.d0*delta2*(s-m2s)))))/delta22)

      hq1 = hq1 / acotlo(1)

      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function hq21(xhat,s)
      implicit double precision (a-z)
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s
      common /powers/ Q4,Q6,Q8,m14,m16,m18,m24,m26,m28
      common /delta1/ delta1 
      common /vecax/ Sp, Sm, Rqp, Rqm
      
      delta2 = delta(m1s,s,-Q2)
      delta12 = delta1  * delta1
      delta14 = delta12 * delta12
      delta22 = delta2  * delta2
      delta24 = delta22 * delta22
      spp = Q2 + m2s + m1s
      smm = Q2 - m2s - m1s
      spm = Q2 + m2s - m1s
      smp = Q2 - m2s + m1s
      spm2 = spm * spm
      smp2 = smp * smp
      diff = m1s+s+Q2-delta2
      sum  = m1s+s+Q2+delta2
      LN = dlog(diff/sum)

       int2 = 
     > ((16.d0*(-2.d0*delta14*Sp*(m2s/(s-m2s)**2+s/(s-m2s)**2+
     > (LN*s*spp)/(delta2*(s-m2s)**2))+
     > 2.d0*m1*m2*Sm*((LN*s*((s-m2s)**3+2.d0*(s-m2s)**2*spm+
     > ((s-m2s)*(-delta12+3.d0*spm2))/2.d0))/
     > (delta2*(s-m2s))-((s-m2s)*(-2.d0*(m1s-m2s+2.d0*Q2)+(s-m2s)+
     > (-delta12+3.d0*smp2)/(2.d0*(s-m2s)))*((s-m2s)+spp))/(2.d0*s)+
     > 2.d0*(delta12-(2.d0*m1s-2.d0*m2s+Q2)*(s-m2s)+
     > (s-m2s)**2-3.d0*Q2*spp))+
     > Sp*(-2.d0*(m1s+m2s)*(s-m2s)**2-9.d0*m2s*spm2-
     > (2.d0*delta12*(delta12+2.d0*m2s*spm))/(s-m2s)+
     > 2.d0*(s-m2s)*(2.d0*delta12+(m1s-5.d0*m2s)*spm)+
     > (s*(-4.d0*delta12*spm+(s-m2s)*(delta12-3.d0*spm2)))/(s-m2s)+
     > ((s-m2s)*(-2.d0*(m1s-m2s+2.d0*Q2)+(s-m2s)+(-delta12+
     > 3.d0*smp2)/(2.d0*(s-m2s)))*spp*
     > ((s-m2s)+spp))/(2.d0*s)+delta12*(-m2s+2.d0*spp)+
     > (LN*s*(2.d0*delta12*(-3.d0*delta12+4.d0*m1s*smp)-
     > 4.d0*(s-m2s)**2*(m1s*(-m1s+m2s)+spm2)-(s-m2s)**3*spp+
     > (s-m2s)*(delta12*(8.d0*m1s-7.d0*spp)+18.d0*m1s*Q2*spp)))/
     > (delta2*(s-m2s)))))/delta24)

       int1 = 
     > ((8.d0*(-(delta12*(m2s/(s-m2s)**2+s/(s-m2s)**2+
     > (LN*s*spp)/(delta2*(s-m2s)**2))*
     > (-2.d0*m1*m2*Sm+Sp*spp))+
     > 2.d0*m1*m2*Sm*((s*((s-m2s)+2.d0*spm))/(s-m2s)+
     > (delta12+(-3.d0*m1s+2.d0*m2s+Q2)*(s-m2s)+2.d0*m2s*spm)/(s-m2s)+
     > (LN*s*(3.d0*delta12+(s-m2s)**2-4.d0*m1s*smp+(s-m2s)*
     > (m1s+3.d0*spm)))/(delta2*(s-m2s))+
     > ((m2s+(s-m2s))*((s-m2s)+spp))/(2.d0*s))+
     > Sp*(4.d0*m14+2.d0*m1s*(s-m2s)-spm*(m2s+spm)-
     > ((delta12+2.d0*m2s*spm)*spp)/(s-m2s)-
     > (s*spm*((s-m2s)+2.d0*spp))/(s-m2s)+
     > (((s-m2s)+spp)*(delta12-4.d0*m1s*(s-m2s)+
     > (s-m2s)**2-2.d0*m2s*spp))/(4.d0*s)+
     > (LN*s*(-(s-m2s)**3-4.d0*(s-m2s)**2*spm+(s-m2s)*
     > (4.d0*m1s*m2s-7.d0*spm*spp)+
     > 2.d0*spp*(-delta12-2.d0*spm*spp)))/
     > (2.d0*delta2*(s-m2s)))))/delta22)

      hq21 = (s-m2s) / (8.d0*s) *
     >     ( delta2**2/2.d0/Sp/delta1           * int2
     >     - 2.d0*delta1/(Sp*spp-2.d0*Sm*m1*m2) * int1 )

      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function hq31(xhat,s)
      implicit double precision (a-z)
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s
      common /powers/ Q4,Q6,Q8,m14,m16,m18,m24,m26,m28
      common /delta1/ delta1 
      common /vecax/ Sp, Sm, Rqp, Rqm
      
      delta2 = delta(m1s,s,-Q2)
      delta12 = delta1 * delta1
      delta22 = delta2 * delta2
      spp = Q2 + m2s + m1s
      smm = Q2 - m2s - m1s
      spm = Q2 + m2s - m1s
      smp = Q2 - m2s + m1s
      diff = m1s+s+Q2-delta2
      sum  = m1s+s+Q2+delta2
      LN = dlog(diff/sum)

      int3 = 
     > ((16.d0*(2.d0*m1*m2*Rqm*(1.d0-smp/(s-m2s)+(LN*s*((s-m2s)+spm))/
     > (delta2*(s-m2s)))-
     > 2.d0*delta12*Rqp*(m2s/(s-m2s)**2+s/(s-m2s)**2+(LN*s*spp)/
     > (delta2*(s-m2s)**2))+
     > Rqp*(2.d0*(m1s-m2s)-(2.d0*s*spm)/(s-m2s)-
     > (2.d0*(delta12+m2s*spm))/(s-m2s)+
     > (LN*s*(-(s-m2s)**2+4.d0*(-delta12+m1s*smp)-
     > 3.d0*(s-m2s)*spm))/(delta2*(s-m2s))-
     > ((s-m2s)*(1.d0-smp/(s-m2s))*((s-m2s)+spp))/(2.d0*s))))/delta22)

      int1 = 
     > ((8.d0*(-(delta12*(m2s/(s-m2s)**2+s/(s-m2s)**2+
     > (LN*s*spp)/(delta2*(s-m2s)**2))*
     > (-2.d0*m1*m2*Sm+Sp*spp))+
     > 2.d0*m1*m2*Sm*((s*((s-m2s)+2.d0*spm))/(s-m2s)+
     > (delta12+(-3.d0*m1s+2.d0*m2s+Q2)*(s-m2s)+2.d0*m2s*spm)/(s-m2s)+
     > (LN*s*(3.d0*delta12+(s-m2s)**2-4.d0*m1s*smp+(s-m2s)*
     > (m1s+3.d0*spm)))/(delta2*(s-m2s))+
     > ((m2s+(s-m2s))*((s-m2s)+spp))/(2.d0*s))+
     > Sp*(4.d0*m14+2.d0*m1s*(s-m2s)-spm*(m2s+spm)-
     > ((delta12+2.d0*m2s*spm)*spp)/(s-m2s)-
     > (s*spm*((s-m2s)+2.d0*spp))/(s-m2s)+
     > (((s-m2s)+spp)*(delta12-4.d0*m1s*(s-m2s)+
     > (s-m2s)**2-2.d0*m2s*spp))/(4.d0*s)+
     > (LN*s*(-(s-m2s)**3-4.d0*(s-m2s)**2*spm+(s-m2s)*
     > (4.d0*m1s*m2s-7.d0*spm*spp)+
     > 2.d0*spp*(-delta12-2.d0*spm*spp)))/
     > (2.d0*delta2*(s-m2s)))))/delta22)

C ***  PATCH: FIO 5/19/99: IF Rqp=0, F3 SHOULD BE ZERO
      hq31 = 0.0
      IF(Rqp.NE.0.0)   hq31 = (s-m2s) / (8.d0*s) *
     >     ( delta2/2.d0/Rqp                    * int3
     >     - 2.d0*delta1/(Sp*spp-2.d0*Sm*m1*m2) * int1 )

      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function fq(i)
      implicit double precision (a-z)
      integer i
      dimension acotlo(3)
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s
      common /delta1/ delta1
      common /acotlo/ acotlo

      fq = 2.d0 /delta1
     >   *  ( (m1s+m2s+Q2)*dlog((m1s+m2s+Q2+delta1)/(m1s+m2s+Q2-delta1))
     >   - 2.d0*delta1  )

      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine get_couplings(a1,a2,v1,v2,Sp,Sm,Rqp,Rqm)
      implicit double precision (a-z)
      Sp  =  v1*v2 + a1*a2
      Sm  =  v1*v2 - a1*a2
      Rqp = (a1*v2 + a2*v1) / 2.d0
      Rqm = (a1*v2 - a2*v1) / 2.d0
      return
      end 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function delta(a,b,c)
      implicit double precision (a-z)
      delta = dsqrt( a**2 + b**2 + c**2 - 2.d0 * ( a*b + a*c + b*c ) )
      return
      end 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function shat(xp)
      implicit double precision (a-z)
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s
      shat = Q2/xp * chi/xb - m1s*xp * xb/chi + m1s-Q2
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function that(xp,zp)
      implicit double precision (a-z)
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s
      that = ( shat(xp)-m1s+Q2 ) 
     >     * ( zp - ( 1.d0 + 2.d0*m1s/(shat(xp)-m1s+Q2) ) ) + m1s
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine svnlo(rterm,vterm,Sp,Sm,Rqp,Rqm)
      implicit double precision (a-z)
      INTEGER iErr, iActL, iActU
      dimension rterm(3), vterm(3)
      external helpf
      Common  / ActInt /  AERR, RERR, iActL, iActU
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s
      common /delta1/ delta1

      pi2 = dacos(-1.d0)*dacos(-1.d0)
      delta1 = delta(m1s,m2s,-Q2)
      sum  = m1s+m2s+Q2+delta1
      diff = m1s+m2s+Q2-delta1
      ln1  = dlog(sum/diff)

c ... arguments of dilogs as lower bounds for integral evaluation
      lb1 = 2.d0*delta1/sum
      lb2 = -diff/2.d0/delta1

C *** FRED: PATCH IN ADZINT  19 MAY 1999
C *** FRED: ensure integration limits a<b are satisfied: 02 may 2005
C     dilog1 = dinteg(helpf,lb1,0.d0,1.d-4)
C     dilog2 = dinteg(helpf,lb2,0.d0,1.d-4)
      dilog1 =-AdzInt(helpf,0.d0,lb1,AERR,RERR,ErrEst,iErr,iActL,iActU)
      dilog2 = AdzInt(helpf,lb2,0.d0,AERR,RERR,ErrEst,iErr,iActL,iActU)

c ... soft real contributions
      rterm(1) = 2.d0 + (m1s+m2s+Q2)/delta1 * ( ln1 - dilog1 - dilog2
     >         - 0.5d0 * (dlog(diff/2.d0/delta1))**2 - pi2/6.d0 )
     >         + dlog(Q2*(chi-chit)**2/m2s/xb**2)
     >         * ( (m1s+m2s+Q2)/delta1*ln1 - 2.d0 )
c ... app is m1->0 limit of rterm    
c      app = (dlog((Q2+m2s)**2/m2s/m1s)-2.d0)*dlog((Q2+m2s)**2/Q2/m2s)
c     >    +2.d0+dlog((Q2+m2s)**2/m1s/m2s)-pi2/3.d0-0.5d0
c     >    *(dlog((Q2+m2s)**2/m1s/m2s))**2
c      write(6,*) rterm(1), app
      rterm(2) = rterm(1)
      rterm(3) = rterm(1)

c ... virtual contributions 
      lb3 = ( delta1 - (Q2-m1s+m2s) ) / 2.d0/delta1
      lb4 = ( delta1 +  Q2+m1s-m2s  ) / 2.d0/delta1
      lb5 = ( delta1 +  Q2-m1s+m2s  ) / 2.d0/delta1
      lb6 = ( delta1 - (Q2+m1s-m2s) ) / 2.d0/delta1

C23456789012345678901234567890123456789012345678901234567890123456789012
C *** FRED: PATCH IN ADZINT  19 MAY 1999
C *** FRED: ensure integration limits a<b are satisfied: 02 may 2005
C      dilog3 = dinteg(helpf,lb3,0.d0,1.d-4)
C      dilog4 = dinteg(helpf,lb4,0.d0,1.d-4)
C      dilog5 = dinteg(helpf,lb5,0.d0,1.d-4)
C      dilog6 = dinteg(helpf,lb6,0.d0,1.d-4)

      dilog3 =-AdzInt(helpf,0.d0,lb3,AERR,RERR,ErrEst,iErr,iActL,iActU)
      dilog4 =-AdzInt(helpf,0.d0,lb4,AERR,RERR,ErrEst,iErr,iActL,iActU)
      dilog5 =-AdzInt(helpf,0.d0,lb5,AERR,RERR,ErrEst,iErr,iActL,iActU)
      dilog6 =-AdzInt(helpf,0.d0,lb6,AERR,RERR,ErrEst,iErr,iActL,iActU)




      lns1 = ( dlog( ( delta1 - (Q2-m1s+m2s) ) / 2.d0/Q2 ) )**2
      lns2 = ( dlog( ( delta1 +  Q2+m1s-m2s  ) / 2.d0/Q2 ) )**2
      lns3 = ( dlog( ( delta1 +  Q2-m1s+m2s  ) / 2.d0/Q2 ) )**2
      lns4 = ( dlog( ( delta1 - (Q2+m1s-m2s) ) / 2.d0/Q2 ) )**2
      c0gl5 = 1/delta1 * ln1*( (Q2+m1s+m2s)*(dlog(Q2/delta1)+1.d0 )
     >      + delta1**2/2.d0/Q2 )
     >      + (Q2+m1s+m2s)/delta1 * ( 0.5d0 * ( lns1-lns2-lns3+lns4 )
     >      - dilog3 + dilog4 + dilog5 - dilog6 )
     >      - 0.5d0 * (m1s-m2s)/Q2 * dlog(m1s/m2s) + dlog(m1s*m2s/Q2/Q2)
     >      - 4.d0
c ... app is m1->0 limit of c0gl5 
c      app = dlog((Q2+m2s)/m1s)*(0.5d0-dlog((Q2+m2s)/Q2))
c     >    + 0.5d0*(dlog((Q2+m2s)/m1s))**2+dlog((Q2+m2s)/m2s)
c     >    *(0.5d0+m2s/Q2-dlog((Q2+m2s)/Q2))+0.5d0
c     >    *(dlog((Q2+m2s)/m2s))**2+2.d0*dlog((Q2+m2s)/Q2)
c     >    +2.d0*dinteg(helpf,Q2/(Q2+m2s),0.d0,1.d-3)-4.d0
c      write(6,*) app, c0gl5

      c0gr5 = 2.d0*m1*m2/delta1 * ln1
  
      c0p1l5 = -1.d0/Q2*( (Q2-m1s+m2s)/delta1 * ln1
     >       +             dlog(m1s/m2s) )

  
      c0p1r5 = -1.d0/Q2*( (Q2+m1s-m2s)/delta1 * ln1
     >       -             dlog(m1s/m2s) )


      vterm(1) = c0gl5 
     >         + ( Sm*(m1s+m2s+Q2) - 2.d0*Sp*m1*m2 )
     >         / ( Sp*(m1s+m2s+Q2) - 2.d0*Sm*m1*m2 ) * c0gr5 

      vterm(2) = c0gl5 + (m1s * c0p1r5 + m2s * c0p1l5)/2.d0
     >         + Sm/Sp * ( c0gr5 + m1*m2/2.d0*(c0p1l5+c0p1r5) )

C ***  PATCH: FIO 5/19/99: IF Rqp=0, F3 SHOULD BE ZERO
      vterm(3) = 0.0
      IF(Rqp.NE.0.0)   vterm(3) = c0gl5 + Rqm/Rqp * c0gr5

      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c ... help-function for evaluation of dilogs
      function helpf(x)
      implicit double precision (a-z)
      helpf = dlog(1.d0-x)/x
      return
      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


