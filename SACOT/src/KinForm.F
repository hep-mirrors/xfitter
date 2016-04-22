      Subroutine SFnConv (SFnIred, g2chi, SFnTen, SFnHel)
      
C                  ------------------------------
C     Given the Ireducible set of structure functions,
C       Parity conserving: 
C              transverse/longitudinal, chirality symmetric/mixed SFnIred(1,2,3,4)
C       Parity violating, transverse, chirality anti-symmetric    SFnIred(5) ;
C     and the chiral coupling combinations,
C       g2chi(1,2,3)             (symmetric/mixed/anti-symmetric) ,
C         1 : gLq(1)*gLq(2) + gRq(1)*gRq(2) 
C         2 : gLq(1)*gRq(2) + gRq(1)*gLq(2)       ! (1,2) refer to the two VecBoson vertices
C         3 : gLq(1)*gLq(2) - gRq(1)*gRq(2)
C     this simple routine returns the Helicity Structure Functions, 
C       SFnHel(-1:1)             (F_Left/Long./Right) 
C     and the Tensor Structure Functions,
C       SFnTen(1,2,3), defined as F_1/2/3:
C       F_1 = F_L + F_R              (2.F_transverse = 2.F_1^conventional)  
C       F_2 = F_L + F_R + 2.F_long   (F_2^conventional / x)
C       F_3 = F_L - F_R
C                  ------------------------------
      Implicit None
C                                    ! Dummy variables
      Real*8, Intent(in) :: SFnIred(5), g2chi(3)
      Real*8, Intent(out):: SfnTen(3), SFnHel(-1:1)
C                                    ! Local variables
      Real*8 TransPC, TransPV
C                  ------------------------------
C                                                 parity-conserving combinations
      TransPC   = g2chi(1) *SFnIred(1) + g2chi(2) *SFnIred(2)      
      SFnHel(0) = g2chi(1) *SFnIred(3) + g2chi(2) *SFnIred(4)
C                                                  parity-violating Piece   
      TransPV   = g2chi(3) *SFnIred(5)   

      SFnHel( 1) = TransPC - TransPV 
      SFnHel(-1) = TransPC + TransPV 
      
      SFnTen(1) = 2d0 * TransPC
      SFnTen(2) = SFnTen(1) + 2d0 * SFnHel(0)    ! A Power-law correction factor omitted here
      SFnTen(3) = 2d0 * TransPV                  ! A Power-law correction factor omitted here
      
      End Subroutine SfnConv                 ! ******************************
 

C                --------------------------------
      Subroutine F2LConvNC(F2o, FLo, SFnTen, SFnHel)
Cpn11 Construct tensor and helicity structure functions SFnTen and
c     SFnHel from input structure functions F2, FL
c
c     Currently works for neutral-current DIS only

      Real*8, Intent(in) :: F2o, FLo
      Real*8, Intent(out):: SfnTen(3), SFnHel(-1:1)

                    ! input variables Flo3,F2o,F3o3 inherited from parent.
      SFnHel(0)= FLo
      SFnTen(2)= F2o
      SFnTen(3)= 0d0 ! No F_3 structure function for the moment

      SFnTen(1) = SFnTen(2)- 2d0 * SFnHel(0)
      SFnHel(1) =(SFnTen(1)- SFnTen(3))/2d0
      SFnHel(-1)=(SFnTen(1)+ SFnTen(3))/2d0

      End Subroutine F2LConvNC

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE HEL2TEN_CT(XBJ,Q,HMASS,XHLR,XH123)
C-----------------------------------------------------------------------------
C      Computes F123 HADRON Helicity Amps from CAOT paper
C  2005-12-18:    wkt: This version uses the "modern" definitions of F_1,2,3
C           while the original ACOT papers uses the "conventional" definitions.
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension XH123(3), XHLR(-1:1)
      PARAMETER (PI=3.14159265359)


      XPLUS=XHLR( 1)
      XZERO=XHLR( 0)
      XMINU=XHLR(-1)

      RHO=1.0
      IF(HMASS.GT.0d0) THEN 
        XNU=Q**2/(2.0*HMASS*XBJ)
        RHO=SQRT(1.0+(Q/XNU)**2)  
      ENDIF

      XH123(1)=(XPLUS+XMINU)
      XH123(2)=(XPLUS+XMINU+2d0*XZERO)/RHO**2
      XH123(3)=(-XPLUS+XMINU)/RHO

      RETURN
      END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE TEN2HEL(XBJ,Q,HMASS,XHLR,XH123)
C-----------------------------------------------------------------------------
C      Computes F123 HADRON Helicity Amps from CAOT paper
C  2005-12-18:    wkt: This version uses the "modern" definitions of F_1,2,3
C           while the original ACOT papers uses the "conventional" definitions.
C---------------------------------------------------------------------

      Implicit Double Precision (A-H, O-Z)
      Dimension XH123(3), XHLR(-1:1)
      PARAMETER (PI=3.14159265359)


      F1=XH123(1) 
      F2=XH123(2) 
      F3=XH123(3) 

      RHO=1.0
      IF(HMASS.GT.0d0) THEN 
        XNU=Q**2/(2.0*HMASS*XBJ)
        RHO=SQRT(1.0+(Q/XNU)**2)  
      ENDIF

      XPLUS= (F1 -RHO*F3) / 2d0
      XMINU= (F1 +RHO*F3) / 2d0
      XZERO= (-F1 +RHO**2 *F2) / 2d0

      XHLR( 1)=XPLUS
      XHLR( 0)=XZERO
      XHLR(-1)=XMINU

      RETURN
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      Function TEN2SIG
     >  (XSIGN,ICHARGED,XBJ,Q,GLLEP,GRLEP,XLEPOL,HMASS,SBIG,F123)
C-----------------------------------------------------------------------------
C     Compute dSig/dx/dy 
C     Arguments modified: FIO 8/15/02
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       Dimension F123(3)
       PARAMETER (PI=3.14159265359)

         TEN2SIG= 0.0

C       ICHARGED     =  0 FOR NC;  NON-ZERO FOR CC
C       XSIGN= +1.0  !*** THIS IS FOR      LEPTON SCATTERING.
C       XSIGN= -1.0  !*** THIS IS FOR ANTI-LEPTON SCATTERING.
C  **** ANTI-LEPTON SCATTERING IS OBTAINED BY GLLEP <=> GRLEP _OR_ BY FLIPPING XSIGN

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
