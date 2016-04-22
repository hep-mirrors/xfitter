C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      SUBROUTINE GotSK1Q(Ihadx, Ipartx,F2Mx, Xx,Qx, aMux, SfnG)
    
    !  NLO contribution to the general SfnG due to initial state quark parton iPrtn
    !  without the EW coupling--the result is flavor-independent;
    !  Also not included is the overall \alpha_s / 2pi factor.
    !                             2006.02.10
    !  Formulas based on Gottschalk, as corrected by Stefan Kretzer etal
    !  Programmed by A. Belyaev, based on the SfZm05 structure (wkt)
    !           ---------------------------------------------- 
      Use mDisPac
      Implicit None
C >>                                       ! dummy variable
      integer, intent(in) :: IHadx, iPartx
      Real*8, intent (in) :: F2Mx,Xx, Qx, aMux
      Real*8, intent(out) :: SfnG(5)
C                                          ! Global variables:
      INTEGER idxi,idxj 
      Common/INDEX/idxi,idxj

      Integer :: iHad, iPart, Ired
      Real*8  :: x, Q, F2M, aMu
      Common / CxNloSK / x, Q, F2M, aMu, iHad,IPart, Ired
      
      REAL*8 DL
      Common/lambda/DL
C --                                     ! Local variables    
      Integer IR
      Real*8 :: ER,  F1,f2,f3, OmgLo(5)
     >        , Qmu1, Qmu2, Smll = 1d-6
C --                                     ! External functions
      Real*8 :: HQ0_N,HQ0_P1,HQ0_P,HQB_P, HQB_P1, 
     >          aPqq, aPqqFxSK, ADZINT_CT, PDF_CT,DKA,
     >        HQA_D,HQ0_D,HQ0_PM,HQB_PM

      EXTERNAL HQ0_N,HQ0_P1,HQ0_P,HQB_P, HQB_P1, 
     >          aPqq, aPqqFxSK, ADZINT_CT, PDF_CT,DKA,
     >        HQA_D,HQ0_D,HQ0_PM,HQB_PM

C                                !    ------- fill integration common
   
      REAL*8  FHQA_D,FHQB_P,FHQI(3),Pdfnx,FHQ0_N,FHQ0_D,FHQ0_P
       
      iHad	= iHadx
      iPart	= iPartx
      x    	= xx
      aMu  	= aMux
      F2M	= F2mx
      Q		= Qx
      
      DL=Q**2/(Q**2+F2M**2) !lambda
      Pdfnx=PDF_CT(iHad, iPart, x, aMu, IR)
      
      FHQ0_N= ADZINT_CT(HQ0_N, X, 1d0, aErr, rErr, ER, IR, 1, 1)			     
      
!  HQ0_D  is the delta term of h^q  
      FHQ0_D= HQ0_D(DL)*Pdfnx			     
      
!   HQ0_PM  is the multiplier of the Plus function in h^q  
!   HQ0_P   is the Plus function in h^q  
      FHQ0_P=ADZINT_CT(HQ0_P1, X, 1d0, aErr, rErr, ER, IR, 1, 1)
     >    -HQ0_PM(1d0)*Pdfnx
     >       *ADZINT_CT(HQ0_P,  0d0, X, aErr, rErr, ER, IR, 1, 1) 
c-----------------      
      FHQA_D= HQA_D(DL)*Pdfnx		     
      DO idxi=1,3
         FHQI(idxi)=FHQ0_N+FHQ0_D+FHQ0_P+FHQA_D

	 DO idxj=1,3
      	
           FHQB_P=ADZINT_CT(HQB_P1, X, 1d0, aErr, rErr, ER, IR, 1, 1)
     >       -HQB_PM(1d0)*Pdfnx
     >           *ADZINT_CT(HQB_P,  0d0, X, aErr, rErr, ER, IR, 1, 1)
	 
	    FHQI(idxi)=FHQI(idxi)+FHQB_P
	  ENDDO
      ENDDO
c-----------------      
      F1 =  FHQI(1)
      F2 =  FHQI(2)
      F3 =  FHQI(3)
      
      SfnG(1) = F1      / 2d0          ! SfnG_1 = F1/2
      SfnG(3) = (F2-F1) / 2d0          ! SfnG_3 = Flong = (F2-F1)/2
      SfnG(2) = 0d0        ! Zero-mass case has no chirality-mixing terms
      SfnG(4) = 0d0        ! Zero-mass case has no chirality-mixing terms
      SfnG(5) = F3 /2d0              ! SfnG_5 = F3/2 
				     
				     
      If (Abs(Q/aMu - 1d0) > Smll) Then
      Qmu1 =  ADZINT_CT(aPqq,  0d0, x, aErr, rErr, ER, IR, 1, 0)
     >     * Pdfnx
      Qmu2 = ADZINT_CT(aPqqFxSK, X, 1d0, aErr, rErr, ER, IR, 1, 1)
      
      OmgLo = (/5d-1, 0d0, 0d0, 0d0, 5d-1/)
      SfnG = SfnG + 2d0 * OmgLo * Log(Q/aMu) * (Qmu2 - Qmu1) 
      End If
              
      SfnG = SfnG * CF               ! QCD color factor

      End Subroutine GotSK1Q            !  *****************************
c================================================

      FUNCTION DKA(DL)
      IMPLICIT NONE
      real*8   DKA,DL
      DKA=0d0
      IF(DL.lt.1d0-1e-09)DKA=1d0/DL*(1d0-DL)*DLOG(1d0-DL)
      END FUNCTION DKA
c----------------------

      FUNCTION HQ0_D(DL)
      IMPLICIT NONE
      real*8   HQ0_D,dl,PI,DKA
      PI=dacos(-1d0)
      HQ0_D=
     > -(4d0+1d0/2d0/DL+PI**2/3d0+(1d0+3d0*DL)/2d0/DL*DKA(DL))
      END FUNCTION HQ0_D
      
c----------------------
      FUNCTION HQ0_N(Z)
      IMPLICIT NONE
c-       				 ! Global variables:
      Integer :: iHad, iPart, Ired
      Real*8  :: x, Q,  F2M, aMu
      Common / CxNloSK / x, Q, F2M, aMu, iHad,IPart, Ired
c-
      Integer IR
      Real*8 :: z, Pdfnxz,HQ0_N,PDF_CT
      
      Pdfnxz = PDF_CT(iHad, iPart, x/z, aMu, IR)
      HQ0_N=-(1d0+Z**2)*DLOG(z)/(1d0-Z)
      HQ0_N=HQ0_N*Pdfnxz/Z
      END FUNCTION HQ0_N
c----------------------
      FUNCTION HQ0_P(Z)
      IMPLICIT NONE
c-
      REAL*8 DL
      Common/lambda/DL
      Real*8  HQ0_P,Z
c-
      HQ0_P=(2d0*log(1d0-z)-log(1d0-dl*z))/(1d0-z)
      END FUNCTION HQ0_P
c----------------------
      FUNCTION HQ0_PM(Z)
      IMPLICIT NONE
      Real*8 :: HQ0_PM,Z
c-
      HQ0_PM=1d0+Z**2
      END FUNCTION HQ0_PM
c----------------------
      FUNCTION HQ0_P1(Z)
      IMPLICIT NONE

      Integer  iHad, iPart, Ired
      Real*8   x, Q, F2M, aMu
      Common / CxNloSK / x, Q, F2M, aMu, iHad,IPart, Ired
c-
      REAL*8 DL
      Common/lambda/DL
      
      Integer IR
      Real*8  z, Pdfnx, Pdfnxz,HQ0_P1,HQ0_P,PDF_CT,HQ0_PM
c-
      Pdfnxz = PDF_CT(iHad, iPart, x/z, aMu, IR)
      Pdfnx  = PDF_CT(iHad, iPart,   x, aMu, IR)
      
      HQ0_P1=HQ0_P(Z)*(HQ0_PM(z)*Pdfnxz/z - HQ0_PM(1d0)*Pdfnx)
      END FUNCTION HQ0_P1
c----------------------
      FUNCTION HQA_D(Z)
      IMPLICIT NONE
c-
      REAL*8 DL
      Common/lambda/DL
      
      INTEGER idxi,idxj 
      Common/INDEX/idxi,idxj
      
      Real*8 :: HQA_D,Z,DKA
       
      HQA_D=0d0
      IF   (idxi==2) HQA_D=DKA(DL)
      END FUNCTION HQA_D
c----------------------
      FUNCTION HQB_P(Z)
      IMPLICIT NONE
c-
      REAL*8 DL
      Common/lambda/DL
      
      INTEGER idxi,idxj 
      Common/INDEX/idxi,idxj
      
      Real*8 :: HQB_P,Z
      
      HQB_P=0
      IF   (idxj==1) THEN
        HQB_P=1d0/(1d0-z)
      ELSEIF(idxj==2) THEN
        HQB_P=1d0/(1d0-dl*z)
      ELSEIF(idxj==3) THEN 
        HQB_P=(1d0-z)/(1d0-dl*z)**2
      ENDIF
      END FUNCTION HQB_P
c----------------------
      FUNCTION HQB_PM(Z)
      IMPLICIT NONE
c-
      REAL*8 DL
      Common/lambda/DL
      
      INTEGER idxi,idxj 
      Common/INDEX/idxi,idxj
      
      Real*8 :: HQB_PM,Z
      HQB_PM=0
      IF   (idxj==1) THEN
	IF(idxi==1) HQB_PM=1d0-4d0*z+z**2
	IF(idxi==2) HQB_PM=2d0-2d0*z**2-2d0/Z
	IF(idxi==3) HQB_PM=-1d0-Z**2
      ELSEIF(idxj==2) THEN
	IF(idxi==1) HQB_PM=Z-Z**2
	IF(idxi==2) HQB_PM=2d0/Z-1d0-Z
	IF(idxi==3) HQB_PM=1d0-Z
      ELSEIF(idxj==3) THEN 
	HQB_PM=1d0/2d0
      ENDIF      
      END FUNCTION HQB_PM     
c----------------------
      FUNCTION HQB_P1(Z)
      IMPLICIT NONE

      Integer :: iHad, iPart, Ired
      Real*8  :: x, Q, F2M, aMu
      Common / CxNloSK / x, Q, F2M, aMu, iHad,IPart, Ired
c-
      REAL*8 DL
      Common/lambda/DL
      
      Integer IR
      Real*8 :: z, Pdfnx, Pdfnxz,HQB_P1,HQB_P,HQB_PM,PDF_CT
      
c-
      Pdfnxz = PDF_CT(iHad, iPart, x/z, aMu, IR)
      Pdfnx  = PDF_CT(iHad, iPart,   x, aMu, IR)

      HQB_P1=HQB_P(Z)*(HQB_PM(z)*Pdfnxz/z - HQB_PM(1d0)*Pdfnx)
      END FUNCTION HQB_P1
C=============================================
      Function aPqqFxSK (z)
      
      Implicit None
C                                          ! Global variables:
      Integer :: iHad, iPart, Ired
      Real*8  :: x, Q, F2M, aMu
      Common / CxNloSK / x, Q, F2M, aMu, iHad,IPart, Ired
C                                          ! Local variable
      Integer IR
      Real*8 :: z, Pdfnx, Pdfnxz, aPqqFxSK
C                                          ! External functions
      Real*8 :: aPqq, PDF_CT
    
      Pdfnx =  PDF_CT(iHad, iPart, x,   aMu, IR)
      Pdfnxz = PDF_CT(iHad, iPart, x/z, aMu, IR)

      aPqqFxSK =   aPqq(z) * (Pdfnxz/z - Pdfnx)

      End Function aPqqFxSK
