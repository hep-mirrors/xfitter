      subroutine SF_SACOTCHI_wrap(idis,XSecType,x,q2,y,StrFn,iLptn1,iLptn2,xsec)
      implicit none
#include "steering.inc"
c call the input file which sets the parameters for the
c S-Acot-chi calculation
      integer Ldtadj,idis
      character*(*) XSecType
      logical UseKFactors
      double precision x, Q2, y, StrFn(0:3),StrFnGm
      double precision Sfn123(3),SfnHel(-1:1)
      double precision F123(1:3,1:5),Fhel(-1:1,1:5)
      double precision Q
      integer iHadn,iCmptn, iFlv, jSf
      double precision Prefactor
      double precision gL,gR,gv,ga 
      double precision AlfEwk, vBnMas, zeta, gEWlH             !, EffXsc
      external StrFnGm, AlfEwk, vBnMas, gEWlH

      Integer :: iOrdr, iSchm,iFschm, mSwtch, iFscle
      double precision  :: Q0, ScalFc, achi     
      Common / cDisCal / Q0, ScalFc, achi, iOrdr, iSchm

c the flag that choose which program to use for alphas and PDF evolution
      integer iAlsEvl, iPDFEvl 
      Common / Evolpak / iAlsEvl, iPDFEvl

c the flag for choosing pole mass of running mass as input
c currently, we concern only polemass input.
      Double precision POLEMASS, MSBMASS, qms(9)
      Integer iMass
      Logical MassesSet
      COMMON / COMQMS/ POLEMASS(9), MSBMASS(9), iMass, MassesSet

      double precision fswitch, al
      integer nf, norder
      logical set
      COMMON / QCDPAR_CT / FSWITCH, AL, NF, NORDER, Set

c
      double precision xSecLh,xsec 
      integer iLptn1, iLptn2, iXsc
      
c test part
c      integer i, j,k
c      integer Ir,iret
c      double precision Pdf_CT,tmp, AMASS,alpi,value
c      external Pdf_CT,amass,alpi
c      integer NGN
c      PARAMETER (NGN = 3)
c      double precision VKM
c      COMMON / KMATRX / VKM (NGN, NGN)

      !use polemass as input
      imass=0 

      ! interface to replcae the polemass from the ewparam.txt
      call polemassinterface

      !interface for reassign the number of flavor for future using
      !call NFinterface

      ! set switching points in alpha_s evolution
      ! preserve for future to use running mass as input.
      !call SetAmhatf            

      ! synchronize the pole and MSbar quark masses
      ! preserve for future to use running mass as input.
      !call MassSync             

      ! set the using of 
      iAlsEvl=1
      iPDFEvl=1

      Q0 = dsqrt(dble(starting_scale))
                ! Q0 is the Qmin below which PDFs are not defined;
                ! value obtained from steering.txt
      ScalFc=dble(scalea1)  !1d0, the value from fitinp
                ! With iFscle = 1, aMu is defined to be always > Q0.
                ! ScalFc is the (dimensionless) factorization scale;
                ! With iFscle = 2,3, aMu ~ ScalFc * the physical scale. 
                ! value obtained from steering.txt
      achi=0d0  ! value obtain from fitinp
                ! achi is the shape parameter in the rescaling variable chi
                ! 0.0 :                   Original rescaling variable chi
                ! 0.0:1.0:       Generalized rescaling variable zeta(lambda) from 
                !                arXiv:0903.2667 by Nadolsky and Tung

      !I_FIT_ORDER  !>  Fit order (LO=1, NLO=2, NNLO=3) 
      iOrdr=I_FIT_ORDER-1
      !iOrdr=0   ! the order of the DIS calculation in PQCD approximation
                ! =0, LO
                ! =1, NLO
                ! =2, NNLO
      iFSchm = 0
                ! flavor # factorization scheme
                ! =1~6, FFN
                ! =0, VFN
      mSwtch = 3
                ! mass treatment switch
                ! =0, Naive ZM
                ! =1, set a special factorization, scale for some realizations of improved ZM
                ! =2, GM,Improved ZM scheme
                ! =3, GM scheme
                ! =8, improved ZM
      iFscle = 2
                ! factorization scale switch,
                ! check explain of ScalFc.
                ! =0,
                ! =1,
                ! =2,
                ! =3 Redefine the fac scale for semi-inclusive HQ prod.
      iSchm = iFschm + mSwtch*10 + iFscle*100

      Fswitch=0 ! rescaling factor
                ! the value from fitinp

      Q=dsqrt(Q2)
c      jsf=0
c      if (jSf==0) then        !CalculateDIS returns F_L = F_2 - 2*x*F_1. 
c                              !StrFnGm returns F_{long} = F_2/(2*x) - F_1.
c                              !Prefactor = F_L/F_{long}
c        Prefactor=2*x
c      else                    !For i=1,2,3 StrFnGm returns F_i/x
c        Prefactor= x
c      endif                   !jsf

      iHadn   = 1  !proton
      if     ((XSecType=='NCDIS').or.(XSecType=='CCDIS')) then
        iFlv    = 0            !inclusive
      else if (XSecType=='CHARMDIS') then
        iFlv    = 14
      end if

      if (idis==1) then !using structure function as output
        StrFn   = 0d0
        F123    = 0d0
        Fhel    = 0d0
        
        do iCmptn = 1,5
          if ((iCmptn==2).or.(iCmptn==3)) cycle
          call SfnGm(iHadn,iCmptn,iFlv, X,Q, Sfn123,SfnHel)
          F123( 1:3,iCmptn) = Sfn123( 1:3)*x
          Fhel(-1:1,iCmptn) = SfnHel(-1:1)*x
          Fhel(   0,iCmptn) = Fhel( 0,iCmptn)*2 
        end do !iCmptn
        
        ! according to Eq(1) and (2) of 1506.06042v3
        zeta = AlfEwk(4) / AlfEwk(1)             ! ratio of gauge couplings
        zeta = zeta * Q2 / (Q2 + vBnMas(4)**2)   ! ratio of propagators
        
        gL=gEWlH ( 1, 4, -1, 1)  !gEWlH ( iLep1, 4, -1, iLep2)
        gR=gEWlH ( 1, 4,  1, 1)  !gEWlH ( iLep1, 4,  1, iLep2)
        gv=gR+gL
        ga=gR-gL
        ! F2
        StrFn(2) = F123(2,1) - gv*zeta * F123(2,5) 
     $           + (gv**2+ga**2)*zeta**2 * F123(2,4)
        ! F_long
        StrFn(0) = Fhel(0,1) - gv*zeta * Fhel(0,5) 
     $           + (gv**2+ga**2)*zeta**2 * Fhel(0,4)     
        ! F3
        StrFn(3) = - ga*zeta * F123(3,5) 
     $           + 2d0*ga*gv*zeta**2 * F123(3,4)
     
      else if (idis==2) then ! using reduce cross section as output

        iXsc=0
        xsec= xSecLh (iLptn1, iLptn2, iHadn, x, Q, y, iFlv, iXsc)

C     iXsc :   0   : DESY Reduced Cross sections
C              1   : My Effective xSec. defined in the DisPac05 tech note
C              2   : d Sig / dx dy
C              3   : d Sig / dx dQ^2
C              4   : .... add as you please!
      end if !idis

      return
      end subroutine

c test
      subroutine t2f(t)
      implicit none
      integer v(8)
      character date*8,time*10,zone*5,t*23

      call date_and_time(date,time,zone,v)
      !t=date(1:4)//"."//date(5:8)//"."//time(1:4)//"."//zone
      t=date(1:4)//"."//date(5:8)//"."//time(1:4)//"."//time(5:6)//"."//zone

      return
      end subroutine

!=====================================================================
! interfaces

      function alpi(m)
      !return alpha_s value
      double precision alpi,m,Q2,pi
      double precision hf_get_alphas
      external hf_get_alphas
      integer iAlsEvl, iPDFEvl,i 
      parameter (pi=3.14159265358979323846263383)
      Common / Evolpak / iAlsEvl, iPDFEvl
      
      if (iAlsEvl==1) then  !using alpha_s from QCDNUM
        Q2=m**2
        alpi=hf_get_alphas(Q2)/pi
      end if

      return
      end function

      subroutine NFinterface
!interface for reassign the number of flavor for future using
      implicit none
      integer iret
      double precision value
      
      value=5
      call PARQCD(1,"NFL",VALUE,IRET)

      return
      end subroutine


      subroutine polemassinterface
c heavy quark mass from ewparam.txt, the xfitter
      implicit none
#include "couplings.inc"
      Double precision POLEMASS, MSBMASS, qms(9)
      Integer iMass
      Logical MassesSet
      COMMON / COMQMS/ POLEMASS(9), MSBMASS(9), iMass, MassesSet

      polemass(1) = mup
      polemass(2) = mdn
      polemass(3) = mst
      polemass(4) = mch
      polemass(5) = mbt
      polemass(6) = mtp

      return
      end subroutine
      

      FUNCTION PDFP (LPRTN, XD, QD)
      double precision PdfP, XD, QD, PDF
      double precision x,Q
      integer LPRTN
      double precision q2,xpdf(-6:6)
      integer iAlsEvl, iPDFEvl 
      Common / Evolpak / iAlsEvl, iPDFEvl

      x=XD
      Q=QD

      if (iPDFEvl==1) then   !using PDF from QCDNUM
        q2=Q**2
        call hf_get_pdfs(x,q2,xpdf)
        PDF=xpdf(LPRTN)/x
        if (LPRTN== 2) PDF=xpdf( 1)/x  /1.00005d0 !agreement test
        if (LPRTN== 1) PDF=xpdf( 2)/x  /1.00005d0 !agreement test
        if (LPRTN==-1) PDF=xpdf(-2)/x
        if (LPRTN==-2) PDF=xpdf(-1)/x

        if (LPRTN== 0) PDF=xpdf( 0)/x  /1.00005d0 !agreement test
        if (LPRTN== 3) PDF=xpdf( 3)/x  /1.00005d0 !agreement test 
        if (LPRTN==-3) PDF=xpdf(-3)/x  /1.00005d0 !agreement test 

        if ((x>=0.999d0).or.(x< 1d-16)) PDF=0d0
      end if

      PDFP=PDF

      RETURN
      END


      Subroutine FILKMM
C                                                   -=-=- filkmm
C                       Given the Mixing parameters for NGN generations,
C                       calculate the KM matrix in the Cabbibo-KM-Maiani-
C                       -Wolfenstein-Chau-Keung...etc scheme.
C                       Used in SetEwk.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "couplings.inc"
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NGN = 3, NANG= NGN*(NGN-1)/2, NPHS=(NGN-1)*(NGN-2)/2)

      COMMON / KMATRX / VKM (NGN, NGN)
      COMMON / MIXPAR / CMX(NANG), DMX(NPHS)

C     (how to handle complex numbers due to phase factors efficiently,
C      bearing in mind that most applications do not need this info??)

C     For the moment, put the absolute values of the matrix elements in by hand

      DIMENSION AK (3, 3)
c CKM values from Particle data group
c W.-M. Yao et al., J. of Phys. G 33, 1 (2006)
c      DATA (AK(I,1), I=1,3) / 0.97383, 0.2272, 0.00396 /
c      DATA (AK(I,2), I=1,3) / 0.2271, 0.97296, 0.04221 /
c      DATA (AK(I,3), I=1,3) / 0.00814, 0.04161, 0.9991 /

c htj xfitter interface
c          Vud, Vus, Vub 
c  AK =    Vcd, Vcs, Vcb
c          Vtd, Vts, Vtb
      AK(1,1) = Vud
      AK(2,1) = Vus
      AK(3,1) = Vub 
      AK(1,2) = Vcd 
      AK(2,2) = Vcs 
      AK(3,2) = Vcb
      AK(1,3) = Vtd 
      AK(2,3) = Vts 
      AK(3,3) = Vtb

      DO 5 I = 1, 3
      DO 6 J = 1, 3
        VKM (I,J) = AK (I,J)
    6 CONTINUE
    5 CONTINUE

      Return
C                        ****************************
      END

C                                                          =-=-= EwCpl0



      Subroutine SetEwk(ID_WZ)
C                                                   -=-=- setewk
ccpy 03/14/2012: ID_WZ=1 for on-shell W and Z oson productions; otherwise 0.

C 05-06-01 : Initialize the EW couplings to Zero first in Data Block.
C            Non-zero values will be filled by the program.
C            Can't assume all compilers will initialize the unset entries
C            to zero.
C 05-05-10 : trivial change: initialize Lprt ->0, so that the coupling
C            don't get printed out each time.  Bug: Lprt can't be set.

C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: SetEwk
C Description: initial setup of the EwkPac.
C ListOfFiles: setewk filcpl filkmm
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2ewk/RCS/SetEwk.f,v 6.1 98/08/16 17:24:43 wkt Exp $
C #Log: SetEwk.f,v $
c Revision 6.1  98/08/16  17:24:43  wkt
c Re-organization; rationalization; initialization for DIS & DY
c 
C   In DisPac: QCPLN ---> EwCplSc
C   In VbpPac: SQCPAN --> EwCpl2An
C              SQCPCM --> EwCpl2Cn
c 
C Revision 6.0  96/12/25  14:38:28  wkt
C Synchronize with version 6 of all the other pac's

c Revision 1.1  94/02/22  11:10:59  lai
c Initial revision
c 
c Revision 2.2  92/03/08  15:21:24  wkt
c Setxxx.f uniformly added to force blockdata linking and 
c to perform other initiation Functions.
c 
C                 Setup the Ewkpac, which contains the following modules

C      Subroutine SetEwk
C      Subroutine FILCPL (CH, GV, GA, GL, GR)
C      Subroutine FILKMM 
C
C      Function ALFEwk (IBOSON)
C      Entry VBNMAS (IBOSON)
C      Entry SWG2F ()
C      Entry GEWLT (IT1, IBS, IBT, IT2)
C      Entry GEWLH (IT1, IBS, IBH, IT2)
C      Entry GEWQT (IQ1, IBS, IBT, IQ2)
C      Entry GEWQH (IQ1, IBS, IBH, IQ2)
C
C      Subroutine SetEwCpl2
C      Function SQCPAN (JP1, JBN, JP2)
C      Entry SQCPCM (JP, JBN)
C      Function EWCPL0 (IBSN, IP1, IP2, IRT)   (Inactive)

C     Set up the basic Electro-weak coupling matrices for the Boson-Fermion
C     Yukawa coupling term in the Effective Lagrangian

C     Boson label: (IBN)   1,   2,   3,   4
C                       gamma   W+   W-   Z    
  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "couplings.inc"
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NSP = 2, NGN = 3, NBN = 4, NPOL = 3)
      PARAMETER (NFL = NSP * NGN)
 
      DIMENSION IP(NSP, NGN), CHL(NSP), CHQ(NSP)
      DIMENSION GGQA(NSP,NGN,NSP,NGN), GGQV(NSP,NGN,NSP,NGN),
     >          GGQR(NSP,NGN,NSP,NGN), GGQL(NSP,NGN,NSP,NGN)
 
      External DatEwk

      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / EwkPAR / WMS, ZMS, SWG2, ALFE, ALFEW(NBN)
      COMMON / EwkPAR2 / WCOUPL, ZCOUPL
      COMMON / KMATRX / VKM (NGN, NGN)

C            Left, Right, Vector, and Axial-vector couplings of Leptons & Quarks

      COMMON / EW1LCP / GLL(NSP,NBN,NSP), GLR(NSP,NBN,NSP), 
     >                  GLV(NSP,NBN,NSP), GLA(NSP,NBN,NSP)
      COMMON / EW1QCP / GQL(NSP,NBN,NSP), GQR(NSP,NBN,NSP), 
     >                  GQV(NSP,NBN,NSP), GQA(NSP,NBN,NSP)
      COMMON / EW2QCP / HQL(NFL,NBN,NFL), HQR(NFL,NBN,NFL), 
     >                  HQV(NFL,NBN,NFL), HQA(NFL,NBN,NFL)
 
      DATA IP, Lprt   / 1, 2, 4, 3, 6, 5, 0 /
      DATA CHL,CHQ / 0d0, -1d0, 0.666666667, -0.333333333 /

c external function from xfitter
      double precision AEMRUN
      external AEMRUN
 
C     Adopt the convention:

C                |  IBSN      coupling       current           coupling**2/4pi
C                |
C   ------------ |------------------------------------------------------------
C      Photon    |    1          e             Gmu            ALFE = Alpha(QED) 
C                |
C      W+/W-     |   2/3      g/2Sqrt2      Gmu(1-G5)          ALFE / SWG2 / 8
C                |
C        Z       |    4        g/2CosW  Gmu((T3-2QSWG2)-T3G5)   ALFE / S2WG2
C                |
C   ------------ |-------------------------------------------------------------
C     where Gmu = Gamma-super-mu; G5 = Gamma-5; Cos W = Cos(Theta-Weinberg);
C     SWG2 = (Sin W)**2; & S2WG2 = (Sin 2W)**2.
C   ---------------------------------------------------------------------------
C                                       Fill the overall coupling**2 array ALFEW
ccpy 03/14/2012
      If (ID_WZ.ne.1) then
CLai 11/14/06 replace the following by CP's note
        ALFEW (1) = ALFE
        ALFEW (2) = ALFE / SWG2 / 8.
        ALFEW (3) = ALFEW (2)
        ALFEW (4) = ALFE / SWG2 / (1.-SWG2) / 4.

      Else
C-----------------------------------------
C the following adopted from CP's note on EW couplings
C for on-shell W or Z productions
C----------------------------------------- 
C Input parameters:

c htj xfitter interface:

C_____W-BOSON MASS
      WMS = mw  !80.39d0
C_____Z BOSON MASS      
      ZMS = mz  !91.187d0
C_____FERMI CONSTANT      
      GMU = gf  !1.16637D-5
       
C_____EFFECTIVE SIN^2(THETA_WEAK) AT Z-POLE
      SW2_EFF_MZ = sin2thw    !0.23143
      SWG2 = SW2_EFF_MZ

c htj 2016.0201
c by now, the running alpha_em function of xfitter has problem. 
c so we set the value of ALFA_EM_MZ value by hardwired value before it fixed. 

C_____RUNNING QED COUPLING AT Z-POLE
      call eprc_init(.true.)
      ALFA_EM_MZ = AEMRUN(MZ**2)
c      ALFA_EM_MZ = 1.0d0/128.937d0


C_____RUNNING QED COUPLING AT ELECTRON MASS SCALE
      ALFA_EM_ME = alphaem !1.0D0/137.0359895D0

C_____PI
      !PI = 3.1415926535898d0 
        
C------------------------------------ 
C Derived variables:

      R2 = SQRT(2.D0)
      XMW2 = WMS**2
      XMZ2 = ZMS**2
      
C_____Weak couplings: g_w
      GWEAK2 = 4.D0*R2*XMW2*GMU
      GWEAK = SQRT(GWEAK2)
      
      
C PHOTON COUPL IS (e**2)/2      
C_____RUNNING QED COUPLING AT THE SCALE AMU (CALCULATED IN 
C_____FUNCTION ALFAEM 
cLai  comment ALFE out since AMU depends on expt. 
c     Use the fine structure constant here from DatEwk
c      ALFE = ALFAEM(AMU) 
CCPY 02/22/12 : Use ALFE at MZ scale
      ALFE = ALFA_EM_MZ
      
C WCOUPL IS (g**2)/8
      WCOUPL = GWEAK2/8.0D0
      
C WE IGNORE THE DIFFERENCE BETWEEN b-quark AND OTHER FERMIONS
C ZCOUPL IS (g/Cos_w)**2/4

C_____Sin^2 and Cos^2 of the weak angle - in the on shell scheme
      CWS = XMW2/XMZ2
      SWS = 1.D0-CWS
      ZCOUPL_MZ = PI*ALFA_EM_MZ/SWS/CWS
      ZCOUPL = ZCOUPL_MZ
      
C------------------------------------
C To a good approximation:

      ALFEW (1) = ALFE
      ALFEW (2) = WCOUPL/4.0/PI
      ALFEW (3) = ALFEW (2)
      ALFEW (4) = ZCOUPL/4.0/PI
  
C------------------------------------

      EndIf        ! ID_WZ
C                                                   ----------------------------
C                                                       For Each Generation:
C                                     Fill the relative coupling in the currents
C                           GXL/R = (GXV -/+ GXA)/2.      (X = Lepton / Quark) ;
C                             hence the currents are:   GXL/R * Gmu * (1-/+G5)

      CALL FILCPL (CHL, GLV, GLA, GLL, GLR)
      CALL FILCPL (CHQ, GQV, GQA, GQL, GQR)
C                                                  ----------------------------
C                          For Quarks, put in Generations and generation-mixing

C                                 fill in the KM matrix array from data in the
C                           MIXing PARameter common block / MIXPAR / consisting
C                                 of Cosines of mixing angles and Phase angles.
      CALL FILKMM

C                     For W+ / W- , Multiply the KM matrix into the Down quarks
      DO 10 IG1 = 1, NGN
      DO 10 IG2 = 1, NGN
      DO 10 IS1 = 1, NSP
      DO 10 IS2 = 1, NSP
        GGQA (IS1, IG1, IS2, IG2) = GQA(IS1, 2, IS2) * VKM(IG1, IG2)
        GGQV (IS1, IG1, IS2, IG2) = GQV(IS1, 2, IS2) * VKM(IG1, IG2)
        GGQR (IS1, IG1, IS2, IG2) = GQR(IS1, 2, IS2) * VKM(IG1, IG2)
        GGQL (IS1, IG1, IS2, IG2) = GQL(IS1, 2, IS2) * VKM(IG1, IG2)
   10 CONTINUE
C                                                  ----------------------------
C                         (Generation, ISPin) labels  ----> parton flavor label
      DO 20 IG1 = 1, NGN
      DO 20 IG2 = 1, NGN
      DO 20 IS1 = 1, NSP
      DO 20 IS2 = 1, NSP
        IP1 = IP(IS1, IG1)
        IP2 = IP(IS2, IG2)

        DO 21 IBN = 1, NBN
        IF (IBN .EQ. 1 .OR. IBN .EQ. 4) THEN
         IF (IG1 .NE. IG2) THEN
            HQA (IP1, IBN, IP2) = 0.0
            HQV (IP1, IBN, IP2) = 0.0
            HQL (IP1, IBN, IP2) = 0.0
            HQR (IP1, IBN, IP2) = 0.0
         ELSE
            HQA (IP1, IBN, IP2) = GQA(IS1, IBN, IS2)
            HQV (IP1, IBN, IP2) = GQV(IS1, IBN, IS2)
            HQL (IP1, IBN, IP2) = GQL(IS1, IBN, IS2)
            HQR (IP1, IBN, IP2) = GQR(IS1, IBN, IS2)
         ENDIF
        ENDIF
   21   CONTINUE

            HQA (IP1, 2, IP2) = GGQA(IS1, IG1, IS2, IG2)
            HQV (IP1, 2, IP2) = GGQV(IS1, IG1, IS2, IG2)
            HQL (IP1, 2, IP2) = GGQL(IS1, IG1, IS2, IG2)
            HQR (IP1, 2, IP2) = GGQR(IS1, IG1, IS2, IG2)

            HQA (IP1, 3, IP2) = GGQA(IS2, IG2, IS1, IG1)
            HQV (IP1, 3, IP2) = GGQV(IS2, IG2, IS1, IG1)
            HQL (IP1, 3, IP2) = GGQL(IS2, IG2, IS1, IG1)
            HQR (IP1, 3, IP2) = GGQR(IS2, IG2, IS1, IG1)

   20 CONTINUE
C                                                  ----------------------------
C                                                                      Finished
C                          Set Lprt>0 to print Quark couplings for confirmation
      If (Lprt .ge. 1) Then
        Print '(A)', 
     >  'Flavor-dependent part of Quark-VectorBoson Vertex:'
        Do 41 Ibn = 1, 4
      Print '(/A, I4)', 'Iboson = ', Ibn
      Print '(/A/(4F12.5))', 'Axial-v: Iflv=1-4:'
     >, ((HqA(Ip1,Ibn,Ip2), Ip1=1,4), Ip2=1,4)
      Print '(/A/(4F12.5))', 'Vector: Iflv=1-4:'
     >, ((HqV(Ip1,Ibn,Ip2), Ip1=1,4), Ip2=1,4)
      Print '(/A/(4F12.5))', 'Left-hand: Iflv=1-4:'
     >, ((HqL(Ip1,Ibn,Ip2), Ip1=1,4), Ip2=1,4)
      Print '(/A/(4F12.5))', 'Right-hand: Iflv=1-4:'
     >, ((HqR(Ip1,Ibn,Ip2), Ip1=1,4), Ip2=1,4)
 41     Continue
      Endif

      Return
C                        ****************************
      End
