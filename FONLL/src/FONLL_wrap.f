************************************************************************
*
*     Initialization routine for the FONLL scheme.
*
************************************************************************
      subroutine FONLL_Set_Input(MassScheme,runm,Mcharm,MBottom,MTop,
     1                           Q_ref,Alphas_ref,PtOrder,Scheme)
*
      implicit none
*
      include 'couplings.inc'
*
      integer PtOrder
      double precision MCharm,MBottom,MTop
      double precision Q_ref,Alphas_ref
      character*7 Scheme
      character*5 MassScheme
      logical runm
*
*     Set EW parameters
*
      call SetZMass(Mz)
      call SetWMass(Mw)
      call SetSin2ThetaW(sin2thw)
      call SetGFermi(gf)
      call SetCKM(Vud, Vus, Vub,
     1            Vcd, Vcs, Vcb,
     2            Vtd, Vts, Vtb)
*
*     APFEL settings
*
c      call EnableTargetMassCorrections(.true.)
      call EnableDynamicalScaleVariations(.true.)
      if(MassScheme(1:4).eq."Pole")then
         call SetPoleMasses(Mcharm,MBottom,MTop)
      elseif(MassScheme.eq."MSbar")then
         call SetMSbarMasses(Mcharm,MBottom,MTop)
         call EnableMassRunning(runm)
      endif
      call SetAlphaQCDRef(Alphas_ref,Q_ref)
      call SetPerturbativeOrder(PtOrder)
      call SetMassScheme(scheme)
      call SetNumberOfGrids(3)
      call SetGridParameters(1,40,3,9.8d-7)
      call SetGridParameters(2,30,3,1d-2)
      call SetGridParameters(3,20,3,7d-1)
*
*     Initialize the APFEL DIS module
*
      call InitializeAPFEL_DIS()
*
      return
      end
*
************************************************************************
*
*     Routine that returns the structure functions
*
************************************************************************
      subroutine sf_fonll_wrap(x,Q2,F2,FL,xF3,F2c,FLc,F2b,FLb)
*
      implicit none
**
*     Input Variables
*
      double precision x,Q2
**
*     Internal Variables
*
      double precision Q
      double precision F2total,FLtotal,F3total
      double precision F2charm,FLcharm
      double precision F2bottom,FLbottom
**
*     Output Variables
*
      double precision F2,FL,xF3,F2c,FLc,F2b,FLb
*
      Q = dsqrt(Q2)
*
      call SetPDFSet("external1")
      call ComputeStructureFunctionsAPFEL(Q,Q)
*
      F2  = F2total(x)
      FL  = FLtotal(x)
      xF3 = F3total(x)
      F2c = F2charm(x)
      FLc = FLcharm(x)
      F2b = F2bottom(x)
      FLb = FLbottom(x)
*
      return
      end
*
************************************************************************
*
*     PDFs to be fet to the DIS module of APFEL
*
************************************************************************
      subroutine ExternalSetAPFEL1(x,Q,xf)
*
      implicit none
**
*     Input Variables
*
      double precision x,Q
**
*     Internal Variables
*
      integer ipdf
      double precision Q2
      double precision pdfsf(-6:6)
**
*     Output Variables
*
      double precision xf(-6:7)
*
      Q2 = Q * Q
      call HF_Get_PDFs(x,Q2,pdfsf)
*
      xf(7) = 0d0
      do ipdf=-6,6
         xf(ipdf) = pdfsf(ipdf)
      enddo
*
      return
      end
