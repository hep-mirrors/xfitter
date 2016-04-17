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
**
*     Input Variables
*
      integer PtOrder
      double precision MCharm,MBottom,MTop
      double precision Q_ref,Alphas_ref
      character*7 Scheme
      character*5 MassScheme
      logical runm
      double precision Q2save
      common / PreviousQ2 / Q2save
**
*     Internal Variables
*
      double precision Qc,Qb,Qt,McQ,MbQ,MtQ
      double precision HeavyQuarkMass
*
*     Initialize Q2save
*
      Q2save = 0d0
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
      call EnableDynamicalScaleVariations(.true.)
      if(MassScheme(1:4).eq."Pole")then
         call SetPoleMasses(MCharm,MBottom,MTop)
      elseif(MassScheme.eq."MSbar")then
         call SetMSbarMasses(MCharm,MBottom,MTop)
         call EnableMassRunning(runm)
      endif
      call SetAlphaQCDRef(Alphas_ref,Q_ref)
      call SetPerturbativeOrder(PtOrder)
      call SetMassScheme(scheme)
      call SetNumberOfGrids(3)
c      call EnableDampingFONLL(.false.)
      call SetGridParameters(1,50,3,9.8d-7)
      call SetGridParameters(2,40,3,1d-2)
      call SetGridParameters(3,40,3,7d-1)
*
      Qc = MCharm
      Qb = MBottom
      Qt = MTop
      if(Qc.ne.MCharm.or.Qb.ne.MBottom.or.Qt.ne.MTop)then
         call InitializeAPFEL()
         McQ = Qc
         if(Qc.ne.MCharm)  McQ = HeavyQuarkMass(4,Qc)
         MbQ = Qb
         if(Qb.ne.MBottom) MbQ = HeavyQuarkMass(5,Qb)
         MtQ = Qt
         if(Qt.ne.MTop)    MtQ = HeavyQuarkMass(6,Qt)
         if(MassScheme(1:4).eq."Pole")then
            call SetPoleMasses(McQ,MbQ,MtQ)
         elseif(MassScheme.eq."MSbar")then
            call SetMSbarMasses(McQ,MbQ,MtQ)
            call SetMassScaleReference(Qc,Qb,Qt)
         endif
      endif
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
      subroutine sf_fonll_wrap(x,Q2,muoQ,F2,FL,xF3,F2c,FLc,F2b,FLb)
*
      implicit none
**
*     Input Variables
*
      double precision x,Q2,muoQ
**
*     Internal Variables
*
      double precision Q
      double precision F2total,FLtotal,F3total
      double precision F2charm,FLcharm
      double precision F2bottom,FLbottom
      double precision Q2save
      common / PreviousQ2 / Q2save
**
*     Output Variables
*
      double precision F2,FL,xF3,F2c,FLc,F2b,FLb
*
      if(Q2.ne.Q2save)then
         Q = dsqrt(Q2)
         call ComputeStructureFunctionsAPFEL(Q*muoQ,Q)
      endif
*
      F2  = F2total(x)
      FL  = FLtotal(x)
      xF3 = F3total(x)
      F2c = F2charm(x)
      FLc = FLcharm(x)
      F2b = F2bottom(x)
      FLb = FLbottom(x)
*
      Q2save = Q2
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
