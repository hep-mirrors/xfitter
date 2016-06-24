      subroutine GetHVQMNRXsection(IDataSet)
C-----------------------------------------------------------------------
C Calculate NLO for heavy quarks in pp (MNR code + fragmentation)
C-----------------------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "fcn.inc"
#include "couplings.inc"

      ! Input dataset index
      integer IDataSet
      ! Data info (to be read from DATASETInfoNames)
      character*(64) Measurement, FinalState, XSecType
      ! Differential binning and normalisation flags
      logical d_pt,d_y,norm_y
      ! Centre-of-mass energy
      double precision sqrt_s
      ! Fragmentation fraction (to be read from DATASETInfo)
      double precision fragfraction
      ! Number of filled bins
      integer nbin
      ! Arrays to fill
      double precision xsec(NTOT)
      double precision ptmin(NTOT),ptmax(NTOT)
      double precision ymin(NTOT),ymax(NTOT)
      ! Bin indices
      integer idx,idxymin,idxymax,idxptmin,idxptmax,idxyminref,idxymaxref
      ! Last iteration counter
      integer IfcnCountLast
      save IfcnCountLast
      ! Heavy-quark masses
      double precision mc,mb
      ! Scale parameters
      double precision mf_A_c,mf_B_c,mf_C_c,mr_A_c,mr_B_c,mr_C_c
      double precision mf_A_b,mf_B_b,mf_C_b,mr_A_b,mr_B_b,mr_C_b
      ! Fragmentation parameters
      double precision fragpar_c,fragpar_b
      ! Called routines
      integer GetBinIndex, GetInfoIndex
      double precision FindXSPtY
      ! Loop counter
      integer i
      ! Debugging flag
      logical debug
      ! HVQMNR parameters (to be passed to C++ part)
      common/HVQMNR_PARS/sqrt_s,mc,mb,
     $                   mf_A_c,mf_B_c,mf_C_c,mr_A_c,mr_B_c,mr_C_c,
     $                   mf_A_b,mf_B_b,mf_C_b,mr_A_b,mr_B_b,mr_C_b,
     $                   fragpar_c,fragpar_b,debug
      save /HVQMNR_PARS/

      ! First call: make some checks, set needed parameters
      if(IfcnCount.eq.1) then
        ! Check HF scheme
        if(HFSCHEME.ne.3.and.HFSCHEME.ne.4) then
           print *,'Error in GetHVQMNRXsection: calculation does not ',
     $     'support HFSCHEME=',HFSCHEME
           call HF_stop
        endif
        ! Check number of bins
        if (NDATAPOINTS(IDataSet).gt.NTOT) then
           print *,'ERROR IN GetHVQMNRXsection'
           print *,'INCREASE NTOT_C to ',NDATAPOINTS(IDataSet)
           call HF_stop
        endif
        ! Set debugging flag
        debug = lDEBUG
      endif
      
      ! Read information from DATASETInfoNames
      if(DATASETInfoDimension(IDataSet).lt.3) then
        print *,'ERROR IN GetHVQMNRXsection: DATASETInfoNames '
        print *,'should start from Measurement, XSecType, FinalState'
        call HF_stop      
      endif
      Measurement = DATASETInfoNames(1,IDataSet)
      XSecType    = DATASETInfoNames(2,IDataSet)
      FinalState  = DATASETInfoNames(3,IDataSet)
      ! Read CMS energy and fragmentation fraction from DATASETInfo
      i = GetInfoIndex(IDataSet,'sqrt(S)')
      if(i.eq.0) then
        print *,'ERROR IN GetHVQMNRXsection: no sqrt(S)'
        call HF_stop      
      endif
      sqrt_s = DATASETInfo(i,IDataSet)
      i = GetInfoIndex(IDataSet,'fragfraction')
      if(i.eq.0.and.XSecType.eq.'abs') then
        print *,'ERROR IN GetHVQMNRXsection: no fragfraction for abs x-sections'
        call HF_stop      
      endif
      if(i.ne.0) then
         fragfraction = DATASETInfo(i,IDataSet)
      else 
         fragfraction = 0d0
      endif

      ! Read needed MINUIT parameters
      ! (enough to do this only once per each iteration)
      if(IfcnCountLast.ne.IfcnCount) then
        ! Update fcn counter
        IfcnCountLast = IfcnCount
        ! Heavy-quark masses
        mc = mch
        mb = mbt
        ! Scale parameters
        call GetMuPar('f','c',mf_A_c,mf_B_c,mf_C_c)
        call GetMuPar('r','c',mr_A_c,mr_B_c,mr_C_c)
        call GetMuPar('f','b',mf_A_b,mf_B_b,mf_C_b)
        call GetMuPar('r','b',mr_A_b,mr_B_b,mr_C_b)
        ! Fragmentation parameters
        call GetFPar('c',fragpar_c)
        call GetFPar('b',fragpar_b)
      endif
      
      ! Call calculation routine
      ! (calculation for any new measurement should be added here)
      if(Measurement.eq.'LHCb_7tev_charm') then
        call hvqmnr_lhcb_7tev_charm(IfcnCount,XSecType,FinalState, 
     $                              nbin,xsec,ymin,ymax,ptmin,ptmax)
        d_pt = .FALSE.
        d_y  = .FALSE.
      else if(Measurement.eq.'LHCb_7tev_beauty') then
        call hvqmnr_lhcb_7tev_beauty(IfcnCount,FinalState, 
     $                               nbin,xsec,ymin,ymax,ptmin,ptmax)
        d_pt = .TRUE.
        d_y  = .TRUE.
      else
        print *,'ERROR IN GetHVQMNRXsection: '
        print *,'unknown Measurement ',Measurement
        call HF_stop      
      endif
      
      ! Prepare bins
      idxymin  = GetBinIndex(IDataSet,'ymin')
      idxymax  = GetBinIndex(IDataSet,'ymax')
      idxptmin = GetBinIndex(IDataSet,'pTmin')
      idxptmax = GetBinIndex(IDataSet,'pTmax')
      if(idxymin.eq.0.or.idxymax.eq.0.or.idxptmin.eq.0.or.idxptmax.eq.0) then
        print *,'ERROR IN GetHVQMNRXsection:'
        print *,'no ymin or ymax or pTmin or pTmax bin'
        call HF_stop
      endif
      
      ! Determine normalisation in y and reference bins if needed
      if(XSecType.eq.'norm_y') then
        norm_y = .TRUE.
      else if(XSecType.eq.'abs') then
        norm_y = .FALSE.
      else
        norm_y = .FALSE.
        print *,'ERROR IN GetHVQMNRXsection: unknown XSecType ',XSecType
        call HF_stop
      endif
      if(norm_y) then
        idxyminref = GetBinIndex(IDataSet,'yminREF')
        idxymaxref = GetBinIndex(IDataSet,'ymaxREF')
        if(idxyminref.eq.0.or.idxymaxref.eq.0) then
           print *,'ERROR IN GetHVQMNRXsection:'
           print *,'no yminREF or ymaxREF bin'
           call HF_stop
        endif
      else
        idxyminref = 0
        idxymaxref = 0
      endif
      
      ! Final filling
      do i=1,NDATAPOINTS(IDataSet)
        idx=DATASETIDX(IDataSet,i)
        THEO(idx) = FindXSPtY(nbin,xsec,ymin,ymax,ptmin,ptmax,
     $    AbstractBins(idxymin,idx),AbstractBins(idxymax,idx),
     $    AbstractBins(idxptmin,idx),AbstractBins(idxptmax,idx),d_pt,d_y)
      ! normalise to another bin if needed
        if(norm_y) then
          THEO(idx) = THEO(idx)/FindXSPtY(nbin,xsec,ymin,ymax,ptmin,ptmax,
     $      AbstractBins(idxyminref,idx),AbstractBins(idxymaxref,idx),
     $      AbstractBins(idxptmin,idx),AbstractBins(idxptmax,idx),d_pt,d_y)
      ! ... or scale to fragmentation fraction
        else
          THEO(idx) = THEO(idx)*fragfraction
        endif
      enddo

      end


      subroutine GetMuPar(mu,q,A,B,C)
C---------------------------------------------------------------------------------------------
C Read parameters for perturbative scales from MINUIT extra parameters
C
C Scales for charm and beauty production are parametrised as:
C mu_f(c)^2 = MNRmf_A_c * pT_c^2 + MNRmf_B_c * m_c^2 + MNRmf_C_c
C mu_r(c)^2 = MNRmr_A_c * pT_c^2 + MNRmr_B_c * m_c^2 + MNRmr_C_c
C mu_f(b)^2 = MNRmf_A_b * pT_b^2 + MNRmf_B_b * m_b^2 + MNRmf_C_b
C mu_r(b)^2 = MNRmr_A_b * pT_b^2 + MNRmr_B_b * m_b^2 + MNRmr_C_b
C where mu_f(c), mu_r(c), mu_f(b), mu_r(b) are factorisation and renormalisation 
C scales for charm and beauty production, respectively, pT is transverse momentum 
C and m_c, m_b are charm and beauty quark masses.
C
C In total, one can provide all 12 parameters (MNRmf_A_c, MNRmf_B_c, MNRmf_C_c,
C MNRmr_A_c, MNRmr_B_c, MNRmr_C_c, MNRmf_A_b, MNRmf_B_b, MNRmf_C_b,
C MNRmr_A_b, MNRmr_B_b, MNRmr_C_b), however there are the foolowing rules:
C 1) if suffix _c (_b) at the end of variable name is omitted, the parameter is applied 
C    for both charm and beauty production
C 2) instead of providing e.g. MNRmr_A_c and MNRmr_B_c the user can provide one 
C    parameter MNRmr_AB_c so then MNRmr_A_c = MNRmr_B_c = MNRmr_AB_c
C 3) if parameters *_C_* are not provided, they are set to 0
C So e.g. for charm production only it is enough to provide just the following two parameters 
C MNRmr_AB and MNRmf_AB.
C---------------------------------------------------------------------------------------------
      implicit none
      ! input variables:
      !   mu should be either 'f' or 'r' (for factorisation or renormalisation scale, respectively)
      !   q should be either 'c' or 'b' (for charm or beauty, respectively)
      character mu,q
      ! output variables
      double precision A,B,C
      ! MINUIT indices
      integer iA,iB,iC
      ! MINUIT parameter variables
      character*256 parname
      double precision unc,ll,ul
      integer st
      ! called routine
      integer GetParameterIndex
      ! default settings
      double precision defmuA,defmuB,defmuC
      data defmuA,defmuB,defmuC /1.0,1.0,0.0/
#include "extrapars.inc"

      iA=GetParameterIndex('MNRm'//mu//'_AB')
      iB=iA
      if(iA.eq.0) then
        iA=GetParameterIndex('MNRm'//mu//'_A')
        iB=GetParameterIndex('MNRm'//mu//'_B')
        if(iA.eq.0.or.iB.eq.0) then
          iA=GetParameterIndex('MNRm'//mu//'_AB_'//q)
          iB=iA
          if(iA.eq.0) then
            iA=GetParameterIndex('MNRm'//mu//'_A_'//q)
            iB=GetParameterIndex('MNRm'//mu//'_B_'//q)
          endif
        endif
      endif
      ! parameters not in ExtraParamMinuit -> using default values
      if(iA.eq.0.or.iB.eq.0) then
        A=defmuA
        B=defmuB
      else
        iA=iExtraParamMinuit(iA)
        iB=iExtraParamMinuit(iB)
        call MNPOUT(iA,parname,A,unc,ll,ul,st)
        ! parameter in ExtraParamMinuit, nut not in MINUIT: this happens, if we are not in 'Fit' mode -> using default value
        if(st.lt.0) then
          A=defmuA
        endif
        call MNPOUT(iB,parname,B,unc,ll,ul,st)
        ! parameter in ExtraParamMinuit, nut not in MINUIT: this happens, if we are not in 'Fit' mode -> using default value
        if(st.lt.0) then
          B=defmuB
        endif
      endif
      iC=GetParameterIndex('MNRm'//mu//'_C')
      if(iC.eq.0) then
        iC=GetParameterIndex('MNRm'//mu//'_C_'//q)
      endif
      if(iC.eq.0) then
        C=defmuC
      else
        iC=iExtraParamMinuit(iC)
        call MNPOUT(iC,parname,C,unc,ll,ul,st)
        ! parameter in ExtraParamMinuit, nut not in MINUIT: this happens, if we are not in 'Fit' mode -> using default value
        if(st.lt.0) then
          C=defmuC
        endif
      endif
      end


      subroutine GetFPar(q,FFpar)
C-----------------------------------------------------------------------
C Read fragmentation parameter
C
C Parameters for non-perturbative fragmentation can be provided 
C as MINUIT extra parameters MNRfrag_c or MNRfrag_b
C for charm and beauty production, respectively.
C-----------------------------------------------------------------------
      implicit none
      ! input variable:
      !   q should be either 'c' or 'b' (for charm or beauty, respectively)
      character q
      ! output variable
      double precision FFpar
      ! MINUIT parameter index
      integer iFFpar
      ! MINUIT parameter variables
      character*256 parname
      double precision unc,ll,ul
      integer st
      ! called routine
      integer GetParameterIndex
      ! default settings
      double precision defFFc,defFFb
      data defFFc,defFFb /4.4,11.0/
#include "extrapars.inc"

      iFFpar=GetParameterIndex('MNRfrag_'//q)
      ! parameter not in ExtraParamMinuit -> using default value
      if(iFFpar.eq.0) then
        if(q.eq.'c') then
          FFpar=defFFc
        else if(q.eq.'b') then
          FFpar=defFFb
        else
          write(*,*)'Warning in GetFPar(): no default value for q = ',q
          call makenan(FFpar)
        endif
      else
        iFFpar=iExtraParamMinuit(iFFpar)
        call MNPOUT(iFFpar,parname,FFpar,unc,ll,ul,st)
        ! parameter in ExtraParamMinuit, nut not in MINUIT: this happens, if we are not in 'Fit' mode -> using default value
        if(st.lt.0) then
          if(q.eq.'c') then
            FFpar=defFFc
          else if(q.eq.'b') then
            FFpar=defFFb
          else
            write(*,*)'Warning in GetFPar(): no default value for q = ',q
            call makenan(FFpar)
          endif
        endif
      endif
      end


      function FindXSPtY(n,XS,ymin,ymax,ptmin,ptmax,
     $                   curymin,curymax,curptmin,curptmax,d_pt,d_y)
C-----------------------------------------------------------------------
C Returns cross section in required bin
C
C Input parameters:
C   n: number of bins
C   XS: cross section array (dimension n)
C   ymin: min rapidity array (dimension n)
C   ymax: max rapidity array (dimension n)
C   ptbin: min pT array (dimension n)
C   ptmax: max pT array (dimension n)
C   curymin: min rapidity of required bin
C   curymax: max rapidity of required bin
C   curptmin: min pT of required bin
C   curptmax: max pT of required bin
C   d_pt: if TRUE, cross section is divided by pT bin width
C   d_y: if TRUE, cross section is divided by rapidity bin width
C-----------------------------------------------------------------------
      implicit none
      double precision FindXSPtY
      integer n
      double precision XS(*),ymin(*),ymax(*),ptmin(*),ptmax(*)
      double precision curymin,curymax,curptmin,curptmax
      logical d_pt,d_y
      integer i
#include "ntot.inc"
#include "indata.inc"

      do i=1,n
        if(ymin(i).eq.curymin.and.ymax(i).eq.curymax.and.
     $     ptmin(i).eq.curptmin.and.ptmax(i).eq.curptmax) then
          FindXSPtY=XS(i)
          ! Check if XS is not nan
          if(FindXSPtY.ne.FindXSPtY) then
            write(*,*)'Warning in FindXSPtY(): nan cross section in bin'
            write(*,*)'y: ',curymin,curymax,' pt: ',curptmin,curptmax
            write(*,*)'Resetting to -1000000000.0'
            FindXSPtY=-1000000000.0
          endif
          ! Divide by bin width if needed
          if(d_pt) then
            FindXSPtY=FindXSPtY/(ptmax(i)-ptmin(i))
          endif
          if(d_y) then
            FindXSPtY=FindXSPtY/(ymax(i)-ymin(i))
          endif
          return
        endif
      enddo
      write(*,*)'ERROR in FindXSPtY(): bin not found'
      write(*,*)'y: ',curymin,curymax,' pt: ',curptmin,curptmax
      call HF_stop
      end
      
      
      subroutine makenan(x)
C-----------------------------------------------------------------------
C Set x to nan
C-----------------------------------------------------------------------
      implicit none
      double precision x
      
      x=0
      x=x/x
      end
