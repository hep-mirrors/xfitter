      Subroutine GetTheoryForDataset(IDataSet)
C---------------------------------------------------------------
C Created  24/05/11
C
C Distribute calculation of theory prediction for a dataset IDataSet
C---------------------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      integer IDataSet,kflag
C-------------------------------------------------------------------
      if (DATASETREACTION(IDataSet).eq.'NC e+-p integrated') then
         Call GetIntegratedNCXsection(IDataSet, HFSCHEME)
      elseif (DATASETREACTION(IDataSet).eq.'NC e+-p') then

        if (DipoleModel.gt.0.and.DipoleModel.lt.2) then
            call DipolePrediction(IDataSet)
         elseif (DipoleModel.eq.2.or.DipoleModel.eq.4) then
            Call GetNCXsection(IDataSet, HFSCHEME)
            call DipolePrediction(IDataSet)
         elseif (DipoleModel.eq.5) then
            Call DipoleBGK(IDataSet)
         else
C Standard DGLAP:
            Call GetNCXsection(IDataSet, HFSCHEME)
         endif

      elseif (DATASETREACTION(IDataSet).eq.'muon p') then
         Call GetNCXsection(IDataSet, HFSCHEME)
      elseif (DATASETREACTION(IDataSet).eq.'NC e+-p charm') then
         Call GetNCCharmXsection(IDataSet, HFSCHEME)
      elseif (DATASETREACTION(IDataSet).eq.'CC e+-p') then
         Call GetCCXsection(IDataSet, HFSCHEME)
      elseif (DATASETREACTION(IDataSet).eq.'CC pp' .or.
     $        DATASETREACTION(IDataSet).eq.'CC ppbar' ) then
         Call GetDYCCXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'NC pp' .or.
     $        DATASETREACTION(IDataSet).eq.'NC ppbar' ) then
         Call GetDYNCXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'pp jets APPLGRID') then
         Call GetJetsPPApplGrid(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'FastNLO ep jets') then
         Call GetJetsFastNLOXsection(IDataSet, .false.)
      elseif (DATASETREACTION(IDataSet).eq.'FastNLO ep jets normalised')then
cv
         call eprc_init(.true.)
         Call GetJetsFastNLOXsectionNormalised(IDataSet)
       elseif (DATASETREACTION(IDataSet).eq.'ttbar') then
         Call GetHathorXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'DDIS') then
         Call GetDiffDisXsection(IDataSet)
      else
CC         print *,'Unknown x-section type',DATASETREACTION(IDataSet)
      endif

      end

      Subroutine GetTheoryIteration
C---------------------------------------------------------------------
C
C Created 24/06/2011. Get theory calculation per iteration, before going into  individual datasets
C
C---------------------------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'endmini.inc'
C--------------------------------------------------------------------
C Drell-Yan:
      if (LFitDY) then
         call dy_do_calc
      endif

      if (LFastAPPLGRID) then
         call Calc_pdf_applgrid_fast
      endif

      if (IPARAM.eq.301) then
         call DDIS_FixModelParams(parminuitsave)
      endif


      end
