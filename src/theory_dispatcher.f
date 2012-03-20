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
      integer IDataSet
C-------------------------------------------------------------------
      write (*,*)' GetTheoryForDataset',IDataSet,
     $     DATASETREACTION(IDataSet)
      if (DATASETREACTION(IDataSet).eq.'NC e+-p integrated') then
         Call GetIntegratedNCXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'NC e+-p') then
         if (DipoleModel.gt.0.and.DipoleModel.le.2) then
            call DipolePrediction(IDataSet)
         elseif (DipoleModel.gt.2) then
            Call GetNCXsection(IDataSet)
            call DipolePrediction(IDataSet)
         else
            Call GetNCXsection(IDataSet)
         endif
      elseif (DATASETREACTION(IDataSet).eq.'muon p') then
         Call GetNCXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'NC e+-p charm') then
         Call GetNCCharmXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'CC e+-p') then
         Call GetCCXsection(IDataSet)
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
      elseif (DATASETREACTION(IDataSet).eq.'FastNLO ep jets normalised') then
         Call GetIntegratedNCXsection(IDataSet)
         Call GetJetsFastNLOXsection(IDataSet, .true.)
      elseif (DATASETREACTION(IDataSet).eq.'ttbar') then
         Call GetHathorXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'DDIS') then
cws         write (*,*)'IDataSet1  ',IDataSet
         Call GetDiffDisXsection(IDataSet)
cws         write (*,*)'IDataSet2  ',IDataSet

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
