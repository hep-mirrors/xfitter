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
      if (DATASETREACTION(IDataSet).eq.'NC e+-p integrated') then
         Call GetIntegratedNCXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'NC e+-p') then
         Call GetReducedNCXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'muon p') then
         Call GetReducedNCXsection(IDataSet)
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
         Call GetJetsFastNLOXsection(IDataSet)
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
C--------------------------------------------------------------------
C Drell-Yan:
      if (LFitDY) then
         call dy_do_calc
      endif

      if (LFastAPPLGRID) then
         call Calc_pdf_applgrid_fast
      endif
      end
