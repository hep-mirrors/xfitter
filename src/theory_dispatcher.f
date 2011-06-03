      Subroutine GetTheoryForDataset(IDataSet)
C---------------------------------------------------------------
C Created by SG 24/05/11
C
C Distribute calculation of theory prediction for a dataset IDataSet
C---------------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      integer IDataSet
C-------------------------------------------------------------------
      if (DATASETREACTION(IDataSet).eq.'NC e+-p integrated') then
         Call GetIntegratedNCXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'NC e+-p') then
         Call GetReducedNCXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'CC e+-p') then
         Call GetCCXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'CC pp') then
         Call GetDYCCXsection(IDataSet)
      elseif (DATASETREACTION(IDataSet).eq.'pp jets APPLGRID') then
         Call GetJetsPPApplGrid(IDataSet)
      else
CC         print *,'Unknown x-section type',DATASETREACTION(IDataSet)
      endif

      end
