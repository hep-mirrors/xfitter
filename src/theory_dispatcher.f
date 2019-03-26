      Subroutine GetTheoryForDataset(IDataSet)
C---------------------------------------------------------------
C Created  24/05/11
C
C Distribute calculation of theory prediction for a dataset IDataSet
C---------------------------------------------------------------
      implicit none
#include "ntot.inc"
C #include "steering.inc"
C #include "for_debug.inc"
#include "datasets.inc"
C #include "scales.inc"
C #include "alphas.inc"
C #include "couplings.inc"
      integer IDataSet
      double precision asref,HF_Get_alphas
C-------------------------------------------------------------------
      if ( UseFixedTheory(IDataSet)) then
         Call UseFixedTheoryXsection(IDataSet) 
      elseif ( DATASETTheoryType(IDataSet).eq.'expression' ) then  ! 6/04/17 This has priority

         call get_theor_eval(IDataSet, 
     $        NDATAPOINTS(IDataSet), DATASETIDX(IDataset,1))

      else
            Call hf_errlog(01110113,'F: theory_disp: 
     $  unknown or old reaction "'
     $                //TRIM(DATASETREACTION(IDataSet)) //'"')
      endif

      end
      !> Copy theo_fix to theory for a dataset
      Subroutine UseFixedTheoryXsection(ISet)
      implicit none
#include "ntot.inc"
#include "datasets.inc"
#include "theo.inc"
      integer i,idx,ISet
C-------------------------------------------------------
      do i=1,NDATAPOINTS(Iset)
         idx = DatasetIdx(Iset,i)
         theo(idx) = theo_fix(idx)
      enddo
C-------------------------------------------------------
      end
