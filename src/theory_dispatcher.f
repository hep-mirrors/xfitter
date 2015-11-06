      Subroutine GetTheoryForDataset(IDataSet)
C---------------------------------------------------------------
C Created  24/05/11
C
C Distribute calculation of theory prediction for a dataset IDataSet
C---------------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
#include "scales.inc"
      integer IDataSet,kflag!
C-------------------------------------------------------------------

      if ( UseFixedTheory(IDataSet)) then
         Call UseFixedTheoryXsection(IDataSet) 
      elseif (DATASETREACTION(IDataSet).eq.'NC e+-p integrated') then         
         if(Itheory.lt.100.) then
            Call GetIntegratedNCXsection(IDataSet, HFSCHEME)
         else
           write(6,*)
     >       'NC e+-p integrated: invalid dataset for itheory>100'
            call hf_stop
         Endif
      elseif (DATASETREACTION(IDataSet).eq.'NC e+-p') then

         if (DipoleModel.eq.1.or.DipoleModel.eq.2) then
            call DipolePrediction(IDataSet)
         elseif (DipoleModel.eq.3.or.DipoleModel.eq.4) then
            Call GetNCXsection(IDataSet, HFSCHEME)
            Call DipolePrediction(IDataSet)
         elseif (DipoleModel.eq.5) then
            Call DipoleBGK(IDataSet)
         else
C Standard DGLAP:& TMDs
            Call GetNCXsection(IDataSet, HFSCHEME)
         endif

      elseif (DATASETREACTION(IDataSet).eq.'muon p') then
         if(Itheory.lt.100) then
            Call GetNCXsection(IDataSet, HFSCHEME)
         else
            write(6,*) ' muon p: invalid dataset for itheory > 100 '
            call hf_stop
         Endif
      elseif (DATASETREACTION(IDataSet).eq.'NC e+-p charm') then
         Call GetNCCharmXsection(IDataSet, HFSCHEME)
      elseif (DATASETREACTION(IDataSet).eq.'NC e+-p beauty') then
         Call GetNCBeautyXsection(IDataSet, HFSCHEME)
       elseif (DATASETREACTION(IDataSet).eq.'NC e+-p FL') then
         Call GetNCFL(IDataSet, HFSCHEME)
       elseif (DATASETREACTION(IDataSet).eq.'NC e+-p F2') then
         Call GetNCF2(IDataSet, HFSCHEME)
c           
      elseif (DATASETREACTION(IDataSet).eq.'CC e+-p integrated') then         
         if(Itheory.lt.100.) then
            Call GetIntegratedCCXsection(IDataSet, HFSCHEME)
         else
           write(6,*)
     >       'CC e+-p integrated: invalid dataset for itheory>100'
            call hf_stop
         Endif
c          
      elseif (DATASETREACTION(IDataSet).eq.'CC e+-p') then
         if(Itheory.lt.100) then
            Call GetCCXsection(IDataSet, HFSCHEME)
         else
            write(6,*) ' CC e+-p: invalid dataset for itheory > 100 '
            call hf_stop
         Endif
      elseif (DATASETREACTION(IDataSet).eq.'CC pp' .or.
     $        DATASETREACTION(IDataSet).eq.'CC ppbar' ) then
         if(Itheory.lt.100) then
           if ( DATASETTheoryType(IDataSet).eq.'expression' ) then
             !call set_theor_CKM(IDataSet,
             call get_theor_eval(IDataSet, 
     $         NDATAPOINTS(IDataSet), DATASETIDX(IDataset,1))
           else
            Call GetDYCCXsection(IDataSet)
           endif
         else
            write(6,*) ' CC pp: invalid dataset for itheory > 100 '
            call hf_stop
         Endif
      elseif (DATASETREACTION(IDataSet).eq.'NC pp' .or.
     $        DATASETREACTION(IDataSet).eq.'NC ppbar' ) then
         if(Itheory.lt.100) then
           if ( DATASETTheoryType(IDataSet).eq.'expression' ) then
             !call set_theor_CKM(IDataSet,
             call get_theor_eval(IDataSet, 
     $         NDATAPOINTS(IDataSet), DATASETIDX(IDataset,1))
           else
             call GetDYNCXsection(IDataSet)
           endif
         else
            write(6,*) ' NC ppbar: invalid dataset for itheory > 100 '
            call hf_stop
         Endif
      elseif (DATASETREACTION(IDataSet).eq.'pp jets APPLGRID') then
         if(Itheory.lt.100) then
           if ( DATASETTheoryType(IDataSet).eq.'expression' ) then
             !call set_theor_CKM(IDataSet,
             call get_theor_eval(IDataSet, 
     $         NDATAPOINTS(IDataSet), DATASETIDX(IDataset,1))
           else
             Call GetJetsPPApplGrid(IDataSet)
	   endif
         else
           write(6,*) 'pp jets APPLGRID: invalid dataset for ithory>100'
            call hf_stop
         Endif
      elseif (DATASETREACTION(IDataSet).eq.'pp jets fastNLO') then
         if(Itheory.lt.100) then
           if ( DATASETTheoryType(IDataSet).eq.'expression' ) then
             !call set_theor_CKM(IDataSet,
             call get_theor_eval(IDataSet, 
     $         NDATAPOINTS(IDataSet), DATASETIDX(IDataset,1))
           else
             Call GetJetsPPApplGrid(IDataSet)
	   endif
         else
           write(6,*) 'pp jets APPLGRID: invalid dataset for ithory>100'
            call hf_stop
         Endif
      elseif (DATASETREACTION(IDataSet).eq.'FastNLO jets' .or.
     $        DATASETREACTION(IDataSet).eq.'FastNLO ep jets') then  ! for backward compatibility
         if(Itheory.lt.100) then
            Call GetJetsFastNLOXsection(IDataSet, .false.)
         else
           write(6,*) 'FastNLO jets: invalid dataset for itheory>100'
            call hf_stop
         Endif
CMK->
      elseif (DATASETREACTION(IDataSet).eq.'FastNLO ttbar') then
         if(Itheory.lt.100) then
            Call GetTopFastNLOXsection(IDataSet, .false.)
         else
           write(6,*)
     >       'FastNLO top pairs: invalid dataset for itheory > 100'
            call hf_stop
         Endif
       elseif (DATASETREACTION(IDataSet).eq.'FastNLO ttbar normalised')
     >     then
         if(Itheory.lt.100) then
            Call GetTopFastNLOXsectionNormalised(IDataSet, .true.)
         else
           write(6,*)
     >       'FastNLO top pairs norm: invalid dataset for itheory>100'
           call hf_stop
         Endif
CMK <-

       elseif (DATASETREACTION(IDataSet).eq.'FastNLO ep jets normalised'
     >     ) then
cv
         if(Itheory.lt.100) then
            call eprc_init(.false.)
            Call GetJetsFastNLOXsectionNormalised(IDataSet)
         else
           write(6,*)
     >       'FastNLO ep jets norm.: invalid dataset for itheory>100'
            call hf_stop
         Endif
       elseif (DATASETREACTION(IDataSet).eq.'ttbar') then
         if(Itheory.lt.100) then
            Call GetHathorXsection(IDataSet)
         else
           write(6,*) 'ttbar: invalid dataset for itheory>100'
            call hf_stop
         Endif
      elseif (DATASETREACTION(IDataSet).eq.'DDIS') then
         if(Itheory.lt.100) then
            Call GetDiffDisXsection(IDataSet)
         else
           write(6,*) 'DDis: invalid dataset for ithory > 100'
            call hf_stop
         Endif

C HVQMNR for heavy-quark production in pp 
      elseif (DATASETREACTION(IDataSet).eq.'HVQMNR pp QQbar') then
         if(Itheory.lt.100) then
            Call GetHVQMNRXsection(IDataSet)
         else
            write(6,*) ' invalid dataset for ithory > 100 '
            call hf_stop
         Endif
         
       elseif ((index(DATASETREACTION(IdataSet), ' Dummy').gt.0).or.
     $     (index(DATASETREACTION(IdataSet), 'Dummy').gt.0)) then
         if(Itheory.lt.100) then
           Call GetDummyXsection(IDataSet)
         else
           write(6,*) ' Dummy: invalid dataset for itheory > 100 '
            call hf_stop
          Endif


      else
            Call hf_errlog(01110113,'F: theory_disp: unknown reaction "'
     $                //TRIM(DATASETREACTION(IDataSet)) //'"')
      endif


      end

      Subroutine GetTheoryIteration
C---------------------------------------------------------------------
C
C Created 24/06/2011. Get theory calculation per iteration, before going into  individual datasets
C
C---------------------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "endmini.inc"
#include "couplings.inc"
C--------------------------------------------------------------------
C Drell-Yan:
      if (LFitDY) then
         call dy_do_calc
      endif

      if (.not.LUseAPPLgridCKM) then
        call update_theor_ckm
      endif

      if (LFastAPPLGRID) then
         call Calc_pdf_applgrid_fast
      endif

      if (PDFStyle.eq.'DDIS') then
         call DDIS_FixModelParams(parminuitsave)
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
