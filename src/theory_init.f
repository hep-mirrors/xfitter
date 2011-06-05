      Subroutine init_theory_datasets
C---------------------------------------------------------------
C
C June 2, 2011, Initialise theory for different data and theory models
C
C---------------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
C-----------------------------------
      integer IDataSet
C---------------------------------------------------------------

C
C Loop over datasets, add information:
C
      print '(''Initialize theory for datasets'')'
      do IDataSet=1,Ndatasets
         if (DATASETREACTION(IDataSet).eq.'NC e+-p integrated') then
            Call InitIntegratedNCXsectionDataset(IDataSet)
         elseif (DATASETREACTION(IDataSet).eq.'NC e+-p') then
            Call InitReducedNCXsectionDataset(IDataSet)
         elseif (DATASETREACTION(IDataSet).eq.'CC e+-p') then
            Call InitCCXsectionDataset(IDataSet)
         elseif (DATASETREACTION(IDataSet).eq.'CC pp') then
            Call InitDYCCXsectionDataset(IDataSet)
         elseif (DATASETREACTION(IDataSet).eq.'pp jets APPLGRID') then
            Call InitJetsPPApplGridDataSet(IDataSet)
         else
C     C         print *,'Unknown x-section type',DATASETREACTION(IDataSet)
         endif
      enddo
C
C Post init for all theory models:
C
      Call InitIntegratedNCXsection
      Call InitReducedNCXsection
      Call InitCCXsection
      Call InitDYCCXsection
      Call InitJetsPPApplGrid
C---------------------------------------------------------
      end

      subroutine InitIntegratedNCXsectionDataset(IDataSet)
      end

      subroutine InitReducedNCXsectionDataset(IDataSet)
      end

      subroutine InitCCXsectionDataset(IDataSet)
      end

      subroutine InitDYCCXsectionDataset(IDataSet)
      end

      subroutine InitJetsPPApplGridDataSet(IDataSet)
C------------------------------------------------------------
C
C Add applgrid file to the list of files
C
C------------------------------------------------------------
      implicit none
      integer IDataSet
C------------------------------------------------------------
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      integer IGridID
      integer i,ibin,idx
      integer idxPt1,idxPt2
      double precision ptLow, ptHigh, Center
      double precision ptLowAP, ptHighAP
C functions:

      integer GetBinIndex
      integer appl_getnbins, appl_getbinnumber
      double precision appl_getbinlowedge,appl_getbinwidth
C------------------------------------------------------------
      
      call appl_readgrid(IGridID,DATASETTheoryFile(IDataSet))

C Store index:
      DATASETTheoryIndex(IDataSet) = IGridID

C Do some checks:
C      print *,'ho',appl_getnbins(IGridID)

C Get low/high Pt indicies
      idxPt1 = GetBinIndex(IDataSet,'pt1')
      idxPt2 = GetBinIndex(IDataSet,'pt2')

      do i=1,NDATAPOINTS(IDataSet)
         idx = DATASETIDX(IDataSet,i)
         ptLow  = AbstractBins(idxPt1,idx)
         ptHigh = AbstractBins(idxPt2,idx)
         center = 0.5*(ptLow + ptHigh)  
         ibin = appl_getbinnumber(IGridID,center)
C Store index:
         IndexTheoryBin(idx) = ibin + 1  ! +1 due to C->Fortran

C Check consistency:
         ptLowAP  = appl_getbinlowedge(IGridID,ibin)
         ptHighAP = ptLowAP + appl_getbinwidth(IGridID,ibin)

C         print *,'hoho',ibin,ptLow, ptHigh, ptLowAP, ptHighAP
      enddo

C------------------------------------------------------------      
      end


      subroutine InitIntegratedNCXsection
      end

      subroutine InitReducedNCXsection
      end

      subroutine  InitCCXsection
      end

      subroutine InitDYCCXsection
      end

      subroutine InitJetsPPApplGrid
C-------------------------------------
      implicit none
      integer n
C-------------------------------------
      call appl_ngrids(n)
      print *, n, ' applgrid grids have been read'
C-------------------------------------
      end


      Subroutine Init_EW_parameters
C-----------------------------------------------------
C
C  Initialise electroweak parameters. Created 5 June 2011
C
C-----------------------------------------------------
      implicit none
      include 'couplings.inc'
C-----------------------------------------------------

C
C Masses:
C
      Mz = 91.187d0
      Mw = 80.41d0
      sin2thw = 1.d0 - Mw**2/Mz**2

cv use mandy's
cv      sin2thw = 0.2315
      cos2thw = 1.d0 - sin2thw


cv electroweak starting values, modified later

      cvu = 0.196
      cau = 0.5
      cvd = -0.346
      cad = -0.5
      
C  DELTA-R AND EFFECTIVE WEAK MIXING ANGLE FROM ZNCV
      call EPRC_INIT(.true.)

C-----------------------------------------------------
      end
