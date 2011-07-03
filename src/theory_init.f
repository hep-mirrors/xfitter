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
         elseif (DATASETREACTION(IDataSet).eq.'NC pp') then
            Call InitDYNCXsectionDataset(IDataSet)
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
C------------------------------------------------------------
C
C Initialise tables for DY process for calculations
C
C------------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      integer IDataSet
C---------------------------------------------------------


      if (DATASETTheoryType(IDataSet).eq.'kfactor') then
         call InitDYCCXsectionDataset_kfactor(IDataSet)
      elseif (DATASETTheoryType(IDataSet).eq.'applgrid') then
         call InitDYXsectionDataset_applgrid(IDataSet)         
      else
         print *,'InitDYCCXsectionDataset: unknown theory type'
     $        ,DATASETTheoryType(IDataSet), ' for set ', IDataSet
         print *,'stop'
         stop
      endif

      end


      subroutine InitDYNCXsectionDataset(IDataSet)
C------------------------------------------------------------
C
C Initialise tables for DY process for calculations
C
C------------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      integer IDataSet
C---------------------------------------------------------


      if (DATASETTheoryType(IDataSet).eq.'kfactor') then
         call InitDYNCXsectionDataset_kfactor(IDataSet)
      elseif (DATASETTheoryType(IDataSet).eq.'applgrid') then
         call InitDYXsectionDataset_applgrid(IDataSet)         
      else
         print *,'InitDYNCXsectionDataset: unknown theory type'
     $        ,DATASETTheoryType(IDataSet), ' for set ', IDataSet
         print *,'stop'
         stop
      endif

      end

      subroutine InitDYNCXsectionDataset_kfactor
      end

      subroutine InitDYXsectionDataset_applgrid(IDataSet)
C------------------------------------------------------------
C
C Initialise tables for DY process for calculations using applgrid
C
C------------------------------------------------------------     
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      integer IDataSet, IGridID
C---------------------------------------------------------------
      call appl_readgrid(IGridID,DATASETTheoryFile(IDataSet))
C Store index:
      DATASETTheoryIndex(IDataSet) = IGridID


      end

      subroutine InitDYCCXsectionDataset_kfactor(IDataSet)
C------------------------------------------------------------
C
C Initialise tables for DY process for calculations using k-factors
C
C------------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      integer IDataSet

      integer GetBinIndex                                                                                                                                    
      double precision dy_mass(2)                                                                                                                            
      double precision dy_y(2)
      double precision pt_cut
      integer NPmax
      parameter(NPmax=100)
      double precision eb(Npmax+1)

      integer idx, idxEta1, idxEta2,i 

C----------------------------------------------------------


      if (NDATAPOINTS(IDataSet).gt.NPmax) then
         print *,'ERROR IN InitDYCCXsectionDataset'
         print *,'INCREASE NPMax to ',NDATAPOINTS(IDataSet)                                                                                                  
         stop
      endif

C Set global parameter:
      LFitDY = .true.

      pt_cut=25.
      dy_mass(1) = 1.
      dy_mass(2) = 7000.
      
      dy_y(1) = -10.
      dy_y(2) =  10.

C Get indicies:
      idxEta1 = GetBinIndex(IDataSet,'eta1')
      idxEta2 = GetBinIndex(IDataSet,'eta2')

      if (idxEta1.eq.0 .or. idxEta2.eq.0) then
         print 
     $        '(''ERROR in GetDYCCXsection, can not find bin index for Eta1, Eta2'',2i6)'
     $        ,idxEta1,idxEta2
         stop
      endif

C Define bins:
      do i=1,NDATAPOINTS(IDataSet) 
         idx =  DATASETIDX(IDataSet,i)
         eb(i)   =  AbstractBins(idxEta1,idx)
         eb(i+1) =  AbstractBins(idxEta2,idx)
      enddo

      print *,'Initialise DY calculations for dataset', IDataSet
      call dy_create_calc(IDataSet, 1, 7000d0, 'W'//char(0), dy_mass, dy_y, pt_cut)
      call dy_set_bins(IDataSet,'eta'//char(0), NDATAPOINTS(IDataSet), eb)

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

      integer NfnloGrids
      parameter (NfnloGrids=100)
      integer IGridIDfnlo(Nfnlogrids)

      integer IGridID
      integer i,ibin,idx,n,n2
      integer idxPt1,idxPt2
      integer idxEta,idxEtaAdd
      double precision ptLow, ptHigh, Center
      double precision ptLowAP, ptHighAP
C functions:

      integer GetBinIndex
      integer appl_getnbins, appl_getbinnumber
      double precision appl_getbinlowedge,appl_getbinwidth

      character *80 ct

C------------------------------------------------------------
      
      
      if (DATASETTheoryType(IDataSet).eq.'FastNLO') then
         call appl_ngrids(n)
         ct = DATASETTheoryFile(IDataSet)
         call appl_readfastnlogrids(IGridIDfnlo,ct(1:Index(ct,' ')-1)//char(0))
         call appl_ngrids(n2)

         if (n2-n+1.gt. Nfnlogrids) then
            print *,'ERROR in InitJetsPPApplGridDataSet'
            print *,'INCREASE Nfnlogrids to',n2-n+1
            stop
         endif

         print *,'ho',n,n2
C Store first index:
         DATASETTheoryIndex(IDataSet) = IGridIDfnlo(1)
         IGridID = IGridIDfnlo(1)

      else
         call appl_readgrid(IGridID,DATASETTheoryFile(IDataSet))
C Store index:
         DATASETTheoryIndex(IDataSet) = IGridID
      endif


C Do some checks:
C      print *,'ho',appl_getnbins(IGridID)

C Get low/high Pt indicies
      idxPt1 = GetBinIndex(IDataSet,'pt1')
      idxPt2 = GetBinIndex(IDataSet,'pt2')
      idxEta = GetBinIndex(IDataSet,'EtaBinNumber')

      do i=1,NDATAPOINTS(IDataSet)
         idx = DATASETIDX(IDataSet,i)
         ptLow  = AbstractBins(idxPt1,idx)
         ptHigh = AbstractBins(idxPt2,idx)
         center = 0.5*(ptLow + ptHigh)  

         if (idxEta.gt.0) then
            idxEtaAdd = int(AbstractBins(idxEta,idx)+0.1)-1 ! Eta index
         else
            idxEtaAdd = 0
         endif

         ibin = appl_getbinnumber(IGridID+idxEtaAdd,center)
C Store index:
         IndexTheoryBin(idx) = ibin + 1 ! +1 due to C->Fortran

C Check consistency:
         ptLowAP  = appl_getbinlowedge(IGridID+idxEtaAdd,ibin)
         ptHighAP = ptLowAP + appl_getbinwidth(IGridID+idxEtaAdd,ibin)
         
c         print *,'hoho',ibin,idxEtaAdd,ptLow, ptHigh, ptLowAP, ptHighAP
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
c      stop
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
