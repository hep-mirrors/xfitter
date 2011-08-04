
      subroutine init_theory_modules
*     ------------------------------------------------

      implicit none 
      include 'steering.inc'
      include 'thresholds.inc'
      include 'couplings.inc'


*     ------------------------------------------------
*     Initialise EW parameters
*     ------------------------------------------------

      call Init_EW_parameters

*     ------------------------------------------------
*     Initialise qcdnum
*     ------------------------------------------------

      if(itheory.eq.0) then
         call qcdnum_ini
      elseif(itheory.eq.1) then       
c          here goes a call to non-DGLAP 
      endif

*     ------------------------------------------------
*     Initialise calculations for each dataset:
*     ------------------------------------------------

      call Init_theory_datasets


      return
      end


*     ------------------------------------------------
*     ------------------------------------------------





      subroutine qcdnum_ini

*     ------------------------------------------------
*     QCDNUM initialisation
*     ------------------------------------------------

      implicit none

      INCLUDE 'steering.inc'
      INCLUDE 'thresholds.inc'
      INCLUDE 'couplings.inc'
      INCLUDE 'datasets.inc'

c set-up of the constants
      integer iord
      integer iosp,nqout
 

C RT parameters:
      double precision alphaS0in,alambdain,flavorin,qsctin,qsdtin
      integer iordin,inullin


      double precision xmin(5)
      integer  iwt(5)
      integer NQGRID, NXGRID
      PARAMETER (NQGRID=11)
      PARAMETER (NXGRID=5)
      

      double precision  QARR(NQGRID),WGT(NQGRID)
      data iosp/2/                   !x grid, lin/qua/spli
      DATA WGT/1.d0,  1.d0,  1.d0,  1.d0, 1.d0,
     $     1.d0, 1.d0, 2.d0, 2.d0, 1.d0, 2.d0/
      DATA QARR/1., 1.1, 1.2,1.6, 1.8, 1.9,
     $     1.96,2.89,22.5625,30625.,500000000./


c      DATA WGT/2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 
c     $         2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 
c     $         2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 
c     $         2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 
c     $         2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 
c     $         2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 
c     $         2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 
c     $         2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 
c     $         2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0, 
c     $         2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 
c     $         2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 
c     $         2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,2.d0, 
c     $         2.d0, 2.d0,2.d0,2.d0, 2.d0,2.d0,4.d0, 4.d0,4.d0,4.d0, 
c     $         4.d0, 4.d0,4.d0,4.d0, 4.d0,4.d0,4.d0, 4.d0/


c      DATA QARR/ 1.000, 1.020, 1.040, 1.060, 1.080, 1.100, 1.120,
c     $           1.140, 1.160, 1.180, 1.200, 1.220, 1.240, 1.260,
c     $           1.280, 1.300, 1.320, 1.340, 1.360, 1.380, 1.400,
c     $           1.420, 1.440, 1.460, 1.480, 1.500, 1.520, 1.540,
c     $           1.560, 1.580, 1.600, 1.620, 1.640, 1.660, 1.680,
c     $           1.700, 1.720, 1.740, 1.760, 1.780, 1.800, 1.8225,
c     $           1.860, 1.880, 1.900, 1.920, 1.940, 1.960,
c     $           1.980, 2.000, 2.020, 2.040, 2.060, 2.080, 2.100,
c     $           2.120, 2.140, 2.160, 2.180, 2.200, 2.220, 
c     $           2.250, 2.280, 2.300, 2.320, 2.340, 2.360, 2.380,
c     $           2.400, 2.420, 2.440, 2.460, 2.480, 2.500, 2.520,
c     $           2.540, 2.560, 2.580, 2.600, 2.620, 2.640, 2.660,
c     $           2.680, 2.700, 2.7225, 2.760, 2.780, 2.800,
c     $           2.820, 2.840, 2.860, 2.880, 2.900, 2.920, 2.940,
c     $           2.960, 2.980, 3.000, 3.020, 3.040, 3.060, 3.080,
c     $           3.100, 3.120, 3.140, 3.160, 3.180, 3.200, 3.220,
c     $           3.240, 3.260, 3.280, 3.300, 3.320, 3.340, 3.360,
c     $           3.380, 3.400, 3.420, 3.440, 3.460, 3.480, 3.500,
c     $           4.000, 6.5,   12.00,   18.49, 20.25, 25.00, 35.,
c     $           60., 120., 200., 400., 1000., 30625., 100000./



c--   linear x grid
      data xmin/9.9d-7,0.01d0,0.10d0,0.40d0,0.70d0/
      data iwt/1,2,4,8,16/






      double precision qmas, hqmass
      dimension qmas(3), hqmass(3)


      double precision as0,r20
      data as0/0.364/, r20/2.D0/!, nfin/0/ !alphas, NNLO, VFNS
      integer I,ndum,ierr

      double precision a,b
      double precision aq2,bq2,qs0,qt

      integer id1,id2,iq0,iqb,iqc,iqt

      integer nflav
      integer nw,nwords,nx

C Functions:
      integer iqfrmq
      double precision asfunc
C---------------------------------------------------------------------------------------

      q0 = starting_scale
      qc = HF_MASS(1)**2
      qb = HF_MASS(2)**2

      call qcinit(6,' ')        !initialize
      call setord(I_FIT_ORDER)         !LO, NLO, NNLO


      call gxmake(xmin,iwt,nxgrid,200,nx,iosp)        !x-grid
      call gqmake(qarr,wgt,nqgrid,120,nqout)          !mu2-grid
      iq0 =iqfrmq(q0)
      iqc =iqfrmq(qc)
      iqb =iqfrmq(qb)



!SG: hardwire top mass

      qt = 200.**2
      iqt =iqfrmq(qt)


      call setcbt(0,iqc,iqb,iqt) !thesholds in the ffns
      call setalf(dble(rtalphas), Mz*Mz) !input alphas

      call readwt(22,'unpolarised.wgt',id1,id2,nw,ierr)
      if(ierr.ne.0) then
         call fillwt(0,id1,id2,nw) !calculate weights
cv         call dmpwgt(1,22,'unpolarised.wgt')
      else 
         print*,' ERRRROR in read unpolarised wieght', ierr
      endif
      write(6,'(/'' weight: words used ='',I10)') nw     

      call zmreadw(22,'zmstf.wgt',nwords,ierr)
      if(ierr.ne.0) then
         call zmfillw(nwords)
cv         call zmdumpw(22,'zmstf.wgt')
      else 
         print*,' ERRRROR in read zmstf weight', ierr
      endif      
      write(6,'(/'' ZMSTF: words used ='',I10)') nwords      

! setting of the evolution parameters


c      call SETABR(1.D0,0.D0)  ! mur scale variation
c      call ZMDEFQ2(1.D0,0.D0) ! muf scale variation

      if ((mod(HFSCHEME,10).eq.2).or.(vfnsINDX.eq.2)) then
         alambdaIn=0.307
         qs0=1.d0
         alphas0in =asfunc(qs0,nflav,ierr)
         qsdtIn = 4.d0 * qc
         qsctIn = 4.d0 * qb
         flavorIn = 3
         iordIn = 1
         inullIn=0
C-
         call RT_Set_Input(
     $        alphaS0in,alambdain,flavorin,qsctin,qsdtin,iordin,inullin)
         call WATE96

      endif
      
      print*,'exit qcdnum_Ini'
      return
      end


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
         elseif (DATASETREACTION(IDataSet).eq.'FastNLO ep jets') then
            Call InitJetsFastNLODataSet(IDataSet)
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
      Call InitDYNCXsection
      Call InitJetsPPApplGrid
      Call InitJetsFastNLO
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

      subroutine InitDYNCXsectionDataset_kfactor(IDataSet)
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
      double precision ranges(7)

      integer NPmax
      parameter(NPmax=100)
      double precision yb(Npmax+1)

      integer idx, idxY1, idxY2,i 

C----------------------------------------------------------


      if (NDATAPOINTS(IDataSet).gt.NPmax) then
         print *,'ERROR IN InitDYNCXsectionDataset'
         print *,'INCREASE NPMax to ',NDATAPOINTS(IDataSet)                                                                                                  
         stop
      endif

C Set global parameter:
      LFitDY = .true.

      ! pt
      ranges(7)=20.d0
      ! mass
      ranges(1) = 66.d0
      ranges(2) = 116.d0
      
      ! rap
      ranges(3) = -10.d0
      ranges(4) =  10.d0
      ! eta
      ranges(5) = -10.d0
      ranges(6) =  10.d0

C Get indicies:
      idxY1 = GetBinIndex(IDataSet,'y1')
      idxY2 = GetBinIndex(IDataSet,'y2')

      if (idxY1.eq.0 .or. idxY2.eq.0) then
         print 
     $        '(''ERROR in GetDYNCXsection, can not find bin index for y1, y2'',2i6)'
     $        ,idxY1,idxY2
         stop
      endif

C Define bins:
      do i=1,NDATAPOINTS(IDataSet) 
         idx =  DATASETIDX(IDataSet,i)
         yb(i)   =  AbstractBins(idxY1,idx)
         yb(i+1) =  AbstractBins(idxY2,idx)
      enddo

      print *,'Initialise DY calculations for dataset', IDataSet
      call dy_create_calc(IDataSet, 1, 7000d0, 'Z'//char(0), ranges, 
     $   'y'//char(0), NDATAPOINTS(IDataSet), yb)
      
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
      double precision ranges(7)

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

      ! pt
      ranges(7)=25.d0
      ! mass
      ranges(1) = 1.d0
      ranges(2) = 7000.d0
      
      ! rap
      ranges(3) = -10.d0
      ranges(4) =  10.d0
      ! eta
      ranges(5) = -10.d0
      ranges(6) =  10.d0

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
      call dy_create_calc(IDataSet, 1, 7000d0, 'W'//char(0), ranges, 
     $   'eta'//char(0), NDATAPOINTS(IDataSet), eb)
      
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

      subroutine InitJetsFastNLODataSet(IDataSet)
C------------------------------------------------------------
C
C Initialize FastNLO reader
C
C------------------------------------------------------------
      implicit none
      integer IDataSet
      include 'steering.inc'
      include 'datasets.inc'
      call fastnloinit(DATASETLABEL(IDataSet),IDataSet);
      end


      subroutine InitIntegratedNCXsection
      end

      subroutine InitReducedNCXsection
      end

      subroutine  InitCCXsection
      end

      subroutine InitDYCCXsection
      end

      subroutine InitDYNCXsection
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

      subroutine InitJetsFastNLO
      end


      Subroutine Init_EW_parameters
C-----------------------------------------------------
C
C  Initialise electroweak parameters. Created 5 June 2011
C
C-----------------------------------------------------
      implicit none
      include 'couplings.inc'
      namelist/EWpars/alphaem, gf, sin2thw, convfac,
     $ Mz, Mw, Mh, wz, ww, wh, wtp,
     $ Vud, Vus, Vub, Vcd, Vcs, Vcb,
     $ men, mel, mmn, mmo, mtn, mta, mup, mdn,
     $ mch, mst, mtp, mbt
C-----------------------------------------------------

      open (51,file='ewparam.txt',status='old')
      read (51,NML=EWpars,END=71,ERR=72)
      close (51)

      goto 73
C 
 71   continue
      print '(''Namelist @EWPars NOT found'')'
      goto 73
 72   continue
      print '(''Error reading namelist @EWPars, STOP'')'
      stop
 73   continue

C
C set derived values
C
      cos2thw = 1.d0 - sin2thw

C
C set same values for DY calculations
C
      call dy_set_ewpars

cv electroweak starting values, modified later
c      cvu = 0.196
c      cau = 0.5
c      cvd = -0.346
c      cad = -0.5
      
C  DELTA-R AND EFFECTIVE WEAK MIXING ANGLE FROM ZNCV
C      call EPRC_INIT(.true.)

C-----------------------------------------------------
      end
