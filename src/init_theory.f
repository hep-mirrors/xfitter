
      subroutine init_theory_modules
*     ------------------------------------------------

      implicit none 
      include 'ntot.inc'
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
         if (ewfit.gt.0) call eprc_init(.true.)
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

      include 'ntot.inc'
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
      PARAMETER (NQGRID=2)
      PARAMETER (NXGRID=5)
      

      integer NQGridHF  !> Increase number of Q2 grid points to ensure that HF thresholds and the starting
                        !> Q2 scale are in.
      parameter (NQGridHF=NQGrid+4)  

      integer NQall     !> Actual number of grid points
      double precision Q2Grid(NQGridHF)
      double precision WQGrid(NQGridHF)

      double precision  QARR(NQGRID),WGT(NQGRID)
      data iosp/2/                   !x grid, lin/qua/spli
      DATA WGT/1.d0,  2.d0/
      DATA QARR/1.,   500000000./


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






      double precision  tmp


      double precision as0,r20
      data as0/0.364/, r20/2.D0/!, nfin/0/ !alphas, NNLO, VFNS
      integer I,ndum,ierr,j

      double precision a,b
      double precision qs0,qt

      integer id1,id2,iq0,iqb,iqc,iqt

      integer nflav
      integer nw,nwords,nx

      double precision hqmass
      dimension hqmass(3)


C Functions:
      integer iqfrmq
      double precision asfunc
C---------------------------------------------------------------------------------------

      q0 = starting_scale
      qc = HF_MASS(1)**2
      qb = HF_MASS(2)**2
      qt = HF_MASS(3)**2

C Check that scales are in proper order:
      if (q0.gt.qc) then
         print *,'Starting scale must be below charm threshold, stop'
         stop
      endif
      
      if (qc.gt.qb) then
         print *,'Charm mass must be below bottom mass, stop'
         stop
      endif
      
      if (qb.gt.qt) then
         print *,'Bottom mass must be below top mass, stop'
      endif

      do i=1,NQGrid 
         Q2Grid(i) = QARR(i)
         WQGrid(i) = WGT(i)
      enddo
C Add extra points:
      Q2Grid(NQGrid+1) = q0
      Q2Grid(NQGrid+2) = qc
      Q2Grid(NQGrid+3) = qb
      Q2Grid(NQGrid+4) = qt
      WQGrid(NQGrid+1) = 1
      WQGrid(NQGrid+2) = 1
      WQGrid(NQGrid+3) = 1
      WQGrid(NQGrid+4) = 1

C Sort the Q2Grid:
      do i=1,NQGridHF
         do j=i+1,NQGridHF
            if ( Q2Grid(j) .lt. Q2GRID(i) ) then
               tmp =  Q2GRID(j)
               Q2Grid(j) = Q2GRID(i)
               Q2GRID(i) = tmp
               tmp =  WQGRID(j)
               WQGrid(j) = WQGRID(i)
               WQGRID(i) = tmp               
            endif
         enddo
      enddo

      NQAll = NQGridHF


C Remove duplicates:
      i = 1
      do while (i.lt.NQAll)
         if (  abs(Q2Grid(i)-Q2Grid(i+1)).lt.1.D-5 ) then
            do j=i+1,NQAll
               Q2Grid(j) = Q2Grid(j+1)
               WQGrid(j) = WQGrid(j+1)
            enddo
            NQAll = NQAll-1
         else
            i = i + 1
         endif
      enddo

      print *,' '
      print *,'Info FROM QCDNUM_INI'
      print '('' Init Q2 grid with number of nodes='',i5)',NQALL      
      print '('' Q2 values at:'',20F14.2)',(Q2grid(i),i=1,NQALL)
      print '('' Weights are :'',20F14.2)',(WQgrid(i),i=1,NQALL)
      print *,' '


      call qcinit(6,' ')        !initialize
      call setord(I_FIT_ORDER)         !LO, NLO, NNLO


      call gxmake(xmin,iwt,nxgrid,200,nx,iosp)        !x-grid
      call gqmake(Q2Grid,WQGrid,NQAll,120,nqout)          !mu2-grid
      iq0 =iqfrmq(q0)
      iqc =iqfrmq(qc)
      iqb =iqfrmq(qb)



 !> Top:
      iqt =iqfrmq(qt)


      if ((mod(HFSCHEME,10).eq.3)) then
         call setcbt(3,iqc,iqb,iqt) !thesholds in the ffns
         print *,'Fixed Flavour Number Scheme set with nf=3'
      else
        call setcbt(0,iqc,iqb,iqt) !thesholds in the ffns
      endif


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
      call zswitch(IPDFSET)


! setting of the evolution parameters


c      call SETABR(1.D0,0.D0)  ! mur scale variation
c      call ZMDEFQ2(1.D0,0.D0) ! muf scale variation

      if ((mod(HFSCHEME,10).eq.2)) then
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


c Fixed Flavour Number Scheme (FFNS)
      elseif ((mod(HFSCHEME,10).eq.3)) then
        if(I_FIT_ORDER.gt.2) then
          print *,'FFN scheme can be used only with NLO, stop'
          stop
        endif

         hqmass(1) = HF_MASS(1)
         hqmass(2) = HF_MASS(2)
         hqmass(3) = HF_MASS(3)

C--- aq=1.0 and bq=0 sets the heavy quarks factorisation scale
C--- Q^2 = aq2*mu_f + bq2  

         if(massh.eq.1) then
             bq2  = bq2 * qc
         elseif(massh.eq.2) then
             bq2  = bq2 * qb
         endif
c         print*,'1 HQ scale (Q^2=a*mu_F^2 + b) a,b,mh', aq2,bq2,massh 

         call hqreadw(22,'hqstf.wgt',nwords,ierr)
         if(ierr.ne.0) then
            call hqfillw(3,hqmass,aq2,bq2,nwords)
cv            call hqdumpw(22,'hqstf.wgt')
         else 
            print*,'ERRRROR in hqreadw!', ierr
         endif      
         write(6,'(/'' HQSTF: words used ='',I10)') nwords      
         call hswitch(IPDFSET)

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
      include 'ntot.inc'
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
         elseif (DATASETREACTION(IDataSet).eq.'CC pp' .or.
     $	         DATASETREACTION(IDataSet).eq.'CC ppbar' ) then
            Call InitDYCCXsectionDataset(IDataSet)
         elseif (DATASETREACTION(IDataSet).eq.'NC pp' .or.
     $           DATASETREACTION(IDataSet).eq.'NC ppbar' ) then
            Call InitDYNCXsectionDataset(IDataSet)
         elseif (DATASETREACTION(IDataSet).eq.'pp jets APPLGRID') then
            Call InitJetsPPApplGridDataSet(IDataSet)
         elseif (DATASETREACTION(IDataSet).eq.'FastNLO ep jets') then
            Call InitJetsFastNLODataSet(IDataSet)
         elseif (DATASETREACTION(IDataSet).eq.'ttbar') then
            Call InitHathorDataSet(IDataSet)
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
      Call InitHathor
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
      include 'ntot.inc'
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
      include 'ntot.inc'
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
      include 'ntot.inc'
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'indata.inc'
      integer IDataSet

      integer GetBinIndex                                                                                                                                    
      integer GetInfoIndex
      double precision ranges(7)

      integer NPmax
      parameter(NPmax=100)
      double precision sqrtS
      double precision yb(Npmax+1)
      double precision y1, y2
      integer bchpr

      integer idx, idxY1, idxY2, idxSqrtS, idxBchpr, i 
      integer idxPte, idxMinv1, idxMinv2

C----------------------------------------------------------


      if (NDATAPOINTS(IDataSet).gt.NPmax) then
         print *,'ERROR IN InitDYNCXsectionDataset'
         print *,'INCREASE NPMax to ',NDATAPOINTS(IDataSet)                                                                                                  
         stop
      endif

C Set global parameter:
      LFitDY = .true.
      
      ! beam 
      idxSqrtS = GetInfoIndex(IDataSet, 'sqrt(S)')
      sqrtS = 7000d0 ! defaults to LHC
      if ( idxSqrtS .ne. 0 ) sqrtS = DATASETInfo(idxSqrtS, IDataSet)
      ! beam charge product
      bchpr = 1 ! defaults to LHC
      if ( DATASETREACTION(IDataSet) .eq. 'NC ppbar' ) 
     $ bchpr = -1

      ! pt
      idxPte = GetInfoIndex(IDataSet,'pte cut')
      ranges(7) = 0.d0
      if ( idxPte .ne. 0 ) 
     $   ranges(7) = DATASETInfo(idxPte, IDataSet)

      ! mass
      idxMinv1 = GetInfoIndex(IDataSet, 'Minv min')
      idxMinv2 = GetInfoIndex(IDataSet, 'Minv max')
      ranges(1) = 1.d0
      ranges(2) = sqrtS ! default to full range
      if ( idxMinv1 .ne. 0 ) 
     $   ranges(1) = DATASETInfo(idxMinv1, IDataSet)
      if ( idxMinv2 .ne. 0 )
     $   ranges(2) = DATASETInfo(idxMinv2, IDataSet)
      
      ! rap
      ! hardcoded so far
      ranges(3) = -10.d0
      ranges(4) =  10.d0

      ! eta
      ! hardcoded so far ( will not work with narrow eta)
      !idxEta1 = GetInfoIndex(IDataSet, 'Eta_el min')
      !idxEta2 = GetInfoIndex(IDataSet, 'Eta_el max')
      ranges(5) = -10.d0
      ranges(6) =  10.d0

C Get indicies:
      idxY1 = GetBinIndex(IDataSet,'y1')
      idxY2 = GetBinIndex(IDataSet,'y2')

      if (idxY1.eq.0 .or. idxY2.eq.0) then
         print 
     $'(''ERROR in GetDYNCXsection, can not find bin index for y1, y2'',2i6)'
     $        ,idxY1,idxY2
         stop
      endif

C Define bins:
      do i=1,NDATAPOINTS(IDataSet) 
         idx =  DATASETIDX(IDataSet,i)
         yb(i)   =  AbstractBins(idxY1,idx)
         yb(i+1) =  AbstractBins(idxY2,idx)
      enddo

      !print *, bchpr, sqrts, ranges
      print *,'Initialise DY calculations for dataset', IDataSet
      call dy_create_calc(IDataSet, bchpr, sqrtS, 'Z'//char(0), ranges, 
     $   'y'//char(0), NDATAPOINTS(IDataSet), yb)
      
      end

      subroutine InitDYXsectionDataset_applgrid(IDataSet)
C------------------------------------------------------------
C
C Initialise tables for DY process for calculations using applgrid
C
C------------------------------------------------------------     
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
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
      include 'ntot.inc'
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'indata.inc'
      integer IDataSet

      integer GetBinIndex                                                                                                                                    
      integer GetInfoIndex                                                                                                                                    
      double precision ranges(7)

      integer NPmax
      parameter(NPmax=100)
      double precision sqrtS, pte, ptnu
      double precision eb(Npmax+1)
      integer bchpr

      integer idx, idxEta1, idxEta2,i 
      integer idxSqrtS, idxBchpr
      integer idxPte, idxPtnu, idxMinv1, idxMinv2


C----------------------------------------------------------


      if (NDATAPOINTS(IDataSet).gt.NPmax) then
         print *,'ERROR IN InitDYCCXsectionDataset'
         print *,'INCREASE NPMax to ',NDATAPOINTS(IDataSet)                                                                                                  
         stop
      endif

C Set global parameter:
      LFitDY = .true.

      ! beam 
      idxSqrtS = GetInfoIndex(IDataSet, 'sqrt(S)')
      sqrtS = 7000d0 ! defaults to LHC
      if ( idxSqrtS .ne. 0 ) sqrtS = DATASETInfo(idxSqrtS, IDataSet)
      ! beam charge product
      bchpr = 1 ! defaults to LHC
      if ( DATASETREACTION(IDataSet) .eq. 'CC ppbar' ) 
     $ bchpr = -1

      ! pt
      idxPte = GetInfoIndex(IDataSet,'pte cut')
      idxPtnu = GetInfoIndex(IDataSet, 'ptnu cut')
      pte = 0.d0
      ptnu = 0.d0
      if ( idxPte .ne. 0 ) 
     $   pte = DATASETInfo(idxPte, IDataSet)
      if ( idxPtnu .ne. 0 )  
     $   ptnu = DATASETInfo(idxPtnu, IDataSet)
      if ( ptnu .le. pte ) then 
         ranges(7) = pte
      else  
         ranges(7) = ptnu
      endif

      ! mass
      idxMinv1 = GetInfoIndex(IDataSet, 'Minv min')
      idxMinv2 = GetInfoIndex(IDataSet, 'Minv max')
      ranges(1) = 1.d0
      ranges(2) = sqrtS ! default to full range
      if ( idxMinv1 .ne. 0 )
     $   ranges(1) = DATASETInfo(idxMinv1, IDataSet)
      if ( idxMinv2 .ne. 0 )
     $   ranges(2) = DATASETInfo(idxMinv2, IDataSet)
      
      ! rap, hardcoded
      ranges(3) = -10.d0
      ranges(4) =  10.d0
      ! eta, hardcoded
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

      !print *, bchpr, sqrts, ranges
      print *,'Initialise DY calculations for dataset', IDataSet
      call dy_create_calc(IDataSet, bchpr, sqrtS, 'W'//char(0), ranges, 
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
      include 'ntot.inc'
C      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
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
      include 'ntot.inc'
c      include 'steering.inc'
      include 'datasets.inc'
      call fastnloinit(DATASETLABEL(IDataSet),IDataSet
     $     ,DATASETTheoryFile(IDataSet)(1:Index(DATASETTheoryFile(IDataSet),' ')-1)//char(0));
      end

      subroutine InitHathorDataSet(IDataSet)
C------------------------------------------------------------
C
C Initialize Hathor reader
C
C------------------------------------------------------------
      implicit none
      integer IDataSet
      include 'ntot.inc'
      include 'datasets.inc'
      include 'steering.inc'
      integer GetInfoIndex
      integer idxSqrtS, idxReaction, idxPrecisionLevel
      double precision sqrtS
      double precision reaction
      double precision mtop
      integer precisionLevel
      logical ppbar

      idxSqrtS = GetInfoIndex(IDataSet, 'sqrt(S)')
      sqrtS = 7000d0 ! defaults to LHC
      if ( idxSqrtS .ne. 0 ) sqrtS = DATASETInfo(idxSqrtS, IDataSet)

      idxReaction = GetInfoIndex(IDataSet, 'ppbar')
      ppbar = .FALSE. ! defaults to LHC
      if ( idxReaction .ne. 0 ) then
         reaction = DATASETInfo(idxReaction, IDataSet)
         if ( reaction .eq. 1 ) ppbar = .TRUE.
      endif

      idxPrecisionLevel = GetInfoIndex(IDataSet, 'precisionLevel')
      precisionLevel = 2 ! defaults to Hathor::MEDIUM
      if ( idxPrecisionLevel .ne. 0 ) precisionLevel = DATASETInfo(idxPrecisionLevel, IDataSet)

      mtop = HF_MASS(3)

      call hathorinit(sqrtS, ppbar, mtop, I_FIT_ORDER, precisionLevel)
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

      subroutine InitHathor
      end

      Subroutine Init_EW_parameters
C-----------------------------------------------------
C
C  Initialise electroweak parameters. Created 5 June 2011
C
C  Comments:
C    12/08/2011: Move reading of the namelist to read_steer
C
C-----------------------------------------------------
      implicit none
      include 'couplings.inc'
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
      Subroutine LHAPDFsubr(x, qmu2, xf)
C-------------------------------------------------------
C
C External PDF reading for QCDNUM
C
C--------------------------------------------------------
      implicit double precision (a-h,o-z)

      double precision x,qmu2
      dimension xf(-6:6)
      call evolvePDF(x, sqrt(qmu2), xf)
      end

