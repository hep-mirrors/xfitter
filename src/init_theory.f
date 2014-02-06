
      subroutine init_theory_modules
*     ------------------------------------------------

      implicit none 
      include 'ntot.inc'
      include 'steering.inc'


*     ------------------------------------------------
*     Initialise EW parameters
*     ------------------------------------------------

      call Init_EW_parameters

*     ------------------------------------------------
*     Initialise qcdnum
*     ------------------------------------------------

      if(itheory.eq.0) then
C Init evolution code:
         call qcdnum_ini

         call Init_heavy_flavours

         if (ewfit.gt.0) call eprc_init(.true.)
      elseif(itheory.ge.100) then       
c          here goes a call to non-DGLAP 
cc           write(6,*) ' in ini_theory for itheory =',itheory
      endif

*     ------------------------------------------------
*     Initialise calculations for each dataset:
*     ------------------------------------------------
      if(Itheory.ge.100) then
ccc        write(6,*) ' ini_theory: no data sets initialised for theory ',itheory
      else
      call Init_theory_datasets
      endif

      return
      end



      Subroutine Init_heavy_flavours()
*-----------------------------------------------------
*
*   Init heavy flavour modules
*
*----------------------- ------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
C----------------------------------------


c Extra init for QCDNUM schemes:
      Call ZMVFNS_init()

C Other schemes:
      if (mod(HFSCHEME,10).eq.3) then
         Call FF_init()
      elseif ( mod(HFSCHEME,10).eq.2) then
         Call RT_Init()
      elseif ( mod(HFSCHEME,10).eq.4) then
         Call ABKM_init()
      endif
C---------------------------------
      end


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

         !> QCDNUM grid definitions
c     xmin_grid defined in steering.inc
      integer  iwt_xgrid(NMXGRID)         !> X-grid population

      integer NQGRID                     !> Min. number of Q2 sub-ranges
      PARAMETER (NQGRID=2)

      integer NQGridHF  !> Increase number of Q2 grid points to ensure that HF thresholds and the starting
                        !> Q2 scale are in.
      parameter (NQGridHF=NQGrid+4)  


c set-up of the constants
      integer iosp,nqout
 
      integer NQall     !> Actual number of grid points
      double precision Q2Grid(NQGridHF)
      double precision WQGrid(NQGridHF)

      double precision  QARR(NQGRID),WGT_Q2(NQGRID)
      data iosp/2/                   !x grid, lin/qua/spli

      double precision  tmp

      double precision as0,r20
      data as0/0.364/, r20/2.D0/!, nfin/0/ !alphas, NNLO, VFNS
      integer I,ndum,ierr,j

      double precision a,b,qt
      
      integer id1,id2
      integer nw,nwords,nx
      integer iqb,iqc,iqt


      integer NQ2bins, NXbins !> requested number of x,q2 bins
      logical ReportXGrid, ReadXGrid

      namelist/qcdnum/xmin_grid, iwt_xgrid, iosp, wgt_q2, QARR,
     $     NQ2bins, NXbins, Read_QCDNUM_Tables, ICheck_QCDNUM,
     $     ReportXGrid, ReadXGrid

      double precision xgrid(NXGridMax)
      integer NXGridAct

C Functions:
      integer iqfrmq

      common/ext_xgrid/xgrid,NXGridAct,ReadXGrid
C---------------------------------------------------------------------------------------

C-----  DEFAULTS -----------

C More detailed at high x:
      xmin_grid(1) = 9.9D-7
      xmin_grid(2) = 0.01D0
      xmin_grid(3) = 0.10D0
      xmin_grid(4) = 0.4D0
      xmin_grid(5) = 0.7D0

C Increasingly more dense grid:
      iwt_xgrid(1) = 1
      iwt_xgrid(2) = 2
      iwt_xgrid(3) = 4
      iwt_xgrid(4) = 8
      iwt_xgrid(5) = 16
C Q2 grid weights 
      WGT_q2(1) = 1.d0
      WGT_q2(2) = 1.d0
C Basic Q2 grid:
      QARR(1) = 1.
      QARR(2) = 2.05D8 ! needed for lhapdf grid  
c      QARR(2) =  64000000.      ! enough for 8 TeV LHC.

C Default sizes
      NQ2bins = 120
      NXbins  = 200

C Default QCDNUM check
      ICheck_QCDNUM = 0

C Default extra info:
      ReportXGrid = .false.
      ReadXGrid= .false.

C---------------------
      q0 = starting_scale
      qc = HF_MASS(1)**2
      qb = HF_MASS(2)**2
      qt = HF_MASS(3)**2


C----- Read grid definitions from the steering.txt
      Read_QCDNUM_Tables = .false.

      open (51,file='steering.txt',status='old')
      read (51,NML=QCDNUM,ERR=7117,END=1771)

 1771 continue

      close (51)

C Check that scales are in proper order:
      if (q0.gt.qc) then
         print *,'Starting scale must be below charm threshold, stop'
         call HF_stop
      endif
      
      if (qc.gt.qb) then
         print *,'Charm mass must be below bottom mass, stop'
         call HF_stop
      endif
      
      if (qb.gt.qt) then
         print *,'Bottom mass must be below top mass, stop'
      endif

      do i=1,NQGrid 
         Q2Grid(i) = QARR(i)
         WQGrid(i) = WGT_Q2(i)
      enddo
C Add extra points:
      Q2Grid(NQGrid+1) = q0
      Q2Grid(NQGrid+2) = qc
      Q2Grid(NQGrid+3) = qb
      Q2Grid(NQGrid+4) = qt
      WQGrid(NQGrid+1) = 4.D0
      WQGrid(NQGrid+2) = 2.D0
      WQGrid(NQGrid+3) = 1.5D0
      WQGrid(NQGrid+4) = 1.D0

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
      print '('' Q2 values at:'',20F11.1)',(Q2grid(i),i=1,NQALL)
      print '('' Weights are :'',20F11.1)',(WQgrid(i),i=1,NQALL)
      print *,' '
      print '('' Init X  grid with number of nodes='',i5)',NMXGRID
      print '('' X  values at:'',20E11.2)',(xmin_grid(i),i=1,NMXGRID),1.0
      print '('' Weights are :'',20I11)',(iwt_xgrid(i),i=1,NMXGRID)
     $     ,iwt_xgrid(NMXGRID)
      print *,' '

      call qcinit(6,' ')        !initialize
      call setord(I_FIT_ORDER)         !LO, NLO, NNLO

      call gxmake(xmin_grid,iwt_xgrid,nmxgrid,NXBINS,nx,iosp)    !x-grid
      call gqmake(Q2Grid,WQGrid,NQAll,NQ2bins,nqout)             !mu2-grid

      print '(''Requested, actual number of x bins are: '',2i5)', 
     $     nXbins,nx
      print '(''Requested, actual number of Q2 bins are: '',2i5)', 
     $     nQ2bins,nqout

      if (ReportXGrid) then
         open (51, file='xgrid.nml',status='unknown')
         call GXCOPY(xgrid, NXGridMax, NXGridAct)
         write (51,'(''&XGrid'')')
         write (51,'(''  NXgrid = '',I5)') NXGridAct
         write (51,'(''  Q20 = '',F10.2)') starting_scale
         write (51,'(''&End'')')
         write (51,'(10E26.18)') (xgrid(j),j=1,NXGridAct)
         close (51)
C         stop
      endif

      if(ReadXGrid) then 
        call ReadXGridNML()
      endif

      iqc =iqfrmq(qc)  !> Charm
      iqb =iqfrmq(qb)  !> Bottom
      iqt =iqfrmq(qt)  !> Top
      if ((mod(HFSCHEME,10).eq.3).or.HFSCHEME.eq.4.or.HFSCHEME.eq.444) then
         call setcbt(3,iqc,iqb,iqt) !thesholds in the ffns
         print *,'Fixed Flavour Number Scheme set with nf=3'
      else
        call setcbt(0,iqc,iqb,iqt) !thesholds in the ffns
      endif

      ierr = 1
      if (Read_QCDNUM_Tables) then
         call readwt(22,'unpolarised.wgt',id1,id2,nw,ierr)
      endif
      if(ierr.ne.0) then
         call fillwt(0,id1,id2,nw) !calculate weights
         call dmpwgt(1,22,'unpolarised.wgt')
      else 
         print*,' Read unpolarised weight file'
      endif
      write(6,'(/'' weight: words used ='',I10)') nw     


! setting of the evolution parameters


c      call SETABR(1.D0,0.D0)  ! mur scale variation
c      call ZMDEFQ2(1.D0,0.D0) ! muf scale variation

      print*,'exit qcdnum_Ini'
      return
 7117 continue
      print '(''Error reading QCDNUM namelist. Stop'')'
      call hf_stop
      end

      subroutine RT_init()
C-------------------------------------------------------------
C
C Created 31/12/12. Split RT code initialisation from QCDNUM
C
C-------------------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'couplings.inc'

C RT parameters
      INTEGER alphaSorderin,alphaSnfmaxin
      DOUBLE PRECISION mCharmin,mBottomin,alphaSQ0in,alphaSMZin
      DOUBLE PRECISION varin(4)
      INTEGER iordin, i
!distancein,tolerancein

      double precision qs0
C Functions:
      double precision hf_get_alphas

      

C------------------------------------
      qs0=1.d0
      alphaSQ0in = hf_get_alphas(qs0)
      mCharmin = mch
      mBottomin  = mbt
      alphaSMZin = hf_get_alphas(mz*mz)


      alphaSorderin =0.d0
      alphaSnfmaxin =3

      iordIn = I_FIT_ORDER-1
            print*,' ---------------------------------------------'
            print*,'Info from RT_init:'


      if ((HFSCHEME.eq.2).or.(HFSCHEME.eq.22)) then
         do i=1,4
            varin(i)=0d0
         enddo

         if (I_FIT_order.eq.1) then
            print*,'YOU SELECTED RT_STD LO in the steering'
         elseif (I_FIT_order.eq.2) then
            print*,'YOU SELECTED RT_STD NLO in the steering'
         elseif (I_FIT_order.eq.3) then
            print*,'YOU SELECTED RT_STD NNLO in the steering'
         endif

      elseif ((HFSCHEME.eq.202).or.(HFSCHEME.eq.222)) then
         
         varin(1)=0d0
         varin(2)=1d0
         varin(3)=-2d0/3d0
         varin(4)=1d0

         if (I_FIT_order.eq.1) then
            print*,'YOU SELECTED RT_OPT LO in the steering'
         elseif (I_FIT_order.eq.2) then
            print*,'YOU SELECTED RT_OPT NLO in the steering'
         elseif (I_FIT_order.eq.3) then
            print*,'YOU SELECTED RT_OPT NNLO in the steering'
         endif


      endif

      print*,' ---------------------------------------------'
C-
         call RT_Set_Input(varin,
     $     mCharmin,mBottomin,alphaSQ0in,alphaSMZin,
     $     alphaSorderin,alphaSnfmaxin,iordin)

cv!     $     distancein,tolerancein,

         call WATE96

C-----------------------------------------------
      end

      subroutine ZMVFNS_init()
C-------------------------------------------------------------
C
C Created 31/12/12. Split FF code initialisation from QCDNUM
C
C-------------------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'couplings.inc'
      integer nwords, ierr
C-------------------------------------------------------------
      ierr = 1
      if (Read_QCDNUM_Tables) then
         call zmreadw(22,'zmstf.wgt',nwords,ierr)
      endif
      if(MASSH.eq.1) then
      hqscale2inmass=bq2*mch*mch
      elseif(MASSH.eq.2) then
      hqscale2inmass=bq2*mbt*mbt
      endif
      if(mod(HFSCHEME,10).eq.0.and.MASSH.eq.1) then
      print*,'factorisation scale for heavy quarks is set to  sqrt(', hqscale1in,'*Q^2 + ',hqscale2in , '* 4m_c^2 )'   
      elseif(mod(HFSCHEME,10).eq.0.and.MASSH.eq.2) then
      print*,'factorisation scale for heavy quarks is set to  sqrt(', hqscale1in,'*Q^2 + ',hqscale2in , '* 4m_b^2 )'   
      endif
c    
      if(ierr.ne.0) then
         call zmfillw(nwords)
      if ( mod(HFSCHEME,10).eq.0) then
         call ZMDEFQ2(aq2,hqscale2inmass) ! muf scale variation     
      endif 
        call zmdumpw(22,'zmstf.wgt')
      else 
         print*,'Read zmstf weight file'
      endif      
      write(6,'(/'' ZMSTF: words used ='',I10)') nwords      
      call zswitch(IPDFSET)

C-------------------------------------------------------------
      end

      subroutine FF_init()
C-------------------------------------------------------------
C
C Created 31/12/12. Split FF code initialisation from QCDNUM
C
C-------------------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'couplings.inc'

      integer nwords,ierr
      double precision hqmass(3)

C-------------------------------------------------------------

      if(I_FIT_ORDER.gt.2) then
         print *,'FFN scheme can be used only with NLO, stop'
         call HF_stop
      endif

      hqmass(1) = HF_MASS(1)
      hqmass(2) = HF_MASS(2)
      hqmass(3) = HF_MASS(3)

C--- scalea1 and scaleb1 sets the heavy quarks factorisation scale

      if(MASSH.eq.1) then
         bq2  = bq2 * hqmass(1)**2
      print*,'factorisation scale for heavy quarks is set to  sqrt(', hqscale1in,'*Q^2 + ',hqscale2in , '* 4m_c^2 )'   
      elseif(MASSH.eq.2) then
         bq2  = bq2 * hqmass(2)**2
      print*,'factorisation scale for heavy quarks is set to  sqrt(', hqscale1in,'*Q^2 + ',hqscale2in , '* 4m_b^2 )'   
      endif
c         print*,'1 HQ scale (Q^2=a*mu_F^2 + b) a,b,mh', aq2,bq2,massh 

      ierr = 1
      if (Read_QCDNUM_Tables) then
         call hqreadw(22,'hqstf.wgt',nwords,ierr)
      endif
      if(ierr.ne.0) then
         call hqfillw(3,hqmass,aq2,bq2,nwords)
         call hqdumpw(22,'hqstf.wgt')
      else 
         print*,'Read hqstf.wgt file'
      endif      
      write(6,'(/'' HQSTF: words used ='',I10)') nwords      
      call hswitch(IPDFSET)

C-------------------------------------------------------------
      end


      subroutine ABKM_init()
C-------------------------------------------------------------
C
C Created 31/12/12. Split ABKM code initialisation from QCDNUM
C
C-------------------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'couplings.inc'

c ABKM parameters:
      double precision rmass8in,rmass10in
      integer kschemepdfin,kordpdfin
      logical msbarmin

C-------------------------------------------------------------
      call initgridconst

!  Take the 3-flavour scheme as a default
      kschemepdfin=0
c  c and b - quark masses         
      rmass8in=HF_MASS(1)
      rmass10in=HF_MASS(2)

! the pole mass definition by default =false (for running mass def in msbar: msbarmin=.true.)
      if(HFSCHEME.eq.444) then
        msbarmin=.true.
      else
        msbarmin=.false.
      endif
      print*,'---------------------------------------------'
      print*,'INFO from ABKM_init:'
      print*,'FF ABM running mass def? T(rue), (F)alse:', msbarmin
      print*,'---------------------------------------------'
      print*,'factorisation scale for heavy quarks  is set to sqrt(', hqscale1in,'*Q^2 + ',hqscale2in , '* 4m_q^2 )'   

c NLO or NNLO: kordpdfin=1 NLO, kordpdfin=2 NNLO
c this flag will set kordhq,kordalps,kordf2,kordfl,kordfl so same order!         
      kordpdfin  = I_FIT_ORDER-1

c set scale for FFNS only         
c      if(HFSCHEME.eq.4.or.HFSCHEME.eq.444) then
!  Set the factorization scale as sqrt(Q2*hqscale1 + 4m^2*hqscale2) for the 
!  pair heavy-quark DIS production and as sqrt(Q2*hqscale1 + m^2*hqscale2) 
!  for the single heavy-quark DIS production
c         hqscale1in=1d0
c         hqscale2in=1d0
! NEEDS TO BE IMPROVED            
c here VFNS (BMSN)           
c      else    
c         hqscale1in=1d0
c         hqscale2in=0d0
c      endif  

! ren.scale=fac.scale as a default
cc        rscale=1d0


      call ABKM_Set_Input(
     $     kschemepdfin,kordpdfin,rmass8in,rmass10in,msbarmin,
     $     hqscale1in,hqscale2in)
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
      include 'theorexpr.inc'
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
         elseif (DATASETREACTION(IDataSet).eq.'FastNLO jets' .or.   
     $           DATASETREACTION(IDataSet).eq.'FastNLO ep jets') then ! for backward compatibility
            Call InitJetsFastNLODataSet(IDataSet)
         elseif (DATASETREACTION(IDataSet)
     $           .eq.'FastNLO ep jets normalised') then
            Call InitIntegratedNCXsectionDataset(IDataSet)
            Call InitJetsFastNLODataSet(IDataSet)
         elseif (DATASETREACTION(IDataSet).eq.'ttbar') then
            Call InitHathorDataSet(IDataSet)
         elseif (DATASETREACTION(IDataSet).eq.'DDIS') then
            Call InitDDisDataSet(IDataSet)            
         elseif (DATASETREACTION(IDataSet).eq.'Dummy') then
            Call InitDummy(IDataSet)            
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

      subroutine InitDummy(IDataSet)
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

!      if (DATASETTheoryType(IDataSet).eq.'kfactor') then
!cv         call InitDummy_kfactor(IDataSet)
!      elseif (DATASETTheoryType(IDataSet).eq.'applgrid') then
!         call InitDummy_applgrid(IDataSet)         
!      else
!         print *,'InitDummy: unknown theory type'
!     $        ,DATASETTheoryType(IDataSet), ' for set ', IDataSet
!         call HF_stop
!      endif
!
!      end

C-------------------------------------------------------------

      end

      subroutine InitIntegratedNCXsectionDataset(IDataSet)
      call eprc_init(.true.)    ! needed for running of alphaem
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
      elseif (DATASETTheoryType(IDataSet).eq.'expression') then
        continue
      else
         print *,'InitDYCCXsectionDataset: unknown theory type'
     $        ,DATASETTheoryType(IDataSet), ' for set ', IDataSet
         call HF_stop
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


      print *, DATASETTheoryType(IDataSet)
      if (DATASETTheoryType(IDataSet).eq.'kfactor') then
         call InitDYNCXsectionDataset_kfactor(IDataSet)
      elseif (DATASETTheoryType(IDataSet).eq.'applgrid') then
         call InitDYXsectionDataset_applgrid(IDataSet)         
      elseif (DATASETTheoryType(IDataSet).eq.'expression') then
        continue
      else
         print *,'InitDYNCXsectionDataset: unknown theory type'
     $        ,DATASETTheoryType(IDataSet), ' for set ', IDataSet
         call HF_stop
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
         call HF_stop
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
         call HF_stop
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
      integer i

C---------------------------------------------------------------
C      call appl_readgrid(IGridID,DATASETTheoryFile(IDataSet))
C Store index:
C      DATASETTheoryIndex(IDataSet) = IGridID


      if (DATASETNapplgrid(IDataSet).gt.1) then
         do i=1,DATASETNapplgrid(IDataSet)
            call appl_readgrid(IGridID,DATASETapplgridNames(i,IDataSet))
            DATASETTheoryIndex(IDataSet) = IGridID
            DATASETApplgridTheoryIndex(i,IDataSet) = IGridID
            print*,'display the Index of Applgrid: ', i, '----',
     $           DATASETApplgridTheoryIndex(i,IDataSet)
         enddo
      endif
      
      if (DATASETNapplgrid(IDataSet).eq.1) then
         call appl_readgrid(IGridID,DATASETTheoryFile(IDataSet))
C     Store index:
         DATASETTheoryIndex(IDataSet) = IGridID
      endif



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
         call HF_stop
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
     $        '(''ERROR in GetDYCCXsection, can not find bin
     $        index for Eta1, Eta2'',2i6)'
     $        ,idxEta1,idxEta2
         call HF_stop
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
            call HF_stop
         endif

C         print *,'ho',n,n2
C Store first index:
         DATASETTheoryIndex(IDataSet) = IGridIDfnlo(1)
         IGridID = IGridIDfnlo(1)

      elseif (DATASETTheoryType(IDataSet).eq.'expression') then
        return
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
      include 'datasets.inc'

      integer GetInfoIndex

      logical PubUnits      
      double precision RealPubUnits, MurDef, MufDef, MurScale, MufScale
      double precision hf_get_mur, hf_get_muf


      RealPubUnits=(DATASETInfo(GetInfoIndex(IDataSet,
     $     'PublicationUnits'),IDataSet))
      if(RealPubUnits .eq. 1.) then
         PubUnits = .True.
      else                  
         PubUnits = .False. 
      endif
      
      MurDef=(DATASETInfo(GetInfoIndex(IDataSet,
     $     'MurDef'),IDataSet))
      MufDef=(DATASETInfo(GetInfoIndex(IDataSet,
     $     'MufDef'),IDataSet))
      MurScale=hf_get_mur(IDataSet)
      MufScale=hf_get_muf(IDataSet)

      call fastnloinit(DATASETLABEL(IDataSet),IDataSet
     $     ,DATASETTheoryFile(IDataSet)(1:Index(DATASETTheoryFile(IDataSet),' ')-1)//char(0)
     $     ,PubUnits, MurDef, MurScale, MufDef, MufScale);
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

      call hathorinit(IDataSet, sqrtS, ppbar, mtop, I_FIT_ORDER, precisionLevel)
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
c      print *, n, ' applgrid grids have been read'
C-------------------------------------
c      call HF_stop
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


c ------------------------------------------------------
c
c Diffraction
c
c-------------------------------------------------------
      Subroutine  InitDDisDataSet(IDataSet)

      implicit none
      integer IDataSet
      include 'ntot.inc'
c      include 'steering.inc'
      include 'datasets.inc'
      double precision sqrtS
      integer idxSqrtS
      integer GetInfoIndex
      
      print *,'Initialising DDIS for dataset', IDataSet
      call flush(6)
      call ddisinit(IDataSet)
      idxSqrtS = GetInfoIndex(IDataSet, 'sqrt(S)')
      sqrtS = 318d0 ! defaults to HERA2
      if ( idxSqrtS .ne. 0 ) sqrtS = DATASETInfo(idxSqrtS, IDataSet)
      call SetECMsq(IDataSet, sqrtS**2)
      
      end

      Subroutine ReadXGridNML
        include 'ntot.inc'
        double precision grid(NXGridMax)
        double precision Q20
        integer NXgrid
        logical ReadXGrid

        common/ext_xgrid/grid,nxgrid,ReadXGrid
        namelist/XGrid/NXgrid,Q20

        open (51,file='xgrid.nml',status='old')
        read (51,NML=XGrid,ERR=8118,END=1881)

 1881   continue
        
        read (51,'(10E26.18)') (grid(idx),idx=1,NXgrid);
        close (51)
        return 
 8118   continue
        print '(''Error reading XGrid namelist. Stop'')'
        call hf_stop
 8128   continue
        print '(''Error reading xgrid. Stop'')'
        call hf_stop
      end
