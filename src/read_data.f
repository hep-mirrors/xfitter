C----------------------------------------------------
!> Read data from data table
C----------------------------------------------------

      subroutine read_data

      implicit none
*     ------------------------------------------------
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "systematics.inc"
#include "indata.inc"
#include "couplings.inc"
#include "for_debug.inc"

      character*10  cdummy
      character   adum 
      character*1   a1dum 

      double precision ECM     
      double precision fac

      double precision UncorNew(NTot),UncorConstNew(NTot),
     $     StatNew(NTot), StatConstNew(NTot),UncorPoissonNew(Ntot)

      logical FIRST             !  true : cov matrix recalculated
      logical GNORM             !  correlated part for the luminosity errors
      

*     ------------------------------------------------
*     study of PDF uncertainties using MC method
*     ------------------------------------------------
     
      double precision alnorm
      external alnorm
      real logshift
      external logshift
      real lsig, lmu, lrunif
      real dummy, dummy_st
      integer vi,icount
      real ranmflat
      double precision rand_shift(NSYSMAX)
      double precision r_sh_fl(NSYSMAX)
      real rndsh, ranflat
      integer num,iseedrand, idate,is,ntime, ndate

      COMMON/SLATE/IS(40)
      real alumlognorm(300)
      real alumierr(300)
      data alumierr/300*0.0/
*     ------------------------------------------------end MC
 

      integer i,j,k,iset,n0,isys,iq2bin,iebin,jsys
      integer NSysSave

*     ------------------------------------------------
*     initilialising
*     ------------------------------------------------

C      NSYS = 0
      DEBUG  = lDEBUG
      GNORM = .true.



      do i=1,nsysMax
         do j=1,ntot
            BETA(i,j) = 0.d0
            betaasym(i,1,j) = 0d0
            betaasym(i,2,j) = 0d0
         enddo
         n_syst_meas(i) = 0  ! Also zero reference table
      enddo


      do j=1,nset
         UseFixedTheory(j) = .false.
      enddo

      do i=1,nset
         NDATAPOINTS(i) = 0
      enddo
      npoints = 0

*     --------------------------------------------end initialising


*     ------------------------------------------------
*     Read data from namelists:
*     ------------------------------------------------

      do i=1,NInputFiles
         call ReadDataFile(InputFileNames(i))
      enddo

      ! Check and read fixed theory predictions (if present)
      call read_theoryfilesNML 

C-----------------------------------------
      print*,'number of points', npoints

      do i=1,nsys
         do k=1,npoints
            beta(i,k) = beta(i,k) / 100.
            betaasym(i,1,k) = betaasym(i,1,k) / 100.
            betaasym(i,2,k) = betaasym(i,2,k) / 100.
C Get omega (quadratic term coefficient):
            omega(i,k) = (betaasym(i,1,k) + betaasym(i,2,k))/2.0
            if (beta(i,k).ne.0 .and. debug) then
               print '(3E14.4,'' omega, gamma:'')'
     $              ,omega(i,k)*100,beta(i,k)*100, 
     $              omega(i,k)/beta(i,k)
            endif
         enddo
      enddo


*     ----------------------------------------------------------------------

* prepare correlations
      call prep_corr


* Save original number of syst. errors:

      NSysSave = NSys
      call covar_to_nui(UncorNew,UncorConstNew,
     $     StatNew,StatConstNew,UncorPoissonNew) ! covariance to nuicance parameters, if needed.


      call reduce_nui(UncorNew,UncorConstNew
     $     ,UncorPoissonNew) ! We can also reduce number of nuisance parameters 


      if (LDebug) then
C
C Dump beta matrix
C
         open (61,file='beta.dat',status='unknown')
         do k=1,npoints
            write (61,'(I5,500(F6.2))') k,(Beta(i,k)*100.0,i=1,NSYSMAX)
         enddo
         close(61)
      endif

C
C Check if MC method is requested and some of the uncertainties given
C using covariance matrix. In this case LConvertCotToNui is set to true
C and warrning message is issued.
C
      if (LRand) then
         if (STATYPE.ne.0.or.SYSTYPE.ne.0) then
            if (NSys .ne. NSysSave) then
               if (.not. LConvertCovToNui ) then
                  LConvertCovToNui = .true.
                  call hf_errlog(14080401,
     $              'W: READ_DATA: MC method requested for cov. info.'
     $              //' Set LConvertCovToNui to true')
               endif
            endif
         endif
      endif
C
C MC method moved to fcn.f
C

* 
      IF (LConvertCovToNui) then
         do k=1,npoints
            is_covariance(k) = .false.
C Also re-set uncorrelated errors:
            e_uncor_mult(k)   = UncorNew(k)
            e_stat_poisson(k) = StatNew(k)
            e_uncor_const(k)  = UncorConstNew(k)
            e_stat_const(k)   = StatConstNew(k)
            e_uncor_poisson(k) = UncorPoissonNew(k)

         enddo
      else
         NSys = NSysSave
      endif

C Calculate alpha:
      do k=1,npoints
         alpha(k) =  sqrt(e_uncor_mult(k)**2
     $        +e_stat_poisson(k)**2
     $        +e_uncor_const(k)**2
     $        +e_stat_const(k)**2
     $        +e_uncor_poisson(k)**2)
     $        *daten(k)
      enddo


!
!  Split control/fit sample:
!
      if (ControlFitSplit) then
         call prepare_control_fit()
      endif

      return
      end


C------------------------------------------------------------------------
C
C  Created 20 May 2011 by SG.
!>  Read data set using namelist format
!> @param CFile input data file 
C------------------------------------------------------------------------
      subroutine ReadDataFile(CFile)

      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "indata.inc"
#include "systematics.inc"
#include "theorexpr.inc"
#include "scales.inc"
#include "for_debug.inc"

      character *(*) CFile
C Namelist  variables:    
      integer ndataMax,nsystMax,ncolumnMax
      parameter (ndataMax=ntot)
      parameter (nsystMax=nsysmax)

      parameter (ncolumnMax = nsystMax+NBinDimensionMax+1)

      character *80 Name
      integer  NData
      integer  NUncert
      integer  NBinDimension

      
      character *80 BinName(NBinDimensionMax)
c      double precision datainfo(ninfoMax)
c      character *80 CInfo(ninfoMax)
      character *80 Reaction

      double precision buffer(ncolumnMax)

C
C Name and type of columns:
C      
      integer   NColumn 
      character *64 ColumnName(ncolumnMax)
      character *64 ColumnType(ncolumnMax)

C Systematics:
      character *64 SystematicType(nsystMax)
      logical Percent(1:nsystMax)

C Reference table
      integer CompressIdx(nsystMax)

      integer IndexDataset
      double precision SystScales(nsystMax)
C Extra info about k-factors, applegrid file(s):
      character*1000 TheoryInfoFile(NKFactMax) !Is this used anymore?  --Ivan
      character*80  TheoryType(2)
      character*80 KFactorNames(NKFactMax)
      integer      NKFactor
C Infomation for open more than 1 applgrid
C     character*80 applgridNames(NapplgridMax)
      integer      NTheoryFiles
      logical ForceAdditive     ! force all errors to be treated as additive
C Variables for plotting
      integer PlotN
      character *64 PlotDefColumn
      double precision PlotDefValue(ncolumnMax)
      character *64 PlotDefTitle(ncolumnMax)
      character *64 PlotVarColumn

      character *256 PlotOptions(ncolumnMax)
      integer PlotDefColIdx, PreviousPlots
      double precision tempD
C Namelist definition:
      namelist/Data/Name,NData
     $     ,NInfo,datainfo,CInfo,Reaction,Percent
     $     ,SystScales, IndexDataset
     $     ,TheoryInfoFile,TheoryType,KFactorNames,NKFactor
     $     ,TermName,TermType,TermInfo, TermSource,TheorExpr
     $     ,ColumnName, ColumnType, NColumn
     $     ,NTheoryFiles, ForceAdditive 

      namelist/PlotDesc/PlotN, PlotDefColumn, PlotDefValue, 
     $     PlotVarColumn, PlotOptions
C--------------------------------------------------------------

      double precision XSections(ndataMax)
      integer          binFlags(ndataMax)
      integer          nDSbins
      double precision AllBins(NBinDimensionMax,ndataMax)
      double precision Syst(nsystmax)

      double precision Akfact(NKFactMax)

      double precision StatError   ! stat
      double precision StatErrorConst ! stat. error to be treated as constant
      double precision UncorError  ! uncorrelated systematics
      double precision UncorConstError  ! uncorrelated systematics
      double precision TotalError  ! total uncertainty

      double precision TotalErrorRead ! total error, provided by the data file

      integer idxSigma

      integer idxUnit
      double precision TheoryUnit  ! scale factor for theory to bring to data units.
      integer GetInfoIndex         ! function thet returns index of an information string.

c     select ppbar reaction for applgrid PDF convolution
      character*80 Msg
      integer idxReaction
      double precision ppbar_reaction

c     Normalise applgrid prediction to 1
      double precision theory_normalised

c     bin-by-bin dynamic scale in applgrid prediction
      double precision theory_dynscale

  
      integer i,j,iBin,iError
      logical LReadKFactor

C Temporary buffer to read the data (allows for comments starting with *)
      character *112768 CTmp

      integer SystematicsExist,iLen
      integer NAsymPlus(NSYSMAX), NAsymMinus(NSYSMAX)
      logical isPlus, isMinus
      
C Functions
      logical FailSelectionCuts
      integer GetBinIndex
      integer fileCheckSum
      integer checksum
      
C-------------------------------------------------------      

C Reset to default:
      NUncert = 0
      NData = 0
      NInfo = 0
      NBinDimension = 0
      Reaction = ' '
      Name     = ' '
      IndexDataSet = 0
      idxSigma = 0
      NKFactor = 0
      TheoryInfoFile = ' '
      LReadKFactor = .false.
      NTheoryFiles = 0

      do i = 1,2
        TheoryInfoFile(i) = ' '
        TheoryType(i) = ' '
      enddo
      NTerms = 0
      do i=1,NTermsMax
        TermName(i) = ' '
        TermType(i) = ' '
        TermInfo(i) = ' '
        TermSource(i) = ' '
      enddo
      TheorExpr = ' '

      nDSbins = 0
      do i = 1, ndataMax
        binFlags(i) = 1
      enddo

C Reset plotting variables
      PlotN=0
      PlotDefColumn='undefined'
c      double precision PlotDefValue(ncolumnMax)
      PlotDefTitle(1)='undefined'
      PlotVarColumn='undefined'

      ForceAdditive = .false.

C Reset scales to 1.0
      do i=1,nsysmax
         SystScales(i) = 1.0
         ColumnType(i) = ' '
         ColumnName(i) = ' '
      enddo


      open(51,file=CFile,status='old',err=99)

cc      if(DEBUG)then
        print *,'Reading data file ',trim(CFile)
cc      endif
      read(51,NML=Data,err=98)

      PlotN = -1
      read(51,NML=PlotDesc,end=96,err=97)
 96   continue

      if(PlotN.eq.-1) then  ! SUPPORT OLD FORMAT WITHOUT PLOTDESC NAMELIST, GET TO END OF DATA NAMELIST
         close(51)
         open(51,file=CFile,status='old',err=99)
         read(51,NML=Data,err=98)
      endif   

      
C
C Check dimensions
C
      if (NColumn.gt. Ncolumnmax) then
         print '(''Error in ReadDataFile for File='',A80)',cfile
         print '(''NColumn = '',i6,'' exeeds NColumnMax='',i6)',ncolumn
     $        ,ncolumnmax
         call HF_stop
      endif

C
C Store 
C
      NDATASETS = NDATASETS + 1

C     IndexDataSet = fileCheckSum(CFile)   ! can use file content
      if (IndexDataSet>=0 .and. .not. UseDataSetIndex) then
         IndexDataSet = checksum(CFile) ! use file name
         IndexDataSet = mod(abs(IndexDataset),10000000)
      else
         IndexDataSet = abs(IndexDataSet)
      endif
      print *,'Dataset index set to = ',IndexDataSet
      DATASETNUMBER(NDATASETS)   = 10000+NDATASETS
      DATASETLABEL(NDATASETS)    = Name
      DATASETNUMBER(NDATASETS)   = IndexDataset

C Reaction info:
      DATASETREACTION(NDATASETS) = Reaction

C Reset bit-masks for error types:
      iStatTypesBitMask(NDATASETS) = 0
      iUncorTypesBitMask(NDATASETS) = 0

C Parse ColumnType, count systematics, etc
      do i=1,NColumn
         if (ColumnType(i).eq.'Flag') then
            continue
         elseif (ColumnType(i).eq.'Bin') then
            NBinDimension = NBinDimension + 1
            BinName(NBinDimension) = ColumnName(i)
         elseif (ColumnType(i).eq.'Sigma') then
            idxSigma = i
         elseif (ColumnType(i).eq.'Error') then
            NUncert = NUncert + 1
            ! Special case: uncorrelated errors (constant or mult)
            if (index(ColumnName(i),'uncor const').gt.0
     $           .or.index(ColumnName(i),'uncor:A').gt.0) then
               SystematicType(NUncert) = 'uncor const'
               iUncorTypesBitMask(NDATASETS) = 
     $              IOR(iUncorTypesBitMask(NDATASETS), ibConst)
            elseif ((ColumnName(i).eq.'stat:A')
     $              .or.(ColumnName(i).eq.'stat const')) then
               SystematicType(NUncert) = 'stat const'
               iStatTypesBitMask(NDATASETS) = 
     $              IOR(iStatTypesBitMask(NDATASETS), ibConst)
            elseif (index(ColumnName(i),'uncor').gt.0) then
               SystematicType(NUncert) = 'uncor'
               iUncorTypesBitMask(NDATASETS) = 
     $              IOR(iUncorTypesBitMask(NDATASETS), ibLinear)
            elseif (ColumnName(i).eq.'stat') then
               SystematicType(NUncert) = 'stat'
               iStatTypesBitMask(NDATASETS) = 
     $              IOR(iStatTypesBitMask(NDATASETS), ibPoisson)
            else
               SystematicType(NUncert) = ColumnName(i)
            endif
         elseif (ColumnType(i).eq.'Dummy') then
! Ignore dummy column
         else
            print '(''Unknown Column type for dataset'',A80)',CFile
            print '(''Column='',i5,'' type='',A32)',i,ColumnType(i)
            print '(''STOP in ReadDataFile'')'
            call HF_stop
         endif
      enddo

C Binning info:
      DATASETBinningDimension(NDATASETS) = NBinDimension
C Filling with 'dummy' first four names for proper formation of fittedresults.txt
      do i=1,4
         DATASETBinNames(i,NDATASETS) = 'dummy'
      enddo
      do i=1,NBinDimension
         DATASETBinNames(i,NDATASETS) = BinName(i)
      enddo

C Extra info:
      DATASETInfoDimension(NDATASETS) = NInfo
      do i=1,NInfo
         DATASETInfoNames(i,NDATASETS) = CInfo(i)
         DATASETInfo(i,NDATASETS) =      DataInfo(i)
      enddo

      dsname = name
      ds_index = IndexDataset 

C Prepare systematics:
      do i=1,NUncert
C--- Statistical: special case
         if (SystematicType(i).eq.'stat') then
         else if (SystematicType(i).eq.'stat const') then
            Call HF_ERRLOG(16020001,'I: Stat Const error used in: '//Name)
C--- Uncorrelated: special case
         else if (SystematicType(i).eq.'uncor const') then
           Call HF_ERRLOG(14030505,'I: Uncor Const error used in: '//Name)
C--- Uncorrelated: special case
         else if (SystematicType(i).eq.'uncor') then
c            Call HF_ERRLOG(14030506,'I: Uncor Error type used')
C--- Total error: special case
         else if (SystematicType(i).eq.'total') then
            Call HF_ERRLOG(14030507,'I: Total error used in: '//Name)
C--- Ignore: special case
         else if (SystematicType(i).eq.'ignore') then

         else
C--- Check if the source already exists:         
            j = SystematicsExist(SystematicType(i))
C Not found:
            if (j.eq.0)  then
C--- Add new source
               Call AddSystematics(SystematicType(i))
               CompressIdx(i) = NSYS               
            else
               CompressIdx(i) = j
            endif
         endif
      enddo

C Count theory expression terms
      CTmp = ' '
      do i = 1,NTermsMax
        if (TermName(i) .eq. ' ' ) goto 88
        NTerms = i
        CTmp = TermName(i)
        TermName(i) = trim(CTmp)
        CTmp = TermType(i)
        TermType(i) = trim(CTmp)
        CTmp = TermInfo(i)
        TermInfo(i) = trim(CTmp)
        CTmp = TermSource(i)
        TermSource(i) = trim(CTmp)
      enddo
 88   continue 

      if (Reaction .eq. ' '
     $     .or.TermType(1).eq.'reaction'
     $     .or.TermSource(1).ne.' ') then
         if (TheoryType(1) .eq. ' ') then
            TheoryType(1) = 'expression'
         endif
      endif

      if ( TheoryType(1) .eq. 'expression') then
         do i =1, NTerms
            if (TermType(i) .eq. ' ') then
               TermType(i) = 'reaction'
            endif
         enddo
      endif

      if (TheoryType(1).ne. 'expression') then
         call hf_errlog(18030710+NDATASETS,
     $        'W: Using obsolete theory calculation for data file '
     $        //trim(CFile) )
      endif

C Theory file if present:
      DATASETTheoryType(NDATASETS) = ' '
      do i=1,2
         if (Index(TheoryInfoFile(i),' ').gt.nchar_theory-1) then
            print *,
     $        'File name too long (>199). ',trim(TheoryInfoFile(i))
            print *,'Increase DATASETTheoryFile'//
     $        ' string length in dataset.inc'
            print *,'or reduce the name (e.g. make symbolic link)'
           call hf_stop
         endif
         if (TheoryType(i).ne.' ' .and. TheoryInfoFile(i).ne.' ') then
            if (TheoryType(i).ne.'kfactor') then
  !   not k-factor, overwrite
               DATASETTheoryFile(NDATASETS) = TheoryInfoFile(i)
               DATASETTheoryType(NDATASETS) = TheoryType(i)
            else
  !   k-factor, depends if nothing else is present
               if (DATASETTheoryType(NDATASETS).eq.' ') then
                  DATASETTheoryFile(NDATASETS) = TheoryInfoFile(i)
                  DATASETTheoryType(NDATASETS) = TheoryType(i)
               endif
               open (53,file=TheoryInfoFile(i),status='old',err=100)
               lreadkfactor = .true.
            endif
         endif
      enddo
     
      DATASETNKfactors(NDATASETS) = NKFactor
      do i=1,NKFactor
         DATASETKFactorNames(i,NDATASETS) = KFactorNames(i)
      enddo
     
     
C     Count applgrids
      do i=1,2
        if(TheoryType(i).EQ.'applgrid') then 
           NTheoryFiles = NTheoryFiles+1   
        endif
      enddo
c      print*,'NTheoryFiles with allpgrids ',NTheoryFiles
     
      DATASETNapplgrid(NDATASETS) = NTheoryFiles
      do i=1,NTheoryFiles
         print*,'Theory files: ', TheoryInfoFile(i)
C     ---> copy the names in a new variable 
         DATASETapplgridNames(i,NDATASETS) = TheoryInfoFile(i)
      enddo

      ! set parameters for general theory interface here instead of
      ! src/init_theory.f.  A.S.
      if ( TheoryType(1).eq.'expression' ) then
        if ( NTerms .eq. 0 ) then
          print *,'Expression theory type selected,but no terms/expression specified'
          call hf_stop
	endif
        DATASETTheoryType(NDATASETS) = TheoryType(1)
        idxReaction = GetInfoIndex(NDATASETS,'ppbar')

        idxReaction = GetInfoIndex(NDATASETS,'Normalised')
        normalised = 0    ! defaults to absolute cross section
        if ( idxReaction .ne. 0 ) then
           normalisation = DATASETInfo(idxReaction, NDATASETS)
           if ( normalisation .ge. 0 ) normalised = 1

       write (Msg,'(''I: Normalise APPLGRID prediction dataset: '',A20,'' '')')
     $        Name
           call HF_errlog(14030401,trim(Msg))
        endif

        idxReaction = GetInfoIndex(NDATASETS,'DynamicScale')
        dynscale = 0               ! defaults to applgrid scale
        if ( idxReaction .ne. 0 ) then
           dynscale = DATASETInfo(idxReaction, NDATASETS)

      write (Msg,'(''I: Emulate dynamic scale dataset: '',A20,'' '')')
     $        Name
           call HF_errlog(14042001,trim(Msg))
        endif

        call set_theor_eval(NDATASETS)
      endif


C Read data info:
      do j=1,NData
C Allow for comments:
 89      read (51,'(A)',err=1017,end=1018) ctmp
         if (ctmp(1:1).eq.'*') then
C     Comment line, read another one
            goto 89
         endif

C Check coherence of the table info
         if (idxSigma.eq.0) then
            print *,
     $'No column contains Sigma keyword for the x-section info!!!'
            call HF_stop
         endif
         do i=1,NColumn
            if (ColumnName(i) .eq. ' ') then
               print *,'Undefined ColumnName !!!'
               print *,'Check name for column number = ',i
               call HF_stop
            endif
         enddo

C Read the colums
         read (ctmp,*,err=1019)(buffer(i),i=1,NColumn)

C Decode the columns
         iBin   = 0
         iError = 0
         do i=1,NColumn
            if (ColumnType(i).eq.'Flag') then
               binFlags(j) = nint(buffer(i))
            elseif (ColumnType(i).eq.'Bin') then
               iBin = iBin + 1
               allbins(iBin,j) = buffer(i)
            elseif (ColumnType(i).eq.'Sigma') then
               XSections(j) = buffer(i)
            elseif (ColumnType(i).eq.'Error') then
               iError = iError + 1
               syst(iError) = buffer(i)
            endif
         enddo

c         read(51,*)(allbins(i,j),i=1,NBinDimension),XSections(j)
c     $        ,(syst(i),i=1,NUncert)

C Scale the syst. erros:
         do i=1,NUncert
            Syst(i) = Syst(i) * SystScales(i)
         enddo

         if (lreadkfactor) then
            read (53,*) (akfact(i),i=1,NKFactor)
         endif

         nDSbins = nDSbins +1

C Apply cuts:
         if (FailSelectionCuts(Reaction,NBinDimension,allbins(1,j),BinName,IndexDataset)) then
	   ! set excluding flag for those bins that were cut
           binFlags(j) = 0
! Since xFitter 2.1 "Reaction" field in dataset is no longer used to select
! reaction, we have a system with terms and reaction modules now
!           if((Reaction.eq.'FastNLO jets').or.
!     $       (Reaction.eq.'FastNLO ep jets').or.
!     $       (Reaction.eq.'FastNLO ep jets normalised')) then
!              call fastnlopointskip(NDataSets, j, NData);
!           endif
           goto 1717
         endif

C skip those bins that have 0 flag
         if ( binFlags(j) .eq. 0 ) then
           goto 1717
         endif

C Add a point:
         npoints = npoints+1

C By default it is fitted:
         FitSample(Npoints) = .true.
         
         if (npoints.ge.NTOT) then
            print 
     $           '(''ReadDataFile Error: exceeding NTOT'')'
            print '(''Current NTOT='',i6)',NTOT
            print '('' Increase NTOT_C in include/dimensions.h'')'
            call HF_stop
         endif

         NDATAPOINTS(NDATASETS) = NDATAPOINTS(NDATASETS) + 1
         DATASETIDX(NDATASETS,NDATAPOINTS(NDATASETS)) = npoints

C Translate errors in %:
         TotalError = 0.
         UncorError = 0.
         UncorConstError = 0.
         StatError = 0.
         StatErrorConst = 0.
         TotalErrorRead = 0.


         do i=1,NUncert
            if (.not.Percent(i)) then
               syst(i) = syst(i)/XSections(j)*100.
            endif



            if (SystematicType(i).eq.'total') then
               TotalErrorRead = Syst(i)
            elseif (SystematicType(i).eq.'ignore') then
C Ignore error source called 'ignore'
            else
               TotalError = TotalError + Syst(i)**2
            endif

c RP handle case when only tot error given (and e.g. full covariance matrix)
c this affects only plots            
            if(TotalError.eq.0.and.TotalErrorRead.ne.0) then
               TotalError = TotalErrorRead**2 
            endif

C Uncor const:            
            if (SystematicType(i).eq.'uncor const') then
               UncorConstError = UncorConstError + Syst(i)**2
            endif

            if (SystematicType(i).eq.'uncor') then
C Uncorrelated error:
               UncorError = UncorError +  Syst(i)**2
            endif
            if (SystematicType(i).eq.'stat') then
C Stat error:
               StatError = StatError +  Syst(i)**2
            endif

            if (SystematicType(i).eq.'stat const') then
C Stat error:
               StatErrorConst = StatErrorConst +  Syst(i)**2
            endif

         enddo

         StatError = sqrt(StatError)
         StatErrorConst = sqrt(StatErrorConst)
         UncorConstError = sqrt(UncorConstError)
         UncorError = sqrt(UncorError)
         TotalError = sqrt(TotalError)

         DATEN(npoints) = XSections(j)

C  XXXXXXXXXXXXXXXXXXXXXXXXX START to become obsolete !!!
         E_UNC(npoints)  = UncorError
         E_UNC_Const(npoints) = UncorConstError
         E_TOT(npoints)  = TotalError
         E_STA(npoints)  = StatError
         E_STA_CONST(npoints) = StatErrorConst
         
C  XXXXXXXXXXXXXXXXXXXXXXXXX END to become obsolete !!!


C XXXXXXXXXXXXXXXXXXXXXXXXX
         Call SetUncorErrors(npoints, StatError,
     $        StatErrorConst,UncorError,UncorConstError)

         LForceAdditiveData(npoints) = ForceAdditive
         
         !  Check total error
         if (TotalErrorRead.ne.0) then
            if ( abs(TotalError -TotalErrorRead)/TotalErrorRead.gt.0.01) then
               print 
     $'(''WARRNING IN READDATA, LARGE DEVIATION FOR TOTAL ERROR'')'
               print '(''Total calculated='',G10.4,'' READ='',G10.4)',
     $              totalError,TotalErrorRead
            endif
         endif

         do i=1,NBinDimension
            AbstractBins(i,npoints) = allbins(i,j)
         enddo


C Reset:
         do i=1,NUncert
            if ( CompressIdx(i).gt.0 ) then
               NAsymPlus(CompressIdx(i))     =  0
               NAsymMinus(CompressIdx(i))    =  0
            endif
         enddo

         do i=1,NUncert
            if (SystematicType(i).ne.'uncor' .and. 
     $           SystematicType(i).ne.'uncor const'.and.
     $           SystematicType(i).ne.'ignore'.and.
     $           SystematicType(i).ne.'stat'.and.
     $           SystematicType(i).ne.'total'.and.
     $           SystematicType(i).ne.'stat const'
     $           ) then


               BETA(CompressIdx(i),npoints) = syst(i)
     $              *SysScaleFactor(CompressIdx(i))

               
C     Store also asymmetric errors:
               iLen   = Len_trim( SystematicType(i))
               isPlus  = SystematicType(i)(iLen:iLen).eq.'+'
               isMinus = SystematicType(i)(iLen:iLen).eq.'-'

               if (isPlus) then
                  NAsymPlus(CompressIdx(i)) = NAsymPlus(CompressIdx(i)) 
     $                 + 1

C Too many pluses and minuses !
                  if (NAsymPlus(CompressIdx(i)).gt.1) then
                     print *,' '
                     print *,'===== ERROR ERROR ERROR ===='
                     print *,' ' 
                     print *,'Problem with systematic source ',
     $                    SystematicType(i)
                     print *,
     $ 'Positive variations defined more than once'
                     print *,'Check the data file, stopping'
                     call hf_errlog(17112012,
     $                    'F: Problem with asymmetric errors')
                     call hf_stop
                  endif
C Store:
                  BetaAsym(CompressIdx(i),1,npoints) = syst(i)
     $                 *SysScaleFactor(CompressIdx(i))                
               endif

               if (isMinus) then
                  NAsymMinus(CompressIdx(i)) = NAsymMinus(CompressIdx(i)) 
     $                 + 1

C  Too many pluses and minuses !
                  if (NAsymMinus(CompressIdx(i)).gt.1) then
                     print *,' '
                     print *,'===== ERROR ERROR ERROR ===='
                     print *,' ' 
                     print *,'Problem with systematic source ',
     $                    SystematicType(i)
                     print *,
     $ 'Negative variations defined more than once'
                     print *,'Check the data file, stopping'
                     call hf_errlog(17112012,
     $                    'F: Problem with asymmetric errors')
                     call hf_stop
                  endif
C  Store:
                  BetaAsym(CompressIdx(i),2,npoints) = syst(i)
     $                 *SysScaleFactor(CompressIdx(i))                
               endif

C  Symmetrise:
               if (NAsymPlus(CompressIdx(i)).eq.1
     $              .and. NAsymMinus(CompressIdx(i)).eq.1 ) then
                  
                  BETA(CompressIdx(i),npoints) = 
     $                 0.5*( BetaAsym(CompressIdx(i),1,npoints)-
     $                        BetaAsym(CompressIdx(i),2,npoints))

                  LAsymSyst(CompressIdx(i)) = .true.
               endif

C     Correct total error shown in plots
               if (NAsymPlus(CompressIdx(i)).eq.1
     $              .and. NAsymMinus(CompressIdx(i)).eq.1 ) then
                  E_TOT(npoints)  = sqrt(E_TOT(npoints)**2
     $                 -(BetaAsym(CompressIdx(i),1,npoints)
     $                 /SysScaleFactor(CompressIdx(i)))**2
     $                 -(BetaAsym(CompressIdx(i),2,npoints)
     $                 /SysScaleFactor(CompressIdx(i)))**2
     $                 +(BETA(CompressIdx(i),npoints)
     $                 /SysScaleFactor(CompressIdx(i)))**2)
               endif

               
               if ( (NAsymPlus(CompressIdx(i)).eq.1
     $              .and. NAsymMinus(CompressIdx(i)).eq.1)
     $              .or. 
     $              ( NAsymPlus(CompressIdx(i)).eq.0
     $              .and.  NAsymMinus(CompressIdx(i)).eq.0)
     $              ) then
                
C--- Add data point to the syst. list (this will help to speedup loops):
                  n_syst_meas(CompressIdx(i)) = n_syst_meas(CompressIdx(i))
     $                 + 1
                  syst_meas_idx(n_syst_meas(CompressIdx(i)),CompressIdx(i)) 
     $                 = npoints

               endif

            endif
         enddo


         JSET(npoints) = NDATASETS ! IndexDataset  
         GPlotVarCol(NDATASETS) = PlotVarColumn
         GNPlots(NDATASETS) = PlotN
         PreviousPlots = 0

         do i=1,NDATASETS-1
c additional check to avoid boundary error (GNPlots getting negative) 
c for data which have no Drawing options defined       
            if(GNPlots(i).gt.0) then
               PreviousPlots = PreviousPlots + GNPlots(i)
            endif   
         enddo
         
         do i=1,PlotN
            GPlotOptions(PreviousPlots + i) = PlotOptions(i)
         enddo

c Find plot numbers         
         i=0
         if(PlotN.gt.0) then
            PlotDefColIdx = GetBinIndex(NDataSets,TRIM(plotdefcolumn))
            if(PlotDefColIdx.eq.0) then
               call HF_Errlog(13012801,
     $              'W:Plotting: Can not find one of the columns')
            endif

            if(PlotDefColIdx.ne.0) then
               tempD = AbstractBins(PlotDefColIdx,npoints)
               do while ((PlotDefValue(i+1).lt.tempD).AND.(i.lt.PlotN))
                  i = i+1
               enddo
               if(PlotDefValue(i+1).lt.tempD) then
                  i=0
               endif
            endif
         endif
         if (i.eq.0) then
            i=999
         endif
         JPLOT(npoints) = i
            

C Store k-factors:
         if (lreadkfactor) then
            do i=1,nkfactor
               kfactors(i,npoints) = akfact(i)
            enddo
         endif
         
 1717  continue
      enddo

C Set data binning information in theory evaluations
c but firest check that there are two columns per each bin dimension
      if ( DATASETTheoryType(NDATASETS).eq.'expression' ) then
c        if ( mod(NBinDimension,2) .ne. 0 ) then
c          print *, 'Problem reading data from ', CFile
c          print *, 'There must be two bin columns per each bin dimension'
c          print *, 'for applgrid based fits.'
C          call hf_stop
c        endif
      
        call set_theor_bins(NDATASETS, NBinDimension, nDSbins, 
     &    binFlags, allbins, binname )

        idxUnit = GetInfoIndex(NDATASETS,'theoryunit')
        if (idxUnit.gt.0) then
          Theoryunit = DATASETInfo(idxUnit,NDATASETS)
        else
          Theoryunit = 1.
        endif

!        call set_theor_units(NDATASETS, Theoryunit)
        call init_theor_eval(NDATASETS)
      endif

      close (51)
      if (lreadkfactor) then
         close (53)
      endif

c      if(DEBUG)then
        print '(''Read'',i8,'' data points for '',A80)',NData,Name
        print '(''Printing first'',i5,'' data points'')',min(Ndata,3)
        print '(20A14)',(BinName(i),i=1,NBinDimension),' sigma'
        do j=1,min(NData,3)
           print '(20E14.4)',(Allbins(i,j),i=1,NBinDimension),XSections(j)
        enddo
c      endif
      return

 97   continue
      print '(''Error reading namelist PlotDesc'')'
      print *,CFile
      call HF_stop

 98   continue
      print '(''Error reading namelist Data'')'
      print *,CFile
      call HF_stop

 99   continue
      print '(''Can not open file '')'
      print *,CFile
      call HF_stop
100   continue
      print '(''Can not open file '')'
      print *,TheoryInfoFile
      call HF_stop

 1017 continue
      print '(''Error reading file'')'
      call HF_stop
 1018 continue
      print '(''End of file while reading file'')'
      call HF_stop
 1019 continue
      print '(''Problem interpreting data line='',i6)',j
      call HF_stop

      end


C---------------------------------------------------------------------------
C
C  Created 16/07/2012. 
!>  Split data into "fit" and "control" samples.
!>  Fit sample is used in chi2 minimization, control sample is used to make
!>  sure that there is no overfitting. 
!>
C----------------------------------------------------------------------------
      subroutine prepare_control_fit()

      implicit none
#include "ntot.inc"
#include "indata.inc"
      integer i
      real ranflat
C----------------------------------------------------------------------------

      NControlPoints = 0
      NFitPoints     = 0
      if (ControlFitSplit) then
         do i=1,Npoints
            call ranlux(ranflat,1)
            FitSample(i) = (ranflat.ge.0.5)
            if (FitSample(i)) then
               NFitPoints = NFitPoints + 1
            else
               NControlPoints = NControlPoints + 1
            endif
         enddo
      else
         do i=1,Npoints
            FitSample(i) = .true.
         enddo
      endif


      end

C------------------------------------------------------------
!> Scale data uncertainties if requiered
!> @param Idx data point index
!> @param StatError, UncorError generic errors which move around
!> @param StatErrorConst 
!> @param UncorConstError
C-------------------------------------------------------------
      Subroutine SetUncorErrors(Idx,StatError,
     $     StatErrorConst,UncorError,UncorConstError)

      implicit none
      integer Idx
      double precision StatError, StatErrorConst,UncorError,UncorConstError
#include "ntot.inc"
#include "steering.inc"
#include "indata.inc"
C---------------------------------------------------------
C 
      if (StatScale.eq.'NoRescale') then
         e_stat_poisson(idx)  = 0.
         e_stat_const(idx)    = sqrt(
     $        StatErrorConst**2+StatError**2)  / 100.        
      elseif (StatScale.eq.'Poisson') then
         e_stat_poisson(idx)  = StatError       / 100.
         e_stat_const(idx)    = StatErrorConst  / 100.
      else
         print *,'ERROR !!!'
         print *,'Unknown StatScale = ',StatScale
         print *,'STOP'
         call hf_stop
      endif

      if (UncorSysScale.eq.'NoRescale') then
         e_uncor_mult(idx)    = 0.
         e_uncor_const(idx)   = sqrt(
     $        UncorError**2+UncorConstError**2) / 100.
         e_uncor_poisson(idx) = 0.
         e_uncor_logNorm(idx) = 0.
      elseif (UncorSysScale.eq.'Poisson') then
         e_uncor_mult(idx)    = 0.
         e_uncor_const(idx)   = UncorConstError / 100.
         e_uncor_poisson(idx) = UncorError / 100.
         e_uncor_logNorm(idx) = 0.      
      elseif (UncorSysScale.eq.'Linear') then
         e_uncor_mult(idx)    = UncorError / 100.
         e_uncor_const(idx)   = UncorConstError / 100.
         e_uncor_poisson(idx) = 0.
         e_uncor_logNorm(idx) = 0.      
      elseif (UncorSysScale.eq.'LogNormal') then
         e_uncor_mult(idx)    = UncorError / 100.
         e_uncor_const(idx)   = UncorConstError / 100.
         e_uncor_poisson(idx) = 0.
         e_uncor_logNorm(idx) = 0.      
         Call Hf_Errlog(1208201201,
     $   'E: Lognormal errors unsuported, use Linear scaling')
      else
         print *,'ERROR !!!'
         print *,'Unknown UncorSysScale = ',UncorSysScale
         print *,'STOP'
         call hf_stop
      endif

c      print *,idx,e_stat_poisson(idx),e_uncor_mult(idx)

C---------------------------------------------------------
      end


! Created 29/05/13. 
!> Check if fixed theory predictions are requested for datasets
      subroutine read_theoryfilesNML
C
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "systematics.inc"
      character*256 InputTheoNames(NSET)
      Namelist/InTheory/InputTheoNames
      integer i
C-----------------------------------------------------------------
      do i=1,NSET ! InputFiles
         InputTheoNames(i) = ''
      enddo      

      open (51,file='steering.txt',status='old')
      read (51,NML=InTheory,END=141,ERR=42)
      

      do i=1,NSET ! InputFiles
         if ( InputTheoNames(i) .ne. '') then
            call hf_errlog(13052901,'I:Use fixed theory predictions') 
            Call read_theory_file(InputTheoNames(i),i)
         endif
      enddo

 141  continue
      close (51)
      return
 42   call hf_errlog(1,'F:Error reading InTheory namelist')
C-----------------------------------------------------------------
      end

      !> Read fixed theory file, associate with dataset idxDataSet
      !> @param FileName input theory file
      !> @param IdxDataSet data set theory is associated with
      Subroutine read_theory_file(FileName,IdxDataSet)
      
      implicit none
      character *(*) FileName
      integer IdxDataSet
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "indata.inc"
#include "systematics.inc"
#include "theo.inc"
C
      integer nsystMax,ncolumnMax
      parameter (nsystMax=nsysmax)
      parameter (ncolumnMax = nsystMax+NBinDimensionMax+1)
      
      character*80 Name
      integer NData,NColumn
      character *64 ColumnName(ncolumnMax)
      character *64 ColumnType(ncolumnMax)
      character *64 SystematicType(nsystMax)
      logical Percent(1:nsystMax)

      integer i,j,idx
      integer ilen
      integer ipoint
      logical isPlus, isMinus
      integer iBin
      integer NBinDimension, idxSigma, NUncert
      character*4096 CTmp
      double precision buffer(ncolumnMax)
      double precision syst(nsystmax)

      integer NAsymPlus(NSYSMAX), NAsymMinus(NSYSMAX)
      integer iError

      double precision 
     $     theo_err2_up(NTOT),
     $     theo_err2_down(NTOT)

C Reference table
      integer CompressIdx(nsystMax)

      namelist/Data/Name, NData, NColumn, ColumnType, ColumnName, Percent
C Function:
      integer SystematicsExist
C---------------------------------------------------------------------
      NBinDimension = 0
      idxSigma = 0 
      NUncert = 0
      NData = 0

      UseFixedTheory(idxdataset) = .True.

      open (52,file=FileName,err=101)
      read (52,nml=Data,err=102,end=103)
C Basic consistency check:
      if (.not. pdfrotate) then
         if (NData.ne.NDATAPOINTS(IdxDataSet)) then
            print *,ndata,NDATAPOINTS(IdxDataSet),IdxDataSet
            call hf_errlog(4,
     $           'F:Mismatch for number of points in theory file '
     $           //trim(FileName))
         endif
      else
         NDataSets = max(NDataSets,idxdataset)
         NDATAPOINTS(IdxDataSet) = NData
         Datasetlabel(idxdataset) = Name
      endif

       do i=1,NColumn
         if (ColumnType(i).eq.'Bin') then
            NBinDimension = NBinDimension + 1
         elseif (ColumnType(i).eq.'Theory') then
            idxSigma = i
         elseif (ColumnType(i).eq.'Error') then
            NUncert = NUncert + 1
            SystematicType(NUncert) = ColumnName(i)
         elseif (ColumnType(i).eq.'Flag') then
            continue
         else
            call hf_errlog(5,'F:Unknown column type in file '
     $           //trim(FileName))
         endif   
      enddo  

C Some more basic checks:
      if (.not. pdfrotate) then
         if (DATASETBinningDimension(IdxDataSet).ne. NBinDimension) then
            call hf_errlog(6,
     $           'F:Binning dimension does not match in file '
     $           //trim(filename))
         endif
      else
         DATASETBinningDimension(IdxDataSet) = NBinDimension
         do i=1,NBinDimension
            DATASETBinNames(i,idxdataset) = ColumnName(i)
         enddo
      endif
      
      if (idxSigma.eq.0) then
         call hf_errlog(7,'F:Did not find theory column in file '
     $        //trim(filename))  
      endif


C Prepare theory systematics:
      do i=1,NUncert
         if (SystematicType(i).eq.'stat') then
            call hf_errlog(13052902,
     $  'I:Theory prediction includes stat. uncertainty')
         else if (SystematicType(i).eq.'stat const') then
            call hf_errlog(14030501,
     $  'I:Theory prediction includes stat const uncertainty')
C--- Uncorrelated const
         else if (SystematicType(i).eq.'uncor const') then
            call hf_errlog(14030502,
     $  'I:Theory prediction includes uncor const uncertainty')
C--- Uncorrelated
         else if (SystematicType(i).eq.'uncor') then
            call hf_errlog(14030503,
     $  'I:Theory prediction includes uncor uncertainty')
C--- Total error: special case
         else if (SystematicType(i).eq.'total') then
            call hf_errlog(14030504,
     $  'I:Theory prediction includes total uncertainty')
C--- Ignore: special case
         else if (SystematicType(i).eq.'ignore') then


         else
C--- Check if the source already exists:         
            if (index(SystematicType(i),':D').eq.0
     $           .and.index(SystematicType(i),':T').eq.0 ) then
               j = len_trim(SystematicType(i))
               if ( SystematicType(i)(j:j).eq.'+') then
                  SystematicType(i) = SystematicType(i)(1:j-1)//':T+'
               elseif ( SystematicType(i)(j:j).eq.'-') then
                  SystematicType(i) = SystematicType(i)(1:j-1)//':T-'
               else
                  SystematicType(i) = SystematicType(i)(1:j)//':T'
               endif
            endif
            j = SystematicsExist(SystematicType(i))


C Not found:
            if (j.eq.0)  then
C--- Add new source
               Call AddSystematics(SystematicType(i))
               CompressIdx(i) = NSYS
            else
               CompressIdx(i) = j
            endif
         endif
      enddo

C Read the predictions:
      ipoint = 0
      do j=1,NDATA
 89      read (52,'(A)',err=1017,end=1018) ctmp
         if (ctmp(1:1).eq.'*') then
C     Comment line, read another one
            goto 89
         endif

C Read the colums
         ipoint = ipoint + 1

         read (ctmp,*,err=1019)(buffer(i),i=1,NColumn)


C Store:
         if (pdfrotate) then
            NPoints = NPoints + 1
            idx = NPoints
            DATASETIDX(idxdataset,ipoint) = idx
         else
            idx = DATASETIDX(idxdataset,ipoint)
         endif


         iError = 0
         iBin   = 0
         do i=1,NColumn
            if (ColumnType(i).eq.'Error') then
               iError = iError + 1
               syst(iError) = buffer(i)    
               if (.not. Percent(iError)) then
                  if (daten(idx).ne.0) then
                     syst(iError) = syst(iError)/daten(idx)*100. ! BUG buffer(idxSigma)*100.
                  else
                     syst(iError) = syst(iError)/buffer(idxSigma)*100.  ! for rotate.
                  endif
               endif
            elseif (ColumnType(i).eq.'Bin') then
               iBin  = iBin + 1
               if (pdfrotate) then                  
                  AbstractBins(iBin,idx) = buffer(i)
               endif
            endif
         enddo


C Reset:
         do i=1,NUncert
            if ( CompressIdx(i).gt.0 ) then
               NAsymPlus(CompressIdx(i))     =  0
               NAsymMinus(CompressIdx(i))     =  0
            endif
         enddo


         theo_fix(idx)  = buffer(idxSigma)

         do i=1,NUncert
            if (SystematicType(i).ne.'uncor' .and. 
     $           SystematicType(i).ne.'uncor const'.and.
     $           SystematicType(i).ne.'ignore'.and.
     $           SystematicType(i).ne.'stat'.and.
     $           SystematicType(i).ne.'total'.and.
     $           SystematicType(i).ne.'stat const'
     $           ) then


               
               BETA(CompressIdx(i),idx) = -syst(i)

C     Store also asymmetric errors:
               iLen   = Len_trim( SystematicType(i))
               isPlus  = SystematicType(i)(iLen:iLen).eq.'+'
               isMinus = SystematicType(i)(iLen:iLen).eq.'-'

               if (isPlus) then
                  NAsymPlus(CompressIdx(i)) = NAsymPlus(CompressIdx(i)) 
     $                 + 1

C ! Too many pluses and minuses !
                  if (NAsymPlus(CompressIdx(i)).gt.1) then
                     print *,' '
                     print *,'===== ERROR ERROR ERROR ===='
                     print *,' ' 
                     print *,'Problem with systematic source ',
     $                    SystematicType(i)
                     print *,
     $ 'Positive variations defined more than once'
                     print *,'Check the data file, stopping'
                     call hf_errlog(17112012,
     $                    'F: Problem with asymmetric errors')
                     call hf_stop
                  endif
C Store:
                  BetaAsym(CompressIdx(i),1,idx) = - syst(i)
     $                 *SysScaleFactor(CompressIdx(i))                
               endif

               if (isMinus) then
                  NAsymMinus(CompressIdx(i)) = NAsymMinus(CompressIdx(i)) 
     $                 + 1

C ! Too many pluses and minuses !
                  if (NAsymMinus(CompressIdx(i)).gt.1) then
                     print *,' '
                     print *,'===== ERROR ERROR ERROR ===='
                     print *,' ' 
                     print *,'Problem with systematic source ',
     $                    SystematicType(i)
                     print *,
     $ 'Negative variations defined more than once'
                     print *,'Check the data file, stopping'
                     call hf_errlog(17112012,
     $                    'F: Problem with asymmetric errors')
                     call hf_stop
                  endif
C ! Store:
                  BetaAsym(CompressIdx(i),2,idx) = - syst(i)
     $                 *SysScaleFactor(CompressIdx(i))                
               endif

C ! Symmetrise:
               if (NAsymPlus(CompressIdx(i)).eq.1
     $              .and. NAsymMinus(CompressIdx(i)).eq.1 ) then
                  
                  BETA(CompressIdx(i),idx) = 
     $                 0.5*( BetaAsym(CompressIdx(i),1,idx)-
     $                        BetaAsym(CompressIdx(i),2,idx))

                  LAsymSyst(CompressIdx(i)) = .true.
               endif


               if ( (NAsymPlus(CompressIdx(i)).eq.1
     $              .and. NAsymMinus(CompressIdx(i)).eq.1)
     $              .or. 
     $              ( NAsymPlus(CompressIdx(i)).eq.0
     $              .and.  NAsymMinus(CompressIdx(i)).eq.0)
     $              ) then
                
C--- Add data point to the syst. list (this will help to speedup loops):
                  n_syst_meas(CompressIdx(i)) = n_syst_meas(CompressIdx(i))
     $                 + 1
                  syst_meas_idx(n_syst_meas(CompressIdx(i)),CompressIdx(i)) 
     $                 = idx

               endif
            elseif (SystematicType(i).eq.'stat') then
               theo_stat(idx) = Syst(i)
            elseif (SystematicType(i).eq.'uncor') then
               theo_unc(idx) =  Syst(i)
            endif
         enddo

      enddo

c      do i=1,NData
      do i=1,Npoints
         idx = DATASETIDX(idxdataset,i)
         THEO_ERR2_UP(idx) = 0
         THEO_ERR2_DOWN(idx) = 0
      enddo

      do i=1,Npoints
c         do j=1,NUncert
         if (NUncert.gt.0) then
               do j=CompressIdx(1),CompressIdx(NUncert)
                  if (j.gt.0) then
                  idx = DATASETIDX(idxdataset,i)
                  
               endif
            enddo
         endif
      enddo

      do i=1,Npoints
         idx = DATASETIDX(idxdataset,i)
         THEO_TOT_UP(idx) = 0.0 !SQRT(THEO_ERR2_UP(idx))*theo_fix(idx)
     +        /100d0
         THEO_TOT_DOWN(idx) = 0.0 ! SQRT(THEO_ERR2_DOWN(idx))*theo_fix(idx)
     +        /100d0
      enddo

      close (52)

      return
 101  Call HF_ErrLog(1,'F:Can not open file '//Trim(FileName))
 102  Call HF_ErrLog(2,'F:Error reading data namelist from the file '
     $     //Trim(FileName))
 103  Call HF_ErrLog(3,'F:Namelist data not found in the file '
     $     //Trim(FileName))
 1017 Call HF_ErrLog(8,'F:Can not read theory file content '
     $     //trim(FileName))
 1018 Call HF_ErrLog(9,
     $     'F:End of theory file while expecting more lines '
     $     //trim(FileName))
 1019 Call HF_ErrLog(10,'F:Problems interpreting content of file '
     $     //trim(FileName))
C---------------------------------------------------------------------
      end


