C---------------------------------------------------
C 
!> Read steering file steer.txt
C
C---------------------------------------------------
      subroutine read_steer

      implicit none

      include 'steering.inc'
C=================================================

      call Set_Defaults  ! global defaults

C Read various namelists:
      call read_hfitternml  ! main steering FIRST
      call read_systematicsnml ! Read (optional) systematics namelist SECOND
      call read_infilesnml   ! Read data file names THIRD
      call read_ewparsnml   ! electroweak parameters
      call read_outputnml   ! output options
      call read_outdirnml   ! output dir 

      if(Itheory.lt.100) then
         call read_lhapdfnml    ! read lhapdf 
C
C Decode PDF type:
C      
         call SetPDFType
C
C Decode PDF style:
C      
         call SetPDFStyle
      endif   ! Itheory < 100

      call read_mcerrorsnml  ! MC uncertainties
      call read_chebnml      ! chebyshev parameterisation extra pars
      call read_polynml
      call read_hqscalesnml  ! read HQ scales

      if (itheory.ge.100) then
         call read_ccfmfilesnml
      endif
      
      call Read_InCorrNml   ! Covariance matrix
      call read_scalesnml   ! Read scales namelist
c WS 2013-01-07 always read CSOffsetNML
      ! if(CorrSystByOffset) then
        call Read_CSOffsetNML   ! Offset method parameters
      ! endif

      if(Itheory.lt.100) then
C
C Also read extra minuit parameters:
C      
         call readextraparam
      endif ! Itheory > 100

C 07/12/2011 Dipole ===>
      call SetDipoleType

C 09/01/2013 Check consistency of the input
      call CheckInputs

      end

!> Set default values for stearable variables.
      subroutine Set_Defaults
C ===========================================
C
C Set default values for stearable variables.
C
C -------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'couplings.inc'
      include 'pdflength.inc'
      include 'pdfparam.inc'
      include 'datasets.inc'
      include 'systematics.inc'
      include 'reweighting.inc'
      include 'scales.inc'
      include 'indata.inc'
      include 'for_debug.inc'
      include 'extrapars.inc'

      integer i
C------------------------------------------------------
*     ------------------------------------------------
*     Initialise basic parameters
*     ------------------------------------------------


      nExtraParam = 0

      Itheory = 0
      EWFIT=0

      iDH_MOD = 0  ! no Dieter Heidt modifications to stat. errros.

      PDFStyle  = 'HERAPDF'
      PDFType  = 'proton'

      H1QCDFUNC= .False.
C=================================================


C  PDF length on/off:
      ILENPDF = 0

C PDF length weight factor:
      do i=1,5
         pdfLenWeight(i) = 0.
      enddo
C Chebyshev param. of the gluon:
      NCHEBGLU = 0

C Chebyshev param. of the Sea:
      NCHEBSEA = 0

C Offset for the Sea chebyshev parameters (default:20)
      IOFFSETCHEBSEA = 20

C Type of Chebyshev parameterization:
      ichebtypeGlu = 0
      ichebtypeSea = 0

      Chi2MaxError = 1.E10  ! turn off.

C     Initialise LHAPDF parameters
      LHAPDFSET = 'cteq65.LHgrid'
      ILHAPDFSET = 0
      IPDFSET = 1

C 25 Jan 2011
C     Pure polinomial param for the valence quarks:
C     by default starting from N=61 for Uv and N=71 for Dv
      NPOLYVAL = 0

C Add option to change Z of valence PDFs at x=1 from  (1-x) to (1-x)^2
      IZPOPOLY = 1

C Square polynom before calculating dv,uv. This forces positivity
      IPOLYSQR = 0
C  Key for W range 
      WMNlen =  20.
      WMXlen = 320.

      chebxmin = 1.E-5

C  Hermes-like strange (off by default):
      ifsttype = 0

C  Cache PDF calls
      CachePDFs     = .false.

! Do not split the data into fit and control sub-samples:
      ControlFitSplit = .false.

C  Fast applgrid:     
      LFastAPPLGRID = .false.
      LUseAPPLgridCKM = .true.
* 
C MC Errors defaults:
      lRAND = .false.
      lRandData = .true.
      iSEEDmc = 0
      STATYPE = 1
      SYSTYPE = 1

C PDF output options:

c 2012-11-08 WS: set default for DoBands
      DoBands = .false.
      outnx = 101
      do i=1,NBANDS
       Q2VAL(i) = -1.
      enddo
      outxrange(1) = 1e-4
      outxrange(2) = 1.0

      strange_frac = 0.31
      charm_frac = 0.00

      mch        = 1.4D0
      mbt        = 4.75D0
      mtp        = 174.D0

      HFSCHEME = 0

      OutDirName  = 'output'
      UseGridLHAPDF5=.false.
      WriteLHAPDF6=.true.

      Debug = .false.

C
C Names of syst. error sources:
C
      do i=1,NSYS
         System(i) = ' '
      enddo
      end


C =============================================
C
!> Read the main steering namelisit 
C----------------------------------------------
      subroutine Read_HFitternml  

      implicit none
      
      include 'ntot.inc'
      include 'datasets.inc'
      include 'steering.inc'
      include 'scales.inc'
      include 'indata.inc'
      include 'for_debug.inc'
C-----------------------------------------------

      character*32 Chi2SettingsName(5)
      character*32 Chi2Settings(5)
      character*32 Chi2ExtraParam(8)

      real*8 Q02       ! Starting scale
      integer IOrder   ! Evolution order
      character*8 Order  ! 
      character*16 TheoryType
      integer i

C Main steering parameters namelist
      namelist/HERAFitter/
     $     ITheory, IOrder,         ! keep for backward compatibility
     $     Q02, HF_SCHEME, PDFStyle, PDFType, 
     $     LDebug, ifsttype,  LFastAPPLGRID, LUseAPPLgridCKM,
     $     Chi2MaxError, EWFIT, iDH_MOD, H1qcdfunc, CachePDFs, 
     $     ControlFitSplit,Order,TheoryType,
     $     Chi2SettingsName, Chi2Settings, Chi2ExtraParam,
     $     AsymErrorsIterations

C--------------------------------------------------------------

C Some defaults
      Order     = ' '
      TheoryType = ' '
      Chi2SettingsName(1) = 'undefined' ! triggering the old style chi2 settings
      do i=1, 8
         Chi2ExtraParam(i) = 'undefined'
      enddo
      AsymErrorsIterations = 0
C
C  Read the main HERAFitter namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=HERAFitter,END=141,ERR=42)
      close (51)

C Backward compatibility for b0.3 release !!!!
      goto 142
 141  continue
      close (51)

!      call HF_ErrLog(13011501,
!     $ 'W:WARNING: Using obsolete h1fitter namelist.'//
!     $     ' Will be deprecated with v1.0!')
!      open (51,file='steering.txt',status='old')
!      read (51,NML=H1Fitter,END=41,ERR=42)
!      close (51)
 142  continue
C  End of backward compatibility !!!

C 
      if (AsymErrorsIterations .gt. 0) then
         call hf_errlog(13080601,'I: Use asymmetric uncertainties')
      else
         call hf_errlog(13080601,
     $        'I: Symmetrise asymmetric uncertainties if present')
      endif


C Decode computation order:
      if (Order.ne.' ') then
         Call DecodeOrder(Order)
      else
         I_FIT_ORDER = IOrder
      endif
      
C Decode theory type:
      if (TheoryType.ne.' ') then
         Call DecodeTheoryType(TheoryType)
      endif 

C     set debug flag used elsewhere according to steering
      Debug = lDebug
C
C Decode Chi2 style:
C

      call SetChi2Style(Chi2SettingsName, Chi2Settings, 
     $     Chi2ExtraParam)

      if (itheory.lt.100) then
C
C Decode HFSCHEME:
C      
         call SetHFSCHEME
      endif

      starting_scale = Q02

      if (LDebug) then
C Print the namelist:
         print HERAFitter
      endif

      return

 41   continue
      print '(''Namelist &HERAFitter NOT found'')'
      call HF_stop
 42   continue
      print '(''Error reading namelist &HERAFitter, STOP'')'
      call HF_stop
      end

C---------------------------------------- 
!> Read electroweak parameters
C-----------------------------------------
      subroutine read_ewparsnml

      implicit none
C Namelist for EW parameters:
      include 'couplings.inc'
      include 'steering.inc'

      namelist/EWpars/alphaem, gf, sin2thw, convfac,
     $ Mz, Mw, Mh, wz, ww, wh, wtp,
     $ Vud, Vus, Vub, Vcd, Vcs, Vcb, Vtd, Vts, Vtb,
     $ men, mel, mmn, mmo, mtn, mta, mup, mdn,
     $ mch, mst, mtp, mbt
C--------------------------------------------------
      open (51,file='ewparam.txt',status='old')
      read (51,NML=EWpars,END=43,ERR=44)
      close (51)

      HF_MASS(1) = mch
      HF_MASS(2) = mbt
      HF_MASS(3) = mtp

* --- Check the consistency of the steering file

      if (HFSCHEME.eq.1.and.HF_MASS(2)**2.lt.starting_scale) then
       write(6,*)
       write(6,*) 'Bottom thres. has to be larger than starting scale'
       write(6,*)
       call HF_stop
      endif

      if (HFSCHEME.eq.1.and.HF_MASS(2).lt.HF_MASS(1)) then
       write(6,*)
       write(6,*) 'Bottom thres. has to be larger than charm thres.'
       write(6,*)
       call HF_stop
      endif

      if (LDebug) then
C Print the namelist:
         print EWpars
      endif

      return

 43   continue
      print '(''Namelist @EWPars NOT found, STOP'')'
      call HF_stop

 44   continue
      print '(''Error reading namelist @EWPars, STOP'')'
      call HF_stop
      end

C-------------------------------------------------------
!> Read InCorr namelist
C-------------------------------------------------------
      subroutine Read_InCorrNml

      implicit none
      include 'ntot.inc'
      include 'datasets.inc'
      include 'steering.inc'
C Namelist for statistical correlations to read
      namelist/InCorr/NCorrFiles,CorrFileNames
C----------------------------------------------------------
C     
C  Read statistical correlations namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=InCorr,END=76,ERR=77)
      print '(''Read '',I4,'' correlation files'')',NCorrFiles
 136  continue
      close (51)
      
      if (LDebug) then
C Print the namelist:
         print InCorr
      endif

      return

 76   print '(''Namelist &InCorr NOT found'')'
      close(51)
      Return
 77   continue
      print '(''Error reading namelist &InCorr, STOP'')'
      call HF_stop

      end

C
!> Read Heavy Quark (HQ) scales
C-------------------------------------------------------
      subroutine read_hqscalesnml

      implicit none
      include 'steering.inc'
      character*32 MassHQ
C-------------------------------------------------------
C (Optional) set HQ scale
        namelist/HQScale/MassHQ,scalea1,scaleb1
C-------------------------------------------------------

C scale for HQ
       scalea1 = 1
       scaleb1 = 0
      MassHQ = 'mc'
C
C  Read the HQScale namelist:
C
      open (51,file='steering.txt',status='old') 
      read (51,NML=HQScale,ERR=70,end=69) 
 69   continue
      close (51)
C
C asign mc or mb to hq scale
C      
      call SetMHSCALE(MassHQ)
       aq2 = 1/scalea1
       bq2 = -4*scaleb1/scalea1
       hqscale1in = scalea1  
       hqscale2in = scaleb1 
       if(mod(HFSCHEME,10).eq.1) then
       if(massh.eq.1) then
       print*,'factorisation scale for heavy quarks is set to  sqrt(', hqscale1in,'*Q^2 + ',hqscale2in , '* 4m_c^2 )'   
       elseif(massh.eq.2) then
       print*,'factorisation scale for heavy quarks is set to  sqrt(', hqscale1in,'*Q^2 + ',hqscale2in , '* 4m_b^2 )'   
        endif   
        endif
      if (LDebug) then
C Print the namelist:
         print HQScale
      endif
         
      return

 70   continue
      print '(''Error reading namelist &HQScale, STOP'')'
      call HF_stop
      end

C
!> Read lhapdf and reweighting namelists
C----------------------------------------
      subroutine read_lhapdfnml

      implicit none
      include 'steering.inc'
      include 'reweighting.inc'
C------------------------------------
C (Optional) LHAPDF steering card
      namelist/lhapdf/LHAPDFSET,ILHAPDFSET,
     $     LHAPDFErrors,Scale68,LHAPDFVARSET,NPARVAR,
     $     WriteAlphaSToMemberPDF

C (Optional) reweighting steering card
      namelist/reweighting/FLAGRW,RWPDFSET,RWDATA
     $     ,RWMETHOD,DORWONLY,RWREPLICAS,RWOUTREPLICAS

C------------------------------------------------------------
C Reweighting defaults

      FLAGRW = .false.
      RWMETHOD = 1
      DORWONLY = .false.
      RWPDFSET = ''
      RWREPLICAS = 0
      RWOUTREPLICAS = 0

C LHAPDFErrors default
      LHAPDFErrors = .false.
      Scale68 = .false.
      NPARVAR = 0
      LHAPDFVARSET = ''

C
      WriteAlphaSToMemberPDF = .false.
C
C  Read the lhapdf namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=lhapdf,ERR=67,end=68)
 68   continue
      close (51)

C
C  Read the reweighting namelist:
C 
      open (51,file='steering.txt',status='old')
      read (51,NML=reweighting,END=75,ERR=74)
 75   continue
      close (51)
C
C check whether RWPDFSET and LHAPDF set are equal
C
      if (FLAGRW) then
         if (TRIM(RWPDFSET) .ne. TRIM(LHAPDFSET)) then
            call HF_ErrLog(12032302,'W:WARNING: Setting LHAPDF set to '
     $           //TRIM(RWPDFSET))
            LHAPDFSET=RWPDFSET
         endif

C  check if the PDFstyle is indeed Ok
         if (PDFStyle.ne.'LHAPDF' .and. PDFStyle.ne.'LHAPDFQ0') then
            call HF_Errlog(12032303,
     $           'W:WARNING: Setting PDF style to LHAPDFQ0')
            PDFStyle = 'LHAPDFQ0'
         endif
      endif

      if (LDebug) then
C Print the namelist:
         print lhapdf
         print reweighting
      endif

      return
C---
 67   continue
      print '(''Error reading namelist &lhapdf, STOP'')'
      call HF_stop

 74   continue
      print '(''Error reading namelist &reweighting, STOP'')'
      call HF_stop
      end

C
!> Read MC errors namelist
C-------------------------------------------------------
      subroutine read_mcerrorsnml

      implicit none
      include 'steering.inc'
C (Optional) MC method namelist
      namelist/MCErrors/LRand, ISeeDMC, StaType, SysType, LRandData 
C------------------------------------------------------
C
C  Read the MC method namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=MCErrors,ERR=62,end=61)
 61   continue
      close (51)

      if (LDebug) then
         print MCErrors
      endif

      return
C-----------------------------------------------
 62   continue
      print '(''Error reading namelist &MCErrors, STOP'')'
      call HF_stop
      end


C
!> Read optional chebyshev namelist
C--------------------------------------------------------
      subroutine read_chebnml

      implicit none
      include 'steering.inc'
      include 'pdflength.inc'
      include 'pdfparam.inc'
C---------------------------------------------
C (Optional) Chebyshev namelist
      namelist/Cheb/ILENPDF,pdfLenWeight,NCHEBGLU,NCHEBSEA
     $     ,IOFFSETCHEBSEA,ichebtypeGlu,ichebtypeSea
     $     ,WMNlen,WMXlen, ChebXMin
C-------------------------------------------------      
C
C  Read the Chebyshev namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=Cheb,ERR=64,end=63)
 
 63   continue
      close (51)

      chebxminlog = log(chebxmin)
      if (NCHEBGLU.ne.0) then
         print *,'Use Chebyshev polynoms for gluon with N=',NCHEBGLU
      endif

      if (NCHEBSEA.ne.0) then
         print *,'Use Chebyshev polynoms for sea with N=',NCHEBSEA
         print *,'Offset for minuit parameters is',IOFFSETCHEBSEA
      endif

      if (LDebug) then
         print Cheb
      endif

      return
C-----------------
 64   continue
      print '(''Error reading namelist &Cheb, STOP'')'
      call HF_stop
      end

      
C
!> Optional polynomial parametrisation for valence quarks
C-------------------------------------------------------------
      subroutine read_polynml

      implicit none
      include 'steering.inc'
C (Optional) Polynomial parameterisation for valence
      namelist/Poly/NPOLYVAL,IZPOPOLY,IPOLYSQR
C-------------------------------------------
C
      open (51,file='steering.txt',status='old')
      read (51,NML=Poly,ERR=66,end=65)    
 65   continue
      close (51)

      if (LDebug) then
         print Poly
      endif

      return
C--------------------------------------------------------
 66   continue
      print '(''Error reading namelist &Poly, STOP'')'
      call HF_stop

      end

C
!> Read InFiles namelist
C-------------------------------------------------------
      subroutine read_infilesnml

      implicit none
      include 'ntot.inc'
      include 'datasets.inc'
      include 'steering.inc'
      include 'scales.inc'
C---
      integer i
C Namelist for datafiles to read
      namelist/InFiles/NInputFiles,InputFileNames
C-------------------------------------------------
C  Read the data namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=InFiles,END=71,ERR=72)
      print '(''Read '',I4,'' data files'')',NInputFiles
      close (51)
C---------------------
C
C  Data-set dependent scales. First set defaults
C
      do i=1,NInputFiles
         DataSetMuR(i)    = 1.0D0
         DataSetMuF(i)    = 1.0D0
         DataSetIOrder(i) = I_Fit_Order
      enddo
C---------------------
      if (LDebug) then
         print InFiles
      endif

      return

 71   continue
      print '(''Namelist &InFiles NOT found'')'
      call HF_stop
 72   continue
      print '(''Error reading namelist &InFiles, STOP'')'
      call HF_stop

      end


C
!> Read ccfm namelist
C------------------------------------------------------
      subroutine read_ccfmfilesnml

      implicit none
C updf stuff
C Namelist for datafiles 
      include 'steering.inc'

      character*132 CCFMfilename  !> Names of input files
      namelist/CCFMFiles/CCFMfilename
      character CCFMfile*132
      Common/CCFMout/CCFMfile
      Integer idx
C-------------------------------------------------
C Read the CCFM data file name
      open (51,file='steering.txt',status='old')
      read (51,NML=CCFMFiles,END=71,ERR=72)
      close (51)
      idx=index(CCFMFilename,' ')-5
      CCFMfile = CCFMFilename(1:idx)

      if (LDebug) then
         print CCFMFilename
      endif

      return
      
 71   continue
      print '(''Namelist &CCFMFiles NOT found'')'
      call HF_stop
 72   continue
      print '(''Error reading namelist &CCFMFiles, STOP'')'
      call HF_stop

      end


C--------------------------------------------------------
!>  Read the scales namelist:
C---------------------------------------------------------
      subroutine read_scalesnml

      implicit none
      include 'ntot.inc'
      include 'scales.inc'
      include 'steering.inc'
C (Optional) Data-set dependent scales
      namelist/Scales/DataSetMuR,DataSetMuF,DataSetIOrder
C---------------------------------------------
      open (51,file='steering.txt',status='old')
      read (51,NML=Scales,END=123,ERR=124)
 123  Continue
      close (51)

      if (LDebug) then
         print Scales
      endif

      return

 124  print '(''Error reading namelist &scales, STOP'')'
      Call HF_stop

      end


C
!> Read output namelist
C------------------------------------------------
      subroutine read_outputnml

      implicit none
      include 'steering.inc'
      integer i, ilastq2

C Output style namelist
      namelist/Output/DoBands, Q2VAL, OutNX, OutXRange,
     $                      UseGridLHAPDF5, WriteLHAPDF6

C--------------------------------------------------------
C  Read the output namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=Output,END=51,ERR=52)
      close (51)

      ilastq2 = NBANDS/2
      do i=NBANDS/2,1,-1
         if (q2val(i).lt.0) ilastq2 = i-1
      enddo
      do i=1,NBANDS/2
       Q2VAL(i+ilastq2) = Q2VAL(i+NBANDS/2)
      enddo
c      print *,'q2val ', (q2val(i),i=1,NBANDS)

      if (LDebug) then
         print Output
      endif

      return
 51   continue
      print '(''Namelist &Output NOT found'')'
      return
 52   continue
      print '(''Error reading namelist &Output, STOP'')'
      call HF_stop

      end


C
!> Read output dir name
C------------------------------------------------
      subroutine read_outdirnml

      implicit none
      include 'steering.inc'
      namelist/OutDir/OutDirName, LHAPDF6OutDir
      LOGICAL ex
C--------------------------------------------------------
C  Read the OutDir namelist:
C
      LHAPDF6OutDir='herapdf'
      open (51,file='steering.txt',status='old')
      read (51,NML=OutDir,END=152,ERR=56)
 152  continue
      close (51)

C check if limit of 22 char is not exceeded:      
      if(LEN(TRIM(OutDirName)).gt.22) then
          call hf_errlog(09092013,
     $   'F: Name of result directory is too long (max is 22 char) ')
          call hf_stop
      endif

      inquire(FILE=TRIM(OutDirName),EXIST=ex)
      if(ex) then
          call hf_errlog(250420131,
     $   'I: Results written to existing directory: '//TRIM(OutDirName))
      else 
          call hf_errlog(250420132,
     $     'I: Creating directory to store results: '//TRIM(OutDirName))
          CALL system('mkdir -p '//TRIM(OutDirName))
      endif


      if (LDebug) then
         print OutDir
      endif

      return
 56   continue
      print '(''Error reading namelist &OutDir, STOP'')'
      call HF_stop

      end

C---------------------------------------
C
!>  Set PDF parameterisation type
C
C---------------------------------------
      Subroutine SetPDFType()

      implicit none
      include 'steering.inc'


      if (PDFType.eq.'proton'.or. PDFType.eq.'PROTON') then
         lead = .false.
         deuteron = .false.
         print *,'Fitting for PROTON PDFs, PDFType=', PDFType
      elseif (PDFType.eq.'lead'.or. PDFType.eq.'LEAD') then
         lead = .true. 
         deuteron = .false.
         print *,'Fitting for LEAD PDFs, PDFType=', PDFType
      elseif (PDFType.eq.'DEUTERON'.or. PDFType.eq.'deuteron') then
         lead = .true. 
         deuteron = .true. 
         print *,'Fitting for DEUTERON PDFs, PDFType=', PDFType
      else
         call hf_errlog(300920131,
     $   'F: Unsupported PDFType used!')
      endif
      end
C---------------------------------


C---------------------------------------
C
!>  Set PDF parameterisation style
C
C---------------------------------------
      Subroutine SetPDFStyle()

      implicit none

      logical lhapdffile_exists
      include 'steering.inc'
C---------------------------------

      FlexibleGluon = .false.
      
      if (PDFStyle.eq.'10p HERAPDF'.or.
     $     PDFStyle.eq.'13p HERAPDF'.or.
     $     PDFStyle.eq. 'HERAPDF'.or.
     $     PDFStyle.eq. 'strange') then
         iparam = 2011
         FlexibleGluon = .true.
         PDF_DECOMPOSITION = 'Dv_Uv_Dbar_Ubar_Str'
 
      elseif (PDFStyle.eq.'CTEQHERA') then
         PDF_DECOMPOSITION = 'Dv_Uv_Dbar_Ubar_Str'

      elseif (PDFStyle.eq.'CTEQ') then
         iparam = 171717
         PDF_DECOMPOSITION = 'Dv_Uv_Dbar_Ubar'


      elseif ((PDFStyle.eq.'AS').or.(PDFStyle.eq.'BiLog')) then
         iparam = 1977
         PDF_DECOMPOSITION = 'Dv_Uv_Dbar_Ubar'

      elseif (PDFStyle.eq.'CHEB') then
         iparam = 4
         PDF_DECOMPOSITION = 'Dv_Uv_Sea_Delta'

      elseif (PDFStyle.eq.'LHAPDFQ0') then
         iparam = 0
         PDF_DECOMPOSITION = 'LHAPDF'

      elseif (PDFStyle.eq.'LHAPDF') then
         iparam = 0
         PDF_DECOMPOSITION = 'LHAPDF'

      elseif (PDFStyle.eq.'DDIS') then
         iparam = 301        
         PDF_DECOMPOSITION = 'Diffractive'
      elseif (PDFStyle.eq.'QCDNUM_GRID') then
         PDF_DECOMPOSITION = 'QCDNUM_GRID'
      else
         print *,'Unsupported PDFStyle =',PDFStyle
         print *,'Check value in steering.txt'
         call HF_stop
      endif

      if ((PDFStyle.eq.'LHAPDF').or.(PDFStyle.eq.'LHAPDFQ0')) then
         INQUIRE(FILE=LHAPDFSET, EXIST=lhapdffile_exists) 
         if(lhapdffile_exists) then
            call InitPDFset(LHAPDFSET)
         else
            call InitPDFsetByName(LHAPDFSET)
         endif

      ! Get number of sets:
         call numberPDF(nLHAPDF_Sets)  
         
         
         call InitPDF(ILHAPDFSET)

         if(PDFStyle.eq.'LHAPDF') then
            IPDFSET = 5
         endif
      endif

      end



C---------------------------------------
C
!>  Set Heavy Flavour Scheme
C
C---------------------------------------
      Subroutine SetHFSCHEME

      implicit none
      include 'steering.inc'
C---------------------------------
      
      if (HF_SCHEME.eq.'ZMVFNS') then
          HFSCHEME = 0
      elseif (HF_SCHEME.eq.'ACOT ZM') then
          HFSCHEME = 1 
      elseif (HF_SCHEME.eq.'ACOT Full') then
          HFSCHEME = 11 
      elseif (HF_SCHEME.eq.'ACOT Chi') then
          HFSCHEME = 111 
      elseif (HF_SCHEME.eq.'RT') then
          HFSCHEME = 2
      elseif (HF_SCHEME.eq.'RT FAST') then
          HFSCHEME = 22 
      elseif (HF_SCHEME.eq.'RT OPT') then
          HFSCHEME = 202
      elseif (HF_SCHEME.eq.'RT OPT FAST') then
          HFSCHEME = 222 
      elseif (HF_SCHEME.eq.'FF') then
          HFSCHEME = 3 
      elseif (HF_SCHEME.eq.'FF ABM') then
         HFSCHEME = 4 
      elseif (HF_SCHEME.eq.'BMSN ABM') then
         HFSCHEME = 44 
      elseif (HF_SCHEME.eq.'FF ABM RUNM') then
         HFSCHEME = 444
      else
         print *,'Unsupported HFSCHEME =',HF_SCHEME
         print *,'Check value in steering.txt'
         call HF_stop
      endif
      end



C---------------------------------------
C
!>  Set HQ scale parameter
!>  @param MassHQ heavy quark mass
C---------------------------------------
      Subroutine SetMHSCALE(MassHQ)

      implicit none
      character*(*) MassHQ
      include 'steering.inc'
C---------------------------------

      if (MassHQ.eq.'mc') then
          MASSH = 1 
      elseif (MassHQ.eq.'mb') then
          MASSH = 2
      else
         print *,'Unsupported MassHQ =',MassHQ
         print *,'Check value in steering.txt'
         call HF_stop
      endif

      end



C---------------------------------------
C
!>  Set Chi2 style
!>  @param Chi2SettingsName bias corrections for uncertainties and treatment of systematics in chi2 
!>  @param Chi2Settings values corresponding to each of Chi2SettingsName parameters
!>  @param Chi2ExtraParam extra corrections in chi2
C---------------------------------------
      Subroutine SetChi2Style(Chi2SettingsName, Chi2Settings, 
     $     Chi2ExtraParam)

      implicit none
      character*32 Chi2SettingsName(5)
      character*32 Chi2Settings(5)
      character*32 Chi2ExtraParam(8) 
      integer i
      include 'steering.inc'
C---------------------------------

      if (Chi2SettingsName(1).eq.'undefined') then
C
C  Reset defaults if Chi2SettingsName parameter is not set.
C 
         CorrSystByOffset=.false.
         CorSysScale = 'Linear'
         StatScale   = 'Poisson'
         UncorSysScale = 'Linear'
         CorChi2Type = 'Hessian'

      else
         ! CorrSystByOffset=.false.
c     $     ,StatScale, UncorSysScale, CorSysScale,UncorChi2Type,CorChi2Type
         do i=1, 5
            if(Chi2SettingsName(i).eq.'StatScale') then
               StatScale = Chi2Settings(i)
            elseif(Chi2SettingsName(i).eq.'UncorSysScale') then
               UncorSysScale = Chi2Settings(i)
            elseif(Chi2SettingsName(i).eq.'CorSysScale') then
               CorSysScale = Chi2Settings(i)
            elseif(Chi2SettingsName(i).eq.'UncorChi2Type') then
               UncorChi2Type = Chi2Settings(i)
            elseif(Chi2SettingsName(i).eq.'CorChi2Type') then
               CorChi2Type = Chi2Settings(i)
            else
               print *,'Unsupported Chi2SettingsName =',Chi2SettingsName(i)
               call HF_stop
            endif
         enddo

C some defaults
         Chi2PoissonCorr = .false.
         Chi2FirstIterationRescale = .false.
         Chi2ExtraSystRescale = .false.
         do i=1, 8
            if(Chi2ExtraParam(i).eq.'PoissonCorr') then
               Chi2PoissonCorr = .true.
            elseif(Chi2ExtraParam(i).eq.'FirstIterationRescale') then
               Chi2FirstIterationRescale = .true.
            elseif(Chi2ExtraParam(i).eq.'ExtraSystRescale') then
               Chi2ExtraSystRescale = .true.
            elseif(Chi2ExtraParam(i).ne.'undefined') then
               print *,'Unsupported Chi2ExtraParam = ',Chi2ExtraParam(i)
               call HF_stop
            endif
         enddo
      endif
      end
      

C
!> Read ExtraMinimisationParameters namelists
C-------------------------------------
      Subroutine ReadExtraParam

      implicit none
      include 'extrapars.inc'
      integer maxExtra
      parameter (maxExtra=50)
      character*32 name(maxExtra)
      double precision Value(maxExtra),Step(maxExtra)
     $     ,Min(maxExtra),Max(maxExtra)

      namelist/ExtraMinimisationParameters/Name,Value,Step,Min,Max
      integer i
C----------------------------------------
      
      open (51,file='steering.txt',status='old')
C
C Read as many instances of the namelist as exists:
C
      do while (.true.)
C
C Reset names
C
         do i=1,maxExtra
            name(i) = ' '
         enddo
         read (51,NML=ExtraMinimisationParameters,END=71,ERR=72)
         
         do i=1,maxExtra
            if (name(i).ne.' ') then
               call AddExternalParam(name(i),value(i), step(i), min(i), max(i))
            endif
         enddo
      enddo
 71   continue
      print '(''Got '',i5,'' extra minuit parameters'')',nExtraParam
      close (51)
      return
 72   continue
      print *,'Problem reading namelist ExtraMinimisationParameters'
      call HF_stop

C----------------------------------------
      end

C
!> Add extra fitting parameters
!> @param name of extra paramet
!> @param value of extra parameter
!> @param step gradient of parameter in case of fitting
!> @param min, max range of allowed values in case of fitting
C-----------------------------------------------
      Subroutine AddExternalParam(name, value, step, min, max)

      implicit none
      include 'extrapars.inc'
      character*(*) name
      double precision value, step, min, max
C---------------------------------------------
C Add extra param
C
      nExtraParam = nExtraParam + 1
      if (nExtraParam.gt. nExtraParamMax) then
         print *,'Number of extra parameters exceeds the limit'
         print *,'nExtraParam=',nExtraParam
         print *,'Check your steering, stopping'
         call HF_stop
      endif
      ExtraParamNames(nExtraParam) = name
      ExtraParamValue(nExtraParam) = value
      ExtraParamStep (nExtraParam) = step
      ExtraParamMin  (nExtraParam) = min
      ExtraParamMax  (nExtraParam) = max
      end


C-----------------------------------------
C
!> Read optional systematics namelist
C
C-----------------------------------------
      Subroutine read_systematicsnml

      implicit none
      include 'ntot.inc'
      include 'systematics.inc'
      include 'steering.inc'
      character*64 ListOfSources(nsysmax),ScaleByNameName(nsysmax)
      double precision ScaleByNameFactor(nsysmax)
      
      namelist/ Systematics/ListOfSources,ScaleByNameName
     $     ,ScaleByNameFactor
      integer i,ii
C----------------------------------------

C Initialisation:
      nsys = 0
      do i=1,nsysmax
         SysScaleFactor(i) = 1.0D0
         ListOfSources(i) = ' '
         ScaleByNameName(i) = ' '

 ! Set default scaling behaviour:
         if (CorSysScale .eq. 'Linear' ) then
            SysScalingType(i)  =  isLinear
         else if (CorSysScale .eq. 'NoRescale') then
            SysScalingType(i)  =  isNoRescale
         else if (CorSysScale .eq. 'LogNorm') then
            SysScalingType(i)  =  isLogNorm
            call hf_errlog(251120122,
     $           'F:LogNormal rescaling not included yet')
         else
            print *,'Unknown correlated systematics scaling behaviour'
            print *,'CorSysScale=',CorSysScale
            print *,'Check your steering'
            call hf_errlog(25112012,
     $           'F:Wrong CorSysScale value from the steering')
            call hf_stop
         endif

   !  Set nuisance parameter behaviour:
         if (CorChi2Type .eq. 'Hessian') then
            SysForm(i)         =  isNuisance
         elseif (CorChi2Type .eq. 'Offset') then
            SysForm(i)         =  isOffset
         elseif (CorChi2Type .eq. 'Matrix') then
            SysForm(i)         =  isMatrix
         else
            print *,'Unknown correlated systatics treatment'
            print *,'CorChi2Type=',CorChi2Type
            print *,'Check your steering'
            call hf_errlog(251120123,
     $           'F:Wrong CorChi2Type value from the steering')
            call hf_stop
         endif

      enddo
      
      open (51,file='steering.txt',status='old')
      read (51,NML=Systematics,END=123,ERR=124)

      if (LDebug) then
         print Systematics
      endif
C----
C Decode:
C
      do i=1,nsysmax
         if (ListOfSources(i).ne.' ') then
            Call AddSystematics(ListOfSources(i))
         else
            goto 77
         endif
      enddo
 77   continue

      do i=1,nsysmax
         if (ScaleByNameName(i).ne.' ') then
            do ii=1,nsys
               if (ScaleByNameName(i) .eq. System(ii)) then
                  SysScaleFactor(ii) = ScaleByNameFactor(i)
               endif
            enddo
         endif
      enddo
 90   continue

 123  Continue
      close (51)
      return

 124  print '(''Error reading namelist &systematics, STOP'')'
      Call HF_stop

C----------------------------------------
      end

C=========================================
C
!> Read optional CSOffset namelist
C
C-----------------------------------------
      Subroutine Read_CSOffsetNML

      implicit none
      include 'ntot.inc'
      include 'systematics.inc'
      include 'steering.inc'
      namelist/CSOffset/ CorSysIndex, UsePrevFit
      ! .................................

      ! --- Initialisation:
      UsePrevFit = 0   ! Do not use previous fit results
      CorSysIndex = NSYSMAX+1  ! trick to calculate all offsets in one job
      
      open (51,file='steering.txt',status='old')
      read (51,NML=CSOffset,END=123,ERR=124)

      if (LDebug) then
         print CSOffset
      endif

 123  Continue
      close (51)
      return

 124  print '(''Error reading namelist &CSOffset, STOP'')'
      Call HF_stop

C----------------------------------------
      end


C
!> Decode computation order
!> @param Order of theoretical calculation
C
      subroutine DecodeOrder(Order)

      implicit none
      character*(*) Order
      include 'steering.inc'
C--------------------------------------------------
      if (Order.eq.'LO') then
         I_Fit_Order = 1
      else if (Order.eq.'NLO') then
         I_Fit_Order = 2
      else if (Order.eq.'NNLO') then
         I_Fit_Order = 3
      else
         print *,'Unknown computation order = ',order
         print *,'Check your steering.txt file'
         call hf_stop
      endif
C--------------------------------------------------
      end

C--------------------------------------------------
!> Decode type of theory to be used
!> @param TheoryType name of theory
C--------------------------------------------------
      subroutine DecodeTheoryType(TheoryType)

      character*(*) TheoryType
      include 'steering.inc'
C------------------------------------------------
      if (TheoryType.eq.'DGLAP') then
         iTheory =  0
      else if ( TheoryType.eq.'DIPOLE') then
      else if ( TheoryType.eq.'FRACTAL') then
         iTheory = 50
      else if ( TheoryType.eq.'uPDF') then
         iTheory = 101
      else if ( TheoryType.eq.'uPDF1') then
         iTheory = 101
      else if ( TheoryType.eq.'uPDF2') then
         iTheory = 102
      else if ( TheoryType.eq.'uPDF3') then
         iTheory = 103
      else if ( TheoryType.eq.'uPDF4') then
         iTheory = 104
      else if ( TheoryType.eq.'uPDF5') then
         iTheory = 105
      else
         print *,'Unknown TheoryType = TheoryType'
         print *,'Check your steering.txt file'
         call hf_stop
      endif
      
      end

C
!> Check if the systematic source is already on the list. 
!> Takes care of asymmetric errors and : modifier.
C     
      integer Function SystematicsExist(SourceName) 

      implicit none
      character*(*) SourceName
      include 'ntot.inc'
      include 'systematics.inc'
      integer j,i,iasym
      character*64 Name
C----------------------------------------------------------------
      SystematicsExist = 0

C Check for +- signs:      
      if ( SourceName( len_trim(Sourcename):len_trim(Sourcename))
     $     .eq.'+' .or.
     $     SourceName( len_trim(Sourcename):len_trim(Sourcename))
     $     .eq.'-')
     $           then
         Name = SourceName(1:len_trim(Sourcename)-1)
      else
         Name = SourceName
      endif

C Check for :
      i = index(Name,':')
      if (i.ne.0) then
         Name = Name(1:i-1)
      endif

      do j=1,NSYS            
         if ( system(j) .eq. Name ) then
            SystematicsExist = j
            Return
         endif
      enddo    
C----------------------------------------------------------------
      end

C-----------------------------------------------------------------------------
!> Add systematic source
C
!>  Detect "+" and "-" signs, at the end of source name, for asymmetric errors
!>
!>  Detect ":" modifiers
!>
!>   :M  - "multiplicative"
!>
!>   :A  - "additive"
!>
!>   :P  - "poisson"
!>
!>   :N  - "nuisance"   -- use nuisance parameters
!>
!>   :C  - "covariance" -- use covariance matrix
!>
!>   :O  - "offset"     -- use offset method for error propagation.
!>
!>   :E  - "external"   -- use minuit to minimise.
C
!> @param SName name of added systematic source. 
C-----------------------------------------------------------------------------
      Subroutine AddSystematics(SName)

      implicit none
      include 'ntot.inc'
      include 'systematics.inc'
      character*(*) SName

      character*64 SourceName
      
      integer ii,iasym
C-----------------------------------------
      
      SourceName = SName

      nsys = nsys + 1
      if (NSYS.gt.NSysMax) then
         print 
     $        '(''ReadDataFile Error: exeeding NSysMax'')'
         print '(''Current NSysMax='',i6)',NSysMax
         print '(''Increase NSysMax in systematics.inc'')'
         call HF_stop
      endif

c     Initialise iasym to 0
      iasym = 0
C
C Detect "+" and "-" signs
C
      if ( SourceName(len_trim(SourceName):len_trim(SourceName)).eq.'+') then
         iasym = len_trim(SourceName)
      endif

      if ( SourceName(len_trim(SourceName):len_trim(SourceName)).eq.'-') then
         iasym = len_trim(SourceName)
      endif

      ii = index(SourceName,':')
      if (ii.eq.0) then
         if (iasym.gt.0) then
            System(nsys) = SourceName(1:iasym-1)
         else
            System(nsys) = SourceName
         endif
      else
         if (iasym.gt.0) then
            System(nsys) = SourceName(1:iasym-1)
         else
            System(nsys) = SourceName(1:ii-1)
         endif
      endif

      do while (ii.gt.0) 
         if ( SourceName(ii+1:ii+1) .eq.'A' ) then
            SysScalingType(nsys) = isNoRescale
            Call HF_errlog(12090001,
     $           'I: Some systematic sources are additive')
         elseif ( SourceName(ii+1:ii+1) .eq.'M' ) then
            SysScalingType(nsys) = isLinear
         elseif ( SourceName(ii+1:ii+1) .eq.'P' ) then
            SysScalingType(nsys) = isPoisson
         elseif ( SourceName(ii+1:ii+1) .eq.'N' ) then
            SysForm(nsys) = isNuisance
         elseif ( SourceName(ii+1:ii+1) .eq.'C' ) then
            SysForm(nsys) = isMatrix
         elseif ( SourceName(ii+1:ii+1) .eq.'O' ) then
            SysForm(nsys) = isOffset
         elseif ( SourceName(ii+1:ii+1) .eq.'E' ) then
            SysForm(nsys) = isExternal
         else
            print *,'WARRNING: Unknown systematics modifier ',
     $            SourceName(ii+1:)
            Call HF_errlog(12090002,
     $'W:WARNING: wrong form or bias correction for a systematic source')
         endif
         
         SourceName = SourceName(ii+2:)
         ii = index(SourceName,':')
      enddo

C Register external systematics:
      if ( SysForm(nsys) .eq. isExternal) then
         call AddExternalParam(System(nsys),0.0D0, 1.0D0, 0.0D0, 0.0D0)
      endif

      end

 !>
 !>  Check consistency of the data input, abort for unsupported combinations
 !>
      Subroutine CheckInputs

      implicit none
      include 'steering.inc'
      character*48 CMess 
C----------------------------------------------------------
!      if ( I_Fit_order .eq. 1 ) then
!         if ( index(HF_SCHEME,'RT').gt.0 ) then
!            CMess = 'RT scheme does not support LO evolution'
!            goto 998
!         endif
!      endif


      if (LHAPDFErrors) then
         if(PDFStyle.ne.'LHAPDF'.and.PDFStyle.ne.'LHAPDFQ0') then
            call HF_Errlog(03062013,
     $ 'W:WARRNING PDFstyle is not LHAPDF, setting PDFErrors to False')
             LHAPDFErrors = .false.
         endif
      endif

      return
 998  continue
      call HF_ERRLOG(13010901,'F: Inconsistent steering: '//CMess)
      call hf_stop
      end
