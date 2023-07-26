C---------------------------------------------------
C
!> Read steering file steering.txt
C
C---------------------------------------------------
      subroutine read_steer

      implicit none



#include "steering.inc"
#include "ntot.inc"
#include "indata.inc"
#include "systematics.inc"
C=================================================

      call Set_Defaults  ! global defaults

C Read various namelists:
      call read_hfitternml  ! main steering FIRST
      call read_systematicsnml ! Read (optional) systematics namelist SECOND

C Special branch for rotation
      if (pdfrotate) then
         call read_theoryfilesNML
         call rediagonalize(NPoints,NSys)
         call hf_stop
      endif

      call read_infilesnml   ! Read data file names THIRD
      call read_outputnml   ! output options

      call read_lhapdfnml    ! read lhapdf
      call read_chi2scan        ! read chi2scan

      call read_mcerrorsnml  ! MC uncertainties
      call read_hqscalesnml  ! read HQ scales

      call Read_InCorrNml   ! Covariance matrix
      call read_scalesnml   ! Read scales namelist
c WS 2013-01-07 always read CSOffsetNML
      ! if(CorrSystByOffset) then
        call Read_CSOffsetNML   ! Offset method parameters
      ! endif

C 30/08/2015  KK - Twist analyses
      call read_HigherTwists

      end

!> Set default values for steerable variables.
      subroutine Set_Defaults
C ===========================================
C
C Set default values for steerable variables.
C
C -------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "couplings.inc"
#include "pdflength.inc"
#include "datasets.inc"
#include "systematics.inc"
#include "scales.inc"
#include "indata.inc"
#include "for_debug.inc"
#include "extrapars.inc"

      integer i
C------------------------------------------------------
*     ------------------------------------------------
*     Initialise basic parameters
*     ------------------------------------------------


      !nExtraParam = 0

      ExtraPdfs = .false.

      iDH_MOD = 0  ! no Dieter Heidt modifications to stat. errros.

C=================================================


C PDF length weight factor:
      do i=1,5
         pdfLenWeight(i) = 0.
      enddo

      Chi2MaxError = 1.E10  ! turn off.

C     Initialise LHAPDF parameters
      LHAPDFSET = 'cteq65.LHgrid'
      ILHAPDFSET = 0
      IPDFSET = 1
      vIPDFSET = IPDFSET

! Do not split the data into fit and control sub-samples:
      ControlFitSplit = .false.
*
C MC Errors defaults:
      lRAND = .false.
      lRandData = .true.
      iSEEDmc = 0
      STATYPE = 1
      SYSTYPE = 1

C PDF output options:
      outnx = 101
      do i=1,NBANDS
       Q2VAL(i) = -1.
      enddo
      outxrange(1) = 1e-4
      outxrange(2) = 1.0

      mch        = 1.4D0
      mbt        = 4.75D0
      mtp        = 174.D0

      !OutDirName  = 'output'

      Debug = .false.

      pdfrotate = .false.
C
C Names of syst. error sources:
C
      do i=1,NSYS
         System(i) = ' '
      enddo

C Check variables for common blocks:
      steering_check = 171717
      call common_check(steering_check)

      end


C =============================================
C
!> Read the main steering namelist
C----------------------------------------------
      subroutine Read_HFitternml

      implicit none

#include "ntot.inc"
#include "steering.inc"
#include "indata.inc"
#include "for_debug.inc"
C-----------------------------------------------

      character*32 Chi2SettingsName(5)
      character*32 Chi2Settings(5)
      character*32 Chi2ExtraParam(8)
      integer i

C Main steering parameters namelist
      namelist/xFitter/
     $     LDebug,
     $     Chi2MaxError, iDH_MOD, 
     $     ControlFitSplit,
     $     Chi2SettingsName, Chi2Settings, Chi2ExtraParam,
     $     AsymErrorsIterations, pdfRotate, UseDataSetIndex

C--------------------------------------------------------------

C     Some defaults
      UseDataSetIndex = .false.
      Chi2SettingsName(1) = 'undefined' ! triggering the old style chi2 settings
      do i=1, 8
         Chi2ExtraParam(i) = 'undefined'
      enddo
      AsymErrorsIterations = 0
C
C  Read the main xFitter namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=xFitter,END=141,ERR=42)
 141  continue
      close (51)

      if (AsymErrorsIterations .gt. 0) then
         call hf_errlog(13080601,'I: Use asymmetric uncertainties')
      else
         call hf_errlog(13080601,
     $        'I: Symmetrise asymmetric uncertainties if present')
      endif

C     set debug flag used elsewhere according to steering
      Debug = lDebug
C
C Decode Chi2 style:
C

      call SetChi2Style(Chi2SettingsName, Chi2Settings,
     $     Chi2ExtraParam)

      if (LDebug) then
C Print the namelist:
         print xFitter
      endif

      return

 41   continue
      print '(''Namelist &xFitter NOT found'')'
      call HF_stop
 42   continue
      print '(''Error reading namelist &xFitter, STOP'')'
      call HF_stop
      end


C-------------------------------------------------------
!> Read InCorr namelist
C-------------------------------------------------------
      subroutine Read_InCorrNml

      implicit none
#include "ntot.inc"
#include "datasets.inc"
#include "steering.inc"
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
#include "steering.inc"
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
!> Read lhapdf namelists
C----------------------------------------
      subroutine read_lhapdfnml

      implicit none
#include "steering.inc"
C------------------------------------
C (Optional) LHAPDF steering card
      namelist/lhapdf/
     $     LHAPDFErrors,Scale68,LHAPDFVARSET,NPARVAR,
     $     DataToTheo,nremovepriors,
     $     lhapdfprofile,lhascaleprofile

      logical lhapdferrors_save

C LHAPDFErrors default

      lhapdferrors_save = lhapdferrors ! may be set by running mode.

      lhapdfprofile = .true.
      lhascaleprofile = .false.

      Scale68 = .false.
      NPARVAR = 0
      LHAPDFVARSET = ''
      NREMOVEPRIORS = 0
C
      DataToTheo = .false.
C
C  Read the lhapdf namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=lhapdf,ERR=67,end=68)
 68   continue
      close (51)

      if ( RunningMode .ne. ' ' ) then
         lhapdferrors = lhapdferrors_save
      endif

      if (LDebug) then
C Print the namelist:
         print lhapdf
      endif

      return
C---
 67   continue
      print '(''Error reading namelist &lhapdf, STOP'')'
      call HF_stop

      end

      subroutine read_chi2scan

      implicit none
#include "steering.inc"
#include "ntot.inc"
#include "theorexpr.inc"
#include "chi2scan.inc"
C------------------------------------
C (Optional) Chi2Scan steering card
      namelist/chi2scan/label,central,values,
     $     dataid,term,TheorySources,scan,pdferrors,
     $     pdfprofile,scaleprofile,
     $	   chi2lhapdfref,chi2lhapdfset,chi2lhapdfvarset,chi2nparvar,
     $     chi2parpoint,datatotheo

C Chi2Scan default
      scan = .false.
      label = ''
      term = ''
      dataid = 0
      central = 0
      pdferrors = .false.
      pdfprofile = .false.
      scaleprofile = .false.
      chi2lhapdfref = ''
      chi2lhapdfset = ''
      chi2lhapdfvarset = ''
      chi2nparvar = 0
      chi2parpoint = 0
      datatotheo =.false.

C
C  Read the chi2scan namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=chi2scan,ERR=70,end=69)
 69   continue
      close (51)

      if (LDebug) then
C Print the namelist:
         print chi2scan
      endif

      return
C---
 70   continue
      print '(''Error reading namelist &chi2scan, STOP'')'
      call HF_stop

      end
C
!> Read MC errors namelist
C-------------------------------------------------------
      subroutine read_mcerrorsnml

      implicit none
#include "steering.inc"
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
!> Read InFiles namelist
C-------------------------------------------------------
      subroutine read_infilesnml

      implicit none
#include "ntot.inc"
#include "datasets.inc"
#include "steering.inc"
#include "scales.inc"
C---
      integer i, nf
C Namelist for datafiles to read
      namelist/InFiles/NInputFiles,InputFileNames

      character*(80) cMsg

C reset defaults:
      NInputFiles = 0
      do i = 1,NSET
         InputFileNames(i) = ''
      enddo
C-------------------------------------------------
C  Read the data namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=InFiles,END=71,ERR=72)
      close (51)

C Determine how many files to process. First count them:
      nf = 0
      do i =1, NSET
         if (InputFileNames(i).ne.'') then
            nf = nf + 1
         endif
      enddo

      if ( NInputFiles.eq.0) then
         NInputFiles = nf       ! by default use all files
      else
         if (NInputFiles.gt.nf) then
            write (cMsg,
     $ '(''W: NInputFiles='',i4
     $ ,'' exceeds actual number of files='',i4,'', reset'')')
     $           NInputFiles, nf
            call hf_errlog(18030601,cMsg)
     $
            NInputFiles = nf
         endif
      endif
      print '(''Will read '',I4,'' data files'')',NInputFiles
C---------------------
C
C  Data-set dependent scales. First set defaults
C
      do i=1,NInputFiles
         DataSetMuR(i)            = 1.0D0
         DataSetMuF(i)            = 1.0D0
         DataSetIOrder(i)         = I_Fit_Order
         DataSetMaxNF(i)          = 0
         DataSetSwitchScales(4,i) = 0d0
         DataSetSwitchScales(5,i) = 0d0
         DataSetSwitchScales(6,i) = 0d0
      enddo
      UseHVFNS = .false.
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
#include "steering.inc"

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
#include "ntot.inc"
#include "scales.inc"
#include "steering.inc"
#include "datasets.inc"
#include "couplings.inc"
C (Optional) Data-set dependent scales
      integer i_fit_order_save,i,ihq
      character*8 DataSetTheoryOrder(NSet)
      namelist/Scales/DataSetMuR,DataSetMuF,DataSetIOrder,
     $     DataSetTheoryOrder,DataSetMaxNF,DataSetSwitchScales
C---------------------------------------------
      do i=1,NSet
         DataSetTheoryOrder(i) = ''
      enddo

      open (51,file='steering.txt',status='old')
      read (51,NML=Scales,END=123,ERR=124)
 123  Continue
      close (51)

C Check datasetorder
      I_Fit_Order_Save = I_Fit_Order
      do i=1,NSet
         if (DataSetTheoryOrder(i).ne.'') then
            call DecodeOrder(DataSetTheoryOrder(i))
            DataSetIOrder(i) = I_Fit_Order
         endif
C Check if the H-VFNS has to be used because some dataset
C has "DataSetMaxNF" different from zero (default).
         if(DataSetMaxNF(i).ne.0)then
C Check that MaxNF is between 3 and 6
            if(DataSetMaxNF(i).lt.3.or.
     1         DataSetMaxNF(i).gt.6) call hf_errlog(2105201601,
     2              'F: DataSetMaxNF must be between 3 and 6')
            UseHVFNS = .true.
         else
            DataSetMaxNF(i) = 6
         endif
C Check if the H-VFNS has to be used because some dataset
C has "DataSetSwitchScales" different from zero (default).
         do ihq=4,6
            if(DataSetSwitchScales(ihq,i).ne.0d0)then
               UseHVFNS = .true.
            else
               if(ihq.eq.4) DataSetSwitchScales(ihq,i) = mch
               if(ihq.eq.5) DataSetSwitchScales(ihq,i) = mbt
               if(ihq.eq.6) DataSetSwitchScales(ihq,i) = mtp
            endif
C Check that the switching scales are ordered
            if(ihq.gt.4)then
               if(DataSetSwitchScales(ihq,i).lt.
     1              DataSetSwitchScales(ihq-1,i))
     2              call hf_errlog(24081601,
     3              'F: DataSetSwitchScales must be ordered')
            endif
         enddo
      enddo
      I_Fit_Order = I_Fit_Order_Save

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
#include "steering.inc"
      integer i, ilastq2

C Output style namelist
      namelist/Output/Q2VAL, OutNX, OutXRange,
     $     WriteLHAPDF5,
     $     ReadParsFromFile, ParsFileName, CovFileName

C--------------------------------------------------------
C  Read the output namelist:
C
      WriteLHAPDF5 = .false.  ! default: no
      ReadParsFromFile = .false.
      ParsFileName = ''
      CovFileName = ''

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
C---------------------------------------
C
!>  Set HQ scale parameter
!>  @param MassHQ heavy quark mass
C---------------------------------------
      Subroutine SetMHSCALE(MassHQ)

      implicit none
      character*(*) MassHQ
#include "steering.inc"
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
#include "steering.inc"
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
c switch on the log poisson correction if ExtraSysRescale was called
               Chi2PoissonCorr = .true.
               call HF_errlog(15012601,
     $    'I: extra log corr (Poisson) activated with ExtraSystRescale')
            elseif(Chi2ExtraParam(i).ne.'undefined') then
               print *,'Unsupported Chi2ExtraParam = ',Chi2ExtraParam(i)
               call HF_stop
            endif
         enddo
      endif
      end

!> Add extra fitting parameters
!> @param name of extra parameter
!> @param value of extra parameter
!> @param step gradient of parameter in case of fitting
!> @param min, max range of allowed values in case of fitting
!> @param constrval constrain to this value in case of fitting
!> @param construnc uncertainty on constrain in case of fitting
!> @param to_gparam send to gParameters or not
C-----------------------------------------------
C As far as I understand, currently ALL parameters are treated as extra
C This routine registers a parameter in some array where MINUIT will
C find them (???)
      Subroutine AddExternalParam(name, value, step, min, max,
     $                            constrval, construnc, to_gParam
     $     ,gParam)

      implicit none
#include "extrapars.inc"
      character*(*) name
      double precision value, step, min, max, constrval, construnc
      double precision gParam
      logical to_gParam
      integer iglobal
      CHARACTER TMP
      INTEGER I
C---------------------------------------------
C Add extra param
C

      nExtraParam = nExtraParam + 1
      if (nExtraParam.gt. nExtraParamMax) then
         print *,'Number of extra parameters exceeds the limit'
         print *,'nExtraParam=',nExtraParam
         print *,'nExtraParamMax=',nExtraParamMax
         print *,'Check your steering'
         print *,'or increase NEXTRAPARAMMAX_C in include/dimensions.h'
         print *,'stopping'
         call HF_stop
      endif
      ExtraParamNames(nExtraParam) = REPEAT(' ',LEN(ExtraParamNames(nExtraParam)))
      DO I = 1, LEN(name)
        TMP = name(I:I)
        ExtraParamNames(nExtraParam)(I:I) = TMP
      ENDDO
      ExtraParamValue(nExtraParam) = value
      ExtraParamStep (nExtraParam) = step
      ExtraParamMin  (nExtraParam) = min
      ExtraParamMax  (nExtraParam) = max
      ExtraParamConstrVal  (nExtraParam) = constrval
      ExtraParamConstrUnc  (nExtraParam) = construnc

      iglobal = 0
      if ( gParam.eq.0.0) then
         iglobal = 1
      endif

C Also add it to c++ map ...
      if (to_gParam) then
         call add_To_Param_Map( gParam, ExtraParamValue(nExtraParam)
     $        ,  iglobal, ExtraParamNames(nExtraParam)//char(0))
      endif

      end


C-----------------------------------------
C
!> Read optional systematics namelist
C
C-----------------------------------------
      Subroutine read_systematicsnml

      implicit none
#include "ntot.inc"
#include "systematics.inc"
#include "steering.inc"
      character*64 ListOfSources(nsysmax),ScaleByNameName(nsysmax),
     $     PriorScaleName(nsysmax)
      double precision ScaleByNameFactor(nsysmax),
     $     PriorScaleFactor(nsysmax)

      namelist/ Systematics/ListOfSources,ScaleByNameName
     $     ,ScaleByNameFactor, PriorScaleName, PriorScaleFactor
      integer i,ii,iType
C----------------------------------------

C Initialisation:
      nsys = 0
      do i=1,nsysmax
         SysScaleFactor(i) = 1.0D0
         ListOfSources(i) = ' '
         ScaleByNameName(i) = ' '
         PriorScaleName(i) = ' '
         PriorScaleFactor(i) = 1.0D0
         SysPriorScale(i) = 1.0D0
! Set default scaling behaviour:
         if (CorSysScale .eq. 'Linear' ) then
            SysScalingType(i)  =  isLinear
         else if (CorSysScale .eq. 'NoRescale') then
            SysScalingType(i)  =  isNoRescale
         else if (CorSysScale .eq. 'Poisson') then
            SysScalingType(i)  =  isPoisson
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
         if (PriorScaleName(i).ne.' ') then
            do ii=1,nsys
               if ( PriorScaleName(i) .eq. System(ii) ) then
                  SysPriorScale(ii) = PriorScaleFactor(i)
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
#include "ntot.inc"
#include "systematics.inc"
#include "steering.inc"
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
#include "steering.inc"
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

C
!> Check if the systematic source is already on the list.
!> Takes care of asymmetric errors and : modifier.
C
      integer Function SystematicsExist(SourceName)

      implicit none
      character*(*) SourceName
#include "ntot.inc"
#include "systematics.inc"
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
!>
!>   :D  - "data" (not theory), default for data files
!>
!>   :T  - "theory" (not data), default for theory files
!>
C
!> @param SName name of added systematic source.
C-----------------------------------------------------------------------------
      Subroutine AddSystematics(SName)

      implicit none
#include "ntot.inc"
#include "systematics.inc"
      character*(*) SName

      character*64 SourceName

      integer ii,iasym
C-----------------------------------------

      SourceName = SName

      nsys = nsys + 1
      if (NSYS.gt.NSysMax) then
         print
     $        '(''ReadDataFile Error: exceeding NSysMax'')'
         print '(''Current NSysMax='',i6)',NSysMax
         print '(''Increase NSYSMAX_C in include/dimensions.h'')'
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
         System(nsys) = SourceName(1:ii-1)
      endif

      ISystType(nsys) = iDataSyst       ! Default

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
         elseif ( SourceName(ii+1:ii+1) .eq.'D' ) then
            ISystType(nsys) = iDataSyst
         elseif ( SourceName(ii+1:ii+1) .eq.'T' ) then
            ISystType(nsys) = iTheorySyst
         else
            print *,'WARNING: Unknown systematics modifier ',
     $            SourceName(ii+1:)
            Call HF_errlog(12090002,
     $'W:WARNING: wrong form or bias correction for a systematic source')
         endif

         SourceName = SourceName(ii+2:)
         ii = index(SourceName,':')
      enddo


C Register external systematics:
      if ( SysForm(nsys) .eq. isExternal) then
         call AddExternalParam(System(nsys),0.0D0, 1.0D0, 0.0D0, 0.0D0
     $                         ,0.0D0,0.0D0,.false.,0.0D0)
      endif

      end

! 30/08/2015 KK - read Higher Twist parameters.
! WS 2015-10-04 - read steering options. Parameters moved to \c ExtraMinimisationParameters.
C--------------------------------------------------------
!>  Read higher twists options from the 'HighTwist' namelist.
      subroutine read_HigherTwists

      implicit none
      include 'steering.inc'

      namelist/HighTwist/doHiTwist,HiTwistType,HiTwistSubType

      doHiTwist = .false.
      HiTwistType = 'Twist4'
      HiTwistSubType = 'lam-sig-x0'
      open (51,file='steering.txt',status='old')
      read (51,NML=HighTwist,ERR=134,end=131)

 131  continue
      close (51)

      if (doHiTwist) then
        if (LDebug) then
          print HighTwist
        endif
      endif

      return
C-----------------
134   continue
      print '(''Error reading namelist &HighTwist, STOP'')'
      call HF_stop
      end
