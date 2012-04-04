
      subroutine read_steer
C---------------------------------------------------
C 
C> Read steering file steer.txt
C
C---------------------------------------------------
      implicit none


      include 'ntot.inc'
      include 'steering.inc'
      include 'couplings.inc'
      include 'pdflength.inc'
      include 'pdfparam.inc'
      include 'datasets.inc'
      include 'systematics.inc'
      include 'nnpdf.inc'
      include 'scales.inc'
      include 'for_debug.inc'
C=================================================


C Define namelists:

      character*32 PDFStyle, Chi2Style, HF_SCHEME, MassHQ

      real*8 Q02       ! Starting scale
      integer IOrder   ! Evolution order
C Main steering parameters namelist
      namelist/H1Fitter/ITheory, IOrder, Q02, HF_SCHEME, PDFStyle, 
     $     Chi2Style, LDebug, ifsttype, ASatur, LSatur, LFastAPPLGRID,
     $     Chi2MaxError, EWFIT, iDH_MOD, H1qcdfunc, CachePDFs


C Output style namelist
      namelist/Output/DoBands, Q2VAL, OutNX, OutXRange

C (Optional) MC method namelist
      namelist/MCErrors/LRand, ISeeDMC, StaType, SysType, LRandData
 
C (Optional) Chebyshev namelist
      namelist/Cheb/ILENPDF,pdfLenWeight,NCHEBGLU,NCHEBSEA
     $     ,IOFFSETCHEBSEA,ichebtypeGlu,ichebtypeSea
     $     ,WMNlen,WMXlen, ChebXMin

C (Optional) Polynomial parameterisation for valence
      namelist/Poly/NPOLYVAL,IZPOPOLY,IPOLYSQR
   
C (Optional) Data-set dependent scales
      namelist/Scales/DataSetMuR,DataSetMuF,DataSetIOrder

      character*128  LHAPDFSET
      integer ILHAPDFSET
      logical lhapdffile_exists

C (Optional) LHAPDF steering card
      namelist/lhapdf/LHAPDFSET,ILHAPDFSET

C (Optional) NNPDF steering card
      namelist/nnpdf/FLAGNNPDF,NNPDFSET,NNPDFRWDATA
     $     ,NNPDFREWEIGHTMETHOD,DONNPDFONLY,NNPDFOUTREPLICAS


      integer i, ilastq2
      integer StdCin,StdCout

C (Optional) set HQ scale
      namelist/HQScale/aq2,bq2,MassHQ


C Namelist for datafiles to read
      namelist/InFiles/NInputFiles,InputFileNames



C Namelist for EW parameters:
      namelist/EWpars/alphaem, gf, sin2thw, convfac,
     $ Mz, Mw, Mh, wz, ww, wh, wtp,
     $ Vud, Vus, Vub, Vcd, Vcs, Vcb,
     $ men, mel, mmn, mmo, mtn, mta, mup, mdn,
     $ mch, mst, mtp, mbt
C-----------------------------------------------------


C---------

*     ------------------------------------------------
*     Initialise basic parameters
*     ------------------------------------------------



      Itheory = 0
      EWFIT=0

      iDH_MOD = 0  ! no Dieter Heidt modifications to stat. errros.

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

C scale for HQ
      aq2 = 1
      bq2 = 0
      MassHQ = 'mc'

C  Hermes-like strange (off by default):
      ifsttype = 0

C  Cache PDF calls
      CachePDFs     = .false.

C  Fast applgrid:     
      LFastAPPLGRID = .false.
* 
      PDFStyle  = '13p HERAPDF'
      Chi2Style = 'HERAPDF'


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


C NNPDF defaults

      FLAGNNPDF = .false.
      NNPDFREWEIGHTMETHOD = 1
      DONNPDFONLY = .false.
      NNPDFSET = ''
      NNPDFOUTREPLICAS = 0


C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      strange_frac = 0.31
      charm_frac = 0.00


      mch        = 1.4D0
      mbt        = 4.75D0
      mtp        = 174.D0

      HFSCHEME = 0


C==== 24/08/2010: Add Saturation inspired cut ====
      ASatur = 0.0
      LSatur = 0.0

      Debug = .false.

C=================================================



C
C  Read the main H1Fitter namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=H1Fitter,END=41,ERR=42)
      close (51)

C     set debug flag used elsewhere according to steering
      Debug = lDebug

      open (51,file='ewparam.txt',status='old')
      read (51,NML=EWpars,END=43,ERR=44)
      close (51)

C
C  Read the output namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=Output,END=51,ERR=52)
      close (51)

C
C  Read the lhapdf namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=lhapdf,ERR=67,end=68)
 68   continue
      close (51)

C  Read the nnpdf namelist:
C 
      open (51,file='steering.txt',status='old')
      read (51,NML=nnpdf,END=75,ERR=74)
 75   continue
      close (51)

C check whether NNPDF and LHAPDF set are equal

      if (FLAGNNPDF) then
         if (TRIM(NNPDFSET) .ne. TRIM(LHAPDFSET)) then
            call HF_ErrLog(12032302,'W:WARNING: Setting LHAPDF set to '
     $           //TRIM(NNPDFSET))
            LHAPDFSET=NNPDFSET
         endif
C  check if the PDFstyle is indeed Ok
         if (PDFStyle.ne.'LHAPDF' .and. PDFStyle.ne.'LHAPDFQ0') then
            call HF_Errlog(12032303,
     $           'W:WARNING: Setting PDF style to LHAPDFQ0')
            PDFStyle = 'LHAPDFQ0'
         endif
      endif


C
C Decode HFSCHEME:
C      
      call SetHFSCHEME(HF_SCHEME)

C
C Decode PDF style:
C      
      call SetPDFStyle(PDFStyle)
      if ((PDFStyle.eq.'LHAPDF').or.(PDFStyle.eq.'LHAPDFQ0')) then
         INQUIRE(FILE=LHAPDFSET, EXIST=lhapdffile_exists) 
         if(lhapdffile_exists) then
            call InitPDFset(LHAPDFSET)
         else
            call InitPDFsetByName(LHAPDFSET)
         endif
         call InitPDF(ILHAPDFSET)
         if(PDFStyle.eq.'LHAPDF') then
            IPDFSET = 5
         endif
      endif

C
C Decode Chi2 style:
C
      call SetChi2Style(Chi2Style)

      HF_MASS(1) = mch
      HF_MASS(2) = mbt
      HF_MASS(3) = mtp


C
C  Read the MC method namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=MCErrors,ERR=62,end=61)
 61   continue
      close (51)

C
C  Read the Chebyshev namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=Cheb,ERR=64,end=63)
 
 63   continue
      close (51)


C
C  Read the Poly namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=Poly,ERR=66,end=65)
    
 65   continue
      close (51)

  

C
C  Read the HQScale namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=HQScale,ERR=70,end=69)
 69   continue
      close (51)

C
C  Read the data namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=InFiles,END=71,ERR=72)
      print '(''Read '',I4,'' data files'')',NInputFiles
      close (51)

C
C  Data-set dependent scales. First set defaults
C
      do i=1,NInputFiles
         DataSetMuR(i)    = 1.0D0
         DataSetMuF(i)    = 1.0D0
         DataSetIOrder(i) = IOrder
      enddo
C
C  Read the scales namelist:
C

      open (51,file='steering.txt',status='old')
      read (51,NML=Scales,END=123,ERR=124)
 123  Continue
      close (51)


C
C asign mc or mb to hq scale
C      
      call SetMHSCALE(MassHQ)


C
C Also read extra minuit parameters:
C      
      call readextraparam

      if (lDebug) then
C Print the namelists:
         print *,'Input Namelists:'
         print H1Fitter
         print Output
         print MCErrors
         print InFiles
         print EWpars
         print HQScale
      endif


      I_FIT_ORDER = IOrder
      starting_scale = Q02

C 07/12/2011 Dipole ===>
      call SetDipoleType
C <==

C
C Names of syst. error sources:
C
      do i=1,NSYS
         System(i) = ' '
      enddo

      goto 73
C 
 41   continue
      print '(''Namelist &H1Fitter NOT found'')'
      goto 73
 42   continue
      print '(''Error reading namelist &H1Fitter, STOP'')'
      call HF_stop
 43   continue
      print '(''Namelist @EWPars NOT found, STOP'')'
      call HF_stop
 44   continue
      print '(''Error reading namelist @EWPars, STOP'')'
      call HF_stop
 51   continue
      print '(''Namelist &Output NOT found'')'
      goto 73
 52   continue
      print '(''Error reading namelist &Output, STOP'')'
      call HF_stop
 62   continue
      print '(''Error reading namelist &MCErrors, STOP'')'
      call HF_stop
 64   continue
      print '(''Error reading namelist &Cheb, STOP'')'
      call HF_stop
 66   continue
      print '(''Error reading namelist &Poly, STOP'')'
      call HF_stop
 67   continue
      print '(''Error reading namelist &lhapdf, STOP'')'
      call HF_stop

 74   continue
      print '(''Error reading namelist &nnpf, STOP'')'
      call HF_stop

 70   continue
      print '(''Error reading namelist &HQScale, STOP'')'
      call HF_stop
 71   continue
      print '(''Namelist &InFiles NOT found'')'
      goto 73
 72   continue
      print '(''Error reading namelist &InFiles, STOP'')'
      call HF_stop
 124  print '(''Error reading namelist &scales, STOP'')'
      Call HF_stop


 73   continue

      chebxminlog = log(chebxmin)


      ilastq2 = NBANDS/2
      do i=NBANDS/2,1,-1
         if (q2val(i).lt.0) ilastq2 = i-1
      enddo
      do i=1,NBANDS/2
       Q2VAL(i+ilastq2) = Q2VAL(i+NBANDS/2)
      enddo
c      print *,'q2val ', (q2val(i),i=1,NBANDS)

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


      if (NCHEBGLU.ne.0) then
         print *,'Use Chebyshev polynoms for gluon with N=',NCHEBGLU
      endif

      if (NCHEBSEA.ne.0) then
         print *,'Use Chebyshev polynoms for sea with N=',NCHEBSEA
         print *,'Offset for minuit parameters is',IOFFSETCHEBSEA
      endif

      return
      end


      Subroutine SetPDFStyle(PDFStyle)
C---------------------------------------
C
C>  Set PDF parameterisation type
C
C---------------------------------------
      implicit none
      character*(*) PDFStyle
      character*32  LHAPDFSET
      integer ILHAPDFSET 
      include 'steering.inc'
C---------------------------------

      FlexibleGluon = .false.
      
      if (PDFStyle.eq.'10p HERAPDF') then
         iparam = 22
         PDF_DECOMPOSITION = 'Dv_Uv_Dbar_Ubar'

      elseif (PDFStyle.eq.'13p HERAPDF') then
         iparam = 229
         FlexibleGluon = .true.
         PDF_DECOMPOSITION = 'Dv_Uv_Dbar_Ubar'

      elseif (PDFStyle.eq.'strange') then
         iparam = 2011
         FlexibleGluon = .true.
         PDF_DECOMPOSITION = 'Dv_Uv_Dbar_Ubar_Str'

      elseif (PDFStyle.eq.'CTEQ') then
         iparam = 171717
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

      else
         print *,'Unsupported PDFStyle =',PDFStyle
         print *,'Check value in steering.txt'
         call HF_stop
      endif

      end



      Subroutine SetHFSCHEME(HF_SCHEME)
C---------------------------------------
C
C>  Set PDF parameterisation type
C
C---------------------------------------
      implicit none
      character*(*) HF_SCHEME
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
      elseif (HF_SCHEME.eq.'FF') then
          HFSCHEME = 3 
      elseif (HF_SCHEME.eq.'ABKM FFNS') then		
         HFSCHEME = 4 		
      elseif (HF_SCHEME.eq.'ABKM BMSN') then		
         HFSCHEME = 44 
      else
         print *,'Unsupported HFSCHEME =',HF_SCHEME
         print *,'Check value in steering.txt'
         call HF_stop
      endif

      end



      Subroutine SetMHSCALE(MassHQ)
C---------------------------------------
C
C>  Set HQ scale mh parameter
C
C---------------------------------------
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


      Subroutine SetChi2Style(Chi2Style)
C---------------------------------------
C
C>  Set Chi2 style
C
C---------------------------------------
      implicit none
      character*(*) Chi2Style
      include 'steering.inc'
C---------------------------------

      if (Chi2Style.eq.'HERAPDF') then
         ICHI2 = 11
      elseif (Chi2Style.eq.'HERAPDF Sqrt') then
         ICHI2 = 31
      elseif (Chi2Style.eq.'HERAPDF Linear') then
         ICHI2 = 21
      elseif (Chi2Style.eq.'CTEQ') then
         ICHI2 = 2
      elseif (Chi2Style.eq.'H12000') then
         ICHI2 = 1        
      elseif (Chi2Style.eq.'H12011') then
         ICHI2 = 41        
      elseif (Chi2Style.eq.'Offset') then
         ICHI2 = 3
      else
         print *,'Unsupported Chi2Style =',Chi2Style
         print *,'Check value in steering.txt'
         call HF_stop
      endif
      
      end


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
      
      nExtraParam = 0

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
C Add extra param
               nExtraParam = nExtraParam + 1
               if (nExtraParam.gt. nExtraParamMax) then
                  print *,'Number of extra parameters exceeds the limit'
                  print *,'nExtraParam=',nExtraParam
                  print *,'Check your steering, stopping'
                  call HF_stop
               endif
               ExtraParamNames(nExtraParam) = name(i)
               ExtraParamValue(nExtraParam) = value(i)
               ExtraParamStep (nExtraParam) = step(i)
               ExtraParamMin  (nExtraParam) = min(i)
               ExtraParamMax  (nExtraParam) = max(i)
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
