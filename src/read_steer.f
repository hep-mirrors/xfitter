
      subroutine read_steer
      implicit none


      include 'steering.inc'
      include 'couplings.inc'
      include 'pdflength.inc'
      include 'pdfparam.inc'
      include 'datasets.inc'

C==== 26/07/2010: ADDED FOR DIPOLE MODEL =========
      INCLUDE 'dip.inc'
C=================================================

      integer i, ilastq2
      integer StdCin,StdCout

      real*4 SPACE
      common/CFREAD/SPACE(6000)
      integer Isch, Iset, Iflg, Ihad
      Common /Ischeme/ Isch, Iset, Iflg, Ihad  !*** pass info out to ACOT

C SG: add namelist for datafiles to read
      namelist/InFiles/NInputFiles,InputFileNames
C---------

*     ------------------------------------------------
*     Initialise basic parameters
*     ------------------------------------------------

      CALL FFINIT(6000)
      StdCin  = 5
      StdCout = 6

      OPEN(UNIT=StdCin,FILE='steering.txt',STATUS='unknown')
      CALL FFSET ( 'LINP' , StdCin )
      CALL FFSET ( 'LOUT' , StdCout)
      CALL FFSET ( 'SIZE' , 10)


      Itheory = 0
      CALL FFKEY('itheory',ITHEORY,1,'INTE')

      ISEED = 2313134
      CALL FFKEY('ISEED',Iseed,1,'INTE')

      ISDRN = 42
      Call FFKEY('ISDRN',IsdRN,1,'INTE')

      I_FIT_ORDER = 2
      CALL FFKEY('IORDER',I_FIT_ORDER,1,'INTE')

C=================================================


C SG: Key for PDF length:
      ILENPDF = 0
      Call FFKEY('ILENPDF',ILENPDF,1,'INTE')

C SG: Key for PDF length weight factor:
      do i=1,5
         pdfLenWeight(i) = 0.
      enddo
      Call FFKEY('PDFLENWEIGHT',pdfLenWeight,5,'REAL')

C SG: Key for Chebyshev param. of the gluon:
C
      NCHEBGLU = 0
      Call FFKEY('NCHGLU',NCHEBGLU,1,'INTE')

C SG: Key for Chebyshev param. of the Sea:
      NCHEBSEA = 0
      Call FFKEY('NCHSEA',NCHEBSEA,1,'INTE')

C 2 Feb 2010 SG: Offset for the Sea chebyshev parameters (default:20)
      IOFFSETCHEBSEA = 20
      Call FFKEY('IOFS',IOFFSETCHEBSEA,1,'INTE')

C SG: Type of Chebyshev parameterization:
      ichebtypeGlu = 0
      Call FFKEY('CHTGLU',ichebtypeGlu,1,'INTE')
      ichebtypeSea = 0
      Call FFKEY('CHTSEA',ichebtypeSea,1,'INTE')


C 25 Jan 2011, SG
C SG: Pure polinomial param for the valence quarks:
C     by default starting from N=61 for Uv and N=71 for Dv
      NPOLYVAL = 0
      Call FFKEY('NPVA',NPOLYVAL,1,'INTE')

C Add option to change Z of valence PDFs at x=1 from  (1-x) to (1-x)^2
      IZPOPOLY = 1
      Call FFKEY('IZPO',IZPOPOLY,1,'INTE')

C Square polynom before calculating dv,uv. This forces positivity
      IPOLYSQR = 0
      Call FFKEY('IPSQ',IPOLYSQR,1,'INTE')

C SG: Key for W range 
      WMNlen =  20.
      WMXlen = 320.
      Call FFKEY('WMNLEN',WMNlen,1,'REAL')
      Call FFKEY('WMXLEN',WMXlen,1,'REAL')

      IPARAM = 1
      CALL FFKEY('IPAR',IPARAM,1,'INTE')

      chebxmin = 1.E-5
      Call FFKEY('CHEBXMIN',chebxmin,1,'REAL')

C SG: Hermes-like strange:
      ifsttype = 0
      Call FFKEY('IFSTTYPE',ifsttype,1,'INTE')

*     ICHI2 = 1 : Pascaud-like, +10/20/30 for scaled error variants
*     ICHI2 = 2 : CTEQ-like
*     ICHI2 = 3 : use full covariant matrix
      ICHI2 = 1
      CALL FFKEY('ICHI2',ICHI2,1,'INTE')

      lfirst = .true. 
      CALL FFKEY('FIRST',lfirst,1,'LOGICAL')

      lcorr = .false.
      CALL FFKEY('CORR',lcorr,1,'LOGICAL')

      lCORWEAK = .false.
      CALL FFKEY('CORWEAK',lCORWEAK,1,'LOGICAL')

      lONLINE = .true.
      CALL FFKEY('ONLINE',lONLINE,1,'LOGICAL')

      lDEBUG = .false.
      CALL FFKEY('DEBUG',lDEBUG,1,'LOGICAL')

      lNORMA = .false.
      CALL FFKEY('NORMA',lNORMA,1,'LOGICAL')

      lTHEO = .false.
      CALL FFKEY('THEO',lTHEO,1,'LOGICAL')

      lRAND = .false.
      CALL FFKEY('RAND',lRAND,1,'LOGICAL')

      iSEEDmc = 0
      Call FFKEY('SEED',iSeeDmc,1,'INTE')

      STATYPE = 1
      CALL FFKEY('STATYPE',STATYPE,1,'INTE')

      SYSTYPE = 1
      CALL FFKEY('SYSTYPE',SYSTYPE,1,'INTE')

      DOBANDS = .false.
      call ffkey('BANDS',DOBANDS,1,'LOGICAL')

      do i=1,NBANDS
       Q2VAL(i) = -1.
      enddo
      call ffkey('Q2VAL',Q2VAL,NBANDS/2,'REAL')
      call ffkey('Q3VAL',Q2VAL(21),NBANDS/2,'REAL')
         
      outform = 0
      CALL FFKEY('OUTFORM',outform,1,'INTE')
      outnx = 101
      CALL FFKEY('OUTNX',outnx,1,'INTE')
      outxrange(1) = 1e-4
      outxrange(2) = 1.0
      CALL FFKEY('OUTXRANGE',outxrange,2,'REAL')

      starting_scale = 4.
      call ffkey('Q02',starting_scale,1,'REAL')

      strange_frac = 0.33
      call ffkey('FSTRANGE',strange_frac,1,'REAL')

      charm_frac = 0.15
      call ffkey('FCHARM',charm_frac,1,'REAL')

      rtalphas = 0.1185
      call ffkey('RTALPHAS',rtalphas,1,'REAL')

      HF_MASS(1) = 1.4
      HF_MASS(2) = 4.5
      call ffkey('HFMAS',HF_MASS,2,'REAL') 

      IFUDGEFL = 0
      CALL FFKEY('FAFL',IFUDGEFL,1,'INTE')
	
      IFUDGEF2 = 0
      CALL FFKEY('FAF2',IFUDGEF2,1,'INTE')

      HFSCHEME = 0
      call ffkey('HFSCHEME',HFSCHEME,1,'INTE')

      vfnsINDX =0 
      call ffkey('vfnsINDX',vfnsINDX,1,'INTE')

      ISCH = 1
      call ffkey('ISCH',ISCH,1,'INTE')

      HF_THRE(1) = 1.4
      HF_THRE(2) = 4.5
      call ffkey('HFTHRES',HF_THRE,2,'REAL')

      pxmax = -1.
      CALL FFKEY('xmax',pxmax,1,'REAL')
      pxmin = -1.
      CALL FFKEY('xmin',pxmin,1,'REAL') 

      pq2min = -1.
      CALL FFKEY('q2min',pq2min,1,'REAL')

C==== 26/07/2010: ADDED Q2MAX CUT ================
      pq2max = -1.
      CALL FFKEY('q2max',pq2max,1,'REAL')


C==== 24/08/2010: Add Saturation inspired cut ====
      ASatur = 0.0
      Call FFKEY('ASat',ASatur,1,'REAL')
      LSatur = 0.0
      Call FFKEY('LSat',LSatur,1,'REAL')


C=================================================


      CALL FFGO

C
C SG 25/05/11
C
C  Read the data namelist:
C
      open (51,file='steering.txt',status='old')
      read (51,NML=InFiles,END=71,ERR=72)
      print '(''Read '',I4,'' data files'')',NInputFiles
      close (51)
      goto 73
C 
 71   continue
      print '(''Namelist @InFiles NOT found'')'
      goto 73
 72   continue
      print '(''Error reading namelist @InFiles, STOP'')'
      stop
 73   continue

      chebxminlog = log(chebxmin)


      ilastq2 = NBANDS/2
      do i=NBANDS/2,1,-1
         if (q2val(i).lt.0) ilastq2 = i-1
      enddo
      do i=1,NBANDS/2
       Q2VAL(i+ilastq2) = Q2VAL(i+NBANDS/2)
      enddo
      print *,'q2val ', (q2val(i),i=1,NBANDS)

* --- Check the consistency of the steering file

      if (HFSCHEME.eq.1.and.HF_THRE(2)**2.lt.starting_scale) then
       write(6,*)
       write(6,*) 'Bottom thres. has to be larger than starting scale'
       write(6,*)
       stop
      endif

      if (HFSCHEME.eq.1.and.HF_THRE(2).lt.HF_THRE(1)) then
       write(6,*)
       write(6,*) 'Bottom thres. has to be larger than charm thres.'
       write(6,*)
       stop
      endif

      if (HFSCHEME.eq.1.and.(LDOFIT1(4).eq.1.or.
     +                       LDOFIT1(6).eq.1.or.
     +                       LDOFIT1(9).eq.1)) then
       write(6,*)
       write(6,*) 'Do not include CC in massive scheme'
       write(6,*)
       stop
      endif

      if (NCHEBGLU.ne.0) then
         print *,'Use Chebyshev polynoms for gluon with N=',NCHEBGLU
      endif


      if (IFUDGEF2.ne.0) then
         print *,'Use Fudge for F2 in the HFSCHEME=2 (HT effects)',IFUDGEF2
      endif
      if (IFUDGEFL.ne.0) then
         print *,'Use Fudge for FL in the HFSCHEME=2 ', IFUDGEFL
      endif



      if (NCHEBSEA.ne.0) then
         print *,'Use Chebyshev polynoms for sea with N=',NCHEBSEA
         print *,'Offset for minuit parameters is',IOFFSETCHEBSEA
      endif

      if ( napplgrids .ne. 0 ) useapplg = .true.
      if ( useapplg ) then
        call getAPPLgrids(napplgrids)
      endif

      return
      end

