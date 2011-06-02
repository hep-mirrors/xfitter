      subroutine read_data
      
c     
c     Read data
c     


      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'

      
      include 'ntot.inc'
      include 'systematics.inc'
      include 'hadcor.inc'
      include 'indata.inc'
      include 'couplings.inc'


      character*10  cdummy
      character   adum 
      character*1   a1dum 

      double precision ECM
c     

      double precision fac

      logical FIRST             !  true : cov matrix recalculated
      logical GNORM             !  correlated part for the luminosity errors




      

cv ///// study of PDF uncertainties
cv
      double precision alnorm
      external alnorm
      real logshift
      external logshift
      real lsig, lmu, lrunif
      real dummy, dummy_st
      integer vi,icount
      double precision  voica, voica_fl, voica_un
      real ranmflat
      double precision rand_shift(NSYS)
      double precision r_sh_fl(NSYS)
      integer numsys
      real rndsh, ranflat
      integer num,iseedrand, idate,is,ntime, ndate

      real f_un
      COMMON/SLATE/IS(40)
      real alumlognorm(300)
      real alumierr(300)
      data alumierr/300*0.0/
cv////////
 

      integer i,j,k,iset,n0,isys,iq2bin,iebin,jsys

      ONLINE = lONLINE
      FIRST = lFIRST

      GNORM = .true.

      ISCHARMANY  = .false.
      ISBOTTOMANY = .false.


      do iset=1,nset
         POLAR(iset) = 0.d0
      enddo

      do iset=1,nset
         ISJET(iset) = .false.
         ISCHARM(iset) = .false.
         ISBOTTOM(iset) = .false.
         ISDYZ(iset) = .false.
         ISDYW(iset) = .false.
      enddo

      do i=1,nsys
         do j=1,ntot
            BETA(i,j) = 0.d0
         enddo
      enddo



      do i=1,nset
         NDATAPOINTS(i) = 0
      enddo
      npoints = 0

cv============
cv RANDOM SHIFTS


      if (lRAND) then
         f_un = 2.
         icount= time()
         print*,' clock = ', icount
         
         call datime(ndate,ntime)
         ntime = ntime*1000000+is(6)
         icount=ntime

cv initialise the random shifts
         do numsys=1,nsys
            rand_shift(numsys)=0.
            r_sh_fl(numsys)=0.
         enddo

cv initialise the random seed gener
         if (iseedmc.ne.0) then
C SG: Overwrite seed to one selected in the steering:
            icount = iseedmc
         endif

         call rmarin(icount,0,0)
         call rluxgo(3,icount,0,0)
         print*,'initialize smeering with a seed isdrn = ',icount
         
c         stop
         do numsys=1,nsys
            call rnorml(rndsh,1)    ! gauss random number
            call ranlux(ranflat,1)   ! uniform random number

            rand_shift(numsys) = rndsh
            r_sh_fl(numsys) = ranflat

            print*,'random numbers: sys, gauss, flat ',
     $           numsys,rand_shift(numsys),
     $           r_sh_fl(numsys)

cv save shifts for lumi uncertainties only (but we are not using it)            
            if ((numsys.eq.1).or.(numsys.eq.7).or. (numsys.eq.14)
     $           .or.(numsys.eq.18).or.(numsys.eq.19)
     $           .or.(numsys.eq.25).or.(numsys.eq.26).or.(numsys.eq.27)
     $           .or.(numsys.eq.32).or.(numsys.eq.8)) then
               alumlognorm(numsys) = alnorm(1.,alumierr(numsys))
            else
               
            endif
         enddo

      endif

cv ============done with initializing=========


C
C Read data from namelists:
C
      do i=1,NInputFiles
         call ReadDataFile(InputFileNames(i))
      enddo

C-----------------------------------------
      print*,'number of points', npoints



C
C MINUIT PARAMETERS FOR NORMALISATIONS !!!!!!!
C


*     ------------------------------------------------------------------
*     -- with this the normalisations are minuit parameters instead...

      if (lNORMA) then

         do jsys=1,nsys
            do k=1,npoints
               if (jsys.eq.1.or.jsys.eq.7.or.jsys.eq.8.or.
     +              jsys.eq.14.or.jsys.eq.18.or.jsys.eq.19.or.
     +              jsys.eq.25.or.jsys.eq.26.or.jsys.eq.27) then
                  beta(jsys,k) = 0.d0
               endif
            enddo
         enddo

      endif
*     ------------------------------------------------------------------




      do i=1,nsys
         do k=1,npoints
            beta(i,k) = beta(i,k) / 100.
            if (.not.lCORR) beta(i,k)=0.d0
         enddo
      enddo


      do i=1,npoints
         ALPHA(i) = ALPHA(i) / 100.
         if (alpha(i).le.0) write(6,*) 'alpha(i) = 0 for point ',i
      enddo





*     ------------------------------------------------------------------
*     -- Calculate or read the full covariance matrix

      if (ICHI2.eq.3) then

         if (FIRST) then
            call CoVarMatrix
         else
            call Read_CoVarMatrix
         endif

      endif


*     ----------------------------------------------------------------------
*     -- CTEQ-like chi2 :
*     -- compute the matrix (A) as given in Eq. (B.4) of JHEP07 (2002) 012
*     -- This is the matrix "sysa" in common/systema/

      if (ICHI2.eq.2) then
         call Systematics
      endif



      return
      end



      subroutine ReadDataFile(CFile)
C------------------------------------------------------------------------
C
C  Created 20 May 2011 by SG.
C   Read data set using namelist format
C
C------------------------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'systematics.inc'

      character *(*) CFile
C Namelist  variables:    
      integer ndataMax,ninfomax,nsystMax
      parameter (ndataMax=1000)
      parameter (ninfoMax=100)
      parameter (nsystMax=500)

      character *80 Name
      integer  NData
      integer  NSyst
      integer  NInfo
      integer  NBinDimension

      
      character *80 BinName(NBinDimensionMax)
      double precision datainfo(ninfoMax)
      character *80 CInfo(ninfoMax)
      character *80 Reaction
      logical Percent(0:nsystMax)
      integer SystematicType(nsystMax)
      integer IndexDataset
      double precision SystScales(nsystMax)
C Extra info about k-factors, applegrid file(s):
      character*80 TheoryInfoFile,TheoryType
      character*80 KFactorNames(NKFactMax)
      integer      NKFactor
C Namelist definition:
      namelist/Data/Name,NData,NSyst,NBinDimension
     $     ,BinName,NInfo,datainfo,CInfo,Reaction,Percent
     $     ,SystematicType, SystScales, IndexDataset
     $     ,TheoryInfoFile,TheoryType,KFactorNames,NKFactor

      double precision XSections(ndataMax)
      double precision StatErrors(ndataMax)
      double precision AllBins(10,ndataMax)
      double precision Syst(nsystmax)

      double precision Akfact(NKFactMax)

      double precision StatError

      double precision UncorError  ! uncorrelated systematics
      double precision TotalError  ! total uncertainty

      integer i,j
      logical LReadKFactor

C Function to check cuts
      logical FailSelectionCuts
      
C-------------------------------------------------------      

C Reset to default:
      NSYST = 0
      NData = 0
      NInfo = 0
      NBinDimension = 0
      Reaction = ' '
      Name     = ' '
      IndexDataSet = 0
      NKFactor = 0
      TheoryInfoFile = ' '
      LReadKFactor = .false.

      open(51,file=CFile,status='old',err=99)

      print *,'Reading data file ...'
      print *,CFile
      read(51,NML=Data,err=98)
C
C Store 
C
      NDATASETS = NDATASETS + 1
      DATASETNUMBER(NDATASETS)   = 10000+NDATASETS
      DATASETLABEL(NDATASETS)    = Name
      DATASETNUMBER(NDATASETS)   = IndexDataset   !!!  
C Reaction info:
      DATASETREACTION(NDATASETS) = Reaction
C Binning info:
      DATASETBinningDimension(NDATASETS) = NBinDimension
      do i=1,NBinDimension
         DATASETBinNames(i,NDATASETS) = BinName(i)
      enddo
C Extra info:
      DATASETInfoDimension(NDATASETS) = NInfo
      do i=1,NInfo
         DATASETInfoNames(i,NDATASETS) = CInfo(i)
         DATASETInfo(i,NDATASETS) =      DataInfo(i)
      enddo
      NQ2BINS(NDATASETS) = 0            ! uff

C Check if we need to read kfactor file:
      if (TheoryInfoFile.ne.' ') then
         if (TheoryType.eq.'kfactor') then
            open (52,file=TheoryInfoFile,status='old',err=100)
            lreadkfactor = .true.
         endif
      endif
      DATASETNKfactors(NDATASETS) = NKFactor
      do i=1,NKFactor
         DATASETKFactorNames(i,NDATASETS) = KFactorNames(i)
      enddo


C Read data info:
      do j=1,NData
         read(51,*)(allbins(i,j),i=1,NBinDimension),XSections(j)
     $        ,StatErrors(j),(syst(i),i=1,NSyst)

         if (lreadkfactor) then
            read (52,*) (akfact(i),i=1,NKFactor)
         endif

C Apply cuts:
         if (FailSelectionCuts(Reaction,NBinDimension,allbins(1,j),BinName)) then
            goto 1717
         endif

C Add a point:
         npoints = npoints+1
         
         if (npoints.ge.NTOT) then
            print '('' ReadDataFile Error, increase NTOT value inside ntot.inc'')'
            stop
         endif

         NDATAPOINTS(NDATASETS) = NDATAPOINTS(NDATASETS) + 1
         DATASETIDX(NDATASETS,NDATAPOINTS(NDATASETS)) = npoints

C Translate errors in %:
         TotalError = 0.
         UncorError = 0.

         if (.not.Percent(0)) then
            StatErrors(j) = StatErrors(j)/XSections(j)*100.
         endif
         TotalError = TotalError + StatErrors(j)**2

         do i=1,NSyst
            if (.not.Percent(i)) then
               syst(i) = syst(i)/XSections(j)*100.
            endif
   
            TotalError = TotalError + Syst(i)**2
            if (SystematicType(i).eq.0) then
C Uncorrelated error:
               UncorError = UncorError +  Syst(i)**2
            endif
         enddo

         UncorError = sqrt(UncorError)
         TotalError = sqrt(TotalError)

         DATEN(npoints) = XSections(j)
         E_UNC(npoints)  = UncorError
         E_TOT(npoints)  = TotalError
         E_STA(npoints)  = StatErrors(j)

         do i=1,NBinDimension
            AbstractBins(i,npoints) = allbins(i,j)
         enddo

         ALPHA(npoints) = sqrt(UncorError**2+StatErrors(j)**2)*DATEN(npoints)
         do i=1,NSyst
            if (SystematicType(i).gt.0) then
               if (SystematicType(i).gt. NSys) then
                  print '(''ReadDataFile Error: requested error source'',i6,'' larger than NSYST='',i6)'
     $                 ,SystematicType(i),NSys
                  print '(''Check SystematicType or increase NSys in systematics.inc'')'
                  stop
               endif

               BETA(SystematicType(i),npoints) = syst(i)
            endif
         enddo


C         print *,'hhhh',alpha(npoints),e_sta(npoints)

         JSET(npoints) = IndexDataset  ! XXXXXXXXXXXX

C Store k-factors:
         if (lreadkfactor) then
            do i=1,nkfactor
               kfactors(i,npoints) = akfact(i)
            enddo
         endif
         
 1717 enddo

      close (51)
      if (lreadkfactor) then
         close (52)
      endif

      print '(''Read'',i8,'' data points for '',A80)',NData,Name
      print '(''Printing first'',i5,'' data points'')',min(Ndata,5)
      print '(20A14)',(BinName(i),i=1,NBinDimension),' sigma'
     $     ,' stat. err'
      do j=1,min(NData,5)
         print '(20E14.4)',(Allbins(i,j),i=1,NBinDimension),XSections(j)
     $        ,StatErrors(j)
      enddo
      return

 98   continue
      print '(''Error reading namelist Data'')'
      print *,CFile
      stop

 99   continue
      print '(''Can not open file '')'
      print *,CFile
      stop
100   continue
      print '(''Can not open file '')'
      print *,TheoryInfoFile
      stop
      end


