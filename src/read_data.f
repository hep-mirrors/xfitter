      subroutine read_data

*     ------------------------------------------------
*     Read data     
*     ------------------------------------------------



      implicit none
*     ------------------------------------------------
      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'      
      include 'systematics.inc'
      include 'indata.inc'
      include 'couplings.inc'
      include 'for_debug.inc'

      character*10  cdummy
      character   adum 
      character*1   a1dum 

      double precision ECM     
      double precision fac

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

*     ------------------------------------------------
*     initilialising
*     ------------------------------------------------

      NSYS = 0
      DEBUG  = lDEBUG
      GNORM = .true.



      do i=1,nsysMax
         do j=1,ntot
            BETA(i,j) = 0.d0
         enddo
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

C-----------------------------------------
      print*,'number of points', npoints

      do i=1,nsysMax
         do k=1,npoints
            beta(i,k) = beta(i,k) / 100.
         enddo
      enddo

      do i=1,npoints
         ALPHA(i) = ALPHA(i) / 100.
         if (alpha(i).le.0) write(6,*) 'alpha(i) = 0 for point ',i
      enddo


*     ----------------------------------------------------------------------
*     -- CTEQ-like chi2 :
*     -- compute the matrix (A) as given in Eq. (B.4) of JHEP07 (2002) 012
*     -- This is the matrix "sysa" in common/systema/

      if (ICHI2.eq.2) then
         call Systematics
      endif

      if (LDebug) then
C
C Dump beta,alpha matricies
C
         open (61,file='beta.dat',status='unknown')
         do k=1,npoints
            write (61,'(I5,500(F6.2))') k,(Beta(i,k)*100.0,i=1,NSYSMAX)
         enddo
         close(61)
      endif

C
C MC method: fluctuate data according to their uncertainteis.
C
      if (lrand) then
         call MC_method
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
      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'systematics.inc'

      character *(*) CFile
C Namelist  variables:    
      integer ndataMax,ninfomax,nsystMax,ncolumnMax
      parameter (ndataMax=1000)
      parameter (ninfoMax=100)
      parameter (nsystMax=500)

      parameter (ncolumnMax = nsystMax+NBinDimensionMax+1)

      character *80 Name
      integer  NData
      integer  NUncert
      integer  NInfo
      integer  NBinDimension

      
      character *80 BinName(NBinDimensionMax)
      double precision datainfo(ninfoMax)
      character *80 CInfo(ninfoMax)
      character *80 Reaction

      double precision buffer(ncolumnMax)

C
C Name and type of columns:
C      
      integer   NColumn 
      character *32 ColumnName(ncolumnMax)
      character *32 ColumnType(ncolumnMax)

C Systematics:
      character *32 SystematicType(nsystMax)
      logical Percent(1:nsystMax)

C Reference table
      integer CompressIdx(nsystMax)

      integer IndexDataset
      double precision SystScales(nsystMax)
C Extra info about k-factors, applegrid file(s):
      character*80 TheoryInfoFile,TheoryType
      character*80 KFactorNames(NKFactMax)
      integer      NKFactor
C Namelist definition:
      namelist/Data/Name,NData
     $     ,NInfo,datainfo,CInfo,Reaction,Percent
     $     ,SystScales, IndexDataset
     $     ,TheoryInfoFile,TheoryType,KFactorNames,NKFactor
     $     ,ColumnName, ColumnType, NColumn

      double precision XSections(ndataMax)
      double precision AllBins(10,ndataMax)
      double precision Syst(nsystmax)

      double precision Akfact(NKFactMax)

      double precision StatError   ! stat
      double precision UncorError  ! uncorrelated systematics
      double precision TotalError  ! total uncertainty

      double precision TotalErrorRead ! total error, provided by the data file

      integer idxSigma

      integer i,j,iBin,iError
      logical LReadKFactor

C Temporary buffer to read the data (allows for comments starting with *)
      character *4096 CTmp

C Function to check cuts
      logical FailSelectionCuts
      
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

C Reset scales to 1.0
      do i=1,nsysmax
         SystScales(i) = 1.0
         ColumnType(i) = ' '
         ColumnName(i) = ' '
      enddo

      open(51,file=CFile,status='old',err=99)

      print *,'Reading data file ...'
      print *,CFile
      read(51,NML=Data,err=98)

C
C Check dimensions
C
      if (NColumn.gt. Ncolumnmax) then
         print '(''Error in ReadDataFile for File='',A80)',cfile
         print '(''NColumn = '',i6,'' exeeds NColumnMax='',i6)',ncolumn,ncolumnmax
         stop
      endif
C
C Store 
C
      NDATASETS = NDATASETS + 1
      DATASETNUMBER(NDATASETS)   = 10000+NDATASETS
      DATASETLABEL(NDATASETS)    = Name
      DATASETNUMBER(NDATASETS)   = IndexDataset   !!!  

C Reaction info:
      DATASETREACTION(NDATASETS) = Reaction

C Parse ColumnType, count systematics, etc
      do i=1,NColumn
         if (ColumnType(i).eq.'Bin') then
            NBinDimension = NBinDimension + 1
            BinName(NBinDimension) = ColumnName(i)
         elseif (ColumnType(i).eq.'Sigma') then
            idxSigma = i
         elseif (ColumnType(i).eq.'Error') then
            NUncert = NUncert + 1
            ! Special case: uncorrelated errors
            if (index(ColumnName(i),'uncor').gt.0) then
               SystematicType(NUncert) = 'uncor'
            else
               SystematicType(NUncert) = ColumnName(i)
            endif
         elseif (ColumnType(i).eq.'Dummy') then
! Ignore dummy column
         else
            print '(''Unknown Column type for dataset'',A80)',CFile
            print '(''Column='',i5,'' type='',A32)',i,ColumnType(i)
            print '(''STOP in ReadDataFile'')'
            stop
         endif
      enddo

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

C Prepare systematics:
      do i=1,NUncert
C--- Statistical: special case
         if (SystematicType(i).eq.'stat') then
C--- Uncorrelated: special case
         else if (SystematicType(i).eq.'uncor') then
C--- Total error: special case
         else if (SystematicType(i).eq.'total') then
C--- Ignore: special case
         else if (SystematicType(i).eq.'ignore') then

         else
C--- Check if the source already exists:         
            do j=1,NSYS            
               if ( system(j).eq.SystematicType(i) ) then
                  CompressIdx(i) = j
                  goto 80
               endif
            enddo
C--- Add new source
            NSYS = NSYS + 1
            if (NSYS.gt.NSysMax) then
               print 
     $              '(''ReadDataFile Error: exeeding NSysMax'')'
               print '(''Current NSysMax='',i6)',NSysMax
               print '(''Increase NSysMax in systematics.inc'')'
               stop
            endif

            System(NSYS) = SystematicType(i)
            CompressIdx(i) = NSYS
 80         continue
         endif
      enddo

C Theory file if present:
      DATASETTheoryFile(NDATASETS) = TheoryInfoFile
      DATASETTheoryType(NDATASETS) = TheoryType

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
            print *,'STOP'
            stop
         endif
         do i=1,NColumn
            if (ColumnName(i) .eq. ' ') then
               print *,'Undefined ColumnName !!!'
               print *,'Check name for column number = ',i
               print *,'STOP'
               stop
            endif
         enddo

C Read the colums
         read (ctmp,*,err=1019)(buffer(i),i=1,NColumn)

C Decode the columns
         iBin   = 0
         iError = 0
         do i=1,NColumn
            if (ColumnType(i).eq.'Bin') then
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
            read (52,*) (akfact(i),i=1,NKFactor)
         endif

C Apply cuts:
         if (FailSelectionCuts(Reaction,NBinDimension,allbins(1,j),BinName)) then
            goto 1717
         endif

C Add a point:
         npoints = npoints+1
         
         if (npoints.ge.NTOT) then
            print 
     $ '('' ReadDataFile Error, increase NTOT value inside ntot.inc'')'
            stop
         endif

         NDATAPOINTS(NDATASETS) = NDATAPOINTS(NDATASETS) + 1
         DATASETIDX(NDATASETS,NDATAPOINTS(NDATASETS)) = npoints

C Translate errors in %:
         TotalError = 0.
         UncorError = 0.
         StatError = 0.
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
            if (SystematicType(i).eq.'uncor') then
C Uncorrelated error:
               UncorError = UncorError +  Syst(i)**2
            endif
            if (SystematicType(i).eq.'stat') then
C Stat error:
               StatError = StatError +  Syst(i)**2
            endif

         enddo

         StatError = sqrt(StatError)
         UncorError = sqrt(UncorError)
         TotalError = sqrt(TotalError)

         DATEN(npoints) = XSections(j)
         E_UNC(npoints)  = UncorError
         E_TOT(npoints)  = TotalError
         E_STA(npoints)  = StatError

         ! > Check total error
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

         ALPHA(npoints) = sqrt(UncorError**2+StatError**2)*DATEN(npoints)

         do i=1,NUncert
            if (SystematicType(i).ne.'uncor' .and. 
     $           SystematicType(i).ne.'ignore'.and.
     $           SystematicType(i).ne.'stat'
     $           ) then
               print *,CompressIdx(i),npoints,SystematicType(i)
               BETA(CompressIdx(i),npoints) = syst(i)
            endif
         enddo


         JSET(npoints) = IndexDataset  

C Store k-factors:
         if (lreadkfactor) then
            do i=1,nkfactor
               kfactors(i,npoints) = akfact(i)
            enddo
         endif
         
 1717  continue
      enddo

      close (51)
      if (lreadkfactor) then
         close (52)
      endif

      print '(''Read'',i8,'' data points for '',A80)',NData,Name
      print '(''Printing first'',i5,'' data points'')',min(Ndata,5)
      print '(20A14)',(BinName(i),i=1,NBinDimension),' sigma'
    
      do j=1,min(NData,5)
         print '(20E14.4)',(Allbins(i,j),i=1,NBinDimension),XSections(j)
    
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

 1017 continue
      print '(''Error reading file'')'
      stop
 1018 continue
      print '(''End of file while reading file'')'
      stop
 1019 continue
      print '(''Problem interpreting data line='',i6)',j
      stop

      end

