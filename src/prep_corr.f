      subroutine prep_corr
*     ------------------------------------------------
*     Prepare statistical correlation matrix corr_stat
*     and given covariance matrix cov
*     ------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'covar.inc'


      integer k,i,j,b,m,n,NCorr,NIdColumns1,NIdColumns2,NIdMax,NCorrMax

      parameter (NIdMax = 3)
      character *80 Name1
      character *80 Name2
      character *80 MatrixType
      character *80 IdColumns1(NIdMax)
      character *80 IdColumns2(NIdMax)


C     Temporary buffer to read the data (allows for comments starting with *)
      character *4096 CTmp
      double precision buffer(2*NIdMax+1)
      integer IdIdx1(NIdMax)
      integer IdIdx2(NIdMax)

c     matrix buffer
      parameter (NCorrMax = 100)
      double precision matrixbuffer(NCorrMax,NCorrMax)
      double precision CovRescale(NCorrMax)

      namelist/StatCorr/Name1,Name2,NIdColumns1,NIdColumns2,IdColumns1,IdColumns2,NCorr,MatrixType

      integer idataset1, idataset2, FindDataSetByName, idx1, idx2
      integer ndataset1, ndataset2, GetNPointsInDataSet
      character *40 ndataTmp
      double precision Values1(NIdMax), Values2(NIdMax)
      logical is_cov_m, is_stat_m, is_syst_m, is_syst_cm, is_cov_resize
      integer tdataset1, tdataset2 
      data is_cov_m /.false./
      data is_stat_m /.false./
      data is_syst_m /.false./
      data is_syst_cm /.false./
      data is_cov_resize /.false./

C Functions:
      integer GetBinIndex, FindIdxForCorrelation
      
C     reset statistical and systematic correlation matrices
      do i=1, npoints
         do j=1, npoints
            corr_stat(i,j) = 0.d0
            corr_syst(i,j) = 0.d0
            cov(i,j) = 0.d0
         enddo
         corr_stat(i,i) = 1.d0
         corr_syst(i,i) = 1.d0
         is_covariance(i) = .false.
      enddo

      do k=1,NCorrFiles
         print*, 'parsing correlation file ', CorrFileNames(k)

         open(51,file=CorrFileNames(k),status='old',err=11)
         read(51,NML=StatCorr,END=12,ERR=13)
         
         print*, 'matrix type: ', MatrixType

c         if (MatrixType.eq.'Statistical correlations') then
c            CovMatrixStyle = 'S'
c         elseif (MatrixType.eq.'Systematic covariance') then
c            
c         endif

         if((NIdColumns1>NIdMax).or.(NIdColumns2>NIdMax)) then
            Call HF_ERRLOG(10040001,
     $           'S: NIdMax parameter in prep_corr.f too small')
         endif

         idataset1 = FindDataSetByName(Name1)
         idataset2 = FindDataSetByName(Name2)


         if((idataset1.lt.1) .or. (idataset2.lt.1)) then
            Call HF_ERRLOG(10040003,
     $           'W: Correlation for data set requested 
     $ but not taken into account (dataset not identified)!')
            goto 1111
         endif

c RP check if new data set and if yet then and reset matrices flags         
         if(k.eq.1) then
           tdataset1 = idataset1
           tdataset2 = idataset2
         endif

         if(tdataset1.ne.idataset1 .or. tdataset2.ne.idataset2) then
            is_cov_m   = .false.
            is_stat_m  = .false.
            is_syst_m  = .false.
            is_syst_cm = .false.
            tdataset1 = idataset1
            tdataset2 = idataset2
            is_cov_resize = .false.
         endif


         do i=1,NIdColumns1
            IdIdx1(i) = GetBinIndex(idataset1,IdColumns1(i))
            if(IdIdx1(i).le.0) then
               Call HF_ERRLOG(10040004,
     $              'W: Correlation for data set requested 
     $ but not taken into account (bin not identified)!')
               goto 1111
            endif
         enddo
         do i=1,NIdColumns2
            IdIdx2(i) = GetBinIndex(idataset2,IdColumns2(i))
            if(IdIdx1(i).le.0) then
               Call HF_ERRLOG(10040004,
     $              'W: Correlation for data set requested 
     $ but not taken into account (bin not identified)!')
               goto 1111
            endif
         enddo

c get number of data points in dataset with cov matrix        
         ndataset1 = GetNPointsInDataSet(idataset1)
         ndataset2 = GetNPointsInDataSet(idataset2)

         if(ndataset1.ne.ndataset2) then
            Call HF_ERRLOG(13112618,
     $           'W: Correlation files for data sets requested 
     $ are not the same!')
            goto 1111
         endif

c check if cov matrix size correspond to number of data points
         write(ndataTmp,'(i3)') idataset1
         if(NCorr.lt.ndataset1**2) then 
            Call HF_ERRLOG(13112510,
     $          'W: Cov matrix is too small for data set '//ndataTmp)
         elseif(NCorr.gt.ndataset1**2) then
            Call HF_ERRLOG(13112511,
     $ 'W: Cov mat too big, will be resized to actual #points in set'//ndataTmp)
            is_cov_resize = .true.
         endif

c     Matrix format for covariance/correlation files
         if (MatrixType.eq.'CovMatrix'.or.
     +        MatrixType.eq.'CorrMatrix'.or.
     +        MatrixType.eq.'CovMatrix'.or.
     +        MatrixType.eq.'CorrStatMatrix'.or.
     +        MatrixType.eq.'CorrSystMatrix') then
            if(NCorr.gt.NCorrMax) then
               Call HF_ERRLOG(13012501,
     $              'S: NCorrMax parameter in prep_corr.f too small')
            endif

            if (MatrixType.eq.'CovMatrix'.or.
     +           MatrixType.eq.'CorrStatMatrix'.or.
     +           MatrixType.eq.'CorrSystMatrix') then
               do i=1,NCorr
                  read (51,*,err=1019) (matrixbuffer(i,j),j=1,NCorr+NIdColumns1)
               enddo
            elseif (MatrixType.eq.'CorrMatrix') then ! convert full correlation matrix to full covariance matrix
               do i=1,NCorr
                  read (51,*,err=1019) (matrixbuffer(i,j),j=1,NCorr+NIdColumns1+1)
               enddo
            endif

c     store delta_i for conversion of correlation matrix to covariance matrix
            if (MatrixType.eq.'CorrMatrix') then
               do i=1,NCorr
                  CovRescale(i)=matrixbuffer(i,NIdColumns1+1)
               enddo
            endif


            do i=1,NCorr
c     bin_i boundaries
               do b=1, NIdColumns1
                  Values1(b) = matrixbuffer(i,b)
               enddo
               do j=1,NCorr
c     bin_j boundaries
                  do b=1, NIdColumns1
                     Values2(b) = matrixbuffer(j,b)
                  enddo

c     get bin_i and bin_j ids                  
                  Idx1 = FindIdxForCorrelation(idataset1, NIdColumns1, IdIdx1, NIdMax, Values1)
                  Idx2 = FindIdxForCorrelation(idataset2, NIdColumns2, IdIdx2, NIdMax, Values2)
                  
                  if (MatrixType.eq.'CorrStatMatrix') then
                     corr_stat(idx1,idx2) = matrixbuffer(i,j+NIdColumns1+1)
                  elseif (MatrixType.eq.'CorrSystMatrix') then
                     corr_syst(idx1,idx2) = matrixbuffer(i,j+NIdColumns1+1)
                  elseif (MatrixType.eq.'CovMatrix') then
                     Cov(idx1,idx2) = matrixbuffer(i,j+NIdColumns1)
                  elseif (MatrixType.eq.'CorrMatrix') then
                     Cov(idx1,idx2) = matrixbuffer(i,j+NIdColumns1+1)
     +                    *CovRescale(i)*CovRescale(j)
                  endif
c                  print *,'Idx1 =', Idx1, 'Idx2 =', Idx2, 'Cov(i,j) =', Cov(idx1,idx2)

C     Mark the points for covariance matrix method:
                  is_covariance(Idx1) = .true.
                  is_covariance(Idx2) = .true.
               enddo
            enddo

            close (51)
            goto 1111
         endif


C Reading correlation values
         do j=1,NCorr
C     Allow for comments:
 89         read (51,'(A)',err=1017,end=1018) ctmp
            if (ctmp(1:1).eq.'*') then
C     Comment line, read another one
               goto 89
            endif
     
C Read the colums
            read (ctmp,*,err=1019)(buffer(i),i=1,NIdColumns1+NIdColumns2+1)

            do i=1, NIdColumns1
               Values1(i) = buffer(i)
            enddo
            do i=1, NIdColumns2
               Values2(i) = buffer(NIdColumns1+i)
            enddo

            Idx1 = FindIdxForCorrelation(idataset1, NIdColumns1, IdIdx1, NIdMax, Values1)
            Idx2 = FindIdxForCorrelation(idataset2, NIdColumns2, IdIdx2, NIdMax, Values2)
c            print *, j, 'Idx1 =', Idx1, 'Idx2 =', Idx2


C SG: Mark the points for covariance matrix method:
            is_covariance(Idx1) = .true.
            is_covariance(Idx2) = .true.
            

            if((Idx1.le.0).or.(Idx2.le.0)) then
               if(is_cov_resize) then
                  cycle
               else 
                 Call HF_ERRLOG(10040005,
     $                'W: Unable to find a proper correlation point')
                 goto 1112
               endif
            endif

            if (MatrixType.eq.'Statistical correlations') then
               corr_stat(Idx1, Idx2) = buffer(NIdColumns1+NIdColumns2+1)
               corr_stat(Idx2, Idx1) = buffer(NIdColumns1+NIdColumns2+1)
               is_stat_m = .true.

            elseif (MatrixType.eq.'Systematic correlations') then
               corr_syst(Idx1, Idx2) = buffer(NIdColumns1+NIdColumns2+1)
               corr_syst(Idx2, Idx1) = buffer(NIdColumns1+NIdColumns2+1)
               is_syst_m = .true.

            elseif (MatrixType.eq.'Systematic covariance matrix') then
               cov(Idx1, Idx2) = buffer(NIdColumns1+NIdColumns2+1)
               cov(Idx2, Idx1) = buffer(NIdColumns1+NIdColumns2+1)
               is_syst_cm = .true.

            elseif (MatrixType.eq.'Full covariance matrix') then
               cov(Idx1, Idx2) = buffer(NIdColumns1+NIdColumns2+1)
               cov(Idx2, Idx1) = buffer(NIdColumns1+NIdColumns2+1)
               is_cov_m = .true.

               if(idataset1.ne.idataset2) then
                  Call HF_ERRLOG(07061113,
     $  'S: Error, cannot use Full Cov Matrix for different datasets')
               endif

            else
               Call HF_ERRLOG(26060000,
     $              'W: Matrix type not recognised! ignore!')
            endif

C RP check that there are no stat or syst corr together with full cov:            
            if(is_stat_m.and.is_cov_m .or. is_syst_m.and.is_cov_m .or.
     $         is_syst_cm.and.is_cov_m) then
               Call HF_ERRLOG(07062013,
     $ 'S: Error, cannot set Stat or Syst matrix with Full Cov Matrix')
            endif
C RP check that there are no syst corr together with syst cov:            
            if(is_syst_m.and.is_syst_cm) then
               Call HF_ERRLOG(07061013,
     $ 'S: Error, cannot set Syst corr and cov matrix for same dataset')
            endif


 1112    enddo
         close (51)

 1111 enddo                     ! loop over correlation files
      return

 11   continue
      print '(''Can not open file vbbb'')'
      print *,CorrFileNames(i)
      call HF_stop
 12   continue
      print '(''Namelist &StatCorr NOT found'')'
      goto 14
 13   continue
      print '(''Error reading namelist &StatCorr, STOP'')'
      call HF_stop

 14   continue

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

      
      integer Function FindDataSetByName(name)
C----------------------------------------------------
C  Return idataset index for a give data set name
C----------------------------------------------------

      implicit none 
      include 'ntot.inc'
      include 'datasets.inc'      

      character *80 name
      integer i
      
      FindDataSetByName = -1

      do i=1,NInputFiles
         if(DATASETLABEL(i).eq.name) then
            if(FindDataSetByName.ne.-1) then
               Call HF_ERRLOG(10040002,
     $              'W: Ambigious data set name, 
     $              needed for statistical correlations')
               FindDataSetByName = -1
            else 
               FindDataSetByName = i
            endif
         endif
      enddo
      end

      integer Function GetNPointsInDataSet(idataset)
C----------------------------------------------------
C     Returns the # of points in a given dataset
C----------------------------------------------------

      implicit none 
      include 'ntot.inc'
      include 'datasets.inc'      

      integer idataset

      GetNPointsInDataSet = -1

      if(NDATAPOINTS(idataset).eq.0) then
         Call HF_ERRLOG(13112619,
     $       'W: Given data set has no data')
      else 
         GetNPointsInDataSet = NDATAPOINTS(idataset)
      endif
      end

      integer Function FindIdxForCorrelation(idataset, NValues, IdIdx, IdIdxMax, Values)
C----------------------------------------------------
C     Returns the idx point that abstract bins given in the IdIdx array
C     assume values Values
C----------------------------------------------------

      implicit none 
      include 'ntot.inc'
      include 'datasets.inc'      
      include 'indata.inc'


      integer m, i, idataset, NValues, idx, IdIdxMax, ngoodpoints
      logical goodpoint(ntot)
      double precision var
      integer IdIdx(IdIdxMax)
      double precision Values(IdIdxMax)

      FindIdxForCorrelation = -1

      do m=1, NDATAPOINTS(idataset)
         goodpoint(m) = .true.
      enddo
      ngoodpoints = 0

      do i=1, NValues
         do m=1, NDATAPOINTS(idataset)
            idx = DATASETIDX(idataset,m)
            var = AbstractBins(IdIdx(i),idx)
            if(var.ne.Values(i)) then
               goodpoint(m) = .false.
            endif
         enddo
      enddo                     ! loop over columns ids
      
      do m=1, NDATAPOINTS(idataset)
         if(goodpoint(m)) then
            ngoodpoints = ngoodpoints+1
            FindIdxForCorrelation = DATASETIDX(idataset,m)
         endif
      enddo
      
      if(ngoodpoints.ne.1) then
         FindIdxForCorrelation = -1
      endif
      end
      
      
      subroutine PrintCorrMatrix
C----------------------------------------------------
C     Dumps correlation matrices to correlation_matrix.txt
C     for debuggins purposes
C----------------------------------------------------
      implicit none 
      include 'ntot.inc'
      include 'indata.inc'
      include 'covar.inc'

      integer i,j
      
      print *, 'Printing correlation_matrix.txt ...'

      OPEN(UNIT=12, FILE="correlation_matrix.txt", ACTION="write", STATUS="replace")
      WRITE(12, *) 'corr_stat: (will be multiplied by the data statistical error)'

      WRITE(12, 1001) '     |', (i, i=1,npoints)

      DO i=1,npoints
         WRITE(12,1000) i, ' | ', (corr_stat(i,j), j=1,npoints)
      END DO
      

      WRITE(12, *) 'corr_syst: (will be multiplied by theory)'

      WRITE(12, 1003) '    |', (i, i=1,npoints)

      DO i=1,npoints
         WRITE(12,1002) i, ' | ', (corr_syst(i,j), j=1,npoints)
      END DO
      

 1000 format (I4,A,500F5.2)
 1002 format (I4,A,500F7.4)
 1001 format (A, 500I5)
 1003 format (A, 500I7)
            
      end
