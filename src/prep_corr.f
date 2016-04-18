      subroutine prep_corr
*     ------------------------------------------------
*     Prepare statistical correlation matrix corr_stat
*     and given covariance matrix cov
*     ------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "indata.inc"
#include "covar.inc"


      integer k,i,j,b,m,n,NCorr,NIdColumns1,NIdColumns2,NIdMax,NCorrMax

      parameter (NIdMax = 5)
      character *80 Name1
      character *80 Name2
      character *80 MatrixType
      character *80 IdColumns1(NIdMax)
      character *80 IdColumns2(NIdMax)
      integer NBins1, NBins2


C     Temporary buffer to read the data (allows for comments starting with *)
      character *4096 CTmp
      double precision buffer(2*NIdMax+1)
      integer IdIdx1(NIdMax)
      integer IdIdx2(NIdMax)

      integer iCov_type_file(NSET**2)

c     matrix buffer
      parameter (NCorrMax = 100*100)
      double precision matrixbuffer(NCorrMax,NCorrMax)
      double precision CovRescale(NCorrMax)
      logical MatrixFormatIsTable

      namelist/StatCorr/Name1,Name2,NIdColumns1,NIdColumns2,
     +     IdColumns1,IdColumns2,NCorr,NBins1,NBins2,
     +     MatrixFormatIsTable,MatrixType
      integer idataset1, idataset2, FindDataSetByName, idx1, idx2
      integer ndataset1, ndataset2, GetNPointsInDataSet
      character *40 ndataTmp
      double precision BinValues1(NIdMax,NCorrMax), BinValues2(NIdMax,NCorrMax)
      double precision Values1(NIdMax), Values2(NIdMax)
      logical has_cov_m(NSET,NSET)
      logical has_corr_m(NSET,NSET)
      logical has_stat_m(NSET,NSET)
      logical has_syst_m(NSET,NSET)
      logical has_syst_cm(NSET,NSET)
      logical cov_resize

C Functions:
      integer GetBinIndex, FindIdxForCorrelation
      
C     Reset statistical and systematic correlation matrices
      do i=1, npoints
         do j=1, npoints
            corr_stat(i,j) = 0.d0
            corr_syst(i,j) = 0.d0
            cov(i,j) = 0.d0
            corr(i,j) = 0.d0
         enddo
         corr_stat(i,i) = 1.d0
         corr_syst(i,i) = 1.d0
         is_covariance(i) = .false.
         iCov_type(i) = 0
      enddo


C     Reset dataset flags for statistical and systematic correlation matrices
      do i=1,NDATASETS
         do j=1,NDATASETS
            has_cov_m(i,j) = .false.
            has_corr_m(i,j) = .false.
            has_stat_m(i,j)  = .false.
            has_syst_m(i,j)  = .false.
            has_syst_cm(i,j) = .false.
         enddo
      enddo


c     First Loop on correlation files
c     to check conflicts of incompatible correlation/covariance files associated 
c     to the same pair of datasets
      do k=1,NCorrFiles
         print*, 'Parsing correlation file: ', CorrFileNames(k)

c     Open file
         open(51,file=CorrFileNames(k),status='old',err=11)

c     Reset the MatrixFormatIsTable  flag for each cov matrix        
         MatrixFormatIsTable = .false.

c     Read namelist StatCorr
         read(51,NML=StatCorr,END=12,ERR=13)

         
         print*, 'Matrix type: ', trim(MatrixType),
     +        '; Is Table format?: ', MatrixFormatIsTable

c     Identify datasets
         idataset1 = FindDataSetByName(Name1)
         idataset2 = FindDataSetByName(Name2)
         write(ndataTmp,'(i3)') idataset1

         if (MatrixType.eq.'Statistical correlations') then
            has_stat_m(idataset1,idataset2) = .true.
            iCov_type_file(k) = iCovStatCorr
         elseif (MatrixType.eq.'Systematic correlations') then
            has_syst_m(idataset1,idataset2) = .true.
            iCov_type_file(k) = iCovSystCorr
         elseif (MatrixType.eq.'Systematic covariance matrix') then
            has_syst_cm(idataset1,idataset2) = .true.
            iCov_type_file(k) = iCovSyst
         elseif (MatrixType.eq.'Full covariance matrix') then
            has_cov_m(idataset1,idataset2) = .true.
            iCov_type_file(k) = iCovTotal
         elseif (MatrixType.eq.'Full correlation matrix') then
            has_corr_m(idataset1,idataset2) = .true.
            iCov_type_file(k) = iCovTotalCorr
         else
            Call HF_ERRLOG(14012917,
     $ 'S: MatrixType not recognised for dataset: '//ndataTmp)
         endif
         
         close (51)
      enddo                     ! End of Loop on correlation files


c     Check that two incompatible correlation/covariance files do not belong to the same pair of datasets
      do i=1,NDATASETS
         do j=1,NDATASETS
C     Check that no statistical or systematic correlation files 
c     are provided together with full covariance files
            if(has_stat_m(i,j).and.has_cov_m(i,j).or.
     $           has_syst_m(i,j).and.has_cov_m(i,j).or.
     $           has_syst_cm(i,j).and.has_cov_m(i,j).or.
     $           has_corr_m(i,j).and.has_cov_m(i,j)) then
               Call HF_ERRLOG(07062013,
     $'S: Error, cannot set Stat, Syst or Full corr matrix with Full Cov Matrix')
            endif
C     check that there are no syst corr together with syst cov:            
            if(has_syst_m(i,j).and.has_syst_cm(i,j)) then
               Call HF_ERRLOG(07061013,
     $      'S: Error, cannot set Syst corr and cov matrix for same dataset')
            endif
         enddo
      enddo

c     Second Loop on correlation files to get bin-to-bin correlations
      do k=1,NCorrFiles

         open(51,file=CorrFileNames(k),status='old',err=11)

c     Reset the MatrixFormatIsTable  flag for each cov matrix        
         MatrixFormatIsTable = .false.

         NCorr = 0
         read(51,NML=StatCorr,END=12,ERR=13)

c     Backward compatibility with old correlation files, in which NCorr is the total number of correlation entries
         if (NCorr.eq.0) then
            NCorr = NBins1 * NBins2
         endif

c     Identify datasets
         idataset1 = FindDataSetByName(Name1)
         idataset2 = FindDataSetByName(Name2)

         write(ndataTmp,'(i3)') idataset1

c     Get number of points in datasets
         ndataset1 = GetNPointsInDataSet(idataset1)
         ndataset2 = GetNPointsInDataSet(idataset2)

         do i=1,NIdColumns1
            IdIdx1(i) = GetBinIndex(idataset1,IdColumns1(i))
         enddo
         do i=1,NIdColumns2
            IdIdx2(i) = GetBinIndex(idataset2,IdColumns2(i))
         enddo

c     ******************************************
c     Begin checks on correlation file format
c     ******************************************

c     Check the number of bin ids does not exceed maximum allowed
         if((NIdColumns1>NIdMax).or.(NIdColumns2>NIdMax)) then
            Call HF_ERRLOG(10040001,
     $           'S: NIdMax parameter in prep_corr.f too small')
         endif

c     Check the number of correlation entries does not exceed maximum allowed
         if(NCorr.gt.NCorrMax) then
            Call HF_ERRLOG(13012501,
     $           'S: NCorrMax parameter in prep_corr.f too small')
         endif

c     Check that the corresponding datasets exist 
         if((idataset1.lt.1).or.(idataset2.lt.1)) then
            Call HF_ERRLOG(10040003,
     $'S: Error in reading correlation file: dataset not identified')
         endif

c     Check bin identifiers
         do i=1,NIdColumns1
            if(IdIdx1(i).le.0) then
               Call HF_ERRLOG(10040004,
     $'S: Error in reading correlation file: bin not identified')
            endif
         enddo
         do i=1,NIdColumns2
            if(IdIdx2(i).le.0) then
               Call HF_ERRLOG(10040004,
     $'S: Error in reading correlation file: bin not identified')
            endif
         enddo

c     Check if cov matrix size correspond to number of data points
         if(NCorr.lt.ndataset1*ndataset2) then 
            Call HF_ERRLOG(13112510,
     $'W: Cov matrix is smaller than bins in dataset: '//ndataTmp)
         endif

         cov_resize = .false.
         if(NCorr.gt.ndataset1*ndataset2) then
            Call HF_ERRLOG(13112511,
     $'W: Cov matrix larger than dataset bins, resize to actual points:'//ndataTmp)
            cov_resize = .true.
         endif

c     ******************************************
c     End of checks on correlation file format
c     ******************************************

c     Read correlation entries

c     ******************************************
c     Table format for covariance/correlation matrix
c     ******************************************
         if (MatrixFormatIsTable) then

            do i=1,NBins1+NIdColumns2
               read (51,*,err=1019) (matrixbuffer(i,j),j=1,NBins2+NIdColumns1)
            enddo


c     Build matrix of identification info for bins of dataset1 (bin boundaries)
            do b=1,NIdColumns1
               do i=NIdColumns2+1, NIdColumns2+NBins1
                  BinValues1(b,i-NIdColumns2) = matrixbuffer(i,b)
               enddo
            enddo
c     Build matrix of identification info for bins of dataset2 (bin boundaries)
            do b=1, NIdColumns2
               do i=NIdColumns1+1,NIdColumns1+NBins2
                  BinValues2(b,i-NIdColumns1) = matrixbuffer(b,i)
               enddo
            enddo

            do i=1,NBins1
               do j=1,NBins2

c     Make arrays of bin identification info for bin_i and bin_j
                  do b=1, NIdColumns1
                     Values1(b) = BinValues1(b,i)
                  enddo
                  do b=1, NIdColumns2
                     Values2(b) = BinValues2(b,j)
                  enddo

c     Get bin_i and bin_j ids
                  Idx1 = FindIdxForCorrelation(idataset1, NIdColumns1, IdIdx1, NIdMax, Values1)
                  Idx2 = FindIdxForCorrelation(idataset2, NIdColumns2, IdIdx2, NIdMax, Values2)

c     Check that bin_i and bin_j are identified                  
                  if ((Idx1.le.0).or.(Idx2.le.0)) then
                     if (cov_resize) then
                        cycle
                     else 
                        Call HF_ERRLOG(10040005,
     $ 'S: Unable to identify points for correlation entry')
                     endif
                  endif

                  if (MatrixType.eq.'Statistical correlations') then
c    Additional check that uncertainty for corresponding correlation matrix is given
                     if(E_STA(idx1).eq.0.d0.and.E_STA_CONST(idx1).eq.0.d0.or.
     $                  E_STA(idx2).eq.0.d0.and.E_STA_CONST(idx2).eq.0.d0) then
                        Call HF_ERRLOG(17030214,
     $'S: Stat correlations without stat error given in dataset: '//ndataTmp)
                     endif
                     corr_stat(idx1,idx2) = matrixbuffer(i+NIdColumns2,j+NIdColumns1)
                  elseif (MatrixType.eq.'Systematic correlations') then
c    Additional check that uncertainty for corresponding correlation matrix is given
                     if(E_UNC(idx1).eq.0.d0.and.E_UNC_CONST(idx1).eq.0.d0.or.
     $                  E_UNC(idx2).eq.0.d0.and.E_UNC_CONST(idx2).eq.0.d0) then
                        Call HF_ERRLOG(18030214,
     $'S: Syst correlations without syst error given in dataset: '//ndataTmp)
                     endif
                     corr_syst(idx1,idx2) = matrixbuffer(i+NIdColumns2,j+NIdColumns1)
                  elseif (MatrixType.eq.'Systematic covariance matrix') then
c    Additional check that no uncertainty for covariance matrix is given
                     if(E_UNC(idx1).ne.0.d0.and.E_UNC_CONST(idx1).ne.0.d0.or.
     $                  E_UNC(idx2).ne.0.d0.and.E_UNC_CONST(idx2).ne.0.d0) then
                        Call HF_ERRLOG(19030214,
     $'S: For sys cov matrix no sys error should be given in dataset: '//ndataTmp)
                     endif
                     cov(idx1,idx2) = matrixbuffer(i+NIdColumns2,j+NIdColumns1)
                  elseif (MatrixType.eq.'Full covariance matrix') then
c    Additional check that no uncertainty for covariance matrix is given
                     if(E_STA(idx1).ne.0.d0.and.E_STA_CONST(idx1).ne.0.d0.or.
     $                  E_STA(idx2).ne.0.d0.and.E_STA_CONST(idx2).ne.0.d0.or.
     $                  E_UNC(idx1).ne.0.d0.and.E_UNC_CONST(idx1).ne.0.d0.or.
     $                  E_UNC(idx2).ne.0.d0.and.E_UNC_CONST(idx2).ne.0.d0) then
                        Call HF_ERRLOG(20030214,
     $'S: For full cov matrix no sys, stat error should be given, set: '//ndataTmp)
                     endif
                     cov(idx1,idx2) = matrixbuffer(i+NIdColumns2,j+NIdColumns1)
                  elseif (MatrixType.eq.'Full correlation matrix') then
c    Additional check that stat and sys uncertainties with full corr matrix are given
                     if(E_STA(idx1).eq.0.d0.and.E_STA_CONST(idx1).eq.0.d0.or.
     $                  E_STA(idx2).eq.0.d0.and.E_STA_CONST(idx2).eq.0.d0.or.
     $                  E_UNC(idx1).eq.0.d0.and.E_UNC_CONST(idx1).eq.0.d0.or.
     $                  E_UNC(idx2).eq.0.d0.and.E_UNC_CONST(idx2).eq.0.d0) then
                        Call HF_ERRLOG(21030214,
     $'S: For full corr matrix sys and stat errors should be given: '//ndataTmp)
                     endif
                     corr(idx1,idx2) = matrixbuffer(i+NIdColumns2,j+NIdColumns1)
c     Reset corr_stat and corr_syst matrices to zero to avoid double counting in chi2 code            
                     corr_stat(idx1,idx2) = 0.d0
                     corr_syst(idx1,idx2) = 0.d0
                  endif
                  
c                  print *,'Idx1 =', Idx1, 'Idx2 =', Idx2, 'Cov(i,j) =', Cov(idx1,idx2)

C     Mark the points for covariance matrix method:
                  is_covariance(Idx1) = .true.
                  is_covariance(Idx2) = .true.
C     Store the type too (it is a bit mask)

                  iCov_type(Idx1) = IOR(icov_type(Idx1),iCov_type_file(k))
                  iCov_type(Idx2) = IOR(icov_type(Idx2),iCov_type_file(k))
               enddo
            enddo

         else ! end of table format for covariance/correlation matrix

c     *******************************
c     Column format for covariance/correlation matrix
c     ******************************************

C     Reading correlation values
            do j=1,NCorr
C     Allow for comments:
 89            read (51,'(A)',err=1017,end=1018) ctmp
               if (ctmp(1:1).eq.'*') then ! Comment line, skip
                  goto 89
               endif
C     Read the colums
               read (ctmp,*,err=1019)(buffer(i),i=1,NIdColumns1+NIdColumns2+1)

               do i=1, NIdColumns1
                  Values1(i) = buffer(i)
               enddo
               do i=1, NIdColumns2
                  Values2(i) = buffer(NIdColumns1+i)
               enddo
               
               Idx1 = FindIdxForCorrelation(idataset1, NIdColumns1, IdIdx1, NIdMax, Values1)
               Idx2 = FindIdxForCorrelation(idataset2, NIdColumns2, IdIdx2, NIdMax, Values2)

c     Check that bin_i and bin_j are identified                  
               if ((Idx1.le.0).or.(Idx2.le.0)) then
                  if (cov_resize) then
                     cycle
                  else 
                     Call HF_ERRLOG(10040005,
     $ 'W: Unable to find a identify points for correlation entry')
                     cycle
                  endif
               endif

               if (MatrixType.eq.'Statistical correlations') then
c    Additional check that uncertainty for corresponding correlation matrix is given
                     if(E_STA(idx1).eq.0.d0.and.E_STA_CONST(idx1).eq.0.d0.or.
     $                  E_STA(idx2).eq.0.d0.and.E_STA_CONST(idx2).eq.0.d0) then
                        Call HF_ERRLOG(17030214,
     $'S: Stat correlations without stat error given in dataset: '//ndataTmp)
                     endif
                  corr_stat(Idx1, Idx2) = buffer(NIdColumns1+NIdColumns2+1)
               elseif (MatrixType.eq.'Systematic correlations') then
c    Additional check that uncertainty for corresponding correlation matrix is given
                     if(E_UNC(idx1).eq.0.d0.and.E_UNC_CONST(idx1).eq.0.d0.or.
     $                  E_UNC(idx2).eq.0.d0.and.E_UNC_CONST(idx2).eq.0.d0) then
                        Call HF_ERRLOG(18030214,
     $'S: Syst correlations without syst error given in dataset: '//ndataTmp)
                     endif
                  corr_syst(Idx1, Idx2) = buffer(NIdColumns1+NIdColumns2+1)
               elseif (MatrixType.eq.'Systematic covariance matrix') then
c    Additional check that no uncertainty for covariance matrix is given
                     if(E_UNC(idx1).ne.0.d0.and.E_UNC_CONST(idx1).ne.0.d0.or.
     $                  E_UNC(idx2).ne.0.d0.and.E_UNC_CONST(idx2).ne.0.d0) then
                        Call HF_ERRLOG(19030214,
     $'S: For sys cov matrix no sys error should be given in dataset: '//ndataTmp)
                     endif
                  cov(Idx1, Idx2) = buffer(NIdColumns1+NIdColumns2+1)
               elseif (MatrixType.eq.'Full covariance matrix') then
c    Additional check that no uncertainty for covariance matrix is given
                     if(E_STA(idx1).ne.0.d0.and.E_STA_CONST(idx1).ne.0.d0.or.
     $                  E_STA(idx2).ne.0.d0.and.E_STA_CONST(idx2).ne.0.d0.or.
     $                  E_UNC(idx1).ne.0.d0.and.E_UNC_CONST(idx1).ne.0.d0.or.
     $                  E_UNC(idx2).ne.0.d0.and.E_UNC_CONST(idx2).ne.0.d0) then
                        Call HF_ERRLOG(20030214,
     $'S: For full cov matrix no sys, stat error should be given, set: '//ndataTmp)
                     endif
                  cov(Idx1, Idx2) = buffer(NIdColumns1+NIdColumns2+1)
               elseif (MatrixType.eq.'Full correlation matrix') then
c    Additional check that stat and sys uncertainties with full corr matrix are given
                     if(E_STA(idx1).eq.0.d0.and.E_STA_CONST(idx1).eq.0.d0.or.
     $                  E_STA(idx2).eq.0.d0.and.E_STA_CONST(idx2).eq.0.d0.or.
     $                  E_UNC(idx1).eq.0.d0.and.E_UNC_CONST(idx1).eq.0.d0.or.
     $                  E_UNC(idx2).eq.0.d0.and.E_UNC_CONST(idx2).eq.0.d0) then
                        Call HF_ERRLOG(21030214,
     $'S: For full corr matrix sys and stat errors should be given: '//ndataTmp)
                     endif
                  corr(Idx1, Idx2) = buffer(NIdColumns1+NIdColumns2+1)
               endif

C     Mark the points for covariance matrix method:
               is_covariance(Idx1) = .true.
               is_covariance(Idx2) = .true.
C Store type too

               iCov_type(Idx1) = IOR(icov_type(Idx1),iCov_type_file(k))
               iCov_type(Idx2) = IOR(icov_type(Idx2),iCov_type_file(k))
            enddo
         endif
         close (51)
      enddo                     ! End of loop over correlation files


c    Check and fill symetric off-diagonal elements
      call CheckCorrMatrix(corr_stat)
      call CheckCorrMatrix(corr_syst)
      call CheckCorrMatrix(cov)
      call CheckCorrMatrix(corr)


      goto 1200

 11   continue
      print '(''Can not open correlation file'')'
      print *,CorrFileNames(k)
      call HF_stop
 12   continue
      print '(''Namelist &StatCorr NOT foundin file'')'
      print *,CorrFileNames(k)
      call HF_stop
 13   continue
      print '(''Error reading namelist &StatCorr, STOP'')'
      print *,CorrFileNames(k)
      call HF_stop

 1017 continue
      print '(''Error reading correlation file'')'
      print *,CorrFileNames(k)
      call HF_stop
 1018 continue
      print '(''Correlation file ended before expected'')'
      print *,CorrFileNames(k)
      call HF_stop
 1019 continue
      print '(''Problem interpreting correlation file data line='',i6)',j
      print *,CorrFileNames(k)
      call HF_stop

 1200 continue
      end

      
      integer Function FindDataSetByName(name)
C----------------------------------------------------
C  Return idataset index for a give data set name
C----------------------------------------------------

      implicit none 
#include "ntot.inc"
#include "datasets.inc"

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
#include "ntot.inc"
#include "datasets.inc"

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
#include "ntot.inc"
#include "datasets.inc"
#include "indata.inc"


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
      

      subroutine CheckCorrMatrix(cmatrix)
C----------------------------------------------------
C     Check and fill symetric off-diagonal elements if only 
c     one side is provided 
C----------------------------------------------------
      implicit none 
#include "ntot.inc"
#include "indata.inc"

      double precision cmatrix(NTOT,NTOT)
      integer i,j

      do i=1,npoints
        do j=1,npoints
      
          if(cmatrix(i,j).ne.cmatrix(j,i)) then
             if(cmatrix(i,j).ne.0.d0) then 
                 cmatrix(j,i) = cmatrix(i,j)
             else if(cmatrix(j,i).ne.0.d0) then
                 cmatrix(i,j) = cmatrix(j,i)
             else 
                Call HF_ERRLOG(14070901,
     $ 'S: Error: cov matrix element for same bins is not equal')
             endif
          endif

        enddo
      enddo
            
      end
      
      subroutine PrintCorrMatrix
C----------------------------------------------------
C     Dumps correlation matrices to correlation_matrix.txt
C     for debuggins purposes
C----------------------------------------------------
      implicit none 
#include "ntot.inc"
#include "indata.inc"
#include "covar.inc"

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

      subroutine PrintCovMatrix
C----------------------------------------------------
C     Dumps covariance matrices to covariance_matrix.txt
C     for debuggins purposes
C----------------------------------------------------
      implicit none 
#include "ntot.inc"
#include "indata.inc"
#include "covar.inc"

      integer i,j
      
      print *, 'Printing covariance_matrix.txt ...'

      OPEN(UNIT=12, FILE="covariance_matrix.txt", ACTION="write", STATUS="replace")

      DO i=1,npoints
         WRITE(12,*) (cov(i,j), j=1,npoints)
      END DO
      
      end
