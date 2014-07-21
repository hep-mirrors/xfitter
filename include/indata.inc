C
C Data information
C
      integer npoints       !> Actual number of data points
      common/CNDataPoints/NPoints


      double precision DATEN(NTOT)   !> Data values
     $     ,E_UNC(NTOT)              !> Uncorelated uncertainty 
     $     ,E_UNC_Const(NTOT)        !> Uncorelated, unscaled error
     $     ,E_TOT(NTOT)              !> Total uncertainty
      double precision E_STA(NTOT)   !> Stat. errors
      double precision E_STA_CONST(NTOT) !> Stat. errors which are to be kept constant
      integer JSET(NTOT)             !> Reference to dataset index

C !> New global definitions to transfer plot information to output
      integer JPLOT(NTOT)            !> Reference to plot number

C !> New definition of stat. and uncor uncertainties:
      double precision e_stat_poisson(NTOT)  
      double precision e_stat_const(NTOT)  
      double precision e_uncor_poisson(NTOT)  
      double precision e_uncor_const(NTOT)  
      double precision e_uncor_mult(NTOT)  
      double precision e_uncor_logNorm(NTOT)  
      common/CUncErrors/e_stat_poisson,e_stat_const,
     $     e_uncor_poisson, e_uncor_const, e_uncor_mult, e_uncor_logNorm

C SG: Abstract binning:
      double precision  AbstractBins(NBinDimensionMax,NTOT)

C SG:  k-factors:
      double precision KFactors(NKFactMax,NTOT)

C SG: index of theory bin, if needed
      integer IndexTheoryBin(NTOT)

      common/indata2/
     + DATEN,E_UNC,E_TOT, E_STA, E_STA_CONST, E_UNC_Const,
     + JSET,IndexTheoryBin,JPLOT

      common/CAbstractBins/AbstractBins,KFactors

C Add covariance type:
      integer NCovTypeMax
      parameter (NCovTypeMax=5)
      integer iCovTotal, iCovTotalCorr, iCovSystCorr, iCovSyst, iCovStatCorr
      parameter (iCovTotal = 1)
      parameter (iCovTotalCorr = 2)
      parameter (iCovSyst = 4)
      parameter (iCovStatCorr = 8)
      parameter (iCovSystCorr = 16)

      logical is_covariance(NTOT)  !> Use covaraince matrix if true
      integer iCov_type(NTOT)      !> type of cov. matrix 
      common/CStatCovar/is_covariance, iCov_type


C !> Add overfitting test by "control" and "fit" sets.
      logical ControlFitSplit          !> Divide data into "control" and "fit"
      logical FitSample(NTOT)          !> true if belongs to fit sample
      integer NFitPoints, NControlPoints !> counts of points
      double precision chi2_fit, chi2_cont
      common/COverFit/NFitPoints,NControlPoints,FitSample
     $     ,ControlFitSplit
      common/COverFit2/ chi2_fit, chi2_cont
