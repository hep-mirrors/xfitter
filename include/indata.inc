C
C Data information
C
      integer npoints       !> Actual number of data points
      common/CNDataPoints/NPoints


      double precision DATEN(NTOT)   !> Data values
     $     ,E_UNC(NTOT)              !> Uncorelated uncertainty 
     $     ,E_TOT(NTOT)              !> Total uncertainty
      double precision E_STA(NTOT)   !> Stat. errors
      double precision E_STA_CONST(NTOT) !> Stat. errors which are to be kept constant
      integer JSET(NTOT)             !> Reference to dataset index


C SG: Abstract binning:
      double precision  AbstractBins(NBinDimensionMax,NTOT)

C SG:  k-factors:
      double precision KFactors(NKFactMax,NTOT)

C SG: index of theory bin, if needed
      integer IndexTheoryBin(NTOT)

      common/indata2/
     + DATEN,E_UNC,E_TOT, E_STA, E_STA_CONST,
     + JSET,IndexTheoryBin

      common/CAbstractBins/AbstractBins,KFactors


C !> Add overfitting test by "control" and "fit" sets.
      logical ControlFitSplit          !> Divide data into "control" and "fit"
      logical FitSample(NTOT)          !> true if belongs to fit sample
      integer NFitPoints, NControlPoints !> counts of points
      double precision chi2_fit, chi2_cont
      common/COverFit/NFitPoints,NControlPoints,FitSample
     $     ,ControlFitSplit
      common/COverFit2/ chi2_fit, chi2_cont
