C
C Data information
C
      integer npoints       !> Actual number of data points
      common/CNDataPoints/NPoints


      double precision DATEN(NTOT)   !> Data values
     $     ,E_UNC(NTOT)              !> Uncorelated uncertainty 
     $     ,E_TOT(NTOT)              !> Total uncertainty
      double precision E_STA(NTOT)   !> Stat. errors
      integer JSET(NTOT)             !> Reference to dataset index


C SG: Abstract binning:
      double precision  AbstractBins(NBinDimensionMax,NTOT)

C SG:  k-factors:
      double precision KFactors(NKFactMax,NTOT)

C SG: index of theory bin, if needed
      integer IndexTheoryBin(NTOT)

      common/indata2/
     + DATEN,E_UNC,E_TOT, E_STA,
     + JSET,IndexTheoryBin

      common/CAbstractBins/AbstractBins,KFactors


