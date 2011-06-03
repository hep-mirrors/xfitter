

      double precision DATEN(NTOT),E_UNC(NTOT),E_TOT(NTOT)
      double precision E_STA(NTOT)  ! stat. errors
      integer JSET(NTOT)


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


