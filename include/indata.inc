      double precision 
     + VAL_Q2,VAL_X,VAL_Y,
     + VAL_S,VAL_E,VAL_STA,
     + VAL_UNC, VAL_K, VAL_KM, VAL_KP

      common/indata/
     + VAL_Q2(NSET,NMAX),VAL_X(NSET,NMAX),VAL_Y(NSET,NMAX),
     + VAL_S(NSET,NMAX),VAL_E(NSET,NMAX),VAL_STA(NSET,NMAX),
     + VAL_UNC(NSET,NMAX), VAL_KP(NSET,NMAX),VAL_KM(NSET,NMAX),
     + VAL_K(NSET,NMAX)


      double precision DATEN(NTOT),E_UNC(NTOT),E_TOT(NTOT)
      double precision E_STA(NTOT)  ! stat. errors
      integer JSET(NTOT)
      double precision VX(NTOT),VQ2(NTOT),VY(NTOT)
      double precision VKP(NTOT),VKM(NTOT),VK(NTOT)

C SG: Add rapidity and mass ranges:
      double precision Y_Low(NTOT), Y_high(NTOT), M_Low(NTOT), M_high(NTOT)

C SG: Abstract binning:
      double precision  AbstractBins(NBinDimensionMax,NTOT)

C SG:  k-factors:
      double precision KFactors(NKFactMax,NTOT)

      common/indata2/
     + VX,VQ2,VY,DATEN,E_UNC,E_TOT, E_STA,
     + Y_low, Y_high, M_low, M_high, VKM, VKP, VK,
     + JSET

      common/CAbstractBins/AbstractBins,KFactors


