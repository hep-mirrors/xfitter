      double precision
     & THEO(NTOT),THEO_MOD(NTOT)             !> Theory predictions, filled for each iteration

      double precision THEO_FIX(NTOT)        !> Fixed theory prediction (if given by &InTheory namelist)

      double precision
     & THEO_UNC(NTOT),             !> Uncorrelated uncertainty on theory predictions
     & THEO_TOT_UP(NTOT),              !> Total up uncertainty on theory predictions
     & THEO_TOT_DOWN(NTOT),            !> Total down uncertainty on theory predictions
     $ Theo_stat(NTOT)                ! Stat. theo uncertainty 
      common/c_theo/theo,theo_mod, theo_fix,theo_unc,theo_tot_up,
     $     theo_tot_down, theo_stat

C Theory uncertainty from the error-band scan (PDF eigenvector variations,
C Bartlett-scaled through bart_scale), filled by GetTheoErrorsAsym/Sym after
C the band scan and zero before it. Kept in its OWN common: c_theo is
C mirrored as a C struct (ftheor_eval.cc, chi2scan.cc), so it must not grow.
C WriteFittedPoints rebuilds THEO_TOT_UP/DOWN from these plus the :T-type
C sources on every call: idempotent across repeated calls (the property the
C 2f839fff reset introduced) without wiping the band errors (the bug the
C 2f839fff reset introduced -- therr+/- in fittedresults.txt were always 0
C for fits without :T sources since May 2023).
      double precision theo_band_up(NTOT), theo_band_down(NTOT)
      common/c_theoband/ theo_band_up, theo_band_down

      double precision TheoFCN3(NTOT)
      double precision TheoModFCN3(NTOT)
      double precision Alphamodfcn3(NTOT)
      common/c_theo3/TheoFCN3, TheoModFCN3, Alphamodfcn3

      double precision residuals(NTOT+1000+1) ! residuals for data vs theory
      common/c_resid/residuals 
