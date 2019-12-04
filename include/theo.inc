      double precision
     & THEO(NTOT),THEO_MOD(NTOT),            !> Theory predictions, filled for each iteration
     & THEOA1(NTOT),THEOA2(NTOT)             !> 30/03/2017 Marina Walt, new parameters for theory predictions
                                             !> filled if experimental data is provided for a ratio sigma(A1)/sigma(A2)

      double precision THEO_FIX(NTOT)        !> Fixed theory prediction (if given by &InTheory namelist)

      double precision
     & THEO_UNC(NTOT),             !> Uncorrelated uncertainty on theory predictions
     & THEO_TOT_UP(NTOT),              !> Total up uncertainty on theory predictions
     & THEO_TOT_DOWN(NTOT),            !> Total down uncertainty on theory predictions
     $ Theo_stat(NTOT)                ! Stat. theo uncertainty 
      common/c_theo/theo,theo_mod, theoa1, theoa2, theo_fix,theo_unc,theo_tot_up,
     $     theo_tot_down, theo_stat

      double precision TheoFCN3(NTOT)
      double precision TheoModFCN3(NTOT)
      double precision Alphamodfcn3(NTOT)
      common/c_theo3/TheoFCN3, TheoModFCN3, Alphamodfcn3
