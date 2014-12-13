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
