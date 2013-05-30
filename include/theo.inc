      double precision
     & THEO(NTOT),THEO_MOD(NTOT)             !> Theory predictions, filled for each iteration

      double precision THEO_FIX(NTOT)        !> Fixed theory prediction (if given by &InTheory namelist)
      
      common/c_theo/theo,theo_mod, theo_fix
