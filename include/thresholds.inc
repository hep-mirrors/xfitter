C
C Starting scale and heavy quark thresholds
C
      double precision 
     $      q0           !> Evolution starting scale
     $     ,qc           !> Charm threshold
     $     ,qb           !> Bottom threshold
     $     ,qt
      common/thresholds/q0,qc,qb,qt

      double precision 
     $      kmuc           !> Charm threshold rescaling factor
     $     ,kmub           !> Bottom threshold rescaling factor
     $     ,kmut           !> Top threshold rescaling factor
      common/thresholdsrescaling/kmuc,kmub,kmut
