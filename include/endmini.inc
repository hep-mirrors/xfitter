C
C Interface for minuit paramters
C
      integer mne          !> Number of external parameters
      parameter(MNE=100)

      double precision pkeep(MNE)   !> Store PDF parameters

      common/endmini/pkeep
