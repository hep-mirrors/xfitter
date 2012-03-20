C
C Interface for minuit paramters
C
      integer mne          !> Number of external parameters
      parameter(MNE=200)

      double precision pkeep(MNE)   !> Store PDF parameters, for iflag=3 only
      double precision parminuitsave(MNE) !> Store for each iteration
      common/endmini/pkeep,parminuitsave
