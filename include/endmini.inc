C
C Interface for minuit paramters
C
      integer mne          !> Number of external parameters
      parameter(MNE=200)

      double precision pkeep(MNE)   !> Store PDF parameters, for iflag=3 only
      double precision parminuitsave(MNE) !> Store for each iteration
      common/endmini/pkeep,parminuitsave

C
C !> Store parameters after each FCN 3 call, useful for NNPDF-like stopping
C
      integer MaxFCN3 !> Max. number of FCN=3 calls in a minuit job
      integer NFCN3   !> Actual number of FCN=3 calls.
      parameter (MaxFCN3=100)  

      double precision chi2cont3(MaxFCN3), chi2fit3(MaxFCN3) !> Store chi2 

      double precision pkeep3(MNE,MaxFCN3)  !> parameter store
      common/ckeep3/pkeep3,chi2cont3,chi2fit3
      common/cfcn3/nfcn3