C
C Theory derivatives w.r.t. the POIs, D_im = d theo_i / d mu_m, obtained by
C central finite differences of the THEORY ONLY (no chi2, no profiling) at the
C MLE. Filled by Bartlett_ComputeD (fcn.f) at iflag=3, before the chi2 call;
C consumed by chi2_calc_syst_shifts on the extPass to build
C     K = D^T W D          (POI block of the quadratic Hessian)
C     G = Gamma^T W D      (POI-NP coupling)
C which give the marginal NP covariance  j = (A - G K^-1 G^T)^-1  of Eq. (34)
C of arXiv:2407.05322, form-independently for external and nuisance sources.
C
C NB: differencing the chi2 instead of the theory would return the Schur
C complement H_mumu - B A^-1 B^T rather than K, collapsing the construction.
C
      integer BartMaxPOI                !> cap on the number of POIs
      parameter (BartMaxPOI = 100)

C BartlettNPOI is the number of FREE MINUIT parameters that are not external
C systematic NPs -- i.e. the true dimension of mu. It is NOT the same as
C nparFCN - nExtSyst, which is (freePOI + freeExt) - allExt = freePOI -
C fixedExt, and therefore under-counts whenever an external systematic is a
C fixed MINUIT parameter. Use BartlettNPOI (when BartlettHaveNPOI) for the
C c_CI denominator and for ndf_bart = npoints - nPOI.
      integer BartlettNPOI              !> number of free POIs
      logical BartlettHaveNPOI          !> .true. once the POI list is built
      logical BartlettHaveD             !> .true. if BartlettD is filled at the MLE
      integer BartlettPOIMinuit(BartMaxPOI)  !> MINUIT external index of each POI
      double precision BartlettD(NTOT, BartMaxPOI)

      common/CBartlettFD/ BartlettD, BartlettPOIMinuit,
     $     BartlettNPOI, BartlettHaveD, BartlettHaveNPOI
