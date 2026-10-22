C----------------------------------------------------------------------
C> @file xfitter_gn_bridge.f
C> @brief Read-only access to the chi2 state for the GaussNewton minimizer
C>
C> All reads of xFitter commons by the GaussNewton module happen here,
C> through the real include files, so that the common-block layouts
C> cannot drift. Nothing is written to the commons. See README.md.
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C> @brief Sizes and chi2 options
C> @param[out] NFit     fitted (FitSample) points with diagonal errors
C> @param[out] NNuis    nuisance-parameter sources (isNuisance)
C> @param[out] NCov     points treated with a covariance matrix
C> @param[out] NMatrix  isMatrix sources
C> @param[out] NOffset  isOffset sources
C> @param[out] NExtern  isExternal sources
C> @param[out] Flags    PoissonCorr, ExtraSystRescale,
C>                      FirstIterationRescale, AsymErrorsIterations > 0
C----------------------------------------------------------------------
      subroutine gn_bridge_sizes(NFit, NNuis, NCov, NMatrix, NOffset,
     $     NExtern, Flags)
      implicit none
#include "ntot.inc"
#include "indata.inc"
#include "systematics.inc"
#include "steering.inc"
      integer NFit, NNuis, NCov, NMatrix, NOffset, NExtern, Flags(4)
      integer i, k
      logical isCov(NTOT)
C--------------------------------------------------------
      NNuis = 0
      NMatrix = 0
      NOffset = 0
      NExtern = 0
      do i=1,NPoints
         isCov(i) = is_covariance(i)
      enddo
      do k=1,NSys
         if (SysForm(k).eq.isNuisance) NNuis = NNuis + 1
         if (SysForm(k).eq.isOffset)   NOffset = NOffset + 1
         if (SysForm(k).eq.isExternal) NExtern = NExtern + 1
         if (SysForm(k).eq.isMatrix) then
            NMatrix = NMatrix + 1
            do i=1,n_syst_meas(k)
               isCov(syst_meas_idx(i,k)) = .true.
            enddo
         endif
      enddo

      NFit = 0
      NCov = 0
      do i=1,NPoints
         if (isCov(i)) then
            NCov = NCov + 1
         else if (FitSample(i)) then
            NFit = NFit + 1
         endif
      enddo

      Flags(1) = 0
      Flags(2) = 0
      Flags(3) = 0
      Flags(4) = 0
      if (Chi2PoissonCorr)           Flags(1) = 1
      if (Chi2ExtraSystRescale)      Flags(2) = 1
      if (Chi2FirstIterationRescale) Flags(3) = 1
      if (AsymErrorsIterations.gt.0) Flags(4) = 1
      end

C----------------------------------------------------------------------
C> @brief Indices (1-based) of the fitted diagonal points and of the
C>        nuisance-parameter sources, in the order used by the module.
C>        Call only when gn_bridge_sizes reports NCov = 0.
C----------------------------------------------------------------------
      subroutine gn_bridge_indices(NFit, IPoint, NNuis, ISource)
      implicit none
#include "ntot.inc"
#include "indata.inc"
#include "systematics.inc"
      integer NFit, IPoint(*), NNuis, ISource(*)
      integer i, k, n
C--------------------------------------------------------
      n = 0
      do i=1,NPoints
         if (FitSample(i)) then
            n = n + 1
            if (n.le.NFit) IPoint(n) = i
         endif
      enddo
      n = 0
      do k=1,NSys
         if (SysForm(k).eq.isNuisance) then
            n = n + 1
            if (n.le.NNuis) ISource(n) = k
         endif
      enddo
      end

C----------------------------------------------------------------------
C> @brief Theory THEO of the selected points (no chi2 is computed)
C----------------------------------------------------------------------
      subroutine gn_bridge_theory(NFit, IPoint, T)
      implicit none
#include "ntot.inc"
#include "theo.inc"
      integer NFit, IPoint(NFit)
      double precision T(NFit)
      integer n
C--------------------------------------------------------
      do n=1,NFit
         T(n) = THEO(IPoint(n))
      enddo
      end

C----------------------------------------------------------------------
C> @brief Chi2 ingredients of the last GetNewChisquare call, README Eqs. (1)-(6)
C>
C> Per point: data, theory, sigma (scerrors), coefficients of |d t| and
C> t^2 in sigma^2, forced-additive flag, log-term flag and alpha.
C> Per source: scaled Gamma (scgamma, zero where the source does not
C> act), power of t in the scaling, shift (sysshift), prior scale.
C----------------------------------------------------------------------
      subroutine gn_bridge_snapshot(NFit, IPoint, NNuis, ISource,
     $     D, T, Sigma, S2DT, S2TT, Additive, LogTerm, Alph,
     $     Gamma, Power, Shift, Prior)
      implicit none
#include "ntot.inc"
#include "indata.inc"
#include "theo.inc"
#include "systematics.inc"
#include "steering.inc"
      integer NFit, IPoint(NFit), NNuis, ISource(NNuis)
      double precision D(NFit), T(NFit), Sigma(NFit), S2DT(NFit),
     $     S2TT(NFit), Alph(NFit), Gamma(NNuis,NFit), Power(NNuis),
     $     Shift(NNuis), Prior(NNuis)
      integer Additive(NFit), LogTerm(NFit)
      integer i, k, n, m, l
      integer column(NTOT)
C--------------------------------------------------------
      do i=1,NPoints
         column(i) = 0
      enddo
      do n=1,NFit
         i = IPoint(n)
         column(i) = n
         D(n) = DATEN(i)
         T(n) = THEO(i)
         Sigma(n) = scerrors(i)
         S2DT(n) = e_stat_poisson(i)**2 + e_uncor_poisson(i)**2
         S2TT(n) = e_uncor_mult(i)**2
         Additive(n) = 0
         if (LForceAdditiveData(i)) Additive(n) = 1
         LogTerm(n) = 0
         if (Chi2PoissonCorr .and. ALPHA(i).gt.0) LogTerm(n) = 1
         Alph(n) = ALPHA(i)
      enddo

      do m=1,NNuis
         k = ISource(m)
         do n=1,NFit
            Gamma(m,n) = 0d0
         enddo
         do l=1,n_syst_meas(k)
            n = column(syst_meas_idx(l,k))
            if (n.gt.0) Gamma(m,n) = scgamma(k,syst_meas_idx(l,k))
         enddo
         if (SysScalingType(k).eq.isLinear) then
            Power(m) = 1d0
         else if (SysScalingType(k).eq.isPoisson) then
            Power(m) = 0.5d0
         else
            Power(m) = 0d0
         endif
         Shift(m) = sysshift(k)
         Prior(m) = SysPriorScale(k)
      enddo
      end

C----------------------------------------------------------------------
C> @brief Sum of squared polarisation shifts added by fcn to the chi2
C----------------------------------------------------------------------
      subroutine gn_bridge_polarisation_chi2(Chi2)
      implicit none
#include "polarity.inc"
      double precision Chi2
C--------------------------------------------------------
      Chi2 = shift_polRHp**2+shift_polRHm**2+
     $     shift_polLHp**2+shift_polLHm**2+
     $     shift_polL**2+shift_polT**2
      end
