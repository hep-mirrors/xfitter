#include "dimensions.h"
* 
      integer NSYSMAX                  !> Maximum number of correlated systematic sources
      parameter (NSYSMAX = NSYSMAX_C)          

      integer nSys                     !> Actual number of correlated systematic sources
      integer nSysData             ! data systematics (not including PDFs)

      character*64 SYSTEM              !> Names of correlated systematic errors
      dimension SYSTEM(NSYSMAX) 
      double precision BETA(NSYSMAX,NTOT)  !> Influence of systematic errors on measurements
      double precision ALPHA(NTOT)         !> Total uncorrelated errors
      double precision ALPHA_MOD(NTOT)         !> Total uncorrelated errors modified !deprecated, to be removed
      double precision sysa(nsysMAX,nsysMAX) !> Correlation matrix of systematics

      double precision omega(NSYSMAX,NTOT)  !> Quadratic term for influence of syst. errors on measurements.

      common/systema/ALPHA,ALPHA_MOD,BETA,sysa,SYSTEM,NSys, Nsysdata

      integer n_syst_meas(NSYSMAX)         ! !> Number of measurements syst. source affects
      integer syst_meas_idx(NTOT,NSYSMAX)  ! !> data points syst. source affects
      Common/Sysmeas/n_syst_meas,syst_meas_idx

      double precision SysScaleFactor(nsysmax)  !> Input scale factors
      !> Scaling behaviour:
      integer isNoRescale, isLinear, isPoisson, isSqrt, isLogNorm
      parameter (isNoRescale = 0)
      parameter (isLinear  = 1)
      parameter (isPoisson = 2)
      parameter (isSqrt = isPoisson)  ! Same  
      parameter (isLogNorm = 3)
      integer n_sys_scaling_max       !> Keep track of different scaling methods
      parameter (n_sys_scaling_max = 3)

      integer SysScalingType(nsysmax)      !> Store scaling behaviour              
      
      !> Type, form: 
      integer isNuisance, isMatrix, isOffset, isExternal, isTNuisance, isTExternal   !> Store type
      parameter (isNuisance = 1)
      parameter (isMatrix   = 2)
      parameter (isOffset   = 3)
      parameter (isExternal = 4)

      integer SysForm(nsysmax)

      logical doOffset
      Logical Chi2OffsFinal, Chi2OffsRecalc
      common/systscal/SysScaleFactor,SysScalingType,SysForm
     & ,doOffset, Chi2OffsFinal, Chi2OffsRecalc

C !> Add also original systematic errors for asymmetric uncertainties
      double precision BetaAsym(NSYSMAX,2,NTOT)
      logical LAsymSyst(NSYSMAX) !> asymmetric uncertainty
      Common/SystAsym/BetaAsym,omega,LAsymSyst

C !> gamma omega, shifts and uncorrelated uncertainties for external usage
      double precision scgamma(NSYSMAX,NTOT)
      double precision scomega(NSYSMAX,NTOT)
      double precision sysshift(NSYSMAX)
      double precision scerrors(NTOT)
      common/systexport/scgamma,scomega,sysshift,scerrors


C !> Special flag to force recalculations
      logical ResetCommonSyst
      common/SystematicsFlags/ ResetCommonSyst

C !> Covariance matrix -> nuisance
      logical LConvertCovToNui
      logical do_reduce
      logical useBlas
      common/CCovNui/LConvertCovToNui, do_reduce, useBlas

C Keep track of data and theory systematics:
      integer iDataSyst, iTheorySyst
      parameter (iDataSyst = 1)
      parameter (iTheorySyst = 2)

      integer ISystType(NSYSMAX)        ! data/theory (something else ?)
      common/CSystType/ISystType

C Penalty term for systematic sources, 1 by default.
      double precision     SysPriorScale(nsysmax)
      common/CSystPrior/SysPriorScale

C  Errors-on-Errors settings
      logical          EoEEnabled
      integer          EoE_n_iterations   !> fixed iteration count (only if EoE_fixed_iter)
      logical          EoE_fixed_iter     !> .true. iff n_iterations was set explicitly
      integer          EoE_max_iterations !> cap for the tolerance-driven Newton loop
      double precision EoE_tolerance      !> convergence tol on max|delta theta|
      logical          EoEActive(NSYSMAX)
      double precision EoEEpsilon(NSYSMAX)

      common/CEonE/ EoEEnabled, EoE_n_iterations, EoE_fixed_iter,
     $     EoE_max_iterations, EoE_tolerance, EoEActive, EoEEpsilon

C --- Bartlett corrections for EoE
      logical BartlettEnabled
      double precision BartlettLRFactor
      double precision BartlettGoFFactor
      double precision BartlettSysFactor(NSysMax)
      double precision BartlettBTilde(NSysMax)   !> b~_theta,s (Eq. 35 of arXiv:2407.05322) per EoE source
      double precision BartlettSigmaU2(NSysMax)  !> sigma_u_hat^2 (Eq. 14) at the MLE
      double precision BartlettRatio(NSysMax)    !> r_s = j^ss / sigma_u_hat^2 (marginal); must be <= 1
C BartlettExtErr: reported error of an EXTERNAL NP in the quadratic convention,
C sqrt([A^-1]_ss)*sqrt(1+b~_theta,s) from the extended NP Hessian of the extPass
C -- the same statistic reported for nuisance NPs. Replaces MINUIT's err (whose
C curvature is the eps-dependent true-GVM one) in every reporting path. Set at
C iflag=3 when EoEEnabled; 0 when unset (consumers must test .gt. 0 and fall
C back to MINUIT's err). Persists across the band scans, unlike ersys.
      double precision BartlettExtErr(NSysMax)
C SysExtFixed: .true. for an external source whose MINUIT parameter is FIXED.
C A fixed parameter is a constant, not a fitted NP: it is excluded from the
C extended NP Hessian of the extPass and from the Bartlett factors (its
C Gamma*theta contribution to the residual stays, through ShiftExternal).
C Determined by Bartlett_CountFreePOI at each iflag=3.
      logical SysExtFixed(NSysMax)
C BartlettFailed: .true. when the factors could not be computed (missing
C theory derivatives or singular POI Hessian) and were set to zero; makes the
C degradation visible in the Results.txt legend instead of silent.
      logical BartlettFailed
      common/CBartlett/ BartlettEnabled, BartlettLRFactor, BartlettGoFFactor, BartlettSysFactor,
     $     BartlettBTilde, BartlettSigmaU2, BartlettRatio, BartlettExtErr,
     $     SysExtFixed, BartlettFailed