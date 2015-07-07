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
      double precision ALPHA_MOD(NTOT)         !> Total uncorrelated errors modified
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
      integer isNuisance, isMatrix, isOffset, isExternal   !> Store type
      parameter (isNuisance = 1)
      parameter (isMatrix   = 2)
      parameter (isOffset   = 3)
      parameter (isExternal = 4)        !> Use minuit to propagate

      integer SysForm(nsysmax)

      logical doOffset
      Logical Chi2OffsFinal, Chi2OffsRecalc
      common/systscal/SysScaleFactor,SysScalingType,SysForm
     & ,doOffset, Chi2OffsFinal, Chi2OffsRecalc

C !> Add also original systematic errors for asymmetric uncertainties
      double precision BetaAsym(NSYSMAX,2,NTOT)
      logical LAsymSyst(NSYSMAX) !> asymmetric uncertainty
      Common/SystAsym/BetaAsym,omega,LAsymSyst

C !> Special flag to force recalculations
      logical ResetCommonSyst
      common/SystematicsFlags/ ResetCommonSyst

C !> Covariance matrix -> nuisance
      logical LConvertCovToNui
      common/CCovNui/LConvertCovToNui

C Keep track of data and theory systematics:
      integer iDataSyst, iTheorySyst
      parameter (iDataSyst = 1)
      parameter (iTheorySyst = 2)

      integer ISystType(NSYSMAX)        ! data/theory (something else ?)
      common/CSystType/ISystType

C Penalty term for systematic sources, 1 by default.
      double precision     SysPriorScale(nsysmax)
      common/CSystPrior/SysPriorScale
