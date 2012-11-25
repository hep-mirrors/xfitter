
* 
      integer NSYSMAX                  !> Maximum number of correlated systematic sources
      parameter (NSYSMAX = 300)          

      integer nSys                     !> Actual number of correlated systematic sources

      character*32 SYSTEM              !> Names of correlated systematic errors
      dimension SYSTEM(NSYSMAX) 
      double precision BETA(NSYSMAX,NTOT)  !> Influence of systematic errors on measurements
      double precision ALPHA(NTOT)         !> Total uncorrelated errors
      double precision ALPHA_MOD(NTOT)         !> Total uncorrelated errors modified
      double precision sysa(nsysMAX,nsysMAX) !> Correlation matrix of systematics
      common/systema/ALPHA,ALPHA_MOD,BETA,sysa,SYSTEM,NSys

      double precision SysScaleFactor(nsysmax)  !> Input scale factors
      !> Scaling behaviour:
      integer isNoRescale, isLinear, isPoisson, isSqrt, isLogNorm
      parameter (isNoRescale = 0)
      parameter (isLinear  = 1)
      parameter (isPoisson = 2)
      parameter (isSqrt = isPoisson)  ! Same  
      parameter (isLogNorm = 3)

      integer SysScalingType(nsysmax)      !> Store scaling behaviour              
      
      !> Type, form: 
      integer isNuisance, isMatrix, isOffset   !> Store type
      parameter (isNuisance = 1)
      parameter (isMatrix   = 2)
      parameter (isOffset   = 3)

      integer SysForm(nsysmax)

      common/systscal/SysScaleFactor,SysScalingType,SysForm

C !> Add also original systematic errors for asymmetric uncertainties
      double precision BetaAsym(NSYSMAX,2,NTOT)
      logical LAsymSyst(NSYSMAX) !> asymmetric uncertainty
      Common/SystAsym/BetaAsym,LAsymSyst

