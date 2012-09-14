
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
      logical SysAdditive(nsysmax)              !> Additive error
      common/systscal/SysScaleFactor,SysAdditive
