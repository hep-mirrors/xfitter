
* 
      integer NSYSMAX                  !> Maximum number of correlated systematic sources
      parameter (NSYSMAX = 300)          

      integer nSys                     !> Actual number of correlated systematic sources

      character*32 SYSTEM              !> Names of correlated systematic errors
      dimension SYSTEM(NSYSMAX) 
      double precision BETA(NSYSMAX,NTOT)  !> Influence of systematic errors on measurements
      double precision ALPHA(NTOT)         !> Total uncorrelated errors
      double precision sysa(nsysMAX,nsysMAX) !> Correlation matrix of systematics
      common/systema/ALPHA,BETA,sysa,SYSTEM,NSys

      
