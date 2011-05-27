      REAL FUNCTION ZEROTH(XX,i)
C-----------------------------------------------------------------------
C-
C-   Purpose and Methods: 
C-
C-   Inputs  :
C-   Outputs :
C-   Controls:
C-
C-   Created  12-JUN-2008   Voica Radescu
C-
C-----------------------------------------------------------------------
      IMPLICIT NONE

      real runif,erf
      integer i, isdrn, iseed
      real mmu, ssig
      real xx, mu, sig,stdlog,amlog
      real zeroth1

      COMMON/PARAM/mu,sig,runif

cv transform the formula from mean, std of x to log(x) 
      stdlog = sqrt(log(1+(sig/mu)**2 ) )
      amlog  = log(mu) - 0.5 * log(1+(sig/mu)**2)

      zeroth1=0.5+0.5*erf((log(xx)-amlog)/stdlog/1.4142)
c      if (sig.gt.0) then
      zeroth=zeroth1-runif
c     else
c         zeroth = (1-zeroth1)-runif
c      endif
C-----------------------------------------------------------------------
 999  RETURN
      END
