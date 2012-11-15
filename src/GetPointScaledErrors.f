*     ---------------------------------------------------------  	 
*     GetPointScaledErrors 
*     calculates the rescaled statistical, uncorrelated and constant
*     errors for point ipoint
*     input: ipoint 
*            fac - factor by which theory should be rescaled (for beta shifts)
*     output: errorsta,errorunc,errorconst
*     ----------------------------------------------------

      subroutine GetPointScaledErrors(ipoint,errorsta,errorunc,errorconst)

      implicit none

      include 'ntot.inc'
      include 'steering.inc'
      include 'indata.inc'
      include 'systematics.inc'
      INCLUDE 'theo.inc'

      double precision d, t
      double precision error, errorunc, errorconst, errorsta
      integer ipoint

      d = daten(ipoint)
      t = theo(ipoint)

      error = alpha(ipoint)
      
*     DECOMPOSE ERROR
      errorunc = E_UNC(ipoint)*d/100.
      errorconst = E_STA_CONST(ipoint)*d/100.
      errorsta = error**2-errorunc**2-errorconst**2
      
      if (errorsta.gt.0) then
         errorsta = sqrt(errorsta)
      else
         errorsta = 0.
      endif
      
      
*     RESCALE ERROR
      if (ICHI2.eq.11 .or. ICHI2.eq.41) then
***   mixed scaling - decompose - scale - recombine
         if (t.gt.0) then
            errorsta = errorsta*dsqrt(abs(t/d))
            errorunc = errorunc*(abs(t/d))
         endif
         
      else if (ICHI2.eq.21) then
***   linear scaling
         errorunc = errorunc*(abs(t/d))
         errorsta = errorsta*(abs(t/d))
         errorconst = errorconst*(abs(t/d))
      else if (ICHI2.eq.31) then
***   sqrt scaling
         errorunc = errorunc*dsqrt(abs(t/d))
         errorsta = errorsta*dsqrt(abs(t/d))
         errorconst = errorconst*dsqrt(abs(t/d))
      endif
      
      return
      end
