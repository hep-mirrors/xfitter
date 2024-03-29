*
* $Id: mncalf.F,v 1.1.1.1 2007/01/03 00:23:21 kluget Exp $
*
* $Log: mncalf.F,v $
* Revision 1.1.1.1  2007/01/03 00:23:21  kluget
* Initial import of the source tree
*
* Revision 1.1.1.1  2007/01/02 14:35:42  kluget
* A test for herafitter
*
* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni
* Minuit
*
*
      SUBROUTINE MNCALF(FCN,PVEC,YCALF,FUTIL)
#include "./d506dp.inc"
CC        Called only from MNIMPR.  Transforms the function FCN
CC        by dividing out the quadratic part in order to find further
CC        minima.    Calculates  ycalf = (f-fmin)/(x-xmin)*v*(x-xmin)
CC
#include "./d506cm.inc"
      EXTERNAL FCN,FUTIL
      DIMENSION PVEC(15)
      NPARX = NPAR
      CALL MNINEX(PVEC)
      CALL FCN(NPARX,GIN,F,U,4,FUTIL)
      NFCN = NFCN + 1
      DO I= 1, NPAR
         GRD(I) = 0.
         DO J= 1, NPAR
            M = MAX(I,J)
            N = MIN(I,J)
            NDEX = M*(M-1)/2 + N
            GRD(I) = GRD(I) + VTHMAT(NDEX) * (XT(J)-PVEC(J))
         END DO
      END DO
      DENOM = 0.
      DO I= 1, NPAR
         DENOM = DENOM + GRD(I) * (XT(I)-PVEC(I))
      END DO
      IF (DENOM .LE. ZERO)  THEN
         DCOVAR = 1.
         ISW(2) = 0
         DENOM = 1.0
      ENDIF
      YCALF = (F-APSI) / DENOM
      RETURN
      END
