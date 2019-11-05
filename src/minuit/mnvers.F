*
* $Id: mnvers.F,v 1.1.1.1 2007/01/03 00:23:21 kluget Exp $
*
* $Log: mnvers.F,v $
* Revision 1.1.1.1  2007/01/03 00:23:21  kluget
* Initial import of the source tree
*
* Revision 1.1.1.1  2007/01/02 14:35:42  kluget
* A test for herafitter
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      SUBROUTINE MNVERS(CV)
#include "./d506dp.inc"
CC         Returns the Minuit version in CV, char*6
CC
#include "./d506cm.inc"
      CHARACTER*(*) CV
      CV = CVRSN
      RETURN
      END
