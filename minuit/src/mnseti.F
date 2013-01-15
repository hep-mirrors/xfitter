*
* $Id: mnseti.F,v 1.1.1.1 2007/01/03 00:23:21 kluget Exp $
*
* $Log: mnseti.F,v $
* Revision 1.1.1.1  2007/01/03 00:23:21  kluget
* Initial import of the source tree
*
* Revision 1.1.1.1  2007/01/02 14:35:42  kluget
* A test for herafitter
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNSETI(TIT)
#include "./d506dp.inc"
CC       Called by user to set or change title of current task.
CC
#include "./d506cm.inc"
      CHARACTER*(*) TIT
      CTITL = TIT
      RETURN
      END
