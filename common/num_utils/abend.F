*
* $Id: abend.F,v 1.1.1.1 2008/10/12 14:18:02 jung Exp $
*
* $Log: abend.F,v $
* Revision 1.1.1.1  2008/10/12 14:18:02  jung
* importing sources
*
* Revision 1.1.1.1  2008/10/12 13:56:30  jung
* importing sources
*
* Revision 1.1  2005/07/21 17:49:15  mkirsano
* substitution of cernlib
*
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
cc#if defined(CERNLIB_QMSUN)
cc#include "sungs/abend.F"
cc#elif defined(CERNLIB_QMVAX)
cc#include "vaxgs/abend.F"
cc#else
      SUBROUTINE ABEND
C
C CERN PROGLIB# Z035    ABEND           .VERSION KERNFOR  4.31  911111
C ORIG.  8/02/88  JZ
C

      STOP  7
      END
