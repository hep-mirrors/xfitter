*
* $Id: datime.F,v 1.1.1.1 1996/02/15 17:50:37 mclareni Exp $
*
* $Log: datime.F,v $
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
* Adopted from CERNLIB by V. Kolesnikov and A. Sapronov (24.07.2014)
*
      SUBROUTINE DATIME (ID,IT)
C
C CERN PROGLIB# Z007    DATIME  DUMMY   .VERSION KERNFOR  4.22  890913
C
C-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE

      COMMON /SLATE/ ISLATE(40)

      DO J=1,6
         ISLATE(J) = 0
      END DO
      
      ID = 790929
      IT = 1200
      RETURN
      END
