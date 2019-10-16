*
* $Id: dfinv.F,v 1.1.1.1 1996/02/15 17:49:03 mclareni Exp $
*
* $Log: dfinv.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
* Adopted from CERNLIB by V. Kolesnikov and A. Sapronov (24.07.2014)
*
          SUBROUTINE          DFINV(N,A,IDIM,IR)
          INTEGER             IR(*)
          DOUBLE PRECISION    A(IDIM,*),ZERO,     X, Y, TI
          DOUBLE PRECISION    S31, S32, S33, S34, DOTF
          CHARACTER*6         HNAME
          DATA      HNAME               /  ' DFINV'  /
          DOTF(X,Y,S31)  =  X*Y + S31
          DATA      ZERO      /  0.D0  /
#include "finv.inc"
          RETURN
          END
