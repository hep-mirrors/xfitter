*
* $Id: dfact.F,v 1.1.1.1 1996/02/15 17:49:03 mclareni Exp $
*
* $Log: dfact.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
* Adopted from CERNLIB by V. Kolesnikov and A. Sapronov (24.07.2014)
*
          SUBROUTINE          DFACT(N,A,IDIM,IR,IFAIL,DET,JFAIL)
          INTEGER             IR(*),    IPAIRF
          DOUBLE PRECISION    A(IDIM,*),DET,      ZERO,     ONE,X,Y,TF
          REAL*8                G1,       G2
          REAL*8                PIVOTF,   P,        Q,        SIZEF,  T
          DOUBLE PRECISION    S11, S12, DOTF
          CHARACTER*6         HNAME
          IPAIRF(J,K)  =  J*2**12 + K
          PIVOTF(X)    =  ABS((X))
          SIZEF(X)     =  ABS((X))
          DOTF(X,Y,S11)  =  X * Y + S11
          DATA      G1, G2              /  1.D-137, 1.D137  /
          DATA      HNAME               /  ' DFACT'  /
          DATA      ZERO, ONE           /  0.D0, 1.D0  /
          DATA      NORMAL, IMPOSS      /  0, -1  /
          DATA      JRANGE, JOVER, JUNDER  /  0, +1, -1  /
#include "fact.inc"
          RETURN
          END
