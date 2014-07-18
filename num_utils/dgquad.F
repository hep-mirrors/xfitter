*
* $Id: dgquad.F,v 1.1.1.1 1996/04/01 15:02:15 mclareni Exp $
*
* $Log: dgquad.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:15  mclareni
* Mathlib gen
*
* Adopted from CERNLIB by V. Kolesnikov and A. Sapronov (14.07.2014)
*
*         Double Prec. version only on 32 bit machines.
      FUNCTION DGQUAD(F,A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(1),W(1)
      EXTERNAL F

      CALL D107D1(1,F,A,B,N,X,W)
      DGQUAD=X(1)
      RETURN
      END
