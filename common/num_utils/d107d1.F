*
* $Id: d107d1.F,v 1.1.1.1 1996/04/01 15:02:14 mclareni Exp $
*
* $Log: d107d1.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:14  mclareni
* Mathlib gen
*
*
      SUBROUTINE D107D1(MODE,F,A,B,N,X,W)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*6 NAME(2)
      DATA NAME /'DGQUAD','DGSET'/
#include "d107rd.inc"

