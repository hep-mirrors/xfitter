      SUBROUTINE MTMUL(L,M,N,A,B,C)
C                                                   -=-=- mtmul
C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: Qcdmsc
C Description: Miscellaneous functions and subroutines for the pac.
C ListOfFiles: mtmul qwarn qzbrnt sud anom evoluf
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Qcdmsc.f,v 1.1 97/12/21 20:34:49 wkt Exp $ 
C #Log:	Qcdmsc.f,v $
c Revision 1.1  97/12/21  20:34:49  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  TO DO MATRIX MULTIPLICATION
      DIMENSION A(11,11),B(11,11),C(11,11)
      DO 20 I=1,L
      DO 20 K=1,N
      C(I,K)=0.
      DO 20 J=1,M
20       C(I,K)=C(I,K)+A(I,J)*B(J,K)
      RETURN
      END
C
C***********************************************************
C
