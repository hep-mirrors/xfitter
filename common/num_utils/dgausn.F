
*
* $Id: gausin64.F,v 1.1.1.1 1996/04/01 15:02:43 mclareni Exp $
*
* $Log: gausin64.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:43  mclareni
* Mathlib gen
*
*
      FUNCTION DGAUSN(P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Computes a "Normal Deviate"
C     Based on G.W. Hill & A.W. Davis, Algorithm 442 Normal Deviate
C     Collected Algorithms from CACM

      CHARACTER NAME*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAME = 'DGAUSN')


      PARAMETER (C = 2.50662 82746 31000 50D0)
      PARAMETER (Z1 = 1, HF = Z1/2, C1 = 3*Z1/4, C2 = 7*Z1/8, C3 = Z1/3)

      IF(P .LE. 0 .OR. P .GE. 1) THEN
       H=0
       WRITE(ERRTXT,101) P
c       CALL MTLPRT(NAME,'G105.1',ERRTXT)
      ELSEIF(P .EQ. HF) THEN
       H=0
      ELSE
       X=P
       IF(P .GT. HF) X=1-P
       X=SQRT(-2*LOG(X))
       X=X-((7.47395*X+494.877)*X+1637.720)/
     1     (((X+117.9407)*X+908.401)*X+659.935)
       IF(P .LT. HF) X=-X
       S=X**2
       Z=C*(P-DFREQ(X))*EXP(HF*S)
       H=(((((C1*S+C2)*Z+X)*X+HF)*C3*Z+HF*X)*Z+1)*Z+X
      ENDIF
      DGAUSN=H
      RETURN
  101 FORMAT('ARGUMENT P =',1P,D15.8,' NOT IN RANGE')
      END
