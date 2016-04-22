      FUNCTION QZBRNT(FUNC, X1, X2, TOLIN, IRT)
C appears to be Wijngaarden-Dekker-Brent algorithm zbrent from Numerical Recipes.
c see also 1utl/zbrnt.f
C                                                   -=-=- qzbrnt
 
C                          Return code  IRT = 1 : limits do not bracket a root;
C                                             2 : function call exceeds maximum
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      PARAMETER (ITMAX = 1000, EPS = 3.E-12)

      external func
 
      IRT = 0
      TOL = ABS(TOLIN)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.)  THEN
        WRITE (NOUT, *) 'Root must be bracketed for QZBRNT. Set = 0'
        IRT = 1
        QZBRNT=0.

	WRITE (NOUT,*) 'QZBRNT fatal error'	!make it a fatal error
	stop

c **    RETURN

      ENDIF
      C=B                       !Previously missing -- would matter if FC=0 -- (jcp 11/12/01)
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          QZBRNT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      WRITE (NOUT, *) 'QZBRNT exceeding maximum iterations.'
      IRT = 2
      QZBRNT=B
      RETURN
C**************************************************
      END
C
