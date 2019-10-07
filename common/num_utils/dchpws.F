*
* $Id: dchpws.F,v 1.1.1.1 1996/04/01 15:02:29 mclareni Exp $
*
* $Log: dchpws.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:29  mclareni
* Mathlib gen
*
* Adopted from CERNLIB by V. Kolesnikov and A. Sapronov (14.07.2014)
* Fixed on 24 Feb 2017 by SG.
*
      SUBROUTINE DCHPWS(N,C,A)

      CHARACTER NAMECP*(*),NAMEPC*(*)
      CHARACTER*80 ERRTXT
      PARAMETER (NAMECP = 'DCHPWS', NAMEPC = 'DPWCHS')

      Double precision C(0:*),A(0:*),QU(0:101),QV(0:101)  ! Double precision (!)

      PARAMETER (R1 = 1, HF = R1/2)

      IF(N .LT. 0 .OR. N .GT. 100) THEN
       call hf_errlog(14071401,
     $ 'I: DCHPWS, NUMBER OF TERMS  N < 0  OR  > 100')
      ELSEIF(N .EQ. 0) THEN
       A(0)=C(0)
      ELSE
       QU(0)=C(N)
       QU(1)=0
       QV(0)=C(N-1)
       QV(1)=2*QU(0)
       DO 2 K = 2,N
       TT=C(N-K)-QU(0)
       QU(K)=0
       DO 1 J = 0,K-1
       T=2*QV(J)-QU(J+1)
       QU(J)=QV(J)
       QV(J)=TT
       TT=T
    1  CONTINUE
       QV(K)=TT
    2  CONTINUE
      ENDIF
      A(0)=QV(0)
      DO 3 J = 1,N
      A(J)=QV(J)-QU(J-1)
    3 CONTINUE
      RETURN

      ENTRY DPWCHS(N,A,C)

      IF(N .LT. 0 .OR. N .GT. 100) THEN
       call hf_errlog(14071401,
     $ 'I: DCHPWS, NUMBER OF TERMS  N < 0  OR  > 100')
      ELSEIF(N .EQ. 0) THEN
       QU(0)=A(0)
      ELSE
       QU(0)=A(N-1)
       QU(1)=A(N)
       IF(N .GE. 2) QU(2)=0
       DO 4 K = 2,N
       TT=QU(0)+HF*QU(2)
       QU(0)=A(N-K)+HF*QU(1)
       DO 5 J = 2,K
       T=HF*QU(J-1)
       IF(J .LE. K-2) T=T+HF*QU(J+1)
       QU(J-1)=TT
       TT=T
    5  CONTINUE
       QU(K)=T
    4  CONTINUE
      ENDIF
      DO 6 J = 0,N
      C(J)=QU(J)
    6 CONTINUE
      RETURN
      END
