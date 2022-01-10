*
* $Id: mnvert.F,v 1.1.1.1 2007/01/03 00:23:21 kluget Exp $
*
* $Log: mnvert.F,v $
* Revision 1.1.1.1  2007/01/03 00:23:21  kluget
* Initial import of the source tree
*
* Revision 1.1.1.1  2007/01/02 14:35:42  kluget
* A test for herafitter
*
* Revision 1.2  1996/03/15 18:02:54  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      SUBROUTINE MNVERT(A,L,M,N,IFAIL)
#include "./d506dp.inc"
C     inverts a symmetric matrix.   matrix is first scaled to
C     have all ones on the diagonal (equivalent to change of units)
C     but no pivoting is done since matrix is positive-definite.
C     
#include "./d506cm.inc"
      DIMENSION A(L,M) ,PP(MNI), Q(MNI),  S(MNI)
      IFAIL=0
      IF (N .LT. 1)  GO TO 100
      IF (N .GT. MAXINT)  GO TO 100
C     scale matrix by sqrt of diag elements
      DO I=1,N
         SI = A(I,I)
         IF (SI .lt. 0) GOTO 100
         S(I) = 1.0/SQRT(SI)
      END DO
      DO I= 1, N
         DO J= 1, N
            A(I,J) = A(I,J) *S(I)*S(J)
         END DO
      END DO
C     . . . start main loop . . . .
      DO I=1,N
         K = I
C     preparation for elimination step1
         IF (A(K,K) .EQ. ZERO)  GO TO 100
         Q(K)=1./A(K,K)
         PP(K) = 1.0
         A(K,K)=0.0
         KP1=K+1
         KM1=K-1
         IF(KM1.lt.0) GO TO 100
         IF(KM1.eq.0) GO TO 50
         DO J=1,KM1
            PP(J)=A(J,K)
            Q(J)=A(J,K)*Q(K)
            A(J,K)=0.0
         END DO
 50      IF( (K-N).eq.0  ) GO TO 60
         IF( (K-N).gt.0  ) GO TO 100
         DO J=KP1,N
            PP(J)=A(K,J)
            Q(J)=-A(K,J)*Q(K)
            A(K,J)=0.0
         END DO
C     elimination proper
 60      DO J=1,N
            DO K=J,N
               A(J,K)=A(J,K)+PP(J)*Q(K)
            END DO
         END DO
      END DO
C     elements of left diagonal and unscaling
      DO J= 1, N
         DO K= 1, J
            A(K,J) = A(K,J) *S(K)*S(J)
            A(J,K) = A(K,J)
         END DO
      END DO
      RETURN
C     failure return
 100  IFAIL=1
      RETURN
      END
