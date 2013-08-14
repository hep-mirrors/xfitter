***********************************************
*
*     NNPDFINTevolveLHA
*
*     Routine  for calculating
*     the value of all xpdfs at x and Q from replica KREP
*     in a (x,Q) point
* 
*     Similar to the evolvePDF routine of the LHAPDF package
*     with an extra input variable for the number of replica
*
************************************************

      subroutine NNPDFINTevolveLHA(X,Q,XPDF,KREP)
      IMPLICIT none
*
      INTEGER MXREP
      PARAMETER(MXREP=1e3)
*
      integer NPX,NPQ2,NPL,IX,IQ2
      parameter(NPX=60,NPQ2=50)
      parameter(NPL=3000)

      double precision Q2MIN,Q2MAX,XPDFMIN,XPDFMAX,XCH,Q2CH,XCH2
      parameter(Q2MAX=1d8,Q2CH=4d0)
      parameter(XPDFMIN=1d-9,XPDFMAX=1d0,XCH=1D-1,XCH2=0.75D0)
*
      double precision XG(NPX),Q2G(NPQ2),XPDFEV(NPX,NPQ2,-6:6,MXREP)
      common/CPDFGR/XPDFEV,XG,Q2G,Q2MIN,IX,IQ2
      
      INTEGER I,IPDF,KREP,JISEARCH,INTERP
      DOUBLE PRECISION X,Q,Q2,XPDF(-6:6)
      DOUBLE PRECISION AXB(NPX,NPQ2,-6:6), BXB(NPX,NPQ2,-6:6)
      DOUBLE PRECISION CXB(NPX,NPQ2,-6:6),TQ,DX
      DOUBLE PRECISION XPDF1(-6:6),XPDF2(-6:6)
*      
      integer m,n,nmax,mmax
      parameter(m=4,n=4)
      parameter(nmax=1e3,mmax=1e3)
      double precision dy,x1,x2,y,x1a(mmax),x2a(nmax),ya(mmax,nmax)
      integer ix1a(m),ix2a(n),J

*     Set correct scale
      Q2=Q**2d0
      
*     Evolved PDF interpolation

      IF ( X.LT.XPDFMIN .OR. X.GT.XPDFMAX ) THEN
         WRITE(6,2000) 
 2000    FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE -- STOP')
         write(6,*) "X= ",X," XMAX, XMIN = ",XPDFMAX,XPDFMIN
         !STOP
      ENDIF
*
      IF ( Q2.LT.Q2MIN .OR. Q2.GT.Q2MAX ) THEN
         WRITE(6,2001) 
 2001    FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE -- STOP')
         write(6,*) "Q2 ,Q2MIN, Q2MAX = ",Q2,Q2MIN,Q2MAX
         !STOP
      ENDIF
*
      INTERP = 1               !POLYN
*
*     CUBIC SPLINE INTERPOLATION, LOG(Q2): LINEAR INTERPOLATION
*
      if(INTERP.eq.0) then 
         
*     spline coefficients
         DO IPDF = -6,6,1
            DO IQ2 = 1,NPQ2
               CALL JSPLINE(AXB,BXB,CXB,IPDF,IQ2,KREP)
            ENDDO
         ENDDO

*     Binary search of points in grid
         IQ2 = JISEARCH(NPQ2,Q2G,Q2)
         IF (IQ2 .EQ. NPQ2) IQ2 = NPQ2-1
         IX = JISEARCH(NPX,XG,X)
         DX = X - XG(IX)     

*     Compute the values of xpdfs for two neighbouring values of Q2
*     using splines of order 3

         DO IPDF = -6,6,1
            XPDF1(IPDF) = XPDFEV(IX,IQ2,IPDF,KREP)
     1           + DX*(AXB(IX,IQ2,IPDF) + DX*(BXB(IX,IQ2,IPDF) 
     2           + DX*CXB(IX,IQ2,IPDF)) )
            XPDF2(IPDF) = XPDFEV(IX,IQ2+1,IPDF,KREP)
     1           + DX*(AXB(IX,IQ2+1,IPDF) + DX*(BXB(IX,IQ2+1,IPDF) 
     2        + DX*CXB(IX,IQ2+1,IPDF)) )
         ENDDO
     
*     Linear interpolation in log Q2

         TQ = (DLOG(Q2)-DLOG(Q2G(IQ2))) 
     1        / (DLOG(Q2G(IQ2+1))-DLOG(Q2G(IQ2)))
*     
         DO IPDF = -6,6,1
            XPDF(IPDF)  = (1.0D0-TQ)*XPDF1(IPDF) + TQ*XPDF2(IPDF)
         ENDDO
*
*     POLYNOMIAL INTERPOLATION
*        
      elseif(INTERP.eq.1) then
*
         IQ2 = JISEARCH(NPQ2,Q2G,Q2)
         IF (IQ2 .EQ. NPQ2) IQ2 = NPQ2-1
         IX = JISEARCH(NPX,XG,X)

*     Assign grid for interpolation. M, N -> order of polyN interpolation

         do I=1,M
            if(IX.ge.M/2.and.IX.le.(NPX-M/2)) IX1A(I) = IX - M/2 + I
            if(IX.lt.M/2) IX1A(I) = I
            if(IX.gt.(NPX-M/2)) IX1A(I) = (NPX - M) + I

*     Check grids
            if(IX1A(I).le.0.or.IX1A(I).gt.NPX) then
               write(6,*) "Error in grids! "
               write(6,*) "I, IXIA(I) = ",I, IX1A(I)
               call exit(-10)
            endif
         enddo

         do J=1,N
            if(IQ2.ge.N/2.and.IQ2.le.(NPQ2-N/2)) IX2A(J) = IQ2 - N/2 + J
            if(IQ2.lt.N/2) IX2A(J) = J
            if(IQ2.gt.(NPQ2-N/2)) IX2A(J) = (NPQ2 - N) + J
*     Check grids
            if(IX2A(J).le.0.or.IX2A(J).gt.NPQ2) then
               write(6,*) "Error in grids! "
               write(6,*) "J, IXIA(J) = ",J,IX2A(J)
               call exit(-10)
            endif
         enddo
            
*     Define points where to evaluate interpolation
*     Choose between linear or logarithmic (x,Q2) interpolation

         IF(X.LT.XCH)THEN
            X1=dlog(X)          
         ELSE
            X1=X
         ENDIF
         X2=dlog(Q2)

         DO IPDF = -6,6,1
            
*     Choose between linear or logarithmic (x,Q2) interpolation

            DO I=1,M
               IF(X.LT.XCH)THEN
                  X1A(I)= dlog(XG(IX1A(I)))
               ELSE
                  X1A(I)= XG(IX1A(I))
               ENDIF
               DO J=1,N
                  X2A(J) = dlog(Q2G(IX2A(J)))
                  YA(I,J) = XPDFEV(IX1A(I),IX2A(J),IPDF,KREP)
               enddo
            enddo
            
*     2D polynomial interpolation
            call polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
            XPDF(IPDF) = y

         enddo                 
      endif
*
      RETURN
*
      END


*****************************************************************
*
* CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.
* INTERPOLATION SUBROUTINES ARE TAKEN FROM
* G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,
* COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977).
*
* SUBROUTINE TAKEN FROM AAC GROUP (KUMANO et al.)
*
********************************************************************
*
      SUBROUTINE JSPLINE(B,C,D,I,J,KREP)
      IMPLICIT none
*
      INTEGER MXREP
      PARAMETER(MXREP=1e3)
*
      integer NPX,NPQ2,NPL,IX,IQ2
      parameter(NPX=60,NPQ2=50)
      parameter(NPL=3000)

      double precision Q2MIN,Q2MAX,XPDFMIN,XPDFMAX,XCH,Q2CH,XCH2
      parameter(Q2MAX=1d8,Q2CH=4d0)
      parameter(XPDFMIN=1d-9,XPDFMAX=1d0,XCH=1D-1,XCH2=0.75D0)
*
      double precision XG(NPX),Q2G(NPQ2),XPDFEV(NPX,NPQ2,-6:6,MXREP)
      common/CPDFGR/XPDFEV,XG,Q2G,Q2MIN,IX,IQ2
*
*
      INTEGER I,J,NM1,K,IB,KREP
      DOUBLE PRECISION B,C,D,T
      DIMENSION B(NPX,NPQ2,-6:6),C(NPX,NPQ2,-6:6),D(NPX,NPQ2,-6:6)
*
      NM1=NPX-1
      IF(NPX.LT.2) RETURN
      IF(NPX.LT.3) GO TO 250
      D(1,J,I)=XG(2)-XG(1)
      C(2,J,I)=(XPDFEV(2,J,I,KREP)-XPDFEV(1,J,I,KREP))/D(1,J,I)
      DO 210 K=2,NM1
         D(K,J,I)=XG(K+1)-XG(K)
         B(K,J,I)=2.0D0*(D(K-1,J,I)+D(K,J,I))
         C(K+1,J,I)=(XPDFEV(K+1,J,I,KREP)-XPDFEV(K,J,I,KREP))/D(K,J,I)
         C(K,J,I)=C(K+1,J,I)-C(K,J,I)
  210 CONTINUE
      B(1,J,I)=-D(1,J,I)
      B(NPX,J,I)=-D(NPX-1,J,I)
      C(1,J,I)=0.0D0
      C(NPX,J,I)=0.0D0
      IF(NPX.EQ.3) GO TO 215
      C(1,J,I)=C(3,J,I)/(XG(4)-XG(2))-C(2,J,I)/(XG(3)-XG(1))
      C(NPX,J,I)=C(NPX-1,J,I)/(XG(NPX)-XG(NPX-2))
     1     -C(NPX-2,J,I)/(XG(NPX-1)-XG(NPX-3))
      C(1,J,I)=C(1,J,I)*D(1,J,I)**2.0D0/(XG(4)-XG(1))
      C(NPX,J,I)=-C(NPX,J,I)*D(NPX-1,J,I)**2.0D0/(XG(NPX)-XG(NPX-3))
 215  CONTINUE
      DO 220 K=2,NPX
         T=D(K-1,J,I)/B(K-1,J,I)
         B(K,J,I)=B(K,J,I)-T*D(K-1,J,I)
         C(K,J,I)=C(K,J,I)-T*C(K-1,J,I)
 220  CONTINUE
      C(NPX,J,I)=C(NPX,J,I)/B(NPX,J,I)
      DO 230 IB=1,NM1
         K=NPX-IB
         C(K,J,I)=(C(K,J,I)-D(K,J,I)*C(K+1,J,I))/B(K,J,I)
 230  CONTINUE
      B(NPX,J,I)=(XPDFEV(NPX,J,I,KREP)-XPDFEV(NM1,J,I,KREP))/D(NM1,J,I)
     1     +D(NM1,J,I)*(C(NM1,J,I)+2.0D0*C(NPX,J,I))
      DO 240 K=1,NM1
         B(K,J,I)=(XPDFEV(K+1,J,I,KREP)-XPDFEV(K,J,I,KREP))/D(K,J,I)
     1        -D(K,J,I)*(C(K+1,J,I)+2.0D0*C(K,J,I))
         D(K,J,I)=(C(K+1,J,I)-C(K,J,I))/D(K,J,I)
         C(K,J,I)=3.0D0*C(K,J,I)
 240  CONTINUE
      C(NPX,J,I)=3.0D0*C(NPX,J,I)
      D(NPX,J,I)=D(NPX-1,J,I)
      RETURN
 250  CONTINUE
      B(1,J,I)=(XPDFEV(2,J,I,KREP)-XPDFEV(1,J,I,KREP))/(XG(2)-XG(1))
      C(1,J,I)=0.0D0
      D(1,J,I)=0.0D0
      B(2,J,I)=B(1,J,I)
      C(2,J,I)=0.0D0
      D(2,J,I)=0.0D0
      RETURN
      END
*
****************************************************************************
*     THIS FUNCTION SEARCHES "I" WHICH SATISFIES THE RELATION
*     X(I) <= Y < X(I+1) BY USING A BINARY SEARCH.
*
*     FUNCTION TAKEN FROM AAC GROUP (KUMANO et al.)
****************************************************************************

      INTEGER FUNCTION JISEARCH(N,X,Y)
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*     Dynamical memory allocation
      DOUBLE PRECISION X(*)
*
      MIN=1
      MAX=N+1
*
   10 CONTINUE
      MID=(MIN+MAX)/2
      IF(Y.LT.X(MID)) THEN
        MAX=MID
      ELSE
        MIN=MID
      END IF
      IF((MAX-MIN).GT.1) GO TO 10
*
      JISEARCH=MIN
*
      RETURN
      END

*****************************************************
*
*     polin2.f
*
*     2D interpolation of arbitrary polinomial order
*     Uses polint
*     Given arrays x1a(1:m) and x2a(1:n) of independent variables,
*     and an m by n array of function values ya(1:m,1:n) tabulated
*     at the grid points defined by x1a,x2a; and given values x1,x2
*     of the independent variable, this routine returns
*     an interpolated function value y with error dy
*
*     Taken from NR fortran
*
*****************************************************

      subroutine polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
      implicit none
*
      integer m,n,nmax,mmax
      integer j,k
      parameter(nmax=1e3,mmax=1e3)

      DOUBLE PRECISION dy,x1,x2,y,x1a(mmax),x2a(nmax),ya(mmax,nmax)
      DOUBLE PRECISION ymtmp(nmax),yntmp(nmax)

      do j=1,m
         do k=1,n
            yntmp(k)=ya(j,k)
         enddo
         call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
      enddo
      call polint(x1a,ymtmp,m,x1,y,dy)
*
      return
      end

***********************************************
*
*     polint.f
*
*     Order N polynomial interpolation using Lagrange's formula
*     as descrived in Numerical Recipees:
*     Given arrays xa and ya each of length n, and given a value
*     x, this routine returns a value y and an error estimate dy.
*     If P(x) is the polynomial of degree N-1 such that 
*     P(xa_i)=ya_i,i=1,...,n, then the returned value is y=P(x)
*     The algorithm used is Neville's algorithm
*
*******************************************************

      subroutine POLINT(xa,ya,n,x,y,dy)
      implicit none
*
      integer n,NMAX
*     Largest anticipated value of n
      parameter(nmax=1e3)
      DOUBLE PRECISION dy,x,y,xa(nmax),ya(nmax)     
      INTEGER i,m,ns
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
         dift=abs(x-xa(i))
         if(dift.lt.dif) then
            ns=i
            dif=dift
         endif
         c(i)=ya(i)
         d(i)=ya(i)
 11   enddo
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
         do i=1,n-m
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            if(den.eq.0) then
               write (*,*) 'failure in polint'
               read(5,*)
            endif
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
         enddo
         if(2*ns.lt.(n-m)) then
            dy=c(ns+1)
         else
            dy=d(ns)
            ns=ns-1
         endif
         y=y+dy
      enddo

      return
      end

***********************************************

