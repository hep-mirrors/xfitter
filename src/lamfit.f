      double precision function lamfit(Q2,X,Y,fitpar)
      
C
C fitpar(1):  C       for Q2=0.2
C fitpar(2):  lambda  for Q2=0.2
C ....
C fitpar(47): C for Q2=150
C fitpar(48): lambda Q2=150
C
      implicit none

      include 'steering.inc'
      include 'fitlam.inc'
C----------------------------------------
      
      double precision  Q2,X,Y,fitpar(100)

      integer iq,iqbin
      
      double precision c,lam, rfd, lam_constd
C---------------------------------------

C From Q2 find bin number:

      do iq=1,nq2h1
         if (abs(q2-q2vh1(iq)).lt.0.01) then
            iqbin = iq
            goto 17
         endif
      enddo
 17   continue

      c   = fitpar( (iqbin-1)*2 + 1)
      lam = fitpar( (iqbin-1)*2 + 2)
      
      rfd = 1.0d0 * rf
      lam_constd = 1.0d0 * lam_const

      if (.not.FITLAMBDA2) then
C --- Ordinary fit formula ------------------------
	lamfit = c * (X**(-1.00d0*lam)) *
     &		( 1.00d0 - Y*Y/(1+(1-Y)**2)*Rfd/(Rfd+1.00d0) )
      else
C --- Modified fit formula ------------------------
	lamfit = c * ( X**(-1.00d0*lam_constd+lam*dlog(X)) ) *
     &		( 1.00d0 - Y*Y/(1+(1-Y)**2)*Rfd/(Rfd+1.00d0) )
      endif

      end

