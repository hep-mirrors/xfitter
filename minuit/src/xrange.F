c ---------------------------------------------------------
      SUBROUTINE XRANGE(FCN,FUTIL)
c ---------------------------------------------------------
c Study the allowed displacements for each coordinate...
c
#include "./d506dp.inc"
#include "./d506cm.inc"
	EXTERNAL FCN,FUTIL
	CHARACTER CBF1*22

	dimension x0(MNI), xvec(MNI), yvec(MNI), x0best(MNI)

	common /umatco/ hquad(MNI,MNI), umat(MNI,MNI), vmat(MNI,MNI), 
     &	                tvec(MNI), errmat(MNI,MNI)

	write(isyswr,100) 
100	format(1x,'XRANGE: find allowed displacements ',
     &	          'according to Hessian, for DeltaChisq=1.0')

c save initial x vector, assumed to be the location of the minimum...
	do j = 1, npar
	   x0(j) = x(j)
	enddo

      CFROM = 'XRANGE '
      NFCNFR = NFCN
      CSTATU= 'OK        '
C                 make sure starting at the right place
      CALL MNINEX(X)
      NPARX = NPAR
      CALL FCN(NPARX,GIN,FS0,U,4,FUTIL)
      NFCN = NFCN + 1
      IF (FS0 .NE. AMIN) THEN
         DF = AMIN - FS0
         WRITE (CBF1(1:12),'(G12.3)') DF
         CALL MNWARN('D','ITERATE',
     +       'function value differs from AMIN by '//CBF1(1:12) )
      ENDIF
      AMIN = FS0

	write(isyswr,110) FS0, NFCN
110	format(1x,'XRANGE: Chisq =',g14.6,' NFCN=',i6)


c keep track of any new best minimum found...
	amin4 = fs0		!error corrected: had fs1 previously
	do i = 1, npar
	   X0BEST(i) = x(i)
	enddo

c loop over the coordinates...
	do ipar = 1, npar

	   sum = 0.d0
	   do j = 1, npar
	      sum = sum + umat(ipar,j)**2
	   enddo

	   do j = 1, npar
	      yvec(j) = umat(ipar,j)/sqrt(sum)
	   enddo
c yvec is now the displacement in y space, which are the approximate diagonal
c coordinates, that gives the maximum increase in search parameter ipar for an
c increase of 1.0 in FCN.

	   do j = 1, npar
	      sum = 0.d0
	      do k = 1, npar
	         sum = sum + umat(j,k)*yvec(k)
	      enddo
	      xvec(j) = sum
	   enddo
c xvec is now the displacement in search parameter space that gives
c maximum change in coordinate ipar for increase of 1.0 in FCN.

c find the actual increases in chisq for displacement that would 
c give UP according to the quadratic approximation...
	   do j = 1, npar
	      x(j) = x0(j) + xvec(j)*sqrt(UP)
	   enddo
	   CALL FCALL(FS1,amin4,x0best,FCN,FUTIL)

	   do j = 1, npar
	      x(j) = x0(j) - xvec(j)*sqrt(UP)
	   enddo
	   CALL FCALL(FS2,amin4,x0best,FCN,FUTIL)

	   write(isyswr,120) ipar, nexofi(ipar), x0(ipar), abs(xvec(ipar))
120	   format(/1x,'Par=',i2,' TruePar=',i2,' at',e13.5,' +/-',
     &	          e13.5,' for DeltaChisq=1.0')

	   write(isyswr,130) ipar, nexofi(ipar), FS1-amin, FS2-amin, UP
130	   format(1x,'Steps in +/- directions along',i3,' TruePar',i3,
     &	' give DeltaChisq =',g12.5,',',g12.5,
     &	' (expected about',g12.5,')')

	   write(isyswr,140) ipar, nexofi(ipar)
140	   format(/5x,'Displacements that give max increase in',
     &	              ' Par=',i2,' (TruePar=',i2,') for DeltaChisq=1.0')

	   do j = 1, npar
	      write(isyswr,150) j, nexofi(j), xvec(j)
150	      format(8x,'x(',i2,' (TruePar=',i2,')) =',g14.6)
	   enddo

	   write(isyswr,160) ipar, nexofi(ipar)
160	   format(/5x,'Eigenvector components that give max increase in',
     &	              ' Par=',i2,' (TruePar=',i2,') for DeltaChisq=1.0')

	   do j = 1, npar
	      write(isyswr,170) j, yvec(j)
170	      format(8x,'y(',i2,') =',g14.6)
	   enddo

	enddo

c restore x vector to initial values, which are assumed to be the location 
c of the minimum...
	do j = 1, npar
	   x(j) = x0(j)
	enddo
	call mninex(x)

	return
	end
