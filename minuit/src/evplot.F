c ---------------------------------------------------------
      SUBROUTINE EVPLOT(FCN,FUTIL)
c ---------------------------------------------------------
c Print a table to plot the eigenvalues
c
#include "./d506dp.inc"
#include "./d506cm.inc"
	EXTERNAL FCN,FUTIL

	common /umatco/ hquad(MNI,MNI), umat(MNI,MNI), vmat(MNI,MNI), 
     &	                tvec(MNI), errmat(MNI,MNI)

	dimension errtmp(MNI,MNI)

	write(isyswr,100) 
100	format(1x,'EVPLOT: eigenvalues of the Hessian')

c =====================================================
c write the eigenvalues and the corresponding distances
c =====================================================
	do i = 1, npar
	   sum = 0.d0
	   do j = 1, npar
	      sum = sum + umat(j,i)**2
	   enddo

	   dist = sqrt(sum)
	   write(isyswr,110) i, dist, 2.d0/dist**2
110	   format(1x,'par=',i2,' dist=',e12.5,1x,'e.v.=2/dist**2=',e12.5)
	enddo

c ==================================================
c write the eigenvectors (normalized to unit length)
c ==================================================
	write(isyswr,150)
150	format(/1x,'Eigenvectors normalized to unit length:')

	do ipar = 1, npar
	   write(isyswr,160) ipar
160	   format(/1x,'eigenvector',i3)
	   if(npar .le. 20) then
	      write(isyswr,170) ipar, (vmat(j,ipar),j=1,npar)
170	      format(1x,i2,'::',20f6.2)
	   elseif(npar .le. 40) then
	      write(isyswr,180) ipar, (vmat(j,ipar),j=1,npar)
180	      format(1x,i2,'::',20f6.2,/5x,20f6.2)
	   else
	      write(isyswr,190) ipar, (vmat(j,ipar),j=1,npar)
190	      format(1x,i2,'::',20f6.2,/5x,20f6.2,/5x,20f6.2)
	   endif

	enddo

c ======================
c write the Error Matrix
c ======================
	write(isyswr,200) 2*UP
200	format(/1x,'Error matrix multiplied by 2*UP =',g14.6,':')

	do i = 1, npar
	do j = 1, npar
	   errtmp(i,j) = 2.d0*UP*errmat(i,j)
	enddo
	enddo

	do i = 1, npar
c **	   write(isyswr,210) 
210	   format(1x)
	   write(isyswr,220) (errtmp(i,j),j=1,i)
220	   format(1x,11e10.3)
	enddo

	return
	end
