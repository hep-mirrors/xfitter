c ---------------------------------------------------------
      SUBROUTINE ITERATE(FCN,FUTIL)
c ---------------------------------------------------------
c Written by J. Pumplin (Aug 2000)
c
c This is an iterative procedure to calculate the hessian matrix
c near a minimum of FCN.  It makes a linear transformation to
c normal coordinates y(1),...,y(npar), which are related to the
c ordinary coordinate displacements from the minimum by
c x(i) = x0(i) + sum_j(umat(i,j)*y(j)) where x0 is the current
c estimate of the location of the minimum and FCN is approximated
c by FCN = FCN0 + (y(1)**2 + ... + y(npar)**2).  For a given
c sum_i(y(i)**2), i.e., a given increase in FCN, y(1) gives the
c smallest displacement sqrt(sum_(x(i)-x0(i))**2)), while y(npar)
c gives the largest displacement.
c
c The matrix umat(i,j) is also written as 
c      umat(i,j) = vmat(i,j)*tvec(j)   [no sums]
c where vmat(i,j) is an orthogonal matrix and
c      tvec(j) = sqrt(2/epsilon(j)) 
c with epsilon(j) the eigenvalue of the j-th eigenvector.
c
c For description of the algorithm, see hep-ph/0008191.
c
c =====================================================================
c
c To use this procedure:
c   (1) Get parameter values that are not too far from
c       the minimum of FCN.  If your initial values are
c       not good enough, you can improve them with one of
c       the standard Minuit commands such as
c               MIGRAD 2000
c       If you don't do this, ITERATE will find the minimum
c       by itself, but it is not very efficient.
c
c   (2) Set the characteristic scale of interest for changes
c       in FCN, using the standard Minuit command
c               SET ERRDEF nn
c       ITERATE will adjust the displacements that are used to
c       calculate derivatives, to make them give an increase 
c	in FCN of approximately nn. (We use nn=5)
c
c   (3) Invoke this subroutine using the new Minuit command
c               ITERATE nn,
c       where nn = desired number of iterations.
c       If nn=0 or is omitted, no iterations are done, but the
c       coordinates are rescaled to give an approximate increase
c       of ERRDEF along each direction. (We use nn=10)
c
c The results from ITERATE are placed in the common block /umatco/.
c In the notation of hep-ph/0008191,
c    hquad(MNI,MNI) is the Hessian matrix called H_{i,j} 
c    vmat(MNI,MNI) is the orthogonal matrix called u_{i,j}
c    tvec(MNI) are the scale factors called t_i = sqrt(2/epsilon_i) 
c              where epsilon_i are the eigenvalues of the Hessian.
c    umat(i,j) = vmat(i,j)*tvec(j)
c    errmat(MNI,MNI) = error matrix = inverse of Hessian.
c
c Factors of 2: 
c       The definition of the Hessian as the matrix of 
c       second derivatives, with no other factors, is the standard 
c       definition in mathematics, and is used here.  However, 
c       the notation in the paper has been changed so that the 
c       Hessian is the inverse of the usual Error Matrix.
c
c Printout from ITERATE is written to the same file as the rest of
c Minuit's output.  The print level is controlled by the standard
c Minuit command, e.g.,
c               SET PRINT 3
c
c =====================================================================
c
c To install this procedure, the original Minuit code has to be
c modified in its command parsing subroutine mnexcm.F as follows:
c   (1) Add ITErate to the list of commands CNAME(30)
c
c   (2) Include label 3000 in the jump table
c
c   (3) Add the function call
c          3000 CALL ITERATE(FCN,FUTIL)
c               GO TO 5000
c
c It is useful to modify mnexcm.F further by similarly adding
c the function call
c          3100 CALL MYSTUFF(FCN,FUTIL)
c               GO TO 5000
c to define the new Minuit command
c               MYSTUFF nn
c which invokes a user-definable subroutine that can contain the 
c code for whatever you want to do with the results from ITERATE.
c (Minuit already has a user-definable subroutine called 
c  STAndard, but it doesn't pass the external FCN and FUTIL to it.)
c
c The MYSTUFF routine distributed here contains two utilities 
c that can be invoked after iterate has been run:  
c mystuff 1000 and mystuff 2000.
c
c =====================================================================
c
c
#include "./d506dp.inc"
#include "./d506cm.inc"
      EXTERNAL FCN,FUTIL
      LOGICAL LDEBUG
      CHARACTER CBF1*22
	dimension x0(MNI), x0new(MNI), dmove(MNI)
	dimension x0best(MNI)
	dimension tvec1(MNI), svec(MNI)
	dimension vold(MNI,MNI)
	dimension g11(MNI), g22(MNI)
	dimension vhvec((MNI*(MNI+1))/2)
	dimension qmat(MNI,MNI)

	dimension yvec(MNI)

        parameter (ncyc=20)     !max. iters to choose step size for diagonal derivatives
        dimension tmvec(2)


	common /umatco/ hquad(MNI,MNI), umat(MNI,MNI), vmat(MNI,MNI), 
     &	                tvec(MNI), errmat(MNI,MNI)

        integer NCPU, GetParamI
c get the requested number of iterations...
      niters = WORD7(1) + 0.01
      niters = max(niters,0)
      WRITE (ISYSWR,20) niters
20	format(1x,'ITERATE BEGINNING:  requested iterations=',i5)

      LDEBUG = (IDBG(3) .GE. 1)

      NCPU = GetParamI("NCPUPumplin")
      
c get the starting point and starting FCN value...
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)

      CFROM = 'ITERATE '
      NFCNFR = NFCN
      CSTATU= 'OK        '
C                 make sure starting at the right place
      CALL MNINEX(X)
      NPARX = NPAR
      CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
      NFCN = NFCN + 1
      IF (FS1 .NE. AMIN) THEN
         DF = AMIN - FS1
         WRITE (CBF1(1:12),'(G12.3)') DF
         CALL MNWARN('D','ITERATE',
     +       'function value differs from AMIN by '//CBF1(1:12) )
      ENDIF
      AMIN = FS1


c keep track of best minimum found during each iteration...
	amin4 = fs1
	do i = 1, npar
	   X0BEST(i) = x(i)
	enddo


	if(ldebug) then
	   write(isyswr,30)
30	   format(1x,'Initial parameter values in ITERATE:',
     &         /3x,'Par',2x,'TruePar',7x,'x')
	   do i = 1, npar
	      write(isyswr,40) i, nexofi(i), x(i)
40	      format(3x,i2,5x,i2,3x,e15.7)
	   enddo
	endif

C AIMSAG is the desired sagitta = average rise in FCN for the two points
c used to calculate diagonal second derivatives.  We set it equal to
c UP, which is the quantity that can be set by the standard Minuit command
c SET ERRDEF nn.  (A fairly large size, e.g., 5 or 10, is used in our
c application, to give a broad picture of the region around the minimum
c that is hopefully insensitive to the fine-scale noise from numerical
c integrations in our FCN.)
      aimsag = abs(UP)

      WRITE (ISYSWR,50) amin, nfcn, aimsag
50    format(1x,'ITERATE STARTING POINT:  FCN=',g16.9,' NFCN=',i6,
     &   ' UP=AIMSAG=',g12.6)

c initialize to x0=current notion of minimum...
	do i = 1, npar
	   x0(i) = x(i)
	   if(gstep(i) .le. 0.d0) then
	      write(ISYSWR,60) gstep(I), I, NEXOFI(I)
60	      format(1x,'ITERATE quits because gstep =',e12.5,
     &  ' not positive for par=',i2,' (real par=',i2,
     &  ') -- perhaps due to constraint')
	      goto 1200
	   endif

c initialize umat=transformation to normal coordinates to identity matrix * scale factors...
	   do j = 1, npar
	      umat(i,j) = 0.d0
	   enddo
	   umat(i,i) = gstep(i)
	enddo

c ---------------------------------------------------------------
c ---------------------------------------------------------------
c read umat(i,j) in from file umat.dat, if that file exists...
c otherwise leave it as identity matrix.
c ***	call reumat
c *** commented out for distribution version ***
c ---------------------------------------------------------------
c ---------------------------------------------------------------

	dstep = sqrt(aimsag)

c ==========================
c Begin loop over iterations
c ==========================
c in the zero-th iteration, only the scales of the new coordinates
c (which are proportional to the old coordinates if no umat was read in)
c can be changed -- no linear combinations are taken.
	do 1000 iter = 0, niters

	if(ldebug) then
	   write(isyswr,70) iter
70	   format(1x,'parameter values beginning iteration',i3,
     &    /3x,'Par',2x,'TruePar',7x,'x')
	   do i = 1, npar
	      write(isyswr,80) i, nexofi(i), x(i)
80	      format(3x,i2,5x,i2,3x,g15.7)
	   enddo
	endif

c write umat as an orthogonal matrix times a diagonal matrix...
	call getvt


c on first pass, only reset step sizes...
	if(iter .eq. 0) goto 12345

	if(ldebug) write(isyswr,90) iter, amin, nfcn
90      format(1x,'Begin calculating diagonal elements of Hessian: ',1x,
     &  'iteration',i3,' FCN=',g16.7,' NFCN=',i6)

C get diagonal elements of second derivative matrix...
      IF (LDEBUG) WRITE (ISYSWR,100) iter
 100  format(1x,'ITERATE iteration',i2,/1x,
     &  ' Par               STEP     ',
     &  '   2nd Deriv.      1st Deriv.        SAGITA  ',
     &  '        F(+)            F(-)')

      amin4s = amin4

C SG: compute all in parallel
      if ( (NCPU.gt.1) .and. (NPAR.gt.1) ) then
         call hf_errlog(2023060901,'I: Compute Hessian using fork()')
         call compute_hessian(NPAR,
     $        NCPU, x0, x0best, dstep, amin,
     $        amin4, vhvec, cormax, cormin, G11, FCN,FUTIL)
      else
c Loop over variable parameters for diagonal second derivatives...
         do I = 1, NPAR
c     take steps in opposite directions...

            do j = 1, npar
               x(j) = x0(j) + umat(j,i)*dstep
            enddo
            CALL FCALL(FS1,amin4,x0best,FCN,FUTIL)

            do j = 1, npar
               x(j) = x0(j) - umat(j,i)*dstep
            enddo
            CALL FCALL(FS2,amin4,x0best,FCN,FUTIL)

            SAG = 0.5d0*(FS1+FS2-2.0*AMIN)
            G22(I) = 2.d0*SAG/dstep**2 !estimate of second derivative
            G11(I) = (FS1-FS2)/(2.d0*dstep) !estimate of first derivative

            sum = 0.d0
            do j = 1, npar
               sum = sum + umat(j,i)**2
            enddo
            step = dstep*sqrt(sum)

            IF (LDEBUG) WRITE (ISYSWR,31)
     &           I,step,G22(I),G11(I),SAG, FS1, FS2
 31         FORMAT (2x,I2,8x,4G16.7,1x,2G16.7)

            NDEX = I*(I+1)/2
            vhvec(NDEX) = G22(I)
         enddo

         IF(ldebug .and. (amin4 .lt. amin4s))
     &        WRITE (ISYSWR,130) amin4,amin4s-amin4
 130     format(1x,'new minimum found while calculating',
     &        ' diagonal elements:',e15.7,
     &        ' lower by',e12.5,' (not used now)')

         if(ldebug) write(isyswr,140) iter, amin, nfcn
 140     format(1x,'diagonal elements done; begin off-diag',1x,
     &        'iteration',i3,' FCN=',g16.9,' NFCN=',i6)

c first derivatives and diagonal second derivatives are finished;
c now compute the off-diagonal second derivatives.
         amin4s = amin4

         IF (NPAR .EQ. 1)  GOTO 230

         DO I= 2, NPAR
            DO J= 1, I-1
c ******************************************************************************
c evaluate the fcn at 4 new points: (a,b), (-a,-b), (a, -b), (-a, b).
c Minuit's Hesse uses just one new point (a,b), and combines with the already
c known values at (0,0), (a,0), (0,b) to get an approximation that is accurate
c to first order.  If take instead two new points, e.g., (a,b) and (-a,-b),
c can combine with the known diagonal points to get a formula that is accurate
c to second order.  Using 4 points is still formally only accurate to second
c order, but the error contributions from the derivatives F13 and F31 are only
c half is big, and the error contribution from the derivative F22 is completely
c cancelled; so it is worth doing this even though it requires 4 function calls
c for each of the npar*(npar-1)/2 off-diagonal elements of the matrix.
c ******************************************************************************

c Steps for off-diagonal elements are taken to be smaller than
c the steps for diagonal elements by factor 1/sqrt(2), so that when
c convergence is achieved, they will correspond to raising the function
c by aimsag, just like the diagonal computations.

               cmult = dstep/sqrt(2.d0)

               do k = 1, npar
                  x(k) = x0(k)+cmult*(umat(k,i) + umat(k,j))
               enddo
               CALL FCALL(FS1,amin4,x0best,FCN,FUTIL)

               do k = 1, npar
                  x(k) = x0(k)-cmult*(umat(k,i) + umat(k,j))
               enddo
               CALL FCALL(FS2,amin4,x0best,FCN,FUTIL)

               do k = 1, npar
                  x(k) = x0(k)+cmult*(umat(k,i) - umat(k,j))
               enddo
               CALL FCALL(FS3,amin4,x0best,FCN,FUTIL)

               do k = 1, npar
                  x(k) = x0(k)-cmult*(umat(k,i) - umat(k,j))
               enddo
               CALL FCALL(FS4,amin4,x0best,FCN,FUTIL)

               slow = (FS1 + FS2 - FS3 - FS4)/(4.d0*cmult**2)
               if(ldebug) write(isyswr,150)
     &              i,j,fs1-amin,fs2-amin,fs3-amin,fs4-amin,slow
 150           format(1x,'i,j=',2i3,' corners=',4e12.4,1x,
     &              'dF/dY(i)dY(j)=',e13.5)

               cmax = max(fs1-amin,fs2-amin,fs3-amin,fs4-amin)
               cmin = min(fs1-amin,fs2-amin,fs3-amin,fs4-amin)


               NDEX = (I*(I-1))/2 + J
               vhvec(NDEX) = slow

               if(NDEX .eq. 2) then
                  cormax = cmax
                  cormin = cmin
               else
                  cormax = max(cormax,cmax)
                  cormin = min(cormin,cmin)
               endif

            ENDDO
         ENDDO

      endif

      if(ldebug) write(isyswr,210) cormin, cormax
 210  format(1x,'off-diagonal corner points:',e13.5,' < F-FMIN <',e13.5)

      IF(ldebug .and. (amin4 .lt. amin4s))
     &     WRITE (ISYSWR,220) amin4, amin4s-amin4
 220  format(1x,'new minimum found while calculating',
     &     ' off-diagonal elements:',e15.7,
     &     ' lower by',e12.5,' (not used now)')

 230  CONTINUE

c get qmat = hessian matrix in y's:  qmat(i,j) = D^2F / Dy_i Dy_j.
	DO I= 1, NPAR
	   DO J= 1, I
	      NDEX = I*(I-1)/2 + J
	      tmp = vhvec(NDEX)
	      qmat(I,J) = tmp
	      qmat(J,I) = tmp
	   ENDDO
	ENDDO

c See how close qmat(i,j) is to the desired 2*(identity matrix)...
	diamin = qmat(1,1)
	diamax = qmat(1,1)
	do i = 1, npar
	   diamin = min(diamin,qmat(i,i))
	   diamax = max(diamax,qmat(i,i))
	enddo
	write(isyswr,240) diamin, diamax, iter
240	format(1x,'    diagonal elements of D^2F / Dy_i Dy_j:',
     &       e13.5,' to',e13.5,'( iter=',i4,')')

	if(npar .gt. 1) then
	   offmin = qmat(1,2)
	   offmax = qmat(1,2)
	   do i = 2, npar
	      do j = 1, i-1
	         offmin = min(offmin,qmat(i,j))
	         offmax = max(offmax,qmat(i,j))
	      enddo
	   enddo
	   write(isyswr,250) offmin, offmax, iter
250	   format(1x,'off-diagonal elements of D^2F / Dy_i Dy_j:',
     &          e13.5,' to',e13.5,'( iter=',i4,')')
	endif

c get P = the hessian matrix in x's:  P(i,j) = D^2F / Dx_i Dx_j...
	do i = 1, npar
	do j = 1, i
	   sum = 0.d0
	   do m = 1, npar
	   do n = 1, npar
	      sum = sum + qmat(m,n)*vmat(i,m)*vmat(j,n)/(tvec(m)*tvec(n))
	   enddo
	   enddo
	   P(i,j) = sum
	   P(j,i) = sum
	enddo
	enddo

c get the eigenvalues (pstar(1),pstar(2),...) of Hessian,
c and their corresponding eigenvectors....

	PRECIS = 1.e-12		!*** maybe can set this to the machine precision ***

	CALL MNEIG(P,MAXINT,NPAR,MAXINT,PSTAR,PRECIS,IFAULT)
c (MNEIG overwrites P with the eigenvectors:  P(i,j) is now the j'th eigenvector.)

	IF (IFAULT .NE. 0)  THEN
	   CALL MNWARN('W','ITERATE', 'error return from first MNEIG call')
	   GO TO 1200
	ENDIF

	PMIN = PSTAR(1)
	PMAX = PSTAR(1)
	DO IP= 1, NPAR
	   IF (PSTAR(IP) .LT. PMIN)  PMIN = PSTAR(IP)
	   IF (PSTAR(IP) .GT. PMAX)  PMAX = PSTAR(IP)
	ENDDO
	IF (LDEBUG) THEN
	   WRITE (ISYSWR,500)
500	   FORMAT(1x,'iterate: eigenvalues of the hessian =')
	   WRITE (ISYSWR,510) (PSTAR(IP),IP=1,NPAR)
510	   FORMAT(1X,5E14.6)
	ENDIF

	write(isyswr,520) pmin, pmax, iter
520	format(1x,'eigenvalue range:',g14.6,' to',g14.6,
     &       ' (iter=',i4,')')

	if(ldebug) then
	   if(min(pmin,pmax) .gt. 0.d0) then
	      write(isyswr,530) pmax/pmin, sqrt(pmax/pmin)
530	      format(1x,'eigenvalue ratio=',e12.5,' distance ratio=',e12.5)
	   endif
	endif

c if notion of hessian is not positive definite, add a constant to all
c of the negative ones to fix it.
c *** Might or might not be better to add the constant to only the negative eigenvalues.***
	IF (PMIN .LE. ZERO) THEN
	   nneg = 0
	   tmp = -2.d0*pmin

	   do i = 1, npar
	      if(pstar(i) .le. 0.d0) then
	         nneg = nneg + 1
	      endif
	      pstar(i) = pstar(i) + tmp	
	   enddo

	   WRITE (ISYSWR,540) nneg, tmp, iter
540	   FORMAT(1x,i2,' negative eigenvalue(s): add',e13.5,
     &      ' to all eigenvalues to force hessian to be',
     &      ' positive definite -- iter=',i4)
	ENDIF

c ------------------------------------------------------------------------------------
c ------------------------------------------------------------------------------------
c *** remove for distribution version ***
c test orthonormality of the eigenvectors...
c the subroutine mneig is supposed to have normalized them for us.
c	if (ldebug) then
c	   tmp = 0.d0
c	   do i = 1, npar
c	   do j = 1, npar
c	      sum = 0.d0
c	      do k = 1, npar
c	         sum = sum + P(i,k)*P(j,k)
c	      enddo
c	      if(i .eq. j) sum = sum - 1.d0
c
c	      tmp = max(tmp,abs(sum))
c	   enddo
c	   enddo
c	   if(abs(tmp) .gt. 1.e-6) then
c	      write(isyswr,550) tmp
c 550	      format(1x,'iterate: warning, first test zero=',e12.5)
c	   endif
c	endif
c ------------------------------------------------------------------------------------
c ------------------------------------------------------------------------------------

c compute new estimate of location of the minimum...
c (g11(j) is the estimate of dF/dy_j)

	do i = 1, npar
	   sum = 0.d0
	   do j = 1, npar
	      sum = sum + vmat(i,j)*g11(j)/tvec(j)
	   enddo
	   svec(i) = sum
	enddo

	do i = 1, npar
	   sum = 0.d0
	   do j = 1, npar
	      sum = sum + P(j,i)*svec(j)
	   enddo
	   tvec1(i) = -sum/pstar(i)
	enddo

	do i = 1, npar
	   sum = 0.d0
	   do j = 1, npar
	      sum = sum + P(i,j)*tvec1(j)
	   enddo
	   dmove(i) = sum
	enddo


c find the displacement in new coordinates...
	do i = 1, npar
	   sum = 0.d0
	   do j = 1, npar
	      sum = sum + vmat(j,i)*dmove(j)
	   enddo
	   yvec(i) = sum/tvec(i)
	enddo

c it may not be safe to move far from the region where the derivatives
c have been estimated, so see how far that is...
	rmax = 0.d0
	do i = 1, npar
	   rmax = max(rmax,abs(yvec(i))/dstep)
	   if(ldebug) write(isyswr,560) i, yvec(i)/dstep
560	   format(1x,'yvec(',i2,')/dstep =',g14.6)
	enddo

	if(ldebug) write(isyswr,570) rmax
570	format(1x,'rmax=',e12.5)


	rmaxc = 1.0d0		!*** set threshold for testing by steps ***
	if(rmax .gt. rmaxc) then
	   if(ldebug) write(isyswr,580) rmax, rmaxc
580	   format(1x,'rmax=',e12.5,' >',e12.5,': take small steps')
	   nsteps = 50
	   stmax = 5.0d0		!largest step to be tested
	else
	   nsteps = 1
	   stmax = 1.d0
	endif

	fbest = amin
	do istep = 1, nsteps
	   step = stmax*float(istep)/float(nsteps)
	   do i = 1, npar
	      x(i) = x0(i) + step*dmove(i)
	   enddo
	   CALL FCALL(FS1,amin4,x0best,FCN,FUTIL)

	   if(ldebug) write(isyswr,590) step, fs1-amin
590	   format(5x,'step fraction=',f8.4,' f-fmin=',e14.6)

c once we have stepped as far as the calculated position of minimum,
c continue stepping only so long as each new step gives a new best value...
	   if((step .gt. 1.d0) .and. (fs1 .gt. fbest)) goto 600

	   if(fs1 .lt. fbest) then
	      fbest = fs1
	      stbest = step
	      do i = 1, npar
	         x0new(i) = x(i)
	      enddo
	   endif
	enddo
600	continue

c install the new minimum point if it's an improvement...
	if(fbest .lt. amin) then
	   if (ldebug) then
	      if(nsteps .eq. 1) then
	         write(isyswr,610)
610	         format(1x,'calculated step was taken')
	      else
	         write(isyswr,620) stbest
620	         format(1x,'step taken was a fraction',f8.4,
     &                ' of the calculated step')
	      endif

	      write(isyswr,630)
630	      format(1x,'par',7x,'moved',12x,'from',12x,'to')
	      do i = 1, npar
	         diff = x0new(i) - x0(i)
	         write(isyswr,640) i, diff, x0(i), x0new(i)
640	         format(2x,i2,3x,e12.5,1x,2e16.8)
	      enddo
	   endif

	   do i = 1, npar
	      x0(i) = x0new(i)
	      x(i) = x0new(i)
	   enddo

	   CALL FCALL(FS1,amin4,x0best,FCN,FUTIL)

	   if (fs1 .lt. amin) then
	      if (ldebug) write(isyswr,650) fs1, amin-fs1
650	      format(1x,'new minimum installed:',e15.7,
     &             ' (',e12.5,' below old)')
	   else
	      write(isyswr,660) fs1, fs1-amin
660	      format(1x,'new "minimum"',e15.7,' (',e12.5,
     &             ' above old: fatal error SHOULD NOT HAPPEN)')
	      stop
	   endif

	   amin = fs1
	else
	   if (ldebug) write(isyswr,670)
670	   format(1x,'step was not taken')
	endif


c get the new version of the transform...
c Reverse the order of the eigenvectors, so the ones with the largest
c eigenvalues of the hessian (hence smallest eigenvalues of the error
c matrix), i.e., the best-determined eigenvectors, are listed first.
	do i = 1, npar
	   tvec(i) = sqrt(2.d0/pstar(npar+1-i))
	enddo

	do i = 1, npar
	do j = 1, npar
	   vold(i,j) = vmat(i,j)	!save copy of old ones to take dot-product with new.
	   vmat(i,j) = P(i,npar+1-j)
	enddo
	enddo

c To make the answer unique, multiply eigenvectors by -1 as
c necessary to make the largest component of each vector positive...
	do i = 1, npar
	   tmp = vmat(1,i)
	   do j = 1, npar
	      if(abs(vmat(j,i)) .gt. abs(tmp)) tmp = vmat(j,i)
	   enddo
	   if(tmp .lt. 0.d0) then
	      do j = 1, npar
	         vmat(j,i) = -vmat(j,i)
	      enddo
	   endif
	enddo

	do i = 1, npar
	do j = 1, npar
	   umat(i,j) = vmat(i,j)*tvec(j)
	enddo
	enddo


c --------------------------------------------------------------------
c --------------------------------------------------------------------
c *** remove for distribution version ***
c another orthogonality test... 
c	if(ldebug) then
c	   errmax = 0.d0
c	   do i = 1, npar
c	   do j = 1, npar
c	      sum1 = 0.d0
c	      sum2 = 0.d0
c	      do k = 1, npar
c	         sum1 = sum1 + vmat(k,i)*vmat(k,j)
c	         sum2 = sum2 + vmat(i,k)*vmat(j,k)
c	      enddo
c
c	      if(i .ne. j) then
c	         errmax = max(errmax,abs(sum1),abs(sum2))
c	      else
c	         errmax = max(errmax,abs(sum1-1.d0),abs(sum2-1.d0))
c	      endif
c
c	   enddo
c	   enddo
c	   if(abs(errmax) .gt. 1.e-6) then
c	      write(isyswr,680) errmax
c 680	      format(1x,'iterate warning:  second test zero=',e12.5)
c	   endif
c	endif
c --------------------------------------------------------------------
c --------------------------------------------------------------------
c
c at this point, the quadratic approximation is supposed to be
c chisqr = amin + y(1)**2 + y(2)**2 + ...
c
c hence, in the next iteration, we expect the first derivatives g11(i)
c to be approximately zero, and the second derivative matrix to be 
c approximately 2*(identity matrix).

c print the lengths...
	if(ldebug) then
	   do i = 1, npar
	      write(isyswr,690) i, aimsag, tvec(i),
     &	                        sqrt(aimsag)*tvec(i)
690	      format(1x,'distances to move along eigenvector direction',i3, 
     &             ' to get DeltaChisqr=1,',e12.5,
     &             ' are',e12.5,',',e12.5)
	   enddo
	endif

c (iteration 0 jumps directly to here)
12345	continue

c adjust step sizes for each direction to make the increases in Chisqr 
c to be approximately equal to the desired value set by errdef...
	amin4s = amin4

        if (ncpu.gt.1) then
           if  (iter .lt. 5) then
              call step_method_one(npar,ncpu,x0,x0best,
     $             AIMSAG,AMIN,AMIN4,ncyc,FCN, FUTIL)
           else
              call step_method_two(npar,ncpu,x0,x0best,
     $             AIMSAG,AMIN,AMIN4,ncyc,FCN, FUTIL)
           endif
        else
           do ipar = 1, npar


              if(ldebug) write(isyswr,710)
 710          format(1x)
           
              if(iter .ge. 5) goto 750
           
c *** BEGINNING WITH ITERATION 5, SWITCH TO THE SLOWER BUT 
c *** MORE ACCURATE "ALTERNATIVE" METHOD ***

              do isize = 1, ncyc
              
                 do j = 1, npar
                    x(j) = x0(j) + umat(j,ipar)*sqrt(aimsag)
                 enddo
                 CALL FCALL(FS1,amin4,x0best,FCN,FUTIL)
              
                 do j = 1, npar
                    x(j) = x0(j) - umat(j,ipar)*sqrt(aimsag)
                 enddo
                 CALL FCALL(FS2,amin4,x0best,FCN,FUTIL)
              
              
                 if(ldebug) write(isyswr,720)
     &                ipar, aimsag, fs1-amin, fs2-amin, isize, iter
 720             format(1x,'Direction #',i2,' expected UP=',g14.6,
     &      ' got',g14.6,' ,',g14.6,' (isize=',i2,' iter=',i4,')')

                 rat1 = (max(fs1,fs2) - amin)/aimsag
                 rat2 = (min(fs1,fs2) - amin)/aimsag

c ideally, rat1 and rat2 would both be 1.0; scale the step size to approximate that...
c (in bizarre case that rat2 is negative, instead increase step size.)
              
                 if(rat2 .le. 0.d0) then
                    tmp = 2.0d0
                    if(isize .gt. 2) then
                       if(tmpst1 .gt. 1.d0) tmp = min(tmp,tmpst1)
                       if(tmpst2 .gt. 1.d0) tmp = min(tmp,tmpst2)
                    endif
                 else
                    tmp = sqrt(2.d0/(rat1 + rat2))
                 endif

c force the scale change to lie between 0.5 and 2.0, so very large
c changes are avoided...
                 if(tmp .gt. 1.5d0) then
                    tmp = 2.d0 - 0.75d0/tmp
                 elseif(tmp .lt. 0.7d0) then
                    tmp = 0.5d0 + (2.d0/7.d0)*tmp
                 endif

c terminate iterations if last few multipliers were close to 1...
                 if(isize .gt. 2) then
                    if(max(abs(tmp-1.d0),
     &                   abs(tmpst1-1.d0),
     &                   abs(tmpst2-1.d0)) .lt. 0.1d0) then
                    
                       if(ldebug) write(isyswr,730)
     &              ipar,min(rat1,rat2),max(rat1,rat2),isize,iter
 730                   format(1x,'dir',i2,' rats',2f15.5,' after',
     &                   i3,' tries (iter=',i4,')')
                    
                       goto 900
                    
                    endif
                 endif

c store the previous multipliers...
                 if(isize .gt. 1) tmpst2 = tmpst1
                 tmpst1 = tmp

                 do j = 1, npar
                    umat(j,ipar) = tmp*umat(j,ipar)
                 enddo
                 print *,'SCALE umat',ipar,tmp
              enddo


              if(ldebug) write(isyswr,740) ipar, iter
 740          format(1x,
     $             'warning -- step size didn''t converge along dir=',
     &          i2,' (iter=',i4,') --- use alternative procedure')

c ===========================
c begin alternative procedure
c ===========================
       
 750          continue

c loop over forward/backward along this direction...
              igood = 0
              do idir = 1, 2

                 tmult = sqrt(aimsag)
                 if(idir .eq. 2) then
                    tmult = -tmult
                 endif

c loop over iterations to find the displacement that 
c gives the desired increase of aimsag in chisqr...
c (This scheme avoids using derivatives -- it is 
c intended to be robust in the face of discontinuities in 
c FCN which can be present due to numerical integration noise.)
                 ntry = 50
                 dbesth = -9.d99
                 dbestl = +9.d99
                 tbesth = 1.d0
                 tbestl = 1.d0
                 do itry = 1, ntry

                    do j = 1, npar
                       x(j) = x0(j) + umat(j,ipar)*tmult
                    enddo
                    CALL FCALL(FS1,amin4,x0best,FCN,FUTIL)

                    del = FS1 - amin

c **************************************************************
c EXTRA PRINTOUT FOR TESTING
	if (ldebug) write(isyswr,830) tmult, del, ipar, idir, itry
830	format(1x,'tmult=',e16.8,' del=',e12.5,
     &       ' ipar=',i2,' idir=',i2,' itry=',i3)
c **************************************************************

c go if convergence good enough...
c *** sets desired accuracy ***
                    accu = 1.e-3

                    if(abs(del-aimsag) .lt. accu*aimsag) then
                       igood = igood + idir
                       goto 880	
                    endif

                    if(itry .eq. 1) then
                       tbest = tmult
                       dbest = del
                    endif

c keep track of best above and best below...
                    if(del .ge. aimsag) then
                       if(abs(del-aimsag) .lt. abs(dbesth-aimsag)) then
                          tbesth = tmult
                          dbesth = del
                       endif
                    endif

                    if(del .le. aimsag) then
c *** error corrected 9/26/00 ***
                       if(abs(del-aimsag) .lt. abs(dbestl-aimsag)) then
                          tbestl = tmult
                          dbestl = del
                       endif
                    endif


c del should always be positive, because amin should be the local minimum;
c bail out if it isn't.
                    if(del .le. 0.d0) then
                       if(ldebug)
     $                      write(isyswr,840) (-del), ipar, idir, itry
 840       format(1x,'ITERATE: warning, FS1 below amin by',e12.5,
     &                      ' ipar=',i3,' idir=',i2,' itry=',i3)
                       goto 880
                    endif

c first time, calculate the multiplier by assuming chisqr is a quadratic 
c function of the displacement, as it is expected to be near a local minimum...
                    if(itry .eq. 1) then
                       xmult = sqrt(aimsag/del)

c force the multiplier to be not too far from 1; but 
c if it is already close to 1, leave it nearly as is...
                       ascale = 0.3d0 !*** parameter choice ***
                       xmult = 1.d0 + ascale*atan((xmult-1.d0)/ascale)
                       goto 850
                    endif

c keep track of best found so far...
                    if(abs(del-aimsag) .lt. abs(dbest-aimsag)) then
                       tbest = tmult
                       dbest = del
                    endif

c reduce step size if have stepped to the other side of the desired value...
                    if(((xmult .gt. 1.d0) .and. (del .ge. aimsag)) .or.
     &               ((xmult .lt. 1.d0) .and. (del .le. aimsag))) then
                       tmult = tbest
                       del = dbest
                       xmult = sqrt(xmult) !makes step factor closer to 1.0, i.e., smaller step size
                    endif

                    if(((xmult .gt. 1.d0) .and. (del .ge. aimsag)) .or.
     &                ((xmult .lt. 1.d0) .and. (del .le. aimsag))) then
                       xmult = 1.d0/xmult
                    endif

 850                continue
                    tmult = xmult*tmult
	            
                 enddo
   
c here it didn't converge.  This appears to happen only if the function is not 
c monotonic.  In that case, go back to the best result found...
                 tmult = tbest

                 do j = 1, npar
                    x(j) = x0(j) + umat(j,ipar)*tmult
                 enddo
                 CALL FCALL(FS1,amin4,x0best,FCN,FUTIL)

                 write(isyswr,860) FS1-amin, ipar, idir
 860             format(1x,'ITERATE: poor convergence warning -- ',
     &                'Delta Chisqr=',e12.5,' ipar=',i2,' idir=',i1)

                 if(ldebug) write(isyswr,870) 
     &                tbestl, dbestl, tbesth, dbesth, tmult, fs1-amin
 870             format(11x,'tmult=',e16.8,' ==> del=',e15.8,
     &                /11x,'tmult=',e16.8,' ==> del=',e15.8,
     &                /7x,'use tmult=',e16.8,' ==> del=',e15.8)

 880             continue


                 if(ldebug) write(isyswr,890) itry, ipar, idir
 890      format(1x,'alternative method:  itry=',i3,' ipar=',i3,' idir=',i2)

                 tmvec(idir) = tmult

              enddo		!end of loop over forward/back

c get an "average" tmult, based on assuming a quadratic + cubic form...
c NOTE, tmvec(2) is negative.
              if(igood .eq. 3) then
                 tmp = 1.d0/tmvec(1)**2
     $                + 1.d0/tmvec(2)**2 + 1.d0/(tmvec(1)*tmvec(2))
                 tmult = 1.d0/sqrt(tmp)

                 if(ldebug) write(isyswr,2230) tmvec(1), tmvec(2),tmult 
 2230            format(1x,'tmvec(1), tmvec(2), tmult=',3e13.5)

              elseif(igood .eq. 1) then
                 tmult = abs(tmvec(1))
                 if(ldebug) write(isyswr,2240) tmult
 2240            format(1x,'only idir=1 converged; use tmult=',e14.6)

              elseif(igood .eq. 2) then
                 tmult = abs(tmvec(2))
                 if(ldebug) write(isyswr,2250) tmult
 2250            format(1x,'only idir=2 converged, use tmult=',e14.6)

              else
                 tmult = 0.001d0
                 if(ldebug) write(isyswr,2260) tmult
 2260   format(1x,'warning: neither direction converged; take tmult=',
     &                f10.5)
        
              endif


c update umat...
              do j = 1, npar
                 umat(j,ipar) = umat(j,ipar)*tmult/sqrt(aimsag)
              enddo

C     print *,'TMUL STD',tmult,tmvec(1)-1,tmvec(2)-1
c ============================
c end of alternative procedure
c ============================
 900          continue
           
           
              if(ldebug) then
                 sum = 0.d0
                 do j = 1, npar
                    sum = sum + umat(j,ipar)**2
                 enddo
                 write(isyswr,920) ipar, sqrt(sum*aimsag), aimsag
 920             format(1x,'Eigenvector #',i2,1x,'has length=',e12.5,
     &                ' to increase fcn by',e12.5)

                 do j = 1, npar
                    JEXT = NEXOFI(j)
                    write(isyswr,930) j, jext, umat(j,ipar)*sqrt(aimsag)
 930                format(5x,'par',i3,' TruePar',i3,1x,g16.8)
                 enddo
              endif


           enddo		!end of loop over eigenvector directions ipar
        endif  ! fortran instead of C
        
	IF(ldebug .and. (amin4 .lt. min(amin,amin4s)))
     &	   WRITE (ISYSWR,950) amin4, amin4s-amin4
950	   format(1x,'new minimum found while adjusting lengths',
     &       e15.7,' lower by',e12.5,' (not used now)')


c See how much have the eigenvector directions have changed from the previous iteration...
	if(iter .ge. 1) then
	   do i = 1, npar
	      sum = 0.d0
	      do j = 1, npar
	         sum = sum + vmat(j,i)*vold(j,i)
	      enddo
	      if(ldebug) write(isyswr,960) i, sum, iter
960	      format(1x,'Dot product between old and new directions',
     &             i3,' =',f10.8,' iter=',i4)
	   enddo
	endif

c replace the location of the minimum by any accidentally discovered improvement...
	if(amin4 .lt. amin) then
	   if(ldebug) write(isyswr,970) amin, amin4, amin-amin4
970	   format(1x,'replace FMIN=',e15.8,
     &          ' by accidentally discovered',e15.8,' lower by',e12.5)
	   amin = amin4
	   do k = 1, npar
	      if(ldebug) write(isyswr,980) k, x0(k), x0best(k)
980	      format(1x,'old x0(',i2,') =',e15.8,2x,'new =',e15.8)
	      x0(k) = x0best(k)
	   enddo
	endif

c set Minuit's vector to our current best notion of the minimum...
	do i = 1, npar
	   x(i) = x0(i)
	enddo
	call mninex(x)

c ---------------------------------------------------------------------------
c ---------------------------------------------------------------------------
c write the current version of umat(i,j)...
c *** commented out for distribution version ***
c 	call wrumat
c ---------------------------------------------------------------------------
c ---------------------------------------------------------------------------

	write(isyswr,990) iter, amin, nfcn
990	format(1x,'End of iteration',i4,' FCN=',g16.9,' NFCN=',i6)

1000	continue

        if(ldebug) then
           write(isyswr,1010)
1010       format(1x,'Final parameter values in ITERATE:',
     &           /3x,'Par',2x,'TruePar',7x,'x')
           do i = 1, npar
              write(isyswr,1020) i, nexofi(i), x(i)
1020          format(3x,i2,5x,i2,3x,e15.7)
           enddo
        endif

c calculate the final notion of the hessian and the error matrix...
	do i = 1, npar
	do j = 1, i
	   sum1 = 0.d0
	   sum2 = 0.d0
	   do k = 1, npar
	      tmp = tvec(k)**2 / 2.d0
	      sum1 = sum1 + vmat(i,k)*vmat(j,k)/tmp
	      sum2 = sum2 + vmat(i,k)*vmat(j,k)*tmp
	   enddo

	   hquad(i,j) = sum1
	   hquad(j,i) = sum1

	   errmat(i,j) = sum2
	   errmat(j,i) = sum2

	enddo
	enddo

	write(isyswr,1100) amin, nfcn
1100	format(1x,'ITERATE finished:  FCN=',g16.9,' NFCN=',i6)
	goto 1300

1200	continue
	CSTATU = 'FAILED    '
	WRITE (ISYSWR,'(A)') '  ITERATE FAILS '

1300	continue

	RETURN
	END



      subroutine set_x_values(Xinput)
#include "./d506dp.inc"
#include "./d506cm.inc"
      double precision xinput(*)
      do i=1, npar
         X(i) = Xinput(i)
      enddo
      end

c---------------------------------------------------------
	SUBROUTINE FCALL(FSHERE,amin4,x0best,FCN,FUTIL)
c ---------------------------------------------------------
c call the chisqr function, recording any new best value and its location.
#include "./d506dp.inc"
#include "./d506cm.inc"
	EXTERNAL FCN,FUTIL
	dimension x0best(MNI)

	CALL MNINEX(X)

	NPARX = NPAR
	CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
	NFCN = NFCN + 1

	IF(FS1 .LT. amin4) THEN
	   amin4 = FS1
	   DO K = 1, NPAR
	      X0BEST(K) = X(K)
	   ENDDO
	ENDIF

	FSHERE = FS1

	RETURN
	END

c ---------------------------------------------------------
	subroutine getvt
c ---------------------------------------------------------
c write umat as an orthogonal matrix times a diagonal matrix...
#include "./d506dp.inc"
#include "./d506cm.inc"
	common /umatco/ hquad(MNI,MNI), umat(MNI,MNI), vmat(MNI,MNI), 
     &	                tvec(MNI), errmat(MNI,MNI)

	LOGICAL LDEBUG
	LDEBUG = (IDBG(3) .GE. 1)

	do i = 1, npar
	   sum = 0.d0
	   do j = 1, npar
	      sum = sum + umat(j,i)**2
	   enddo
	   tvec(i) = sqrt(sum)
	enddo

	do i = 1, npar
	do j = 1, npar
	   vmat(i,j) = umat(i,j) / tvec(j)
	enddo
	enddo

c check that vmat is orthogonal...
	errmax = 0.d0
	do i = 1, npar
	do j = 1, npar
	   sum1 = 0.d0
	   sum2 = 0.d0
	   do k = 1, npar
	      sum1 = sum1 + vmat(k,i)*vmat(k,j)
	      sum2 = sum2 + vmat(i,k)*vmat(j,k)
	   enddo

	   if(i .ne. j) then
	      errmax = max(errmax,abs(sum1),abs(sum2))
	   else
	      errmax = max(errmax,abs(sum1-1.d0),abs(sum2-1.d0))
	   endif

	enddo
	enddo

	if(abs(errmax) .gt. 1.e-5) then
	   write(isyswr,100) errmax
100	   format(1x,'getvt warning: vmat not orthogonal -- zero=',
     &          e12.5)
	endif

	return
	end


      double precision function GetUmat(i,j)
C
C  Wraper around Umat common.
C  Return values from UMAT matrix
C

#include "./d506dp.inc"
#include "./d506cm.inc"
      integer i,j
      double precision hquad,umat,vmat,tvec,errmat
C--------------------------------------------------------------------
      common /umatco/ hquad(MNI,MNI), umat(MNI,MNI), vmat(MNI,MNI), 
     &	                tvec(MNI), errmat(MNI,MNI)
C--------------------------------------------------------------------

      GetUmat = umat(i,j)

      end

      subroutine ScaleUmat(i,j,scale)
C
C  Wraper around Umat common.
C  Return values from UMAT matrix
C

#include "./d506dp.inc"
#include "./d506cm.inc"
      integer i,j
      double precision scale
      double precision hquad,umat,vmat,tvec,errmat
C--------------------------------------------------------------------
      common /umatco/ hquad(MNI,MNI), umat(MNI,MNI), vmat(MNI,MNI), 
     &	                tvec(MNI), errmat(MNI,MNI)
C--------------------------------------------------------------------

      umat(i,j) = umat(i,j)*scale

      end


      
