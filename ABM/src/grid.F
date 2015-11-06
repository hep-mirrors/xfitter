C-----------------
      SUBROUTINE PDFFILLGRID
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'
      include 'CONSTCOM.'
      INCLUDE 'PDFCOM.'

      real*8 fsp(nxtot),bs(nxtot),cs(nxtot),ds(nxtot),xx(nxtot)
      real*8 pdfs(0:8),dpdfs(0:8,27)
c RP      
      dimension pdfsff(-6:6)

      real*8 f(-6:6)

!  Setup of the PDF grid

      smax=log(log(q2max/0.04)/log(q2ini(0)/0.04))
      smin=-log(log(q2min/0.04)/log(q2ini(0)/0.04))
      dels1(kschemepdf)=smax/nspgrid
      dels2(kschemepdf)=smin/nsmgrid
      xlog0=log(xbmin)
      delx2=-(xlog0-xlog1)/nxmgrid

      delx3=(log(1-xbmax)-log(1-x1))/(nxpgrid-1)

! Mapping with log(1-x) at large x
      do i=0,nxpgrid-1
        xgrid(i)=1-exp(log(1-xbmax)+delx3*(i-nxpgrid+1))
      end do
c 
      xgrid(nxpgrid)=1.
! Mapping with log(x) at small x
      do i=-nxmgrid,0
        xgrid(i)=exp(xlog1+delx2*i)
      end do

      DO I=0,nspgrid
        SGRID(I,kschemepdf)=DELS1(kschemepdf)*I
      END DO

      DO I=-nsmgrid,0
        SGRID(I,kschemepdf)=DELS2(kschemepdf)*I
      END DO


! ...then fill the PDF grid....

      DO is=-nsmgrid,nspgrid
        q2=exp(exp(sgrid(is,kschemepdf))
     *        *log(q2ini(kschemepdf)/0.04))*0.04
        DO IX=-nxmgrid,nxpgrid
          xb=xgrid(ix)


c get PDFs ...    
          CALL HF_Get_PDFs(xb,q2,PDFSFF)

c ... and alpha_s
          an = HF_Get_alphas(q2)

!  \alpha_s  
          Y(kschemepdf,0,IX,is)=an

!  gluon
          Y(kschemepdf,1,IX,is)=pdfsff(0)
!      write(6,*)'RP as,q2,g',an,q2,Y(kschemepdf,1,IX,is)
!          write(6,*)'gluon:',Y(kschemepdf,1,IX,is)
          
!         if(q2.lt.5.) then
!      print*,'x,q2',xb,q2     
!      write(6,*)'RP alpha,g',Y(kschemepdf,0,IX,is),Y(kschemepdf,1,IX,is) 
!         endif
!  d
          Y(kschemepdf,2,IX,is)=pdfsff(1)
!  dbar
          Y(kschemepdf,3,IX,is)=pdfsff(-1)
!  u
          Y(kschemepdf,4,IX,is)=pdfsff(2)
!  ubar
          Y(kschemepdf,5,IX,is)=pdfsff(-2)
!  s,sbar
          Y(kschemepdf,6,IX,is)=pdfsff(3)
          Y(kschemepdf,7,IX,is)=pdfsff(-3)
          Y(kschemepdf,8,IX,is)=0.
          Y(kschemepdf,9,IX,is)=0.
          Y(kschemepdf,10,IX,is)=0.
          Y(kschemepdf,11,IX,is)=0.
!  c,cbar
          if (kschemepdf.ge.1) then
            Y(kschemepdf,8,IX,is)=pdfsff(-4)          
            Y(kschemepdf,9,IX,is)=pdfsff(-4)          
            Y(kschemepdf,10,IX,is)=0.
            Y(kschemepdf,11,IX,is)=0.
          end if
!  b,bbar
          if (kschemepdf.ge.2) then
            Y(kschemepdf,10,IX,is)=pdfsff(-5)          
            Y(kschemepdf,11,IX,is)=pdfsff(-5)          
          end if

        end do

! ...and finally, construct the spline interpolation coefficients for the grid.
        do iq=0,7+kschemepdf*2
          do ix=-nxmgrid,nxpgrid
            xx(ix+nxmgrid+1)=xgrid(ix)
            fsp(ix+nxmgrid+1)=y(kschemepdf,iq,ix,is)
          end do
          call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
          do ix=-nxmgrid,nxpgrid
              bcoeff(kschemepdf,iq,ix,is)=bs(ix+nxmgrid+1)
              ccoeff(kschemepdf,iq,ix,is)=cs(ix+nxmgrid+1)
              dcoeff(kschemepdf,iq,ix,is)=ds(ix+nxmgrid+1)
          end do
        end do
      end do

      RETURN
      END
C------------------
      REAL*8 FUNCTION XQG(IQ,xx,Q2,kp)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'

      DATA Q2SAV/-1./
      data xbsav/-1./
      save xbsav,q2sav,is,q,ix,p

      xb=xx               

      if(q2.gt.q2max.or.q2.lt.q2min)
     -            print*,'q2 out of range in xqg:',q2,q2min,q2max
      s=log(log(q2/0.04)/log(q2ini(kp)/0.04))
      if (s.ge.0.) then
        is=int(s/dels1(kp))
        dss=(s-sgrid(is,kp))/dels1(kp)
      else
        is=int(s/dels2(kp))-1
        dss=(s-sgrid(is,kp))/dels2(kp)
      end if
      q=(s-sgrid(is,kp))/(sgrid(is+1,kp)-sgrid(is,kp))

      if(xb.ne.xbsav) then 
        if (xb.gt.1d0+1d-8.or.xb.lt.xbmin)
     -              print *,'xb out of range in xqg:',xb

        do ix=-nxmgrid,nxpgrid-1
          if (xb.lt.xgrid(ix+1)) goto 300
        end do
 300  xbsav=xb
      end if

      aa=xb-xgrid(ix)

      if(is.ne.-nsmgrid.and.ix.ne.-nxmgrid) then 
        f0=y(kp,iq,ix,is)+aa*bcoeff(kp,iq,ix,is)
     +    +aa**2*ccoeff(kp,iq,ix,is)+aa**3*dcoeff(kp,iq,ix,is)
        fp=y(kp,iq,ix,is+1)+aa*bcoeff(kp,iq,ix,is+1)
     +    +aa**2*ccoeff(kp,iq,ix,is+1)+aa**3*dcoeff(kp,iq,ix,is+1)
        fm=y(kp,iq,ix,is-1)+aa*bcoeff(kp,iq,ix,is-1)
     +    +aa**2*ccoeff(kp,iq,ix,is-1)+aa**3*dcoeff(kp,iq,ix,is-1)
        xqg=fm*q*(q-1d0)/2d0 + f0*(1d0-q**2) + fp*q*(q+1d0)/2d0
      else
        f0=y(kp,iq,ix,is)+aa*bcoeff(kp,iq,ix,is)
     +    +aa**2*ccoeff(kp,iq,ix,is)+aa**3*dcoeff(kp,iq,ix,is)
        fp=y(kp,iq,ix,is+1)+aa*bcoeff(kp,iq,ix,is+1)
     +    +aa**2*ccoeff(kp,iq,ix,is+1)+aa**3*dcoeff(kp,iq,ix,is+1)
        xqg=f0*(1d0-q) + fp*q
      end if

      RETURN
      END
!------------------
      real*8 function xqgn(iq,xb,q2,kp)
      implicit double precision (a-h,o-z)

!  The neutron PDFs

      k=iq
      if (iq.eq.2) k=4
      if (iq.eq.3) k=5
      if (iq.eq.4) k=2
      if (iq.eq.5) k=3

      xqgn=xqg(k,xb,q2,kp)

      RETURN
      END
