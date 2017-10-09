c-------------
      subroutine fillvfngrid
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'
      INCLUDE 'CONSTCOM.' 
      INCLUDE 'PDFCOM.' 
 
      common  /forqgspl/ xb0,q2,xlog,an,an2,an3,kn,iqn,ixn,isn,ihqn

!  Set up the boundary conditions for the strong coupling evolution 
!  using the 3-flavour strong coupling at the scale of 4-flavour matching
!  stored in the grid   
        alphas0=xqg(0,0.1d0,vfnth(4)**2,0)
        q20alphas=vfnth(4)**2

      do is=-nsmgrid,nspgrid
!  Take the 3-flavour PDFs as an input for the matching conditions
        isch0=0
!  filling the LO 4-flavour PDFs 
        isch1=-2
        Q2=0.04*exp(exp(sgrid(is,isch1))*log(q2ini(isch1)/0.04))
        AN=ALPHAS_ffn4(q2)/4./PI 
        an2=an**2
        call fillvfx(is,8,0,isch0,isch1)
!  filling the NLO 4-flavour PDFs 
        if (kordhq.ge.1) then 
          isch1=-3
          call fillvfx(is,8,1,isch0,isch1)
        end if 
!  Take the LO 4-flavour PDFs as an input for the matching conditions
        isch0=-2
!  filling the LO 5-flavour PDFs 
        isch1=-4
        Q2=0.04*exp(exp(sgrid(is,isch1))*log(q2ini(isch1)/0.04))
        AN=ALPHAS_ffn5(q2)/4./PI 
        an2=an**2
        call fillvfx(is,10,0,isch0,isch1)
        if (kordhq.ge.1) then 
!  Take the NLO 4-flavour PDFs as an input for the matching conditions
          isch0=-3
!  filling the NLO 5-flavour PDFs 
          isch1=-5
          call fillvfx(is,10,1,isch0,isch1)
        end if 
      end do
 
      return 
      end
c-------------
      subroutine fillvfx(is,nhq,kome,isch0,isch1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'
      INCLUDE 'CONSTCOM.'
      INCLUDE 'PDFCOM.' 

      common  /forqgspl/ xb0,q2,xlog,an,an2,an3,kn,iqn,ixn,isn,ihqn

      real*8 fsp(nxtot),bs(nxtot),cs(nxtot),ds(nxtot),xx(nxtot)

      an2=an**2
      ischem=isch1

      do ix=-nxmgrid,nxpgrid
        xx(ix+nxmgrid+1)=xgrid(ix)
      end do

      do IX=-nxmgrid,nxpgrid-1
! FOPT matching 
        Y(ischem,0,IX,is)=an*4*pi
        do iq=1,nhq
          Y(ischem,iq,IX,is)=hqpdf(XGRID(IX),is,iq,nhq,kome,isch0)
        end do 
        Y(ischem,nhq+1,IX,is)=Y(ischem,nhq,IX,is)
      end do
      Y(ischem,0,nxpgrid,is)=Y(ischem,0,nxpgrid-1,is)

      do iq=1,nhq+1
        Y(ischem,iq,nxpgrid,is)=0D0
        do ix=-nxmgrid,nxpgrid 
          fsp(ix+nxmgrid+1)=y(ischem,iq,ix,is)
        end do
        call spline (nxmgrid+nxpgrid+1,xx,fsp,bs,cs,ds)
        do ix=-nxmgrid,nxpgrid
          bcoeff(ischem,iq,ix,is)=bs(ix+nxmgrid+1)
          ccoeff(ischem,iq,ix,is)=cs(ix+nxmgrid+1)
          dcoeff(ischem,iq,ix,is)=ds(ix+nxmgrid+1)
        end do
      end do

      return 
      end
C------------------
      real*8 function hqpdf(xb,is,iq,ihq,kome,isch0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'
      include 'CONSTCOM.'
      include 'PRECCOM.'

      COMMON  /FORQGSPL/ XB0,Q2,XLOG,AN,AN2,an3,kn,iqn,ixn,isn,ihqn
      common /forhqpdf/ r,kome0
      real*8 q(nflim)

      external hqpdfi1,hqpdfi2

      xb0=xb
      iqn=iq
      isn=is
      ihqn=ihq
      kn=isch0

      kome0=kome
      r=q2/rmass(ihq)**2  

      hqpdf=0.

c Local terms for the gluon and quark pieces

      if (iq.eq.1) then 
        hqpdf=1 + an*ome_gg_1_local(xb,r)
        if (kome.ge.1) hqpdf=hqpdf + an2*ome_gg_2_local(xb,r)
        hqpdf=hqpdf*XQGX1(kn,1,XB,IS) 
      end if

      if (iq.ge.2.and.iq.le.ihq-1) then 
        hqpdf=1
        if (kome.ge.1) hqpdf=hqpdf + an2*ome_qqns_2_local(xb,r)
        hqpdf=hqpdf*XQGX1(kn,iq,XB,IS)  
      end if

c  integration of the total regular term 

      if (xb.ge.0.1) then  
        CALL GAUSS1(hqpdfi1
     ,            ,log(1d-8),log(1.-xb),nmthq,res,EPS)
      else 
        CALL GAUSS1(hqpdfi1
     ,            ,log(1d-8),log(0.9d0),nmthq,df1,EPS1)
        CALL GAUSS1(hqpdfi2,log(xb)
     ,            ,log(0.1d0),nmthq,df2,EPS2)
        res=df1+df2
        eps=sqrt(eps1**2+eps2**2)
      end if 

      hqpdf=hqpdf+res

      return 
      end
C------------------
      real*8 function hqpdfi1(t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      y=1.-exp(t)

      hqpdfi1=hqpdfi(y)*(1-y)

      return 
      end
C------------------
      real*8 function hqpdfi2(t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      y=exp(t)

      hqpdfi2=hqpdfi(y)*y

      return 
      end
C------------------
      real*8 function hqpdfi(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      real*8 q(13)

      common  /forqgspl/ xb0,q2,xlog,an,an2,an3,kn,iqn,ixn,isn,ihqn
      common /forhqpdf/ r,kome0

      y=xb0/z
      hqpdfi=0.

c   gluon distribution 

      if (iqn.eq.1) then 
        if (kome0.ge.1) then 
          glu0=xqgx1(kn,1,xb0,isn)   
          glu=xqgx1(kn,1,y,isn)      
          hqpdfi=hqpdfi+ome_gg_2_singular(z,r)*(glu-glu0)
          hqpdfi=hqpdfi+ome_gg_2(z,r)*glu
          qps=0.
          do k=2,ihqn-1
            qps=qps+xqgx1(kn,k,y,isn)   
          end do
          hqpdfi=hqpdfi+ome_gq_2(z,r)*qps
          hqpdfi=hqpdfi*an2
        end if
      end if

c  light quark disributions

      if (iqn.ge.2.and.iqn.le.ihqn-1) then 
        if (kome0.ge.1) then 
          pdf0=xqgx1(kn,iqn,xb0,isn)   
          pdfc=xqgx1(kn,iqn,y,isn)      
          hqpdfi=hqpdfi+ome_qqns_2_singular(z,r)*(pdfc-pdf0)
          hqpdfi=hqpdfi+ome_qqns_2(z,r)*pdfc
          hqpdfi=hqpdfi*an2
        end if
      end if

c  heavy quark distributions

      if (iqn.ge.ihqn) then 
        glu=xqgx1(kn,1,y,isn)   
        hqpdfi=hqpdfi+ome_g_1(z,r)*an*glu
        if (kome0.ge.1) then 
          qps=0.
          do k=2,ihqn-1
            qps=qps+xqgx1(kn,k,y,isn)   
          end do
          hqpdfi=hqpdfi + an2*(ome_g_2(z,r)*glu + ome_q_2(z,r)*qps)
        end if
        hqpdfi=hqpdfi/2.
      end if

      return 
      end
C------------------
      real*8 function XQGX1(k,iq,XB,IS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'

      data xbsave /-1d0/

      save xbsave,ix,aa
      
      if (xb.ne.xbsave) then 
        do ix=-nxmgrid,nxpgrid-1
          if (xb.lt.xgrid(ix+1)) goto 300
        end do
 300    xbsave=xb
      end if

      aa=xb-xgrid(ix)

      if (khalf.eq.0) then 
        xqgx1=y(k,iq,ix,is)+aa*bcoeff(k,iq,ix,is)
     +     +aa**2*ccoeff(k,iq,ix,is)+aa**3*dcoeff(k,iq,ix,is)
      else 
        xqgx1=yhalf(k,iq,ix)+aa*bcoeffh(iq,ix)
     +     +aa**2*ccoeffh(iq,ix)+aa**3*dcoeffh(iq,ix)
      end if

      RETURN
      END
