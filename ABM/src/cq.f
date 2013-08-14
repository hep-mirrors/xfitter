c-------------
      subroutine fillvfngrid
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'
      INCLUDE 'CONSTCOM.'

      common  /forqgspl/ xb0,q2,xlog,an,an2,an3,kn,iqn,ixn,isn

!  Set up the boundary conditions for the strong copuling evolution 
!  using the 3-flavour strong coupling at the scale of the c-quark mass
!  stored in the grid   
      alphas0=xqg(0,0.1d0,rmass(8)**2,0)
      q20alphas=rmass(8)**2

      do is=-nsmgrid,nspgrid
        Q2=0.04*exp(exp(sgrid(is))*log(q2ini/0.04))
!  generation of the LO 4-flavour PDFs
        AN=ALPHAS_ffn4(q2)/4./PI 
        an2=an**2
        call fillvfx(is,8,0)
!  generation of the NLO 4-flavour PDFs
        if (kordhq.ge.1) then 
          call fillvfx(is,8,1)
        end if 
!  generation of the LO 5-flavour PDFs
        AN=ALPHAS_ffn5(q2)/4./PI 
        an2=an**2
        call fillvfx(is,10,0)
!  generation of the NLO 5-flavour PDFs
        if (kordhq.ge.1) then 
          call fillvfx(is,10,1)
        end if 
      end do

      return 
      end
c-------------
      subroutine fillvfx(is,nhq,kome)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'
      INCLUDE 'CONSTCOM.'

      common  /forqgspl/ xb0,q2,xlog,an,an2,an3,kn,iqn,ixn,isn

      real*8 fsp(nxtot),bs(nxtot),cs(nxtot),ds(nxtot),xx(nxtot)

!  Set ischem=-2 for the LO 4-flavour 
!      ischem=-3 for the NLO 4-flavour
!      ischem=-4 for the LO 5-flavour
!      ischem=-5 for the NLO 5-flavour PDFs

      ischem=-(nhq-8)-kome-2
      do ix=-nxmgrid,nxpgrid
        xx(ix+nxmgrid+1)=xgrid(ix)
      end do

      DO IX=-nxmgrid,nxpgrid-1
        Y(ischem,0,IX,is)=an*4*pi
        do iq=1,nhq
          Y(ischem,iq,IX,is)=hqpdf(XGRID(IX),is,iq,nhq,kome)
        end do 
        Y(ischem,nhq+1,IX,is)=Y(ischem,nhq,IX,is)
      end do
        Y(ischem,IQ,nxpgrid,is)=an*4*pi
      do iq=1,nhq+1
        Y(ischem,IQ,nxpgrid,is)=0D0
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
      real*8 function hqpdf(xb,is,iq,ihq,kome)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'
      include 'CONSTCOM.'
      include 'PRECCOM.'

      COMMON  /FORQGSPL/ XB0,Q2,XLOG,AN,AN2,an3,kn,iqn,ixn,isn
      common /forhqpdf/ r,kome0
      real*8 q(nflim)

      external hqpdfi1,hqpdfi2

      xb0=xb
      iqn=iq
      isn=is
!  Take the 3-flavour PDFs as an input for the matching conditions
!  both for 4- and 5-flavour PDFs generation
      kn=0
      kome0=kome
      r=q2/rmass(ihq)**2  

      hqpdf=0.

c Local terms for the gluon and quark pieces

      if (iq.eq.1) then 
        hqpdf=1 + an*ome_gg_1_local(xb,r)
        if (kome.ge.1) hqpdf=hqpdf + an2*ome_gg_2_local(xb,r)
        hqpdf=hqpdf*XQGX1(kn,1,XB,IS) 
      end if

      if (iq.ge.2.and.iq.le.7+2*kn) then 
        hqpdf=1
        if (kome.ge.1) hqpdf=hqpdf + an2*ome_qqns_2_local(xb,r)
        hqpdf=hqpdf*XQGX1(kn,iq,XB,IS)  
      end if

c  integration of the total regular term 

        CALL GAUSS1(hqpdfi2,log(xb),0d0,nmthq,res,EPS2)

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

      common  /forqgspl/ xb0,q2,xlog,an,an2,an3,kn,iqn,ixn,isn
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
          do k=2,7+2*kn
            qps=qps+xqgx1(kn,k,y,isn)   
          end do
          hqpdfi=hqpdfi+ome_gq_2(z,r)*qps
          hqpdfi=hqpdfi*an2
        end if
      end if

c  light quark disributions

      if (iqn.ge.2.and.iqn.le.7+2*kn) then 
        if (kome0.ge.1) then 
          pdf0=xqgx1(kn,iqn,xb0,isn)   
          pdfc=xqgx1(kn,iqn,y,isn)      
          hqpdfi=hqpdfi+ome_qqns_2_singular(z,r)*(pdfc-pdf0)
          hqpdfi=hqpdfi+ome_qqns_2(z,r)*pdfc
          hqpdfi=hqpdfi*an2
        end if
      end if

c  heavy quark distributions

      if (iqn.ge.8+2*kn) then 
        glu=xqgx1(kn,1,y,isn)   
        hqpdfi=hqpdfi+ome_g_1(z,r)*an*glu
        if (kome0.ge.1) then 
          qps=0.
          do k=2,7+2*kn
            qps=qps+xqgx1(kn,k,y,isn)   
          end do
          hqpdfi=hqpdfi + an2*(ome_g_2(z,r)*glu + ome_q_2(z,r)*qps)
        end if
        hqpdfi=hqpdfi/2.
      end if

      return 
      end
