!------------------
      real*8 function flcharm_ffn(xb,q2,nq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'
      include 'PRECCOM.'

      external flcharmi
      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0

      if (kordhq.ge.3) 
     -   print *,'FLCHARM_FFN IS NOT READY FOR KORDHQ=',kordhq

!  Save the current scheme status
      kschemepdfs=kschemepdf
      kschemepdf=0   

!  The mass and charge of heavy quark
      rm2=rmass(nq)**2
      qq=rcharge(nq)

!  The factorization scale
      rmu2=q2*hqscale1+4*rm2*hqscale2

!  Take strong coupling constant in the 3-flavour scheme
      an=xqg(0,0.1d0,rmu2,kschemepdf)
      xb0=xb
      q2ss=q2
      qqs=qq

      a=1./(1.+4.*rm2/q2)
      if (xb.lt.a) then
        CALL GAUSS1(flcharmi,dlog(xb),dlog(a),nflhq,flc,EPS)
        flcharm_ffn=flc
      else
        flcharm_ffn=0.
      end if

!  Restore the initial scheme
      kschemepdf=kschemepdfs

      return
      end
!------------------
      real*8 function flcharmi(t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'
      include 'PRECCOM.'

      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0
!  The decoupling constants for nf=3 
      data d1dec,d1dec2, d2dec, beta0
     -    /1.33333333d0, 1.77777778d0, 10.3192955234169d0, 2.25d0/

      z=dexp(t)
      y=xb0/z

      xi=q2ss/rm2
      sp=q2ss*(1./z-1.)
      eta=sp/4./rm2-1.

!  The number of fermions for the strong coupling constant.
      nfe=nfloops(q2ss,kschemepdf)

      rfac0=1.
      rfac1=1.
      alr=-log(rscale)
      bet0=(11-2./3.*nfe)*alr
      bet1=(102-38/3.*nfe)*alr

      if (kordhq.eq.1) then
        rfac0=1-an/4./pi*bet0
      end if
      if (kordhq.eq.2) then
        rfac0=1-an/4./pi*bet0-(an/4./pi)**2*(bet1-bet0**2)
        rfac1=1-an/4./pi*2*bet0
      end if

      glu0=xqg(1,y,rmu2,kschemepdf)
      flcharmi=cllog(eta,xi)*glu0*an/4./pi**2*q2ss/rm2*qqs**2*rfac0 

      if (kordhq.ge.1) then 
        flcharmnloi=clnlog(eta,xi)+clnlobarg(eta,xi)*log(rmu2/rm2)
!  The term appearing for the MS-bar definition of the heavy-quark mass
        if (msbarm) then 
          glubornd=dborn_l(eta,xi)
          flcharmnloi=flcharmnloi+glubornd/xi*d1dec/2./pi**2
        end if

        flcharmnloi=flcharmnloi*qqs**2*glu0

        flcharmnloi=flcharmnloi+qqs**2*xpscqqpm(y,rmu2)
     *   *(clnloq(eta,xi)+clnlobarq(eta,xi)*log(rmu2/rm2))

!  The non-singlet contribution 
        if (hqnons) then 
          flcharmnloi=flcharmnloi+dlnloq(eta,xi)*f2cqqpm(3,1,22,y,rmu2)
        end if 

        flcharmi=flcharmi+flcharmnloi*an**2/pi*q2ss/rm2*rfac1 
      end if

      if (kordhq.ge.2) then 
        call clhqg21(eta,xi,cc21)
        call clhqg22(eta,xi,cc22)

        flcharmnnloi=qqs**2*glu0
     *   *(cc21*log(rmu2/rm2) + cc22*(log(rmu2/rm2))**2)
!  The terms appearing for the MS-bar definition of the heavy-quark mass
        if (msbarm) then 
          call clhqg20m(eta,xi,cc20m)
          call clhqg21m(eta,xi,cc21m)
          call clhqg20p(eta,xi,cc20p)
          call clhqg21p(eta,xi,cc21p)
          ypl=min(y*(1+delder),1d0)
          delder0=ypl/y-1.
          ymn=y*(1-delder0)
          glup=xqg(1,ypl,rmu2,kschemepdf)
          glum=xqg(1,ymn,rmu2,kschemepdf)
          glud=(glup-glum)/(2.*delder0)

          flcharmnnloi=flcharmnnloi + qqs**2
     *   * (glu0*(cc20m + cc21m*log(rmu2/rm2)
     +   + (glubornd*(2*d2dec - d1dec2
     +   + 4*beta0*d1dec*log(rscale*rmu2/rm2)
     -   - 2*d1dec*beta0*log(rmu2/rm2))
     +   + 2*d1dec*beta0*xi*cllog(eta,xi))/(4*pi**2)**2/xi)
     +   + glud*(cc20p + cc21p*log(rmu2/rm2)))
        end if

        flcharmi=flcharmi + flcharmnnloi*an**3*4*xi
      end if

      return
      end
C------------------
      real*8 function f2charm_ffn(xb,q2,nq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'
      include 'PRECCOM.'

      external f2charmi
      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0

!  Save the current scheme status
      kschemepdfs=kschemepdf
      kschemepdf=0

      if (kordhq.ge.3) 
     -   print *,'F2CHARM_FFN IS NOT READY FOR KORDHQ=',kordhq

!  The mass and charge of heavy quark
      rm2=rmass(nq)**2
      qq=rcharge(nq)
!  The factorization scale
      rmu2=q2*hqscale1+4*rm2*hqscale2

!  Take strong coupling constant
      an=xqg(0,0.1d0,rmu2,kschemepdf)
      xb0=xb
      q2ss=q2
      qqs=qq
 
      a=1./(1.+4.*rm2/q2)
      if (xb.lt.a) then
        CALL GAUSS1(f2charmi,log(xb),log(a),nf2hq,f2c,EPS)
        f2charm_ffn=f2c
      else
        f2charm_ffn=0.
      end if

!  Restore the initial scheme
      kschemepdf=kschemepdfs

      return
      end
C------------------
      real*8 function f2charmi(t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'
      include 'PRECCOM.'

      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0
!  The decoupling constants for nf=3 
      data d1dec, d1dec2, d2dec, beta0
     -    /1.33333333d0, 1.77777778d0, 10.3192955234169d0, 2.25d0/

      z=dexp(t)
      y=xb0/z

      xi=q2ss/rm2
      sp=q2ss*(1./z-1.)
      eta=sp/4./rm2-1.

      glu0=xqg(1,y,rmu2,kschemepdf)
      f2charmi=c2log(eta,xi)*glu0*an/4./pi**2*q2ss/rm2*qqs**2

      if (kordhq.ge.1) then 
        f2charmnloi=c2nlog(eta,xi)+c2nlobarg(eta,xi)*log(rmu2/rm2)

!  The term appearing for the MS-bar definition of the heavy-quark mass
        if (msbarm) then 
          glubornd=dborn_l(eta,xi) + dborn_t(eta,xi)
          f2charmnloi=f2charmnloi+glubornd/xi*d1dec/2./pi**2
        end if

        f2charmnloi=f2charmnloi*qqs**2*glu0

        f2charmnloi=f2charmnloi + qqs**2*xpscqqpm(y,rmu2)
     *   *(c2nloq(eta,xi)+c2nlobarq(eta,xi)*log(rmu2/rm2))

!  The non-singlet contribution 
        if (hqnons) then 
          f2charmnloi=f2charmnloi+d2nloq(eta,xi)*f2cqqpm(3,1,22,y,rmu2) !tblt'd
!          f2charmnloi=f2charmnloi                                      !exact 
!     +    +ch2qns10_exact(z,xi)/(16*pi*xi/z)*f2cqqpm(3,1,22,y,rmu2)
        end if 

        f2charmi=f2charmi + f2charmnloi*an**2/pi*q2ss/rm2
      end if

      if (kordhq.ge.2) then 
        call c2hqg20(eta,xi,cc20)
        call c2hqg21(eta,xi,cc21)
        call c2hqg22(eta,xi,cc22)
        f2charmnnloi=qqs**2*glu0
     *   *(cc20 + cc21*log(rmu2/rm2) + cc22*(log(rmu2/rm2))**2)

        if (msbarm) then 
          call c2hqg20m(eta,xi,cc20m)
          call c2hqg21m(eta,xi,cc21m)
          call c2hqg20p(eta,xi,cc20p)
          call c2hqg21p(eta,xi,cc21p)
          ypl=min(y*(1+delder),1d0)
          delder0=ypl/y-1.
          ymn=y*(1-delder0)
          glup=xqg(1,ypl,rmu2,kschemepdf)
          glum=xqg(1,ymn,rmu2,kschemepdf)
          glud=(glup-glum)/(2.*delder0)

          f2charmnnloi=f2charmnnloi + qqs**2
     *   * (glu0*(cc20m + cc21m*log(rmu2/rm2)
     +   + (glubornd*(2*d2dec - d1dec2
     +   + 4*beta0*d1dec*log(rscale*rmu2/rm2)
     -   - 2*d1dec*beta0*log(rmu2/rm2))
     +   + 2*d1dec*beta0*xi*c2log(eta,xi))/(4*pi**2)**2/xi)
     +   + glud*(cc20p + cc21p*log(rmu2/rm2)))
        end if

        f2charmi=f2charmi + f2charmnnloi*an**3*4*q2ss/rm2
      end if

      return
      end
C------------------
      real*8 function ftnucharm(nb,nt,ni,xb,q2,nq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!  The seminclusive neutrino-nucleon charm CC structure function F_T=2xF_1 
!  with account of the QCD corrections up to NLO (cf. PLB 380, 171 (1996))

      include 'CONSTCOM.'
      include 'PDFCOM.'
      include 'PRECCOM.'

      external ftnucharmi
      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0
!  The decoupling constant for nf=3 
      data d1dec /1.33333333d0/

!  Set the 3-flavour scheme
      kschemepdfs=kschemepdf
      kschemepdf=0

      rm2=rmass(8)**2
!  The factorization scale
      rmu2=q2*hqscale1+rm2*hqscale2
      rq=q2/(q2+rm2)
      xi=min(xb/rq,1d0)

!  The number of fermions for the strong coupling constant.
      nfe=nfloops(q2ss,kschemepdf)

!  The factors appearing if the renormalization scale is not equal
!  to the factorization scale rmu2. 

      rfac0=1.
      rfac1=1.
      alr=-log(rscale)
      bet0=(11-2./3.*nfe)*alr
      bet1=(102-38/3.*nfe)*alr

      if (kordhq.ge.1) then
!  The change in the strong coupling scale is taken into account in the grid
        an=xqg(0,0.1d0,rmu2,kschemepdf)/4./pi
        rfac0=1-an*bet0
      end if
      if (kordhq.eq.2) then
        rfac0=1-an*bet0-an**2*(bet1-bet0**2)
        rfac1=1-an*2*bet0
      end if

!  The LO term 
      ftnucharm=2*xb*f2cqqpm(nb,nt,ni,xi,rmu2)/xi

!  The NLO contribution 
      if (kordhq.ge.1.and.xi.lt.1d0) then 
        xb0=xb
        q2ss=q2
        nb0=nb
        nt0=nt
        ni0=ni      
        CALL GAUSS1(ftnucharmi,xi,1d0,nf2hq,f2c,EPS)
        fac=2*xb*2*an
        ak=1./rq*(1-rq)*log(1-rq)      
        dc=cf*(-4-0.5/rq-pi**2/3.-(1+3*rq)/2./rq*ak 
     + +(xi+xi**2/2.+2*log(1-xi))*log((q2+rm2)/rmu2)
     + +0.5*(rq*xi-rq**2*xi+log(1-rq*xi)*(1-rq*xi))/rq**2/(1-rq*xi)
     - -2*(-(log(1-xi))**2+ddilog((1-rq*xi)/(1-rq))-ddilog(1./(1-rq))
     + +log(1-rq*xi)*log(rq*(1-xi)/(1-rq)))  
     - -2*log(1-xi))*f2cqqpm(nb,nt,ni,xi,rmu2)/xi
        ftnucharm=ftnucharm+(f2c+dc)*fac

!  The term of Ref.[Phys.Lett. B699, 345 (2011)] appearing for the 
!  MS-bar definition of the heavy-quark mass 

        if(msbarm) then
          xip=min(xi*(1+delder),1d0)
          delder0=xip/xi-1.
          xim=xi*(1-delder0)
          ftxp=f2cqqpm(nb0,nt0,ni0,xip,rmu2)/xip
          ftxm=f2cqqpm(nb0,nt0,ni0,xim,rmu2)/xim
          dftx=(ftxp-ftxm)/2./delder0/xi
          ftnucharm=ftnucharm + 4*d1dec*dftx*(xi-xb)*fac
        end if 
      end if
      ftnucharm=ftnucharm*rfac0

!  Restore the initial scheme
      kschemepdf=kschemepdfs

      return
      end
!------------------
      real*8 function ftnucharmi(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'

      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0
      data d1dec /1.33333333d0/

      xi=xb0/rq
      z=xi/y

      xin=q2ss/rm2
      eta=xin*(1./(rq*y)-1.)-1.

      rl=log((1-rq*y)/(1-rq)/y)
      rlmf=log((rm2+q2ss)/rmu2)

!  The number of fermions for the strong coupling constant.
      nfe=nfloops(q2ss,kschemepdf)

!  The factors appearing if the renormalization scale is not equal
!  to the factorization scale rmu2. 
      rfac0=1.
      rfac1=1.
      alr=-log(rscale)
      bet0=(11-2./3.*nfe)*alr
      bet1=(102-38/3.*nfe)*alr

      if (kordhq.ge.1) then
        rfac0=1-an*bet0
      end if
      if (kordhq.eq.2) then
        rfac0=1-an*bet0-an**2*(bet1-bet0**2)
        rfac1=1-an*2*bet0
      end if

      glu0=xqg(1,z,rmu2,kschemepdf)

!  The gluon-initiated contribution 
      pqg0=(y**2+(1-y)**2)/2.
      ftg=pqg0*(rl+log((q2ss+rm2)/rmu2))
      ftg=ftg+pqg0*(2*log(1-y)-log(1-rq*y)-log(y))
      ftg=ftg+y*(1-y)*(4-4*(1-rq))
      ftg=ftg+(1-rq)*y/(1-rq*y)-1
      ftg=ftg+(1-rq)*y*rl*(2+rq*y*(-4))

      ftg=ftg*glu0

!  The quark-initiated contribution 
      ftz=f2cqqpm(nb0,nt0,ni0,z,rmu2)
      ftx=f2cqqpm(nb0,nt0,ni0,xi,rmu2)

      ftq=-(1+y**2)*log(y)/(1-y)
      ftq=cf*ftq*ftz

      ftq0=(1+y**2)/(1-y)*log((q2ss+rm2)/rmu2)*(ftz-ftx)
      ftq0=ftq0+(2*log(1-y)-log(1-rq*y))/(1-y)*((1+y**2)*ftz-2*ftx)
      ftq0=ftq0+1./(1-y)*((1-4*y+y**2)*ftz+2*ftx)
      ftq0=ftq0+1./(1-rq*y)*(y-y**2)*ftz
      ftq0=ftq0+0.5*(1-y)/(1-rq*y)**2*(ftz-ftx)
      ftq0=cf*ftq0

!  The total NLO term
      ftnucharmi=(ftq0+ftq+ftg)/xi

!  The asymptotic NNLO terms 
      if (kordhq.ge.2.and.q2ss.ge.50d0) then 
!  Gluon contribution 
        call cthqg20cc_asy(eta,xin,cc20)
        call cthqg21cc_asy(eta,xin,cc21)
        call chqg22cc_asy(eta,xin,cc22)

        ftcharmnnloi=glu0*(cc20 + cc21*rlmf + cc22*rlmf**2)/16.

!  Pure-singlet contribution 
        ps0=xpscqqpm(z,rmu2)

        call cthqps20cc_asy(eta,xin,cc20)
        call cthqps21cc_asy(eta,xin,cc21)
        call chqps22cc_asy(eta,xin,cc22)

        ftcharmnnloi=ftcharmnnloi + ps0*(cc20 
     +   + cc21*rlmf + cc22*rlmf**2)/16.

        if(msbarm) then
          dhg10 = 2-4*y+4*y**2
          ftcharmnnloi=ftcharmnnloi - d1dec/4.*glu0
        end if

        ftnucharmi=ftnucharmi + ftcharmnnloi*8*an/xi*rfac1
      end if

      return
      end
!------------------
      real*8 function f3nucharm(nb,nt,ni,xb,q2,nq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!  The seminclusive neutrino-nucleon charm CC structure function F_3 
!  with account of the QCD corrections up to NLO (cf. PLB 380, 171 (1996))

      include 'CONSTCOM.'
      include 'PDFCOM.'
      include 'PRECCOM.'

      external f3nucharmi
      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0
!  The decoupling constant for nf=3 
      data d1dec /1.33333333d0/

!  Set the 3-flavour scheme
      kschemepdfs=kschemepdf
      kschemepdf=0

      rm2=rmass(8)**2
!  The factorization scale
      rmu2=q2*hqscale1+rm2*hqscale2
      rq=q2/(q2+rm2)
      xi=min(xb/rq,1d0)

!  The number of fermions for the strong coupling constant.
      nfe=nfloops(q2ss,kschemepdf)

!  The factors appearing if the renormalization scale is not equal
!  to the factorization scale rmu2. 
      rfac0=1.
      rfac1=1.
      alr=-log(rscale)
      bet0=(11-2./3.*nfe)*alr
      bet1=(102-38/3.*nfe)*alr

      if (kordhq.ge.1) then
!  The change in the strong coupling scale is taken into account in the grid
        an=xqg(0,0.1d0,rmu2,kschemepdf)/4./pi
        rfac0=1-an*bet0
      end if
      if (kordhq.eq.2) then
        rfac0=1-an*bet0-an**2*(bet1-bet0**2)
        rfac1=1-an*2*bet0
      end if

!  The LO term 
      f3nucharm=2*f3cqqpm(nb,nt,ni,xi,rmu2)/xi
      
!  The NLO contribution 
      if (kordhq.ge.1.and.xi.lt.1d0) then 
        xb0=xb
        q2ss=q2
        nb0=nb
        nt0=nt
        ni0=ni      
        CALL GAUSS1(f3nucharmi,xi,1d0,nf3hq,f3c,EPS)
        fac=2*2*an
        ak=1./rq*(1-rq)*log(1-rq)      
        dc=cf*(-4-0.5/rq-pi**2/3.-(1+3*rq)/2./rq*ak 
     + +(xi+xi**2/2.+2*log(1-xi))*log((q2+rm2)/rmu2)
     - -2*(-(log(1-xi))**2+ddilog((1-rq*xi)/(1-rq))-ddilog(1./(1-rq))
     + +log(1-rq*xi)*log(rq*(1-xi)/(1-rq)))  
     + +0.5*(rq*xi-rq**2*xi+log(1-rq*xi)*(1-rq*xi))/rq**2/(1-rq*xi)
     - -2*log(1-xi))*f3cqqpm(nb,nt,ni,xi,rmu2)/xi
        f3nucharm=f3nucharm+(f3c+dc)*fac

!  The term of Ref.[Phys.Lett. B699, 345 (2011)] appearing for the 
!  MS-bar definition of the heavy-quark mass 

        if(msbarm) then
          xip=min(xi*(1+delder),1d0)
          delder0=xip/xi-1.
          xim=xi*(1-delder0)
          f3xp=f3cqqpm(nb0,nt0,ni0,xip,rmu2)/xip
          f3xm=f3cqqpm(nb0,nt0,ni0,xim,rmu2)/xim
          df3x=(f3xp-f3xm)/2./delder0/xi
          f3nucharm=f3nucharm + 4*d1dec*df3x*(xi-xb)*fac
        end if
      end if
      f3nucharm=f3nucharm*rfac0

!  Restore the initial scheme
      kschemepdf=kschemepdfs

      return
      end
!------------------
      real*8 function f3nucharmi(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'

      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0
      data d1dec /1.33333333d0/

      xi=xb0/rq
      z=xi/y

      xin=q2ss/rm2
      eta=xin*(1./(rq*y)-1.)-1.

      rl=log((1-rq*y)/(1-rq)/y)
      rlmf=log((rm2+q2ss)/rmu2)
      rlm=log((rm2+q2ss)/rm2)

!  The number of fermions for the strong coupling constant.
      nfe=nfloops(q2ss,kschemepdf)

!  The factors appearing if the renormalization scale is not equal
!  to the factorization scale rmu2. 
      rfac0=1.
      rfac1=1.
      alr=-log(rscale)
      bet0=(11-2./3.*nfe)*alr
      bet1=(102-38/3.*nfe)*alr

      if (kordhq.ge.1) then
        rfac0=1-an*bet0
      end if
      if (kordhq.eq.2) then
        rfac0=1-an*bet0-an**2*(bet1-bet0**2)
        rfac1=1-an*2*bet0
      end if

      glu0=xqg(1,z,rmu2,kschemepdf)

!  The gluon-initiated contribution 
      pqg0=(y**2+(1-y)**2)/2.
      f3g=pqg0*(-rl+log((q2ss+rm2)/rmu2))
      f3g=f3g+pqg0*(2*log(1-y)-log(1-rq*y)-log(y))
      f3g=f3g+y*(1-y)*2*(1-rq)
      f3g=f3g+(1-rq)*y*rl*(-2*(1-y)+rq*y*2)
      f3g=f3g*glu0

!  The quark-initiated contribution 
      f3z=f3cqqpm(nb0,nt0,ni0,z,rmu2)
      f3x=f3cqqpm(nb0,nt0,ni0,xi,rmu2)

      f3q=-(1+y**2)*log(y)/(1-y)
      f3q=cf*f3q*f3z

      f3q0=(1+y**2)/(1-y)*log((q2ss+rm2)/rmu2)*(f3z-f3x)
      f3q0=f3q0+(2*log(1-y)-log(1-rq*y))/(1-y)*((1+y**2)*f3z-2*f3x)
      f3q0=f3q0+1./(1-y)*((-1-y**2)*f3z+2*f3x)
      f3q0=f3q0+1./(1-rq*y)*(1-y)*f3z
      f3q0=f3q0+0.5*(1-y)/(1-rq*y)**2*(f3z-f3x)
      f3q0=cf*f3q0

!  The total NLO contribution 
      if (nb0.eq.6) f3nucharmi=(f3q0+f3q+f3g)/xi
      if (nb0.eq.7) f3nucharmi=(f3q0+f3q-f3g)/xi

!  The asymptotic NNLO terms 
      if (kordhq.ge.2.and.q2ss.ge.50d0) then 
!  Gluon contribution 
        call c3hqg20cc_asy(eta,xin,cc20)
        call c3hqg21cc_asy(eta,xin,cc21)
        call chqg22cc_asy(eta,xin,cc22)

        f3charmnnloi=glu0*(cc20 + cc21*rlmf + cc22*rlmf**2)/16.

!  Pure-singlet contribution 
        ps0=xpscqqpm(z,rmu2)

        call c3hqps20cc_asy(eta,xin,cc20)
        call c3hqps21cc_asy(eta,xin,cc21)
        call chqps22cc_asy(eta,xin,cc22)

        f3charmnnloi=f3charmnnloi + ps0*(cc20 
     +   + cc21*rlmf + cc22*rlmf**2)/16.

        if(msbarm) then
          dhg10 = 2-4*y+4*y**2
          f3charmnnloi=f3charmnnloi + d1dec/4.*glu0
        end if

        f3charmnnloi=f3charmnnloi*8*an/xi*rfac1

        if (nb0.eq.6) f3nucharmi=f3nucharmi + f3charmnnloi
        if (nb0.eq.7) f3nucharmi=f3nucharmi - f3charmnnloi


      end if

      return
      end
!------------------
      real*8 function f2nucharm(nb,nt,ni,xb,q2,nq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!  The seminclusive neutrino-nucleon charm CC structure function F_2 
!  with account of the QCD corrections up to NLO (cf. PLB 380, 171 (1996))

      include 'CONSTCOM.'
      include 'PDFCOM.'
      include 'PRECCOM.'

      external f2nucharmi
      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0
c  The decoupling constant for nf=3 
      data d1dec /1.33333333d0/

!  Set the 3-flavour scheme
      kschemepdfs=kschemepdf
      kschemepdf=0

      rm2=rmass(8)**2
!  The factorization scale
      rmu2=q2*hqscale1+rm2*hqscale2
      rq=q2/(q2+rm2)
      xi=min(xb/rq,1d0)

!  The number of fermions for the strong coupling constant.
      nfe=nfloops(q2ss,kschemepdf)

!  The factors appearing if the renormalization scale is not equal
!  to the factorization scale rmu2. 

      rfac0=1.
      rfac1=1.
      alr=-log(rscale)
      bet0=(11-2./3.*nfe)*alr
      bet1=(102-38/3.*nfe)*alr

      if (kordhq.ge.1) then
!  The change in the strong coupling scale is taken into account in the grid
        an=xqg(0,0.1d0,rmu2,kschemepdf)/4./pi
        rfac0=1-an*bet0
      end if
      if (kordhq.eq.2) then
        rfac0=1-an*bet0-an**2*(bet1-bet0**2)
        rfac1=1-an*2*bet0
      end if

!  The LO term 
      f2nucharm=2*f2cqqpm(nb,nt,ni,xi,rmu2)

!  The NLO contribution 
      if (kordhq.ge.1.and.xi.lt.1d0) then 
        xb0=xb
        q2ss=q2
        nb0=nb
        nt0=nt
        ni0=ni      
        CALL GAUSS1(f2nucharmi,xi,1d0,nf2hq,f2c,EPS)
        fac=2*xi*2*an
        ak=1./rq*(1-rq)*log(1-rq)      
        dc=cf*(-4-0.5/rq-pi**2/3.-(1+3*rq)/2./rq*ak + ak
     + +(xi+xi**2/2.+2*log(1-xi))*log((q2+rm2)/rmu2)
     + +0.5*(rq*xi-rq**2*xi+log(1-rq*xi)*(1-rq*xi))/rq**2/(1-rq*xi)
     - -2*(-(log(1-xi))**2+ddilog((1-rq*xi)/(1-rq))-ddilog(1./(1-rq))
     + +log(1-rq*xi)*log(rq*(1-xi)/(1-rq)))  
     - -2*log(1-rq*xi)/rq
     + +2*(1-rq)/rq*log(1-rq))*f2cqqpm(nb,nt,ni,xi,rmu2)/xi
        f2nucharm=f2nucharm+(f2c+dc)*fac

!  The term of Ref.[Phys.Lett. B699, 345 (2011)] appearing for the 
!  MS-bar definition of the heavy-quark mass 

        if(msbarm) then
          xip=min(xi*(1+delder),1d0)
          delder0=xip/xi-1.
          xim=xi*(1-delder0)
          f2x=f2cqqpm(nb0,nt0,ni0,xi,rmu2)/xi
          f2xp=f2cqqpm(nb0,nt0,ni0,xip,rmu2)/xip
          f2xm=f2cqqpm(nb0,nt0,ni0,xim,rmu2)/xim
          df2x=(f2xp-f2xm)/2./delder0/xi
          f2nucharm=f2nucharm + 4*d1dec*(f2x+xi*df2x)*(1-xb/xi)*fac
        end if 
        f2nucharm=f2nucharm*rfac0
      end if

!  Restore the initial scheme
      kschemepdf=kschemepdfs

      return
      end
!------------------
      real*8 function f2nucharmi(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'

      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0
      data d1dec /1.33333333d0/

      xi=xb0/rq
      z=xi/y

      xin=q2ss/rm2
      eta=xin*(1./(rq*y)-1.)-1.

      rl=log((1-rq*y)/(1-rq)/y)
      rlmf=log((rm2+q2ss)/rmu2)

!  The number of fermions for the strong coupling constant.
      nfe=nfloops(q2ss,kschemepdf)

!  The factors appearing if the renormalization scale is not equal
!  to the factorization scale rmu2. 
      rfac0=1.
      rfac1=1.
      alr=-log(rscale)
      bet0=(11-2./3.*nfe)*alr
      bet1=(102-38/3.*nfe)*alr

      if (kordhq.ge.1) then
        rfac0=1-an*bet0
      end if
      if (kordhq.eq.2) then
        rfac0=1-an*bet0-an**2*(bet1-bet0**2)
        rfac1=1-an*2*bet0
      end if

      glu0=xqg(1,z,rmu2,kschemepdf)

!  The gluon-initiated contribution 
      pqg0=(y**2+(1-y)**2)/2.
      f2g=pqg0*(rl+rlmf)
      f2g=f2g+pqg0*(2*log(1-y)-log(1-rq*y)-log(y))
      f2g=f2g+y*(1-y)*(8-18*(1-rq)+12*(1-rq)**2)
      f2g=f2g+(1-rq)/(1-rq*y)-1
      f2g=f2g+(1-rq)*y*rl*(6*rq+rq*y*(-12*rq))
      f2g=f2g*glu0

!  The quark-initiated contribution 
      f2z=f2cqqpm(nb0,nt0,ni0,z,rmu2)
      f2x=f2cqqpm(nb0,nt0,ni0,xi,rmu2)

      f2q=-(1+y**2)*log(y)/(1-y)
      f2q=cf*f2q*f2z

      f2q0=(1+y**2)/(1-y)*rlmf*(f2z-f2x)
      f2q0=f2q0+(2*log(1-y)-log(1-rq*y))/(1-y)*((1+y**2)*f2z-2*f2x)
      f2q0=f2q0+2.*(1+y)*f2z
      f2q0=f2q0+1./(1-rq*y)*((-1-y)*f2z+2*f2x)
      f2q0=f2q0+0.5*(1-y)/(1-rq*y)**2*(f2z-f2x)
      f2q0=f2q0+2*(rq-1)/(1-y)/(1-rq*y)*(f2z-f2x)
      f2q0=cf*f2q0

!  The total NLO term 
      f2nucharmi=(f2q0+f2q+f2g)/xi*rfac0

!  The asymptotic NNLO terms 
      if (kordhq.ge.2.and.q2ss.ge.50d0) then 
!  Gluon contribution 
        call c2hqg20cc_asy(eta,xin,cc20)
        call c2hqg21cc_asy(eta,xin,cc21)
        call chqg22cc_asy(eta,xin,cc22)
 
        f2charmnnloi=glu0*(cc20 + cc21*rlmf + cc21*rlmf**2)/16.

!  Pure-singlet contribution 
        ps0=xpscqqpm(z,rmu2)

        call c2hqps20cc_asy(eta,xin,cc20)
        call c2hqps21cc_asy(eta,xin,cc21)
        call chqps22cc_asy(eta,xin,cc22)

        f2charmnnloi=f2charmnnloi + ps0*(cc20 
     +   + cc21*rlmf + cc22*rlmf**2)/16.

        if(msbarm) then
          dhg10 = 2-4*y+4*y**2
          f2charmnnloi=f2charmnnloi - d1dec/4.*glu0
        end if

        f2nucharmi=f2nucharmi + f2charmnnloi*8*an/xi*rfac1
      end if

      return
      end


