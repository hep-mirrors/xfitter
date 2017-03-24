c The functions OME...(Z,R) return the values of the heavy-quark OMEs
c coefficients calculated by Buza-Matiounine-Smith-van Neerven 
c in [Eur.Phys.J.C1:301-320,1998].
c The argument Z iz longitudinal variable, R=(\mu/m_H)^2.
c A common in 'CONSTCOM.' must contain variables ZETA2 and ZETA3, values of 
c the Rieman zeta-function at argument equal to 2 and 3, respectively,
c and the QCD structure constants TR=1/2, CG=3, CF=4/3.
c The dilog function DDILOG and polylog function WGPLG 
c from the CERN library are used.

C------------------
      real*8 function ome_q_2(z,r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{Hq}^{PS}
      include 'PRECCOM.'

      real*8 cc(2)

      d0l=-log(r)

      if (omeint) then 
c  grid interpolation 
        call omeintx(2,z,cc)
        ome_q_2=ome_q_2_2(z)*d0l**2 + cc(2)*d0l + cc(1)
      else
c  exact expression 
        ome_q_2=ome_q_2_2(z)*d0l**2 + ome_q_2_1(z)*d0l + ome_q_2_0(z)
      end if

      return
      end
C------------------
      real*8 function ome_q_2_2(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{Hq}^{PS}, the ln(\mu/m)^2 term

      include 'CONSTCOM.'
   
      dlz=log(z)

      a2=-8*(1+z)*dlz-16/z/3d0-4+4*z+16*z*z/3d0

      ome_q_2_2=cf*tr*a2

      return
      end
C------------------
      real*8 function ome_q_2_1(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{Hq}^{PS}, the ln(\mu/m)^1 term

      include 'CONSTCOM.'
   
      dlz=log(z)
      dlz2=dlz*dlz

      a1=8*(1+z)*dlz2-(8+40*z+64*z*z/3d0)*dlz
     -   - 160/z/9d0+16-48*z+448*z*z/9d0

      ome_q_2_1=cf*tr*a1

      return
      end
C------------------
      real*8 function ome_q_2_0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{Hq}^{PS}, the ln(\mu/m)^0 term
   
      include 'CONSTCOM.'
      complex*16 WGPLG

      s121mz=dreal(WGPLG(1,2,1d0-z))
      s111mz=ddilog(1-z)
      dlz=log(z)
      dlz2=dlz*dlz
      dlz3=dlz2*dlz

      a0=(1+z)*(32*s121mz+16*dlz*s111mz-16*zeta2*dlz
     -   - 4.*dlz3/3d0)+(32./z/3d0+8-8*z-32.*z*z/3d0)*s111mz
     +   + (-32./z/3d0-8+8*z+32.*z*z/3d0)*zeta2
     +   + (2+10*z+16*z*z/3d0)*dlz2-(56./3d0+88.*z/3d0+448.*z*z/9d0)*dlz
     -   - 448./z/27d0-4./3d0-124.*z/3d0+1600.*z*z/27d0

      ome_q_2_0=cf*tr*a0

      return
      end
C------------------
      real*8 function ome_g_1(z,r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s) coefficient in A_{Hg}^{PS}  

      include 'CONSTCOM.'

      ome_g_1=4*tr*(z**2+(1-z)**2)*log(r) 

      return
      end
C------------------
csm new beg
csm new function; derivate of ome_g_1(z,r) wrt to mass
      real*8 function dome_g_1(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

csm  the derivative d/dm of the O(\alpha_s) coefficient in A_{Hg}^{PS}  
csm  1/2*m d/dm log(r) = - 1
      include 'CONSTCOM.'

      dome_g_1= - 4*tr*(z**2+(1-z)**2)

      return
      end
csm new end
C------------------
      real*8 function ome_g_2(z,r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{Hg}^{PS}  
      include 'PRECCOM.'

      real*8 cc(2)
      d0l=-log(r)

      if (omeint) then 
c  grid interpolation
        call omeintx(1,z,cc)
        ome_g_2=ome_g_2_2(z)*d0l**2 + cc(2)*d0l + cc(1)
      else
c  exact expression 
        ome_g_2=ome_g_2_2(z)*d0l**2 + ome_g_2_1(z)*d0l + ome_g_2_0(z)
      end if

      return
      end
C------------------
      real*8 function ome_g_2_2(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{Hg}^{S}, ln(\mu/m)^2 term

      include 'CONSTCOM.'

      dlz=log(z)
      dlm=log(1d0-z)

      a2=(8-16*z+16*z*z)*dlm-(4-8*z+16*z*z)*dlz-2.+8*z
      b2=-(8-16*z+16*z*z)*dlm-(8+32*z)*dlz-16./z/3.
     -      -4-32*z+124.*z*z/3.
 
      ome_g_2_2=a2*cf*tr + b2*cg*tr 

      return
      end
C------------------
      real*8 function ome_g_2_1(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{Hg}^{S}, ln(\mu/m)^1 term

      include 'CONSTCOM.'

      s11mz=ddilog(-z)
      dlz=log(z)
      dlz2=dlz*dlz
      dlm=log(1d0-z)
      dlm2=dlm*dlm
      dlp=log(1d0+z)

      a1=(8-16.0*z+16*z*z)*(2*dlz*dlm-dlm2+2*zeta2)
     -    -  (4-8.0*z+16*z*z)*dlz2-32*z*(1-z)*dlm
     -    -  (12-16.0*z+32*z*z)*dlz-56+116*z-80*z*z

      b1=(16+32*z+32*z*z)*(s11mz+dlz*dlp)+(8-16*z+16*z*z)*dlm2
     +    +  (8+16*z)*dlz2+32*z*zeta2+32*z*(1-z)*dlm
     -    -  (8+64*z+352*z*z/3.)*dlz-160/z/9.+16-200*z+1744*z*z/9.

      ome_g_2_1=a1*cf*tr + b1*cg*tr

      return
      end
C------------------
      real*8 function ome_g_2_0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{Hg}^{S}, ln(\mu/m)^0 term

      include 'CONSTCOM.'
      include 'PDFCOM.'
      complex*16 WGPLG

csm decoupling coefficients to one loop (expansion in alphas/(4*pi)
csm see e.g. hep-ph/0004189 eq.(17)
      d1dec = 4.d0*(4.d0/3.d0)

      s121mz=dreal(WGPLG(1,2,1d0-z))
      s12mz=dreal(WGPLG(1,2,-z))
      s211mz=dreal(WGPLG(2,1,1d0-z))
      s21mz=dreal(WGPLG(2,1,-z))
      s111mz=ddilog(1-z)
      s11mz=ddilog(-z)
      dlz=log(z)
      dlz2=dlz*dlz
      dlz3=dlz2*dlz
      dlm=log(1d0-z)
      dlm2=dlm*dlm
      dlm3=dlm2*dlm
      dlp=log(1d0+z)
      dlp2=dlp*dlp

      a01=(1-2*z+2*z*z)*(8*zeta3+4*dlm3/3
     -    -  8*dlm*s111mz+8*zeta2*dlz-4*dlz*dlm2+2*dlz3/3
     -    -  8*dlz*s111mz+8*s211mz-24*s121mz)
      
      a02=-(4+96*z-64*z*z)*s111mz-(4-48*z+40*z*z)*zeta2
     -    -  (8+48*z-24*z*z)*dlz*dlm
     +    +  (4+8*z-12*z*z)*dlm2-(1+12*z-20*z*z)*dlz2
     -    -  (52*z-48*z*z)*dlm-(16+18*z+48*z*z)*dlz+26
     -    -  82*z+80*z*z+z*z*(-16*zeta2*dlz+4*dlz3/3.0
     +    +  16*dlz*s111mz+ 32*s121mz)

      a0=a01+a02

      b01=(1-2*z+2*z*z)*(-4*dlm3/3.+8*dlm*s111mz
     -    -  8*s211mz)+(1+2*z+2*z*z)*(-8*zeta2*dlp
     -    -  16*dlp*s11mz-8*dlz*dlp2+4*dlz2*dlp+8*dlz*s11mz
     -    -  8*s21mz-16*s12mz)+(16+64*z)*(2*s121mz
     +    +  dlz*s111mz)-(4+8*z)*dlz3/3+(8-32*z
     +    +  16*z*z)*zeta3-(16+64*z)*zeta2*dlz
 
      b02=(16*z+16*z*z)*(s11mz+dlz*dlp)+(32/z/3.+12.
     +    +  64*z-272*z*z/3.)*s111mz-(12+48*z
     -    -  260*z*z/3.+32/z/3.)*zeta2-4*z*z*dlz*dlm
     -    -  (2+8*z-10*z*z)*dlm2+(2+8*z+46*z*z/3.)*dlz2
     +    +  (4+16*z-16*z*z)*dlm-(56/3.+172*z/3.+1600*z*z/9.)*dlz
     -    -  448/z/27.-4/3.-628*z/3.+6352*z*z/27.

      b0=b01+b02
 
      ome_g_2_0=a0*cf*tr + b0*cg*tr

csm addition for MSbar running mass
      if (msbarm) then 
        ome_g_2_0=ome_g_2_0 + 2.d0*d1dec * dome_g_1(z)
      end if 

      return
      end
C------------------
      real*8 function ome_qqns_2(z,r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{qq}^{NS}, regular piece  

      include 'CONSTCOM.'

      dlz=log(z)
      dlz2=dlz*dlz

      d0l=-log(r)
      d0l2=d0l*d0l

      a2=-4./3.*(1+z)

      a1=8./3.*(1+z**2)/(1-z)*dlz+8./9.-88./9.*z

      a0=(1+z**2)/(1-z)*(2./3.*dlz2+20./9.*dlz)
     +  +8./3.*(1-z)*dlz+44./27.-268./27.*z

      ome_qqns_2=cf*tr*(a2*d0l2 + a1*d0l + a0)

      return
      end
C------------------
      real*8 function ome_qqns_2_singular(z,r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{qq}^{NS}, singular piece  

      include 'CONSTCOM.'

      d0l=-log(r)
      d0l2=d0l*d0l

      ome_qqns_2_singular=cf*tr*(8./3.*d0l2 + 80./9.*d0l + 224./27.)
     /    /(1-z)

      return
      end
C------------------
      real*8 function ome_qqns_2_local(z,r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{qq}^{NS}, local piece  

      include 'CONSTCOM.'

      dlm=log(1d0-z)
      d0l=-log(r)
      d0l2=d0l*d0l

      ome_qqns_2_local=cf*tr*((8./3.*dlm+2)*d0l2
     +             +(80./9.*dlm+16./3.*zeta2+2./3.)*d0l
     +             +(224./27.*dlm-8./3.*zeta3+40./9.*zeta2+73./18.))

      return
      end
C------------------
      real*8 function ome_gq_2(z,r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{gq}^{S}  

      include 'CONSTCOM.'

      dlz=log(z)
      dlm=log(1d0-z)
      dlm2=dlm**2

      d0l=-log(r)
      d0l2=d0l*d0l

      a2=16./3./z-16./3.+8./3.*z

      a1=160./9./z-160./9.+128/9.*z+(32./3./z-32./3.+16./3.*z)*dlm

      a0=4./3.*(2./z-2+z)*dlm2+8./9.*(10./z-10+8*z)*dlm
     +  +(448./z-448+344*z)/27.

      ome_gq_2=cf*tr*(a2*d0l2 + a1*d0l + a0)

      return
      end
C------------------
      real*8 function ome_gg_1_local(z,r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s) coefficient in A_{gg}^{S}, local piece 

      include 'CONSTCOM.'

      d0l=-log(r)

      ome_gg_1_local=4./3.*tr*d0l

      return 
      end
C------------------
csm new beg
csm new function; derivate of ome_gg_1_local wrt to mass
      real*8 function dome_gg_1_local(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

csm  the derivative d/dm of the O(\alpha_s) coefficient in A_{gg}^{S}, local piece
csm  1/2*m d/dm log(r) = - 1
      include 'CONSTCOM.'

      dome_gg_1_local=4./3.*tr

      return 
      end
csm new end
C------------------
      real*8 function ome_gg_2(z,r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{gg}^{S}, regular piece

      include 'CONSTCOM.'

      dlz=log(z)
      dlz2=dlz*dlz
      dlz3=dlz2*dlz
      dlm=log(1d0-z)
      d0l=-log(r)
      d0l2=d0l*d0l

      a2=8*(1+z)*dlz+16./3./z+4-4*z-16./3.*z**2
      b2=8./3.*(1./z-2+z-z**2)

      a1=8*(1+z)*dlz2+(24+40*z)*dlz-16./3./z+64-32*z-80./3.*z**2
      b1=16./3.*(1+z)*dlz+184./9./z-232./9.+152./9.*z-184./9.*z**2

      a0=4./3.*(1+z)*dlz3+(6+10*z)*dlz2+(32+48*z)*dlz
     -  -8./z+80-48*z-24*z**2
      b0=4./3.*(1+z)*dlz2+(52+88*z)*dlz/9.-4./3.*z*dlm
     +   +(556./z-628+548*z-700*z**2)/27.

      ome_gg_2=(a2*cf*tr + b2*cg*tr)*d0l2
     +  + (a1*cf*tr + b1*cg*tr)*d0l
     +  +  a0*cf*tr + b0*cg*tr

      return
      end
C------------------
      real*8 function ome_gg_2_singular(z,r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{gg}^{S}, singular piece

      include 'CONSTCOM.'

      d0l=-log(r)
      d0l2=d0l*d0l

      ome_gg_2_singular=cg*tr*(8./3.*d0l2 + 80./9.*d0l + 224./27.)
     /   /(1-z)

      return 
      end
C------------------
      real*8 function ome_gg_2_local(z,r)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c  The O(\alpha_s^2) coefficient in A_{gg}^{S}, local piece  

      include 'CONSTCOM.'
      include 'PDFCOM.'

csm decoupling coefficients to one loop (expansion in alphas/(4*pi)
csm see e.g. hep-ph/0004189 eq.(17)
      d1dec = 4.d0*(4.d0/3.d0)

      dlm=log(1d0-z)
      d0l=-log(r)
      d0l2=d0l*d0l

      ome_gg_2_local=(cg*tr*8./3.*dlm)*d0l2
     +                  +(cg*tr*(16./3.+80./9.*dlm) + 4*cf*tr)*d0l
     +                  +(cg*tr*(10./9.+224./27.*dlm) - cf*tr*15)

csm addition for MSbar running mass
      if (msbarm) then
        ome_gg_2_local= ome_gg_2_local + 2d0*d1dec*dome_gg_1_local(z)
      end if

      return 
      end
C-----------------
      SUBROUTINE omegridini
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'APSCOM6.'

c  Initialization of the grid -- it is called once at the beginning

      xo1=0.2d0
      xbomax=0.999d0

      XLOG0=LOG(Xbmin)
      xolog1=log(xo1)
      xolog2=log(1-xo1)

      DELop=(xolog2-log(1-xbomax))/(nxpgrid-1)
      DELom=-(XLOG0-xolog1)/nxmgrid

      DO I=0,nxpgrid-1
        xogrid(I)=1-exp(xolog2-delop*i)
      end do
      xogrid(nxpgrid)=1.

      DO I=-nxmgrid,0
        xogrid(I)=EXP(xolog1+delom*I)
      end do

      DO IX=-nxmgrid,nxpgrid-1
        z=xogrid(ix)
        omegrid(1,1,ix)=ome_g_2_0(z)*z*(1-z)
        omegrid(1,2,ix)=ome_g_2_1(z)*z*(1-z)
        omegrid(2,1,ix)=ome_q_2_0(z)*z*(1-z)
        omegrid(2,2,ix)=ome_q_2_1(z)*z*(1-z)
      end do

      RETURN
      END
C------------------
      subroutine omeintx(kint,xx,cc)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'APSCOM6.'

      real*8 cc(2)

c  returns the O(\alpha_s^2) OMEs coefficients 
c  for A_Hg(XX) at KINT=1 and A_Hq(XX) at KINT=2
c  in CC(1) (constant term) and CC(2) (linear-log term)

      xb=min(xx,xbomax)
      xl1=log(xb)
      xl2=log(1-xb)

      if (xb.ge.xo1) then
        IX=int((xolog2-xl2)/delop)
        DXX=(xolog2-xl2)/delop-ix
        if (ix.eq.0) then 
          do iq=1,2
            cc(iq)=omegrid(kint,iq,0)*(1-dxx)
     +            +omegrid(kint,iq,1)*dxx
          end do
        else 
          do iq=1,2
            cc(iq)=omegrid(kint,iq,ix-1)*(dxx-1)*dxx/2.
     +            +omegrid(kint,iq,ix)*(1-dxx**2)
     +            +omegrid(kint,iq,ix+1)*(1+dxx)*dxx/2.
          end do
        end if          
      else
        IX=int((XL1-XoLOG1)/DELom)-1
        DXX=(xl1-xolog1)/delom-ix
        do iq=1,2
          cc(iq)=omegrid(kint,iq,ix-1)*(dxx-1)*dxx/2.
     +          +omegrid(kint,iq,ix)*(1-dxx**2)
     +          +omegrid(kint,iq,ix+1)*(1+dxx)*dxx/2.
        end do
      end if

      do iq=1,2
        cc(iq)=cc(iq)/xb/(1-xb)
      end do

      RETURN
      END


