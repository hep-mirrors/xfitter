C------------------
      real*8 function ch2g00_asymp(z,xi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      ch2g00_asymp=ome_g_1(z,xi)+c2g_1_0(z)

      return
      end
C------------------
      real*8 function ch2g11_asymp(z,xi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'

      bb1=(16-32*z+32*z**2)*log(1-z) + (16+64*z)*log(z) 
     +  + 32./3./z+8+64*z-248./3.*z**2
      bb0=(16+64*z)*ddilog(1-z) - (16-32*z+32*z**2)*zeta2
     +  + (96*z-32*z**2)*log(z)*log(1-z) + (16-32*z+32*z**2)*log(1-z)**2
     -  - (8+32*z)*log(z)**2 + (32./3./z-8+192*z-632./3.*z**2)*log(1-z)
     -  - (8+256*z-248./3.*z**2)*log(z) 
     +  + 16./3./z-172./3.-968./3.*z+1124./3*z**2

      ch2g11_asymp=-cg*tr*(bb0+bb1*log(xi))

      return
      end
C------------------
      real*8 function ch2g10_asymp(z,xi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      complex*16 WGPLG

      s121mz=realpart(WGPLG(1,2,1d0-z))
      s12mz=realpart(WGPLG(1,2,-z))
      s211mz=realpart(WGPLG(2,1,1d0-z))
      s21mz=realpart(WGPLG(2,1,-z))
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

      a2=(8-16*z+16*z*z)*dlm-(4-8*z+16*z*z)*dlz-2.+8*z
      b2=-(8-16*z+16*z*z)*dlm-(8+32*z)*dlz-16/z/3.
     -      -4-32*z+124*z*z/3.

      a1=(8-16*z)*s111mz - (32-64*z+64*z**2)*zeta2
     -    - (24-48*z+64*z**2)*log(z)*log(1-z)
     +    + (16-32*z+32*z**2)*log(1-z)**2 + (8-16*z+32*z**2)*log(z)**2
     -    - (28-96*z+80*z**2)*log(1-z) + (8-48*z+80*z**2)*log(z)
     +    + 36-68*z+16*z**2

      b1=-(16+32*z+32*z**2)*(s11mz+log(z)*log(1+z)) + (16+64*z)*s111mz
     -    - (16+32*z**2)*zeta2 + (96*z-32*z**2)*log(z)*log(1-z)
     +    + (8-16*z+16*z**2)*log(1-z)**2 - (16+48*z)*log(z)**2
     +    + (32./3./z-8+160*z-536./3.*z**2)*log(1-z)
     -    - (192*z-200*z**2)*log(z)+208/9./z-220/3.-368./3.*z
     +    + 1628/9.*z**2

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

      ch2g10_asymp=(a2*cf*tr - b2*cg*tr)*log(xi)**2
     +  + (a1*cf*tr + b1*cg*tr)*log(xi)
     +  +  a0*cf*tr + b0*cg*tr
     +  +  c2g_2_0(z)

      return
      end
C------------------
      real*8 function ch2qps10_asymp(z,xi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      complex*16 WGPLG

      s121mz=realpart(WGPLG(1,2,1d0-z))
      s111mz=ddilog(1-z)
      dlz=log(z)
      dlz2=dlz*dlz
      dlz3=dlz2*dlz
      dlm=log(1d0-z)

      a2=-8*(1+z)*dlz-16./z/3d0-4+4*z+16.*z*z/3d0

      a1=16*(1+z)*(s111mz+dlz*dlm-dlz2) 
     +  + (32./z/3d0+8-8*z-32./3d0*z**2)*dlm 
     +  + 32*z**2*dlz+208./9d0/z-208./3d0+160./3d0*z-64./9d0*z**2

      a0=(1+z)*(32*s121mz+16*dlz*s111mz-16*zeta2*dlz
     -   - 4*dlz3/3d0)+(32/z/3d0+8-8*z-32*z*z/3d0)*s111mz
     +   + (-32/z/3d0-8+8*z+32*z*z/3d0)*zeta2
     +   + (2+10*z+16*z*z/3d0)*dlz2 - (56/3d0+88*z/3d0+448*z*z/9d0)*dlz
     -   - 448/z/27d0-4/3d0-124*z/3d0+1600*z*z/27d0

      ch2qps10_asymp=cf*tr*(-a2*log(xi)**2 + a1*log(xi) + a0)
     +  +  c2ps_2_0(z)

      return
      end
C------------------
      real*8 function ch2qps11_asymp(z,xi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      complex*16 WGPLG

      s111mz=ddilog(1-z)
      dlz=log(z)
      dlz2=dlz*dlz
      dlz3=dlz2*dlz
      dlm=log(1d0-z)

      b1=16*(1+z)*dlz+32./3d0/z+8-8*z-32./3d0*z**2      

      b0=8*(1+z)*(2*s111mz+2*dlz*dlm-dlz2)
     +  + (32./3d0/z+8-8*z-32./3d0*z**2)*dlm-(8+40*z-32./3d0*z**2)*dlz
     +  + 16./3d0/z-160./3d0+16./3d0*z+128./3d0*z**2 

      ch2qps11_asymp=-cf*tr*(b1*log(xi) +b0) 

      return
      end
C------------------
      real*8 function ch2qns10_asymp(z,xi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      complex*16 WGPLG

      s111mz=ddilog(1-z)
      dlz=log(z)
      dlz2=dlz*dlz
      dlz3=dlz2*dlz
      dlm=log(1d0-z)

      a2=4./3.*(1+z**2)/(1-z)

      a1=(1+z**2)/(1-z)*(8./3.*dlm-16./3.*dlz-58./9.)
     +  +2./3.+26./3.*z

      a0=(1+z**2)/(1-z)*(-8./3.*s111mz-8./3.*zeta2
     -  -16./3.*dlz*dlm+4./3.*dlm**2+4*dlz2-58./9.*dlm
     +  +134./9.*dlz+359./27.)
     +  +(2./3.+26./3.*z)*dlm-(2+46./3.*z)*dlz+29./9.-295./9.*z

      ch2qns10_asymp=cf*tr*(a2*log(xi)**2 + a1*log(xi) + a0)

      return
      end
C------------------
      real*8 function chlqns10_asymp(z,xi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      complex*16 WGPLG

      dlz=log(z)
      dlm=log(1d0-z)

      chlqns10_asymp=cf*tr*(16./3.*z*(log(xi)+dlm-2*dlz)
     +  +16./3.-400./18*z)

      return
      end
C------------------
      real*8 function ch2qns10_exact(z,xi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'

      sq1=sqrt(1-4*z/(1-z)/xi)
      sq2=sqrt(1-4*z/xi)
      rl1=log((1+sq1)/(1-sq1))
      rl2=log((1+sq2)/(1-sq2))
      rl3=log((sq1+sq2)/(sq2-sq1))
      dil1=ddilog((1-z)*(1+sq1)/(1+sq2))
      dil2=ddilog((1-sq2)/(1+sq1))
      dil3=ddilog((1-sq1)/(1+sq2))
      dil4=ddilog((1+sq1)/(1+sq2))

      ch2qns10_exact=cf*tr*((4./3.*(1+z**2)/(1-z)
     -  -16./(1-z)*z**2/xi**2*(1-9*z+9*z**2))
     *  *(log((1-z)/z**2)*rl1+rl1*rl2+2*(-dil1+dil2+dil3-dil4))    
     +  +(-8./3.+4./(1-z)
     +  +(z/(1-z)/xi)**2*(128-432*z+288*z**2-8./(1-z)))*rl1
     +  +rl3/sq2*(88./9.+136./9.*z-152/9./(1-z)
     +  +(z/(1-z)/xi)*(464./9.-512./3.*z+2048./9.*z**2)
     +  +(z/(1-z)/xi)**2
     *  *(-832./9.+6208./9.*z-11392./9.*z**2+6016./9.*z**3))
     +  +sq1*(-272./27.-1244./27.*z+718./27./(1-z)
     +  +(z/(1-z)/xi)*(-3424./27.+15608./27.*z
     -  -4304./9.*z**2+20./27./(1-z))))

      return 
      end
C------------------
      real*8 function chlqns10_exact(z,xi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'

      sq1=sqrt(1-4*z/(1-z)/xi)
      sq2=sqrt(1-4*z/xi)
      rl1=log((1+sq1)/(1-sq1))
      rl2=log((1+sq2)/(1-sq2))
      rl3=log((sq1+sq2)/(sq2-sq1))
      dil1=ddilog((1-z)*(1+sq1)/(1+sq2))
      dil2=ddilog((1-sq2)/(1+sq1))
      dil3=ddilog((1-sq1)/(1+sq2))
      dil4=ddilog((1+sq1)/(1+sq2))

      chlqns10_exact=cf*tr*(96*z**3/xi**2
     *  *(log((1-z)/z**2)*rl1+rl1*rl2+2*(-dil1+dil2+dil3-dil4))
     +  +(z/(1-z)/xi)**2*(64-288*z+192*z**2)*dl1
     +  +z*(16./3.-416./3.*z/xi+1408./3.*z**2/xi**2)*rl3/sq2
     +  +(16./3.-400./18.*z
     +  +(z/(1-z)/xi)*(-160./3.+3824/9.*z-992./3.*z**2))*sq1)

      return 
      end
