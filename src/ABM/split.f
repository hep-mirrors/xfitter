C-----------------
      subroutine pgridini
      implicit double precision (a-h,o-z)

      include 'APSCOM6.'
      include 'CONSTCOM.'
      real*8 p(nflim)

c set up the steps for the x-grid

c log(x) transform for the low-x region x 
      xlog0=log(xbmin)
      delx2=-(xlog0-xlog1)/nxmgrid
c log(1-x) transform for the high-x region 
      delxp=(xlog2-log(1-xbmax))/(nxpgrid-1)

c fill the x-grid

c large-x region 
      do i=0,nxpgrid-1
        xpgrid(i)=1-exp(xlog2-delxp*i)
      end do
      xpgrid(nxpgrid)=1.

c small-x region 
      do i=-nxmgrid,0
        xpgrid(i)=exp(xlog1+delx2*i)
      end do

c  Select the upper margin for the initial fermion number, nfc, 
c  in order to have the grid ready for the 
c  cases when the initial number of fermion and the number of 
c  fermions in the loops, nfe, are not equal.

      nfc=(nflim-1)/2
      do i=3,6

        nfe=i
c  The regular pieces of the splitting functions do not depend
c  on the initial fermion number, nfc, up to the NNLO therefore this 
c  splitting function interpolation can be used even if the 
c  initial fermion number is not equal to the number of fermions 
c  in the loops, nfe. 

        do iq=1,2*nfc+1
          DO IX=-nxmgrid,nxpgrid-1
            z=xpgrid(ix)
            call QGSPLIT(IQ,p,Z)
            do jq=1,2*nfc+1
              pgrid(0,0,nfe,iq,jq,ix)=p(jq)
            end do
            call QGSPLITnlo(IQ,p,Z)
            do jq=1,2*nfc+1
              pgrid(0,1,nfe,iq,jq,ix)=p(jq)
            end do
            call QGSPLITnnlo(IQ,p,Z,kernelmod)
            do jq=1,2*nfc+1
              pgrid(0,2,nfe,iq,jq,ix)=p(jq)
            end do
          end do 
        end do
      end do

      return 
      end
C------------------
      subroutine pintx(k,iq,xx,p)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'APSCOM6.'
      include 'CONSTCOM.'
      INCLUDE 'PDFCOM.'

      real*8 p(0:2,nflim)

      xb=min(xx,xbmax)
      xl1=log(xb)
      xl2=log(1-xb)

      if (xb.ge.x1) then
        IX=int((xlog2-xl2)/delxp)
        DXX=(xlog2-xl2)/delxp-ix
        if (ix.eq.0) then 
          do jq=1,2*nfc+1
            do ko=0,2
              p(ko,jq)=pgrid(k,ko,nfe,IQ,jq,IX)*(1-dxx)
     +               +pgrid(k,ko,nfe,IQ,jq,IX+1)*dxx
            end do
          end do
        else 
          do jq=1,2*nfc+1
            do ko=0,2
              p(ko,jq)=PGRID(k,ko,nfe,IQ,jq,IX-1)*(dxx-1)*dxx/2.
     +          +PGRID(k,ko,nfe,IQ,jq,IX)*(1-dxx**2)
     +          +PGRID(k,ko,nfe,IQ,jq,IX+1)*(1+dxx)*dxx/2.
            end do
          end do
        end if          
      else
        IX=int((XL1-XLOG1)/DELX2)-1
        DXX=(xl1-xlog1)/delx2-ix
          do jq=1,2*nfc+1
            do ko=0,2
              p(ko,jq)=PGRID(k,ko,nfe,IQ,jq,IX-1)*(dxx-1)*dxx/2.
     +          +PGRID(k,ko,nfe,IQ,jq,IX)*(1-dxx**2)
     +          +PGRID(k,ko,nfe,IQ,jq,IX+1)*(1+dxx)*dxx/2.
            end do
          end do
      end if

      RETURN
      END
C------------------
      real*8 function pgg_0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The LO gluon-gluon splitting function of [Nucl. Phys. B126, 298 (1977)]
c  (regular piece) 
      include 'CONSTCOM.'

      pgg_0=2*cg*(1./Z +Z*(1-Z) - 2)

      return 
      end
C------------------
      real*8 function pqq_0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The LO quark-quark splitting function of [Nucl. Phys. B126, 298 (1977)]
c  (regular piece) 
      include 'CONSTCOM.'

      pqq_0=-cf*(Z+1)

      return 
      end
C------------------
      real*8 function pqg_0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The LO quark-gluon splitting function of [Nucl. Phys. B126, 298 (1977)]
c  (regular piece) 
      include 'CONSTCOM.'

      pqg_0=tr*(Z**2+(1-Z)**2)

      return 
      end
C------------------
      real*8 function pgq_0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The LO gluon-quark splitting function of [Nucl. Phys. B126, 298 (1977)]
c  (regular piece) 
      include 'CONSTCOM.'

      pgq_0=cf*(1+(1-Z)**2)/Z

      return 
      end
C------------------
      real*8 function pgg_0_singular(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The LO gluon-gluon splitting function of [Nucl. Phys. B126, 298 (1977)]
c  (singular piece) 
      include 'CONSTCOM.'

      pgg_0_singular=-2.*cg/(1.-Z)

      return 
      end
C------------------
      real*8 function pqq_0_singular(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The LO quark-quark splitting function of [Nucl. Phys. B126, 298 (1977)]
c  (singular piece) 

      include 'CONSTCOM.'

      pqq_0_singular=-2.*cf/(1.-Z)

      return 
      end
C------------------
      real*8 function pgg_0_local_nf0(x)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The LO gluon-gluon splitting function of [Nucl. Phys. B126, 298 (1977)]
c  (NF-independent, local piece) 

      include 'CONSTCOM.'

      pgg_0_local_nf0=2*cg*(11./12.+LOG(1-X))

      return 
      end
C------------------
      real*8 function pqq_0_local(x)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The LO quark-quark splitting function of [Nucl. Phys. B126, 298 (1977)]
c  (local piece) 

      include 'CONSTCOM.'

      pqq_0_local=cf*(2*LOG(1-X)+1.5)

      return 
      end
C------------------
      real*8 function pgg_1_nf0_singular(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO gluon-gluon splitting function of [Phys. Lett. B97, 437 (1980)]
c  (NF-independent, singular piece) 
      include 'CONSTCOM.'

      pgg_1_nf0_singular=-cg**2*(67./9.-pi**2/3.)/(1.-z)

      return 
      end
C------------------
      real*8 function pgg_1_nf1_singular(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO gluon-gluon splitting function of [Phys. Lett. B97, 437 (1980)]
c  (NF^1, singular piece) 
      include 'CONSTCOM.'

      pgg_1_nf1_singular=20./9.*cg*tr/(1-z)

      return 
      end
C------------------
      real*8 function pqq_1_nf0_singular(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO quark-quark splitting function of [Nucl. Phys. B175, 27 (1980)]
c  (NF-independent, singular piece) 
      include 'CONSTCOM.'

      pqq_1_nf0_singular=-2*(cf*cg*(67./18.-pi**2/6.))/(1-z)

      return 
      end
C------------------
      real*8 function pqq_1_nf1_singular(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO quark-quark splitting function of [Nucl. Phys. B175, 27 (1980)]
c  (NF^1, singular piece) 
      include 'CONSTCOM.'

      pqq_1_nf1_singular=cf*tr*20./9./(1.-z)

      return 
      end
C------------------
      real*8 function pgg_1_nf0_local(x)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO gluon-gluon splitting function of [Phys. Lett. B97, 437 (1980)]
c  (NF-independent, local piece) 
      include 'CONSTCOM.'
c The local term in (-x*pgg_1_nf0)_+
      data gg3 /-2.11826d0/

      pgg_1_nf0_local=-cg**2*gg3+(log(1-x)+1)*(67./9.-pi**2/3.)*cg**2

      return 
      end
C------------------
      real*8 function pgg_1_nf1_local(x)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO gluon-gluon splitting function of [Phys. Lett. B97, 437 (1980)]
c  (NF^1, local piece) 
      include 'CONSTCOM.'
c The local terms in (-x*pgg_1_nf1)_+
      data gg1,gg2 /-0.37037d0,-1.53704d0/

      pgg_1_nf1_local=-cf*tr*gg1-tr*cg*gg2-(log(1-x)+1)*20./9.*cg*tr

      return 
      end
C------------------
      real*8 function pqq_1_nf0_local(x)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO quark-quark splitting function of [Nucl. Phys. B175, 27 (1980)]
c  (NF-independent, local piece) 
      include 'CONSTCOM.'
c The local terms in (-pqqns_1_nf0)_+
      data ff1,ff2 /-2.65254d0,-0.01766d0/

      pqq_1_nf0_local=-(cf**2*ff1+cf*cg*ff2)
     +   +(cf*cg*(67./18.-pi**2/6.))*(2*LOG(1-X)+1.5)

      return 
      end
C------------------
      real*8 function pqq_1_nf1_local(x)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO quark-quark splitting function of [Nucl. Phys. B175, 27 (1980)]
c  (NF^1, local piece) 
      include 'CONSTCOM.'
c The local term in (-pqqns_1_nf1)_+, equal
c to the sum of terms in (-x*pqqns_1_nf1)_+, (-pgq_1_nf1)_+, (-2x*pqqps_1)_+
c                          (1.322877+20./9.,    -1.92593,     2*0.37037)
      data ff3 /2.35991d0/
c The local term in (-2x*pqqps_1)_+
      data ffps /0.37037d0/ 

      pqq_1_nf1_local=cf*tr*(-ff3 - 20./9.*LOG(1-X)
c the pure-singlet local term subtracted
     +  + 2*ffps)
      return 
      end
C------------------
      real*8 function pgg_1_nf0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO gluon-gluon splitting function of [Phys. Lett. B97, 437 (1980)]
c  (NF-independent, regular piece) 
      include 'CONSTCOM.'

      pgg3=27./2*(1.-z)+67./9.*(z**2-1./z)
     ++(-25./3.+11./3.*z-44./3.*z**2)*dlog(z)+4.*(1.+z)*dlog(z)**2
     ++(67./9.-4.*dlog(z)*dlog(1.-z)+dlog(z)**2-pi**2/3.)
     **(1./(1.-z)+1./z-2.+z-z**2)+2*(1./(1.+z)-1./z-2.-z-z**2)*s2nlo(z)

      pgg_1_nf0=cg**2*(pgg3-(67./9.-pi**2/3.)/(1-z))

      return 
      end
C------------------
      real*8 function pgg_1_nf1(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO gluon-gluon splitting function of [Phys. Lett. B97, 437 (1980)]
c  (NF^1, regular piece) 
      include 'CONSTCOM.'

      pgg1=-16.+8.*z+20./3.*z**2+4./3./z+(-6.-10.*z)*dlog(z)
     ++(-2.-2.*z)*dlog(z)**2
      pgg2=2.-2.*z+26./9*z**2-26./9./z-4./3.*(1.+z)*dlog(z)
     --20./9*(1./(1.-z)+1./z-2.+z-z**2) + 20./9./(1-z)

      pgg_1_nf1=cf*tr*pgg1+cg*tr*pgg2

      return 
      end
C------------------
      real*8 function pgq_1_nf0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO gluon-quark splitting function of [Phys. Lett. B97, 437 (1980)]
c  (NF-independent, regular piece) 
      include 'CONSTCOM.'

      pfg1=-5./2.-7./2*z+(2.+7./2*z)*dlog(z)+(-1.+z/2.)*dlog(z)**2
     --2.*z*dlog(1.-z)+(-3.*dlog(1.-z)-dlog(1.-z)**2)*(1.+(1.-z)**2)/z
      pfg2=28./9.+65./18.*z+44./9.*z**2+(-12.-5.*z-8./3.*z**2)*dlog(z)
     ++(4.+z)*dlog(z)**2+2.*z*dlog(1.-z)+(-2.*dlog(z)*dlog(1.-z)
     ++0.5*dlog(z)**2+11./3.*dlog(1.-z)+dlog(1.-z)**2-pi**2/6.+0.5)
     **(1.+(1.-z)**2)/z-(1.+(1.+z)**2)/z*s2nlo(z)

      pgq_1_nf0=cf**2*pfg1+cf*cg*pfg2

      return 
      end
C------------------
      real*8 function pgq_1_nf1(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO gluon-quark splitting function of [Phys. Lett. B97, 437 (1980)]
c  (NF^1, regular piece) 
      include 'CONSTCOM.'

      pfg3=-4./3.*z-(20./9.+4./3.*dlog(1.-z))*(1.+(1.-z)**2)/z

      pgq_1_nf1=cf*tr*pfg3

      return 
      end
C------------------
      real*8 function pqg_1(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO quark-gluon splitting function of [Phys. Lett. B97, 437 (1980)]
c  (regular piece) 
      include 'CONSTCOM.'

      pgf1=4.-9.*z+(-1.+4.*z)*dlog(z)+(-1.+2.*z)*dlog(z)**2
     ++4.*dlog(1.-z)+(-4.*dlog(z)*dlog(1.-z)+4.*dlog(z)+2.*dlog(z)**2
     --4.*dlog(1.-z)+2.*dlog(1.-z)**2-2./3.*pi**2+10.)*(z**2+(1.-z)**2)
      pgf2=182./9.+14./9.*z+40./9./z+(136./3*z-38./3.)*dlog(z)
     --4.*dlog(1.-z)-(2.+8.*z)*dlog(z)**2+(-dlog(z)**2+44./3.*dlog(z)
     --2.*dlog(1.-z)**2+4.*dlog(1.-z)+pi**2/3.-218./9.)*(z**2+(1.-z)**2)
     ++2.*(z**2+(1.+z)**2)*s2nlo(z)

      pqg_1=cf*tr/2.*pgf1+cg*tr/2.*pgf2

      return 
      end
C------------------
      real*8 function pqqns_1_nf0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO quark-quark splitting function of [Nucl. Phys. B175, 27 (1980)]
c  (NF-independent, non-singlet, regular piece) 
      include 'CONSTCOM.'

      pf=-2.*(1.+z**2)/(1.-z)*log(z)*log(1.-z)
     --(3./(1.-z)+2.*z)*log(z)-0.5*(1.+z)*log(z)**2-5.*(1.-z)
      pg=(1.+z**2)/(1.-z)*(log(z)**2+11./3.*log(z)+67./9.-pi**2/3.)
     ++2.*(1.+z)*log(z)+40./3.*(1.-z) - 2*(67./9.-pi**2/3.)/(1-z)
      pqqns_1_nf0=cf**2*pf+cf*cg/2.*pg

      return 
      end
C------------------
      real*8 function pqqns_1_nf1(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO quark-quark splitting function of [Nucl. Phys. B175, 27 (1980)]
c  (NF^1, non-singlet, regular piece) 
      include 'CONSTCOM.'

      pqqns_1_nf1=cf*tr*(2./3.*((1+z**2)/(1-z)*(-log(z)-5./3.)-2*(1-z))
     +   +  20./9./(1.-z))

      return 
      end
C------------------
      real*8 function pqaqns_1(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO quark-antiquark splitting function of [Nucl. Phys. B175, 27 (1980)]
c  (non-singlet, regular piece) 
      include 'CONSTCOM.'

      pa=2*(1+z**2)/(1+z)*s2nlo(z)+2*(1+z)*log(z)+4*(1-z)
      pqaqns_1=(cf**2-cf*cg/2.)*pa

      return 
      end
C------------------
      real*8 function pqqps_1(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NLO quark-quark splitting function of [Phys. Lett. B97, 437 (1980)]
c  (pure-singlet, regular piece) 
      include 'CONSTCOM.'

      pqqps_1=cf*tr*(-4+12*z+(10*z+16./3.*z**2+2)*log(z)-112./9.*z**2
     +  +40./9./z-2*(1+z)*log(z)**2)/2.

      return 
      end
C------------------
      real*8 function pgg_2_nf0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO gluon-gluon splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF-independent, regular piece) 

      DL  = LOG (y)
      DL1 = LOG (1-y)

      pgg_2_nf0= 3589*dl1 - 20852 + 3968*y - 3363*y**2 + 4848*y**3
     +  + dl*dl1*(7305 + 8757*dl) + 274.4*dl - 7471*dl**2 + 72*dl**3
     -  - 144*dl**4 + 14214./y + 2675.8/y*dl

      return 
      end
C------------------
      real*8 function pgg_2_nf1(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO gluon-gluon splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF^1, regular piece) 

      DL  = LOG (y)
      DL1 = LOG (1-y)

      pgg_2_nf1=-320*dl1 - 350.2 + 755.7*y - 713.8*y**2 + 559.3*y**3
     +   + dl*dl1*(26.15 - 808.7*dl) + 1541*dl + 491.3*dl**2 
     +   + 832./9.*dl**3 + 512./27.*dl**4 + 182.96/y + 157.27/y*dl

      return 
      end
C------------------
      real*8 function pgg_2_nf2(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO gluon-gluon splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF^2, regular piece) 

      DL  = LOG (y)
      DL1 = LOG (1-y)

      pgg_2_nf2=-13.878 + 153.4*y - 187.7*y**2 + 52.75*y**3
     -   - dl*dl1*(115.6 - 85.25*y + 63.23*dl) - 3.422*dl + 9.68*dl**2
     -   - 32./27.*dl**3-680./243./y

      return 
      end
C------------------
      real*8 function pgq_2_nf0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO gluon-quark splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF-independent, regular piece) 

      DL  = LOG (Y)
      DL1 = LOG (1-Y)

      pgq_2_nf0=400./81.*dl1**4 + 2200./27.*dl1**3 + 606.3*dl1**2
     +  + 2193*dl1 - 4307 + 489.3*y + 1452*y**2 + 146*y**3 
     -  - 447.3*dl**2*dl1 - 972.9*y*dl**2 + 4033*dl - 1794*dl**2
     +  + 1568./9.*dl**3 - 4288./81.*dl**4 + 6163.1/y + 1189.3/y*dl

      return 
      end
C------------------
      real*8 function pgq_2_nf1(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO gluon-quark splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF^1, regular piece) 

      DL  = LOG (Y)
      DL1 = LOG (1-Y)

      pgq_2_nf1=-400./81*dl1**3 - 68.069*dl1**2 - 296.7*dl1 - 183.8
     +   + 33.5*y - 277.9*y**2 + 108.6*y*dl**2 - 49.68*dl*dl1 
     +   + 174.8*dl + 20.39*dl**2 + 704./81.*dl**3 + 128./27.*dl**4
     -   - 46.41/y + 71.082/y*dl

      return 
      end
C------------------
      real*8 function pgq_2_nf2(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO gluon-quark splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF^2, regular piece) 

      DL  = LOG (Y)
      DL1 = LOG (1-Y)

      pgq_2_nf2=96./27.*dl1**2 * (1./y-1+y/2.) 
     +   + 320./27.*dl1 * (1./y-1+4./5.*y) - 64./27. * (1./y-1-2*y)

      return 
      end
C------------------
      real*8 function pqg_2_nf0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-gluon splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF-independent, regular piece) 

      DL  = LOG (Y)
      DL1 = LOG (1-Y)

      pqg_2_nf0=(100./27.*dl1**4 - 70./9.*dl1**3 - 120.5*dl1**2
     +   + 104.42*dl1 + 2522 - 3316*y + 2126*y**2 
     +   + dl*dl1*(1823 - 25.22*dl) - 252.5*y*dl**3 + 424.9*dl 
     +   + 881.5*dl**2 - 44./3.*dl**3 + 536./27.*dl**4 - 1268.3/y
     -   - 896./3./y*dl)/2.

      return 
      end
C------------------
      real*8 function pqg_2_nf1(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-gluon splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF^1, regular piece) 

      DL  = LOG (Y)
      DL1 = LOG (1-Y)

      pqg_2_nf1= (20./27.*dl1**3 + 200./27*dl1**2 
     -   - 5.496*dl1 - 252 + 158*y + 145.4*y**2 - 139.28*y**3 
     -   - dl*dl1*(53.09 + 80.616*dl) - 98.07*y*dl**2 + 11.7*y*dl**3
     -   - 254*dl - 90.8*dl**2 - 376./27.*dl**3 - 16./9*dl**4 
     +   + 1112./243./y)/2.

      return 
      end
C------------------
      real*8 function pnsp_2_nf0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl Phys B688, 101 (2004)]
c  (NF-independent, non-singlet, C-even, regular piece) 

      DL  = LOG (Y)
      DL1 = LOG (1-Y)

      pnsp_2_nf0= 714.1*DL1 + 1641.1 - 3135*y + 243.6*Y**2 
     -  - 522.1*y**3 + dl*dl1*(563.9 + 256.8*dl) + 1258*dl
     +  + 294.9*dl**2 + 800./27.*dl**3 + 128./81.*dl**4

      return 
      end
C------------------
      real*8 function pnsp_2_nf1(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl Phys B688, 101 (2004)]
c  (NF^1, non-singlet, C-even, regular piece) 

      DL  = LOG (Y)
      DL1 = LOG (1-Y)

      pnsp_2_nf1=-5120./81.*dl1-197.+381.1*y+72.94*y**2+44.79*y**3
     -            -1.497*y*dl**3-56.66*dl*dl1
     -            -152.6*dl-2608./81.*dl**2-64./27.*dl**3 

      return 
      end
C------------------
      real*8 function pnsp_2_nf2(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl Phys B688, 101 (2004)]
c  (NF^2, non-singlet, C-even, regular piece) 

      DL  = LOG (Y)

      pnsp_2_nf2=-(- 32.* Y*DL/(1.-Y) * (3.* DL + 10.) - 64.
     -           - (48.* DL**2 + 352.* DL + 384.) * (1.-Y) )/81.D0

      return 
      end
C------------------
      real*8 function pnsm_2_nf0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl Phys B688, 101 (2004)]
c  (NF-independent, non-singlet, C-odd, regular piece) 

      DL  = LOG (Y)
      DL1 = LOG (1-Y)

      pnsm_2_nf0=714.1*DL1 + 1860.2 - 3505*y + 297*Y**2 
     -  - 433.2*y**3 + dl*dl1*(684 + 251.2*dl) + 1465*dl
     +  + 399.2*dl**2 + 320./9.*dl**3 + 116./81.*dl**4 

      return 
      end
C------------------
      real*8 function pnsm_2_nf1(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl Phys B688, 101 (2004)]
c  (NF^1, non-singlet, C-odd, regular piece) 

      DL  = LOG (Y)
      DL1 = LOG (1-Y)

      pnsm_2_nf1=-5120./81.*dl1-216.62+406.5*y+77.89*y**2+34.76*y**3
     -            -1.136*y*dl**3-65.43*dl*dl1
     -            -172.69*dl-3216./81.*dl**2-256./81.*dl**3

      return 
      end
C------------------
      real*8 function pnsm_2_nf2(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl Phys B688, 101 (2004)]
c  (NF^2, non-singlet, C-odd, regular piece) 

      DL  = LOG (Y)

      pnsm_2_nf2=-(- 32.* Y*DL/(1.-Y) * (3.* DL + 10.) - 64.
     -               - (48.* DL**2 + 352.* DL + 384.) * (1.-Y) )/81.D0

      return 
      end
C------------------
      real*8 function ppsp_2_nf0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF-independent, pure-singlet, C-even, regular piece) 

      DL  = LOG (Y)
      DL1 = LOG (1-Y)

      ppsp_2_nf0=(-5.926*dl1**3 - 9.751*dl1**2 - 72.11*dl1 + 177.4
     + + 392.9*y - 101.4*y**2 - 57.04*dl*dl1 - 661.6*dl + 131.4*dl**2
     - - 400./9.*dl**3 + 160./27.*dl**4 - 506./y - 3584./27./y*dl)*(1-y)

      return 
      end
C------------------
      real*8 function ppsp_2_nf1(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF^1, pure-singlet, C-even, regular piece) 

      DL  = LOG (Y)
      DL1 = LOG (1-Y)

      ppsp_2_nf1=(1.778*dl1**2 + 5.944*dl1 + 100.1 - 125.2*y
     +   + 49.26*y**2 - 12.59*y**3 - 1.889*dl*dl1 + 61.75*dl 
     +   + 17.89*dl**2 + 32./27.*dl**3 + 256./81./y)*(1-y) 

      return 
      end
C------------------
      real*8 function ppsm_2(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (pure-singlet, C-odd, regular piece) 

      DL  = LOG (Y)
      DL1 = LOG (1-Y)

      ppsm_2=(dl1*(-163.9/y - 7.208*y) + 151.49 + 44.51*y
     -   - 43.12*y**2 + 4.82*y**3)*(1-y) 
     +   + dl*dl1*(-173.1 + 46.18*dl) + 178.04*dl + 6.892*dl**2
     +   + 40./27.*(dl**4-2*dl**3)

      return 
      end
C------------------
      real*8 function pgg_2_nf0_local(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO gluon-gluon splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF-independent, local piece) 

      DL1 = LOG (1-y)

      pgg_2_nf0_local=2643.521*DL1 + 4425.894 

      return 
      end
C------------------
      real*8 function pgg_2_nf1_local(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO gluon-gluon splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF^1, local piece) 
c The local term in (-x*pqg_2_nf0)_+
      data gf0/18.80866d0/

      DL1 = LOG (1-y)

      pgg_2_nf1_local=-412.172*DL1 - 528.723 
c the quark-gluon local term subtracted
     +  + gf0

      return 
      end
C------------------
      real*8 function pgg_2_nf2_local(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO gluon-gluon splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF^2, local piece) 
c The local term in (-x*pqg_2_nf1)_+
      data gf1 /-6.056172d0/

      DL1 = LOG (1-y)

      pgg_2_nf2_local=-16./9.*DL1+6.4630
c the quark-gluon local term subtracted
     +  + gf1

      return 
      end
C------------------
      real*8 function pns_2_nf0_local(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl Phys B688, 101 (2004)]
c  (NF-independent, non-singlet, local piece) 

      DL1 = LOG (1-Y)

      pns_2_nf0_local=1174.898*DL1 + 1295.384 

      return 
      end
C------------------
      real*8 function pns_2_nf1_local(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl Phys B688, 101 (2004)]
c  (NF^1, non-singlet, local piece) 
c The local term in (-2x*ppsp_2_nf0)_+ 
      data ffps0 /42.21255d0/  

      DL1 = LOG (1-Y)

      pns_2_nf1_local=-(183.187*DL1 + 173.927)
c the pure-singlet local term subtracted
     +  + ffps0

      return 
      end
C------------------
      real*8 function pns_2_nf2_local(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl Phys B688, 101 (2004)]
c  (NF^2, non-singlet, local piece) 
      include 'CONSTCOM.'
c The local term in (-2x*ppsp_2_nf1)_+ 
      data ffps1 /3.445817d0/

      DL1 = LOG (1-Y)

      pns_2_nf2_local=-(16*DL1 + (51+48*zeta3-80*zeta2))*4./81.D0
c the pure-singlet local term subtracted
     +  + ffps1

      return 
      end
C------------------
      real*8 function pgg_2_nf0_singular(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO gluon-gluon splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF-independent, singular piece) 

      pgg_2_nf0_singular=-2643.521/(1-y)

      return 
      end
C------------------
      real*8 function pgg_2_nf1_singular(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO gluon-gluon splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF-independent, singular piece) 

      pgg_2_nf1_singular=412.172/(1-y)

      return 
      end
C------------------
      real*8 function pgg_2_nf2_singular(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO gluon-gluon splitting function of [Nucl. Phys. B691, 129 (2004)]
c  (NF-independent, singular piece) 

      pgg_2_nf2_singular=16./9./(1-y)

      return 
      end
C------------------
      real*8 function pns_2_nf0_singular(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl Phys B688, 101 (2004)]
c  (NF-independent, non-singlet, singular piece) 

      pns_2_nf0_singular=-1174.898/(1-y)

      return 
      end
C------------------
      real*8 function pns_2_nf1_singular(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl Phys B688, 101 (2004)]
c  (NF^1, non-singlet, singular piece) 

      pns_2_nf1_singular=183.187/(1-y)

      return 
      end
C------------------
      real*8 function pns_2_nf2_singular(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The NNLO quark-quark splitting function of [Nucl Phys B688, 101 (2004)]
c  (NF^2, non-singlet, singular piece) 

      pns_2_nf2_singular=64./81.D0/(1-y)

      return 
      end
C------------------
      subroutine QGSPLIT(IQ,p,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'APSCOM6.'
      include 'CONSTCOM.'

      real*8 p(nflim)

      do jq=1,nflim-1
        p(jq)=0.
      end do

      IF (IQ.EQ.1) THEN
        p(1)=pgg_0(z)
        do jq=2,2*nfc+1
          p(jq)=pgq_0(z)
        end do
      END IF

      IF (IQ.GE.2.and.iq.le.2*nfc+1) THEN
        p(1)=pqg_0(z)
        p(iq)=pqq_0(z)
      END IF

      RETURN
      END
C------------------
      subroutine QGSPLIT0(IQ,p,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'APSCOM6.'
      include 'CONSTCOM.'
c The local terms in (-x*pqg_0)_+
      data gf2 /0.666667d0/

      real*8 p(nflim)

      do jq=1,nflim-1
        p(jq)=0.
      end do

      IF (IQ.EQ.1) THEN
        p(1)=pgg_0_local_nf0(x)
c The piece proportional to the number of initial quarks
     -    - nfc*tr*gf2
      END IF

      IF (IQ.GE.2.and.iq.le.2*nfc+1) THEN
        p(iq)=pqq_0_local(x)
      END IF

      RETURN
      END
C------------------
      subroutine QGSPLITX(IQ,p,Z)
      implicit double precision (A-H,O-Z)

      include 'APSCOM6.'
      include 'CONSTCOM.'

      real*8 p(nflim)

      do jq=1,nflim-1
        p(jq)=0.
      end do

      IF (IQ.EQ.1) THEN
        p(1)=pgg_0_singular(z)
      END IF

      IF (IQ.GE.2.and.iq.le.2*nfc+1) THEN
        p(iq)=pqq_0_singular(z)
      END IF

      RETURN
      END
C------------------
      subroutine QGSPLITNLO(IQ,P,Z)
      implicit double precision (a-h,o-z)

      include 'APSCOM6.'
      include 'CONSTCOM.'

      real*8 p(nflim)

      do jq=1,nflim-1
        p(jq)=0.
      end do

      IF (IQ.EQ.1) THEN
        p(1)=pgg_1_nf0(z)+nfe*pgg_1_nf1(z)
        do jq=2,2*nfc+1
          p(jq)=pgq_1_nf0(z)+nfe*pgq_1_nf1(z)
        end do
      END IF

      IF (IQ.ge.2.and.iq.le.2*nfc+1) THEN
        p(1)=pqg_1(z)

        pvqq=pqqns_1_nf0(z)+nfe*pqqns_1_nf1(z)
        pvqaq=pqaqns_1(z)
        ps=pqqps_1(z)

        do jq=2,2*nfc+1
          p(jq)=ps
        end do
        p(iq)=p(iq)+pvqq

        if (iq.eq.2.or.iq.eq.4.or.iq.eq.6.or.iq.eq.8.or.iq.eq.10) then 
          p(iq+1)=p(iq+1)+pvqaq
        end if
        if (iq.eq.3.or.iq.eq.5.or.iq.eq.7.or.iq.eq.9.or.iq.eq.11) then 
          p(iq-1)=p(iq-1)+pvqaq
        end if
      END IF

      RETURN
      END
C------------------
      subroutine QGSPLITNLO0(IQ,p,X)
      implicit double precision (A-H,O-Z)

      include 'APSCOM6.'
      include 'CONSTCOM.'
c The local terms in (-x*pqg_1)_+
      data gf2,gf3 /1.37037d0, 0.64815d0/
c The local term in (-2x*pqqps_1)_+
      data ffps /0.37037d0/ 

      real*8 p(nflim)

      do jq=1,nflim-1
        p(jq)=0.
      end do

      IF (IQ.EQ.1) THEN
        p(1)=pgg_1_nf0_local(x)+nfe*pgg_1_nf1_local(x)
c The piece proportional to the number of initial quarks
     -    -nfc*(cf*tr*gf2+tr*cg*gf3)
      END IF

      IF ( IQ.GE.2.and.iq.le.2*nfc+1) THEN
        p(iq)=pqq_1_nf0_local(x)+nfe*pqq_1_nf1_local(x)
c The piece proportional to the number of initial quarks
     -    -nfc*cf*tr*2*ffps
      END IF

      RETURN
      END
C------------------
      subroutine QGSPLITNLOX(IQ,p,Z)
      implicit double precision (a-h,o-z)

      include 'APSCOM6.'
      include 'CONSTCOM.'

      real*8 p(nflim)

      do jq=1,nflim-1
        p(jq)=0.
      end do

      IF (IQ.EQ.1) THEN
        p(1)=pgg_1_nf0_singular(z)+nfe*pgg_1_nf1_singular(z)
      END IF

      IF ( IQ.GE.2.and.iq.le.2*nfc+1) THEN
        p(iq)=pqq_1_nf0_singular(z)+nfe*pqq_1_nf1_singular(z)
      END IF

      RETURN
      END
C------------------
      subroutine QGSPLITNNLO(IQ,p,Z,kermod)
      implicit double precision (a-h,o-z)

      include 'APSCOM6.'
      include 'CONSTCOM.'

      real*8 p(nflim)

      do jq=1,nflim-1
        p(jq)=0.
      end do

      IF (IQ.EQ.1) THEN
        p(1)=pgg_2_nf0(z) + nfe*pgg_2_nf1(z) + nfe**2*pgg_2_nf2(z)
        do jq=2,2*nfc+1
          p(jq)=pgq_2_nf0(z) + nfe*pgq_2_nf1(z) + nfe**2*pgq_2_nf2(z)
        end do
      END IF

      IF (IQ.GE.2.and.iq.le.2*nfc+1) THEN
        p(1)=pqg_2_nf0(z) + nfe*pqg_2_nf1(z)
        pplus=pnsp_2_nf0(z) + nfe*pnsp_2_nf1(z) + nfe**2*pnsp_2_nf2(z)
        pminus=pnsm_2_nf0(z) + nfe*pnsm_2_nf1(z) + nfe**2*pnsm_2_nf2(z)
        psplus=ppsp_2_nf0(z) + nfe*ppsp_2_nf1(z)
        psminus=ppsm_2(z)

        pvqq=(pplus+pminus)/2.
        pvqaq=(pplus-pminus)/2.
        psqq=(psplus+psminus)/2.
        psqaq=(psplus-psminus)/2.

        p(iq)=p(iq)+pvqq

        if (iq.eq.2.or.iq.eq.4.or.iq.eq.6.or.iq.eq.8.or.iq.eq.10) then 
          p(iq+1)=p(iq+1)+pvqaq
          do jq=2,2*nfc,2
            p(jq)=p(jq)+psqq
            p(jq+1)=p(jq+1)+psqaq
          end do
        end if
        if (iq.eq.3.or.iq.eq.5.or.iq.eq.7.or.iq.eq.9.or.iq.eq.11) then 
          p(iq-1)=p(iq-1)+pvqaq
          do jq=2,2*nfc,2
            p(jq)=p(jq)+psqaq
            p(jq+1)=p(jq+1)+psqq
          end do
        end if
      END IF

      RETURN
      END
C------------------
      subroutine QGSPLITNNLO0(IQ,p,z,kermod)
      implicit double precision (a-h,o-z)

      include 'APSCOM6.'
      include 'CONSTCOM.'
c The local terms in (-x*pqg_2_nf0)_+ and (-x*pqg_2_nf1)_+
      data gf0,gf1 /18.80866d0, -6.056172d0/
c The local terms in (-2x*ppsp_2_nf0)_+ and (-2x*ppsp_2_nf1)_+
      data ffps0,ffps1 /42.21255d0, 3.445817d0/  

      real*8 p(nflim)

      do jq=1,nflim-1
        p(jq)=0.
      end do

      IF (IQ.EQ.1) THEN
        p(1)=pgg_2_nf0_local(z) + nfe*pgg_2_nf1_local(z)
     +      + nfe**2*pgg_2_nf2_local(z)
c The piece proportional to the number of initial quarks
     -      - nfc*(gf0 + nfe*gf1)
      END IF

      IF ( IQ.GE.2.and.iq.le.2*nfc+1) THEN
        p(iq)=pns_2_nf0_local(z) + nfe*pns_2_nf1_local(z)
     +      + nfe**2*pns_2_nf2_local(z)
c The piece proportional to the number of initial quarks
     -      - nfc*(ffps0 + nfe*ffps1)
      END IF

      RETURN
      END
C------------------
      subroutine QGSPLITNNLOX(IQ,p,Z,kermod)
      implicit double precision (A-H,O-Z)

      include 'APSCOM6.'
      include 'CONSTCOM.'

      real*8 p(nflim)

      do jq=1,nflim-1
        p(jq)=0.
      end do

      IF (IQ.EQ.1) THEN
        p(1)=pgg_2_nf0_singular(z) + nfe*pgg_2_nf1_singular(z)
     +      + nfe**2*pgg_2_nf2_singular(z)
      END IF

      IF (IQ.GE.2.and.iq.le.2*nfc+1) THEN
        p(iq)=pns_2_nf0_singular(z) + nfe*pns_2_nf1_singular(z)
     +      + nfe**2*pns_2_nf2_singular(z)
      END IF

      RETURN
      end
c-------------------------
      real*8 FUNCTION P0GGA(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'

      p0gga=4*cg*(1./Z+Z*(1.-Z)-2)

      return 
      end
c------------------------
      real*8 FUNCTION P0GQA(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'

      p0gqa=2*cf*(1.+(1.-Z)**2)/Z

      return 
      end
c------------------------
      real*8 FUNCTION P0QQA(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'

      p0qqa=-2*cf*(Z+1)

      return 
      end
c------------------------
      real*8 FUNCTION P0QGA(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'

      p0qga=2*tr*(Z**2+(1.-Z)**2)

      return 
      end
c------------------------
      real*8 FUNCTION P0GGB(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'

      p0ggb=4*cg/(1.-Z)


      return 
      end
c------------------------
      real*8 FUNCTION P0QQB(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'

      p0qqb=4*cf/(1.-Z)

      return 
      end
c------------------------
      real*8 FUNCTION P0GGC(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'

      p0ggc=2*(2.*cg*(11./12.+LOG(1.-z))-tr*2./3.*nf)

      return 
      end
c------------------------
      real*8 FUNCTION P0QQC(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'

      p0qqc=2*cf*(2.*LOG(1.-z)+3./2.)

      return 
      end
C------------------
      subroutine QGAMSPLIT(IQ,p,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'APSCOM6.'
      include 'CONSTCOM.'

      real*8 p(nflim)

      do jq=1,nflim-1
        p(jq)=0.
      end do

      IF (IQ.EQ.nflim-1) THEN  
        do jq=2,2*nfc+1
          p(jq)=pgq_0(z)/cf*rcharge(jq)**2
        end do
      END IF

      IF (IQ.GE.2.and.iq.le.2*nfc+1) THEN
        p(nflim-1)=pqg_0(z)/tr*rcharge(iq)**2
        p(iq)=pqq_0(z)/cf*rcharge(iq)**2
      END IF

      RETURN
      END
C------------------
      subroutine QGAMSPLIT0(IQ,p,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'APSCOM6.'
      include 'CONSTCOM.'

      data gf2 /0.666667d0/

      real*8 p(nflim)

      do jq=1,nflim-1
        p(jq)=0.
      end do

      IF (IQ.EQ.nflim-1) THEN
c The piece proportional to the number of initial quarks
        p(nflim-1)=-gf2*qsum(nfc) 
      END IF

      IF (IQ.GE.2.and.iq.le.2*nfc+1) THEN
        p(iq)=pqq_0_local(x)/cf*rcharge(iq)**2
      END IF

      RETURN
      END
C------------------
      subroutine QGAMSPLITX(IQ,p,Z)
      implicit double precision (A-H,O-Z)

      include 'APSCOM6.'
      include 'CONSTCOM.'

      real*8 p(nflim)

      do jq=1,nflim-1
        p(jq)=0.
      end do

      IF (IQ.GE.2.and.iq.le.2*nfc+1) THEN
        p(iq)=pqq_0_singular(z)/cf*rcharge(iq)**2
      END IF

      RETURN
      END

