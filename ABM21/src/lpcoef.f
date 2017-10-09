C------------------
      real*8 function c3nspm_2_0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet C-asymmetric constant term for F_3 
c [Nucl. Phys. B568, 263 (2000)]

      dl0=log(y)
      dl1=log(1-y)

      c3nspm_2_0 = - 396.1*dl1 - 92.43*dl1**2*dl0 - 3.049*dl0**3 
     -  - 30.14*dl0**2 - 79.14*dl0 - 467.2*y - 242.9
     -  - (- 409.6*dl1 - 147.9*dl1**2*dl0 - 3.922*dl0**3 
     -  - 33.31*dl0**2 - 67.6*dl0 - 576.8*y - 206.1)

      return 
      end
C------------------
      real*8 function clnspm_2_0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet C-asymmetric constant term for F_L 
c [Nucl. Phys. B568, 263 (2000)]

      dl0=log(y)
      dl1=log(1-y)

      clnspm_2_0 = (13.62*dl1**2 - 55.79*dl1 - 150.5*dl1*dl0
     +  +(26.56*y - 0.031)*dl0**2 - 14.85*dl0 + 97.48*y - 40.41)
     -  -(13.30*dl1**2 - 59.12*dl1 - 141.7*dl1*dl0
     +  +(23.29*y - 0.043)*dl0**2 - 22.21*dl0 + 100.8*y - 52.27)

      return 
      end
C------------------
      real*8 function c2nspm_2_0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet C-asymmetric constant term for F_2
c [Nucl. Phys. B568, 263 (2000)]

      dl0=log(y)
      dl1=log(1-y)
      c2nspm_2_0 = (- 69.59 - 1008*y - 660.7 * dL1
     -          - 2.835 * dL0**3 - 17.08 * dL0**2 + 5.986 * dL0 
     -          - 174.8 * dL0 * dL1**2 + 95.09 * dL0**2 * dL1)
     -          -(- 84.18 - 1010*y - 663.0 * dL1
     -          - 3.748 * dL0**3 - 19.56 * dL0**2 - 1.235 * dL0 
     -          - 192.4 * dL0 * dL1**2 + 80.41 * dL0**2 * dL1)

      return 
      end
C------------------
      real*8 function clps_3_0(y,nfi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) pure-singlet constant term for F_L
c [Nucl. Phys. B724, 3 (2005)]

      include 'CONSTCOM.'

      clps_3_0=clps_3_0_nf0(y) + nfi*clps_3_0_nf1(y)
     +  + (qsum0(nfi)**2/qsum(nfi)-3*qsum0(nfi))/nfi*clps_3_0_nfch(y)

      return 
      end
C------------------
      real*8 function clg_3_0(y,nfi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) gluonic constant term for F_L
c [Nucl. Phys. B724, 3 (2005)]

      include 'CONSTCOM.'

      clg_3_0=clg_3_0_nf0(y) + nfi*clg_3_0_nf1(y)
     +  + qsum0(nfi)**2/qsum(nfi)*clg_3_0_nfch(y)

      return 
      end
C------------------
      real*8 function clns_3_0(y,nfi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) non-singlet constant term for F_L
c [Nucl. Phys. B724, 3 (2005)]

      include 'CONSTCOM.'

      clns_3_0=clns_3_0_nf0(y) + nfi*clns_3_0_nf1(y)
     +  + nfi**2*clns_3_0_nf2(y) + 3*qsum0(nfi)*clns_3_0_nfch(y)

      return 
      end
C------------------
      real*8 function c2g_1_0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S) gluonic constant term of [Nucl. Phys. B383, 525 (1992)]  

      include 'CONSTCOM.'

      c2g_1_0=tr*4*((y**2+(1.-y)**2)*log((1.-y)/y)-1.+8.*y*(1.-y))

      return 
      end
C------------------
      real*8 function c2ns_1_0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S) non-singlet constant term of [Nucl. Phys. B383, 525 (1992)]  

      include 'CONSTCOM.'

      c2ns_1_0=cf*(-2*(1+y)*log(1-y)-2*(1+y**2)/(1-y)*log(y)+6+4*y)

      return 
      end
C------------------
      real*8 function c2ns_1_0_singular(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S) non-singlet constant term of [Nucl. Phys. B383, 525 (1992)]  
c (singular piece)

      include 'CONSTCOM.'

      c2ns_1_0_singular=cf*(4.*log(1-y)-3)/(1-y)

      return 
      end
C------------------
      real*8 function c2ns_1_0_local(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S) non-singlet constant term of [Nucl. Phys. B383, 525 (1992)]  
c (local piece)

      include 'CONSTCOM.'

      dl1=log(1-y)
      c2ns_1_0_local=cf*(2*dl1**2-3*dl1-4*zeta2-9)

      return 
      end
C------------------
      real*8 function c3ns_1_0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S) non-singlet constant term of [ZPC 11, 293 (1982)]  

      include 'CONSTCOM.'

      c3ns_1_0=c2ns_1_0(y)-2*cf*(1+y)

      return 
      end
C------------------
      real*8 function c3ns_1_0_singular(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S) non-singlet constant term of [ZPC 11, 293 (1982)]
c (singular piece)

      c3ns_1_0_singular=c2ns_1_0_singular(y)

      return 
      end
C------------------
      real*8 function c3ns_1_0_local(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S) non-singlet constant term of [ZPC 11, 293 (1982)]
c (local piece)

      c3ns_1_0_local=c2ns_1_0_local(y)

      return 
      end
C------------------
      real*8 function c2ps_2_0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) pure-singlet constant term of [Nucl. Phys. B588, 345 (2000)]

      dl0=log(z)
      dl1=log(1-z)

      c2ps_2_0=5.290 * (1./z-1.) + 4.310 * dL0**3
     -         - 2.086 * dL0**2 + 39.78 * dL0 - 0.101 * (1.-z) * dL1**3
     -         - (24.75 - 13.80 * z) * dL0**2 * dL1 + 30.23 * dL0 * dL1

      return 
      end
C------------------
      real*8 function c2ps_2_axial_0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) pure-singlet constant axial term of Zijlstra's thesis

      include 'CONSTCOM.'

      dl0=log(z)
      dl1=log(1-z)

      c2ps_2_axial_0=cf*(32*(3*z**3 - 3*z**2 + z)*ddilog(1.-z) 
     -  -32*(3./5.*z**3+1./3.*z - 1./15./z**2)
     *  *(ddilog(-z)+dl0*log(1+z))
     +  +16*(18./5.*z**3 - 3*z**2 + 4./3.*z)*dl0**2
     +  +32*(12./5.*z**3 - 3*z**2 + 2./3.*z)*zeta2
     +  +16*(1./15. - 24./5.*z**2 + 12./5.*z - 2./15./z
     +  +1./(1-z) - 1./(1+z))*dl0
     -  -8./5.*(1./3. + 48*z**2 - 42*z - 4./3./z)) 

      return 
      end
C------------------
      real*8 function c2g_2_0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) gluonic constant term of [Nucl. Phys. B588, 345 (2000)]

      dl0=log(z)
      dl1=log(1-z)

      c2g_2_0=(1./z * (11.90 + 1494.* dL1) + 5.319 * dL0**3
     -         - 59.48 * dL0**2 - 284.8 * dL0 + 392.4 - 1483.* dL1
     +         + (6.445 + 209.4 * (1.-z)) * dL1**3 - 24.00 * dL1**2
     -         - 724.1 * dL0**2 * dL1 - 871.8 * dL0 * dL1**2)

      return 
      end
C------------------
      real*8 function c2g_2_0_local(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) gluonic constant term of [Nucl. Phys. B588, 345 (2000)]
c  (local piece)

      c2g_2_0_local=-0.28

      return 
      end
C------------------
      real*8 function c2nsp_2_0_nf0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Nucl. Phys. B568, 263 (2000)]
c  (NF-independent piece)

      dl0=log(z)
      dl1=log(1-z)

      c2nsp_2_0_nf0=- 69.59 - 1008.* z
     -          - 2.835 * dL0**3 - 17.08 * dL0**2 + 5.986 * dL0 
     -          - 17.19 * dL1**3 + 71.08 * dL1**2 - 660.7 * dL1
     -          - 174.8 * dL0 * dL1**2 + 95.09 * dL0**2 * dL1

      return 
      end 
C------------------
      real*8 function c2nsp_2_0_nf0_singular(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Nucl. Phys. B568, 263 (2000)]
c  (NF-independent, singular piece)

      dl1=log(1-z)
      c2nsp_2_0_nf0_singular=(14.2222*dl1**3-61.3333*dl1**2
     -   -31.105*dl1+188.64)/(1-z)

      return 
      end
C------------------
      real*8 function c2nsp_2_0_nf0_local(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Nucl. Phys. B568, 263 (2000)]
c  (NF-independent, local piece)

      dl1=log(1-z)
      c2nsp_2_0_nf0_local=-338.046+14.2222*dl1**4/4.
     -       -61.3333*dl1**3/3.-31.105*dl1**2/2.+188.64*dl1

      return 
      end
C------------------
      real*8 function c2nsp_2_0_nf1(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Nucl. Phys. B568, 263 (2000)]
c  (NF^1 piece)

      dl0=log(z)
      dl1=log(1-z)
      c2nsp_2_0_nf1= - 5.691 - 37.91 * z 
     +          + 2.244 * dL0**2 + 5.770 * dL0 
     -          - 1.707 * dL1**2  + 22.95 * dL1
     +          + 3.036 * dL0**2 * dL1 + 17.97 * dL0 * dL1 

      return 
      end
C------------------
      real*8 function c2nsp_2_0_nf1_singular(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Nucl. Phys. B568, 263 (2000)]
c  (NF^1, singular piece)

      dl1=log(1-z)
      c2nsp_2_0_nf1_singular=(1.77778*dl1**2-8.5926*dl1+6.3489)/(1-z) 

      return 
      end
C------------------
      real*8 function c2nsp_2_0_nf1_local(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Nucl. Phys. B568, 263 (2000)]
c  (NF^1, local piece)

      dl1=log(1-z)
      c2nsp_2_0_nf1_local=46.8405+1.77778*dl1**3/3.
     -     -8.5926*dl1**2/2.+6.3489*dl1

      return 
      end
C------------------
      real*8 function c3nsm_2_0_nf0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Nucl. Phys. B568, 263 (2000)]
c  (NF-independent piece)

      dl0=log(z)
      dl1=log(1-z)

      c3nsm_2_0_nf0=-15.20*dl1**3 + 94.61*dl1**2 - 409.6*dl1 
     -      - 147.9*dl1**2*dl0 - 3.922*dl0**3 - 33.31*dl0**2
     -      - 67.6*dl0 - 576.8*z - 206.1

      return 
      end 
C------------------
      real*8 function c3nsm_2_0_nf0_singular(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Nucl. Phys. B568, 263 (2000)]
c  (NF-independent, singular piece)

      dl1=log(1-z)
      c3nsm_2_0_nf0_singular=(14.2222*dl1**3-61.3333*dl1**2
     -   -31.105*dl1+188.64)/(1-z)

      return 
      end
C------------------
      real*8 function c3nsm_2_0_nf0_local(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Nucl. Phys. B568, 263 (2000)]
c  (NF-independent, local piece)

      dl1=log(1-z)
      c3nsm_2_0_nf0_local=-338.635+14.2222*dl1**4/4.
     -       -61.3333*dl1**3/3.-31.105*dl1*2/2.+188.64*dl1

      return 
      end
C------------------
      real*8 function c3nsm_2_0_nf1(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Nucl. Phys. B568, 263 (2000)]
c  (NF-independent piece)

      dl0=log(z)
      dl1=log(1-z)

      c3nsm_2_0_nf1=0.042*dl1**3 - 0.808*dl1**2 + 25*dl1 
     +  + 9.684*dl1*dl0 + 2.207*dl0**2 + 8.683*dl0 - 14.97*z - 6.337

      return 
      end 
C------------------
      real*8 function c3nsm_2_0_nf1_singular(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Nucl. Phys. B568, 263 (2000)]
c  (NF^1, singular piece)

      dl1=log(1-z)
      c3nsm_2_0_nf1_singular=(1.77778*dl1**2-8.5926*dl1+6.3489)/(1-z) 

      return 
      end
C------------------
      real*8 function c3nsm_2_0_nf1_local(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Nucl. Phys. B568, 263 (2000)]
c  (NF^1, local piece)

      dl1=log(1-z)
      c3nsm_2_0_nf1_local=46.8405+1.77778*dl1**3/3.
     -     -8.5926*dl1**2/2.+6.3489*dl1

      return 
      end
C------------------
      real*8 function clns_1_0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S) non-singlet constant term of [Nucl. Phys. B383, 525 (1992)]  

      include 'CONSTCOM.'

      clns_1_0=cf*4*y

      return 
      end
C------------------
      real*8 function clg_1_0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S) gluonic constant term of [Nucl. Phys. B383, 525 (1992)]  

      include 'CONSTCOM.'

      clg_1_0=tr*16*y*(1-y)

      return 
      end
C------------------
      real*8 function clns_2_0_nf0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Phys. Lett. B606, 123 (2005)]
c  (NF-independent, regular piece)

      include 'CONSTCOM.'

      dl0=log(z)
      dl1=log(1-z)

      clns_2_0_nf0=128./9.*z*dl1**2-46.5*z*dl1-84.094*dl0*dl1
     -   -37.338+89.53*z+33.82*z**2+z*dl0*(32.9+18.41*dl0)-128./9.*dl0

      return 
      end
C------------------
      real*8 function clns_2_0_nf0_local(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Phys. Lett. B606, 123 (2005)]
c  (NF-independent, local piece)

      include 'CONSTCOM.'

      clns_2_0_nf0_local=-0.012

      return 
      end
C------------------
      real*8 function clns_2_0_nf1(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) non-singlet constant term of [Phys. Lett. B606, 123 (2005)]
c  (NF^1 piece)

      include 'CONSTCOM.'

      dl0=log(z)
      dl1=log(1-z)

      clns_2_0_nf1=16./27.*(6*z*dl1-12*z*dl0-25*z+6)

      return 
      end
C------------------
      real*8 function clg_2_0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) gluon constant term of [Phys. Lett. B606, 123 (2005)]

      include 'CONSTCOM.'

      dl0=log(z)
      dl1=log(1-z)

      clg_2_0=(94.74-49.2*z)*(1-z)*dl1**2+864.8*(1-z)*dl1
     + +1161*z*dl0*dl1+60.06*z*dl0**2+39.66*(1-z)*dl0-5.333*(1./z-1)

      return
      end
C------------------
      real*8 function clps_2_0(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^2) pure-singlet constant term of [Phys. Lett. B606, 123 (2005)]

      include 'CONSTCOM.'

      dl0=log(z)
      dl1=log(1-z)

      clps_2_0=(15.94 - 5.212 * z) * (1-z)**2 * dL1
     +         + (0.421 + 1.520 * z) * dL0**2 + 28.09 * (1-z) * dL0
     -         - (2.370/z - 19.27) * (1-z)**3

      return
      end
C------------------
      real*8 function clns_3_0_nf0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) non-singlet constant term of [Phys. Lett. B606, 123 (2005)]
c  (NF-independent, regular piece)

      include 'CONSTCOM.'

      dl0=log(y)
      dl1=log(1-y)

      clns_3_0_nf0=512./27.*dl1**4 - 177.4*dl1**3 + 650.6*dl1**2 
     -  - 2729*dl1 - 2220.5 - 7884*y + 4168*y**2 
     -  - (844.7*dl0 + 517.3*dl1)*dl0*dl1
     +  +(195.6*dl1 - 125.3)*(1-y)*dl1**3 + 208.3*y*dl0**3 - 1355.7*dl0
     -  - 7456./27.*dl0**2 - 1280./81.*dl0**3
 
      return 
      end
C------------------
      real*8 function clns_3_0_nf0_local(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) non-singlet constant term of [Phys. Lett. B606, 123 (2005)]
c  (NF-independent, local piece)

      include 'CONSTCOM.'

      clns_3_0_nf0_local=0.113

      return 
      end
C------------------
      real*8 function clns_3_0_nf1(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) non-singlet constant term of [Phys. Lett. B606, 123 (2005)]
c  (NF^1, regular piece)

      include 'CONSTCOM.'

      dl0=log(y)
      dl1=log(1-y)

      clns_3_0_nf1=1024./81.*dl1**3 - 112.35*dl1**2 + 344.1*dl1 
     +  + 408.4 - 9.345*y - 919.3*y**2  
     +  + (239.7 + 20.63*dl1)*(1-y)*dl1**2
     +  + dl0*dl1*(887.3 + 294.5*dl0 - 59.14*dl1)
     -  - 1792./81.*y*dl0**3 + 200.73*dl0 + 64./3.*dl0**2
 
      return 
      end
C------------------
      real*8 function clns_3_0_nf1_local(z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) non-singlet constant term of [Phys. Lett. B606, 123 (2005)]
c  (NF^1, local piece)

      include 'CONSTCOM.'

      clns_3_0_nf1_local=0.006

      return 
      end
C------------------
      real*8 function clns_3_0_nf2(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) non-singlet constant term of [Phys. Lett. B606, 123 (2005)]
c  (NF^2 piece)

      include 'CONSTCOM.'

      dl0=log(y)
      dl1=log(1-y)

      clns_3_0_nf2=(3*y*dl1**2 + (6-25*y)*dl1 - 19 
     +  + (317./6.-12*zeta2)*y 
     -  - 6*y*dl0*dl1 + 6*y*ddilog(y) + 9*y*dl0**2 
     -  - (6-50*y)*dl0)*64./81.

      return 
      end
C------------------
      real*8 function clns_3_0_nfch(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) non-singlet constant term of [Phys. Lett. B606, 123 (2005)]
c  (charge-factor piece)

      include 'CONSTCOM.'

      dl0=log(y)
      dl1=log(1-y)

      clns_3_0_nfch=y*((107 + 321.05*y - 54.62*y**2)*(1-y) 
     -  - 26.717 + 9.773*dl0
     +  + (363.8 + 68.32*dl0)*y*dl0 - 320./81.*dl0**2*(2+dl0))

      return 
      end
C------------------
      real*8 function clg_3_0_nf0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) gluon constant term of [Phys. Lett. B606, 123 (2005)]
c  (NF-independent piece)

      include 'CONSTCOM.'

      dl0=log(y)
      dl1=log(1-y)

      clg_3_0_nf0=(144*dl1**4 - 47024./27.*dl1**3 + 6319*dl1**2 
     +  + 53160*dl1)*(1-y) + 72549*dl0*dl1
     +  + 88238*dl0**2*dl1 + (3709 - 33514*y - 9533*y**2)*(1-y)
     +  + 66773*y*dl0**2 - 1117*dl0 + 45.37*dl0**2 - 5360./27.*dl0**3
     -  - (2044.7*(1-y) + 409.506*dl0)/y
 
      return 
      end
C------------------
      real*8 function clg_3_0_nf1(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) gluon constant term of [Phys. Lett. B606, 123 (2005)]
c  (NF^1 piece)

      include 'CONSTCOM.'

      dl0=log(y)
      dl1=log(1-y)

      clg_3_0_nf1=(32./3.*dl1**3 - 1216./9.*dl1**2 - 592.3*dl1 
     +  + 1511*y*dl1)*(1-y) + 311.3*dl0*dl1 + 14.24*dl0**2*dl1 
     +  + (577.3 - 729*y)*(1-y) + 30.78*y*dl0**3 + 366*dl0
     +  + 1000./9.*dl0**2 + 160./9.*dl0**3 + 88.5037/y*(1-y)
 
      return 
      end
C------------------
      real*8 function clg_3_0_nfch(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) gluon constant term of [Phys. Lett. B606, 123 (2005)]
c  (charge-factor piece)

      include 'CONSTCOM.'

      dl0=log(y)
      dl1=log(1-y)

      clg_3_0_nfch=(-0.0105*dl1**3 + 1.55*dl1**2 + 19.72*y*dl1 
     -  - 66.745*y + 0.615*y**2)*(1-y) + 20./27.*y*dl0**4 
     +  + (280./81. + 2.26*y)*y*dl0**3 - (15.4 - 2.201*y)*y*dl0**2
     -  - (71.66 - 0.121*y)*y*dl0
 
      return 
      end
C------------------
      real*8 function clps_3_0_nf0(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) pure-singlet constant term of [Phys. Lett. B606, 123 (2005)]
c  (NF-independent piece)

      include 'CONSTCOM.'

      dl0=log(y)
      dl1=log(1-y)

      clps_3_0_nf0=(1568./27.*dl1**3 - 3968./9.*dl1**2 + 5124*dl1)
     *  * (1-y)**2
     +  + (2184*dl0 + 6059*(1-y))*dl0*dl1 - (795.6 + 1036*y)*(1-y)**2
     -  - 143.6*(1-y)*dl0 + 2848./9.*dl0**2 - 1600./27.*dl0**3
     -  - (885.53*(1-y) + 182*dl0)/y*(1-y)
 
      return 
      end
C------------------
      real*8 function clps_3_0_nf1(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) pure-singlet constant term of [Phys. Lett. B606, 123 (2005)]
c  (NF^1 piece)

      include 'CONSTCOM.'

      dl0=log(y)
      dl1=log(1-y)

      clps_3_0_nf1=(-32./9.*dl1**2 + 29.52*dl1)*(1-y)**2 
     +  + (35.18*dl0 + 73.06*(1-y))*dl0*dl1 - 35.24*y*dl0**2
     -  - (14.16 - 69.84*y)*(1-y)**2 - 69.41*(1-y)*dl0 
     -  - 128./9.*dl0**2 + 40.239/y*(1-y)**2
 
      return 
      end

C------------------
      real*8 function clps_3_0_nfch(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c The O(\ALPHA_S^3) pure-singlet constant term of [Phys. Lett. B606, 123 (2005)]
c  (charge-factor piece)

      include 'CONSTCOM.'

      dl0=log(y)
      dl1=log(1-y)

      clps_3_0_nfch=y*(107 + 321.05*y - 54.62*y**2)*(1-y) 
     -  - 26.717 + 9.773*dl0
     +  + (363.8 + 68.32*dl0)*y*dl0 - 320./81.*dl0**2*(2+dl0)
 
      return 
      end
C------------------
      subroutine cintx(kint,iord,ifun,inf,xx,cc)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'APSCOM6.'
      include 'CONSTCOM.'

      real*8 cc(3)

      xb=min(xx,xbcmax)
      xl1=log(xb)
      xl2=log(1-xb)

      if (xb.ge.xc1) then
        IX=int((xclog2-xl2)/delcp)
        DXX=(xclog2-xl2)/delcp-ix
        if (ix.eq.0) then 
          do iq=1,3
            cc(iq)=cgrid(kint,iord,ifun,iq,inf,0)*(1-dxx)
     +            +cgrid(kint,iord,ifun,iq,inf,1)*dxx
          end do
        else 
          do iq=1,3
            cc(iq)=cgrid(kint,iord,ifun,iq,inf,ix-1)*(dxx-1)*dxx/2.
     +            +cgrid(kint,iord,ifun,iq,inf,ix)*(1-dxx**2)
     +            +cgrid(kint,iord,ifun,iq,inf,ix+1)*(1+dxx)*dxx/2.
          end do
        end if          
      else
        IX=int((XL1-XcLOG1)/DELcm)-1
        DXX=(xl1-xclog1)/delcm-ix
        do iq=1,3
          cc(iq)=cgrid(kint,iord,ifun,iq,inf,ix-1)*(dxx-1)*dxx/2.
     +          +cgrid(kint,iord,ifun,iq,inf,ix)*(1-dxx**2)
     +          +cgrid(kint,iord,ifun,iq,inf,ix+1)*(1+dxx)*dxx/2.
        end do
      end if

      RETURN
      END
C-----------------
      SUBROUTINE cgridini
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'APSCOM6.'
      include 'CONSTCOM.'

      xc1=0.2d0
      xbcmax=0.999d0

      XLOG0=LOG(XBmin)
      xclog1=log(xc1)
      xclog2=log(1-xc1)

      DELcp=(xclog2-log(1-xbcmax))/(nxpgrid-1)
      DELcm=-(XLOG0-xclog1)/nxmgrid

      DO I=0,nxpgrid-1
        xcgrid(I)=1-exp(xclog2-delcp*i)
      end do
      xcgrid(nxpgrid)=1.

      DO I=-nxmgrid,0
        XcGRID(I)=EXP(XcLOG1+DELcm*I)
      end do

      do nfi=3,6
        DO IX=-nxmgrid,nxpgrid-1
          z=xcgrid(ix)
          cgrid(1,2,2,1,nfi,ix)=clns_3_0(z,nfi)
          cgrid(1,2,2,2,nfi,ix)=clg_3_0(z,nfi) 
          cgrid(1,2,2,3,nfi,ix)=clps_3_0(z,nfi)
        end do
      end do

      RETURN
      END

