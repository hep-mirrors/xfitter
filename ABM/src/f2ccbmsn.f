c------------------
      real*8 function f2h_bmsn(ni,nb,nt,xb,q2,nq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'CONSTCOM.'
      include 'PDFCOM.'

      if (nq.eq.8) f2ffn=f2charm_ffn(xb,q2,nq)
      if (nq.eq.10) f2ffn=f2charm_ffn4(xb,q2,nq)

! Save current status of KORDHQ  
      kordhqs=kordhq 

! BMSN prescription is ready for KORDHQ=1 only 
      kordhq=1
      f2h_bmsn=f2h_vfn(ni,nb,nt,xb,q2,nq)
     -   - f2h_asymp(ni,nb,nt,xb,q2,nq) + f2ffn

! Restore status of KORDHQ  
      kordhq=kordhqs 

      return 
      end
C------------------
      real*8 function f2h_vfn(nb,nt,ni,xb,q2,nq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      external f2h_vfni

      include 'CONSTCOM.'
      include 'PDFCOM.'
      include 'PRECCOM.'

      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0

!  The charge of heavy quark
      qq=rcharge(nq)

      if (kordhq.ne.1) 
     -    print *,'F2H_VFN IS NOT READY FOR KORDHQ=',kordhq
      
!  Save the current scheme status
      kschemepdfs=kschemepdf
      kschemepdf=0           ! Take 4-, 5- flavour quark distribution 
                             ! without O(\alpha_s) corrections. 
      ischem0=-(nq-8)-2      ! LO 4-, 5-flavour PDFs
      ischem1=-(nq-8)-3      ! NLO 4-, 5-flavour PDFs

      an=xqg(0,0.1d0,q2,ischem0)/4./pi

      f20=2*xqg(nq,xb,q2,ischem1)   
      f2h_vfn=f20*qq**2

      q2ss=q2
      nq0=nq
      qqs=qq
      xb0=xb
      nb0=nb
      nt0=nt
      ni0=ni


      CALL GAUSS1(f2h_vfni,log(xb),0d0,nf2hq,df,EPS)
      f2h_vfn=f2h_vfn+df

      if (kordhq.ge.1) then 
        f2h_vfn=f2h_vfn+an*c2ns_1_0_local(xb)
     *  *xqg(nq,xb,q2,ischem0)*qq**2  
        f2h_vfn=f2h_vfn+an**2*c2g_2_0_local(y)
     *  *xqg(1,xb,q2,ischem0)*qq**2  
        r=q2/rmass(nq)**2  
        f2h_vfn=f2h_vfn + an**2*(ome_qqns_2_local(xb,r)
     +     + c2nsp_2_0_nf1_local(xb))*f2cqqpm(nb,nt,ni,xb,q2)   
      end if

!  Restore the initial scheme
      kschemepdf=kschemepdfs

      return
      end
C------------------
      real*8 function f2h_vfni(t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      INCLUDE 'PDFCOM.' 

      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0

      ischem0=-(nq0-8)-2      ! LO 4-, 5-flavour PDFs
      ischem1=-(nq0-8)-3      ! NLO 4-, 5-flavour PDFs

      y=exp(t)
      z=xb0/y

      f2ns0=0.
      f2ns=0.
      f2g=an*c2g_1_0(y)*xqg(1,z,q2ss,ischem0)     
      f2ps=0.

      f2ns=c2ns_1_0(y)*xqg(nq0,z,q2ss,ischem0)*an*qqs**2   
      f2ns0=c2ns_1_0_singular(y)*an*qqs**2
     *     *(xqg(nq0,z,q2ss,ischem0)-xqg(nq0,xb0,q2ss,ischem0))     

      if (kordhq.ge.1) then 
        f2g=f2g+c2g_2_0(y)*an**2*xqg(1,z,q2ss,kschemepdf)   
        f2ps=c2ps_2_0(y)*an**2*xpscqqpm(z,q2ss)*qqs**2  

        r=q2ss/rmass(nq0)**2  
        f2ns=f2ns+(ome_qqns_2(y,r)+c2nsp_2_0_nf1(y))
     *      *f2cqqpm(nb0,nt0,ni0,z,q2ss)*an**2
        f2ns0=f2ns0+(ome_qqns_2_singular(y,r)+c2nsp_2_0_nf1_singular(y))
     * *(f2cqqpm(nb0,nt0,ni0,z,q2ss)
     - - f2cqqpm(nb0,nt0,ni0,xb0,q2ss))*an**2
      end if

      f2g=f2g*qqs**2 
      f2h_vfni=(f2g+f2ns0+f2ps+f2ns)*y

      return 
      end
c------------------
      real*8 function f2h_asymp(nb,nt,ni,xb,q2,nq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'
      include 'PRECCOM.'

      external f2h_asympi1,f2h_asympi2
      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0

!  The mass and charge of heavy quark
      rm2=rmass(nq)**2
      qq=rcharge(nq)

      if (kordhq.ge.2) 
     -   print *,'F2H_ASYMP IS NOT READY FOR KORDHQ=',kordhq

!  Save the current scheme status
      kschemepdfs=kschemepdf
      kschemepdf=0           ! take 4-, 5- flavour quark distribution 
                             ! without O(\alpha_s) corrections. 
      ischem0=-(nq-8)        ! 3-, and 4-flavour PDFs w.o. O(\alpha_s) corr. 

      rmu2=q2
      an=xqg(0,0.1d0,rmu2,ischem0)/4./pi

      xb0=xb
      q2ss=q2
      nq0=nq
      qqs=qq
      nb0=nb
      nt0=nt
      ni0=ni

      a=1./(1.+4.*rm2/q2)

      CALL GAUSS1(f2h_asympi2,log(xb),0d0,nf2hq,f2c,EPS)

      f2h_asymp=f2c
      if (kordhq.ge.1) then 
        f2h_asymp=f2h_asymp-an**2*0.28*xqg(1,xb,q2,ischem0)*qq**2     

        xqsumc0=f2cqqpm(nb,nt,ni,xb,rmu2)    

        dl1=log(1-xb)

        f2h_asymp=f2h_asymp+an**2*cf*tr*xqsumc0
     *   * ( log(q2/rm2)**2*2*(4./3.*dl1+1)
     +   + log(q2/rm2)*(2*(8./3.*dl1**2/2.-58./9.*dl1)
     -   - (32./3.*zeta2+38./3))
     +   + 2*(-8./3.*zeta2*dl1 + 4./3.*dl1**3/3. - 58./9.*dl1**2/2. 
     +   + 359./27.*dl1) + 268./9.*zeta2 + 265./9.)

c The second-order coefficinet for the non-singlet OME  
c BMSN EPJ C1, 301 (1998)
c        f2h_asymp=f2h_asymp+cf*tr*xqsumc0
c     *   *(224./27.*dl1 - 8./3.*zeta3 + 40./9.*zeta2 
c     +   + 73./18.)
c Approximate formula for the N_f-coefficient in C_2^NS,NNLO
c van Neerven-Vogt Nucl.Phys. B568, 263 (2000). 
c        f2h_asymp=f2h_asymp+xqsumc0
c     *   *(1.77778*dl1**3/3.-8.5926*dl1**2/2.+6.3489*dl1+46.8405)
c Exact formula for the N_f-coefficient in C_2^NS,NNLO
c van Neerven-Zijlstra Phys.Lett. B272, 127 (1991).
c        f2h_asymp=f2h_asymp+cf*xqsumc0
c     *    * (2*(2./3.*dl1**3/3.-29./9.*dl1**2/2.
c     -    - 4./3.*zeta2*dl1+247./54.*dl1)
c     +    + 4./3.*zeta3+38./3.*zeta2+457./36.)
      end if

!  Restore the initial scheme
      kschemepdf=kschemepdfs

      return
      end
C------------------
      real*8 function f2h_asympi1(t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      y=1.-exp(t)

      f2h_asympi1=f2h_asympi(y)*(1-y)

      return 
      end
C------------------
      real*8 function f2h_asympi2(t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      y=exp(t)

      f2h_asympi2=f2h_asympi(y)*y
 
      return 
      end
C------------------
      real*8 function f2h_asympi(t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      INCLUDE 'PDFCOM.' 

      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0
      complex*16 WGPLG

      z=(t)
      y=xb0/z
      dlz=log(z)
      dlm=log(1d0-z)

      ischem0=-(nq0-8)      ! 3-, and 4-flavour PDFs w.o. O(\alpha_s) corr. 
     
      xi=q2ss/rm2

      f2h_asympi=an*ch2g00_asymp(z,xi)*xqg(1,y,rmu2,ischem0)*qqs**2   

      if (kordhq.ge.1) then 
        xqsum=xpscqqpm(y,rmu2)                     
        xqsumc=f2cqqpm(nb0,nt0,ni0,y,rmu2)         
        xqsumc0=f2cqqpm(nb0,nt0,ni0,xb0,rmu2)      

        f2h_asympi=f2h_asympi+an**2*((ch2g10_asymp(z,xi)
     +   +  ch2g11_asymp(z,xi)*log(rmu2/rm2))*xqg(1,y,rmu2,0)  
     +   +  (ch2qps10_asymp(z,xi)
     +   +  ch2qps11_asymp(z,xi)*log(rmu2/rm2))*xqsum)*qqs**2

        f2h_asympi=f2h_asympi+an**2*cf*tr*(-4./3.*(1+z)*log(xi)**2
     +   + ((1+z**2)/(1-z)*(-16./3.*dlz)
     -   - (1+z)*(8./3.*dlm-58./9.) + 2./3.+26./3.*z)*log(xi)
     +   + (1+z**2)/(1-z)*(-8./3.*ddilog(1-z)
     -   - 16./3.*dlz*dlm + 4*dlz**2 + 134./9.*dlz)
     -   - (1+z)*(-8./3.*zeta2 + 4./3.*dlm**2 - 58./9.*dlm + 359./27.)
     +   + (2./3.+26./3.*z)*dlm - (2+46./3.*z)*dlz + 29./9. 
     -   - 295./9.*z)*xqsumc
 
c The second-order coefficinet for the non-singlet OME  
c BMSN EPJ C1, 301 (1998)
c        f2h_asympi=cf*tr*((1+z**2)/(1-z)*(2./3.*dlz**2+20./9.*dlz) 
c     +  + 8./3.*(1-z)*dlz + 44./27.-268./27.*z)*xqsumc
c Approximate formula for the N_f-coefficient in C_2^NS,NNLO
c van Neerven-Vogt Nucl.Phys. B568, 263 (2000). 
c        f2h_asympi=( - 5.691 - 37.91 * z 
c     +          + 2.244 * dLz**2 + 5.770 * dLz 
c     -          - 1.707 * dLm**2  + 22.95 * dLm
c     +          + 3.036 * dLz**2 * dLm + 17.97 * dLz * dLm )*xqsumc
c Exact formula for the N_f-coefficient in C_2^NS,NNLO
c van Neerven-Zijlstra Phys.Lett. B272, 127 (1991).
c        f2h_asympi=cf*((1+z**2)/(1-z)*(-8./3.*dlz*dlm 
c     -   - 4./3.*ddilog(1-z) + 5./3.*dlz**2 + 19./3*dlz)
c     -   - (1+z)*(2./3.*dlm**2-29./9.*dlm-4./3.*zeta2+247./54.)
c     +   + 1./3.*(1+13*z)*dlm-1./3.*(7+19*z)*dlz
c     -   - 23./18.-27./2.*z)*xqsumc

        f2ns0=(8./3.*log(xi)**2 
     +   + 2*(8./3.*dlm-58./9.)*log(xi)
     +   + 2*(-8./3.*zeta2 + 4./3.*dlm**2 - 58./9.*dlm 
     +   + 359./27.))/(1-z)

c The second-order coefficinet for the non-singlet OME  
c BMSN EPJ C1, 301 (1998)
c        f2ns0=224./27./(1-z)
c Approximate formula for the N_f-coefficient in C_2^NS,NNLO
c van Neerven-Vogt Nucl.Phys. B568, 263 (2000). 
c        f2ns0=(1.77778*dlm**2-8.5926*dlm+6.3489)/(1-z)/cf/tr
c Exact formula for the N_f-coefficient in C_2^NS,NNLO
c van Neerven-Zijlstra Phys.Lett. B272, 127 (1991).
c        f2ns0=2*(2./3.*dlm**2-29./9.*dlm-4./3.*zeta2+247./54.)/(1-z)/tr

        f2h_asympi=f2h_asympi+an**2*f2ns0*cf*tr*(xqsumc-xqsumc0)
      end if

      return
      end
C------------------
      real*8 function f2charm_ffn4(xb,q2,nq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'
      include 'PRECCOM.'

      external f2charmi
      common /forf2charm/ xb0,q2ss,rm2,qqs,rmu2,an,rq,nb0,nt0,ni0,nq0

!  Save the current scheme status
      kschemepdfs=kschemepdf
      kschemepdf=-2     !  4-flavour PDFs 

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
        f2charm_ffn4=f2c
      else
        f2charm_ffn4=0.
      end if

!  Restore the initial scheme
      kschemepdf=kschemepdfs

      return
      end
