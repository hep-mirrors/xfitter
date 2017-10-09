C------------------
      real*8 function f3qcd(nb,nt,ni,xb,q2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      external f3nloifg

      include 'PDFCOM.'
      include 'PRECCOM.'
      include 'CONSTCOM.'

      common /forconvol/ xb0,q2ss,an,nb0,nt0,ni0

!  Set the number of fermions in the final state, nfc, and the 
!  number of fermions in the loops, nfe
      nfc=max(kschemepdf+3,3)
      nfe=nfloops(q2,kschemepdf)

      f30=f3qpm(nb,nt,ni,xb,q2)
      f3qcd=F30

      if (kordf3.eq.0) return 

!  Take coupling constant from the grid filled in the PDF evolution  
      an=xqg(0,0.1d0,q2,kschemepdf)/4./pi 

      xb0=xb
      q2ss=q2
      nb0=nb
      nt0=nt
      ni0=ni
      rfac=1.
      if (kordf3.ge.2) then
        alr=-log(rscale)
        bet0=11-2./3.*nfe
        rfac=1-an*bet0*alr
      end if

      CALL GAUSS1(f3nloifg,xb,1D0,nf3qcd,df,EPS)

      f3qcd=f3qcd+df+an*f30*c3ns_1_0_local(xb)*rfac

      if (kordf3.ge.2) then 
        f3qcd=f3qcd+an**2*f30
     *   *(c3nsm_2_0_nf0_local(xb) + nfc*c3nsm_2_0_nf1_local(xb))

      end if

      return
      end
C------------------
      real*8 function f3nloifg(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'

      common /forconvol/ xb0,q2ss,an,nb0,nt0,ni0

      z=xb0/y

      rfac=1.
      if (kordf3.ge.2) then 
        alr=-log(rscale)
        bet0=11-2./3.*nfe
        rfac=1-an*bet0*alr
      end if

      f3ns0=an*c3ns_1_0_singular(y)*rfac
      f3ns=an*c3ns_1_0(y)*f3qpm(nb0,nt0,ni0,z,q2ss)/y*rfac

      if (kordf3.ge.2) then 
        f3ns0=f3ns0+(c3nsm_2_0_nf0_singular(y)
     +          + nfc*c3nsm_2_0_nf1_singular(y))*an**2

        f3ns=f3ns+(c3nsm_2_0_nf0(y) + nfc*c3nsm_2_0_nf1(y))
     *  *an**2*f3qpm(nb0,nt0,ni0,z,q2ss)/y

        if (ni0.eq.24) then                                 !CC 
          cdel= c3nspm_2_0(y)/2.
          if (nb0.eq.6) then                                !neutrino beam 
            f3ns=f3ns+an**2*cdel*(f3qpm(6,nt0,ni0,z,q2ss)
     -                           -f3qpm(7,nt0,ni0,z,q2ss))/y
          end if
          if (nb0.eq.7) then                                !antineutrino beam
            f3ns=f3ns+an**2*cdel*(f3qpm(7,nt0,ni0,z,q2ss)
     -                         -f3qpm(6,nt0,ni0,z,q2ss))/y
          end if
        end if
      end if

      f3ns0=f3ns0*(f3qpm(nb0,nt0,ni0,z,q2ss)/y
     -            -f3qpm(nb0,nt0,ni0,xb0,q2ss))

      f3nloifg=f3ns0+f3ns

      return 
      end
C------------------
      real*8 function flqcd(nb,nt,ni,xb,q2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      external flqcdi

      include 'PDFCOM.'
      include 'CONSTCOM.'
      include 'PRECCOM.'

      common /forconvol/ xb0,q2ss,an,nb0,nt0,ni0

!  Set the number of fermions in the final state, nfc, and the 
!  number of fermions in the loops, nfe
      nfc=max(kschemepdf+3,3)
      nfe=nfloops(q2,kschemepdf)
     
      flqcd=0.

      if (kordfl.eq.0) return   

!  Take coupling constant from the grid filled in the PDF evolution  
      an=xqg(0,0.1d0,q2,kschemepdf)/4./pi 

      xb0=xb
      q2ss=q2
      nb0=nb
      nt0=nt
      ni0=ni

      alr=-log(rscale)
      bet0=11-2./3.*nfe

      CALL GAUSS1(flqcdi,dlog(xb),0D0,nflqcd,fl,EPS)
      flqcd=fl

      if (kordfl.ge.2) then 
        dd=an**2*clns_2_0_nf0_local(xb)*f2qpm(nb,nt,ni,xb,q2ss) 
        flqcd=flqcd*(1-an*bet0*alr)+dd
      end if
      if (kordfl.ge.3) then 
        dd=an**3*(clns_3_0_nf0_local(xb) + nfc*clns_3_0_nf1_local(xb))
     *    *f2qpm(nb,nt,ni,xb,q2ss)
        flqcd=flqcd+dd
      end if

      return
      end
C------------------
      real*8 function flqcdi(t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'CONSTCOM.'
      include 'PDFCOM.'
      include 'PRECCOM.'
      real*8 cc(3)

      common /forconvol/ xb0,q2ss,an,nb0,nt0,ni0

      y=dexp(t)
      z=xb0/y

      flns=an*clns_1_0(y)
      flg=an*clg_1_0(y)
      flps=0.

      if (kordfl.ge.2) then   
        flns=flns+an**2*(clns_2_0_nf0(y) + nfc*clns_2_0_nf1(y))
        flg=flg+an**2*clg_2_0(y)
        flps=an**2*clps_2_0(y)
!  C-odd term for the neutrino-induced case
        if (ni0.eq.24) then                          ! CC
          if (nb0.eq.6) then                         !neutrino beam 
            flns=flns+clnspm_2_0(y)/2.
     *    *(f2qpm(7,nt0,ni0,z,q2ss)-f2qpm(6,nt0,ni0,z,q2ss))*an**2
          end if
          if (nb0.eq.7) then                         !antineutrino beam
            flns=flns+clnspm_2_0(y)/2.
     *    *(f2qpm(6,nt0,ni0,z,q2ss)-f2qpm(7,nt0,ni0,z,q2ss))*an**2
          end if
        end if
      end if

      if (kordfl.ge.3) then   
        if (lpcint) then 
          call cintx(1,2,2,nfc,y,cc)
          flns=flns+an**3*cc(1)
          flg=flg+an**3*cc(2)
          flps=flps+an**3*cc(3)
        else
          flns=flns+an**3*clns_3_0(y,nfc)
          flg=flg+an**3*clg_3_0(y,nfc)
          flps=flps+an**3*clps_3_0(y,nfc)
        end if
      end if

      flns=flns*f2qpm(nb0,nt0,ni0,z,q2ss) 
      flg=flg*xqg(1,z,q2ss,kschemepdf)

      if (ni0.eq.22) then
        flg=qsum(nfc)*flg   
        flps=qsum(nfc)*flps*xpsqpm(z,q2ss,kschemepdf) 
      end if

      if (ni0.eq.23) then
        flg=vaq2sum(nfc)*flg   
        flps=vaq2sum(nfc)*flps*xpsqpm(z,q2ss,kschemepdf) 
      end if

      if (ni0.eq.25) then
        flg=vaqsum(nfc)*flg   
        flps=vaqsum(nfc)*flps*xpsqpm(z,q2ss,kschemepdf) 
      end if

      if (ni0.eq.24) then
        flps=flps*xpsqpm2(z,q2ss,kschemepdf)
      end if

      flqcdi=(flns+flg+flps)*y

      return
      end
C------------------
      real*8 function f2qcd(nb,nt,ni,xb,q2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      external f2nloifg,f2nloifglog,f2nloifglog1

      include 'PDFCOM.'
      include 'PRECCOM.'
      include 'CONSTCOM.'

      common /forconvol/ xb0,q2ss,an,nb0,nt0,ni0

!  Set the number of fermions in the final state, nfc, and the 
!  number of fermions in the loops, nfe
      nfc=max(kschemepdf+3,3)
      nfe=nfloops(q2,kschemepdf)
     
      f20=f2qpm(nb,nt,ni,xb,q2)
      f2qcd=f20

      if (kordf2.eq.0) return 
!  Take coupling constant from the grid filled in the PDF evolution  
      an=xqg(0,0.1d0,q2,kschemepdf)/4./pi 

      q2ss=q2
      xb0=xb
      nb0=nb
      nt0=nt
      ni0=ni
      rfac=1.

      if (kordf2.ge.2) then
        alr=-log(rscale)
        bet0=11-2./3.*nfe
        rfac=1-an*bet0*alr
      end if

      if (xb.ge.0.1) then  
        CALL GAUSS1(f2nloifglog1
     ,            ,log(1d-8),log(1.-xb),nf2qcd2,df,EPS)
      else 
        CALL GAUSS1(f2nloifglog1
     ,            ,log(1d-8),log(0.9d0),nf2qcd2,df1,EPS1)
        CALL GAUSS1(f2nloifglog,log(xb)
     ,            ,log(0.1d0),nf2qcd1,df2,EPS2)
        df=df1+df2
      end if 

      f2qcd=f2qcd+df+an*f20*c2ns_1_0_local(xb)*rfac

      if (kordf2.ge.2) then 
        f2qcd=f2qcd+an**2*(c2nsp_2_0_nf0_local(xb)
     +         + nfc*c2nsp_2_0_nf1_local(xb))*f20
        if (ni.eq.22) then 
          f2qcd=f2qcd+an**2*c2g_2_0_local(xb)
     *                *xqg(1,xb,q2,kschemepdf)*qsum(nfc)  
        end if
        if (ni.eq.23) then 
          f2qcd=f2qcd+an**2*c2g_2_0_local(xb)
     *                *xqg(1,xb,q2,kschemepdf)*vaq2sum(nfc)  
        end if
        if (ni.eq.25) then 
          f2qcd=f2qcd+an**2*c2g_2_0_local(xb)
     *                *xqg(1,xb,q2,kschemepdf)*vaqsum(nfc)  
        end if
      end if

      return
      end
C------------------
      real*8 function f2nloifg(y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'PDFCOM.'
      include 'CONSTCOM.'

      common /forconvol/ xb0,q2ss,an,nb0,nt0,ni0

      z=xb0/y

      rfac=1.
      if (kordf2.ge.2) then 
        alr=-log(rscale)
        bet0=11-2./3.*nfe
        rfac=1-an*bet0*alr
      end if

      f2ns=c2ns_1_0(y)*f2qpm(nb0,nt0,ni0,z,q2ss)*an*rfac
      f2ns0=c2ns_1_0_singular(y)*an*rfac
      f2g=c2g_1_0(y)*an*rfac
      f2ps=0.

      if (kordf2.ge.2) then 
        f2ns0=f2ns0+(c2nsp_2_0_nf0_singular(y)
     +         + nfc*c2nsp_2_0_nf1_singular(y))*an**2
        f2ns=f2ns+(c2nsp_2_0_nf0(y) + nfc*c2nsp_2_0_nf1(y))
     *     *f2qpm(nb0,nt0,ni0,z,q2ss)*an**2
        f2g=f2g+c2g_2_0(y)*an**2
        f2ps=c2ps_2_0(y)*an**2
!  Axial term term for Z-exchange neutral current
        if (ni0.eq.23) then                          
          f2ps=f2ps + c2ps_2_axial_0(y)*an**2
        end if
!  C-odd term for the neutrino-induced case
        if (ni0.eq.24) then                          ! CC
          if (nb0.eq.6) then                         ! neutrino beam 
            f2ns=f2ns+c2nspm_2_0(y)/2.
     *    *(f2qpm(7,nt0,ni0,z,q2ss)-f2qpm(6,nt0,ni0,z,q2ss))*an**2
          end if
          if (nb0.eq.7) then                         ! antineutrino beam
            f2ns=f2ns+c2nspm_2_0(y)/2.
     *    *(f2qpm(6,nt0,ni0,z,q2ss)-f2qpm(6,nt0,ni0,z,q2ss))*an**2
          end if
        end if
      end if

      f2ns0=f2ns0*(f2qpm(nb0,nt0,ni0,z,q2ss)
     -            -f2qpm(nb0,nt0,ni0,xb0,q2ss))
      f2g=f2g*xqg(1,z,q2ss,kschemepdf)

      if (ni0.eq.22) then
        f2g=qsum(nfc)*f2g  
        f2ps=qsum(nfc)*f2ps*xpsqpm(z,q2ss,kschemepdf)
      end if
      if (ni0.eq.23) then
        f2g=vaq2sum(nfc)*f2g  
        f2ps=vaq2sum(nfc)*f2ps*xpsqpm(z,q2ss,kschemepdf)
      end if
      if (ni0.eq.25) then
        f2g=vaqsum(nfc)*f2g  
        f2ps=vaqsum(nfc)*f2ps*xpsqpm(z,q2ss,kschemepdf)
      end if
      if (ni0.eq.24) then
        f2ps=f2ps*xpsqpm2(z,q2ss,kschemepdf)
      end if

      f2nloifg=(f2g+f2ns0+f2ps+f2ns)

      return 
      end
C------------------
      real*8 function f2nloifglog(t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'PDFCOM.'
      include 'CONSTCOM.'

      common /forconvol/ xb0,q2ss,an,nb0,nt0,ni0

      y=exp(t)

      f2nloifglog=f2nloifg(y)*y

      return 
      end
C------------------
      real*8 function f2nloifglog1(t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'PDFCOM.'
      include 'CONSTCOM.'

      common /forconvol/ xb0,q2ss,an,nb0,nt0,ni0

      y=1.-exp(t)

      f2nloifglog1=f2nloifg(y)*(1-y)

      return 
      end
