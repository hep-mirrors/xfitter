c The codes for calculation of the strong coupling constant alpha_s
c at the scale of Q2 for different number of active flavours 
c  
c   ALPHAS_FFN   -- 3 flavours 
c   ALPHAS_FFN4  -- 4 flavours
c   ALPHAS_FFN5  -- 5 flavours
c
c The input parameters are transferred by /COMMON/ stored in CONSTCOM.
c
c The order of QCD is defined by KORDALPS  
c
c   KORDALPS=0   --   1-loop  
c   KORDALPS=1   --   2-loop  
c   KORDALPS=2   --   3-loop  
c
c The value of starting evolution scale is defined by Q20ALPHAS, 
c the value of alphas_s for 3 flavours at this scale is defined by ALPHAS0.
c
c The value of alphas_s for 4(5) flavours is matched with one for 3(4) flavours
c at the scale equal to mass of c-(b-)quark. The mass of c-quarks is defined 
c by RMASS(8) and one of b-quark is defined by RMASS(10).
c
c The service variables NFEFF, Q2REP and the constant PI are aslo stored in 
c /COMMON/ of CONSTCOM.
c
c Subroutine DSNLEQ from the CERN MATHLIB library is used to solve the 
c trancedental equations 
c
C----------------------
      DOUBLE PRECISION FUNCTION ALPHAS_ffn(Q2)

c \alpha_s in the 3-flavour scheme

      implicit double precision (a-h,o-z)

      include 'CONSTCOM.'
      include 'PRECCOM.'

      external alphastt
      real*8 tt(1),ww(20),res(1)

c  transfer the variables by the /COMMON/ of 'CONSTCOM.'

      q2s=q20alphas
      alpss=alphas0
      nfeff=3
      q2rep=q2

c  initial approximation for the solution 
      tt(1)=0.1d0
      CALL DSNLEQ(1,tt,res,alphastol,alphastol,200,0,INFO,alphastt,WW)

      alphas_ffn=tt(1)

      RETURN
      END
C----------------------
      DOUBLE PRECISION FUNCTION ALPHAS_ffn4(Q2)

c   \alpha_s in the 4-flavour scheme

      implicit double precision (a-h,o-z)

      include 'CONSTCOM.'
      include 'PRECCOM.'
      include 'PDFCOM.'

      external alphastt
      real*8 tt(1),ww(20),res(1)

      alpss=alphas_ffn(rmass(8)**2)
      q2s=rmass(8)**2
      nfeff=4
      q2rep=q2

      if (kordalps.eq.2) then 
csm NNLO matching with MSbar mass at scale mu = mbar
csm see e.g. hep-ph/0004189 eq.(23)
        if (msbarm) then 
          alpss=alpss*(1.d0+(alpss/pi)**2*(-11.d0/72.d0))
        else
csm NNLO matching with on-shell mass at scale mu = M_h
csm see e.g. hep-ph/0004189 eq.(25)
          alpss=alpss*(1.d0+(alpss/pi)**2*(7.d0/72.d0))
        end if
      end if

c  initial approximation for the solution 
      tt(1)=0.1d0
      CALL DSNLEQ(1,tt,res,alphastol,alphastol,200,0,INFO,alphastt,WW)

      alphas_ffn4=tt(1)

      RETURN
      END
C----------------------
      DOUBLE PRECISION FUNCTION ALPHAS_ffn5(Q2)

c   \alpha_s in the 5-flavour scheme

      implicit double precision (a-h,o-z)

      include 'CONSTCOM.'
      include 'PRECCOM.'
      include 'PDFCOM.'

      external alphastt
      real*8 tt(1),ww(20),res(1)

c  transfer the variables by the /COMMON/ of 'CONSTCOM.'
c  \alpha_s in the 5-flavour scheme is mathched with the 4-flavour 
c  one at the scale of b-quark mass

      alpss=alphas_ffn4(rmass(10)**2)
      q2s=rmass(10)**2
      nfeff=5
      q2rep=q2

      if (kordalps.eq.2) then 
csm NNLO matching with MSbar mass at scale mu = mbar
csm see e.g. hep-ph/0004189 eq.(23)
        if (msbarm) then
          alpss=alpss*(1.d0+(alpss/pi)**2*(-11.d0/72.d0))
        else
csm NNLO matching with on-shell mass at scale mu = M_h
csm see e.g. hep-ph/0004189 eq.(25)
          alpss=alpss*(1.d0+(alpss/pi)**2*(7.d0/72.d0))
        end if
      end if

c  initial approximation for the solution 
      tt(1)=0.1d0
      CALL DSNLEQ(1,tt,res,alphastol,alphastol,200,0,INFO,alphastt,WW)

      alphas_ffn5=tt(1)

      RETURN
      END
C--------------
      real*8 FUNCTION alphast(XX)

c  solution of the renormgroup differential equation for \alpha_s

      implicit double precision (a-h,o-z)

      include 'CONSTCOM.'

c  set up of the QCD \beta-function constants

      BET0=(11.-2./3.*nfeff)/2./pi
      BET1=(102.-38./3.*nfeff)/8./pi**2
      bet2=(2857.-5033./9.*nfeff+325./27.*nfeff**2)/64./pi**3
      BET=BET1/BET0

c  1-loop approximation  

      if (kordalps.eq.0) then 
        alphast=1./XX-1./alpss-BET0/2.*log(Q2rep/Q2s)
      end if

c  2-loop approximation 

      if (kordalps.eq.1) then 
        alphast=1./XX-1./alpss-BET*log((1./xx+BET)/(1./alpss+BET))
     - -BET0/2.*log(Q2rep/Q2s)
      end if

c  3-loop approximation 

      if (kordalps.eq.2) then
        if (bet2.ge.0.) then 
          det=sqrt(4*bet0*bet2-bet1**2)
          alphast=1./XX-1./alpss+BET*log(xx/alpss)
     -   -bet/2.*log((bet0+bet1*xx+bet2*xx**2)
     /   /(bet0+bet1*alpss+bet2*alpss**2)) 
     -   -(bet1**2-2*bet0*bet2)/bet0/det
     *   *(atan((bet1+2*bet2*xx)/det)-atan((bet1+2*bet2*alpss)/det))
     -   -BET0/2.*log(Q2rep/Q2s)
        else
          alphast=1./XX-1./alpss-BET*log((1./xx+BET)/(1./alpss+BET))
     -   -BET0/2.*log(Q2rep/Q2s)
        end if
      end if

      RETURN
      END
c-----------------------
      SUBROUTINE alphastt(N,X,F,K)  

c  interface for the subroutine DSNLEQ

      implicit double precision (a-h,o-z)

      dimension X(*),F(*) 

crp      f(k)=alphast(x(1))    
      f(k)=alphast(x)    

      RETURN 
      END 
!----------------------
      real*8 FUNCTION alphas_vfn(q2)

!  Combination of \alpha_s value at different scales with matching at 
!  the heavy quark mass values

      implicit double precision (a-h,o-z)

!  The number of loops changes at the scales equal to the quark masses 
      nfa=nfloops(q2,-1)

      if (nfa.eq.3) alphas_vfn=alphas_ffn(q2)
      if (nfa.eq.4) alphas_vfn=alphas_ffn4(q2)
      if (nfa.eq.5) alphas_vfn=alphas_ffn5(q2)

      return
      end
!----------------------
      integer FUNCTION nfloops(q2,kpdf)

!  The number of fermions in the loops of the strong coupling evolution 
!  For negative kpdf values it changes with a scale, otherwise it is 
!  defined as kpdf+3
   
      implicit double precision (a-h,o-z)

      include 'CONSTCOM.'

      if (kpdf.ge.0) then 
        nfloops=kpdf+3
      else
!  the scales below the c-quark mass: 3 fermions 
        if (q2.le.rmass(8)**2) nfloops=3

!  the scales above the c-quark mass and below the c-quark mass: 4 fermions
        if (q2.gt.rmass(8)**2.and.q2.le.rmass(10)**2) nfloops=4

!  the scales above the b-quark mass: 5 fermions
        if (q2.gt.rmass(10)**2) nfloops=5
      end if

      return
      end
!----------------------
      real*8 FUNCTION alphas_abkm(q2)

!  The same as alphas_ffn, used for the compatibility purposes

      implicit double precision (a-h,o-z)

      alphas_abkm=alphas_ffn(q2)

      return 
      end
