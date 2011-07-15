      Subroutine PDF_param_iteration(p,iflag)
C-------------------------------------------------------
C
C Created 5 June 2011. Move PDF parameterisation setting from FCN 
C
C  Input:  p(*)  -- input minuit parameters
C          iflag -- minuit flag 
C
C--------------------------------------------------------
      implicit none
      double precision p(*)
      integer iflag
      include 'pdfparam.inc'
      include 'steering.inc'
      include 'alphas.inc'
      integer i
C-------------------------------------------------------

C  25 Jan 2011: Poly params for valence:
      if (NPOLYVAL.gt.0) then
         Call StorePoly(p,iflag)
      endif

C  22 Apr 2011: CT parameterisation:
      if (IPARAM.eq.171717) then
         Call DecodeCtPara(p)
         alphas = p(14)
      endif

      if (iparam.eq.1) then     !  H1PDF2k like
         
         Bg  = p(1)
         Cg  = p(2)
         Dg  = p(3)
         
         Bu  = p(4)
         Cu  = p(5)
         Fu  = p(6)
         
         Ad  = p(7)
         Cd  = p(8)

         Cubar = p(9)
         Cdbar = p(10)
         
         Bd = Bu
         Bubar = Bu
         Bdbar = Bu
         Adbar = Ad
         aU = aD * (1.-fstrange)/(1.-fcharm)
         aUbar = aU
         
* fixed for H1param
         alphas = p(11)
         Eu = p(12)

      elseif (iparam.eq.21) then
         
         Bg  = p(1)
         Cg  = p(2)

         Bu  = p(3)
         Cu  = p(4)
         Eu = p(5)
         Fu  = p(6)
         
         Ad  = p(7)
         Cd  = p(8)
         
         Cubar = p(9)
         Cdbar = p(10)
         
         Bd = Bu
         Bubar = Bu
         Bdbar = Bu
         Adbar = Ad
         aU = aD * (1.-fstrange)/(1.-fcharm)
         aUbar = aU

* fixed for optimized H1param
         Dg = p(11)
         Alphas = p(12)
C
C  Chebyshev param. for the gluon:
C
         if (NCHEBGLU.gt.0) then
            do i=1,NCHEBGLU
               ChebPars(i) = p(20+i)
            enddo
            call ChebToPoly
         endif
         
      elseif (iparam.eq.2) then

         Bg = p(1)
         Cg = p(2)
         Dg = p(3)

         Buv = p(4)
         Bdv = Buv
         Cuv = p(5)
         Duv = p(6)
         
         Cdv = p(7)
         Ddv = p(8)

         Adbar = p(9)
         Aubar = Adbar * (1.-fstrange)/(1.-fcharm)

         Bdbar = p(10)
         Bubar = Bdbar
         
         Cdbar = p(11)
         Cubar = p(12)

* fixed for inbetween
         Alphas = p(13)
         Euv = p(14)

C
C  Chebyshev param. for the gluon:
C
         if (NCHEBGLU.gt.0) then
            do i=1,NCHEBGLU
               ChebPars(i) = p(20+i)
            enddo
            call ChebToPoly
         endif



      elseif ((iparam.eq.22).or.(iparam.eq.221).or.(iparam.eq.222)) then
         
         Bg = p(1)
         Cg = p(2)

         Buv = p(3)
               
         if (iparam.eq.221) then               
            Bdv=p(17)
         else
            Bdv = Buv
         endif
               
         Cuv = p(4)
         Duv = p(5)
         
         Cdv = p(6)
         
         Adbar = p(7)
         Aubar = Adbar * (1.-fstrange)/(1.-fcharm)
         
         Bdbar = p(8)
         Bubar = Bdbar
         
         Cdbar = p(9)
         Cubar = p(10)
         
         Euv = p(11)
               
* fixed for optimized in between
         Dg = p(12)
         Ddv = p(13)
         Alphas = p(14)
         
         Ddbar=p(15)
         Dubar=p(16)	
         
         Apg=p(18)
         Bpg=p(19)
         
         Rfudge=p(20)
         afudge=p(21)
         f2ht1=p(22)
         
         f2ht2=p(23)
         if (iparam.eq.222) then
            Cpg=25.
         endif

C
C  Chebyshev param. for the gluon:
C     
         if (NCHEBGLU.gt.0) then
            do i=1,NCHEBGLU
               ChebPars(i) = p(20+i)
            enddo
            call ChebToPoly
         endif


      elseif (iparam.eq.225) then

         Bg = p(1)
         Cg = p(2)

         Buv = p(3)
         Bdv = Buv
         Cuv = p(4)
         Duv = p(5)

         Cdv = p(6)

         Adbar = p(7)
C  fstrange -> fs0
c               Aubar = Adbar * (1.-fs0)/(1.-fcharm)
         Aubar = p(15)

         Bdbar = p(8)
c               Bubar = Bdbar
         Bubar = p(16)

         Cdbar = p(9)
         Cubar = p(10)

         Euv = p(11)

* fixed for optimized in between
         Dg = p(12)
         Ddv = p(13)
         Alphas = p(14)

C
C  Chebyshev param. for the gluon:
C
         if (NCHEBGLU.gt.0) then
            do i=1,NCHEBGLU
               ChebPars(i) = p(20+i)
            enddo
            call ChebToPoly
         endif



      elseif (iparam.eq.229) then
         
         Bg = p(1)
         Cg = p(2)
         
         Buv = p(3)
               

         Bdv=p(17)
               
         Cuv = p(4)
         Duv = p(5)

         Cdv = p(6)
               
         Adbar = p(7)
         Aubar = Adbar * (1.-fstrange)/(1.-fcharm)

         Bdbar = p(8)
         Bubar = Bdbar
               
         Cdbar = p(9)
         Cubar = p(10)
               
         Euv = p(11)
               
* fixed for optimized in between
         Dg = p(12)
         Ddv = p(13)
         Alphas = p(14)
               
         Ddbar=p(15)
         Dubar=p(16)	

         Apg=p(18)
         Bpg=p(19)
               
         Rfudge=p(20)
         afudge=p(21)
         Cpg=25.

         if (NCHEBGLU.gt.0) then
            do i=1,NCHEBGLU
               ChebPars(i) = p(20+i)
            enddo
            call ChebToPoly
         endif

      elseif (iparam.eq.3) then ! g,uval,dval,sea as in ZEUS-S 2002 fit
         
         Bg = p(1)
         Cg = p(2)
         Dg = 0.                ! Different from H1PDF2k
         
         Cuv = p(3)
         Buv = 0.5
         Duv = p(4)
         Fuv=0.
               
         Cdv = p(5)
         Bdv = 0.5
         Ddv = p(6)
         Fdv=0.
         
         Asea = p(7)
         Bsea = p(8)
         Csea = p(9)
         Dsea = p(10)
         
         Adel = p(11)
         Bdel = 0.5
         Cdel = Csea +2.

      elseif (iparam.eq.4) then ! g,uval,dval,sea as in ZEUS-JET fit
         
         Bg = p(1)
         Cg = p(2)
         Dg = p(3)         
         
         Buv = p(4)
         Cuv = p(5)
         Duv = p(6)
         
         Bdv = Buv
         Cdv = p(7)
         Ddv = p(8)
         
         Asea = p(9)
         Bsea = p(10)
         Csea = p(11)

*  dbar-ubar (not Ubar - Dbar), Adel fixed to output of ZEUS-S fit   
 
         Adel = 0.27          
         Bdel = 0.5
         Cdel = Csea +2.

* fixed for ZEUS-JETS
         Alphas = p(12)
         Euv = p(13)
C
C  Chebyshev param. for the gluon:
C
         if (NCHEBGLU.gt.0) then
            do i=1,NCHEBGLU
               ChebPars(i) = p(20+i)
            enddo
            call ChebToPoly
         endif
C
C  Chebyshev param. for the sea:
C
         if (NCHEBSea.gt.0) then
            do i=1,NCHEBSea
C  Offset is now steering parameter (default = 20, params start from 41)
               ChebParsSea(i) = p(20+IOFFSETCHEBSEA+i)
            enddo
         endif

         if (NChebGlu.gt.0 .or. NChebSea.gt.0) then
            call ChebToPoly
         endif


      elseif (iparam.eq.24) then ! g,uval,dval,sea as in ZEUS-JET fit

         Bg = p(1)
         Cg = p(2)
         Dg = p(3)         
         
         Buv = p(4)
         Cuv = p(5)
         Duv = p(6)
         
         Bdv = Buv
         Cdv = p(7)

         Euv = p(8)
               
         Asea = p(9)
         Bsea = p(10)
         Csea = p(11)

*  dbar-ubar (not Ubar - Dbar), Adel fixed to output of ZEUS-S fit   
         
         Adel = 0.27          
         Bdel = 0.5
         Cdel = Csea +2.

* fixed for ZEUS-JETS optimized
              
         Ddv = p(12)
         Alphas = p(13)
         
      endif         

      
      end


* -------------------------------------------------------
      double precision function flav_number(q)
* -------------------------------------------------------

      implicit none
      include 'thresholds.inc'

      double precision q
      
      flav_number = 3.d0
      if (q.ge.qc) flav_number = 4.d0
      if (q.ge.qb) flav_number = 5.d0

      return
      end

* -------------------------------------------------------
      double precision function gluon(x)
* -------------------------------------------------------
* x *g(x,Q2)

      implicit none
      include 'pdfparam.inc'
      include 'steering.inc'
      double precision x

C External function:
      double precision PolyParam,ctpara
C-------------------------------------------------

C    22 Apr 11, SG, Add CTEQ-like
      if (iparam.eq.171717) then
         gluon = ctpara(x,ctglue)
         return
      endif


      if (nchebglu.eq.0) then
C
C new2 jf to test term fg
C
         if ((iparam.eq.222.).or.(iparam.eq.229)) then
            gluon = ag * x**bg * (1.-x)**cg * (1. + dg * x + fg * x**3)-
     $           apg*x**bpg*(1.-x)**cpg

         else
            gluon = ag * x**bg * (1.-x)**cg * (1. + dg * x + fg * x**3)
         endif
      else
C
C SG: Use polynomial representation of cheb. 
C
         gluon = ag * PolyParam(x,nchebGlu,polyPars,chebxminlog)
         if (ichebtypeGlu.eq.0) then
C Do nothing
         else if (ichebtypeGlu.eq.1) then
            gluon = gluon * (1 - x) ! force PDFs=0 for x=1
         endif
         
      endif
      end

      Subroutine ChebToPoly()
C
C Utility to convert chebyshev to standard polynomial expansion.
C
      implicit none
      include 'pdfparam.inc'
      include 'steering.inc'
      integer i
C---------------------------------
      if (nchebGlu.gt.0) then
         do i=1,nchebmax
            polyPars(i) = 0
         enddo
         call DCHPWS(nchebGlu,ChebPars,polyPars)
      endif
      if (nchebSea.gt.0) then
         do i=1,nchebmax
            polyParsSea(i) = 0
         enddo
         call DCHPWS(nchebSea,ChebParsSea,polyParsSea)
      endif
      end      


* -------------------------------------------------------
      double precision function H1U(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x

      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.21) then
       H1U = au * x**bu * (1.-x)**cu * (1. + du*x +  eu *x**2+ fu*x**3)
      endif

      return
      end

* -------------------------------------------------------
      double precision function H1D(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x

      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.21) then
       H1D = ad * x**bd * (1.-x)**cd * (1. + dd*x)
      endif

      return
      end

* -------------------------------------------------------
      double precision function Uval(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,x23
      double precision PolyVal,ctpara
C---------------------------------------------------

C    22 Apr 11, SG, Add CTEQ-like
      if (iparam.eq.171717) then
         UVal = ctpara(x,ctuval)
         return
      endif

      if (iparam.eq.2.or.iparam.eq.3.or.iparam.eq.4
     $     .or.iparam.eq.22.or.iparam.eq.24.or.iparam.eq.225
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then

C
C 25 Jan 2011: add polynomial param 
C
         if (NPOLYVAL.eq.0) then
            Uval = aUv * x**bUv * (1.-x)**cUv
     +           * (1. + dUv*x + eUv *x**2+  fUV *x**3)
         else
C 
C PDFs are parameterised as a function of x23 = x^{2/3}
C
            x23 = x**(2.D0/3.D0)
            Uval = aUv * PolyVal(x23,NPOLYVALINT,PolyUval)
         endif
      endif

      return
      end


* -------------------------------------------------------
      double precision function Dval(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,x23
      double precision PolyVal,ctpara
C--------------------------------------------------------

C    22 Apr 11, SG, Add CTEQ-like
      if (iparam.eq.171717) then
         DVal = ctpara(x,ctdval)
         return
      endif


      if (iparam.eq.2.or.iparam.eq.3.or.iparam.eq.4
     $     .or.iparam.eq.22.or.iparam.eq.24.or.iparam.eq.225
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
C
C 25 Jan 2011: add polynomial param 
C
         if (NPOLYVAL.eq.0) then
            Dval = aDv * x**bDv * (1.-x)**cDv
     +           * (1. + dDv*x + fDv * x**3)
         else
C
C PDFs are parameterised as a function of x23 = x^{2/3}
C
            x23 = x**(2.D0/3.D0)
            Dval = aDv * PolyVal(x23,NPOLYVALINT,PolyDval)
         endif
      endif

      return
      end


      double precision function PolyVal(x,NPOLY,Poly)
C-----------------------------------------------------
C 25 Jan 2011
C
C Evaluate polynomial sum fast
C
      implicit none
      integer NPOLY
      double precision x,Poly(NPOLY)
      integer i
      double precision sum
C------------------------------------------------
      sum = 0
      do i=NPoly,1,-1
         sum = x*(sum + Poly(i))
      enddo
      PolyVal = sum
      end


* -------------------------------------------------------
      double precision function sea(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,Ubar,Dbar
C External function:
      double precision PolyParam
C--------------------------------------------------

      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.2
     $     .or.iparam.eq.21.or.iparam.eq.22.or.iparam.eq.225
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
         sea = Ubar(x) + Dbar(x)
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then 
* warning for iparam = 3 or 4, the sea is 2 * sum (ubar +dbar + sbar + cbar)

         if (nchebSea.eq.0) then
            sea = Asea * x**Bsea * (1.-x)**Csea
     +           * (1. + Dsea*x)
         else

            sea = PolyParam(x,nchebSea,polyParsSea,chebxminlog)
            if (ichebtypeSea.eq.0) then
C Do nothing
            else if (ichebtypeSea.eq.1) then
               sea = sea * (1 - x)   ! force PDFs=0 for x=1
            endif
         endif
      endif
 
      return
      end

      double precision Function PolyParam(x,ncheb,poly,xminlog)
C
C  SG: Use polynomial representation of cheb. 
C
      implicit none
      integer ncheb
      double precision x,poly(ncheb),xminlog
      double precision xx,sum
      integer i
C-------------------------------------------------
      xx = (2.D0*log(x)-xminlog)/(-xminlog)
      sum = poly(ncheb)
      do i=ncheb-1,1,-1
         sum = sum*xx + poly(i)
      enddo
      

      PolyParam = sum
      end


* -------------------------------------------------------
      double precision function dbmub(x)
* -------------------------------------------------------
* new jf , added to fit a la ZEUS, dbmub = dbar-ubar (not Dbar - Ubar)
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x
      
      if (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then 
         dbmub = Adel * x**Bdel * (1.-x)**Cdel
      endif
 
      return
      end
* -------------------------------------------------------
      double precision function qstrange (x)
* -------------------------------------------------------
* new jf , added to fit a la ZEUS, qstrange = 0.1 (i.e. fstrange *.5) * sea
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,sea,Dbar

C SG: x-dependent fs:
      double precision fs
      double precision fshermes
C----------------------------------------------------
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(x)
      endif
      
      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.2
     $     .or.iparam.eq.21.or.iparam.eq.22.or.iparam.eq.225
     $     .or.iparam.eq.171717
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then

         qstrange = fs * Dbar(x)
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then 
         qstrange = 0.5 * fs * sea(x)
      endif
 

      return
      end
* -------------------------------------------------------
      double precision function cbar(x)
* -------------------------------------------------------
* new2 jf    
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
      double precision x,pdf,q2
      dimension pdf(-6:6)

      double precision sing,flav_number,QPDFXQ,vcplus,vicplus,cplus
      integer iflag,iq0,iqc,iqfromq

      if (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
         if (x.eq.1) goto 999
         call fpdfxq(1,x,q2,pdf,0)
* charm a la ZEUS 
         if (q0.lt.qc) then
            cbar = 0.
         else

            cbar=pdf(-4)
         endif

      endif
 999  continue
      return
      end
* -------------------------------------------------------
      double precision function Ubar(x)
* -------------------------------------------------------
* new2 jf    
* corrected for iparam=2 and iparam = 3 or 4
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
      double precision x,sea,dbmub,qstrange,cbar
      double precision sing,flav_number,QPDFXQ
      integer iflag,iq0,iqb,iqc,iqfromq,jtest
      double precision ctpara
C----------------------------------------------
* new2 jf SPECIAL TEST with dubar


C    22 Apr 11, SG, Add CTEQ-like
      if (iparam.eq.171717) then
         Ubar = ctpara(x,ctubar)
         return
      endif


      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.2
     $     .or.iparam.eq.21.or.iparam.eq.22.or.iparam.eq.225
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
         Ubar = aubar * x**bubar * (1.-x)**cubar * (1. + dubar *x)
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
 
         Ubar = (0.5d0 * sea(x) - dbmub(x) - qstrange (x) + cbar(x))/2.d0

      endif
 
      return
      end

* -------------------------------------------------------
      double precision function Dbar(x)
* -------------------------------------------------------
* new2 jf
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,sea,Ubar
      double precision ctpara
C------------------------------------------------------------

C    22 Apr 11, SG, Add CTEQ-like
      if (iparam.eq.171717) then
         Dbar = ctpara(x,ctdbar)
         return
      endif


* SPECIAL TEST with ddbar      
      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.2
     $     .or.iparam.eq.21.or.iparam.eq.22.or.iparam.eq.225
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
         Dbar = adbar * x**bdbar * (1.-x)**cdbar * (1. + ddbar *x)
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
         Dbar = sea(x) * 0.5d0 - Ubar(x)
      endif

      return
      end

* -------------------------------------------------------
      double precision function singlet(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x,H1U,H1D,Ubar,Dbar,sea,Uval,Dval

      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.21) then
         singlet = H1U(x)+H1D(x)+Ubar(x)+Dbar(x)

      elseif (iparam.eq.2.or.iparam.eq.22.or.iparam.eq.225
     $        .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
         singlet = 2.d0 * sea(x) + Uval(x) + Dval(x)
* new jf
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
         singlet = sea(x) + Uval(x) + Dval(x)
      endif

      return
      end

* -------------------------------------------------------
      double precision function Uminus(x)
* -------------------------------------------------------
      implicit none
      include 'pdfparam.inc'
      include 'steering.inc'
      double precision x,H1U,Ubar,Uval

      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.21) then
       Uminus = H1U(x) - Ubar(x)
* new jf
      elseif (iparam.eq.2.or.iparam.eq.3.or.iparam.eq.4
     $      .or.iparam.eq.22.or.iparam.eq.225.or.iparam.eq.24
     $      .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
       Uminus = Uval(x)
      endif

      return
      end

* -------------------------------------------------------
      double precision function Dminus(x)
* -------------------------------------------------------
      implicit none
      include 'pdfparam.inc'
      include 'steering.inc'
      double precision x,H1D,Dbar,Dval

      if (iparam.eq.1.or.iparam.eq.11.or.iparam.eq.21) then
       Dminus = H1D(x) - Dbar(x)
* new jf
      elseif (iparam.eq.2.or.iparam.eq.3.or.iparam.eq.4
     $      .or.iparam.eq.22.or.iparam.eq.225.or.iparam.eq.24
     $      .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
       Dminus = Dval(x)
      endif

      return
      end

* -------------------------------------------------------
      double precision function Deltaud(x)
* -------------------------------------------------------
      implicit none
      include 'pdfparam.inc'
      include 'steering.inc'
      double precision x,H1U,Ubar,H1D,Dbar
      double precision Uval,Dval,qstrange,cbar
C SG: x-dependent fs:
      double precision fs
      double precision fshermes
C----------------------------------------------------
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(x)
      endif

      if (iparam.eq.1.or.iparam.eq.21) then
        Deltaud = H1U(x) - fcharm * Ubar(x)
     +          + Ubar(x) * (1. - fcharm)
     +          - (H1D(x) - fs * Dbar(x))
     +          - (Dbar(x) * (1. - fs))
      elseif (iparam.eq.2.or.iparam.eq.22.or.iparam.eq.225
     $       .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
* doubious factor two removed -- REVERT CHANGE (Sasha Glazov)
         Deltaud =  Uval(x) + 2.*(1.-fcharm)*Ubar(x)
     +        - (Dval(x) + 2.*(1.-fs)*Dbar(x))
CSG        Deltaud =  Uval(x) + (1.-fcharm)*Ubar(x)
CSG     +          - (Dval(x) + (1.-fs)*Dbar(x))
* new jf
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
*  doubious factor two removed  -- revert change
         Deltaud =  Uval(x) + 2.*(Ubar(x)-cbar(x))
     +        - (Dval(x) + 2.*(Dbar(x)-qstrange(x)))
CSG         Deltaud =  Uval(x) + (Ubar(x)-cbar(x))
CSG     +          - (Dval(x) + (Dbar(x)-qstrange(x)))
      elseif (iparam.eq.11) then
        Deltaud = H1U(x) + Ubar(x)
     +          - (H1D(x) - fs * Dbar(x))
     +          - (Dbar(x) * (1. - fs))
      endif

      return
      end

* -------------------------------------------------------
      double precision function Splus(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
      double precision x,Dbar,singlet,flav_number,qstrange
C SG: x-dependent fs:
      double precision fs
      double precision fshermes
C----------------------------------------------------
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(x)
      endif


      Splus = 2.d0 * fs * Dbar(x) 
     +      - singlet(x)/flav_number(q0)
*new jf
      If (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
         Splus = 2.d0 * qstrange(x) 
     +      - singlet(x)/flav_number(q0)
      endif

      return
      end

* -------------------------------------------------------
      double precision function Uplus(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
      double precision x,Ubar,H1U,singlet,flav_number
      double precision Uval,cbar


* to be modified for iparam = 2 or 1 when q0.lt.qc

      if (iparam.eq.1.or.iparam.eq.21) then
        Uplus = Ubar(x)*(1.-fcharm)
     +          + H1U(x) - fcharm*Ubar(x) 
     +          - singlet(x)/flav_number(q0)

* new jf
      elseif (iparam.eq.2.or.iparam.eq.22.or.iparam.eq.225
     $       .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
        Uplus = Uval(x) + 2.*(1.-fcharm)*Ubar(x) 
     +        - singlet(x)/flav_number(q0)
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
         Uplus = Uval(x) + 2.d0 * (Ubar(x)-cbar(x))
     +    - singlet(x)/flav_number(q0)
      elseif (iparam.eq.11) then
        Uplus = Ubar(x)
     +          + H1U(x) 
     +          - singlet(x)/flav_number(q0)
      endif

      return
      end

* -------------------------------------------------------
      double precision function Dplus(x)
* -------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
      double precision x,Dbar,H1D,singlet,flav_number,Dval,qstrange
C SG: x-dependent fs:
      double precision fs
      double precision fshermes
C----------------------------------------------------
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(x)
      endif
*mis print corrected 30/04/2008 (jf)
      if (iparam.eq.1.or.iparam.eq.21) then
        Dplus = Dbar(x)*(1.-fs)
     +          + H1D(x) - fs*Dbar(x) 
     +          - singlet(x)/flav_number(q0)
      elseif (iparam.eq.2.or.iparam.eq.22.or.iparam.eq.225
     $       .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
        Dplus = Dval(x) + 2.*(1.-fs)*Dbar(x) 
     +        - singlet(x)/flav_number(q0)
* new jf
      elseif (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
        Dplus = Dval(x) + 2.*(Dbar(x)-qstrange(x)) 
     +        - singlet(x)/flav_number(q0)
      endif

      return
      end


* -------------------------------------------------------
      double precision function Cplus(x)
* -------------------------------------------------------
* given at starting scale if q0 > qc (allows to write
* that c+cbar is a fraction fcharm of Ubar for charm H1 style);
* else defined at qc (new2)

      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
      double precision x,singlet,flav_number,Ubar,factor
      integer iflag,jtest,iqfrmq,iq0,iqc,iqb

* Charm H1 style 

      if (iparam.eq.1.or.iparam.eq.2.or.iparam.eq.21.or.iparam.eq.22.or.iparam.eq.225
     $     .or.iparam.eq.221.or.iparam.eq.222.or.iparam.eq.229) then
         if (q0.ge.qc) then
            Cplus = 2.d0 * fcharm * Ubar(x) 
     +       - singlet(x)/flav_number(q0)
         else

            Cplus = - singlet(x)/flav_number(qc)
           
         endif
      endif

* Charm ZEUS style

      if (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.11.or.iparam.eq.24) then

         Cplus = - singlet(x)/flav_number(qc)
         if(q0.ge.qc) then
            IQ0 = iqfrmq(q0)
            IQC = iqfrmq(qc)
            IQB = iqfrmq(qb)
            factor = -1. /flav_number(qc)
         endif
      endif

      return
      end

* -------------------------------------------------------
      double precision function Bplus(x)
* -------------------------------------------------------
* given at scale = b threshold
      implicit none

      include 'thresholds.inc'
      double precision x,singlet,flav_number
      integer iflag



      Bplus = - singlet(x)/flav_number(qb)
      return
      end

* -------------------------------------------------------
      double precision function Deltacs(x)
* -------------------------------------------------------
      implicit none
      include 'pdfparam.inc'
      include 'steering.inc'
      double precision x,Ubar,Dbar,qstrange,cbar
C SG: x-dependent fs:
      double precision fs
      double precision fshermes
C----------------------------------------------------
      if (ifsttype.eq.0) then
         fs = fstrange
      else
         fs = fshermes(x)
      endif

      Deltacs = 2.*fcharm*Ubar(x) 
     +        - 2.*fs*Dbar(x)
*new jf
      If (iparam.eq.3.or.iparam.eq.4.or.iparam.eq.24) then
         Deltacs = 2.*cbar(x) 
     +        - 2.*qstrange(x)
      endif

      return
      end

      
      double precision function fshermes(x)
C
C X-dependent strange fraction, inspired by HERMES data
C Created 31 Oct 2009 by SG.
C
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision x
      double precision hermes_xcent,hermes_xrise
      parameter (hermes_xcent= 0.07)
      parameter (hermes_xrise=20.0)
      logical lfirstlocal
      data lfirstlocal/.true./
C---------------------------------------------------
      if (lfirstlocal) then
         lfirstlocal = .false.
         print *,'----------------------------------------------------'
         print *,' Called FSHERMES. Hermes-inspired strange density '
         print *,' Amp,Xmean,Xslope=',fstrange,hermes_xcent
     $        ,hermes_xrise
         print *,'----------------------------------------------------'
      endif

      if (x.lt.1.0D-8) then
         fshermes = fstrange
      else
         fshermes = fstrange*( 
     $        0.5D0*(1.D0+tanh(-(x-hermes_xcent)*hermes_xrise)))
      endif
C      print *,'DEBUG:',x,fshermes,fstrange
C---------------------------------------------------
      end

      subroutine StorePoly(p,iflag)
C
C Created 27 Jan 2011 by SG. Transfer parameters from MINUIT array p to
C internal arrays PolyUval and PolyDval
C
      implicit none
      include 'steering.inc'
      include 'pdfparam.inc'
      double precision p(*)
      integer iflag,i,nn
      double precision ptmp(100)

C-----------------------------------      
      if (iflag.eq.1) then
         print *,
     $        'INFO: USING POLY parameterisation for valence quarks'
         print '(''Require PDFs at x=1 to vanish as (1-x)**'',I1)'
     $        ,IZPOPOLY
         if (IPOLYSQR.eq.1) then
            print *,
     $      'Square valence parameterisation, enforce positivity'
         else if (IPOLYSQR.eq.0) then
         else
            print *,'Invalid IPOLYSQR=',IPOLYSQR
            print *,'Can be 1 or 0, stop'
            stop
         endif
         
      endif
      
      if (IPOLYSQR.eq.0) then
         Call DecodePoly(p(21),PolyUVal,NPOLYVAL,IZPOPOLY)
         Call DecodePoly(p(31),PolyDVal,NPOLYVAL,IZPOPOLY)
         NPOLYVALINT = NPOLYVAL+IZPOPOLY
      else if  (IPOLYSQR.eq.1) then
         Call DecodePoly(p(21),Ptmp,NPOLYVAL,IZPOPOLY)
         Call SquarePoly(NPOLYVAL+IZPOPOLY,PTmp,NPOLYVALINT,PolyUVal)
         Call DecodePoly(p(31),Ptmp,NPOLYVAL,IZPOPOLY)
         Call SquarePoly(NPOLYVAL+IZPOPOLY,PTmp,NPOLYVALINT,PolyDVal)
      else
         print *,'Invalid IPOLYSQR=',IPOLYSQR
         print *,'Can be 1 or 0, stop'
         stop
      endif

      end      


      subroutine DecodePoly(pars,poly,np,iz)
C
C Created 29 Jan 2011. Decode input minuit parameters to internal PolyVal
C arrays for different order of (1-x) 
C     
      implicit none
      integer np,iz
      double precision pars(*),poly(*)
      integer i
C----------------------------------------------
      if (iz.eq.1) then
* \times (1-x)
         Poly(1) = pars(1)
         do i=2,np
            Poly(i) = pars(i)-pars(i-1)
         enddo
         Poly(np+1) = -pars(np)
      elseif (iz.eq.2) then
* \times (1-x)^2
         Poly(1) = pars(1)
         Poly(2) = pars(2)-2*pars(1)
         do i=3,np
            Poly(i) = pars(i)-2*pars(i-1)+pars(i-2)
         enddo
         Poly(np+1) = -2*pars(np)+pars(np-1)
         Poly(np+2) = pars(np)
      endif
      
      end

      subroutine squarepoly(Npar,PolyIn,Npar2,PolyOut)
C-----------------------------------------
C
C     1 Feb 2011 by SG. Square a polynom. 
C
C-----------------------------------------
      implicit none
      integer Npar, Npar2
      double precision PolyIn(Npar), PolyOut(*)
      integer j,k,i
C--------------------------------------------------
      NPar2 = NPar*2-1

      do k=1,NPar2
         PolyOut(k) = 0.
         do j=1,k
            i = k-j+1
            if (i.le.NPar.and.j.le.NPar) then
               PolyOut(k) = PolyOut(k) + PolyIn(i)*PolyIn(j)
            endif
         enddo
      enddo

      end

      subroutine DecodeCtPara(pars)
C-------------------------------------------------------
C Created 22 Apr 11 by SG. Decode minuit input for CTEQ-like param.
C   pars(20-25)  - Uv
C   pars(30-35)  - Dv
C   pars(40-45)  - Ubar
C   pars(50-55)  - Dbar
C   pars(60-65)  - gluon
C------------------------------------------------------
      implicit none 
      include 'pdfparam.inc'
      include 'steering.inc'
      double precision pars(*)
      integer i
      logical lfirstt
      data lfirstt /.true./

C---------------------------------------------------------
      if (lfirstt) then
         lfirstt = .false.
         print *,'DecodeCtPara INFO: First time call'        
      endif
C simple copy first:
      do i=0,5
         ctuval(i) = pars(20+i)
         ctdval(i) = pars(30+i)
         ctubar(i) = pars(40+i)
         ctdbar(i) = pars(50+i)
         ctglue(i) = pars(60+i)
      enddo

C Extra constrains:
      if (pars(40).eq.0) then
         ctubar(0) = ctdbar(0) * (1.-strange_frac) ! normalization ubar = dbar 
         ctubar(1) = ctdbar(1)  ! Bubar = Bdbar
      endif

C Impose Buv = Bdv if parameter for Buv = 0.
      if (pars(20+1).eq.0) then
         ctuval(1) = ctdval(1)  ! Buv = Bdv
      endif
C (other constraints from sum-rules)


C---------------------------------------------------------
      end

      double precision function ctpara(x,a)
C----------------------------------------------------
C
C cteq-like parameterisation: 
C  UF = a0*E**(a3*x)*(1 - x)**a2*x**(a1 + n)*(1 + E**a4*x + E**a5*x**2)
C
C-----------------------------------------------------
      implicit none
      double precision x,a(0:5)
      double precision UF
      UF = a(0)*exp(a(3)*x)*(1 - x)**a(2)*x**(a(1))*(1 + exp(a(4))*x 
     $     + exp(a(5))*x**2)
      
      ctpara = UF

      end
