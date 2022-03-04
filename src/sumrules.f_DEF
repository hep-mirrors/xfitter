
* --------------------------------------------------------------
      subroutine SumRules(kflag)
* --------------------------------------------------------------

      implicit none

#include "steering.inc"
#include "pdfparam.inc"
#include "for_debug.inc"

      integer kflag
      double precision t1,t2,t3,t4,term
      double precision t1mt3,t1mt4
      double precision tu,tg,td,tubar,tdbar
      double precision tsbar, tcbar,tsmalldb, tsmallub
      double precision CalcIntegral,CalcIntegralCheb

      double precision CalcIntXpdf,CalcIntXpdfFixN,CalcIntPdf

      double precision SSDINT
      double precision ToInteg
      external ToInteg
      double precision zero,tgMRST

      double precision btmp,c1,c2
      common/For_Integ/btmp,c1,c2
*new joel feltesse
      double precision tuv,tdv,tub,tdb,tsea,tdel
      double precision polyvalint,polyvalint0
      double precision para,x
      double precision ubar,dbar,uval,dval,gluon
      integer i

      double precision fs
      double precision fshermes
      double precision tstr,tNoGlue,tPho
*add for mixed CTEQHERA
      double precision SumRuleCTEQ, SumRuleCTEQhera

C-----------------------------------------
      kflag=0
      zero = 1d-10

C=========================================================
C Nothing to do for LHAPDF or Diffractive:
C
      if (PDF_DECOMPOSITION.eq.'LHAPDF' 
     $     .or. PDF_DECOMPOSITION.eq.'Diffractive'
     $     .or. PDF_DECOMPOSITION.eq.'QCDNUM_GRID' ) then
         Return
      endif

C==========================================================
C CTEQ-like parameterisation:
C
      if (PDFStyle.eq.'CTEQ') then
         Call SumRulesCTeq
         Return
      endif

      if (PDFStyle.eq.'CTEQHERA') then
         Call SumRulesCTEQHera
         Return
      endif

C    22 Sep 11, VR, Add AS parametrisation
      if ((PDFStyle.eq.'AS').or.(PDFStyle.eq.'BiLog')) then
         Call SumRulesAS
         return
      endif



C==========================================================
C Standard parameterisation. 
C

C--------------
C Valence:
      if (Index(PDF_DECOMPOSITION,'Dv_Uv').gt.0) then

C**********************************************************
C*     -- sum rule : D - Dbar = 1   :  gives ADval
C*

         if (pardval(1).eq.0) then
            pardval(1) = 1.0d0/CalcIntPdf(pardval)
         else
            dv_sum = pardval(1)*CalcIntPdf(pardval)
         endif
            
C**********************************************************
C*     -- sum rule : U - Ubar = 2   :  gives AUval
C*
         if (paruval(1).eq.0) then
            paruval(1) = 2.0D0/CalcIntPdf(paruval)
         else
            uv_sum = paruval(1)*CalcIntPdf(paruval)/2.
         endif
            
C Also integrate momenta, for momentum sum rule:
         tUv = paruval(1)*CalcIntXpdf(paruval)
         tDv = pardval(1)*CalcIntXpdf(pardval) 
cv         print*,'sumrules......', tuv, tdv

      else

         print *,'Un-implemented valence decomposition '//PDF_DECOMPOSITION
         print *,'Stop in sumrules'
         call HF_STOP
      endif


C**********************************************************
C*     -- sum rule : x ( gluon + Sigma) = 1  :  gives Ag
C*

C----------------
C Gluon:

C Check chebyshev and flexible gluon:
      if (nchebglu.eq.0) then
         if (FlexibleGluon) then
            tg = CalcIntXpdfFixN(parglue,6)             
            tgMRST=CalcIntegral(parglue(8),parglue(9))
         else
            tg = CalcIntXpdf(parglue)
            tgMRST=0.0 
         endif
      else
         tg = CalcIntegralCheb(nchebglu,
     $        polyPars,chebxminlog, ichebtypeGlu)
      endif

C----------------
C Sea:
      if (Index(PDF_DECOMPOSITION,'Dbar_Ubar').gt.0) then

         tUb = parubar(1)*CalcIntXpdf(parubar)
         tDb = pardbar(1)*CalcIntXpdf(pardbar)
         if (iTheory.eq.11.or.iTheory.eq.35) then
            tPho = parphoton(1)*CalcIntXpdf(parphoton)
         else
            tPho = 0 !> activate it only when QED is needed
         endif

         if (Index(PDF_DECOMPOSITION,'Str').gt.0) then
            tStr = parstr(1)*CalcIntXpdf(parstr)
         else
            tStr = 0   !> Strange already included in Dbar
         endif

C Total sea integral:
         tsea = 2.0d0 * (tUb + tDb + tStr)

      elseif (Index(PDF_DECOMPOSITION,'Sea').gt.0) then

         if (nchebsea.eq.0) then                   
            tsea = CalcIntXpdf(parsea)
         else
            tsea = CalcIntegralCheb(nchebsea,
     $           polyParsSea,chebxminlog,ichebtypeSea )
            Parsea(1) = 1.d0
         endif

      else
         print *,'Un-implemented sea decomposition '//PDF_DECOMPOSITION
         print *,'Stop in sumrules'
         call HF_STOP
      endif

C
C 1 - (valence + sea momentum):
C      
      tNoGlue = 1.D0 - ( tUv + tDv + tSea + tPho)

C*******************************************************************

      if (tg.le.0) then
         tg=0.001
      endif

C Calculate gluon normalisation, taking into account flexible piece:
      if (parglue(1).eq.0) then  !> Impose sum rule
         parglue(1)=(tNoGlue+parglue(7)*tgMRST)/tg
      else
         p_sum = parglue(1)*tg + tUv + tDv + tSea 
     $        + tPho - parglue(7)*tgMRST
      endif

C*******************************************************************

C     propagate the normalizations and other parameters to
C     standard parametrisation

cVR      print*,'........................................tphoton', tPho

      if (NCHEBGLU.eq.0) then         
      if (lprint) then
         print '(''uv:'',11F10.4)',(paruval(i),i=1,10)
         print '(''dv:'',11F10.4)',(pardval(i),i=1,10)
         print '(''Ub:'',11F10.4)',(parubar(i),i=1,10)
         print '(''Db:'',11F10.4)',(pardbar(i),i=1,10)
         print '(''GL:'',11F10.4)',(parglue(i),i=1,10)
         print '(''ST:'',11F10.4)',(parstr(i),i=1,10)
         if (iTheory.eq.11.or.iTheory.eq.35) then
            print '(''PH:'',11F10.4)',(parphoton(i),i=1,10)
         endif
         if (uv_sum.ne.0.or. dv_sum.ne.0 .or. p_sum.ne.0) then
            print '(''Sum rules, uv, dv, p:'',3F10.4)'
     $           ,uv_sum, dv_sum, p_sum
         endif
      endif
      endif
      
 999  continue
      return
      end


* --------------------------------------------------------------
      double precision function ToInteg(x)
* --------------------------------------------------------------

      implicit none

      double precision x
      double precision b,c1,c2
      common/For_Integ/b,c1,c2
      double precision xloc

      xloc = x

      ToInteg = xloc**(b-1) * ( (1.-xloc)**c1 - (1.-xloc)**c2)

      return
      end


* --------------------------------------------------------------
      double precision function CalcIntXpdf(pdfpars)
C---------------------------------------------------------------
C Calculated  \int xpdf(x) dx using the standard PDF
C parameterisation
C---------------------------------------------------------------
      implicit none
      double precision pdfpars(10)
      integer i
      double precision sum
      
      double precision CalcIntegral
C---------------------------------------------------------------
      sum =  CalcIntegral(pdfpars(2),pdfpars(3)) 
      ! WS 2015-10-08
      ! do i=1,7
      do i=1,3
         if ( pdfpars(3+i).ne.0 ) then
            sum = sum + pdfpars(3+i) 
     $           * CalcIntegral(pdfpars(2)+i,pdfpars(3))
         endif 
      enddo

C Also espsilon times sqrt x:
      if ( pdfpars(10).ne.0) then
         sum = sum + pdfpars(10)*CalcIntegral(pdfpars(2)+0.5,pdfpars(3)) 
      endif

      CalcIntXpdf = sum
C---------------------------------------------------------------
      end


* --------------------------------------------------------------
      double precision function CalcIntPdf(pdfpars)
C---------------------------------------------------------------
C Calculated  \int pdf(x) dx using the standard PDF
C parameterisation
C---------------------------------------------------------------
      implicit none
      double precision pdfpars(10)
      integer i
      double precision sum
      
      double precision CalcIntegral
C---------------------------------------------------------------
      sum =  CalcIntegral(pdfpars(2)-1,pdfpars(3)) 
      ! WS 2015-10-08
      ! do i=1,7
      do i=1,3
         if ( pdfpars(3+i).ne.0 ) then
            sum = sum + pdfpars(3+i) 
     $           * CalcIntegral(pdfpars(2)+i-1,pdfpars(3))
         endif 
      enddo

C Also espsilon times sqrt x:
      if ( pdfpars(10).ne.0) then
         sum = sum + pdfpars(10)*CalcIntegral(pdfpars(2)-0.5,pdfpars(3)) 
      endif

      CalcIntPdf = sum
C---------------------------------------------------------------
      end

* --------------------------------------------------------------
      double precision function CalcIntXpdfFixN(pdfpars,n)
C---------------------------------------------------------------
C Calculated  \int xpdf(x) dx using the standard PDF
C parameterisation. Sum up to N-th term (max N=6)
C---------------------------------------------------------------
      implicit none
      double precision pdfpars(10)
      integer N
C---
      integer i
      double precision sum
      
      double precision CalcIntegral
C---------------------------------------------------------------
      sum =  CalcIntegral(pdfpars(2),pdfpars(3)) 
      do i=1,N-3
         if ( pdfpars(3+i).ne.0 ) then
            sum = sum + pdfpars(3+i) 
     $           * CalcIntegral(pdfpars(2)+i,pdfpars(3))
         endif 
      enddo
      CalcIntXpdfFixN = sum
C---------------------------------------------------------------
      end

* --------------------------------------------------------------
      double precision function CalcIntegral(alpha,beta)
* --------------------------------------------------------------

* Calculates int_0^1 dx x^(alpha) (1-x)^(beta)
* Requires alpha > -1 and beta > -1
* Note: DGamma(x) infinite for x above ~ 170

      implicit none

      double precision alpha,beta
      double precision eps,aa,bb,u,v,uv,DGAMMF

      eps = 1d-5
      aa = alpha+1.d0
      bb = beta+1.d0
      if (aa.le.0d0) aa =eps
      if (bb.le.0d0) bb =eps
      u = DGAMMF(aa)
      v = DGAMMF(bb)
      uv = DGAMMF(aa+bb)
      CalcIntegral = u*v / uv

      return
      end

      double precision function CalcIntegralCheb(ncheb,poly,xminlog
     $     ,iflag)
C
C Created by SG 15 July 2009.
C
C Modified 30 Oct 2009: add a new argument "iflag"
C
C iflag = 0: default param.
C iflag = 1: default * (1-x) param.
C
      implicit none
C
C Keep both simple numerical and analytic formula.
C
      integer ncheb
      double precision poly(ncheb),xminlog
      integer iflag
      double precision zero,one,result,xx,temp
      external gluon
      integer i
      double precision PolyParam,chebint,chebint2,chebint2big
C-----------------------------------------------------
      if (ncheb.le.15) then
         if (iflag.eq.0) then
            result = ChebInt(-2.D0/xminlog,poly)
         elseif (iflag.eq.1) then
            result = ChebInt2(-2.D0/xminlog,poly)
         endif
      elseif (ncheb.le.30 .and. iflag.eq.1 ) then
         result = ChebInt2big(-2.D0/xminlog,poly)
      else
         zero = 1.0D-6
         one  = 1.0D0 - zero
      
         result = 0.

C ssdint(zero,gluon,one)
         do i=1,100
            xx = i/100.D0 - 1./200D0
            temp = PolyParam(xx,ncheb,Poly,xminlog)
            if (iflag.eq.0) then
            else if (iflag.eq.1) then
               temp = temp * ( 1 - xx )  ! new 30 oct 2009, SG.
            endif
            result = result + temp
         enddo
         result = result/100.
         
      endif

      print *,'Cheb integral=',result,' flag=',iflag

      CalcIntegralCheb = result

C-----------------------------------------------------      

      end


      double precision function PolyValInt(NPOLY,Poly)
C-----------------------------------------------------
C 25 Jan 2011
C
C Evaluate integral  \int_0^1 dx ( xf(x)/x )   for PolyVal param of xf(x)
C
      implicit none
      integer NPOLY
      double precision Poly(NPOLY)
      integer i
      double precision sum
C------------------------------------------------
      sum = 0
      do i=NPoly,1,-1
         sum = sum + Poly(i)/i
      enddo

C Jackbian gives factor 3./2. :

      PolyValInt = 3.D0/2.D0 * sum

      end


      double precision function PolyValInt0(NPOLY,Poly)
C-----------------------------------------------------
C 25 Jan 2011
C
C Evaluate integral  \int_0^1 dx ( xf(x) )   for PolyVal param of xf(x)
C
      implicit none
      integer NPOLY
      double precision Poly(NPOLY)
      integer i
      double precision sum
C------------------------------------------------
      sum = 0
      do i=NPoly,1,-1
         sum = sum + Poly(i)/(i+1+0.5)
      enddo

C Jackbian gives factor 3./2. :

      PolyValInt0 = 3.D0/2.D0 * sum

      end



      double precision function ChebInt(a,c)
C--------------------------------------
C
C Automaticaly generated by Maple. Up to 15 polynomials
C
C--------------------------------------
      implicit double precision (t)
      double precision c(15),a
      t1 = c(1)
      t2 = c(12)
      t3 = a ** 2
      t4 = t3 ** 2
      t5 = t4 * a
      t8 = c(2)
      t9 = c(5)
      t10 = c(4)
      t11 = c(3)
      t12 = c(7)
      t13 = c(6)
      t14 = c(8)
      t15 = c(15)
      t16 = t4 ** 2
      t17 = t16 * a
      t20 = c(13)
      t21 = t16 * t4
      t24 = t3 * a
      t25 = t4 * t24
      t30 = t1 - 0.55440D5 * t2 * t5 + t8 + t9 + t10 + t11 + t12 + t13 +
     # t14 + t15 - 0.726485760D9 * t15 * t17 + 0.479001600D9 * t20 * t21
     # - 0.1663200D7 * t2 * t25 + t20 + 0.182D3 * t15 * t3
      t33 = c(10)
      t50 = c(14)
      t57 = 0.121080960D9 * t15 * t16 + t33 + 0.110D3 * t2 * t3 + 0.3024
     #D4 * t33 * t4 - 0.15120D5 * t33 * t5 + 0.30D2 * t12 * t3 - 0.4D1 *
     # t9 * a + t2 - 0.504D3 * t33 * t24 - 0.181440D6 * t33 * t25 - 0.99
     #0D3 * t2 * t24 + t50 + 0.51891840D8 * t50 * t16 + 0.360D3 * t12 * 
     #t4 - 0.3991680D7 * t20 * t25
      t61 = t16 * t3
      t70 = t16 * t5
      t73 = t4 * t3
      t80 = c(9)
      t87 = c(11)
      t92 = 0.6227020800D10 * t50 * t21 + 0.1037836800D10 * t50 * t61 - 
     #0.8648640D7 * t50 * t25 - 0.6D1 * t12 * a - 0.259459200D9 * t50 * 
     #t17 - 0.8717829120D11 * t15 * t70 + 0.1235520D7 * t50 * t73 - 0.12
     #0D3 * t12 * t24 + 0.720D3 * t12 * t73 + t80 + 0.24D2 * t9 * t4 - 0
     #.120D3 * t13 * t5 - 0.5D1 * t13 * a + 0.151200D6 * t87 * t73 + 0.2
     #0D2 * t13 * t3
      t123 = -0.720D3 * t12 * t5 - 0.2D1 * t11 * a + 0.17160D5 * t50 * t
     #4 + 0.90D2 * t87 * t3 + 0.42D2 * t14 * t3 + 0.665280D6 * t20 * t73
     # - 0.154440D6 * t50 * t5 - 0.8D1 * t80 * a + 0.56D2 * t80 * t3 + 0
     #.40320D5 * t80 * t16 - 0.336D3 * t80 * t24 + 0.19958400D8 * t20 * 
     #t16 + 0.5040D4 * t14 * t73 + 0.1680D4 * t80 * t4 + 0.2162160D7 * t
     #15 * t73
      t144 = t16 * t24
      t157 = -0.9D1 * t33 * a - 0.6720D4 * t80 * t5 + 0.20160D5 * t80 * 
     #t73 - 0.11D2 * t2 * a - 0.2520D4 * t14 * t5 + 0.11880D5 * t20 * t4
     # - 0.10D2 * t87 * a + 0.4358914560D11 * t15 * t21 - 0.240240D6 * t
     #15 * t5 - 0.1452971520D11 * t15 * t144 - 0.24D2 * t9 * t24 - 0.210
     #D3 * t14 * t24 + 0.156D3 * t50 * t3 - 0.1716D4 * t50 * t24 - 0.622
     #7020800D10 * t50 * t70
      t188 = 0.12D2 * t9 * t3 - 0.39916800D8 * t2 * t144 - 0.6D1 * t10 *
     # t24 + 0.6D1 * t10 * t3 - 0.14D2 * t15 * a - 0.3D1 * t10 * a + 0.6
     #652800D7 * t2 * t16 - 0.95040D5 * t20 * t5 + 0.2D1 * t11 * t3 - 0.
     #3628800D7 * t87 * t17 + 0.7920D4 * t2 * t4 + 0.39916800D8 * t2 * t
     #61 - 0.5040D4 * t14 * t25 - 0.7D1 * t14 * a - 0.40320D5 * t80 * t2
     #5
      t219 = -0.1D1 * t8 * a - 0.604800D6 * t87 * t25 + 0.24024D5 * t15 
     #* t4 - 0.79833600D8 * t20 * t17 - 0.362880D6 * t33 * t17 - 0.12D2 
     #* t20 * a + 0.60480D5 * t33 * t73 - 0.1320D4 * t20 * t24 + 0.36288
     #0D6 * t33 * t16 - 0.60D2 * t13 * t24 + 0.120D3 * t13 * t4 + 0.5040
     #D4 * t87 * t4 + 0.239500800D9 * t20 * t61 + t87 + 0.8717829120D11 
     #* t15 * t16 * t73
      t250 = -0.479001600D9 * t20 * t144 + 0.3628800D7 * t87 * t61 + 0.3
     #632428800D10 * t15 * t61 + 0.72D2 * t33 * t3 + 0.132D3 * t20 * t3 
     #- 0.30240D5 * t87 * t5 - 0.13D2 * t50 * a - 0.17297280D8 * t15 * t
     #25 + 0.1814400D7 * t87 * t16 + 0.840D3 * t14 * t4 + 0.332640D6 * t
     #2 * t73 - 0.3113510400D10 * t50 * t144 - 0.720D3 * t87 * t24 - 0.2
     #184D4 * t15 * t24 - 0.19958400D8 * t2 * t17
      t253 = t30 + t57 + t92 + t123 + t157 + t188 + t219 + t250

      ChebInt = t253

      end

      double precision function chebint2(a,c)
      implicit double precision (t)
      double precision c(15),a
C
C Output of 
C
C with(CodeGeneration);
C Fortran(int((1-x)*sum(c[i]*(a*log(x)+1)^(i-1),i=1..15),x=0..1.0),optimize); 
C
C--------------------------
      t1 = a ** 2
      t2 = t1 ** 2
      t3 = t2 ** 2
      t4 = t3 * t2
      t7 = t1 * a
      t14 = t2 * a
      t15 = t3 * t14
      t20 = t2 * t1
      t23 = t2 * t7
      t40 = 0.6226260666D10 * c(14) * t4 - 0.675D3 * c(11) * t7 + 0.1991
     #941875D8 * c(13) * t3 - 0.1050000000D2 * c(15) * a - 0.8717297026D
     #11 * c(15) * t15 + 0.1208444738D9 * c(15) * t3 + 0.7143750000D3 * 
     #c(7) * t20 - 0.8614856250D7 * c(14) * t23 + 0.1225867500D7 * c(14)
     # * t20 + 0.1050000000D2 * c(5) * t1 - 0.9750000000D1 * c(14) * a -
     # 0.1520268750D6 * c(14) * t14 + 0.1662375000D5 * c(14) * t2 - 0.22
     #50000000D2 * c(5) * t7 - 0.1656703125D7 * c(12) * t23
      t51 = t3 * t1
      t62 = t3 * a
      t73 = -0.5625000000D2 * c(6) * t7 - 0.2976750000D5 * c(11) * t14 +
     # 0.3675000000D2 * c(8) * t1 + 0.1592500000D3 * c(15) * t1 - 0.1237
     #500000D4 * c(13) * t7 + 0.3989730938D8 * c(12) * t51 - 0.9D1 * c(1
     #3) * a + 0.2000250000D5 * c(9) * t20 - 0.93555D5 * c(13) * t14 - 0
     #.3976087500D7 * c(13) * t23 - 0.3625256250D6 * c(10) * t62 - 0.112
     #5000000D3 * c(7) * t7 + 0.2325000000D2 * c(5) * t2 + 0.2145268125D
     #7 * c(15) * t20 + 0.3627028125D7 * c(11) * t51
      t83 = t3 * t7
      t106 = -0.7975563750D8 * c(13) * t62 - 0.1181250000D3 * c(6) * t14
     # + 0.9625000000D2 * c(12) * t1 + 0.4024125000D5 * c(9) * t3 - 0.39
     #90705469D8 * c(12) * t83 - 0.1452616791D11 * c(15) * t83 - 0.75000
     #00000D1 * c(11) * a + 0.4882500000D4 * c(11) * t2 - 0.5625000000D1
     # * c(4) * t7 + 0.3621712500D6 * c(10) * t3 + 0.7875000000D2 * c(11
     #) * t1 - 0.7257763012D9 * c(15) * t62 - 0.2480625000D4 * c(8) * t1
     #4 - 0.1488375000D5 * c(10) * t14 - 0.4788846562D9 * c(13) * t83
      t130 = 0.2929500000D4 * c(10) * t2 - 0.6D1 * c(9) * a + 0.49D2 * c
     #(9) * t1 - 0.1500000000D1 * c(3) * a - 0.9281250000D3 * c(12) * t7
     # - 0.315D3 * c(9) * t7 + 0.5000000000D0 * c(9) + 0.5000000000D0 * 
     #c(15) - 0.4016250000D5 * c(9) * t23 + 0.5000000000D0 * c(8) + 0.50
     #00000000D0 * c(13) + 0.5000000000D0 * c(10) + 0.5000000000D0 * c(4
     #) + 0.5000000000D0 * c(11) + 0.3630655153D10 * c(15) * t51
      t158 = 0.5000000000D0 * c(6) + 0.5250000000D1 * c(4) * t1 + 0.5000
     #000000D0 * c(14) + 0.5000000000D0 * c(12) - 0.6226640733D10 * c(14
     #) * t15 + 0.3300412500D6 * c(12) * t20 - 0.3625256250D7 * c(11) * 
     #t62 - 0.1608750000D4 * c(14) * t7 + 0.1162500000D3 * c(6) * t2 + 0
     #.6600825000D6 * c(13) * t20 + 0.5000000000D0 * c(2) - 0.1807312500
     #D6 * c(10) * t23 - 0.1722971250D8 * c(15) * t23 + 0.2327325000D5 *
     # c(15) * t2 + 0.5000000000D0 * c(7)
      t187 = 0.2625000000D2 * c(7) * t1 - 0.3750000000D1 * c(6) * a - 0.
     #4500000000D1 * c(7) * a - 0.5457375000D5 * c(12) * t14 + 0.2393838
     #562D9 * c(13) * t51 - 0.1993890938D8 * c(12) * t62 + 0.1750000000D
     #2 * c(6) * t1 + 0.5000000000D0 * c(5) + 0.5000000000D0 * c(3) - 0.
     #3112750266D10 * c(14) * t83 + 0.5179048875D8 * c(14) * t3 + 0.4789
     #431281D9 * c(13) * t4 - 0.4725000000D3 * c(10) * t7 - 0.2047500000
     #D4 * c(15) * t7 - 0.2364862500D6 * c(15) * t14
      t218 = -0.7500000000D0 * c(2) * a - 0.1968750000D3 * c(8) * t7 - 0
     #.5020312500D4 * c(8) * t23 + 0.1810856250D7 * c(11) * t3 - 0.52500
     #00000D1 * c(8) * a + 0.1037330044D10 * c(14) * t51 - 0.2250000000D
     #1 * c(4) * a - 0.3D1 * c(5) * a + 0.5000625000D4 * c(8) * t20 + 0.
     #8137500000D3 * c(8) * t2 + 0.5000000000D0 * c(1) - 0.7087500000D3 
     #* c(7) * t14 + 0.3487500000D3 * c(7) * t2 + 0.1365000000D3 * c(14)
     # * t1 - 0.6615D4 * c(9) * t14
      t250 = 0.63D2 * c(10) * t1 - 0.6750000000D1 * c(10) * a + 0.435838
     #2466D11 * c(15) * t4 + 0.6000750000D5 * c(10) * t20 + 0.6639806250
     #D7 * c(12) * t3 - 0.2592058219D9 * c(14) * t62 + 0.8717563073D11 *
     # c(15) * t3 * t20 + 0.1500187500D6 * c(11) * t20 + 0.1155000000D3 
     #* c(13) * t1 + 0.7672500000D4 * c(12) * t2 + 0.1627500000D4 * c(9)
     # * t2 + 0.1150875000D5 * c(13) * t2 - 0.6024375000D6 * c(11) * t23
     # - 0.8250000000D1 * c(12) * a + 0.1750000000D1 * c(3) * t1
      t253 = t40 + t73 + t106 + t130 + t158 + t187 + t218 + t250
      ChebInt2 = t253
      end


      double precision function chebint2big(a,c)
      implicit double precision (t)
      double precision c(30),a
C
C Output of 
C
C with(CodeGeneration);
C Fortran(int((1-x)*sum(c[i]*(a*log(x)+1)^(i-1),i=1..30),x=0..1.0),optimize); 
C
C--------------------------

     
      t1 = a ** 2
      t2 = t1 ** 2
      t3 = t2 ** 2
      t4 = t3 ** 2
      t5 = t4 * a
      t8 = t1 * a
      t9 = t3 * t8
      t12 = t2 * a
      t17 = t3 * a
      t22 = t4 * t8
      t25 = t2 * t1
      t28 = t4 * t1
      t31 = t2 * t8
      t32 = t3 * t31
      t43 = -0.2128781136D19 * c(22) * t5 - 0.4788846562D9 * c(13) * t9 
     #- 0.6275981250D7 * c(26) * t12 + 0.1050000000D2 * c(5) * t1 - 0.18
     #14440753D10 * c(16) * t17 + 0.8137500000D3 * c(8) * t2 - 0.1077166
     #337D22 * c(24) * t22 + 0.1326165750D8 * c(19) * t25 + 0.8515140787
     #D19 * c(22) * t28 - 0.1014200549D21 * c(30) * t32 + 0.2625000000D2
     # * c(7) * t1 - 0.5625000000D2 * c(6) * t8 - 0.1968750000D3 * c(8) 
     #* t8 - 0.5250000000D1 * a * c(8)
      t44 = t4 * t3
      t57 = t3 * t2
      t58 = t4 * t57
      t67 = t4 * t2
      t74 = t3 * t1
      t79 = 0.6204483832D24 * c(25) * t44 + 0.1662375000D5 * c(14) * t2 
     #+ 0.1592500000D3 * c(15) * t1 + 0.1225867500D7 * c(14) * t25 - 0.1
     #050000000D2 * a * c(15) - 0.8614856250D7 * c(14) * t31 + 0.8841761
     #977D31 * c(30) * t58 - 0.3974788125D7 * c(24) * t12 + 0.5069190262
     #D10 * c(21) * t3 + 0.238D3 * c(18) * t1 + 0.2432900848D19 * c(21) 
     #* t67 + 0.8325880287D16 * c(28) * t57 - 0.9281250000D3 * c(12) * t
     #8 + 0.4149559559D13 * c(24) * t74 + 0.4308667402D22 * c(24) * t67
      t85 = t3 * t25
      t94 = t4 * t74
      t101 = t3 * t12
      t112 = -0.2529635062D9 * c(20) * t31 - 0.4896140581D20 * c(29) * t
     #32 + 0.7123905368D17 * c(24) * t85 - 0.6615D4 * c(9) * t12 - 0.140
     #7585670D14 * c(22) * t9 - 0.5457375000D5 * c(12) * t12 + 0.1088886
     #937D29 * c(28) * t94 - 0.6082231818D17 * c(20) * t5 - 0.2047500000
     #D4 * c(15) * t8 - 0.1248958278D18 * c(28) * t101 - 0.1125000000D2 
     #* a * c(16) - 0.1722971250D8 * c(15) * t31 + 0.2904524122D11 * c(1
     #7) * t74 - 0.4938897088D12 * c(18) * t9
      t121 = t4 * t25
      t126 = t4 * t17
      t145 = 0.3630655153D10 * c(15) * t74 + 0.90117D5 * c(20) * t2 + 0.
     #2179191233D12 * c(16) * t57 + 0.2490647949D16 * c(26) * t57 + 0.16
     #80380888D26 * c(27) * t121 - 0.18D2 * a * c(25) - 0.1551120981D26 
     #* c(26) * t126 - 0.2976750000D5 * c(11) * t12 + 0.4882500000D4 * c
     #(11) * t2 + 0.7368134775D29 * c(30) * t44 + 0.1837500000D3 * c(16)
     # * t1 + 0.3575446875D7 * c(16) * t25 - 0.7638066718D22 * c(29) * t
     #5 - 0.2230126145D18 * c(23) * t32 - 0.1762599589D11 * c(19) * t17
      t160 = t4 * t31
      t177 = 0.1216448684D19 * c(21) * t28 - 0.2027387404D17 * c(21) * t
     #32 - 0.2055375000D5 * c(30) * t8 - 0.9961875000D4 * c(24) * t8 - 0
     #.3230571094D8 * c(16) * t31 - 0.11385D5 * c(25) * t8 - 0.775560455
     #9D25 * c(26) * t160 - 0.1216449844D18 * c(20) * t22 - 0.4274408444
     #D19 * c(26) * t32 - 0.4225642174D18 * c(30) * t101 - 0.3990705469D
     #8 * c(12) * t9 - 0.7309575000D6 * c(18) * t12 + 0.1754317647D28 * 
     #c(30) * t121 + 0.2470545000D6 * c(25) * t2
      t194 = t4 * t9
      t197 = t4 * t12
      t210 = 0.7113530672D13 * c(25) * t74 - 0.1181250000D3 * c(6) * t12
     # + 0.6639806250D7 * c(12) * t3 + 0.5179048875D8 * c(14) * t3 - 0.7
     #123687957D16 * c(24) * t101 - 0.4740037048D12 * c(25) * t17 + 0.20
     #16457246D27 * c(27) * t44 + 0.7143750000D3 * c(7) * t25 - 0.442088
     #0980D31 * c(30) * t194 - 0.2192896797D27 * c(30) * t197 - 0.362525
     #6250D6 * c(10) * t17 - 0.6463002644D24 * c(26) * t197 - 0.21543326
     #74D23 * c(26) * t22 - 0.7095855915D17 * c(22) * t32 - 0.2273208127
     #D20 * c(28) * t32
      t240 = 0.1644705562D9 * c(27) * t25 + 0.5129329267D19 * c(24) * t4
     # + 0.3201162430D16 * c(19) * t4 + 0.1473626988D31 * c(30) * t94 + 
     #0.2787615144D17 * c(23) * t85 - 0.2503928239D13 * c(29) * t17 + 0.
     #1208444738D9 * c(15) * t3 + 0.483D3 * c(25) * t1 + 0.4080375000D6 
     #* c(28) * t2 + 0.525D3 * c(26) * t1 + 0.2325000000D2 * c(5) * t2 -
     # 0.1292600529D23 * c(24) * t197 - 0.5447312965D11 * c(16) * t9 + 0
     #.1587114967D12 * c(19) * t74
      t271 = 0.2727870565D21 * c(28) * t4 + 0.1013701436D18 * c(21) * t4
     # + 0.210D3 * c(17) * t1 - 0.4537028667D27 * c(28) * t160 - 0.30162
     #55007D13 * c(20) * t9 - 0.2403725625D7 * c(22) * t12 - 0.137355750
     #0D7 * c(20) * t12 + 0.1810856250D7 * c(11) * t3 + 0.9074056794D26 
     #* c(28) * t121 + 0.1279310852D13 * c(22) * t74 + 0.6226260666D10 *
     # c(14) * t57 - 0.3150D4 * c(17) * t8 - 0.2025000000D2 * a * c(28) 
     #+ 0.38764845D8 * c(22) * t25 + 0.1295136934D16 * c(25) * t57
      t303 = 0.9625000000D2 * c(12) * t1 + 0.1365000000D3 * c(14) * t1 +
     # 0.3300412500D6 * c(12) * t25 + 0.1155000000D3 * c(13) * t1 + 0.17
     #50000000D1 * c(3) * t1 - 0.1500000000D1 * a * c(3) + 0.7672500000D
     #4 * c(12) * t2 + 0.7264940961D14 * c(30) * t74 + 0.1407757536D15 *
     # c(22) * t57 - 0.93555D5 * c(13) * t12 - 0.1402793438D8 * c(30) * 
     #t12 - 0.2554544672D20 * c(22) * t22 + 0.5601267623D24 * c(27) * t6
     #7 + 0.2145268125D7 * c(15) * t25
      t332 = 0.2677500000D3 * c(19) * t1 - 0.15D2 * a * c(21) + 0.332500
     #0000D3 * c(21) * t1 + 0.3675000000D3 * c(22) * t1 - 0.1575000000D2
     # * a * c(22) - 0.1350000000D2 * a * c(19) + 0.6142500000D3 * c(28)
     # * t1 + 0.4427500000D3 * c(24) * t1 - 0.5020312500D4 * c(8) * t31 
     #- 0.315D3 * c(9) * t8 - 0.12D2 * a * c(17) + 0.1750000000D2 * c(6)
     # * t1 - 0.3750000000D1 * a * c(6) + 0.5000000000D0 * c(12) + 0.500
     #0000000D0 * c(13)
      t348 = 0.5000000000D0 * c(14) + 0.5000000000D0 * c(15) + 0.5000000
     #000D0 * c(16) + 0.5000000000D0 * c(17) + 0.5000000000D0 * c(18) + 
     #0.5000000000D0 * c(19) + 0.5000000000D0 * c(20) + 0.5000000000D0 *
     # c(21) + 0.5000000000D0 * c(22) + 0.5000000000D0 * c(23) + 0.50000
     #00000D0 * c(24) + 0.5000000000D0 * c(25) + 0.5000000000D0 * c(26) 
     #+ 0.5000000000D0 * c(27)
      t375 = 0.5000000000D0 * c(1) + 0.5000000000D0 * c(2) + 0.500000000
     #0D0 * c(28) + 0.5000000000D0 * c(29) + 0.5000000000D0 * c(30) - 0.
     #5944250812D10 * c(29) * t31 - 0.1650000000D2 * a * c(23) + 0.40241
     #25000D5 * c(9) * t3 - 0.7500000000D0 * c(2) * a + 0.4683327433D20 
     #* c(23) * t28 - 0.9961375512D14 * c(25) * t9 - 0.8841761986D31 * c
     #(30) * t4 * t101 + 0.3675000000D2 * c(8) * t1 + 0.3885766564D18 * 
     #c(26) * t85 + 0.6402361494D16 * c(19) * t28
      t406 = -0.856184175D9 * c(23) * t31 + 0.4759788906D14 * c(29) * t7
     #4 - 0.1034080423D24 * c(25) * t197 - 0.1488375000D5 * c(10) * t12 
     #+ 0.6000750000D5 * c(10) * t25 + 0.5109091781D20 * c(22) * t67 + 0
     #.6761233822D19 * c(30) * t85 - 0.6702788905D13 * c(21) * t9 + 0.63
     #65031317D21 * c(29) * t4 - 0.5837619375D9 * c(22) * t31 - 0.936663
     #7000D19 * c(23) * t5 - 0.7087500000D3 * c(7) * t12 - 0.5625000000D
     #1 * c(4) * t8 + 0.63D2 * c(10) * t1
      t437 = 0.1500187500D6 * c(11) * t25 - 0.1293750000D5 * c(26) * t8 
     #+ 0.3041514158D10 * c(20) * t3 + 0.3059864297D14 * c(28) * t74 - 0
     #.6537972769D12 * c(16) * t101 + 0.1124000594D22 * c(23) * t121 + 0
     #.5522107500D6 * c(30) * t2 - 0.11609325D8 * c(29) * t12 - 0.544443
     #4644D28 * c(28) * t126 - 0.389174625D9 * c(21) * t31 + 0.170073750
     #0D6 * c(23) * t2 + 0.1524441712D30 * c(29) * t94 + 0.1037330044D10
     # * c(14) * t74 + 0.1419891602D22 * c(30) * t4 - 0.1875000000D2 * a
     # * c(26)
      t467 = -0.1656703125D7 * c(12) * t31 + 0.4042500000D3 * c(23) * t1
     # + 0.1046107569D14 * c(17) * t85 - 0.3486918810D13 * c(17) * t101 
     #- 0.3048883435D30 * c(29) * t194 - 0.1111358914D22 * c(27) * t5 + 
     #0.1089196546D11 * c(16) * t74 + 0.3097066580D15 * c(23) * t57 - 0.
     #1270002108D13 * c(19) * t9 - 0.3630695947D13 * c(30) * t17 - 0.160
     #8750000D4 * c(14) * t8 + 0.3627028125D7 * c(11) * t74 - 0.31127502
     #66D10 * c(14) * t9 + 0.1265158125D9 * c(26) * t25
      t498 = 0.4760437500D6 * c(29) * t2 + 0.2691336375D9 * c(29) * t25 
     #- 0.9763503750D8 * c(18) * t31 + 0.7875000000D2 * c(11) * t1 + 0.7
     #053844298D11 * c(18) * t74 - 0.6411612666D18 * c(24) * t32 - 0.825
     #0000000D1 * a * c(12) - 0.1125000000D3 * c(7) * t8 - 0.1512342619D
     #26 * c(28) * t197 - 0.1520268750D6 * c(14) * t12 - 0.8662500000D4 
     #* c(23) * t8 - 0.7257763012D9 * c(15) * t17 - 0.1689412164D15 * c(
     #20) * t101 - 0.3684067442D30 * c(30) * t126 - 0.1950000000D2 * a *
     # c(27)
      t531 = 0.6286789884D11 * c(27) * t3 - 0.2250000000D2 * c(5) * t8 +
     # 0.27689175D8 * c(21) * t25 + 0.7105000000D3 * c(30) * t1 - 0.3D1 
     #* a * c(5) + 0.3000674791D23 * c(28) * t28 + 0.2585201366D23 * c(2
     #4) * t121 + 0.1727212800D12 * c(30) * t3 + 0.1627500000D4 * c(9) *
     # t2 - 0.4725000000D3 * c(10) * t8 + 0.1709737288D18 * c(25) * t85 
     #- 0.2480625000D4 * c(8) * t12 - 0.8401897427D24 * c(29) * t22 + 0.
     #1286794451D11 * c(23) * t3
      t562 = -0.4500000000D1 * a * c(7) - 0.2559375000D4 * c(16) * t8 - 
     #0.5170398417D22 * c(25) * t22 + 0.2485520145D17 * c(30) * t57 - 0.
     #1228022426D29 * c(30) * t160 + 0.5720715D7 * c(17) * t25 + 0.42575
     #46032D18 * c(22) * t4 + 0.1551120958D26 * c(26) * t44 - 0.33023615
     #62D10 * c(27) * t31 + 0.2959627238D11 * c(25) * t3 + 0.3487500000D
     #3 * c(7) * t2 + 0.9782647875D9 * c(18) * t3 + 0.3621712500D6 * c(1
     #0) * t3 + 0.1126462500D6 * c(21) * t2 - 0.2540736054D28 * c(29) * 
     #t160
      t592 = 0.2160488940D25 * c(28) * t67 + 0.4234559837D27 * c(29) * t
     #121 + 0.2589524438D9 * c(16) * t3 - 0.5334985780D14 * c(19) * t101
     # - 0.4147293150D10 * c(17) * t17 - 0.8569712756D15 * c(29) * t9 + 
     #0.9615201750D8 * c(25) * t25 + 0.2000250000D5 * c(9) * t25 - 0.296
     #2523155D12 * c(24) * t17 + 0.1000224930D23 * c(27) * t28 + 0.42315
     #D5 * c(17) * t2 + 0.5250000000D1 * c(4) * t1 - 0.2250000000D1 * a 
     #* c(4) - 0.3556860713D15 * c(18) * t5
      t623 = -0.159766425D9 * c(19) * t31 + 0.8891100231D13 * c(19) * t5
     #7 - 0.2092247063D14 * c(17) * t32 - 0.1831410D7 * c(21) * t12 + 0.
     #4032914581D27 * c(27) * t94 - 0.1554259191D17 * c(25) * t101 + 0.3
     #378927447D16 * c(21) * t85 - 0.2413264219D10 * c(26) * t31 + 0.111
     #1354674D21 * c(27) * t4 - 0.1452616791D11 * c(15) * t9 - 0.2364862
     #500D6 * c(15) * t12 + 0.4352392997D11 * c(26) * t3 + 0.2585200441D
     #23 * c(25) * t67 - 0.1065571570D12 * c(22) * t17 + 0.5620000959D21
     # * c(23) * t67
      t654 = 0.8419160889D18 * c(27) * t85 - 0.5068468510D16 * c(20) * t
     #32 - 0.1307654414D13 * c(16) * t32 + 0.1013678234D16 * c(20) * t85
     # - 0.7406307887D12 * c(26) * t17 + 0.55335D5 * c(18) * t2 - 0.5450
     #625000D4 * c(20) * t8 + 0.2992500000D3 * c(20) * t1 - 0.1425000000
     #D2 * a * c(20) + 0.1991941875D8 * c(13) * t3 + 0.3173625000D5 * c(
     #16) * t2 - 0.7770262500D7 * c(27) * t12 + 0.2345403229D13 * c(23) 
     #* t74 + 0.1926581224D14 * c(27) * t74
      t685 = 0.3556847144D15 * c(18) * t4 + 0.5000625000D4 * c(8) * t25 
     #+ 0.8188691962D10 * c(22) * t3 - 0.3238039980D17 * c(26) * t101 + 
     #0.2092263026D14 * c(17) * t4 + 0.4358382466D11 * c(15) * t57 + 0.5
     #179048875D9 * c(17) * t3 - 0.3097255633D16 * c(23) * t101 - 0.5743
     #2375D8 * c(17) * t31 + 0.8716764932D12 * c(17) * t57 - 0.243655025
     #4D25 * c(30) * t22 - 0.1845866123D23 * c(30) * t5 - 0.7975563750D8
     # * c(13) * t17 - 0.1088886941D29 * c(28) * t194 - 0.2700609887D24 
     #* c(28) * t22
      t715 = -0.1807312500D6 * c(10) * t31 - 0.675D3 * c(11) * t8 + 0.12
     #70368065D29 * c(29) * t44 + 0.1973084825D11 * c(24) * t3 - 0.21D2 
     #* a * c(29) + 0.2929500000D4 * c(10) * t2 - 0.3000669068D22 * c(28
     #) * t5 - 0.1709763378D19 * c(25) * t32 + 0.2215043573D24 * c(30) *
     # t28 - 0.1067046002D16 * c(19) * t32 + 0.6475684668D15 * c(24) * t
     #57 - 0.2175000000D2 * a * c(30) - 0.2432899688D19 * c(21) * t22 - 
     #0.2592058219D9 * c(14) * t17
      t746 = -0.2331388786D18 * c(29) * t101 + 0.1938242250D8 * c(20) * 
     #t25 + 0.1162500000D3 * c(6) * t2 + 0.4789431281D9 * c(13) * t57 - 
     #0.9D1 * a * c(13) + 0.8617322477D21 * c(25) * t28 - 0.4458188109D1
     #0 * c(28) * t31 - 0.3825D4 * c(18) * t8 - 0.9536231250D7 * c(28) *
     # t12 - 0.1699094162D13 * c(28) * t17 + 0.6600825000D6 * c(13) * t2
     #5 + 0.1391512500D6 * c(22) * t2 - 0.7481250000D4 * c(22) * t8 - 0.
     #6402349283D16 * c(19) * t5 - 0.6049370475D26 * c(29) * t197
      t778 = 0.1457029050D17 * c(29) * t57 - 0.5081472334D29 * c(29) * t
     #126 + 0.6615000000D3 * c(29) * t1 - 0.1842750000D5 * c(29) * t8 + 
     #0.1760876618D10 * c(19) * t3 + 0.1013678234D17 * c(22) * t85 - 0.3
     #348939219D11 * c(20) * t17 + 0.8717563073D11 * c(15) * t85 - 0.520
     #3039888D15 * c(28) * t9 - 0.3976087500D7 * c(13) * t31 + 0.6701152
     #083D12 * c(21) * t74 + 0.2667574300D15 * c(19) * t85 - 0.336076137
     #5D25 * c(27) * t197 - 0.7500000000D1 * a * c(11)
      t809 = 0.1748594954D19 * c(28) * t85 + 0.2114621438D9 * c(28) * t2
     #5 + 0.8933859309D11 * c(28) * t3 - 0.3625256250D7 * c(11) * t17 - 
     #0.6088980398D11 * c(21) * t17 - 0.1743140149D12 * c(17) * t9 - 0.5
     #15970D6 * c(17) * t12 + 0.5330166188D8 * c(23) * t25 - 0.602437500
     #0D6 * c(11) * t31 + 0.1292600221D24 * c(26) * t67 - 0.6750000000D1
     # * a * c(10) - 0.3110703750D7 * c(23) * t12 + 0.2154330619D21 * c(
     #24) * t28 + 0.1561100212D19 * c(23) * t4 - 0.2815171340D14 * c(23)
     # * t9
      t839 = -0.1267059123D16 * c(22) * t101 + 0.3102241639D24 * c(25) *
     # t121 - 0.4016250000D5 * c(9) * t31 + 0.4625489048D16 * c(27) * t5
     #7 + 0.8401889415D23 * c(29) * t28 + 0.2413298634D14 * c(20) * t57 
     #+ 0.6033246585D14 * c(21) * t57 - 0.1481940494D14 * c(18) * t101 +
     # 0.3350576041D12 * c(20) * t74 + 0.1814811521D28 * c(28) * t44 - 0
     #.4032914551D27 * c(27) * t126 + 0.7561711290D25 * c(29) * t67 + 0.
     #8841105D7 * c(18) * t25 - 0.6226640733D10 * c(14) * t101
      t861 = -0.1230764752D10 * c(24) * t31 + 0.5000000000D0 * c(3) + 0.
     #5000000000D0 * c(4) + 0.5000000000D0 * c(5) + 0.5000000000D0 * c(6
     #) + 0.5000000000D0 * c(7) + 0.5000000000D0 * c(8) + 0.5000000000D0
     # * c(9) + 0.5000000000D0 * c(10) + 0.5000000000D0 * c(11) + 0.3077
     #615170D22 * c(26) * t28 - 0.6D1 * a * c(9) - 0.1778817056D15 * c(2
     #6) * t9 + 0.2327325000D5 * c(15) * t2 - 0.1993890938D8 * c(12) * t
     #17
      t892 = -0.14625D5 * c(27) * t8 + 0.4274441055D20 * c(26) * t4 + 0.
     #2941125000D6 * c(26) * t2 - 0.2585201520D23 * c(24) * t160 - 0.127
     #5000000D2 * a * c(18) + 0.71145D5 * c(19) * t2 + 0.2027402872D17 *
     # c(20) * t4 - 0.1132729442D13 * c(27) * t17 - 0.4054821212D18 * c(
     #21) * t5 + 0.1307634461D13 * c(16) * t85 - 0.3547293750D6 * c(16) 
     #* t12 - 0.1737550238D10 * c(25) * t31 + 0.1250740303D12 * c(29) * 
     #t3 + 0.3048883440D30 * c(29) * t58
      t923 = -0.5109092999D20 * c(22) * t197 - 0.5395745069D14 * c(24) *
     # t9 - 0.8001807074D23 * c(27) * t22 + 0.2963700077D13 * c(18) * t5
     #7 + 0.3393424125D9 * c(30) * t25 - 0.1645312500D5 * c(28) * t8 - 0
     #.3083282896D15 * c(27) * t9 - 0.1725000000D2 * a * c(24) + 0.49D2 
     #* c(9) * t1 - 0.8812997944D10 * c(18) * t17 - 0.3590544184D20 * c(
     #24) * t5 - 0.6204483648D24 * c(25) * t160 - 0.5020785D7 * c(25) * 
     #t12 - 0.1873332760D21 * c(23) * t22 + 0.1185588445D14 * c(26) * t7
     #4
      t955 = -0.1124000460D22 * c(23) * t197 + 0.5687500000D3 * c(27) * 
     #t1 + 0.3475875000D6 * c(27) * t2 - 0.6721523951D26 * c(27) * t160 
     #- 0.8717297026D11 * c(15) * t101 - 0.9750000000D1 * a * c(14) + 0.
     #1538798780D20 * c(25) * t4 + 0.2585201366D25 * c(26) * t121 + 0.24
     #36551416D26 * c(30) * t67 - 0.6476079961D17 * c(27) * t101 + 0.592
     #7942890D14 * c(18) * t85 + 0.2393838562D9 * c(13) * t74 - 0.138067
     #5944D16 * c(30) * t9 - 0.3847011625D21 * c(26) * t5 + 0.2058787500
     #D6 * c(24) * t2
      t986 = 0.7211401312D8 * c(24) * t25 - 0.1231043720D21 * c(25) * t5
     # + 0.3989730938D8 * c(12) * t74 + 0.1150875000D5 * c(13) * t2 - 0.
     #1237500000D4 * c(13) * t8 - 0.1012095D7 * c(19) * t12 - 0.18032749
     #64D12 * c(23) * t17 - 0.4590D4 * c(19) * t8 - 0.1010314723D20 * c(
     #27) * t32 - 0.1778410004D15 * c(18) * t32 - 0.6412500000D4 * c(21)
     # * t8 - 0.4826891896D15 * c(21) * t101 + 0.3497189908D19 * c(29) *
     # t85 - 0.7835603344D10 * c(30) * t31 + 0.1216448684D18 * c(20) * t
     #28
      t991 = t240 + t177 + t145 + t839 + t375 + t592 + t778 + t685 + t46
     #7 + t892 + t303 + t623 + t79 + t562 + t746 + t271 + t861 + t654 + 
     #t955 + t986 + t112 + t43 + t437 + t923 + t332 + t809 + t406 + t348
     # + t210 + t531 + t715 + t498

      chebint2big = t991
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      


      Subroutine SumRulesAS
C---------------------------------------------------------------
C 
C  Created 20 Jul 2011 by VR. Add sum-rules for A.Schoening parameterisation 
C
C---------------------------------------------------------------
      implicit none
#include "pdfparam.inc"

      double precision sumUv, sumDv
      double precision sumMom, sumGlue,x
      integer i
      double precision SumRuleASpar,splogn
      double precision ubar,dbar,uval,dval,gluon
      integer IDebug
      data IDebug/1/

C---------------------------------------------------------------

C Counting sum-rule for uv:
      sumUv = SumRuleASpar(-1,asuval)
      asuval(1) = 2.0D0 / sumUv

C Counting sum-rule for dv:
      sumDv = SumRuleASpar(-1,asdval)
      asdval(1) = 1.0D0 / sumDv

C Momentum sum rule:
      sumMom = 2.D0*asubar(1)*SumRuleASpar(0,asubar) +
     $         2.D0*asdbar(1)*SumRuleASpar(0,asdbar) +
     $         asuval(1)*SumRuleASpar(0,asuval) +
     $         asdval(1)*SumRuleASpar(0,asdval) 
      sumGlue = SumRuleASpar(0,asglue)
      asglue(1) = (1.0 - SumMom)/sumGlue

      if (IDebug.eq.1) then
         print '(''uv:'',5F10.4)',(asuval(i),i=1,5)
         print '(''dv:'',5F10.4)',(asdval(i),i=1,5)
         print '(''Ub:'',5F10.4)',(asubar(i),i=1,5)
         print '(''Db:'',5F10.4)',(asdbar(i),i=1,5)
         print '(''GL:'',5F10.4)',(asglue(i),i=1,5)
      endif


      if (IDebug.eq.10) then
         do i=1,8
            x = 10**(-i/2.)
            print '(6F12.5)',x,splogn(x,asuval),splogn(x,asdval),
     $           splogn(x,asubar),splogn(x,asdbar),splogn(x,asglue)
            print '(6F12.5)',x,uval(x),dval(x),ubar(x),dbar(x),gluon(x)
         enddo
      endif

      return
      end

      double precision function SumRuleASpar(n,aas)
C---------------------------------------------------------------
C A wrapper for A. Schoening sum-rule integral.
C 
C parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
C---------------------------------------------------------------
      implicit none 
      integer n
      integer i
      double precision aas(1:5), aass(1:5)
      double precision splognni

      do i=1,5
         aass(i)=aas(i)
      enddo
      aass(1)=1.d0
      aass(2)=aass(2)+n

      SumRuleASpar=splognni(aass)

      if (SumRuleASpar.eq.0) then
         print*, 'sum rule is ZERO---- ERROR'
         STOP
      endif
      return
      end
      double precision function splognni(as)
C---------------------------------------------------------------
*     Numerical Integration of the 
*     Special lognormal function 
*     using the Simpson method
*
*
*     A.Schoening, University Heidelberg, Physikalisches Institut
*     Creation: 12.6.2011
*  
*  ASPDF = A1*x**(A2-A3*log(x))*(1-x)**(A4-A5*log(1-x))
C---------------------------------------------------------------
      implicit none
 
      integer np,i
      data np /1000/
      
      double precision as(1:5)
      double precision xas,xnas
      double precision splogn

      double precision peak,h,sum,xlmin
      data xlmin /-10.0/
      
      logical logflag
      data logflag /.true./

      logical falling
      data falling /.false./
c      data falling /.true./
      double precision eps
      data eps /0.001d0/
      double precision f1,f2,f3


C-----------------------------------------------------

      f1=splogn(eps,as)
      f2=splogn(0.5d0,as)
      f3=splogn((1.d0-eps),as)

      if (f1.gt.f2+f3) then
         logflag=.true.
         falling=.true.
      elseif (f3.gt.f2+f1) then
         logflag=.true.
         falling=.false.
      else
         logflag=.false.
      endif
c      print *,logflag,falling

c linear integration
      if (.not.logflag) then
         h=1.d0/float(np)
         sum=0.d0
c weight 4
         do i=1,np-1,2
            xas=i*h
            sum=sum+splogn(xas,as)
         enddo
         sum=2.d0*sum
c weight 2
         do i=2,np-2,2
            xas=i*h
            sum=sum+splogn(xas,as)
         enddo
         sum=2.d0*sum
         sum=sum+splogn(0.d0,as)+splogn(1.d0,as)
         splognni=h/3.d0*sum
      else
c steeply falling distribution

         if (falling) then

            h=-xlmin/float(np)
            sum=0.d0
c weight 4
            do i=1,np-1,2
               xas=10.d0**(xlmin+i*h)
               sum=sum+xas*splogn(xas,as)
            enddo
            sum=2.d0*sum
c weight 2
            do i=0,np-2,2
               xas=10.d0**(xlmin+i*h)
               sum=sum+xas*splogn(xas,as)
            enddo
            sum=2.d0*sum
            xas=1.d0
            sum=sum+splogn(xas,as)
            splognni=h/3.d0*sum*log(10.d0)


         else
c steeply rising distribution
            h=-xlmin/float(np)
            sum=0.d0
c     weight 4
            do i=1,np-1,2
               xas=10.d0**(xlmin+i*h)
               xnas=1.d0-xas
               sum=sum+xas*splogn(xnas,as)
*     print *,x,x
            enddo
            sum=2.d0*sum
c weight 2
            do i=0,np-2,2
               xas=10.d0**(xlmin+i*h)
               xnas=1.d0-xas
               sum=sum+xas*splogn(xnas,as)
            enddo
            sum=2.d0*sum
            xas=0.d0
            sum=sum+splogn(xas,as)
            splognni=h/3.d0*sum*log(10.d0)
            
         endif

      endif
      if (splognni.eq.0) then
         print*, 'sum rule is ZERO---- Warning', sum, logflag, falling
         splognni=1d-10
      endif
      return
      end




ccccccccccccccccccccccccccccccccccccccccc

      Subroutine SumRulesCTeq
C---------------------------------------------------------------
C 
C  Created 22 Apr 2011 by SG. Add sum-rules for CTEQ-like parameterisation 
C
C---------------------------------------------------------------
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      double precision sumUv, sumDv
      double precision sumMom, sumGlue,x
      integer i
      double precision SumRuleCTEQ,ctpara, tStr
      double precision ubar,dbar,uval,dval,gluon,str
      double precision CalcIntegral
      
      integer IDebug
      data IDebug/1/

C---------------------------------------------------------------

C Counting sum-rule for uv:
      sumUv = SumRuleCTEQ(-1,ctuval)
      ctuval(1) = 2.0D0 / sumUv

C Counting sum-rule for dv:
      sumDv = SumRuleCTEQ(-1,ctdval)
      ctdval(1) = 1.0D0 / sumDv

C Momentum sum rule:
C----------------
C Sea:

      
      if (Index(PDF_DECOMPOSITION,'Str').gt.0) then
         tStr = ctstr(1)*SumRuleCTEQ(0,ctstr)
      else
         tStr = 0               ! Strange already included in Dbar
      endif


      sumMom = 2.D0*ctubar(1)*SumRuleCTEQ(0,ctubar) +
     $     2.D0*ctdbar(1)*SumRuleCTEQ(0,ctdbar) +
     $     ctuval(1)*SumRuleCTEQ(0,ctuval) +
     $     ctdval(1)*SumRuleCTEQ(0,ctdval)+ 2.D0*tStr 
      sumGlue = SumRuleCTEQ(0,ctglue)
      sumMom = sumMom - ctglue(7)*CalcIntegral(ctglue(8),ctglue(9))



      ctglue(1) = (1.0 - SumMom)/sumGlue



      if (IDebug.eq.1) then
         print '(''uv:'',9F10.4)',(ctuval(i),i=1,9)
         print '(''dv:'',9F10.4)',(ctdval(i),i=1,9)
         print '(''Ub:'',9F10.4)',(ctubar(i),i=1,9)
         print '(''Db:'',9F10.4)',(ctdbar(i),i=1,9)
         print '(''GL:'',9F10.4)',(ctglue(i),i=1,9)
         print '(''ST:'',x9F10.4)',(ctstr(i),i=1,9)
      endif
      if (IDebug.eq.10) then
         do i=1,8
            x = 10**(-i/2.)
            print '(7F12.5)',x,ctpara(x,ctuval),ctpara(x,ctdval),
     $           ctpara(x,ctubar),ctpara(x,ctdbar),ctpara(x,ctglue), 
     $           ctpara(x,ctstr)
            print '(7F12.5)',x,uval(x),dval(x),ubar(x),dbar(x),
     $           gluon(x),str(x)
         enddo
      endif

      end


ccccccccccccccccccccccccccccccccccccccccc

      Subroutine SumRulesCTEQHera
C---------------------------------------------------------------
C 
C  Sum-rules for CTEQ-HERA hybrid parameterisation 
C
C---------------------------------------------------------------
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
      double precision sumUv, sumDv
      double precision sumMom, sumGlue,x
      integer i
      double precision SumRuleCTEQhera,ctherapara, tStr
      double precision ubar,dbar,uval,dval,gluon,str
      double precision CalcIntegral
      
      integer IDebug
      data IDebug/1/

C---------------------------------------------------------------

C Counting sum-rule for uv:
      sumUv = SumRuleCTEQhera(-1,ctuval)
      ctuval(1) = 2.0D0 / sumUv

C Counting sum-rule for dv:
      sumDv = SumRuleCTEQhera(-1,ctdval)
      ctdval(1) = 1.0D0 / sumDv

C Momentum sum rule:
C----------------
C Sea:

      
      if (Index(PDF_DECOMPOSITION,'Str').gt.0) then
         tStr = ctstr(1)*SumRuleCTEQhera(0,ctstr)
      else
         tStr = 0               ! Strange already included in Dbar
      endif


      sumMom = 2.D0*ctubar(1)*SumRuleCTEQhera(0,ctubar) +
     $     2.D0*ctdbar(1)*SumRuleCTEQhera(0,ctdbar) +
     $     ctuval(1)*SumRuleCTEQhera(0,ctuval) +
     $     ctdval(1)*SumRuleCTEQhera(0,ctdval)+ 2.D0*tStr 
      sumGlue = SumRuleCTEQhera(0,ctglue)
      sumMom = sumMom - ctglue(7)*CalcIntegral(ctglue(8),ctglue(9))



      ctglue(1) = (1.0 - SumMom)/sumGlue



      if (IDebug.eq.1) then
         print '(''uv:'',9F10.4)',(ctuval(i),i=1,9)
         print '(''dv:'',9F10.4)',(ctdval(i),i=1,9)
         print '(''Ub:'',9F10.4)',(ctubar(i),i=1,9)
         print '(''Db:'',9F10.4)',(ctdbar(i),i=1,9)
         print '(''GL:'',9F10.4)',(ctglue(i),i=1,9)
         print '(''ST:'',x9F10.4)',(ctstr(i),i=1,9)
      endif
      if (IDebug.eq.10) then
         do i=1,8
            x = 10**(-i/2.)
            print '(7F12.5)',x,ctherapara(x,ctuval),ctherapara(x,ctdval),
     $           ctherapara(x,ctubar),ctherapara(x,ctdbar),ctherapara(x,ctglue), 
     $           ctherapara(x,ctstr)
            print '(7F12.5)',x,uval(x),dval(x),ubar(x),dbar(x),
     $           gluon(x),str(x)
         enddo
      endif

      end



      double precision function SumRuleCTEQ(n,acteq)
C---------------------------------------------------------------
C Sum-rule integral for 
C  
C UF = a0*E**(a3*x)*(1 - x)**a2*x**(a1 + n)*(1 + E**a4*x + E**a5*x**2)
C
C parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
C---------------------------------------------------------------
      implicit none 
      integer n
      double precision acteq(1:6)
      double precision YF
      double precision DGammF,HypG1F1r, HypG1F1
C-----------------------------------------------------
      YF = (
     &  DGammF(1 + acteq(3) )*DGammF(1 + acteq(2) + n)*(
     &   Hypg1F1(1 + acteq(2) + n,2 + acteq(2) + acteq(3) + n,
     $     acteq(4)) + 
     &   (1 + acteq(2) + n)*DGammF(2 + acteq(2) + acteq(3) + n)*
     &   ( exp(acteq(5))*
     &    Hypg1F1R(2 + acteq(2) + n,3 + acteq(2)
     $     + acteq(3) + n,acteq(4)) + 
     &     exp(acteq(6))*(2 + acteq(2) + n)*
     &    Hypg1F1R(3 + acteq(2) + n,4 + acteq(2)
     $     + acteq(3) + n,acteq(4))
     &   )
     &  )
     & )/DGammF(2 + acteq(2) + acteq(3) + n)

      SumRuleCTEQ = YF

      end

      double precision function SumRuleCTEQhera(n,acteq)
C---------------------------------------------------------------
C Sum-rule integral for 
C  
C UF = a1*exp(a6*x)*(1 - x)**a3*x**(a2 + n)*(1 + a4*x + a5*x**2)
C
C parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
C---------------------------------------------------------------
      implicit none 
      integer n
      double precision acteq(1:6)
      double precision YF
      double precision DGammF,HypG1F1r, HypG1F1
C-----------------------------------------------------
      YF = (
     &  DGammF(1 + acteq(3) )*DGammF(1 + acteq(2) + n)*(
     &   Hypg1F1(1 + acteq(2) + n,2 + acteq(2) + acteq(3) + n,
     $     acteq(6)) + 
     &   (1 + acteq(2) + n)*DGammF(2 + acteq(2) + acteq(3) + n)*
     &   ( acteq(4)*
     &    Hypg1F1R(2 + acteq(2) + n,3 + acteq(2)
     $     + acteq(3) + n,acteq(6)) + 
     &     acteq(5)*(2 + acteq(2) + n)*
     &    Hypg1F1R(3 + acteq(2) + n,4 + acteq(2)
     $     + acteq(3) + n,acteq(6))
     &   )
     &  )
     & )/DGammF(2 + acteq(2) + acteq(3) + n)

      SumRuleCTEQhera = YF

      end

      double precision function HYPG1F1R(a,b,z)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Regularized Confluent Hypergeometric function 1F1
c
c Implementation: A.Schoening, University Heidelberg
c Date: 08.04.2011
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision HYPG1F1
      double precision a,b,z
      double precision DGammF

      if (b.lt.0) then
         print *,'HYPG1F1R Warning, function not defined for b negative'
         HYPG1F1R=0.0
         return
      endif
      HYPG1F1R=hypg1f1(a,b,z)/DGammF(b)
      return
      end

      double precision function HYPG1F1(a0,b,z0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Confluent Hypergeometric function 1F1
c
c Implementation: A.Schoening, University Heidelberg
c Date: 08.04.2011
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision a0,a,b,z,z0

      integer n
      double precision prec
      data prec /1e-10/

      double precision abnzfact
      double precision factor
      double precision old

c Kummer's transformation
      if (z0.lt.0) then
         a=b-a0
         z=-z0
         factor=exp(z0)
      else
         a=a0
         z=z0
         factor=1.0
      endif

      n=1
      abnzfact=1.0
      hypg1f1=1.0

 10   continue
      abnzfact=abnzfact*(a+n-1.0)/(b+n-1.0)/n*z
      hypg1f1=hypg1f1+abnzfact

c precision reached
      if (abs(abnzfact).lt.abs(hypg1f1)*prec) goto 99
      n=n+1
c      print *,'HYPG1F1: n=',n,' hpyg1f1=',hypg1f1,abnzfact,z
      if (n.lt.1000) goto 10
      print *,'HYPG1F1: no convergence for ',a,b,z

 99   continue
      hypg1f1=hypg1f1*factor

      return
      end

      REAL*8 FUNCTION SSDINT(XL,F,XR)
C-----------------------------------------------------------------------
C     Integrate REAL*8 F over REAL*8 (XL,XR)
C     Note quadrature constants R and W have been converted to explicit
C     REAL*8 (.xxxxxDxx) form.
C
C     Bisset's XINTH
C-----------------------------------------------------------------------
      IMPLICIT NONE
*KEEP,SSLUN.
*KEND.
      EXTERNAL F
      INTEGER NMAX
      REAL*8 TOLABS,TOLREL,XLIMS(200)
      REAL*8 R(93),W(93)
      INTEGER PTR(4),NORD(4)
      INTEGER ICOUNT
      REAL*8 XL,XR,F
      REAL*8 AA,BB,TVAL,VAL,TOL
      INTEGER NLIMS,I,J,KKK
C
      DATA PTR,NORD/4,10,22,46,  6,12,24,48/
      DATA (R(KKK),KKK=1,48)/
     + .2386191860D0,.6612093865D0,.9324695142D0,.1252334085D0,
     + .3678314990D0,.5873179543D0,.7699026742D0,.9041172563D0,
     + .9815606342D0,.0640568929D0,.1911188675D0,.3150426797D0,
     + .4337935076D0,.5454214714D0,.6480936519D0,.7401241916D0,
     + .8200019860D0,.8864155270D0,.9382745520D0,.9747285560D0,
     + .9951872200D0,.0323801710D0,.0970046992D0,.1612223561D0,
     + .2247637903D0,.2873624873D0,.3487558863D0,.4086864820D0,
     + .4669029048D0,.5231609747D0,.5772247261D0,.6288673968D0,
     + .6778723796D0,.7240341309D0,.7671590325D0,.8070662040D0,
     + .8435882616D0,.8765720203D0,.9058791367D0,.9313866907D0,
     + .9529877032D0,.9705915925D0,.9841245837D0,.9935301723D0,
     + .9987710073D0,.0162767488D0,.0488129851D0,.0812974955D0/
      DATA (R(KKK),KKK=49,93)/
     + .1136958501D0,.1459737146D0,.1780968824D0,.2100313105D0,
     + .2417431561D0,.2731988126D0,.3043649444D0,.3352085229D0,
     + .3656968614D0,.3957976498D0,.4254789884D0,.4547094222D0,
     + .4834579739D0,.5116941772D0,.5393881083D0,.5665104186D0,
     + .5930323648D0,.6189258401D0,.6441634037D0,.6687183100D0,
     + .6925645366D0,.7156768123D0,.7380306437D0,.7596023411D0,
     + .7803690438D0,.8003087441D0,.8194003107D0,.8376235112D0,
     + .8549590334D0,.8713885059D0,.8868945174D0,.9014606353D0,
     + .9150714231D0,.9277124567D0,.9393703398D0,.9500327178D0,
     + .9596882914D0,.9683268285D0,.9759391746D0,.9825172636D0,
     + .9880541263D0,.9925439003D0,.9959818430D0,.9983643759D0,
     + .9996895039/
      DATA (W(KKK),KKK=1,48)/ .4679139346D0,.3607615730D0,
     +.1713244924D0,.2491470458D0, .2334925365D0,.2031674267D0,
     +.1600783285D0,.1069393260D0, .0471753364D0,.1279381953D0,
     +.1258374563D0,.1216704729D0, .1155056681D0,.1074442701D0,
     +.0976186521D0,.0861901615D0, .0733464814D0,.0592985849D0,
     +.0442774388D0,.0285313886D0, .0123412298D0,.0647376968D0,
     +.0644661644D0,.0639242386D0, .0631141923D0,.0620394232D0,
     +.0607044392D0,.0591148397D0, .0572772921D0,.0551995037D0,
     +.0528901894D0,.0503590356D0, .0476166585D0,.0446745609D0,
     +.0415450829D0,.0382413511D0, .0347772226D0,.0311672278D0,
     +.0274265097D0,.0235707608D0, .0196161605D0,.0155793157D0,
     +.0114772346D0,.0073275539D0, .0031533461D0,.0325506145D0,
     +.0325161187D0,.0324471637D0/
      DATA (W(KKK),KKK=49,93)/
     + .0323438226D0,.0322062048D0,.0320344562D0,.0318287589D0,
     + .0315893308D0,.0313164256D0,.0310103326D0,.0306713761D0,
     + .0302999154D0,.0298963441D0,.0294610900D0,.0289946142D0,
     + .0284974111D0,.0279700076D0,.0274129627D0,.0268268667D0,
     + .0262123407D0,.0255700360D0,.0249006332D0,.0242048418D0,
     + .0234833991D0,.0227370697D0,.0219666444D0,.0211729399D0,
     + .0203567972D0,.0195190811D0,.0186606796D0,.0177825023D0,
     + .0168854799D0,.0159705629D0,.0150387210D0,.0140909418D0,
     + .0131282296D0,.0121516047D0,.0111621020D0,.0101607705D0,
     + .0091486712D0,.0081268769D0,.0070964708D0,.0060585455D0,
     + .0050142027D0,.0039645543D0,.0029107318D0,.0018539608D0,
     + .0007967921/
C
C      DATA TOLABS,TOLREL,NMAX/1.D-35,5.D-5,100/
C      DATA TOLABS,TOLREL,NMAX/1.D-30,5.D-4,200/

      DATA TOLABS,TOLREL,NMAX/1.D-15,2.D-2,100/

C
      SSDINT=0
      NLIMS=2
      XLIMS(1)=XL
      XLIMS(2)=XR
      ICOUNT=0
C
   10 AA=(XLIMS(NLIMS)-XLIMS(NLIMS-1))/2
      BB=(XLIMS(NLIMS)+XLIMS(NLIMS-1))/2
      TVAL=0
      DO 20 I=1,3
   20 TVAL=TVAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
      TVAL=TVAL*AA
      DO 40 J=1,4
       VAL=0
       DO 30 I=PTR(J),PTR(J)-1+NORD(J)
        ICOUNT=ICOUNT+1
        IF(ICOUNT.GT.1E5) THEN
         WRITE(1,*) 'WARNING IN SSDINT: SET SSDINT TO ZERO'
         WRITE(6,*) 'WARNING IN SSDINT: SET SSDINT TO ZERO'
         SSDINT=0.
         RETURN
        ENDIF
   30  VAL=VAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
       VAL=VAL*AA
       TOL=MAX(TOLABS,TOLREL*ABS(VAL))
       IF(ABS(TVAL-VAL).LT.TOL) THEN
        SSDINT=SSDINT+VAL
        NLIMS=NLIMS-2
        IF (NLIMS.NE.0) GO TO 10
        RETURN
       ENDIF
   40 TVAL=VAL
      IF(NMAX.EQ.2) THEN
       SSDINT=VAL
       RETURN
      END IF
      IF(NLIMS.GT.(NMAX-2)) THEN
       WRITE(1,10000) SSDINT,NMAX,BB-AA,BB+AA
       WRITE(6,10000) SSDINT,NMAX,BB-AA,BB+AA
       RETURN
      ENDIF
      XLIMS(NLIMS+1)=BB
      XLIMS(NLIMS+2)=BB+AA
      XLIMS(NLIMS)=BB
      NLIMS=NLIMS+2
      GO TO 10
C
10000 FORMAT (' SSDINT FAILS, SSDINT,NMAX,XL,XR=',G15.7,I5,2G15.7)
      END

      
