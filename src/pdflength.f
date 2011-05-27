      Subroutine PDFLength(DeltaChi2)
C
C Created 27 June 2009 by SG
C Add extra constraint for the PDF "length" = int_wmin^wmax \sqrt{1+pdf'(W)**2} dw
C
      implicit none
      
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'pdflength.inc'
      double precision DeltaChi2
      double precision pdflen(5)
      double precision zero
      integer i
      integer ngrid
      parameter(ngrid=500)
      double precision grid(ngrid+1),d0,dst,val,wmin,wmax
C External functions
      double precision glulen,sealen
      double precision ssdint
      external glulen,sealen
      logical LFirstIn
      data LFirstIn /.true./
C----------------------------------------------------
C
C
      if (LFirstIn) then
        Wmin = WMNLEN
        Wmax = WMXLEN
        LFirstIn = .false.
        print '(''PDFLENGTH INITIALIZATION: Set Wmin,Wmax to '',2F8.2)',
     $        Wmin,Wmax
      endif

      DeltaChi2 = 0

      if (pdfLenWeight(1).gt.0) then
         pdflen(1) = ssdint(Wmin,glulen,Wmax)
     $     -(Wmax-Wmin)
         
         DeltaChi2 = DeltaChi2 + pdflen(1)*pdfLenWeight(1)
      else
         pdflen(1) = 0.
      endif

      if (pdfLenWeight(2).gt.0) then
         pdflen(2) = ssdint(Wmin,sealen,Wmax)
     $     -(Wmax-Wmin)
         
         DeltaChi2 = DeltaChi2 + pdflen(2)*pdfLenWeight(2)
      else
         pdflen(2) = 0.
      endif
 

      print *,'Gluon length=',pdflen(1),ag,bg,cg,dg,fg
      print *,'Sea length=',pdflen(2)


C---------------------------------------------------------
      end

      double precision function powerLen(W,a,b,c,d,f)
C
C Utility to calculate pdf length element in W for power parameterization.
C
      implicit none
      double precision W,a,b,c,d,f,q2,x,der,derw,p
C----------------------------------------------------
C Assume Q2=4
      Q2 = 4.D0
      X = Q2/(Q2 + W*W)

      p   = (1.D0+d*x+f*x*x*x)
      der = a*x**b*(1.D0-x)**c*p*
     $     (b/x 
     $     - c/(1.0D0-x) 
     $     + (d+3.D0*f*x*x)/p)
C W derrivative:
      derw = - der* (2*W*Q2)/((W*W+Q2)*(W*W+Q2))

      PowerLen = sqrt(1.D0+derw*derw)      
C----------------------------------------------------
      end

      
      double precision function ChebLen(W,ncheb,poly,a,xminlog,iType)
C
C Utility to calculate pdf length element in W for chebyshev parameterization.
C
      implicit none
      integer ncheb,iType
      double precision W,poly(ncheb),a,xminlog
      double precision Q2,X,XX,Sum,der,derw,sum2
      integer i
      logical LFirst
      data LFirst /.true./
C------------------------------------------------------
      if (LFirst) then
         print *,'First time in ChebLen. IType=',itype
         LFirst = .false.
      endif
C Assume Q2=4
      Q2 = 4.D0
      X = Q2/(Q2 + W*W)
C
C get derrivative:
C
      xx =  (2.D0*log(x)-xminlog)/(-xminlog)
      sum = poly(ncheb)*(ncheb-1)
         
      do i=ncheb-1,2,-1
         sum = sum*xx + poly(i)*(i-1.D0)
      enddo

C SG: Fix for (1-x) dumping:
      if (iType.eq.0) then
C do nothing
      else if (iType.eq.1) then        
C subract term corresponding to -x:
         sum2 = poly(ncheb)*ncheb
         do i=ncheb-1,1,-1
            sum2 = sum2*xx + poly(i)*i
         enddo
         sum = sum - sum2
      endif

      der = sum * a * 2.D0/(-xminlog) / x
C W derrivative:
      derw = - der* (2*W*Q2)/((W*W+Q2)*(W*W+Q2))

      ChebLen = sqrt(1.D0 + derw*derw)
C-------------------------------------------------------
      end
      


      double precision function glulen(W)
      implicit none
      double precision W
      include 'pdfparam.inc'
      include 'steering.inc'
C
      double precision PowerLen,ChebLen
C----------------------------------------------------


      if (nchebglu.eq.0) then
         glulen = powerlen(W,ag,bg,cg,dg,fg)
      else
         glulen = cheblen(W,nchebGlu,polyPars,ag,chebxminlog,
     $        ichebtypeGlu)
      endif
      end


      double precision function Sealen(W)
      implicit none
      double precision W
      include 'pdfparam.inc'
      include 'steering.inc'
C
      double precision PowerLen,ChebLen
C----------------------------------------------------


      if (nchebsea.eq.0) then
         Sealen = powerlen(W,aSea,bSea,cSea,dSea,0.0D0)
      else
         Sealen = cheblen(W,nchebSea,polyParsSea,1.D0,chebxminlog,
     $        ichebtypeSea)

      endif
      end
