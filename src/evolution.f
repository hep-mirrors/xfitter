
      subroutine Evolution

*
* DOALL = true : also evolves the Uplus and Dplus distribution.
* This is not needed to calculate the DIS cross-sections,
* but needed when one store all pdfs.
*
     
      implicit none
c double precision (a-h,o-z)
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
c      common/thresholds/q0,qc,qb

      double precision func1,func24,func22

      external func1            !input parton dists: iparam=1
      external func24           !input parton dists: iparam=24
      external func22           !input parton dists: iparam=22
      
      double precision def1, def24, def22,pdfv,glu,glu1,x

      dimension pdfv(-6:6)
      dimension def22(-6:6,12)    !flavor composition
      dimension def1(-6:6,12)    !flavor composition
      dimension def24(-6:6,12)    !flavor composition

      integer iq0, iqfrmq
      double precision eps
cv======
cv Remark:
cv need to make it working for h12k parametrisation 
cv ----
      data def1  /
C--     tb  bb  cb  sb  ub  db g   d   u   s   c   b   t
C--     -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     +     0., 0., 0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 0., !D
     +     0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 1., 0., 0., !U
     +     0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., !Ubar
     +     0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., !Dbar
     +     0., 0., 0., -1., 0., 0., 0., 0., 0., 1., 0., 0., 0., !s-sbar
     +     0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., !s+sbar
     +     78*0.    /


      data def24  /
cccvC--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
cvC--     -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     +     0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0., !dval
     +     0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0., !uval
     +     0., 2., 2., 2., 2., 2., 0., 0., 0., 0., 0., 0., 0., !sea
     +     0., 0., 0., 0., -1., 1., 0., 0., 0., 0., 0., 0., 0., !delta
     +     0., 0., 0., -1., 0., 0., 0., 0., 0., 1., 0., 0., 0., !s-sbar
     +     0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., !s+sbar
     +     78*0.    /

      data def22  /
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--     -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     +     0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0., !dval
     +     0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0., !uval
     +     0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., !s+sbar
     +     0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., !Ubar
     +     0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., ! Dbar
     +     0., 0., 0., -1., 0., 0., 0., 0., 0., 1., 0., 0., 0., !s-sbar
     +     78*0.    /


c      call grpars(nx,xmi,xma,nq,qmi,qma,nord)


      q0=starting_scale

c      iqc  = iqfrmq(qc)         !charm threshold
c      iqb  = iqfrmq(qb)         !bottom threshold
c      call setcbt(nfin,iqc,iqb,999) !thesholds in the vfns

      iq0  = iqfrmq(q0)         !starting scale

cv ===
      if (iparam.eq.1)  call evolfg(1,func1,def1,iq0,eps) !evolve all pdf's: H1
      if (iparam.eq.4)  call evolfg(1,func24,def24,iq0,eps) !evolve all pdf's: ZEUS
cv ===

      if (iparam.eq.2011) call evolfg(1,func22,def22,iq0,eps)

      if (iparam.eq.222222) call evolfg(1,func22,def22,iq0,eps)
      if (iparam.eq.222223) call evolfg(1,func22,def22,iq0,eps)
      if (iparam.eq.22)  call evolfg(1,func22,def22,iq0,eps) !evolve all pdf's: H1ZEUS
      if (iparam.eq.171717)  call evolfg(1,func22,def22,iq0,eps) !evolve all pdf's: H1ZEUS
      if (iparam.eq.225)  call evolfg(1,func22,def22,iq0,eps) !evolve all pdf's: H1ZEUS
      if (iparam.eq.222)  call evolfg(1,func22,def22,iq0,eps) !evolve all pdf's: H1ZEUS
      if (iparam.eq.229)  call evolfg(1,func22,def22,iq0,eps) !evolve all pdf's: H1ZEUS
      if (iparam.eq.221)  call evolfg(1,func22,def22,iq0,eps) !evolve all pdf's: H1ZEUS

c      glu1=pdfval(0,0.000103523178d0,1.95,0)
cv      call allpdf(0.000103523178d0,1.95d0,pdfv,0)


      return
      end

*     ----------------------------------------------------
      double precision function func1(id,x)
*     ----------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'pdfparam.inc'

      if (id.eq.0) func1=gluon(x)
      if (id.eq.1) func1=H1D(x)
      if (id.eq.2) func1=H1U(x)
      if (id.eq.3) func1=Ubar(x)
      if (id.eq.4) func1=Dbar(x)
      if (id.eq.6) func1=2*qstrange(x)
      if (id.eq.5) func1=0.d0


      return
      end
      

*     ----------------------------------------------------
      double precision function func22(id,x)
*     ----------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'pdfparam.inc'
      include 'steering.inc'

      func22 = 0.D0
      if (id.eq.0) func22=gluon(x)
      if (id.eq.1) func22=dval(x)
      if (id.eq.2) func22=uval(x)
      
      if (id.eq.3) func22=2*qstrange(x)
      
      if (id.eq.4) func22=ubar(x)
      if (id.eq.5) func22=dbar(x)
      if (id.eq.6) func22=0.d0

      return
      end

*     ----------------------------------------------------
      double precision function func24(id,x)
*     ----------------------------------------------------

      implicit double precision (a-h,o-z)
      include 'pdfparam.inc'


      if (id.eq.0) func24=gluon(x)
      if (id.eq.1) func24=dval(x)
      if (id.eq.2) func24=uval(x)
      if (id.eq.3) func24=sea(x)
      if (id.eq.4) func24=dbmub(x)
      if (id.eq.5) func24=0.d0
      if (id.eq.6) func24=0.d0

      return
      end
