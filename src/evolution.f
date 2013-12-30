
      subroutine Evolution

*
* DOALL = true : also evolves the Uplus and Dplus distribution.
* This is not needed to calculate the DIS cross-sections,
* but needed when one store all pdfs.
*
     
      implicit double precision (a-h,o-z)
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'thresholds.inc'
c      common/thresholds/q0,qc,qb

      double precision func0,func1,func24,func22

      external func0            !input parton dists: iparam=0
      external func1            !input parton dists: iparam=1
      external func24           !input parton dists: iparam=24
      external func22           !input parton dists: iparam=22
      external func22text       ! text input
      external func30

      double precision def0, def1, def24, def22,pdfv,glu,glu1,x
      double precision def30


      dimension pdfv(-6:6)
      dimension def22(-6:6,12)    !flavor composition
      dimension def1(-6:6,12)    !flavor composition
      dimension def0(-6:6,12)    !flavor composition
      dimension def24(-6:6,12)    !flavor composition
      dimension def30(-6:6,12)  !flavor composition

      integer iq0, iqfrmq
      double precision eps
cjt test
      
      data nfin/0/
      data q2c/3.D0/, q2b/25.D0/, q2b/200.D0/            !thresh and mu20
       
cjt test 
cv======
cv Remark:
cv need to make it working for h12k parametrisation 
cv ----

      data def0  /                                             ! just a copy of def22
C--       tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--       -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     +     0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0., !dval
     +     0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0., !uval
     +     0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., !s+sbar
     +     0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., !Ubar
     +     0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., ! Dbar
     +     0., 0., 0., -1., 0., 0., 0., 0., 0., 1., 0., 0., 0., !s-sbar
     +     78*0.    /




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
C--       tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--       -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     +     0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0., !dval
     +     0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0., !uval
     +     0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., !s+sbar
     +     0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., !Ubar
     +     0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 0., ! Dbar
     +     0., 0., 0., -1., 0., 0., 0., 0., 0., 1., 0., 0., 0., !s-sbar
     +     78*0.    /



      data def30  /
C--       tb  bb  cb   sb   ub   db   g   d   u   s   c   b   t
C--       -6  -5  -4   -3   -2   -1   0   1   2   3   4   5   6  
     +     0., 0., 0.,  0.,  0.,  1., 0., 1., 0., 0., 0., 0., 0., !d+
     +     0., 0., 0.,  0.,  1.,  0., 0., 0., 1., 0., 0., 0., 0., !u+
     +     0., 0., 0.,  1.,  0.,  0., 0., 0., 0., 1., 0., 0., 0., !s+
     +     0., 0., 0.,  0.,  0., -1., 0., 1., 0., 0., 0., 0., 0., !d-
     +     0., 0., 0.,  0., -1.,  0., 0., 0., 1., 0., 0., 0., 0., !u-
     +     0., 0., 0., -1.,  0.,  0., 0., 0., 0., 1., 0., 0., 0., !s-
     +     78*0.    /



c      call grpars(nx,xmi,xma,nq,qmi,qma,nord)


      q0=starting_scale

c      iqc  = iqfrmq(qc)         !charm threshold
c      iqb  = iqfrmq(qb)         !bottom threshold
c      call setcbt(nfin,iqc,iqb,999) !thesholds in the vfns

      iq0  = iqfrmq(q0)         !starting scale

      if (IPDFSET.eq.5) return  ! for external pdf evolution not needed!

cv ===
      if (PDF_DECOMPOSITION.eq.'LHAPDF')  then
         call evolfg(1,func0,def0,iq0,eps) !evolve all pdf's: LHAPDF

      elseif (PDF_DECOMPOSITION.eq.'QCDNUM_GRID') then
         call evolfg(1,func22text,def22,iq0,eps)

      elseif (Index(PDF_DECOMPOSITION,'D_U_Dbar_Ubar').gt.0) then   ! D,U,Dbar,Ubar 
         call evolfg(1,func1,def1,iq0,eps) !evolve all pdf's: H1

      elseif (Index(PDF_DECOMPOSITION,'Sea').gt.0) then
         call evolfg(1,func24,def24,iq0,eps) !evolve all pdf's: ZEUS

      elseif (PDF_DECOMPOSITION.eq.'Diffractive') then
         call evolfg(1,func30,def30,iq0,eps) !evolve all pdf's: ZEUS diffractive (hard Pomeron)

      elseif (Index(PDF_DECOMPOSITION,'Dbar_Ubar').gt.0) then
         call evolfg(1,func22,def22,iq0,eps)  ! uv, dv, Ubar, Dbar (and also strange)

      else
         print *,'Unknown PDF Decomposition: '//PDF_DECOMPOSITION
         print *,'Stop in evolution'
         call HF_Stop
      endif

      return
      end

*     ----------------------------------------------------
      double precision function func0(id,x)
*     ----------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'steering.inc'

      double precision pdfval, q0
      dimension pdfval(-6:6)

      q0=sqrt(starting_scale)
      
      call evolvePDF(x, q0, pdfval)

      if (id.eq.0) func0=pdfval(0)
      if (id.eq.1) func0=pdfval(1)-pdfval(-1)
      if (id.eq.2) func0=pdfval(2)-pdfval(-2)
      if (id.eq.3) func0=pdfval(3)+pdfval(-3)
      if (id.eq.4) func0=pdfval(-2)+pdfval(-4)
      if (id.eq.5) func0=pdfval(-3)+pdfval(-1)
      if (id.eq.6) func0=pdfval(3)-pdfval(-3)

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
      double precision function func22text(id,x)
*     ----------------------------------------------------
      implicit none
      integer id
      double precision x
      double precision pdf_from_text
C----------------------------
      func22text = pdf_from_text(x,id)

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



*     ----------------------------------------------------
      double precision function func30(id,x)
*     ----------------------------------------------------

      implicit double precision (a-h,o-z)
      include 'pdfparam.inc'
      
      PARAMETER(ParDumpFactor=1.d-3)
      
      dfac = dexp(-ParDumpFactor/(1.00001d0-x))
      func30 = 0.D0
      if (id.eq.0) then
*ws: nchebglu in pdf_param.f must be 0
*ws: initialized to 0 in read_steer.f
         func30 = gluon(x)*dfac
      elseif (id.eq.1.or.id.eq.2.or.id.eq.3) then
*ws: NPOLYVAL in pdf_param.f must be 0
*ws: initialized to 0 in read_steer.f
         func30 = 2*Uval(x)*dfac
      else
      endif
      return
      end

