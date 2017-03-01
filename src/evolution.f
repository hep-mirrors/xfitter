      subroutine Evolution

*
* DOALL = true : also evolves the Uplus and Dplus distribution.
* This is not needed to calculate the DIS cross-sections,
* but needed when one store all pdfs.
*
     
c      implicit double precision (a-h,o-z)
      implicit none
#include "steering.inc"
#include "pdfparam.inc"
#include "thresholds.inc"
#include "alphas.inc"
#include "couplings.inc"
#include "ntot.inc"
#include "datasets.inc"
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
      integer nfin
      double precision q2c,q2b

      data nfin/0/
      data q2c/3.D0/, q2b/25.D0/, q2b/200.D0/            !thresh and mu20

*     ---------------------------------------------------------
*     declaration related to alphas
*     for RT code, transfer alpha S
*     --------------------------------------------------------- 
      double precision alphaszero

*     ---------------------------------------------------------
*     Subroutine of the external evolution codes
*     --------------------------------------------------------- 
      external LHAPDFsubr
      external APFELsubr
      external APFELsubrPhoton
      external QEDEVOLsubr

      double precision epsi
      double precision hf_get_alphas,asRef

*     ---------------------------------------------------------
*     Save scale in a common to avoid APFEL to evolve if the
*     scale does not change.
*     --------------------------------------------------------- 
      double precision q2p
      common / PrevoiusQ / q2p

      integer NextraSets

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

      integer i
*
*     Set parameters of the initial scale PDFs to be used by MELA
*
      if (mod(hfscheme,10).eq.6) then
         call SetxFitterParametersMELA(parubar,pardbar,
     1                                 paruval,pardval,
     2                                 parglue,
     3                                 fstrange,fcharm)
      endif

c      call grpars(nx,xmi,xma,nq,qmi,qma,nord)

      q0=starting_scale

c      iqc  = iqfrmq(qc)         !charm threshold
c      iqb  = iqfrmq(qb)         !bottom threshold
c      call setcbt(nfin,iqc,iqb,999) !thesholds in the vfns

      iq0  = iqfrmq(q0)         !starting scale

*
*     Initialize alphas
*
      if (itheory.eq.0.or.itheory.eq.11) then
         call setalf(dble(alphas),Mz*Mz)
      else
*     Make sure that alphas is correctly set at the Z mass
         call SetMaxFlavourPDFs(5)
         call SetMaxFlavourAlpha(5)
         call SetAlphaQCDRef(dble(alphas),dble(Mz))
         asref = HF_Get_alphas(dble(starting_scale))
         call SetAlphaQCDRef(asref,dsqrt(dble(starting_scale)))
      endif
      alphaSzero = hf_get_alphas(1D0)
      call RT_SetAlphaS(alphaSzero)

      NextraSets = 0
      if (ExtraPdfs) then
         NextraSets = 1
      endif

C ---- LHAPDF ----
      if(IPDFSET.eq.5) then
         call PDFEXT(LHAPDFsubr,IPDFSET,NextraSets,dble(0.001),epsi)
         return
C ---- APFEL ----
      elseif (IPDFSET.eq.7) then
         if(itheory.eq.35)then
            q2p = starting_scale
            call SetPDFSet("external")
            call PDFEXT(APFELsubrPhoton,IPDFSET,1,dble(0.001),epsi)
         else
*     6 flavours (default)
            q2p = starting_scale
            call SetPDFSet("external")
            call SetMaxFlavourPDFs(6)
            call SetMaxFlavourAlpha(6)
            call PDFEXT(APFELsubr,IPDFSET,0,dble(0.001),epsi)
C
C     If the use of the H-VFNS is required, enable the evolutions
C     with different values of NFmax
C
            if(UseHVFNS)then
*     5 flavours
               q2p = starting_scale
               call SetPDFSet("external")
               call SetMaxFlavourPDFs(5)
               call SetMaxFlavourAlpha(5)
               call PDFEXT(APFELsubr,IPDFSET+1,0,dble(0.001),epsi)
*     4 flavours (redefine alphas)
               q2p   = starting_scale
               call SetPDFSet("external")
               call SetMaxFlavourPDFs(4)
               call SetMaxFlavourAlpha(4)
               call PDFEXT(APFELsubr,IPDFSET+2,0,dble(0.001),epsi)
*     3 flavours (redefine alphas)
               q2p   = starting_scale
               call SetPDFSet("external")
               call SetMaxFlavourPDFs(3)
               call SetMaxFlavourAlpha(3)
               call PDFEXT(APFELsubr,IPDFSET+3,0,dble(0.001),epsi)
            endif
         endif
         return
C ---- QEDEVOL ----
      elseif (IPDFSET.eq.8) then
         call qedevol_main
CC         call PDFEXT(QEDEVOLsubr,IPDFSET,1,dble(0.001),epsi)
         return
      endif

cv ===
      if (PDF_DECOMPOSITION.eq.'LHAPDF')  then
         call evolfg(1,func0,def0,iq0,eps) !evolve all pdf's: LHAPDF

      elseif (PDF_DECOMPOSITION.eq.'QCDNUM_GRID') then
         call evolfg(1,func22text,def22,iq0,eps)

      elseif (Index(PDF_DECOMPOSITION,'D_U_Dbar_Ubar').gt.0) then ! D,U,Dbar,Ubar 
         call evolfg(1,func1,def1,iq0,eps) !evolve all pdf's: H1

      elseif (Index(PDF_DECOMPOSITION,'Sea').gt.0) then
         call evolfg(1,func24,def24,iq0,eps) !evolve all pdf's: ZEUS

      elseif (PDF_DECOMPOSITION.eq.'Diffractive') then
         call evolfg(1,func30,def30,iq0,eps) !evolve all pdf's: ZEUS diffractive (hard Pomeron)

      elseif (Index(PDF_DECOMPOSITION,'Dbar_Ubar').gt.0) then
         call evolfg(1,func22,def22,iq0,eps) ! uv, dv, Ubar, Dbar (and also strange)

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
#include "steering.inc"

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
#include "pdfparam.inc"

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
#include "pdfparam.inc"
#include "steering.inc"

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
#include "pdfparam.inc"


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
#include "pdfparam.inc"
      
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
*
************************************************************************
*
*     Subroutine that defines the PDFs to be evolved with APFEL.
*     (predefined name)
*
************************************************************************
      subroutine ExternalSetAPFEL(x,q0,xf)
*
      implicit none
#include "steering.inc"
**
*     Input Variables
*
      double precision x
      double precision q0
**
*     Internal Variables
*
      integer ipdf
      double precision gluon
      double precision pdf_from_text
      double precision qstrange,Ubar,Dbar,H1U,H1D
      double precision sea,dbmub,dval,uval
      double precision photon
      double precision dfac,ParDumpFactor
      parameter(ParDumpFactor=1.d-3)
**
*     Output Variables
*
      double precision xf(-6:7)
*
*     Set PDFs to zero
*
      do ipdf=-6,7
         xf(ipdf) = 0d0
      enddo
      if(x.gt.1d0) x = 1d0
*
*     Construct PDFs addording to the PDF decomposition
*
      if(PDF_DECOMPOSITION.eq.'LHAPDF')then
c         q0 = sqrt(starting_scale)
         call evolvePDF(x, q0, xf)

      elseif(PDF_DECOMPOSITION.eq.'QCDNUM_GRID')then
         xf(-3) = ( pdf_from_text(x,3) - pdf_from_text(x,6) ) / 2d0
         xf(-2) = pdf_from_text(x,4)
         xf(-1) = pdf_from_text(x,5)
         xf(0)  = pdf_from_text(x,0)
         xf(1)  = pdf_from_text(x,1) - pdf_from_text(x,5)
         xf(2)  = pdf_from_text(x,2) - pdf_from_text(x,4)
         xf(3)  = ( pdf_from_text(x,3) + pdf_from_text(x,6) ) / 2d0

      elseif(Index(PDF_DECOMPOSITION,'D_U_Dbar_Ubar').gt.0)then ! D,U,Dbar,Ubar 
         xf(-3) = qstrange(x)
         xf(-2) = Ubar(x)
         xf(-1) = Dbar(x)
         xf(0)  = gluon(x)
         xf(1)  = H1D(x) - xf(-3)
         xf(2)  = H1U(x)
         xf(3)  = xf(-3)

      elseif(Index(PDF_DECOMPOSITION,'Sea').gt.0)then
         xf(-2) = sea(x) / 4d0 - dbmub(x) / 2d0
         xf(-1) = sea(x) / 4d0 + dbmub(x) / 2d0
         xf(0)  = gluon(x)
         xf(1)  = dval(x) + xf(-1)
         xf(2)  = uval(x) + xf(-2)

      elseif(PDF_DECOMPOSITION.eq.'Diffractive')then
         dfac = dexp(-ParDumpFactor/(1.00001d0-x))
*
         xf(-3) = dfac * Uval(x)
         xf(-2) = xf(-3)
         xf(-1) = xf(-3)
         xf(0)  = dfac * gluon(x)
         xf(1)  = xf(-3)
         xf(2)  = xf(-3)
         xf(3)  = xf(-3)

      elseif(Index(PDF_DECOMPOSITION,'Dbar_Ubar').gt.0)then
         xf(-3) = qstrange(x)
         xf(-2) = ubar(x)
         xf(-1) = dbar(x) - xf(-3)
         xf(0)  = gluon(x)
         xf(1)  = dval(x) + xf(-1)
         xf(2)  = uval(x) + xf(-2)
         xf(3)  = xf(-3)

      else
         print *,'Unknown PDF Decomposition: '//PDF_DECOMPOSITION
         print *,'Stop in evolution'
         call HF_Stop
      endif
*
*     Photon PDF
*
      if(itheory.eq.35) xf(7) = photon(x)
*
      return
      end

      Subroutine LHAPDFsubr(x, qmu2, xf)
C-------------------------------------------------------
C
C External PDF reading for QCDNUM
C
C--------------------------------------------------------
      implicit none
#include "steering.inc"

      double precision x,qmu2
      double precision xf(-6:7)
      if ( ExtraPdfs ) then
         call evolvePDFphoton(x, sqrt(qmu2), xf, xf(7))
      else
         call evolvePDF(x, sqrt(qmu2), xf)
      endif      
      end

      Subroutine APFELsubr(x, qmu2, xf)
C-------------------------------------------------------
C
C External PDF reading for APFEL
C
C--------------------------------------------------------
      implicit none
*
#include "steering.inc"
*
      integer i
      double precision x,qmu2
      double precision xf(-6:6)

      double precision q2p
      common / PrevoiusQ / q2p
*
*     Perform evolution with APFEL only if the final scale has changed
*
      if(qmu2.ne.q2p)then
         call EvolveAPFEL(dsqrt(q2p),dsqrt(qmu2))
         call SetPDFSet("apfel")
      endif
*
      call xPDFall(x,xf)
      q2p = qmu2
*
      return
      end
*
      Subroutine APFELsubrPhoton(x, qmu2, xf)
C-------------------------------------------------------
C
C External PDF reading for APFEL (including the photon)
C
C--------------------------------------------------------
      implicit none
*
#include "steering.inc"
*
      double precision x,qmu2
      double precision xf(-6:7)

      double precision q2p
      common / PrevoiusQ / q2p
*
*     Perform evolution with APFEL only if the final scale has changed
*
      if(qmu2.ne.q2p)then
c         call EvolveAPFEL(dsqrt(q2p),dsqrt(qmu2))
c         call SetPDFSet("apfel")
         call EvolveAPFEL(dsqrt(dble(starting_scale)),dsqrt(qmu2))
      endif
*
      call xPDFallPhoton(x,xf)
      q2p = qmu2
*
      return
      end
