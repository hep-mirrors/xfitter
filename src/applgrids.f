c----------------------------------------------------------
c     routines that call  the pdf and alphas routines
c     from QCDNUM
c----------------------------------------------------------
      double precision function appl_fnalphas(Q)
C-     
      implicit none
C-
      double precision Q, Q2      ! fact scale
      double precision mur2       ! renorm scale
      double precision alphaspdf
      double precision RFROMF,hf_get_alphas
      integer nfl, ialerr
      Q2 = Q**2
      mur2 = RFROMF(Q2) ! get renorm scale from factscale in QCDNUM


      ialerr = 0
      
      appl_fnalphas = HF_Get_alphas(mur2)

      return
      end

      subroutine appl_fnpdf(x, Q, xf)
C
C> @brief Interface to QCDNUM PDF call. 
C
      implicit none
C------------------------------------
#include "fcn.inc"
#include "steering.inc"
      double precision x, Q, Q2
      double precision xf(-N_CHARGE_PDF:N_CHARGE_PDF)
      double precision xfe(-N_CHARGE_PDF:N_CHARGE_PDF+N_NEUTRAL_PDF)
      integer iqnset, iqnchk,ifl
      
C---------------------------------------
      iqnset = IPDFSET
      iqnchk = 0
      do ifl=-N_CHARGE_PDF,N_CHARGE_PDF
        xf(ifl)=0.d0
        xfe(ifl)=0.d0
      enddo
      xfe(N_CHARGE_PDF+N_NEUTRAL_PDF)=0.d0

c     Return zero if x range falls below qcdnum grid xmin values (to avoid large weights)
      if (PDFStyle.ne.'LHAPDFNATIVE') then
         if ( x .lt. xmin_grid(1) .or. x .gt. 1d0-1d-7 ) then 
            Call HF_ERRLOG(14060218,
     $'W: x value below xmin in qcdnum grid, applgrid weight set to 0')
            return
         endif
      endif

      Q2 = Q*Q

      if (IFlagFCN.eq.1 .and. LFastAPPLGRID) then

         Call Register_pdf_applgrid(x,q2)

      endif


      if (LFastAPPLGRID) then
         if (IFlagFCN.eq.1) then
            call hf_get_pdfs(x, Q2, xfe)
            do ifl=-N_CHARGE_PDF,N_CHARGE_PDF
              xf(ifl)=xfe(ifl)
            enddo
         else
            call Retrive_pdf_applgrid_fast(xf)
         endif
      else
         call hf_get_pdfs( x, Q2, xfe)
         do ifl=-N_CHARGE_PDF,N_CHARGE_PDF
           xf(ifl)=xfe(ifl)
         enddo
      endif

      return
      end

C> @brief PDF for ppbar process
      subroutine appl_fnpdf_bar(x, Q, xf)
      implicit none
#include "steering.inc"
      double precision x, Q, xft
      double precision xf(-N_CHARGE_PDF:N_CHARGE_PDF+N_NEUTRAL_PDF)
      integer ifl

      call appl_fnpdf(x, Q, xf)
      do ifl=0,6
        xft=xf(ifl)
        xf(ifl) = xf(-ifl)
        xf(-ifl)=xft
      enddo

      return
      end

C> @brief PDF for pn process
      subroutine appl_fnpdf_neut(x, Q, xf)
      implicit none
#include "steering.inc"
      double precision x, Q, xft
      double precision xf(-N_CHARGE_PDF:N_CHARGE_PDF+N_NEUTRAL_PDF)
      integer ifl

      call appl_fnpdf(x, Q, xf)
      ! switch up and down
      xft=xf(1)
      xf(1) = xf(2)
      xf(2) = xft
      ! switch anti-up and anti-down
      xft=xf(-1)
      xf(-1) = xf(-2)
      xf(-2) = xft

      return
      end

C> @brief Register x,Q2 point on a grid
      Subroutine register_pdf_applgrid(x,Q2)
C----------------------------------------------
      implicit none
      double precision x,Q2
#include "applgrid_fastpdf.inc"
      logical LFirstTime
      data LFirstTime/.true./
      integer ILoc
C      
      if (lFirstTime) then
C Reset NAPPLPDF
         LFirstTime = .false.
         NAPPLPDF = 0
         NAPPLPDFINT = 0
         print *,'Cashing PDF calls for APPLGRID, this may take a while'
      endif
C

C try to locate X,Q2 point
      Call locateXQ2(XAPPLPDF,QAPPLPDF,NAPPLPDF,X,Q2,ILoc)


C increase ref.
      NAPPLPDF = NAPPLPDF + 1

      if (mod(NAPPLPDF,10000).eq.0) then
         print '(''... cashing '',I7,'' calls to '',I7,'' calls'')',
     $        NAPPLPDF,NAPPLPDFINT
      endif

      if (ILoc.gt.0) then
         IRefApp(NAPPLPDF) = ILoc
      else

         NAPPLPDFINT = NAPPLPDFINT + 1
         IRefApp(NAPPLPDF) = NAPPLPDFINT
         if (NAPPLPDFINT.gt.NAPPLPDFMAX) then
            print *,'Increase NAPPLPDFMAX in applgrid_fastpdf'
            stop
         endif

      
         XAPPLPDF(NAPPLPDFINT) = x
         QAPPLPDF(NAPPLPDFINT) = q2
      endif
      end


C Fast code to cacl. PDFs for applgrid
      Subroutine Calc_pdf_applgrid_fast
      implicit none
#include "applgrid_fastpdf.inc"
#include "steering.inc"
      integer ICheck,i,ibuf,inbuf,iset,isel

      double precision fdef(13,13)
      data fdef/
     $     1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., ! tb
     $     0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., ! bb
     $     0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., ! cb
     $     0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., ! sb
     $     0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., ! ub
     $     0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., ! db
     $     0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., ! g
     $     0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., ! d
     $     0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., ! u
     $     0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., ! s
     $     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., ! c
     $     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., ! b
     $     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1. ! t
     $       /
C----------------------------------------
      if (NAPPLPDF.eq.0) Return
 
      ICheck = 1 ! force QCDNUM checks
      iset = IPDFSET


      ! First call gluon:
      call fflist(iset,fdef(1,7),0,XAPPLPDF,QAPPLPDF, 
     $     APPLPDF(1,0),NAPPLPDFINT,ICHECK) 
      do i=1,13
         if (i.ne.7) then
            call fflist(iset,fdef(1,i),1,XAPPLPDF,QAPPLPDF, 
     $           APPLPDF(1,i-7),NAPPLPDFINT,0)  ! no check          
         endif
      enddo

      end

      subroutine Retrive_pdf_applgrid_fast(PDFs)
C
C Retrieve next PDF value from the buffer
C
      implicit none
#include "steering.inc"
      double precision pdfs(-N_CHARGE_PDF:N_CHARGE_PDF+N_NEUTRAL_PDF)
#include "applgrid_fastpdf.inc"
      integer IRound,i
      data IRound/0/
C----------------------------------------
      if (IRound .eq. NAPPLPDF) then
C Reached the last, start over
         IRound = 0
      endif
      IRound = IRound + 1
      do i=-N_CHARGE_PDF,N_CHARGE_PDF+N_NEUTRAL_PDF
         PDFs(i) = APPLPDF(IRefApp(IRound),i)
      enddo
      end

c----------------------------------------------------------

      subroutine ag_ngrids( ng )
c      integer ng
      call appl_ngrids(ng)

      return 
      end


      subroutine ag_gridids(igrids)
      implicit none
      integer igrids(100)
      call appl_gridids(igrids)

      return 
      end

      integer function ag_getnbins(igrid)
      implicit none
C--
      integer igrid
      integer appl_getnbins
      ag_getnbins = appl_getnbins(igrid)
      return 
      end

      subroutine ag_convolute(igrid, iorder, muR, muF, xsec)
C-----------------
      implicit none
C---------
      double precision xsec(100)
      integer igrid, iorder
      double precision muR, muF

      call appl_convoluteorder(igrid,iorder-1,muR,muF,xsec)
      return 
      end

      subroutine ag_releasegrids
      implicit none
      call appl_releasegrids
      return 
      end



      Subroutine  locateXQ2(XArr,QArr,NAPPLPDF,X,Q,ILoc)
C-------------------------------------------------------------
C
C  Locate X,Q point in  XArr,QArr array
C
C-------------------------------------------------------------
      implicit none
      integer  ILoc,NAPPLPDF
      double precision QArr(NAPPLPDF), XArr(NAPPLPDF),X,Q
      integer I
      double precision epsilonX
      double precision epsilonQ
      parameter (epsilonX = 1.0D-8)
      parameter (epsilonQ = 1.0D-8)

C---------------------------------------------------------

      do I=1,NAPPLPDF
         if ( ( abs(QArr(I)-Q)/Q.lt.epsilonQ) 
     $        .and. (abs (XArr(I)-X)/X.lt. epsilonX)) then
            ILoc = I
            Return
         endif
      enddo
      
      ILoc = -1

      end
*
************************************************************************
*
*     Subroutine that returns the PDFs at the initial scale in the
*     physical basis. Needed by APFELgrid
*
************************************************************************
      subroutine apfel_fnpdf(x, Q0, xf)
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
      double precision xf(-6:6)
*
*     Set PDFs to zero
*
      do ipdf=4,6
         xf(ipdf)  = 0d0
         xf(-ipdf) = 0d0
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
      return
      end
