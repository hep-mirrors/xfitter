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
      double precision RFROMF,ASFUNC
      integer nfl, ialerr
      Q2 = Q**2
      mur2 = RFROMF(Q2) ! get renorm scale from factscale in QCDNUM


      ialerr = 0
c     get alpha_s interpolation
      appl_fnalphas = ASFUNC(mur2, nfl, ialerr)
      if ( ialerr .eq. 1) then
        print *, 'alpha_s convolution error. Renorm scale too low:'
	print *, 'rens^2 = ',mur2, ' < 0.1GeV^2'
	print *, 'setting alpha_s at rens^2 = 0.1GeV^2'
	mur2 = 0.1
	appl_fnalphas = ASFUNC(mur2, nfl, ialerr)
      endif

      return
      end

      subroutine appl_fnpdf(x, Q, xf)
      implicit none
C------------------------------------
      include 'applgrid_fastpdf.inc'
      include 'fcn.inc'
      include 'steering.inc'
      double precision x, Q, Q2
      double precision xf(-6:6)
      integer iqnset, iqnchk,ifl
      logical LFirstTime
      data LFirstTime/.true./
      
C---------------------------------------
      if (lFirstTime) then
C Reset NAPPLPDF
         LFirstTime = .false.
         NAPPLPDF = 0
      endif

      iqnset = 1
      iqnchk = 0
      do ifl=-6,6
        xf(ifl)=0.d0
      enddo
      if ( x .lt. 1.d-7 .or. x .gt. 1d0-1d-7 ) return

      Q2 = Q*Q

      if (IFlagFCN.eq.1 .and. LFastAPPLGRID) then
         NAPPLPDF = NAPPLPDF + 1
         if (NAPPLPDF.gt.NAPPLPDFMAX) then
            print *,'Increase NAPPLPDFMAX in applgrid_fastpdf'
            stop
         endif
         XAPPLPDF(NAPPLPDF) = x
         QAPPLPDF(NAPPLPDF) = q2
      endif



      if (LFastAPPLGRID) then
         if (IFlagFCN.eq.1) then
            call fpdfxq(iqnset, x, Q2, xf, iqnchk)
         else
            call Retrive_pdf_applgrid_fast(xf)
         endif
      else
         call fpdfxq(iqnset, x, Q2, xf, iqnchk)
      endif

      return
      end

C Fast code to cacl. PDFs for applgrid
      Subroutine Calc_pdf_applgrid_fast
      implicit none
      include 'applgrid_fastpdf.inc'
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
      call fastini(XAPPLPDF,QAPPLPDF,NAPPLPDF,ICHECK)
      
      iset = 1
      isel = 7
      ibuf = -1
      inbuf = -ibuf
      do i=2,12
         if (i.eq.7) then
C            call fastsns(iset,fdef(1,i),0,inbuf)
            call fastepm(iset,0,ibuf)
         else
            call fastsns(iset,fdef(1,i),isel,ibuf)
         endif
c         print *,'ho',ibuf,inbuf
         call fastfxq(inbuf,APPLPDF(1,i-7),NAPPLPDF)
      enddo

c      print *,'xxx',xapplpdf(1),qapplpdf(1)
c      print '(''hahaha''13F10.2)',(applpdf(1,i),i=-6,6)
      end

      subroutine Retrive_pdf_applgrid_fast(PDFs)
C
C Retrieve next PDF value from the buffer
C
      implicit none
      double precision pdfs(-6:6)
      include 'applgrid_fastpdf.inc'
      integer IRound,i
      data IRound/0/
C----------------------------------------
      if (IRound .eq. NAPPLPDF) then
C Reached the last, start over
         IRound = 0
      endif
      IRound = IRound + 1
      do i=-6,6
         PDFs(i) = APPLPDF(IRound,i)
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

      subroutine ag_convolute(igrid, xsec)
C-----------------
      implicit none
C---------
      double precision xsec(100)
      integer igrid

      call appl_convolute(igrid,xsec)
      return 
      end

      subroutine ag_releasegrids
      implicit none
      call appl_releasegrids
      return 
      end
