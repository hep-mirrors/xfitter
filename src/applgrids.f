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
      include 'fcn.inc'
      include 'steering.inc'
      double precision x, Q, Q2
      double precision xf(-6:6)
      integer iqnset, iqnchk,ifl
      
C---------------------------------------

      iqnset = IPDFSET
      iqnchk = 0
      do ifl=-6,6
        xf(ifl)=0.d0
      enddo
      if ( x .lt. 1.d-7 .or. x .gt. 1d0-1d-7 ) return

      Q2 = Q*Q

      if (IFlagFCN.eq.1 .and. LFastAPPLGRID) then

         Call Register_pdf_applgrid(x,q2)

      endif



      if (LFastAPPLGRID) then
         if (IFlagFCN.eq.1) then
            call hf_get_pdfs(x, Q2, xf)
         else
            call Retrive_pdf_applgrid_fast(xf)
         endif
      else
         call hf_get_pdfs( x, Q2, xf)
      endif

      return
      end

C> @brief Register x,Q2 point on a grid
      Subroutine register_pdf_applgrid(x,Q2)
C----------------------------------------------
      implicit none
      double precision x,Q2
      include 'applgrid_fastpdf.inc'
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
      include 'applgrid_fastpdf.inc'
      include 'steering.inc'
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

      do i=1,13
         call mypdflist(iset,fdef(1,i),XAPPLPDF,QAPPLPDF,
     $        APPLPDF(1,i-7),NAPPLPDFINT,ICHECK) 
      enddo

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

C     ===========================================      
      subroutine mypdflist(iset,def,x,q,f,n,ichk)
C     ===========================================

C--   Interpolation of linear combination of pdfs using fast engine.
C--   The number of interpolations can be larger than mpt0 since
C--   this routine autmatically buffers in chunks of mpt0 words.
C--
C--   iset       (in)   pdf set [1-9]
C--   def(-6:6)  (in)   coefficients, ..., sb, ub, db, g, d, u, s, ...
C--                     for gluon  set def(0) = non-zero
C--                     for quarks set def(0) = 0.D0   
C--   x          (in)   list of x-points
C--   q          (in)   list of mu2 points
C--   f          (out)  list of interpolated pdfs
C--   n          (in)   number of items in the list
C--   ichk       (in)   0/1  no/yes check grid boundary  
      
      implicit double precision (a-h,o-z)      
      
      parameter(nmax = 5000)    !should be set to mpt0, or less
      
      dimension xx(nmax),qq(nmax)
      dimension def(-6:6), coef(0:12,3:6)
      dimension x(*), q(*), f(*)
      
C--   Fatal error that cannot be caught by QCDNUM
      if(n.le.0) stop 'MYPDFLIST: n.le.0 --> STOP'      
      
      call setUmsg('MYPDFLIST') !QCDNUM will catch all further errors

C--   Translate the coefficients def(-6:6) from flavour space to 
C--   singlet/non-siglet space. These coefs are n_f dependent.
       
      if(def(0).eq.0.D0) then                !quarks
        do nf = 3,6
          coef(0,nf) = 0.D0                  
          call efromqq(def, coef(1,nf), nf) 
        enddo
      else                                   !gluon
        do nf = 3,6
          coef(0,nf) = def(0)
          do i = 1,12
            coef(i,nf) = 0.D0 
          enddo  
        enddo
      endif  
      
C--   Fill output array f in batches of nmax words
      ipt = 0
      jj  = 0      

      qmax =0
      qmin =10000
      xmax =0
      xmin =1
      do i = 1,n
        ipt     = ipt+1
        xx(ipt) = x(i)
        qq(ipt) = q(i)

        xmin = min(x(i),xmin)
        xmax = max(x(i),xmax)
        qmin = min(q(i),qmin)
        qmax = max(q(i),qmax)

        if(ipt.eq.nmax) then
c           print *,xmin,xmax,qmin,qmax

          call fastini(xx,qq,nmax,ichk)
          call fastsum(iset,coef,-1)         !fill sparse buffer 1
          call fastfxq(1,f(jj*nmax+1),nmax)
          ipt = 0
          jj  = jj+1
        endif
      enddo
C--   Flush remaining ipt points
      if(ipt.ne.0) then
        call fastini(xx,qq,ipt,ichk)
        call fastsum(iset,coef,-1)           !fill sparse buffer 1
        call fastfxq(1,f(jj*nmax+1),ipt)          
      endif
        
      call clrUmsg
      
      return
      end
      
