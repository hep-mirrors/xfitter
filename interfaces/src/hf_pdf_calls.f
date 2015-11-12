      subroutine HF_Get_PDFs(x,q2,PDFSF)
C----------------------------------------------------------------------
C Interface to PDF 
C
C  Input:  x, Q2 values
C  Output: 13 PDF values
C----------------------------------------------------------------------
      implicit none
C----------------------------------------------------------------------
#include "steering.inc"
#include "fcn.inc"

      double precision x,q2
     $     ,pdfsf(-N_CHARGE_PDF:N_CHARGE_PDF+N_NEUTRAL_PDF)

      integer i
      double precision A,Z, tmpU,tmpD,tmpUb,tmpDb
      data A,Z /207,82/

      double precision FSNSXQ
C----------------------------------------------------------------------
      if (PDFStyle.eq.'LHAPDFNATIVE') then
         if (ExtraPdfs) then
C photon is present !
            call evolvePDFphoton(x, sqrt(q2), pdfsf, pdfsf(7))
         else
            call evolvePDF(x, sqrt(q2), pdfsf)
         endif
         return
      endif

      if ( ExtraPdfs ) then
C QED evolution:
         call FPDFXQ(iPDFSET,x,q2,PDFSF,ICheck_QCDNUM)         
         PDFSF(N_CHARGE_PDF+N_NEUTRAL_PDF) = FSNSXQ( iPDFSET,13,x,q2,ICheck_QCDNUM)
c         print *,ipdfset,x,q2, PDFSF(N_CHARGE_PDF+1), PDFSF(0)
         return
      endif

      if (x.ge.1.D0) then
         do i=-6,6
            PDFSF(i) = 0
         enddo

         if (ICheck_QCDNUM.gt.0) then
            Call HF_errlog(12042201,'W:HF_GET_PDFS X value >= 1.0') 
         endif
         Return
      endif

      if ( CachePDFs ) then
C Cache PDF calls:
         call GetCachedPDFs(ifcncount,x,q2,PDFSF)
      else
C Get PDFs directly:
         Call HF_Get_PDFs_UnCached(x,q2,PDFSF)
      endif
  
C---- switch for lead PDF: Combine to form nuclear pdf; scale by A
C----     For full cross setion on lead, multiply by A

      if(lead) then
         if (deuteron) then
            A=2
            Z=1
         endif
         tmpU  = (Z*PDFSF( 1) + (A-Z)*PDFSF( 2) )/A
         tmpD  = (Z*PDFSF( 2) + (A-Z)*PDFSF( 1) )/A
         tmpUb = (Z*PDFSF(-1) + (A-Z)*PDFSF(-2) )/A
         tmpDb = (Z*PDFSF(-2) + (A-Z)*PDFSF(-1) )/A
         PDFSF( 1) = tmpU
         PDFSF( 2) = tmpD
         PDFSF(-1) = tmpUb
         PDFSF(-2) = tmpDb

          

      endif


C----------------------------------------------------------------------
      end


      subroutine HF_PDFFAST(pdfdef,Xarray,Q2array,PDFout,Npoints)
C----------------------------------------------------------------------
C
C Interface to QCDNUM fast PDF calls
C
C----------------------------------------------------------------------
      implicit none
C----------------------------------------------------------------------
#include "steering.inc"
      integer NPoints
      double precision pdfdef(*),Xarray(*),Q2array(*),PDFout(*)
      integer icheck
      data icheck/1/  ! Force PDF checks
C----------------------------------------------------------------------
      call mypdflist
     $     (IPDFSET,pdfdef,Xarray,Q2ARRAY,PDFout,Npoints,Icheck)
      end

      double precision Function HF_Get_alphas(q2)
C----------------------------------------------------
C
C  Return alpha_S value for a given Q2
C
C----------------------------------------------------
      implicit none
#include "steering.inc" 
      double precision Q2
      integer nf,ierr
      double precision ASFUNC,alphaQCD
      double precision alphaspdf

C----------------------------------------------------
      if (PDFStyle.eq.'LHAPDFNATIVE') then
         HF_Get_alphas = alphaspdf(dsqrt(q2))
      else
      if (itheory.eq.10) then
         HF_Get_alphas = alphaQCD(dsqrt(q2))
      else
         HF_Get_alphas = ASFUNC(q2,nf,ierr) 
      endif
      endif
C----------------------------------------------------
      end


      double precision function hf_get_mur(iDataSet)
C------------------------------------------------------
C
C Get data-set depdendent renormalisation scale
C
C------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "scales.inc"
      integer iDataSet
C------------------------------------------------------

      hf_get_mur = DataSetMur(iDataSet)

      end


      double precision function hf_get_muf(iDataSet)
C------------------------------------------------------
C
C Get data-set depdendent renormalisation scale
C
C------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "scales.inc"
      integer iDataSet
C------------------------------------------------------

      hf_get_muf = DataSetMuf(iDataSet)

      end

Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

C----------------------------------------------------------------------
C-------------  Internal ----------------------------------------------
C----------------------------------------------------------------------


Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


      subroutine HF_Get_PDFs_UnCached(x,q2,PDFSF)
C----------------------------------------------------------------------
C Interface to PDF 
C
C  Input:  x, Q2 values
C  Output: 13 PDF values
C----------------------------------------------------------------------
      implicit none
C----------------------------------------------------------------------
#include "steering.inc"
      double precision x,q2,pdfsf(-6:6)
C----------------------------------------------------------------------
!$OMP CRITICAL
      call FPDFXQ(iPDFSET,x,q2,PDFSF,ICheck_QCDNUM)
!$OMP END CRITICAL
C----------------------------------------------------------------------
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
      
