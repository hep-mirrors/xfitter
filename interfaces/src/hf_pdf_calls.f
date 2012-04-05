      subroutine HF_Get_PDFs(x,q2,PDFSF)
C----------------------------------------------------------------------
C Interface to PDF 
C
C  Input:  x, Q2 values
C  Output: 13 PDF values
C----------------------------------------------------------------------
      implicit none
C----------------------------------------------------------------------
      include 'steering.inc'
      include 'fcn.inc'
      double precision x,q2,pdfsf(-6:6)
C----------------------------------------------------------------------
      if ( CachePDFs ) then
C Cache PDF calls:
         call GetCachedPDFs(ifcncount,x,q2,PDFSF)
      else
C Get PDFs directly:
         Call HF_Get_PDFs_UnCached(x,q2,PDFSF)
      endif

C----------------------------------------------------------------------
      end



      double precision Function HF_Get_alphas(q2)
C----------------------------------------------------
C
C  Return alpha_S value for a given Q2
C
C----------------------------------------------------
      implicit none 
      double precision Q2
      integer nf,ierr
      double precision ASFUNC

C----------------------------------------------------
      HF_Get_alphas = ASFUNC(q2,nf,ierr) 
C----------------------------------------------------
      end


      double precision function hf_get_mur(iDataSet)
C------------------------------------------------------
C
C Get data-set depdendent renormalisation scale
C
C------------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'scales.inc'
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
      include 'ntot.inc'
      include 'scales.inc'
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
      include 'steering.inc'
      double precision x,q2,pdfsf(-6:6)
C----------------------------------------------------------------------
      call FPDFXQ(iPDFSET,x,q2,PDFSF,ICheck_QCDNUM)
C----------------------------------------------------------------------
      end
