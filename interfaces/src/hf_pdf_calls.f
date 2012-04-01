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
      double precision x,q2,pdfsf(-6:6)
C----------------------------------------------------------------------
      call FPDFXQ(iPDFSET,x,q2,PDFSF,ICheck_QCDNUM)
C----------------------------------------------------------------------
      end
