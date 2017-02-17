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

      integer i,Q2vIPDFSET
      double precision A,Z, tmpU,tmpD,tmpUb,tmpDb
      data A,Z /207,82/


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
C         call FPDFXQ(Q2viPDFSET(Q2),x,q2,PDFSF,ICheck_QCDNUM)         
         call ALLFXQ(Q2viPDFSET(Q2),x,q2,PDFSF,N_NEUTRAL_PDF,ICheck_QCDNUM)         
C         PDFSF(N_CHARGE_PDF+N_NEUTRAL_PDF) = FSNSXQ(Q2viPDFSET(Q2),13,x,q2,ICheck_QCDNUM)
c         print *,vipdfset,x,q2, PDFSF(N_CHARGE_PDF+1), PDFSF(0)
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

      if (pdfdef(7).eq.1.0) then  ! gluon
         call fflist(vIPDFSET,pdfdef,0,XArray,Q2ARRAY, 
     $        PDFout,Npoints,Icheck)
      else
         call fflist(vIPDFSET,pdfdef,1,XArray,Q2ARRAY, 
     $        PDFout,Npoints,Icheck)
      endif

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
      integer Q2vIPDFSET
      double precision x,q2,pdfsf(-6:6)
C----------------------------------------------------------------------
!$OMP CRITICAL
      call ALLFXQ(Q2viPDFSET(q2),x,q2,PDFSF,0,ICheck_QCDNUM)
!$OMP END CRITICAL
C----------------------------------------------------------------------
      end

*
************************************************************************
*
*     Function needed to pick the PDF table with the correct number of
*     active flavours.
*
************************************************************************
      function Q2vIPDFSET(Q2)
*
      implicit none
*
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
**
*     Input variables
*
      double precision Q2
**
*     Internal variables
*
      double precision Q
      double precision HF_Get_alphas
**
*     Output variables
*
      integer Q2viPDFSET
*
      Q2viPDFSET = vIPDFSET
      if(.not.UseHVFNS) return
*
      Q = dsqrt(Q2)
      if(Q.gt.DataSetSwitchScales(6,cIDataSet))then
         Q2viPDFSET = max(IPDFSET,vIPDFSET)
         call SetMaxFlavourPDFs(min(6,DataSetMaxNF(cIDataSet)))
         call SetMaxFlavourAlpha(min(6,DataSetMaxNF(cIDataSet)))
      elseif(Q.gt.DataSetSwitchScales(5,cIDataSet))then
         Q2viPDFSET = max(IPDFSET+1,vIPDFSET)
         call SetMaxFlavourPDFs(min(5,DataSetMaxNF(cIDataSet)))
         call SetMaxFlavourAlpha(min(5,DataSetMaxNF(cIDataSet)))
      elseif(Q.gt.DataSetSwitchScales(4,cIDataSet))then
         Q2viPDFSET = max(IPDFSET+2,vIPDFSET)
         call SetMaxFlavourPDFs(min(4,DataSetMaxNF(cIDataSet)))
         call SetMaxFlavourAlpha(min(4,DataSetMaxNF(cIDataSet)))
      else
         Q2viPDFSET = max(IPDFSET+3,vIPDFSET)
         call SetMaxFlavourPDFs(min(3,DataSetMaxNF(cIDataSet)))
         call SetMaxFlavourAlpha(min(3,DataSetMaxNF(cIDataSet)))
      endif
*
      return
      end
