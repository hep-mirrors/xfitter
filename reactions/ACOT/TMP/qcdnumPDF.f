C============================================================================= 
      FUNCTION PDF(iSet, IHADRONin, IPARTONin, xIN, q, iRet)
C============================================================================= 
C      Front end program for  cteq4 pdf's: Fred Olness  6/20/96
C
C============================================================================= 
      Implicit Double Precision (A-H, O-Z)
C      Common  / pdftx / IHADRON2  !*** OVERRIDE IHADRON SWITCH (CAREFUL!!!)
C     2/3/2005  Fixed IHad; no longer need to override

      dimension xpdf(-6:6)                             !pdfout
      data small /1.0d-16/

      DATA IPARM,IPRINT /1,1/
      SAVE IPARM,IPRINT,IHADRON
      LOGICAL IFIRST 
      data IFIRST,ICTEQ /.TRUE.,4/
      SAVE IFIRST,ICTEQ
C============================================================================= 
C     NOTE: FOR NOW IGNORE: ISET, IHADRONIN, IRET          fio 18 Jan 2011

      pdf=0.0d0

      x=xIN
c     if(x.gt.0.99d0) x=0.99d0
      if(x.gt.0.99d0) return !*** FIO 19 NOV 2012: pdf->0 for x->1
      q2=q**2   !*** QCDNUM USES Q^2 **********

      call hf_get_pdfs(x,q2,xpdf)                  !interpolate all pdf's

cv      print*,'pdfs',x,q2, xpdf(-1), xpdf(-2), xpdf(0), xpdf(1)-xpdf(-1),xpdf(2)-xpdf(-2)
cv      stop
      tmp= xpdf(IPARTONin)
      if(IPARTONin.eq.+1) tmp = xpdf(+2)
      if(IPARTONin.eq.+2) tmp = xpdf(+1)
      if(IPARTONin.eq.-1) tmp = xpdf(-2)
      if(IPARTONin.eq.-2) tmp = xpdf(-1)
!
      pdf=tmp/x
c      if(pdf.le.0.0d0) pdf=0.0d0   !*** Protect negative PDFs: (true only to NLO order)
      if(pdf.le.small) pdf=0.0d0   !*** Round small PDFs to zero to avoid noise; especially for heavy quarks below threshold

      RETURN
      END
C============================================================================= 
C============================================================================= 
