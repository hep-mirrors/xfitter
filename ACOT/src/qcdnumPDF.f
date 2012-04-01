C============================================================================= 
      FUNCTION PDF(iSet, IHADRONin, IPARTONin, x, q, iRet)
C============================================================================= 
C      Front end program for  cteq4 pdf's: Fred Olness  6/20/96
C
C============================================================================= 
      Implicit Double Precision (A-H, O-Z)
C      Common  / pdftx / IHADRON2  !*** OVERRIDE IHADRON SWITCH (CAREFUL!!!)
C     2/3/2005  Fixed IHad; no longer need to override

      dimension xpdf(-6:6)                             !pdfout


      DATA IPARM,IPRINT /1,1/
      SAVE IPARM,IPRINT,IHADRON
      LOGICAL IFIRST 
      data IFIRST,ICTEQ /.TRUE.,4/
      SAVE IFIRST,ICTEQ
C============================================================================= 
C     NOTE: FOR NOW IGNORE: ISET, IHADRONIN, IRET          fio 18 Jan 2011
               
      if(x.gt.0.99d0) x=0.99d0
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

      RETURN
      END
C============================================================================= 
C============================================================================= 
