      FUNCTION Pdf_CT (Ihadron, Iparton, X, Q, Ir)
Cpn11 Returns either the regular PDF for iparton=-6:6, or their
c     singlet and non-singlet combinations arising in NNLO evolution
c
c      iparton         PDF type
c        20             Singlet^+ and
c       -20             Singlet^- combinations:
c        Sigma_{+-} = Sum_{i=1}^{Nf} (q_i +- qbar_i)
c
c        21     : 20+nf     non-singlet q_i^+ and 
c       -(20+nf):-21       non-singlet q_i^- combinations:
c        q_i^{+-} = q_i +- qbar_i - Sigma^{+-}/Nf, i = 1:Nf
c
c        30     : Sum_{i=1}^{Nf} e_i^2 (q_i + qbar_i)
c

C========================================================================
C GroupName: Pdfs
C Description: user callable parton distribution functions
C ListOfFiles: pdf pdfh pdfp pdftst
C========================================================================
C $Header: /Net/u52/wkt/h22/2prz/RCS/Pdfs.f,v 7.1 99/08/20 22:07:38 wkt Exp $
C $Log:	Pdfs.f,v $
c Revision 7.2: 2010/12/02 pn + marco
c Added the option to return singlet and non-singlet PDFs
c
c Note: jset is removed as an argument. It is very important to make
c       for variable switching point sure that Setpdf(Iset) 
c        is called before a specific set of pdf is to be used.
c
c Revision 7.1  99/08/20  22:07:38  wkt
c common PdfSwh fixed
c 
c Revision 7.0  99/08/20  10:56:38  wkt
c Renamed Version 7.0 since common block /PdfSwh/ has changed;
c and in order to synchronize with EvlPac7. Iset = 902 (used .ini) added.
c 
c Revision 6.5  1999/08/19  21:23:24  wkt
c Add Iset = 902 option to evolve from .ini file in non-interactive mode
c and with user setable grids.
c
c Revision 6.4  1999/08/13  11:19:14  wkt
c pdfp.f modified to synchronize with new setpdf concerning the added
c option of Iset=Isetin0=900 -- cf. comments in Setpdf.f
c
c Revision 6.3  97/12/17  23:03:41  wkt
c Iset 1410 extended to 1414
c 
c Revision 6.2  97/11/16  00:18:14  wkt
c Revised pdf, pdfh, pdfp : added IpdMod,Iptn0 switches in / PdfSwh /
c See comments in pdf.f module.
c 
c Revision 6.1  97/11/15  18:07:17  wkt
c OverAll pdf functions + parametrized pdf's (with evolution package excluded)
c 
C========================================================================

C These callable parton distributions are, in order of generality:

C      FUNCTION Pdf (Ihadron, Iparton, X, Q, Ir)
C         It is now only a shell to make it compatible to existing 
C         programs calling the PDF function in this format.

C      FUNCTION PdfH (Ihadron, Iparton, X, Q)
C         This is mainly for Hadron targets other than the proton.

C      FUNCTION PdfP (LPRTN, XD, QD)
C         Parton distribution in proton.

C ========================================================================

      IMPLICIT none
      
      double precision pdf_CT, pdfh, x, q, sigma, ei2(6)
      integer ihadron, iparton, 
     >  jparton,aparton, ir, nfl, nf, ipsign, i, iw, nin, nout, nwrt

      external nfl
      COMMON / IOUNIT / NIN, NOUT, NWRT      
      data iw /0/
      data ei2 /4d0,1d0,1d0,4d0,1d0,4d0/ !9*(frac. quark charge)^2

      Nf = nfl(Q)               ! the current number of active flavors
      jparton = iparton

cmarco 
      if (abs(jparton).le.nf) then !return a "physical-flavor" PDF

        Pdf_CT = PdfH (Ihadron, jparton, X, Q) !for iparton=-6:6
      else if (jparton.eq.30) then ! Sum_{i=1}^{Nf} e_i^2 (q_i + qbar_i)
        
        sigma = 0d0         
        do i=1,nf
          sigma = sigma + 
     >      (PdfH (Ihadron, I, X, Q) + PdfH (Ihadron, -I, X, Q) ) !*ei2 we multiply by the charge outside
        enddo                   !i
        PDF_CT = Sigma/9d0

      else if (abs(jparton).ge.20.and.abs(jparton).le.20+nf) then 
                                           !return a (non-)singlet PDF
        ipsign = jparton/abs(jparton)
        
        sigma = 0d0         
        do i=1,nf
          sigma = sigma + PdfH (Ihadron, I, X, Q) 
     >      + ipsign*PdfH (Ihadron, -I, X, Q)
        enddo                   !i

        if (abs(jparton).eq.20) then
          PDF_CT = Sigma
          return
        endif
         
        jparton = abs(jparton)-20 !map jparton on 1:Nfl
        
        PDF_CT = PdfH (Ihadron, jparton, X, Q) 
     >    + ipsign*PdfH (Ihadron, -I, X, Q) - Sigma/Nf
      else
        Call WARNI(IW, NWRT,
     >    'Requested iparton exceeds Nfl: setting PDF=0d0',
     >    'Iparton', Iparton, 0,0,0)
        Pdf_CT=0.D0
      endif                     !iparton 
      Return
C                  *****************************
      End

