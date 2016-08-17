C------------------
      DOUBLE PRECISION FUNCTION XQG0(k,iq,xb,ix)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      real*8 p_gluon(1:10),p_uval(1:10),p_ubar(1:10)
c      real*8 p_dval(1:10),p_dbar(1:10),p_qstrange(1:10)
!  The PDF input at the evolution starting scale

c      p_gluon=(/6.0610460227606477,0.14928,
c     .     8.6694,0.,0.,0.,0.,0.,0.,0./)
c      p_uval=(/3.6376751243653960,0.64779,
c     .     4.8742999999999999,0.,10.228999999999999,0.,0.,0.,0.,0./)
c      p_ubar=(/0.1120215,-0.16753999999999999,
c     .     1.8717999999999999,0.,0.,0.,0.,0.,0.,0./)
c      p_dval=(/2.0915735564585449,0.64779,4.3596,
c     .     0.,0.,0.,0.,0.,0.,0./)
c      p_dbar=(/0.16235,-0.16754,4.111,
c     .     0.,0.,0.,0.,0.,0.,0./)

      xqg0=0.
c      fs=0.31

c      print*,'gluon'
      IF (IQ.EQ.1) XQG0=gluon(xb) ! 2.37*xb**(-0.3)*(1-xb)**12.                   ! gluon
c      IF (IQ.EQ.1) XQG0=para_st(xb,p_gluon) ! 2.37*xb**(-0.3)*(1-xb)**12.                   ! gluon
c      print*, 'dval, dbar, qstrange'
      IF (IQ.EQ.2) XQG0=dval(xb)+dbar(xb)-qstrange(xb)  !(0.145*xb**(-0.27)+1.6*xb**0.6)*(1-xb)**4.5   ! d
c      IF (IQ.EQ.2) XQG0=para_st(xb,p_dval)+(1-fs)*para_st(xb,p_dbar)  !(0.145*xb**(-0.27)+1.6*xb**0.6)*(1-xb)**4.5   ! d
c      print*,'uval, ubar'
      IF (IQ.EQ.4) XQG0=uval(xb)+ubar(xb)  !(0.16*xb**(-0.26)+3.5*xb**0.7)*(1-xb)**3.7    ! u
c      IF (IQ.EQ.4) XQG0=para_st(xb,p_uval)+para_st(xb,p_ubar)  !(0.16*xb**(-0.26)+3.5*xb**0.7)*(1-xb)**3.7    ! u
c      print*,'qstrange'
      IF (IQ.EQ.6) XQG0=qstrange(xb) !0.108*xb**(-0.29)*(1-xb)**10.                 ! s
c      IF (IQ.EQ.6) XQG0=fs*para_st(xb,p_dbar) !0.108*xb**(-0.29)*(1-xb)**10.                 ! s
c      print*,'dbar,qstrange'
      IF (IQ.EQ.3) XQG0=dbar(xb)-qstrange(xb)  !0.14*xb**(-0.275)*(1-xb)**7                   ! dbar
c      IF (IQ.EQ.3) XQG0=(1-fs)*para_st(xb,p_dbar)  !0.14*xb**(-0.275)*(1-xb)**7                   ! dbar
c      print*,'ubar'
      IF (IQ.EQ.5) XQG0=ubar(xb)  !0.14*xb**(-0.275)*(1-xb)**9                   ! ubar
c      IF (IQ.EQ.5) XQG0=para_st(xb,p_ubar)  !0.14*xb**(-0.275)*(1-xb)**9                   ! ubar
c      print*,'qstrange'
      IF (IQ.EQ.7) XQG0=qstrange(xb) !0.108*xb**(-0.29)*(1-xb)**10.                 ! sbar
c      IF (IQ.EQ.7) XQG0=fs*para_st(xb,p_dbar) !0.108*xb**(-0.29)*(1-xb)**10.                 ! sbar


c      if(ix.eq.1)    print*,'iq = ',iq,', xqg0 = ',XQG0

      RETURN
      END

      double precision function para_st(x,a)
C----------------------------------------------------
C                                               
C standard-like parameterisation:               
C  AF = (a*x**b)*(1 - x)**c*(1 + d*x + e*x**2+f*x**3)-
C     - (ap*x**bp)*(1-x)**cp                          
C                                                     
C-----------------------------------------------------

      implicit none
      double precision x,a(1:10)
      double precision AF
      integer i

c      print*,'in para_st'
c      do i=1,10
c         print*,'a(',i,') = ',a(i)
c      enddo

      AF = a(1)*x**a(2)*(1 - x)**a(3)*(1 + a(4)*x
     $     + a(5)*x**2+a(6)*x**3+a(10)*x**0.5)-a(7)*x**a(8)*(1-x)**a(9)

      para_st = AF

      end

c-----------------------
      real*8 function useralphas(q2,kschemepdf,kordpdf,kpdfset)
      implicit real*8 (a-h,o-z)

!  Must return the value of strong coupling at the scale of q2
!  for the 3-flavour scheme (kschemepdf=0) or 
!  for the 4-flavour scheme (kschemepdf=1) or 
!  for the 5-flavour scheme (kschemepdf=2) 

!  kordpdf=0  --  LO
!  kordpdf=1  --  NLO
!  kordpdf=2  --  NNLO

      character name*80
      real*8 f(-6:6)
      common /foruserpdfs/ inituser
      data inituser /-1/ 

      if (inituser.ne.kpdfset) then 

!  Select for the 5-flavour NNLO PDFs the preliminary NNLO HERAPDFs 
!  (small-alpha_s variant) obtained by courtesy of A.Glasov. In order to make 
!  use of this option one has to copy the file 
!
!  $(top_builddir)/user/HERAPDF1.0_NNLO_1145.LHgrid, 
!
!  which contain the PDF grid in the LHAPDF format, to the standard LHAPDF 
!  grid place. 

        if (kschemepdf.eq.2) then 
         if (kordpdf.eq.2) name='HERAPDF1.0_NNLO_1145.LHgrid'
         if (kordpdf.eq.1) name='abkm09_5_nlo_alp.LHgrid'
        end if

!  Select for the 3-flavour NNLO PDFs ones obtained in the variant of ABKM09 fit
!  with the MSbar definition of the heavy-quark masses (cf. [arXiv:1011.6259]).
!  In order to make use of this option one has to copy the file 
!
!  $(top_builddir)/user/abkm09_3_nnlo_msbar.LHgrid, 
!
!  which contain the PDF grid in the LHAPDF format, to the standard LHAPDF 
!  grid place. 

        if (kschemepdf.eq.0) then 
          if (kordpdf.eq.2) name='abkm09_3_nnlo_msbar.LHgrid'
        end if
        call initPDFSetByName(name)
        call InitPDF(kpdfset)
        inituser=kpdfset
      end if

      useralphas=alphasPDF(sqrt(q2))

      return 
      end      
!------------------------
      real*8 function userpdfs(xb,q2,iq,kschemepdf,kordpdf,kpdfset)
      implicit real*8 (a-h,o-z)

!  Must return the value of PDFs times xb at the scale of q2
!  for the 3-flavour scheme (kschemepdf=0) or
!  for the 4-flavour scheme (kschemepdf=1) or 
!  for the 5-flavour scheme (kschemepdf=2) 

!  kordpdf=0  --  LO
!  kordpdf=1  --  NLO
!  kordpdf=2  --  NNLO
 
!  iq=0 -- gluon
!  iq=1 -- d-quark
!  iq=-1 -- anti-d-quark
!  iq=2 -- u-quark
!  iq=-2 -- anti-u-quark
!  iq=3 -- s-quark
!  iq=-3 -- anti-s-quark
!  iq=4 -- c-quark
!  iq=-4 -- anti-c-quark
!  iq=5 -- b-quark
!  iq=-5 -- anti-b-quark

      character name*80
      real*8 f(-6:6)
      common /foruserpdfs/ inituser
      data inituser /-1/ 

      if (inituser.ne.kpdfset) then 

!  Select for the 5-flavour NNLO PDFs the preliminary NNLO HERAPDFs 
!  (small-alpha_s variant) obtained by courtesy of A.Glasov. In order to make 
!  use of this option one has to copy the file 
!
!  $(top_builddir)/user/HERAPDF1.0_NNLO_1145.LHgrid, 
!
!  which contain the PDF grid in the LHAPDF format, to the standard LHAPDF 
!  grid place. 

        if (kschemepdf.eq.2) then 
          if (kordpdf.eq.2) name='HERAPDF1.0_NNLO_1145.LHgrid'
          if (kordpdf.eq.1) name='abkm09_5_nlo_alp.LHgrid'
        end if

!  Select for the 3-flavour NNLO PDFs ones obtained in the variant of ABKM09 fit
!  with the MSbar definition of the heavy-quark masses (cf. [arXiv:1011.6259]).
!  In order to make use of this option one has to copy the file 
!
!  $(top_builddir)/user/abkm09_3_nnlo_msbar.LHgrid, 
!
!  which contain the PDF grid in the LHAPDF format, to the standard LHAPDF 
!  grid place. 

        if (kschemepdf.eq.0) then 
          if (kordpdf.eq.2) name='abkm09_3_nnlo_msbar.LHgrid'
        end if
        call initPDFSetByName(name)
        call InitPDF(kpdfset)
        inituser=kpdfset
      end if

      call evolvePDF(xb,sqrt(q2),f)
      userpdfs=f(iq)

      return 
      end
