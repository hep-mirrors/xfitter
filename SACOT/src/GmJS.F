C      Module GM-JS :  Interface to FFNS programs of Jack Smith etal.

C      Light quark initiated pieces 10/25/2010 by MG 
C      F77 built by HLL.  07.07 - 07.08


      Real*8 Function F2hg2(X, Q, aMu, aMs, Ihad)
      Use mDisPac
      implicit none
	real*8 :: tmp1, smpsnf, a5

      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad

      Real*8 :: ADZINT_CT, alpi, F2HQgint !functions
      external F2HQgint
      Real*8, Parameter :: D1=1D0 
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST, F2g
      Integer :: IER
      Real*8 :: XX, QQ, bmu, bms, eps, aMu0, als4pi, zmax, Xi
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      bms=aMs
      QQ=Q
      XX=X
      Xi=X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      If(Xi.ge.1d0) then
         F2hg2=0d0
         return
      Endif
      Ihadn=Ihad

      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

	a5 = zmax - 0.00000000001d0*(zmax-xx)

      F2g = 
     >   SMPSNF(F2HQgint,xx + 0.00000000001d0*(a5-xx), 
     >	                 xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(F2HQgint,xx + 0.000001d0*(a5-xx),
     >                   xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQgint,xx + 0.00001d0*(a5-xx),
     >                   xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQgint,xx + 0.0001d0*(a5-xx),
     >                   xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQgint,xx + 0.001d0*(a5-xx),
     >                   xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQgint,xx + 0.01d0*(a5-xx),
     >                   xx + 0.1d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(F2HQgint,xx + 0.1d0*(a5-xx),
     >                   xx + 0.9d0*(a5-xx),16, ErrEst, Ier) 
     > + SMPSNF(F2HQgint,xx + 0.9d0*(a5-xx),
     >                   xx + 0.99d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(F2HQgint,xx + 0.99d0*(a5-xx),
     >                   xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQgint,xx + 0.999d0*(a5-xx),
     >                   xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQgint,xx + 0.9999d0*(a5-xx),
     >                   xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQgint,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c      F2g = ADZINT_CT
c     >  (F2HQgint, XX, zmax, AERR,RERR,
c     >	ERREST,IER, IACTA, IACTB)



c -----------------------------------------------------------------

c -----------------------------------------------------------------

      F2hg2= (als4pi)**2 * F2g  

      end


      Real*8 Function F2HQgint(z)
      implicit none
      Real*8, intent(in):: z

      Real*8 :: HQcoef, PDFH   !functions
      Real*8 :: Coef, gluon
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn

      Coef= HQcoef(2,z,QQ,bmu,bms)
      gluon=PDFH(Ihadn, 0, XX/z, bmu)
      F2HQgint=gluon*Coef/z !*XX the multiplication by XX is made outside in 4fit/chidis.f
      End


cmarco F2hlPS2: contributions prop. to the charge of heavy quarks////////////////////
      Real*8 Function F2hlPS2(X, Q, aMu, aMs, Ihad)
      implicit none
	real*8 :: tmp1, smpsnf, a5
	real*8 :: tmp
      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad

      Real*8 :: ADZINT_CT, alpi !functions
      Real*8 :: F2hint
      external F2hint 

      Real*8, Parameter :: AERR=1D-16, RERR=1D-4,D1=1D0
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST, tF2h
      Integer :: IER
      Real*8 :: Xi, aMu0, als4pi
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      bms=aMs
      QQ=Q
      XX=X
      Xi=X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      If(Xi.ge.1d0) then
         F2hlPS2=0d0
         return
      Endif
      Ihadn=Ihad

      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

	a5 = zmax - 0.00000000001d0*(zmax-xx)

        tF2h = 
     >   SMPSNF(F2hint,xx + 0.00000000001d0*(a5-xx), 
     >	               xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(F2hint,xx + 0.000001d0*(a5-xx),
     >                 xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2hint,xx + 0.00001d0*(a5-xx),
     >                 xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2hint,xx + 0.0001d0*(a5-xx),
     >                 xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2hint,xx + 0.001d0*(a5-xx),
     >                 xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2hint,xx + 0.01d0*(a5-xx),
     >                 xx + 0.1d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(F2hint,xx + 0.1d0*(a5-xx),
     >                 xx + 0.9d0*(a5-xx),16, ErrEst, Ier) 
     > + SMPSNF(F2hint,xx + 0.9d0*(a5-xx),
     >                 xx + 0.99d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(F2hint,xx + 0.99d0*(a5-xx),
     >                 xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2hint,xx + 0.999d0*(a5-xx),
     >                 xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2hint,xx + 0.9999d0*(a5-xx),
     >                 xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2hint,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx        tF2h = ADZINT_CT
c xx     >  (F2hint, XX, zmax, AERR,RERR,
c xx     >	ERREST,IER, IACTA, IACTB)

c ---------------------------------------------------------------

c ---------------------------------------------------------------

        F2hlPS2= (als4pi)**2 * tF2h 

      end                       !F2hlPS2
cmarco F2hlPS2: contributions prop. to the charge of heavy quarks//////////////



cmarco F2hint quark initiated -> photon-heavy qrk 
      Real*8 Function F2hint(z)
      implicit none
      Real*8, intent(in):: z
      Real*8 :: HQhqlcoef,PDF_CT !functions
      Real*8 :: Coef, sum
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn, Ir
      Common /F2LHQhad/ Ihadn
      external PDF_CT

      Coef= HQhqlcoef(2,z,QQ,bmu,bms)
      sum=PDF_CT(ihadn,20, XX/Z, bmu, Ir) !(u+ub+d+db+s+sb+...)
                                       !Nf active flavor is 
                                       !read from subroutine ./flib/2prz/pdf.f
      F2hint=sum*Coef/z !*XX multiplication by XX is made outside in 4fit/chidis.f
      End
cmarco F2hint quark initiated -> photon-heavy qrk




cmarco===================================================================
c Semi-inclusive F2l: contributions prop. to the charge of light quarks
c This corresponds to the d(eta,xi) terms in Eq 5 PLB 347 1995 Riemersma
cmarco===================================================================
      Real*8 Function SiF2l(X, Q, aMu, aMs, Ihad, iPrtn)
      implicit none
	real*8 :: tmp1, smpsnf, a5

      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad, iPrtn

      Real*8 :: ADZINT_CT, alpi !functions
      Real*8 :: F2lint
      external F2lint
      Real*8, Parameter ::  AERR=1D-16, RERR=1D-4,D1=1D0
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST, tF2l
      Integer :: IER
      Real*8 :: aMu0, als4pi, Xi
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn
      Data aMu0 /-1d0/
      save aMu0

      jPrtn=iPrtn

      bms=aMs
      QQ=Q
      XX=X
      Xi=X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      If(Xi.ge.1d0) then
         SiF2l=0d0
         return
      Endif
      Ihadn=Ihad

      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

	a5 = zmax - 0.00000000001d0*(zmax-xx)

      tF2l = 
     >   SMPSNF(F2lint,xx + 0.00000000001d0*(a5-xx), 
     >	               xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(F2lint,xx + 0.000001d0*(a5-xx),
     >                 xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lint,xx + 0.00001d0*(a5-xx),
     >                 xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lint,xx + 0.0001d0*(a5-xx),
     >                 xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lint,xx + 0.001d0*(a5-xx),
     >                 xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lint,xx + 0.01d0*(a5-xx),
     >                 xx + 0.1d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(F2lint,xx + 0.1d0*(a5-xx),
     >                 xx + 0.9d0*(a5-xx),16, ErrEst, Ier) 
     > + SMPSNF(F2lint,xx + 0.9d0*(a5-xx),
     >                 xx + 0.99d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(F2lint,xx + 0.99d0*(a5-xx),
     >                 xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lint,xx + 0.999d0*(a5-xx),
     >                 xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lint,xx + 0.9999d0*(a5-xx),
     >                 xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lint,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx      tF2l= ADZINT_CT
c xx     >  (F2lint, XX, zmax, AERR,RERR,
c xx     >	ERREST,IER, IACTA, IACTB)


    
        SiF2l= (als4pi)**2 * tF2l

      end                       !SiF2l
cmarco Semi-inclusive F2l: contributions prop. to the charge of light quarks



cmarco Semi-inclusive F2lint light quark initiated -> photon-light qrk
      Real*8 Function F2lint(z)
      implicit none
      Real*8, intent(in):: z

      Real*8 :: HQlqcoef, PDF_CT !functions
      Real*8 :: Coef, suml
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Integer :: jNf, jPrtn, Ir
      Common /F2L_NF / jNf, jPrtn
      external PDF_CT

      Coef= HQlqcoef(2,z,QQ,bmu,bms)

      suml = PDF_CT(ihadn, jPrtn, XX/Z, bmu, Ir) 
     >     + PDF_CT(ihadn,-jPrtn, XX/Z, bmu, Ir)

      F2lint=suml*Coef/z !*XX the multiplication by XX is made outside in 4fit/chidis.f

      End
cmarco Semi-inclusive F2lint light quark initiated -> photon-light qrk





cmarco FLhlPS2 contributions prop. to the charge of heavy quarks////////////////////
      Real*8 Function FLhlPS2(X,Q,aMu,aMs,Ihad)
!     NNLO singlet contribution to F_long=F_2/(2x)-F_1 (0 helicity S.F.),
!     as defined by CTEQ

      Use mDisPac
      implicit none
	real*8 :: tmp1, smpsnf, a5

      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad

      Real*8 :: ADZINT_CT, alpi !functions
      Real*8 :: FLhint
      external FLhint 

      Real*8, Parameter :: D1=1D0
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST, tFLh
      Integer :: IER
      Real*8 :: XX, QQ, bmu, bms, eps, aMu0, als4pi, zmax, Xi
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      bms=aMs
      QQ=Q
      XX=X
      Xi=X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      If(Xi.ge.1d0) then
         FLhlPS2=0d0
         return
      Endif
      Ihadn=Ihad

      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

	a5 = zmax - 0.00000000001d0*(zmax-xx)

        tFLh = 
     >   SMPSNF(FLhint,xx + 0.00000000001d0*(a5-xx), 
     >	               xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(FLhint,xx + 0.000001d0*(a5-xx),
     >                 xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLhint,xx + 0.00001d0*(a5-xx),
     >                 xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLhint,xx + 0.0001d0*(a5-xx),
     >                 xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLhint,xx + 0.001d0*(a5-xx),
     >                 xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLhint,xx + 0.01d0*(a5-xx),
     >                 xx + 0.1d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(FLhint,xx + 0.1d0*(a5-xx),
     >                 xx + 0.9d0*(a5-xx),16, ErrEst, Ier) 
     > + SMPSNF(FLhint,xx + 0.9d0*(a5-xx),
     >                 xx + 0.99d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLhint,xx + 0.99d0*(a5-xx),
     >                 xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLhint,xx + 0.999d0*(a5-xx),
     >                 xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLhint,xx + 0.9999d0*(a5-xx),
     >                 xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLhint,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx        tFLh= ADZINT_CT
c xx     >    (FLhint, XX, zmax, AERR,RERR,
c xx     >	  ERREST,IER, IACTA, IACTB)


        FLhlPS2= (als4pi)**2 * tFLh 

cmarco Divide by 2, to match the normalization of CTEQ's F_{long}
        FLhlPS2 =  FLhlPS2/2d0

      return
      End                       !FLh
cmarco FLhlPS2 contributions prop. to the charge of heavy quarks////////////////////





cmarco FLhint quark initiated -> photon-heavy qrk
      Real*8 Function FLhint(z)
      implicit none
      Real*8, intent(in):: z

      Real*8 :: HQhqlcoef, PDF_CT !functions
      Real*8 :: Coef, sum
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn, Ir
      Common /F2LHQhad/ Ihadn
      external PDF_CT

      Coef= HQhqlcoef(0,z,QQ,bmu,bms)
      sum=PDF_CT(ihadn,20, XX/Z, bmu, Ir) !(u+ub+d+db+s+sb+...)
                                       !Nf active flavor is 
                                       !read from subroutine ./flib/2prz/pdf.f
      FLhint=sum*Coef/z !*XX multiplication by XX is made outside in 4fit/chidis.f

      End
cmarco FLhint quark initiated -> photon-heavy qrk






cmarco================================================================
c Semi-inclusive FLl contributions prop. to the charge of light quarks
c NNLO singlet contribution to F_long=F_2/(2x)-F_1 (0 helicity S.F.),
c as defined by CTEQ
cmarco================================================================
      Real*8 Function SiFLl(X, Q, aMu, aMs, Ihad, iPrtn)
      Use mDisPac
      implicit none
	real*8 :: tmp1, smpsnf, a5
      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad, iPrtn

      Real*8 :: ADZINT_CT, alpi !functions
      Real*8 :: FLlint
      external FLlint

      Real*8, Parameter :: D1=1D0
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST, tFLl
      Integer :: IER
      Real*8 :: XX, QQ, bmu, bms, eps, aMu0, als4pi, zmax, Xi
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn
      Data aMu0 /-1d0/
      save aMu0

      jPrtn=iPrtn

      bms=aMs
      QQ=Q
      XX=X
      Xi=X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      If(Xi.ge.1d0) then
         SiFLl=0d0
         return
      Endif
      Ihadn=Ihad

      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

	a5 = zmax - 0.00000000001d0*(zmax-xx)

        tFLl = 
     >   SMPSNF(FLlint,xx + 0.00000000001d0*(a5-xx), 
     >	               xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(FLlint,xx + 0.000001d0*(a5-xx),
     >                 xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLlint,xx + 0.00001d0*(a5-xx),
     >                 xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLlint,xx + 0.0001d0*(a5-xx),
     >                 xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLlint,xx + 0.001d0*(a5-xx),
     >                 xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLlint,xx + 0.01d0*(a5-xx),
     >                 xx + 0.1d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLlint,xx + 0.1d0*(a5-xx),
     >                 xx + 0.9d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(FLlint,xx + 0.9d0*(a5-xx),
     >                 xx + 0.99d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLlint,xx + 0.99d0*(a5-xx),
     >                 xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLlint,xx + 0.999d0*(a5-xx),
     >                 xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLlint,xx + 0.9999d0*(a5-xx),
     >                 xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLlint,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx        tFLl= ADZINT_CT
c xx     >    (FLlint, XX, zmax, AERR,RERR,
c xx     >	    ERREST,IER, IACTA, IACTB)

    
        SiFLl= (als4pi)**2 * tFLl

cmarco Divide by 2, to match the normalization of CTEQ's F_{long}
      SiFLl =  SiFLl/2d0

      return
      End                       !FLl
cmarco FLl contributions prop. to the charge of light quarks////////////////////



cmarco FLlint light quark initiated -> photon-light qrk
      Real*8 Function FLlint(z)
      implicit none
      Real*8, intent(in):: z

      Real*8 :: HQlqcoef, PDF_CT !functions
      Real*8 :: Coef, suml
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Integer :: jNf, jPrtn, Ir
      Common /F2L_NF / jNf, jPrtn
      external PDF_CT

      Coef= HQlqcoef(0,z,QQ,bmu,bms)

      suml = PDF_CT(ihadn, jPrtn, XX/Z, bmu, Ir) 
     >     + PDF_CT(ihadn,-jPrtn, XX/Z, bmu, Ir)

      FLlint=suml*Coef/z !*XX the multiplication by XX is made outside in 4fit/chidis.f

      End
cmarco FLlint light quark initiated -> photon-light qrk


c=======================================================
c New Subtraction term introduced in Eq 19 short paper
c Taken from Eq.B.5 of EPJC1 301 1998 Buza et al
c This is for F2 
c=======================================================


      Real*8 Function F2lNSSub(X, Q, aMu, aMs, Ihad, iPrtn)

      Use mDisPac
      implicit none
	real*8 :: tmp1, smpsnf, a5

      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad, iPrtn

      Real*8 :: ADZINT_CT, alpi, F2lNSSubint !functions
      Real*8 :: fchi, Delta, Remn, PDF_CT!functions
      Real*8 :: suml0
      external F2lNSSubint
      Real*8, Parameter :: D1=1D0
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST, tF2lNSSub
      Integer :: IER
      Real*8, parameter:: Zeta2=1.6449340668482264365D0
      Real*8, parameter:: Zeta3=1.2020569031595942854D0
      Real*8 :: XX, QQ, bmu, bms, eps, aMu0, als4pi, zmax, Xi
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn, Ir
      Common /F2LHQhad/ Ihadn
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn
      Data aMu0 /-1d0/
      save aMu0

      jPrtn = iPrtn

      bms=aMs
      QQ=Q
      XX=X
      Xi=fchi(X,Q,bms,bms)
      zmax=X/Xi
      If(Xi.ge.1d0) then
         F2lNSSub=0d0
         return
      Endif
      Ihadn=Ihad
      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

      suml0=PDF_CT(ihadn, jprtn, XX, bmu, Ir) + 
     >  PDF_CT(ihadn, -jprtn, XX, bmu, Ir)

                                !Dirac delta contribution
      Delta=Cf*Tf*suml0*(2d0*dlog(eps)**2 + (16d0/3d0*Zeta2 + 2d0/3d0)
     >  *dlog(eps) + (-8d0/3d0*Zeta3 + 40d0/9d0*Zeta2 +73d0/18d0))

      Remn=Cf*Tf*suml0*dlog(1d0 - XX)*(8d0/3d0*dlog(eps)**2 + 80d0/9d0
     >  *dlog(eps) + 224d0/27d0) ! This is the f(x)*log(1-x) after the Plus integral

c C^0 (*) A^(2)NSqq,H (*) (u+ub)+(d+db)+...

	a5 = zmax - 0.00000000001d0*(zmax-xx)

      tF2lNSSub = 
     >   SMPSNF(F2lNSSubint,xx + 0.00000000001d0*(a5-xx), 
     >	                    xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(F2lNSSubint,xx + 0.000001d0*(a5-xx),
     >                      xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lNSSubint,xx + 0.00001d0*(a5-xx),
     >                      xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lNSSubint,xx + 0.0001d0*(a5-xx),
     >                      xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lNSSubint,xx + 0.001d0*(a5-xx),
     >                      xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lNSSubint,xx + 0.01d0*(a5-xx),
     >                      xx + 0.1d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lNSSubint,xx + 0.1d0*(a5-xx),
     >                      xx + 0.9d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(F2lNSSubint,xx + 0.9d0*(a5-xx),
     >                      xx + 0.99d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lNSSubint,xx + 0.99d0*(a5-xx),
     >                      xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lNSSubint,xx + 0.999d0*(a5-xx),
     >                      xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lNSSubint,xx + 0.9999d0*(a5-xx),
     >                      xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2lNSSubint,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx      tF2lNSSub= ADZINT_CT
c xx     >  (F2lNSSubint, XX, zmax, AERR,RERR,
c xx     >	ERREST,IER, IACTA, IACTB)
                                                                   
      F2lNSSub = (als4pi)**2 * (tF2lNSSub + Delta + Remn) 

      End


      Real*8 Function F2lNSSubint(z) 
      implicit none
      Real*8, intent(in):: z

      Real*8 :: A2NSqqPlus, A2NSqqFinite, PDF_CT !functions
      Real*8 :: suml, suml0, Sub1, Sub2
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Real*8 :: xz, A1Hgz
      common /Cf2sub/ xz, A1Hgz
      Integer Ihadn, ir
      Common /F2LHQhad/ Ihadn
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn

      Sub1= A2NSqqPlus(z/zmax,eps)
      Sub2= A2NSqqFinite(z/zmax,eps)

      suml= PDF_CT(ihadn, jprtn, XX/z, bmu, Ir) +
     >  PDF_CT(ihadn,-jprtn, XX/z, bmu, Ir)

      suml0= PDF_CT(ihadn, jprtn, XX/zmax, bmu, Ir) +
     >  PDF_CT(ihadn,-jprtn, XX/zmax, bmu, Ir)

      F2lNSSubint = (suml/z - suml0/zmax)*Sub1+(suml*Sub2/z)

      End



      Real*8 Function F2HQSub2(X, Q, aMu, aMs, Ihad)
      Use mDisPac
      implicit none
	real*8 :: tmp1, smpsnf, a5
	   real *8 :: xx1, xx2, xii, fii
	   integer itst, ii

      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad

      Real*8 :: ADZINT_CT, alpi, F2HQS2, F2HQS3, F2HQS4 !functions
      Real*8 :: AnalyticF2HQS1
      Real*8 :: fchi !functions
      external F2HQS2, F2HQS3, F2HQS4
c$$$      real*8 :: adz2nt, F2HQS1
c$$$      external F2HQS1
      external AnalyticF2HQS1
      Real*8, Parameter :: D1=1D0
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST, F2gSub1, F2gSub2, F2gSub3, F2qSub4
      Integer :: IER

      Real*8 :: XX, QQ, bmu, bms, eps, aMu0, als4pi, zmax, Xi
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      bms=aMs
      QQ=Q
      XX=X

cpn2010 Xi=fchi(x,Q,bms,bms)
      Xi=fchi(X,Q,bms,bms)!X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      If(Xi.ge.1d0) then
         F2HQSub2=0d0
         return
      Endif
      Ihadn=Ihad
      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

cmarco11 Numerical double convolution replaced by a faster semi-
c        analytical version
c      F2gSub1= Adz2nt
c     >  (F2HQS1, XX, zmax, AERR,RERR,ERREST,IER, IACTA, IACTB)! C^1 (*) A^(1)Hg (*) g

c C^1 (*) A^(1)Hg (*) g

	a5 = zmax - 0.00000000001d0*(zmax-xx)

        F2gSub1 = 
     >   SMPSNF(AnalyticF2HQS1,xx + 0.00000000001d0*(a5-xx), 
     >	                       xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(AnalyticF2HQS1,xx + 0.000001d0*(a5-xx),
     >                         xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(AnalyticF2HQS1,xx + 0.00001d0*(a5-xx),
     >                         xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(AnalyticF2HQS1,xx + 0.0001d0*(a5-xx),
     >                         xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(AnalyticF2HQS1,xx + 0.001d0*(a5-xx),
     >                         xx + 0.01d0*(a5-xx),4, ErrEst, Ier)  
     > + SMPSNF(AnalyticF2HQS1,xx + 0.01d0*(a5-xx),
     >                         xx + 0.10d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(AnalyticF2HQS1,xx + 0.1d0*(a5-xx),
     >                         xx + 0.9d0*(a5-xx),8, ErrEst, Ier)
     > + SMPSNF(AnalyticF2HQS1,xx + 0.9d0*(a5-xx),
     >                         xx + 0.99d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(AnalyticF2HQS1,xx + 0.99d0*(a5-xx),
     >                         xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(AnalyticF2HQS1,xx + 0.999d0*(a5-xx),
     >                         xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(AnalyticF2HQS1,xx + 0.9999d0*(a5-xx),
     >                         xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(AnalyticF2HQS1,xx + 0.99999d0*(a5-xx),a5,4, 
     >	ErrEst, Ier)

c xx        F2gSub1 = ADZINT_CT
c xx     >  (AnalyticF2HQS1, XX, zmax, AERR,RERR,
c xx     >	ERREST,IER, IACTA, IACTB)



	a5 = zmax - 0.00000000001d0*(zmax-xx)

      F2gSub2 = 
     >   SMPSNF(F2HQS2,xx + 0.00000000001d0*(a5-xx), 
     >	               xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(F2HQS2,xx + 0.000001d0*(a5-xx),
     >                 xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS2,xx + 0.00001d0*(a5-xx),
     >                 xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS2,xx + 0.0001d0*(a5-xx),
     >                 xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS2,xx + 0.001d0*(a5-xx),
     >                 xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS2,xx + 0.01d0*(a5-xx),
     >                 xx + 0.1d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS2,xx + 0.1d0*(a5-xx),
     >                 xx + 0.9d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(F2HQS2,xx + 0.9d0*(a5-xx),
     >                 xx + 0.99d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS2,xx + 0.99d0*(a5-xx),
     >                 xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS2,xx + 0.999d0*(a5-xx),
     >                 xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS2,xx + 0.9999d0*(a5-xx),
     >                 xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS2,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx      F2gSub2 = ADZINT_CT
c xx     >  (F2HQS2, XX, zmax, AERR,RERR,ERREST,IER, IACTA, IACTB)! C0 (*) A^(2)Hg (*) g


cmarco The condition for matching Nf alpha_s to (Nf-1) virtual
cmarco corrections in J. Smith's code is included in SFnGm.f
c      F2gSub3= ADZINT_CT
c     >  (F2HQS3, XX, zmax, AERR,RERR,ERREST,IER, IACTA, IACTB) 

      F2gSub3=0d0

	a5 = zmax - 0.00000000001d0*(zmax-xx)

      F2qSub4 = 
     >   SMPSNF(F2HQS4,xx + 0.00000000001d0*(a5-xx), 
     >	               xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(F2HQS4,xx + 0.000001d0*(a5-xx),
     >                 xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS4,xx + 0.00001d0*(a5-xx),
     >                 xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS4,xx + 0.0001d0*(a5-xx),
     >                 xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS4,xx + 0.001d0*(a5-xx),
     >                 xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS4,xx + 0.01d0*(a5-xx),
     >                 xx + 0.1d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS4,xx + 0.1d0*(a5-xx),
     >                 xx + 0.9d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(F2HQS4,xx + 0.9d0*(a5-xx),
     >                 xx + 0.99d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS4,xx + 0.99d0*(a5-xx),
     >                 xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS4,xx + 0.999d0*(a5-xx),
     >                 xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS4,xx + 0.9999d0*(a5-xx),
     >                 xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2HQS4,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx      F2qSub4= ADZINT_CT
c xx     >  (F2HQS4, XX, zmax, AERR,RERR,ERREST,IER, IACTA, IACTB) !C0 (*) A^(2)PSHq (*) Sigma

      F2HQSub2= (als4pi)**2 *(F2gSub1 + F2gSub2 + F2gSub3 +F2qSub4) 
                                                                    
      end


      Real*8 Function F2HQS1(z) ! A^1 (*) C^1 (*) g
      implicit none
      Real*8, intent(in):: z

      Real*8 :: A1Hg !function
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Real*8 :: A1Hgxz, C1Hg, G1Q(5)

      A1Hgxz=A1Hg(XX/z,eps)
      Call SfZm1Q  (iHadn, 0, z/zmax, QQ, bmu, G1Q)
      C1Hg=(G1Q(1)+G1Q(3))*4d0 ! ~ F2 taken out x, charge, & als4pi factors

      F2HQS1=A1Hgxz*C1Hg/z
      End


      Real*8 Function AnalyticF2HQS1(z) ! A^1 (*) C^1 (*) g
cmarco11**********************************************************
c Here we implemented analytic convolutions for A^1 (*) C^1 (*)
c to speed up the computation of A^1 (*) C^1 (*) g
cmarco*********************************************************

      Use mDisPac
      implicit none
      external xLi
      Real*8, intent(in):: z
      Real*8 :: ConvC1A1, gluon, xLi, PDFH !function
      Real*8, parameter:: Zeta2=1.6449340668482264365D0
      Real*8 :: z1
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn
      
      gluon=PDFH(Ihadn, 0, XX/z, bmu)
      z1=z/zmax

      ConvC1A1=Cf*Tf*(20d0 - 48d0*z1 + 64d0*z1**2 + 24d0*Zeta2 
     >  - 48d0*z1*Zeta2 + 64d0*z1**2*Zeta2 + 28d0*dlog(1d0 - z1) 
     >  - 64d0*z1*dlog(1d0 - z1) 
     >  + 48d0*z1**2*dlog(1d0 - z1) - 8d0*dlog(1d0 - z1)**2 
     >  + 16d0*z1*dlog(1d0 - z1)**2 
     >  - 16d0*z1**2*dlog(1d0 - z1)**2 
     >  + 4d0*dlog(z1) + 32d0*z1*dlog(z1) 
     >  - 48d0*z1**2*dlog(z1)- 4d0*dlog(z1)**2 + 8d0*z1*dlog(z1)**2 
     >  - 16d0*z1**2*dlog(z1)**2 - 16d0*xLi(2,1d0 - z1) 
     >  + 32d0*z1*xLi(2,1d0 - z1) - 32d0*z1**2*xLi(2,1d0 - z1) 
     >  - 8d0*xLi(2,z1)+ 16d0*z1*xLi(2,z1) - 32d0*z1**2*xLi(2,z1))
     >  *dlog(eps)
   
      if(QQ.ne.bmu) then
      ConvC1A1 = ConvC1A1 + Cf*Tf*(4d0 - 16d0*z1 - 16d0*dlog(1d0 - z1) 
     >    + 32d0*z1*dlog(1d0 - z1) - 32d0*z1**2*dlog(1d0 - z1) 
     >    + 8d0*dlog(z1) - 16d0*z1*dlog(z1) 
     >    + 32d0*z1**2*dlog(z1))*dlog(eps)*dlog(QQ**2/bmu**2)
      endif
      
      AnalyticF2HQS1=ConvC1A1*gluon/z
      End


      Real*8 Function F2HQS2(z) ! A^2 (*) C^0 (*)g
      implicit none
      Real*8, intent(in):: z

      Real*8 :: A2Hg, PDFH !functions
      Real*8 :: gluon, Sub2
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Real*8 :: xz, A1Hgz
      common /Cf2sub/ xz, A1Hgz
      Integer Ihadn
      Common /F2LHQhad/ Ihadn

      Sub2= A2Hg(z/zmax,eps)
      gluon=PDFH(Ihadn, 0, XX/z, bmu)

      F2HQS2=gluon*Sub2/z
      End

      Real*8 Function F2HQS3(z) ! substraction due to the difference of alpha_S from FFNS to VFNS
      implicit none
      Real*8, intent(in):: z

      Real*8 :: A1Hg, PDFH !functions
      Real*8 :: gluon, Sub3
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Real*8 :: xz, A1Hgz
      common /Cf2sub/ xz, A1Hgz
      Integer Ihadn
      Common /F2LHQhad/ Ihadn

      A1Hgz = A1Hg(z/zmax,eps)

      Sub3= dlog(1d0/eps)*A1Hgz* 2d0/3d0  

      gluon = PDFH(Ihadn, 0, XX/z, bmu)

      F2HQS3 = gluon *Sub3 /z

      End

cmarco
      Real*8 Function F2HQS4(z) ! A2PSHq (*) C^0 (*)(u+ub)+(d+db)+...
      implicit none
      Real*8, intent(in):: z

      Real*8 :: A2PSHq, PDF_CT !functions
      Real*8 :: Sub, sum
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Real*8 :: xz, A1Hgz
      common /Cf2sub/ xz, A1Hgz 
      Integer Ihadn, ir
      Common /F2LHQhad/ Ihadn

      sum=PDF_CT(ihadn,20, XX/Z, bmu, Ir) !(u+ub+d+db+s+sb+...)

      Sub= A2PSHq(z/zmax,eps)
      F2HQS4=sum*Sub/z

      End
cmarco






c     Real*8 Function F2sub1(z) !A^1 (*) C^1
c     Real*8, intent(in):: z
c
c      Real*8 :: A1Hg, c2np1a0, c2np1Dk !functions
c      Real*8 :: XX, QQ, bmu, bms, eps, zmax
c     common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
c      Real*8 :: xz, A1Hgz, A1Hgzx, CL, Dk
c      common /CF2sub/ xz, A1Hgz
c
c      A1Hgzx=A1Hg(xz/z,eps)
c      CL= c2np1a0(z,0)
c      Dk= c2np1Dk(A1Hgzx,A1Hgz,z,0)

c      F2sub1=CL*A1Hgzx/z + Dk

c      End


c      Function c2np1a0(y,Nf)
c NLO non-singlet coefficient functions for F_2
c References:
c   S. Moch, J. Vermaseren and A. Vogt, hep-ph/0504242

c      Implicit none
c      Real*8 :: c2np1a0, c2np1Dk, c2s1aL, c2g1aL, Dkfun, c2np1Dk0
c      Real*8 :: y, Dkk(0:1)
c      Integer :: Nf

c      Real*8, Parameter:: Cf=4d0/3D0
c      Real*8 :: DL, y1, DL1, QA, QA0

c      DL = Log(y)
c      y1 = 1d0 - y
c      DL1 = Log(y1)

c      c2np1a0= Cf*(-2.*(1.+y)*(DL1-DL) -4./y1*DL +6. +4.*y)

c      Return

c      Entry c2np1Dk(QA, QA0, y, Nf)

c      c2np1Dk = Cf*(4.*Dkfun(1,y,QA,QA0) -3.*Dkfun(0,y,QA,QA0))
c      Dkk(1)=Dkfun(1,y,QA,QA0)
c      Dkk(0)=Dkfun(0,y,QA,QA0)
c      c2np1Dk = c2np1Dk0(Dkk, Cf)

c      Return

c      Entry c2s1aL(y, Nf)
c NLO pure-singlet coefficient functions for F_2
c      c2s1aL=0d0
c      Return

c      Entry c2g1aL(y, Nf)
c NLO gluon coefficient functions for F_2

c      DL = Log(y)
c      y1 = 1d0 - y
c      DL1 = Log(y1)

c      c2g1aL= Nf*((2. -4.*y*y1)*(DL1-DL) -2. +16.*y*y1)

c      Return

c      End



      Real*8 Function F2HQg0(X, Q, aMu, aMs, Ihad, Nptn)
      Use mDisPac
      implicit none
	real*8 :: tmp, tmp1, smpsnf, a5

      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad, Nptn

      Real*8 :: ADZINT_CT, alpi, F2HQg0int !functions
      Real*8 :: efac(6)=(/4d0,1d0,1D0,4D0,1D0,4D0/)  ! e^2*9 for each quark
      external F2HQg0int
      Real*8, Parameter :: D1=1D0
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST
      Integer :: IER

      Real*8 :: XX, QQ, bmu, bms, eps, aMu0, als4pi, zmax, Xi
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      bms=aMs
      QQ=Q
      XX=X
      Xi=X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      If(Xi.ge.1d0) then
         F2HQg0=0d0
         return
      Endif
      Ihadn=Ihad

      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

      
c ******************************************************************
c force an error to reexamine this -- so far, this code is not used.
c ******************************************************************
	if(QQ*QQ .ge. 0.d0) then
	   print *,'fatal error -- F2HQg0int not yet optimized'
	   stop
	endif

      tmp = ADZINT_CT
     >  (F2HQg0int, XX, zmax, AERR/1000.d0,RERR/1000.d0,
     >	ERREST,IER, IACTA, IACTB)

	a5 = zmax - 0.00000000001d0*(zmax-xx)

	tmp1 = 
     >   SMPSNF(F2HQg0int,xx + 0.00000000001d0*(a5-xx), 
     >	       xx + 0.00000000001d0*(a5-xx),64, ErrEst, Ier)
     > + SMPSNF(F2HQg0int,xx + 0.00000000001d0*(a5-xx),
     >         xx + 0.0000000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.0000000001d0*(a5-xx),
     >         xx + 0.000000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.000000001d0*(a5-xx),
     >         xx + 0.00000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.00000001d0*(a5-xx),
     >         xx + 0.0000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.0000001d0*(a5-xx),
     >         xx + 0.000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.000001d0*(a5-xx),
     >         xx + 0.00001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.00001d0*(a5-xx),
     >         xx + 0.0001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.0001d0*(a5-xx),
     >                 xx + 0.001d0*(a5-xx),128, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.001d0*(a5-xx),
     >                   xx + 0.01d0*(a5-xx),256, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.01d0*(a5-xx),
     >                   xx + 0.03d0*(a5-xx),256, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.03d0*(a5-xx),
     >                   xx + 0.1d0*(a5-xx),256, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.1d0*(a5-xx),
     >                   xx + 0.9d0*(a5-xx),4096, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.9d0*(a5-xx),
     >             xx + 0.99d0*(a5-xx),1024, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.99d0*(a5-xx),
     >             xx + 0.999d0*(a5-xx),256, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.999d0*(a5-xx),
     >             xx + 0.9999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.9999d0*(a5-xx),
     >             xx + 0.99999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.99999d0*(a5-xx),
     >             xx + 0.999999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.999999d0*(a5-xx),
     >             xx + 0.9999999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.9999999d0*(a5-xx),
     >             xx + 0.99999999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.99999999d0*(a5-xx),
     >             xx + 0.999999999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQg0int,xx + 0.999999999d0*(a5-xx),a5,32, 
     >	ErrEst, Ier)

        print *,'F2HQg0int tmp candidate 11',
     >  100.d0*abs((tmp1 - tmp)/tmp), 'percent',tmp


      tmp = ADZINT_CT
     >  (F2HQg0int, XX, zmax, AERR,RERR,ERREST,IER, IACTA, IACTB)


      F2HQg0= als4pi * efac(Nptn)/9d0 * tmp
      end

      Real*8 Function F2HQg0int(z)
      implicit none
      Real*8, intent(in):: z

      Real*8 :: HQcoef0 !functions
      Real*8 :: PDFH !functions
      Real*8 :: Coef, gluon
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer :: Ihadn
      Common /F2LHQhad/ Ihadn

      Coef=HQcoef0(2,z,QQ,bmu,bms)
      gluon=PDFH(Ihadn, 0, XX/z, bmu)
cmarco the multiplication by XX is made outside in 4fit/chidis.f
      F2HQg0int=gluon*Coef/z !*XX
cmarco this has been checked against JS code
      End

      Real*8 Function F2HQSub0(X, Q, aMu, aMs, Ihad, Nptn)
      Use mDisPac
      implicit none
	real*8 :: tmp, tmp1, smpsnf, a5

      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad, Nptn

      Real*8 :: ADZINT_CT, alpi, F2HQS0int, fchi !functions
      external F2HQS0int
      Real*8, Parameter :: D1=1D0
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST
      Integer :: IER

      Real*8 :: efac(6)=(/4d0,1d0,1D0,4D0,1D0,4D0/)  ! e^2*9 for each quark
      Real*8 :: XX, QQ, bmu, bms, eps, aMu0, als4pi, zmax, Xi
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      bms=aMs
      QQ=Q
      XX=X
cpn2010 Xi=fchi(x,Q,bms,bms)
      Xi=fchi(x,Q,bms,bms)!X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      If(Xi.ge.1d0) then
         F2HQSub0=0d0
         return
      Endif
      Ihadn=Ihad

      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

c ******************************************************************
c force an error to reexamine this -- so far, this code is not used.
c ******************************************************************
c force an error, to reexamine this if necessary
	if(bmu*bmu .ge. 0.d0) then
	   print *,'fatal error -- F2HQS0int not yet optimized'
	   stop
	endif


      tmp = ADZINT_CT
     >  (F2HQS0int, XX, zmax, AERR/1000.d0,RERR/1000.d0,
     >	ERREST,IER, IACTA, IACTB)

	a5 = zmax - 0.00000000001d0*(zmax-xx)

	tmp1 = 
     >   SMPSNF(F2HQS0int,xx + 0.00000000001d0*(a5-xx), 
     >	       xx + 0.00000000001d0*(a5-xx),64, ErrEst, Ier)
     > + SMPSNF(F2HQS0int,xx + 0.00000000001d0*(a5-xx),
     >         xx + 0.0000000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.0000000001d0*(a5-xx),
     >         xx + 0.000000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.000000001d0*(a5-xx),
     >         xx + 0.00000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.00000001d0*(a5-xx),
     >         xx + 0.0000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.0000001d0*(a5-xx),
     >         xx + 0.000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.000001d0*(a5-xx),
     >         xx + 0.00001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.00001d0*(a5-xx),
     >         xx + 0.0001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.0001d0*(a5-xx),
     >                 xx + 0.001d0*(a5-xx),128, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.001d0*(a5-xx),
     >                   xx + 0.01d0*(a5-xx),256, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.01d0*(a5-xx),
     >                   xx + 0.03d0*(a5-xx),256, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.03d0*(a5-xx),
     >                   xx + 0.1d0*(a5-xx),256, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.1d0*(a5-xx),
     >                   xx + 0.9d0*(a5-xx),4096, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.9d0*(a5-xx),
     >             xx + 0.99d0*(a5-xx),1024, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.99d0*(a5-xx),
     >             xx + 0.999d0*(a5-xx),256, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.999d0*(a5-xx),
     >             xx + 0.9999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.9999d0*(a5-xx),
     >             xx + 0.99999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.99999d0*(a5-xx),
     >             xx + 0.999999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.999999d0*(a5-xx),
     >             xx + 0.9999999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.9999999d0*(a5-xx),
     >             xx + 0.99999999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.99999999d0*(a5-xx),
     >             xx + 0.999999999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(F2HQS0int,xx + 0.999999999d0*(a5-xx),a5,32, 
     >	ErrEst, Ier)

        print *,'F2HQS0int tmp candidate 11',
     >  100.d0*abs((tmp1 - tmp)/tmp), 'percent',tmp


      tmp = ADZINT_CT
     >  (F2HQS0int, XX, zmax, AERR,RERR,ERREST,IER, IACTA, IACTB)

      F2HQSub0= als4pi * efac(Nptn)/9d0 * tmp

      end

      Real*8 Function F2HQS0int(z)
      implicit none
      Real*8, intent(in):: z

      Real*8 :: A1Hg, PDFH !functions
      Real*8 :: Sub, gluon
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer :: Ihadn
      Common /F2LHQhad/ Ihadn

      Sub=A1Hg(z/zmax,eps)
      gluon=PDFH(Ihadn, 0, XX/z, bmu)
      F2HQS0int=gluon*Sub/z

      End

      Real*8 Function FLhg2(X, Q, aMu, aMs, Ihad)
      Use mDisPac
      implicit none
	real*8 :: tmp, tmp1, smpsnf, a5

      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad

      Real*8 :: ADZINT_CT, alpi, FLHQgint !functions
      external FLHQgint
      Real*8, Parameter :: D1=1D0
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST
      Integer :: IER
      Real*8 :: aMu0, als4pi, Xi
      Real*8 :: XX, QQ, bmu, bms, eps, zmax,tFLhg2
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      bms=aMs
      QQ=Q
      XX=X
      Xi=X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      If(Xi.ge.1d0) then
         FLhg2=0d0
         return
      Endif
      Ihadn=Ihad

      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

	a5 = zmax - 0.00000000001d0*(zmax-xx)

      tFLhg2 = 
     >   SMPSNF(FLHQgint,xx + 0.00000000001d0*(a5-xx), 
     >	                 xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(FLHQgint,xx + 0.000001d0*(a5-xx),
     >                   xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQgint,xx + 0.00001d0*(a5-xx),
     >                   xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQgint,xx + 0.0001d0*(a5-xx),
     >                   xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQgint,xx + 0.001d0*(a5-xx),
     >                   xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQgint,xx + 0.01d0*(a5-xx),
     >                   xx + 0.1d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQgint,xx + 0.1d0*(a5-xx),
     >                   xx + 0.9d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(FLHQgint,xx + 0.9d0*(a5-xx),
     >                   xx + 0.99d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQgint,xx + 0.99d0*(a5-xx),
     >                   xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQgint,xx + 0.999d0*(a5-xx),
     >                   xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQgint,xx + 0.9999d0*(a5-xx),
     >                   xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQgint,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx      tFLhg2 = ADZINT_CT
c xx     >  (FLHQgint, XX, zmax, AERR,RERR,ERREST ,IER, IACTA, IACTB)


cmarco Divide by 2, to match the normalization of CTEQ's F_{long}
      FLhg2= (als4pi)**2 * tFLhg2/2d0
                                                                 
      return 
      End                       !FLHQg

      Real*8 Function FLHQgint(z)
      Use mDisPac
      implicit none
      Real*8, intent(in):: z

      Real*8 :: HQcoef, PDFH !functions
      Real*8, Parameter :: D1=1D0
      Integer, Parameter :: IACTA=1, IACTB=1
      Real*8 :: Coef, gluon
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn

      Coef= HQcoef(0,z,QQ,bmu,bms)
      gluon=PDFH(Ihadn, 0, XX/z, bmu)

      FLHQgint=gluon*Coef/z !XX the multiplication by XX is made outside in 4fit/chidis.f
      End




      Real*8 Function FLHQSub2(X, Q, aMu, aMs, Ihad)
      Use mDisPac
      implicit none
	   double precision xx1, xx2, xii, fii, tmp2, tmp3
	   double precision xx3, sum
	   integer :: ieasy
	   integer :: itst, ii
	double precision tmp, tmp1, smpsnf, a5

      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad
      Real*8 :: ADZINT_CT, alpi, FLHQSint, fchi !functions
      external FLHQSint
      Real*8, Parameter :: D1=1D0
      Integer, Parameter :: IACTA=0, IACTB=1
      Real*8 :: ERREST
      Integer :: IER

      Real*8 :: XX, QQ, bmu, bms, eps, aMu0, als4pi, 
     >  zmax, Xi, tFLSub2
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      bms=aMs
      QQ=Q
      XX=X

      Xi=fchi(x,Q,bms,bms)!X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      If(Xi.ge.1d0) then
         FLHQSub2=0d0
         return
      Endif
      Ihadn=Ihad

      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

	a5 = zmax - 0.00000000001d0*(zmax-xx)

      tFLSub2 = 
     >   SMPSNF(FLHQSint,xx + 0.00000000001d0*(a5-xx), 
     >	                 xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(FLHQSint,xx + 0.000001d0*(a5-xx),
     >                   xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQSint,xx + 0.00001d0*(a5-xx),
     >                   xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQSint,xx + 0.0001d0*(a5-xx),
     >                   xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQSint,xx + 0.001d0*(a5-xx),
     >                   xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQSint,xx + 0.01d0*(a5-xx),
     >                   xx + 0.1d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQSint,xx + 0.1d0*(a5-xx),
     >                   xx + 0.9d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(FLHQSint,xx + 0.9d0*(a5-xx),
     >                   xx + 0.99d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQSint,xx + 0.99d0*(a5-xx),
     >                   xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQSint,xx + 0.999d0*(a5-xx),
     >                   xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQSint,xx + 0.9999d0*(a5-xx),
     >                   xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLHQSint,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx      tFLSub2= ADZINT_CT
c xx  >  (FLHQSint, XX,zmax, AERR,RERR,ERREST,IER, IACTA, IACTB)

cmarco Divide by 2, to match the normalization of CTEQ's F_{long}
      FLHQSub2= (als4pi)**2 *tFLSub2/2d0

      return 
      End !FLHQsub

      Real*8 Function FLHQSint(z)

      Use mDisPac
      implicit none
	real*8 :: tmp, tmp1, smpsn2, a5
      Real*8, intent(in):: z

      Real*8 :: adz2nt, FLsub1, PDFH !functions
      external FLsub1
      Real*8, Parameter :: D1=1D0
      Integer, Parameter :: IACTA=1, IACTB=1
      Integer :: IER
      Real*8 :: Sub, XX, QQ, bmu, bms, eps, zmax, errest
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Real*8 :: xz, gluon
      common /FLsub/ xz
      Integer Ihadn
      Common /F2LHQhad/ Ihadn

      xz=z

	a5 = D1 - 0.00000000001d0*(D1-z)

      Sub = 
     >   SMPSN2(FLsub1,z + 0.00000000001d0*(a5-z),a5,16, ErrEst, Ier)

c xx      Sub= Adz2nt
c xx     >        (FLsub1, z,D1, AERR,RERR,ERREST,IER, IACTA, IACTB)

      gluon=PDFH(Ihadn, 0, XX/z, bmu)

      FLHQSint=gluon*Sub/z

      End
      Real*8 Function FLsub1(z) !A^1 (*) C^1
      implicit none

      Real*8, intent(in):: z

      Real*8 :: A1Hg, clnp1a !functions
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Real*8 :: xz
      common /FLsub/ xz

      FLsub1=A1Hg(xz/z,eps)*clnp1a(z,0)/z



      End

      Function clnp1a(y,Nf)
c LO non-singlet coefficient functions for F_L
      Implicit none
      Real*8 :: clnp1a, cls1aL, clg1aL
      Real*8 :: y
      Integer :: Nf

      Real*8, parameter:: Cf=4d0/3D0

      clnp1a= 4d0* Cf * y
      Return

      Entry cls1aL(y, Nf)
c LO pure-singlet coefficient functions for F_L
      cls1aL= 0d0
      Return

      Entry clg1aL(y, Nf)
c LO gluon coefficient functions for F_L
      clg1aL= 8d0 * Nf * y * (1d0 - y)
      Return

      End




      Real*8 Function FLHQg0(X, Q, aMu, aMs, Ihad, Nptn)
      Use mDisPac
      implicit none
	real*8 :: tmp, tmp1, smpsnf, a5

      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad, Nptn

      Real*8 :: ADZINT_CT, alpi, FLHQg0int !functions
      external FLHQg0int
      Real*8, Parameter :: D1=1D0
      Integer, Parameter :: IACTA=1, IACTB=1
      Real*8 :: ERREST
      Integer :: IER

      Real*8 :: efac(6)=(/4d0,1d0,1D0,4D0,1D0,4D0/)  ! e^2*9 for each quark
      Real*8 :: XX, QQ, bmu, bms, eps, aMu0, als4pi, zmax, Xi
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      bms=aMs
      QQ=Q
      XX=X
      Xi=X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      If(Xi.ge.1d0) then
         FLHQg0=0d0
         return
      Endif
      Ihadn=Ihad

      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

      tmp = ADZINT_CT
     >  (FLHQg0int, XX,zmax, AERR/1000.d0,RERR/1000.d0,
     >  ERREST,IER, IACTA, IACTB)

	a5 = zmax - 0.00000000001d0*(zmax-xx)

	tmp1 = 
     >   SMPSNF(FLHQg0int,xx + 0.00000000001d0*(a5-xx), 
     >	       xx + 0.00000000001d0*(a5-xx),64, ErrEst, Ier)
     > + SMPSNF(FLHQg0int,xx + 0.00000000001d0*(a5-xx),
     >         xx + 0.0000000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.0000000001d0*(a5-xx),
     >         xx + 0.000000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.000000001d0*(a5-xx),
     >         xx + 0.00000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.00000001d0*(a5-xx),
     >         xx + 0.0000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.0000001d0*(a5-xx),
     >         xx + 0.000001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.000001d0*(a5-xx),
     >         xx + 0.00001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.00001d0*(a5-xx),
     >         xx + 0.0001d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.0001d0*(a5-xx),
     >                 xx + 0.001d0*(a5-xx),128, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.001d0*(a5-xx),
     >                   xx + 0.01d0*(a5-xx),256, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.01d0*(a5-xx),
     >                   xx + 0.03d0*(a5-xx),256, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.03d0*(a5-xx),
     >                   xx + 0.1d0*(a5-xx),256, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.1d0*(a5-xx),
     >                   xx + 0.9d0*(a5-xx),4096, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.9d0*(a5-xx),
     >             xx + 0.99d0*(a5-xx),1024, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.99d0*(a5-xx),
     >             xx + 0.999d0*(a5-xx),256, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.999d0*(a5-xx),
     >             xx + 0.9999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.9999d0*(a5-xx),
     >             xx + 0.99999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.99999d0*(a5-xx),
     >             xx + 0.999999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.999999d0*(a5-xx),
     >             xx + 0.9999999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.9999999d0*(a5-xx),
     >             xx + 0.99999999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.99999999d0*(a5-xx),
     >             xx + 0.999999999d0*(a5-xx),64, ErrEst, Ier) 
     > + SMPSNF(FLHQg0int,xx + 0.999999999d0*(a5-xx),a5,32, 
     >	ErrEst, Ier)

        print *,'FLHQg0int tmp candidate 11',
     >  100.d0*abs((tmp1 - tmp)/tmp), 'percent',tmp

c force an error, to reexamine this if necessary
	if(tmp*tmp .ge. 0.d0) then
	   print *,'fatal error -- FLHQg0int not yet optimized'
	   stop
	endif


      tmp = ADZINT_CT
     >  (FLHQg0int, XX,zmax, AERR,RERR,ERREST,IER, IACTA, IACTB)

      FLHQg0= als4pi * efac(Nptn)/9d0 * tmp

      end

      Real*8 Function FLHQg0int(z)
      implicit none
      Real*8, intent(in):: z

      Real*8 :: HQcoef0, PDFH !functions
      Real*8 :: Coef
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Real*8 :: gluon
      Integer Ihadn
      Common /F2LHQhad/ Ihadn

      Coef= HQcoef0(0,z,QQ,bmu,bms)
      gluon=PDFH(Ihadn, 0, XX/z, bmu)
cmarco the multiplication by XX is made outside in 4fit/chidis.f
      FLHQg0int=gluon*Coef/z !*XX

      End



cmarco=================================================================
c Implementation of Eq 23  for FL with no plus distrb. 
cmarco=================================================================


      Real*8 Function FLLNS2(X, Q, aMu, aMs, Ihad, iPrtn)
      Use mDisPac
      Implicit none
	real*8 :: tmp, tmp1, smpsnf, a5

      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad, iPrtn
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn
 
      Real*8 :: ADZINT_CT, alpi, FLLNS2int, cL1llint !functions
      external FLLNS2int, cL1llint
      Real*8 :: FLLNS2int0
      external FLLNS2int0
      Real*8, Parameter :: D1=1D0
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST
      Integer :: IER
      Real*8 :: Xi, tFLLNS2, tcL1ll, xxi, aMu0, als4pi
      Real*8 :: XX, QQ, bmu, bms, eps, zmax,  small
      parameter (small=1d-7)
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      jPrtn=iPrtn

      bms=aMs
      QQ=Q
      XX=X
      Xi=X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      xxi=QQ**2/bms**2

      If(Xi.ge.1d0) then
         FLLNS2=0d0
         return
      Endif
      Ihadn=Ihad

      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

	a5 = zmax - 0.00000000001d0*(zmax-xx)

      tFLLNS2 = 
     >   SMPSNF(FLLNS2int,xx + 0.00000000001d0*(a5-xx), 
     >	                  xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(FLLNS2int,xx + 0.000001d0*(a5-xx),
     >                    xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLLNS2int,xx + 0.00001d0*(a5-xx),
     >                    xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLLNS2int,xx + 0.0001d0*(a5-xx),
     >                    xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLLNS2int,xx + 0.001d0*(a5-xx),
     >                    xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLLNS2int,xx + 0.01d0*(a5-xx),
     >                    xx + 0.1d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLLNS2int,xx + 0.1d0*(a5-xx),
     >                    xx + 0.9d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(FLLNS2int,xx + 0.9d0*(a5-xx),
     >                    xx + 0.99d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLLNS2int,xx + 0.99d0*(a5-xx),
     >                    xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLLNS2int,xx + 0.999d0*(a5-xx),
     >                    xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLLNS2int,xx + 0.9999d0*(a5-xx),
     >                    xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(FLLNS2int,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx      tFLLNS2= ADZINT_CT
c xx     >  (FLLNS2int, XX, zmax, AERR,RERR,ERREST,IER, IACTA, IACTB) !we are integrating up to zmax


c      t1FLLNS2= ADZINT_CT          !piece coming from the plus distribution
c     >  (FLLNS2int0, small, XX, AERR,RERR,ERREST,IER, IACTA, IACTB) 
      
	a5 = D1 - 0.00000000001d0*(D1-xx)

      tcL1ll = 
     >   SMPSNF(cL1llint,xx + 0.00000000001d0*(a5-xx), 
     >	                 xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(cL1llint,xx + 0.000001d0*(a5-xx),
     >                   xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(cL1llint,xx + 0.00001d0*(a5-xx),
     >                   xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(cL1llint,xx + 0.0001d0*(a5-xx),
     >                   xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(cL1llint,xx + 0.001d0*(a5-xx),
     >                   xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(cL1llint,xx + 0.01d0*(a5-xx),
     >                   xx + 0.1d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(cL1llint,xx + 0.1d0*(a5-xx),
     >                   xx + 0.9d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(cL1llint,xx + 0.9d0*(a5-xx),
     >                   xx + 0.99d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(cL1llint,xx + 0.99d0*(a5-xx),
     >                   xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(cL1llint,xx + 0.999d0*(a5-xx),
     >                   xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(cL1llint,xx + 0.9999d0*(a5-xx),
     >                   xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(cL1llint,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx      tcL1ll = ADZINT_CT
c xx     >  (cL1llint, XX, D1, AERR,RERR,ERREST,IER, IACTA, IACTB) 

      tcL1ll = Tf*2d0 * 2d0/3d0 * dlog(xxi) * tcL1ll !the extra factor of 2 comes from
                                                      !our notebook computation

      FLLNS2= (als4pi)**2 * (tFLLNS2 - tcL1ll)

      end



      Real*8 Function FLLNS2int(z)
      implicit none
      Real*8, intent(in):: z
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn
 
      Real*8 :: LLNSqr, PDF_CT !functions
      Real*8 :: Coef, xxi
      Real*8 :: XX, QQ, bmu, bms, eps, zmax, suml
c      real*8:: suml0
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn, Ir
      Common /F2LHQhad/ Ihadn
      external PDF_CT

      xxi=QQ**2/bms**2
      
      Coef= LLNSqr(z,xxi)
      
      suml=PDF_CT(ihadn, jprtn, XX/z, bmu, Ir) +
     >     PDF_CT(ihadn,-jprtn, XX/z, bmu, Ir)

c      suml0=PDF_CT(ihadn, jprtn, XX, bmu, Ir)+
c     >      PDF_CT(ihadn,-jprtn, XX, bmu, Ir)
c      FLLNS2int=(suml - z*suml0)*Coef/z !

      FLLNS2int=suml*Coef/z !*XX the multiplication by XX is made outside in 4fit/chidis.f

      End


      Real*8 Function FLLNS2int0(z)! coming from the plus 
      implicit none
      Real*8, intent(in):: z
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn

      Real*8 :: LLNSqr, PDF_CT !functions
      Real*8 :: Coef, xxi
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Real*8 :: suml0
      Integer Ihadn, Ir
      Common /F2LHQhad/ Ihadn
      external PDF_CT

      xxi=QQ**2/bms**2
      
      Coef= LLNSqr(z,xxi)
      
      suml0=PDF_CT(ihadn, jprtn, XX, bmu,  Ir) + 
     >      PDF_CT(ihadn,-jprtn, XX, bmu, Ir)

      FLLNS2int0= -suml0*Coef !*XX the multiplication by XX is made outside in 4fit/chidis.f

      End


      Real*8 Function cL1llint(z)
      implicit none
      Real*8, intent(in):: z

      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn

      Real*8, Parameter:: Cf=4D0/3d0, Tf=0.5d0
      Real*8 :: PDF_CT !functions
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Real*8 :: suml
      Integer Ihadn, Ir
      Common /F2LHQhad/ Ihadn
      
      suml=PDF_CT(ihadn, jprtn, XX/z, bmu, Ir) + 
     >     PDF_CT(ihadn,-jprtn, XX/z, bmu,Ir)

      cL1llint=4d0*Cf*z*suml/z

      End


cmarco=================================================================
c Implementation of the plus distributuion in Eq 23  for F2
cmarco=================================================================


      Real*8 Function F2LNS2(X, Q, aMu, aMs, Ihad, iPrtn)
      Use mDisPac
      Implicit none
	real*8 :: tmp, tmp1, smpsnf, a5

c      external c1ll
      Real*8, intent(in) :: X, Q, aMu, aMs
      Integer, intent(in):: Ihad, iPrtn
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn
 
      Real*8 :: ADZINT_CT, alpi, PDF_CT, F2LNS2int, F2LNS2int0, c1llint!functions
      external F2LNS2int, F2LNS2int0, c1llint
      Real*8, Parameter :: D1=1D0 
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST
      Integer :: IER
      Real*8 :: Xi, t0, t1, t2, tc1ll, aMu0, als4pi, suml0, 
     >  Finite, xxi, small
      parameter (small=1d-7)
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      jPrtn=iPrtn
      bms=aMs
      QQ=Q
      XX=X
      Xi=X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
      zmax=X/Xi
      xxi=QQ**2/bms**2

      If(Xi.ge.1d0) then
         F2LNS2=0d0
         return
      Endif
      Ihadn=Ihad

      bmu=amu
      eps=(bms/bmu)**2 !m^2/mu^2
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

      suml0=PDF_CT(ihadn, jprtn, XX, bmu, Ier) + 
     >  PDF_CT(ihadn,-jprtn, XX, bmu,Ier)

	a5 = zmax - 0.00000000001d0*(zmax-xx)

      t0 = 
     >   SMPSNF(F2LNS2int,xx + 0.00000000001d0*(a5-xx), 
     >	                  xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(F2LNS2int,xx + 0.000001d0*(a5-xx),
     >                    xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2LNS2int,xx + 0.00001d0*(a5-xx),
     >                    xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2LNS2int,xx + 0.0001d0*(a5-xx),
     >                    xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2LNS2int,xx + 0.001d0*(a5-xx),
     >                    xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2LNS2int,xx + 0.01d0*(a5-xx),
     >                    xx + 0.1d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2LNS2int,xx + 0.1d0*(a5-xx),
     >                    xx + 0.9d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(F2LNS2int,xx + 0.9d0*(a5-xx),
     >                    xx + 0.99d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2LNS2int,xx + 0.99d0*(a5-xx),
     >                    xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2LNS2int,xx + 0.999d0*(a5-xx),
     >                    xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2LNS2int,xx + 0.9999d0*(a5-xx),
     >                    xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(F2LNS2int,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx      t0= ADZINT_CT
c xx     >  (F2LNS2int, XX, zmax, AERR,RERR,ERREST,IER, IACTA, IACTB) !we are integrating up to zmax


cmarco An integral from 0 to x, coming from the plus distribution for L2NSq.
c To improve convergence of this integral, we subtract the Log(z) terms
c arising in the asymptotic limit z->0 and evaluate them analytically as
c t2

	a5 = xx - 0.00000000001d0*(xx-small)
      t1= 
     > + SMPSNF(F2LNS2int0,small,a5,16, ErrEst, Ier)

c xx      t1 = ADZINT_CT 
c xx     >  (F2LNS2int0, small, XX, AERR,RERR,ERREST,IER, IACTA, IACTB) 

      t1 = suml0*t1  !multiplication by the PDF(XX) is made outside of the subroutine F2LNS2int0
      
C      t2=-suml0*(8d0/27d0*XX*(-11D0 + 9D0*dlog(XX)**2  
C     >  + dlog(XX)*(11D0 - 12d0*dlog(xxi))-dlog(xxi)+3d0*dlog(xxi)**2)) 

      t2=-suml0*(-4d0*XX*(-157d0 + 12d0*Pi**2 - 54d0*dlog(XX)**2 +
     >  6d0*dlog(xxi) - 18d0*dlog(xxi)**2 
     >  + dlog(XX)*(-66d0 + 72d0*dlog(xxi))))/81d0

cmarco tc1ll + Finite -> contributions from alpha_s piece in Eq.B2 NPB383 1992 Van Neerven et. al.
c finite -> Analytic piece coming by tc1ll plus distrib.
      Finite=Cf*suml0*(2d0*dlog(1d0-XX)**2 - 3d0*dlog(1d0-XX) - 9d0  
     >  -4d0*Pi**2/6d0)

	a5 = D1 - 0.00000000001d0*(D1-xx)

      tc1ll = 
     >   SMPSNF(c1llint,xx + 0.00000000001d0*(a5-xx), 
     >	                xx + 0.000001d0*(a5-xx),4, ErrEst, Ier)
     > + SMPSNF(c1llint,xx + 0.000001d0*(a5-xx),
     >                  xx + 0.00001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(c1llint,xx + 0.00001d0*(a5-xx),
     >                  xx + 0.0001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(c1llint,xx + 0.0001d0*(a5-xx),
     >                  xx + 0.001d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(c1llint,xx + 0.001d0*(a5-xx),
     >                  xx + 0.01d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(c1llint,xx + 0.01d0*(a5-xx),
     >                  xx + 0.1d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(c1llint,xx + 0.1d0*(a5-xx),
     >                  xx + 0.9d0*(a5-xx),8, ErrEst, Ier) 
     > + SMPSNF(c1llint,xx + 0.9d0*(a5-xx),
     >                  xx + 0.99d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(c1llint,xx + 0.99d0*(a5-xx),
     >                  xx + 0.999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(c1llint,xx + 0.999d0*(a5-xx),
     >                  xx + 0.9999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(c1llint,xx + 0.9999d0*(a5-xx),
     >                  xx + 0.99999d0*(a5-xx),4, ErrEst, Ier) 
     > + SMPSNF(c1llint,xx + 0.99999d0*(a5-xx),a5,4, ErrEst, Ier)

c xx      tc1ll = ADZINT_CT
c xx     >  (c1llint, XX, D1, AERR,RERR,ERREST,IER, IACTA, IACTB) 

      tc1ll = Tf*2d0/3d0*dlog(QQ**2/bms**2)*(tc1ll + Finite)*2d0 !the extra factor of 2 comes from
                                                                 !(a2=4)/2 from our notebook computation

      F2LNS2= (als4pi)**2 * (t0 + t1 + t2 - tc1ll)

      end


      Real*8 Function F2LNS2int(z)
      implicit none
      Real*8, intent(in):: z
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn
 
      Real*8 :: L2NSqr, PDF_CT !functions
      Real*8 :: Coef, xxi
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Real*8 :: suml, suml0
      Integer Ihadn, Ir
      Common /F2LHQhad/ Ihadn
      
      xxi=QQ**2/bms**2

      Coef= L2NSqr(z,xxi)

      suml=PDF_CT(ihadn,jprtn, XX/z, bmu, Ir) 
     >  + PDF_CT(ihadn,-jprtn, XX/z, bmu, Ir)

      suml0=PDF_CT(ihadn,jprtn, XX, bmu, Ir) 
     >  + PDF_CT(ihadn,-jprtn, XX, bmu, Ir)

      F2LNS2int=(suml/z - suml0)*Coef !*XX the multiplication by XX is made outside in 4fit/chidis.f

      End



      Real*8 Function F2LNS2int0(z)
      Use mDisPac
      implicit none
      Real*8, intent(in):: z
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn

      Real*8 :: L2NSqr0, PDF_CT !functions
      Real*8 :: Coef, xxi
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Real*8 :: suml0, asym
      Integer Ihadn, Ir
      Common /F2LHQhad/ Ihadn
      
      xxi=QQ**2/bms**2

      Coef= L2NSqr0(z,xxi) !the logarithmic asympototic behavior
                           !is subtracted from L2NSqr0(z,xxi) in GmJSa.F ;    
                           !it is then added after analytic integration in F2LNS2.
 
      F2LNS2int0= -Coef !*XX the multiplication by XX is made outside in 4fit/chidis.f
      End

      Real*8 Function c1llint(z)
      implicit none
      Real*8, intent(in):: z

      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn

      Real*8, Parameter:: Cf=4D0/3d0, Tf=0.5d0
      Real*8 :: PDF_CT !functions
      Real*8 :: P1, P2, P3
      Real*8 :: XX, QQ, bmu, bms, eps, zmax
      common /F2LHQ/ XX, QQ, bmu, bms, eps, zmax
      Real*8 :: suml, suml0
      Integer Ihadn, Ir
      Common /F2LHQhad/ Ihadn
      
      suml=PDF_CT(ihadn, jprtn, XX/z, bmu, Ir) + 
     >     PDF_CT(ihadn,-jprtn, XX/z, bmu, Ir)

      suml0=PDF_CT(ihadn, jprtn, XX, bmu, Ir)+ 
     >      PDF_CT(ihadn,-jprtn, XX, bmu, Ir)

      P1 = 4d0*dlog(1d0-z)/(1d0-z)*(suml/z -suml0)
      P2 = -3d0/(1d0-z)*(suml/z -suml0)
      P3 = (- 2*(1d0+z)*dlog(1d0-z) - 2d0*(1d0 + z**2)/(1d0 - z)*dlog(z)
     >  + 6d0 + 4d0*z)*suml/z

      c1llint=Cf*(P1+P2+P3) 

      End
