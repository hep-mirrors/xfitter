cmarco=================================================================
c Separate implementation of ZM c2NS2 from A. Vogt hep-ph/0504242 
c This has been extrapolated from Liang's ZmVVA.f file
cmarco=================================================================


      Real*8 Function F2ZMNS2(X, Q, amu, Ihad, iPrtn, Nptn1)
      Use mDisPac
      Implicit None
	double precision smpsnf, err
	double precision a5, tmp, sum
	double precision tmp1, tmp2, ran1, hwrgen
	integer inx

      Real*8, intent(in) :: X, Q, aMu
      Integer, intent(in):: Ihad, Nptn1, iPrtn
      external F2CnsQns2
      Real*8 :: aLm
      Real*8 :: ADZINT_CT, alpi, F2CnsQns2, PDF_CT !functions
      Real*8 :: c2mnp2Dk1a, c2np2Dk0, c2mnp2Dk2a   !functions
      Real*8, Parameter :: D1=1D0 
      Integer, Parameter :: IACTA=2, IACTB=2
C      Real*8, Parameter:: Cf=4d0/3D0, Tf=.5d0, 
      Real*8, Parameter::Ca=3d0
      Real*8, parameter:: Zeta2=1.6449340668482264365D0
      Real*8, parameter:: Zeta3=1.2020569031595942854D0
      Real*8 :: ERREST, tF2ZMNS2, Delta
      Integer :: IER, Nff, Ip
      Real*8 :: aMu0, als4pi, XX, QQ, bmu, Qns0, DL1
      integer :: Ihadn, Ir
      common /F2LZM/ XX, QQ, bmu
      Common /F2LHQhad/ Ihadn

      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn
      Real*8 :: Dkfac(0:6)
      Common /DkX/ Dkfac

      Data aMu0 /-1d0/
      save aMu0

      jNf = Nptn1 ; jPrtn = iPrtn

      QQ=Q
      XX=X
      If(XX.ge.1d0) then
        tF2ZMNS2=0d0
        return
      Endif
      
      Ihadn=Ihad
      
      bmu=amu
      
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      
      aMu0=aMu
      
      aLm=dlog(QQ**2/bmu**2)
      
      Nff=jNf
      
      Delta=0d0
      
      DL1=dlog(1d0-XX)
      Do Ip=0,6
         Dkfac(Ip)=DL1**(Ip+1)/(Ip+1)
      Enddo

      Qns0 = PDF_CT(ihadn, jprtn, XX, bmu, Ir) +
     >  PDF_CT(ihadn, -jprtn, XX, bmu, Ir)
      
      Delta = (-338.513 + 46.8531*Nff
     >  + c2np2Dk0(Dkfac,Nff))*Qns0 !it recollects all the 
                                    !Dk terms coming from delta functions 
                                    !and plus distributions
      
      If(QQ.ne.bmu) then
      Delta = Delta + (
     >  aLm*aLm*(Cf*Cf*(-8d0*Zeta2+4.5d0)-5.5d0*Ca*Cf+Nff*Cf
     >  + c2mnp2Dk2a(Dkfac,Nff))
     >  +aLm*(Cf*Cf*(40d0*Zeta3-12d0*Zeta2-25.5d0)
     >  +Ca*Cf*(-12d0*Zeta3+88d0/3d0*Zeta2+215d0/6d0)
     >  +Nff*Cf*(-16d0*Zeta2-19d0)/3d0
     >  +c2mnp2Dk1a(Dkfac,Nff)))*Qns0
      Endif

	a5 = 1.d0 - 0.00000000001d0*(1.d0-XX)

	sum = 
     >   SMPSNF(F2CnsQns2,XX + 0.00000000001d0*(a5-XX), 
     >	                  XX + 0.0000001d0*(a5-XX),4, ERR, IER)
     > + SMPSNF(F2CnsQns2,XX + 0.0000001d0*(a5-XX),
     >                    XX + 0.000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.000001d0*(a5-XX),
     >                    XX + 0.00001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.00001d0*(a5-XX),
     >                    XX + 0.0001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.0001d0*(a5-XX),
     >                    XX + 0.001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.001d0*(a5-XX),
     >                    XX + 0.01d0*(a5-XX),8, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.01d0*(a5-XX),
     >                    XX + 0.1d0*(a5-XX),16, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.1d0*(a5-XX),
     >                    XX + 0.9d0*(a5-XX),32, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.9d0*(a5-XX),
     >                    XX + 0.99d0*(a5-XX),16, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.99d0*(a5-XX),
     >                    XX + 0.999d0*(a5-XX),8, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.999d0*(a5-XX),
     >                    XX + 0.9999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.9999d0*(a5-XX),
     >                    XX + 0.99999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.99999d0*(a5-XX),
     >                    XX + 0.999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.999999d0*(a5-XX),
     >                    XX + 0.9999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CnsQns2,XX + 0.9999999d0*(a5-XX),a5,4, ERR, IER)


      F2ZMNS2= (als4pi)**2 * (sum + Delta)

      end


      Function F2CnsQns2(y)
      Use mDisPac
      Implicit None
	double precision smpsnf, err
	double precision tmp, sum
	integer inx
      Real*8, intent(in) :: y
      Real*8 :: F2CnsQns2

      Real*8 :: XX, QQ, bmu, PDF_CT
      common /F2LZM/ XX, QQ, bmu
      Real*8 :: aLm
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn
      Real*8 :: Qns, Qns0

      Integer :: Nff
      Integer Ihadn, Ir
      Common /F2LHQhad/ Ihadn

      Real*8 :: CL, Dk
      Real*8 :: QGtype, c2np1a0, c2np2a0, c2np3a0  ! function
      Real*8 :: c2np1Dk, c2np2Dk, c2np3Dk
      Real*8 :: c2mnp1a0, c2mnp1Dk
      Real*8 :: c2mnp2a2, c2mnp2a1, c2mnp2Dk2, c2mnp2Dk1

      Qns = PDF_CT(ihadn, jprtn, XX/y, bmu, Ir) +
     >      PDF_CT(ihadn,-jprtn, XX/y, bmu, Ir)

      Qns0 = PDF_CT(ihadn, jprtn, XX, bmu, Ir) +
     >       PDF_CT(ihadn, -jprtn, XX, bmu, Ir)

      Nff=jNf

      CL= c2np2a0(y,Nff)
      Dk= c2np2Dk(Qns,Qns0,y,Nff)

      If(abs(QQ-bmu).gt.1e-6) then
         aLm=dlog(QQ**2/bmu**2)
         CL=CL+(aLm*(aLm*c2mnp2a2(y,NFf)+c2mnp2a1(y,Nff)))
         Dk=Dk+(aLm*(aLm*c2mnp2Dk2(Qns,Qns0,y,Nff)
     >        +c2mnp2Dk1(Qns,Qns0,y,Nff)))
      Endif

      F2CnsQns2 = CL * Qns / y + Dk

      Return
      End


      Function Dkfun(I, yy, QA, QA0)
      Implicit none
      Real*8 :: Dkfun, yy, QA, QA0
      Integer :: I

      Dkfun=(dlog(1d0-yy))**I/(1d0-yy)*(QA/yy-QA0)

      Return
      End



      Function c2np2a0(y,Nf)
c NNLO non-singlet coefficient functions for F_2
c References:
c   S. Moch, J. Vermaseren and A. Vogt, hep-ph/0504242

      Implicit none
      Real*8 :: c2np2a0, c2np2Dk, c2s2aL, c2g2aL, Dkfun,c2np2Dk0
      Real*8 :: y,Dkk(0:3)
      Integer :: Nf

      Real*8, Parameter:: Cf=4D0/3d0
      Real*8 :: DL, y1, DL1, QA, QA0

      DL = Log(y)
      y1 = 1d0 - y
      DL1 = Log(y1)

      c2np2a0= -17.74*DL1**3 + 72.24*DL1*DL1 -628.8*DL1 -181.
     >       -806.7*y + .719*y*DL**4 +DL*DL1*(37.75*DL-147.1*DL1)
     >       -28.384*DL -20.7*DL*DL -80./27.*DL**3
     >   +Nf*(-1.5*DL1*DL1 +24.87*DL1 -7.8109 -17.82*y -12.97*y*y
     >        -.185*y*DL**3 +8.113*DL*DL1 +16./3.*DL +20./9.*DL*DL)


      Return

      Entry c2np2Dk(QA, QA0, y, Nf)

      Dkk(0)=Dkfun(0, y ,QA, QA0)
      Dkk(1)=Dkfun(1, y ,QA, QA0)
      Dkk(2)=Dkfun(2, y ,QA, QA0)
      Dkk(3)=Dkfun(3, y ,QA, QA0)

      c2np2Dk = c2np2Dk0(Dkk,Nf)

      Return

      Entry c2s2aL(y, Nf)
c NNLO pure-singlet coefficient functions for F_2

      DL = Log(y)
      y1 = 1d0 - y
      DL1 = Log(y1)

      c2s2aL= Nf*((8./3.*DL1*DL1 -32./3.*DL1 +9.8937)*y1
     >          +(9.57 -13.41*y +.08*DL1**3)*y1*y1 +5.667*y*DL**3
     >          -DL*DL*DL1*(20.26-33.93*y) +43.36*y1*DL
     >          -1.053*DL*DL +40./9.*DL**3 +5.2903/y*y1*y1)
      Return

      Entry c2g2aL(y, Nf)
c NNLO gluon coefficient functions for F_2

      DL = Log(y)
      y1 = 1d0 - y
      DL1 = Log(y1)

      c2g2aL= Nf*(58./9.*DL1**3 -24.*DL1*DL1 -34.88*DL1 +30.586
     >          -(25.08 +760.3*y+29.65*DL1**3)*y1 +1204.*y*DL*DL
     >          +DL*DL1*(293.8 +711.2*y +1043.*DL) +115.6*DL
     >          -7.109*DL*DL +70./9.*DL**3 +11.9033/y*y1)
      Return

      End


      Real*8 function c2np2Dk0(Dk, Nf)
      Implicit none
      Real*8 Dk(0:*)
      Integer Nf

      c2np2Dk0 = 128./9.*Dk(3)
     >    +(-184./3. +Nf*16./9.)*Dk(2)
     >    +(-31.1052 -Nf*232./27.)*Dk(1)
     >    +(188.641 +Nf*6.34888)*Dk(0)

      Return
      end






      Real*8 Function c2mnp2a2(y,Nf)
      Implicit none
      Real*8 y
      Integer Nf
      Real*8, Parameter:: Cf=4d0/3D0, Tf=.5d0, Ca=3d0
      Real*8, parameter:: Zeta2=1.6449340668482264365D0
      Real*8 c2mnp2a1, c2ms2a2, c2ms2a1, c2mg2a2, c2mg2a1, xLi
      Real*8 DL, y1, DL1, y1p, DL1p, y2p, Li2y, Li2y1,tmp

      DL = Log(y)
      y1 = 1d0 - y
      DL1 = Log(y1)
      y1p= 1d0 + y
      c2mnp2a2=Cf*Cf*(-8d0*y1p*DL1-4d0*(1d0+y*y)/y1*DL+2d0*y1p*DL
     >                -10d0-2d0*y)
     >     +Ca*Cf*11d0/3d0*y1p + Nf*Cf*(-2d0/3d0*y1p)
      return

      Entry c2mnp2a1(y,Nf)
      DL = Log(y)
      y1 = 1d0 - y
      DL1 = Log(y1)
      y1p= 1d0 + y
      DL1p=Log(y1p)
      y2p=1+y*y
      Li2y=xLi(2,-y)               !dilogarithm
      Li2y1=xLi(2,y1)

      c2mnp2a1=Cf*Cf*(y2p/y1*(4d0*DL*DL-24d0*DL*DL1-6d0*DL)
     >  +y2p/y1p*(4d0*DL*DL-16d0*Li2y-16d0*DL*DL1p-8d0*Zeta2)
     >  +y1p*(4d0*Li2y1+4d0*DL*DL1-12d0*DL1*DL1-4d0*DL*DL+16d0*Zeta2)
     >  +8d0*(2d0+3d0*y)*DL1-2d0*(3d0+11d0*y)*DL+2d0*(19d0+14d0*y))
     > +Ca*Cf*(y2p/y1*(44d0/3d0*DL+2d0*DL*DL)
     >  +y2p/y1p*(8d0*Li2y-2d0*DL*DL+8d0*DL*DL1p+4d0*Zeta2)
     >  +y1p*(22d0/3d0*DL1+4d0*Zeta2)-(164d0+434d0*y)/9d0)
     > +Nf*Cf*(-8d0/3d0*y2p/y1*DL-4d0/3d0*y1p*DL1+(32d0+68d0*y)/9d0)

      Return
      
      Entry c2ms2a2(y,Nf)
      DL = Log(y)
      y1p= 1d0 + y
      c2ms2a2=Nf*Cf*Tf*(8d0*y1p*DL+4d0/3d0*(3d0-4d0*y*y-3d0*y+4d0/y))
      return

      Entry c2ms2a1(y,Nf)
      DL = Log(y)
      y1 = 1d0 - y
      DL1 = Log(y1)
      y1p= 1d0 + y
      Li2y1=xLi(2,y1)               !dilogarithm

      c2ms2a1=Nf*Cf*Tf*(16d0*y1p*(Li2y1+DL*DL1-DL*DL)+32d0*y*y*DL
     >  +8d0/3d0*(3d0-4d0*y*y-3d0*y+4d0/y)*DL1
     >  -16d0/9d0*(39d0+4d0*y*y-30d0*y-13d0/y))
      return

      Entry c2mg2a2(y,Nf)
      DL = Log(y)
      y1 = 1d0 - y
      DL1 = Log(y1)
      y1p= 1d0 + y
      tmp=8d0*(1d0+2d0*y*y-2d0*y)*DL1
cmarco
      c2mg2a2=Nf*(Cf*Tf*(-4d0*(1d0+4d0*y*y-2d0*y)*DL+tmp-2d0+8d0*y)
     >  +Ca*Tf*(8d0*(1d0+4d0*y)*DL+tmp !it should be Ca*Tf according to NPB383 1992 Eq.B.6 and not Ca*Cf
     >  +4d0/3d0*(3d0-31d0*y*y+24d0*y+4d0/y)))
      return

      Entry c2mg2a1(y,Nf)
      DL = Log(y)
      y1 = 1d0 - y
      DL1 = Log(y1)
      y1p= 1d0 + y
      DL1p= Log(y1p)
      Li2y1=xLi(2,y1)               !dilogarithm
      Li2y=xLi(2,-y)
cmarco
      c2mg2a1=Nf*(Cf*Tf*(16d0*(1d0+2d0*y*y-2d0*y)*(DL1*DL1-2d0*Zeta2)
     >  +8d0*(1d0+4d0*y*y-2d0*y)*DL*DL-8d0*(3d0+8d0*y*y-6d0*y)*DL*DL1
     >  +8d0*(1d0+10d0*y*y-6d0*y)*DL+8d0*(1d0-2d0*y)*Li2y1
     >  -4d0*(7d0+20d0*y*y-24d0*y)*DL1+4d0*(9d0+4d0*y*y-17d0*y))
     >  +Ca*Tf*(32d0*y*(3d0-y)*DL*DL1 !it should be Ca*Tf according to NPB383 1992 Eq.B.6 and not Ca*Cf
     >   -16d0*(1d0+2d0*y*y+2d0*y)*(DL*DL1p+Li2y)
     >   -16d0*(1d0+3d0*y)*DL*DL+16d0*(1d0+4d0*y)*Li2y1
     >   -16d0*(1d0+2d0*y*y)*Zeta2+8d0*y*(25d0*y-24d0)*DL
     >   +8d0*(1d0+2d0*y*y-2d0*y)*DL1*DL1
     >   -8d0/3d0*(3d0+67d0*y*y-60d0*y-4d0/y)*DL1
     >   -4d0/9d0*(165d0-407d0*y*y+276d0*y-52d0/y)))
      return

      end




      Real*8 Function c2mnp2Dk2(QA, QA0, y, Nf)
      Implicit none
      Real*8 QA, QA0, y
      Integer Nf
      Real*8 Dkk(0:2)
      Real*8 Dkfun, c2mnp2Dk1, c2mnp2Dk2a, c2mnp2Dk1a !function

      Dkk(0)=Dkfun(0,y,QA,QA0)
      Dkk(1)=Dkfun(1,y,QA,QA0)
      c2mnp2Dk2 = c2mnp2Dk2a(Dkk, Nf)
      return

      Entry c2mnp2Dk1(QA, QA0, y, Nf)
      Dkk(0)=Dkfun(0,y,QA,QA0)
      Dkk(1)=Dkfun(1,y,QA,QA0)
      Dkk(2)=Dkfun(2,y,QA,QA0)
      c2mnp2Dk1 = c2mnp2Dk1a(Dkk, Nf)
      Return
      End





      Real*8 function c2mnp2Dk2a(Dk, Nf)
      Implicit none
      Real*8 Dk(0:*)
      Integer Nf
      Real*8, Parameter:: Cf=4d0/3D0, Ca=3D0
      Real*8, parameter:: Zeta2=1.6449340668482264365D0
      Real*8 c2mnp2Dk1a

      c2mnp2Dk2a = Cf*Cf*(16d0*Dk(1)+12d0*Dk(0)) -Ca*Cf*22d0/3d0*Dk(0)
     >     +Nf*Cf*4d0/3d0*Dk(0)
      return

      Entry c2mnp2Dk1a(Dk, Nf)
      c2mnp2Dk1a=Cf*Cf*(24d0*Dk(2)-12d0*Dk(1)-(32d0*Zeta2+45d0)*Dk(0))
     >     +Ca*Cf*(-44d0/3d0*Dk(1)+(367d0/9d0-8d0*Zeta2)*Dk(0))
     >     +Nf*Cf*(8d0/3d0*Dk(1)-58d0/9d0*Dk(0))
      Return
      end


cmarco=================================================================
c Separate implementation of ZM cLNS2 from A. Vogt hep-ph/0504242 
c This has been extrapolated from Liang's ZmVVA.f file
cmarco=================================================================

      Real*8 Function FLZMNS2(X, Q, amu, Ihad, iPrtn, Nptn1)
      Use mDisPac
      Implicit None
	double precision smpsnf, err
	double precision a5, tmp, sum
	double precision tmp1, tmp2, ran1, hwrgen
	integer inx

      Real*8, intent(in) :: X, Q, aMu
      Integer, intent(in):: Ihad, Nptn1, iPrtn

      Real*8 :: ADZINT_CT, alpi, FLCnsQns2, PDF_CT !functions
      external FLCnsQns2
      Real*8, Parameter :: D1=1D0 
      Integer, Parameter :: IACTA=0, IACTB=1
      Real*8 :: ERREST, tFLZMNS2, Delta
      Integer :: IER
      Real*8 :: XX, QQ, bmu, aMu0, als4pi, zmax, Xi, Qns0
      common /F2LZM/ XX, QQ, bmu
      Integer Ihadn, Ir
      Common /F2LHQhad/ Ihadn
 
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn
      Data aMu0 /-1d0/
      save aMu0

      jNf = Nptn1 ; jPrtn = iPrtn

      QQ=Q
      XX=X
c      Xi=X*(1d0 + 4d0*bms*bms/(Q*Q)) !chi; rescaling variable
c      zmax=X/Xi
      If(XX.ge.1d0) then
        tFLZMNS2=0d0
        return
      Endif
      Ihadn=Ihad

      bmu=amu
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

      Qns0 = PDF_CT(ihadn, jprtn, XX, bmu, Ir) +
     >      PDF_CT(ihadn,-jprtn, XX, bmu, Ir)

      Delta= -0.012d0*Qns0      !it recollects all the 
                                !Dk terms coming from delta functions 
                                !and plus distributions

	a5 = 1.d0 - 0.00000000001d0*(1.d0-XX)

	sum = 
     >   SMPSNF(FLCnsQns2,XX + 0.00000000001d0*(a5-XX), 
     >	                  XX + 0.00000001d0*(a5-XX),4, ERR, IER)
     > + SMPSNF(FLCnsQns2,XX + 0.00000001d0*(a5-XX),
     >                    XX + 0.0000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.0000001d0*(a5-XX),
     >                    XX + 0.000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.000001d0*(a5-XX),
     >                    XX + 0.00001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.00001d0*(a5-XX),
     >                    XX + 0.0001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.0001d0*(a5-XX),
     >                    XX + 0.001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.001d0*(a5-XX),
     >                    XX + 0.01d0*(a5-XX),8, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.01d0*(a5-XX),
     >                    XX + 0.1d0*(a5-XX),16, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.1d0*(a5-XX),
     >                    XX + 0.9d0*(a5-XX),32, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.9d0*(a5-XX),
     >                    XX + 0.99d0*(a5-XX),16, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.99d0*(a5-XX),
     >                    XX + 0.999d0*(a5-XX),8, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.999d0*(a5-XX),
     >                    XX + 0.9999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.9999d0*(a5-XX),
     >                    XX + 0.99999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.99999d0*(a5-XX),
     >                    XX + 0.999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.999999d0*(a5-XX),
     >                    XX + 0.9999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCnsQns2,XX + 0.9999999d0*(a5-XX),a5,4, ERR, IER)


c include marco's  Divide by 2, to match with the normalization of F_{long}
      FLZMNS2= (als4pi)**2 * (sum + Delta) / 2.d0

      end


      Function FLCnsQns2(y)
      Implicit None

      Real*8, intent(in) :: y
      Real*8 :: FLCnsQns2, PDF_CT, XX, QQ, bmu, alm, Qns
      common /F2LZM/ XX, QQ, bmu
      Integer :: jNf, jPrtn
      Common /F2L_NF / jNf, jPrtn
      Integer :: Nff, Ihadn, Ir
      Common /F2LHQhad/ Ihadn 
      Real*8 :: CL, clnp2a_CT, clmnp2a

      Nff=jNf
      Qns = PDF_CT(ihadn, jprtn, XX/y, bmu, Ir) +
     >      PDF_CT(ihadn,-jprtn, XX/y, bmu, Ir)

      CL= clnp2a_CT(y,Nff)
      If(abs(QQ-bmu).gt.1e-6) then
         aLm=dlog(QQ**2/bmu**2)
         CL= CL + clmnp2a(y,Nff)*aLm
      endif

      FLCnsQns2 = CL * Qns / y
      Return
      End

      Function clnp2a_CT(y,Nf)
c (alpha_S^2) non-singlet coefficient functions for F_L
c References:
c   S. Moch, J. Vermaseren and A. Vogt, PLB606(2005)123-129
c   W.L. van Neerven, A. Vogt, NPB588(2000)345
c   W.L. van Neerven, A. Vogt, NPB568(2000)263
c   E.B. Zijlstra, W.L. van Neerven, PLB272(1991)127
c   S. Moch, J.A.M. Vermaseren, NPB573(2000)853
c   J.S. Guillen, et al., NPB353(1991)337

      Implicit none
      Real*8 :: clnp2a_CT, cls2aL, clg2aL
      Real*8 :: y
      Integer :: Nf

      Real*8, Parameter:: Cf=4D0/3d0, Ca=3d0
      Real*8, parameter:: Zeta2=1.6449340668482264365D0

      Real*8 :: DL, y1, DL1, yy, yp
      Real*8 :: DL2,DL3,DLm2,DLm3,DLp
      Real*8 :: snp12, xLi

cmarco11 Approximate expression
      yy=y*y
      DL = dlog(y)
      y1 = 1d0 - y
      DL1 = dlog(y1)

      clnp2a_CT = 128./9.*y*DL1*DL1 -46.5*y*DL1-84.094*DL*DL1
     >       -37.338 +89.53*y +33.82*yy +y*DL*(32.9+18.41*DL)
     >       -128./9.*DL +16./27.*Nf*(6.*y*DL1 -12.*y*DL -25.*y +6.)

cmarco11 Exact expression
c$$$      yp = 1d0 + y
c$$$      DLp = dlog(yp)
c$$$      DL2 = xLi(2,y)
c$$$      DL3 = xLi(3,y)
c$$$      DLm2= xLi(2,-y)
c$$$      DLm3= xLi(3,-y)
c$$$
c$$$      clnp2a=(4d0*Cf*(Ca - 2d0*Cf)*yy*(
c$$$     >  4d0*DL*(6d0 - 3d0*y + 47d0*yy - 9d0*yy*y)/15d0/yy
c$$$     >  - 4d0*DLm2*(DL-2d0*DLp) - 8d0*Zeta3 - 2d0*DL**2*dlog(1d0-yy)
c$$$     >  + 4d0*DL*DLp**2 - 4d0*DL*DL2 + 2d0/5d0*(5d0-3d0*yy)*DL**2
c$$$     >  - 4d0*(2d0+10d0*yy+5d0*yy*y-3d0*yy*yy*y)/5d0/yy/y*(DLm2+DL*DLp)
c$$$     >  + 4d0*Zeta2*(dlog(1d0-yy) -(5d0-3d0*yy)/5d0)+8d0*snp12(-y)
c$$$     >  + 4d0*DL3-4d0*DLm3-23d0/3d0*DL1-(144d0+294d0*y -1729*yy+216*yy*y
c$$$     >  )/90d0/yy)
c$$$     >  + 8d0*Cf*Cf*yy*(DL2 + (dlog(y/y1))**2 - 3d0*Zeta2 
c$$$     >  - (3d0-22*y)/3d0/y*DL + (6d0-25d0*y)/6d0/y*DL1
c$$$     >  - (78d0-355d0*y)/36d0/y) 
c$$$     >  - 8d0/3d0*Cf*Nf*yy*(dlog(yy/y1)-(6d0-25d0*y)/6d0/y) )/y

      Return

      Entry cls2aL(y, Nf)
c NLO pure-singlet coefficient functions for F_L

      DL = Log(y)
      y1 = 1d0 - y
      DL1 = Log(y1)

      cls2aL= Nf*((15.94 -5.212*y)*y1*y1*DL1 +(.421 +1.52*y)*DL*DL
     >      +28.09*y1*DL -(2.37/y -19.27)* y1**3)
      Return

      Entry clg2aL(y, Nf)
c NLO gluon coefficient functions for F_L

      DL = Log(y)
      y1 = 1d0 - y
      DL1 = Log(y1)

      clg2aL= Nf*((94.74 -49.2*y)*y1*DL1*DL1 +864.8*y1*DL1
     >       +1161*y*DL1*DL +60.06*y*DL*DL +39.66*y1*DL
     >       -5.333*(1./y-1.))
      Return

      End


      Function clmnp2a(y,Nf) !factorization scale dependent term
      Implicit none
      Real*8 :: y
      Integer :: Nf
      Real*8 :: clmnp2a, clms2a, clmg2a

      Real*8, Parameter:: Cf=4D0/3d0, Ca=3d0, Tf=.5d0
      Real*8 :: DL, y1, DL1

      DL = Log(y)
      y1 = 1d0 - y
      DL1 = Log(y1)

      clmnp2a=Cf*Cf*(8d0*y*(2d0*DL1-DL)+4d0*(2d0+y))
     >     +Ca*Cf*(-44d0/3d0*y)+Nf*Cf*8d0/3d0*y
      return

      Entry clms2a(y,Nf)
      DL = Log(y)
      clms2a= Nf*Cf*Tf*(-32d0*y*DL -32d0/3d0*(3d0-2d0*y*y-1d0/y))
      return


      Entry clmg2a(y,Nf)
      DL = Log(y)
      y1 = 1d0 - y
      DL1 = Log(y1)
      clms2a= Nf*(Cf*Tf*(32d0*y*DL +16d0*(1d0-2d0*y*y+y))
     >     +Ca*Tf*(64d0*y*y1*DL1-128d0*y*DL
     >             -32d0/3d0*(3d0-17d0*y*y+15d0*y-1d0/y)))
      return
      end

cmarco=================================================================
c Separate implementation of ZM c2PS2 from A. Vogt hep-ph/0504242 
c This has been extrapolated from Liang's ZmVVA.f file
cmarco=================================================================


      Real*8 Function F2ZMPS2(X, Q, amu, Ihad)
C     The zero-mass NNLO pure singlet quark F2 for 1 quark flavor 
c     in the coefficient function (and singlet PDF for the global Nf flavors)
      Use mDisPac
      Implicit None
	double precision smpsnf, err
	double precision a5, tmp, sum
	double precision tmp1, tmp2, ran1, hwrgen
	integer inx

      Real*8, intent(in) :: X, Q, aMu
      Integer, intent(in):: Ihad

      Real*8 :: ADZINT_CT, alpi, F2CqQPS !functions
      external F2CqQPS
      Real*8, Parameter :: D1=1D0 
      Integer, Parameter :: IACTA=0, IACTB=1
      Real*8 :: ERREST, tF2ZMPS2
      Integer :: IER
      Real*8 :: XX, QQ, bmu, aMu0, als4pi
      common /F2LZM/ XX, QQ, bmu
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      QQ=Q
      XX=X
      If(XX.ge.1d0) then
        tF2ZMPS2=0d0
        return
      Endif
      Ihadn=Ihad

      bmu=amu
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

	a5 = 1.d0 - 0.00000000001d0*(1.d0-XX)

	sum = 
     >   SMPSNF(F2CqQPS,XX + 0.00000000001d0*(a5-XX), 
     >	                XX + 0.000000001d0*(a5-XX),4, ERR, IER)
     > + SMPSNF(F2CqQPS,XX + 0.000000001d0*(a5-XX),
     >                  XX + 0.00000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.00000001d0*(a5-XX),
     >                  XX + 0.0000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.0000001d0*(a5-XX),
     >                  XX + 0.000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.000001d0*(a5-XX),
     >                  XX + 0.00001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.00001d0*(a5-XX),
     >                  XX + 0.0001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.0001d0*(a5-XX),
     >                  XX + 0.001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.001d0*(a5-XX),
     >                  XX + 0.01d0*(a5-XX),8, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.01d0*(a5-XX),
     >                  XX + 0.1d0*(a5-XX),16, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.1d0*(a5-XX),
     >                  XX + 0.9d0*(a5-XX),32, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.9d0*(a5-XX),
     >                  XX + 0.99d0*(a5-XX),32, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.99d0*(a5-XX),
     >                  XX + 0.999d0*(a5-XX),8, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.999d0*(a5-XX),
     >                  XX + 0.9999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.9999d0*(a5-XX),
     >                  XX + 0.99999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.99999d0*(a5-XX),
     >                  XX + 0.999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.999999d0*(a5-XX),
     >                  XX + 0.9999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.9999999d0*(a5-XX),
     >                  XX + 0.99999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.99999999d0*(a5-XX),
     >                  XX + 0.999999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CqQPS,XX + 0.999999999d0*(a5-XX),a5,4, ERR, IER)


      F2ZMPS2= (als4pi)**2 * sum 

      end

      Function F2CqQPS(y)
      Use mDisPac
      Implicit None

      Real*8 :: F2CqQPS
      Real*8, intent(in) :: y
      Real*8 :: XX, QQ, bmu
      common /F2LZM/ XX, QQ, bmu
      Real*8 :: aLm, PDF_CT, Qs, Qs0
      Integer :: Nff=1
      Integer Ihadn, Ir
      Common /F2LHQhad/ Ihadn

      Real*8 :: CL, Dk
      Real*8 :: QGtype, c2np1a0, c2np2a0, c2np3a0  ! function
      Real*8 :: c2np1Dk, c2np2Dk, c2np3Dk
      Real*8 :: c2mnp1a0, c2mnp1Dk,c2s2aL
      Real*8 :: c2mnp2a2, c2mnp2a1, c2mnp2Dk2, c2mnp2Dk1
      Real*8 :: c2ms2a2, c2ms2a1

      aLm=dlog(QQ**2/bmu**2)

      Qs = PDF_CT(ihadn, 20, XX/y, bmu, Ir)

      Qs0 = PDF_CT(ihadn, 20, XX, bmu, Ir)


      CL= c2s2aL(y,Nff)         !only pure-singlet

      If(abs(QQ-bmu).gt.1e-6) then
! See Eq B.4 in Nucl Phys B383 by Zijlstra and Van Neerven, 
! there is no aLm dependence in the pure singlet case.  
! See also Nucl. Phys. B588 345 (2000) A. Vogt
         CL=CL + 0d0                   
      Endif

      F2CqQPS = CL * Qs / y 
      Return
      End





cmarco=================================================================
c Separate implementation of ZM c2g2 from A. Vogt hep-ph/0504242 
c This has been extrapolated from Liang's ZmVVA.f file
cmarco=================================================================


      Real*8 Function F2ZMG2 (X,Q,aMu,iHad)
C     The zero-mass NNLO gluon contribution for 1 quark flavor 
c     in the coefficient function (and singlet PDF for the global Nf flavors)
      Use mDisPac
      Implicit None
	double precision smpsnf, err
	double precision a5, tmp, sum
	double precision tmp1, tmp2, ran1, hwrgen
	integer inx

      Real*8, intent(in) :: X, Q, aMu
      Integer, intent(in):: Ihad

      Real*8 :: ADZINT_CT, alpi, F2CgG !functions
      external F2CgG
      Real*8, Parameter :: D1=1D0 
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST, tF2ZMG2
      Integer :: IER
      Real*8 :: XX, QQ, bmu, aMu0, als4pi
      common /F2LZM/ XX, QQ, bmu
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      QQ=Q
      XX=X
      If(XX.ge.1d0) then
        tF2ZMG2=0d0
        return
      Endif
      Ihadn=Ihad

      bmu=amu
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

	a5 = 1.d0 - 0.00000000001d0*(1.d0-XX)

	sum = 
     >   SMPSNF(F2CgG,XX + 0.00000000001d0*(a5-XX), 
     >	              XX + 0.0000000001d0*(a5-XX),4, ERR, IER)
     > + SMPSNF(F2CgG,XX + 0.0000000001d0*(a5-XX),
     >                XX + 0.000000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.000000001d0*(a5-XX),
     >                XX + 0.00000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.00000001d0*(a5-XX),
     >                XX + 0.0000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.0000001d0*(a5-XX),
     >                XX + 0.000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.000001d0*(a5-XX),
     >                XX + 0.00001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.00001d0*(a5-XX),
     >                XX + 0.0001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.0001d0*(a5-XX),
     >                XX + 0.001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.001d0*(a5-XX),
     >                XX + 0.01d0*(a5-XX),8, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.01d0*(a5-XX),
     >                XX + 0.1d0*(a5-XX),16, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.1d0*(a5-XX),
     >                XX + 0.9d0*(a5-XX),32, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.9d0*(a5-XX),
     >                XX + 0.99d0*(a5-XX),16, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.99d0*(a5-XX),
     >                XX + 0.999d0*(a5-XX),8, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.999d0*(a5-XX),
     >                XX + 0.9999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.9999d0*(a5-XX),
     >                XX + 0.99999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.99999d0*(a5-XX),
     >                XX + 0.999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.999999d0*(a5-XX),
     >                XX + 0.9999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.9999999d0*(a5-XX),
     >                XX + 0.99999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.99999999d0*(a5-XX),
     >                XX + 0.999999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(F2CgG,XX + 0.999999999d0*(a5-XX),a5,4, ERR, IER)



      F2ZMG2= (als4pi)**2 * sum

      end


      Function F2CgG(y)
      Use mDisPac
      Implicit None

      Real*8, intent(in) :: y
      Real*8 :: F2CgG
      Real*8 :: PDFH !functions
      Real*8 :: XX, QQ, bmu
      common /F2LZM/ XX, QQ, bmu
      Integer Ihadn
      Common /F2LHQhad/ Ihadn

      Integer :: Nff=1

      Real*8 :: CL, GG, alm
      Real*8 :: QGtype, c2g2aL  ! functions
      Real*8 :: c2mg1a
      Real*8 :: c2mg2a2, c2mg2a1

      CL= c2g2aL(y,NFf)
      If(abs(QQ-bmu).gt.1e-6) then
         aLm=dlog(QQ**2/bmu**2)
         CL=CL + aLm*(aLm*c2mg2a2(y,Nff)+c2mg2a1(y,Nff))
      endif
        
      GG = PDFH(Ihadn, 0, XX/y, bmu)
      
      F2CgG = CL * GG / y
      Return
      End


cmarco=================================================================
c Separate implementation of ZM cLPS2 from A. Vogt hep-ph/0504242 
c This has been extrapolated from Liang's ZmVVA.f file
cmarco=================================================================


      Real*8 Function FLZMPS2(X, Q, amu, Ihad)
C     The zero-mass NNLO pure singlet quark F2 for 1 quark flavor 
c     in the coefficient function (and singlet PDF for the global Nf flavors)
      Use mDisPac
      Implicit None
	double precision smpsnf, err
	double precision a5, tmp, sum
	double precision tmp1, tmp2, ran1, hwrgen
	integer inx

      Real*8, intent(in) :: X, Q, aMu
      Integer, intent(in):: Ihad

      Real*8 :: ADZINT_CT, alpi, FLCqQPS !functions
      external FLCqQPS
      Real*8, Parameter :: D1=1D0 
      Integer, Parameter :: IACTA=0, IACTB=1
      Real*8 :: ERREST, tFLZMPS2
      Integer :: IER
      Real*8 :: XX, QQ, bmu, aMu0, als4pi
      common /F2LZM/ XX, QQ, bmu
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      QQ=Q
      XX=X
      If(XX.ge.1d0) then
        tFLZMPS2=0d0
        return
      Endif
      Ihadn=Ihad

      bmu=amu
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

	a5 = 1.d0 - 0.00000000001d0*(1.d0-XX)

	sum = 
     >   SMPSNF(FLCqQPS,XX + 0.00000000001d0*(a5-XX), 
     >	                XX + 0.0000000001d0*(a5-XX),4, ERR, IER)
     > + SMPSNF(FLCqQPS,XX + 0.0000000001d0*(a5-XX),
     >                  XX + 0.000000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.000000001d0*(a5-XX),
     >                  XX + 0.00000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.00000001d0*(a5-XX),
     >                  XX + 0.0000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.0000001d0*(a5-XX),
     >                  XX + 0.000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.000001d0*(a5-XX),
     >                  XX + 0.00001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.00001d0*(a5-XX),
     >                  XX + 0.0001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.0001d0*(a5-XX),
     >                  XX + 0.001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.001d0*(a5-XX),
     >                  XX + 0.01d0*(a5-XX),8, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.01d0*(a5-XX),
     >                  XX + 0.1d0*(a5-XX),16, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.1d0*(a5-XX),
     >                  XX + 0.9d0*(a5-XX),32, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.9d0*(a5-XX),
     >                  XX + 0.99d0*(a5-XX),16, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.99d0*(a5-XX),
     >                  XX + 0.999d0*(a5-XX),8, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.999d0*(a5-XX),
     >                  XX + 0.9999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.9999d0*(a5-XX),
     >                  XX + 0.99999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.99999d0*(a5-XX),
     >                  XX + 0.999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.999999d0*(a5-XX),
     >                  XX + 0.9999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.9999999d0*(a5-XX),
     >                  XX + 0.99999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.99999999d0*(a5-XX),
     >                  XX + 0.999999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCqQPS,XX + 0.999999999d0*(a5-XX),a5,4, ERR, IER)


c includes marco's Divide by 2, to match with the normalization of F_{long}
      FLZMPS2= (als4pi)**2 * sum /2.d0

      end

      Function FLCqQPS(y)
      Use mDisPac
      Implicit None
      Real*8 :: FLCqQPS
      Real*8, intent(in) :: y
      Real*8 :: XX, QQ, bmu, PDF_CT
      common /F2LZM/ XX, QQ, bmu
      Integer Ihadn, ir
      Common /F2LHQhad/ Ihadn     
      Integer :: Nff=1

      Real*8 :: CL, aLm, Qs
      Real*8 :: QGtype, c2np1a0, c2np2a0, c2np3a0  ! function
      Real*8 :: c2np1Dk, c2np2Dk, c2np3Dk
      Real*8 :: c2mnp1a0, c2mnp1Dk,c2s2aL
      Real*8 :: c2mnp2a2, c2mnp2a1, c2mnp2Dk2, c2mnp2Dk1
      Real*8 :: c2ms2a2, c2ms2a1, cls2aL, clms2a

      Qs = PDF_CT(ihadn, 20, XX/y, bmu, Ir)

      CL= cls2aL(y,Nff)    ! pure-singlet from Eq.5 in PLB606 123 2005 A. Vogt
      If(abs(QQ-bmu).gt.1e-6) then
         aLm=dlog(QQ**2/bmu**2)
         CL= CL + clms2a(y,Nff)*aLm !Eq B.3 NPB383 525 1992, Zijlstra, van Neeven
      endif

      FLCqQPS = CL * Qs / y
      Return
      End




cmarco=================================================================
c Separate implementation of ZM cLg2 from A. Vogt hep-ph/0504242 
c This has been extrapolated from Liang's ZmVVA.f file
cmarco=================================================================


      Real*8 Function FLZMG2 (X,Q,aMu,iHad)
C     The zero-mass NNLO gluon FL for 1 quark flavor 
c     in the coefficient function 
      Use mDisPac
      Implicit None
	double precision smpsnf, err
	double precision a5, tmp, sum
	double precision tmp1, tmp2, ran1, hwrgen
	integer inx

      Real*8, intent(in) :: X, Q, aMu
      Integer, intent(in):: Ihad

      Real*8 :: ADZINT_CT, alpi, FLCgG !functions
      external FLCgG
      Real*8, Parameter :: D1=1D0 
      Integer, Parameter :: IACTA=2, IACTB=2
      Real*8 :: ERREST, tFLZMG2
      Integer :: IER
      Real*8 :: XX, QQ, bmu, aMu0, als4pi
      common /F2LZM/ XX, QQ, bmu
      Integer Ihadn
      Common /F2LHQhad/ Ihadn
      Data aMu0 /-1d0/
      save aMu0

      QQ=Q
      XX=X
      If(XX.ge.1d0) then
        tFLZMG2=0d0
        return
      Endif
      Ihadn=Ihad

      bmu=amu
      If(aMu/=aMu0) als4pi=alpi(aMu)/4D0
      aMu0=aMu

	a5 = 1.d0 - 0.00000000001d0*(1.d0-XX)

	sum = 
     >   SMPSNF(FLCgG,XX + 0.00000000001d0*(a5-XX), 
     >	              XX + 0.000000001d0*(a5-XX),4, ERR, IER)
     > + SMPSNF(FLCgG,XX + 0.000000001d0*(a5-XX),
     >                XX + 0.00000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.00000001d0*(a5-XX),
     >                XX + 0.0000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.0000001d0*(a5-XX),
     >                XX + 0.000001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.000001d0*(a5-XX),
     >                XX + 0.00001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.00001d0*(a5-XX),
     >                XX + 0.0001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.0001d0*(a5-XX),
     >                XX + 0.001d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.001d0*(a5-XX),
     >                XX + 0.01d0*(a5-XX),8, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.01d0*(a5-XX),
     >                XX + 0.1d0*(a5-XX),16, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.1d0*(a5-XX),
     >                XX + 0.9d0*(a5-XX),32, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.9d0*(a5-XX),
     >                XX + 0.99d0*(a5-XX),16, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.99d0*(a5-XX),
     >                XX + 0.999d0*(a5-XX),8, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.999d0*(a5-XX),
     >                XX + 0.9999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.9999d0*(a5-XX),
     >                XX + 0.99999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.99999d0*(a5-XX),
     >                XX + 0.999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.999999d0*(a5-XX),
     >                XX + 0.9999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.9999999d0*(a5-XX),
     >                XX + 0.99999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.99999999d0*(a5-XX),
     >                XX + 0.999999999d0*(a5-XX),4, ERR, IER) 
     > + SMPSNF(FLCgG,XX + 0.999999999d0*(a5-XX),a5,4, ERR, IER)


c include marco's Divide by 2, to match with the normalization of F_{long}
      FLZMG2 = (als4pi)**2 * sum  / 2.d0

      end

      Function FLCgG(y)
      Use mDisPac
      Implicit None

      Real*8, intent(in) :: y
      Real*8 :: FLCgG
      Real*8 :: PDFH !functions
      Real*8 :: XX, QQ, bmu
      common /F2LZM/ XX, QQ, bmu
      Integer Ihadn
      Common /F2LHQhad/ Ihadn

      Real*8 :: CL, GG,alM 
      Real*8 :: QGtype, c2g2aL, c2g3aL  ! functions
      Real*8 :: c2mg1a
      Real*8 :: c2mg2a2, c2mg2a1, clmg2a, clg2aL 
      Integer :: Nff=1

      CL= clg2aL(y,Nff)
      If(abs(QQ-bmu).gt.1e-6) then
         aLm=dlog(QQ**2/bmu**2)
         CL = CL + clmg2a(y,Nff)*aLm
      endif
         
      GG = PDFH(Ihadn, 0, XX/y, bmu)
      
      FLCgG = CL * GG / y
      Return
      End
