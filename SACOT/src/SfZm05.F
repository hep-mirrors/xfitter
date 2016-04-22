      SUBROUTINE SfZm0Q (Ihad, Iprtn, X, Q, aMu, SfnG)
C-----------------------------------------------------------------------------
C      Computes ZM LO Hadron Structure Functions SfnG 
C          -- individual contribution from Iprtn in the initial state
C
C        (x, Q) are input kinematic variables;
C        aMu is the factorization scale, which does not have to be equal to Q.
C-----------------------------------------------------------------------------
       Implicit None
C                                              Dummy Arguments       
       Integer, Intent(In) :: Ihad, Iprtn 
c htj
       integer :: i,nx,nfl
       real*8 :: Zmin, Xmin,zz,xx,qq,dz,wght,pdfh
       real*8::glu,uqk,ubr,dqk,dbr,sqk,sbr,cqk,bqk,uvl,dvl,svl,xsv,tqk

       Real*8,  Intent(In) :: X, Q, aMu
       Real*8,  Intent(Out):: SfnG(5)
C                                              Local variables
       Integer:: iRet
       Real*8 :: OmgLo(5), PdfTmp
C                                              External function
       Real*8 :: PDF_CT
C                          -------------
C            At LO, the non-zero partonic irreducible amplitudes 
C            (sans EW couplings) are transverse and = 1/2 (so that F_1=F_2=F_3=1.)
       OmgLo = (/5d-1, 0d0, 0d0, 0d0, 5d-1/)          
C                                              PDF factor
       PdfTmp = PDF_CT(Ihad, Iprtn, X, aMu, iRet)
C                                              Result       
       SfnG = OmgLo * PdfTmp          

c       print*,PDFTmp
c       stop
c htj test  
c      open(1001,file="PDFplt.xfitter_Q100.0.txt",status="unknown")
c      Write (1001, '(A1,A15,20A16)') '#','Z',' x '
c     > ,'xglu','xuqk', 'xubr','xdqk','xdbr','xsqk','xsbr','xcqk'
c     > ,'xbqk','uvl','dvl','svl','xsv','tqk'

c      nx=100
c      nfl=5
c      qq=100.0d0
c      xmin=1d-6
c      Zmin = Xmin**(1d0/3d0)
c      Dz = (1d0-Zmin) / Nx
c      Do I = 1, Nx-1
c      zz = Zmin + I * Dz
c      xx = zz **3d0
c      wght = (3d0 * xx / zz)

c      glu = pdfh(1, 0, xx, qq) * wght * xx
c      uqk = pdfh(1, 1, xx, qq) * wght * xx
c      ubr = pdfh(1,-1, xx, qq) * wght * xx
c      dqk = pdfh(1, 2, xx, qq) * wght * xx
c      dbr = pdfh(1,-2, xx, qq) * wght * xx
c      sqk = pdfh(1, 3, xx, qq) * wght * xx
c      sbr = pdfh(1,-3, xx, qq) * wght * xx
c      if (nfl.ge.4) cqk = pdfh(1, 4, xx, Qq) * wght * xx
c      if (nfl.ge.5) bqk = pdfh(1, 5, xx, Qq) * wght * xx
c      if (nfl.ge.6) tqk = pdfh(1, 6, xx, Qq) * wght * xx
c      uvl = (pdfh(1, 1, xx, Qq) - pdfh(1,-1, xx, Qq)) * wght
c      dvl = (pdfh(1, 2, xx, Qq) - pdfh(1,-2, xx, Qq)) * wght
c      svl = (pdfh(1, 3, xx, Qq) - pdfh(1,-3, xx, Qq)) * wght
c      xsv = svl * xx

c      Write (1001, '(20ES16.8)')
c     >zz, xx, glu, uqk, ubr, dqk, dbr, sqk, sbr, cqk, bqk,
c     >uvl,dvl,svl,xsv,tqk

c      end do
c      close(1001)
c      print*,"pdfplt xfitter done"
c      stop


       END Subroutine SfZm0Q                ! *******************
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      Subroutine SfZm1G (Ihad, XX, Q, vMu, SFnG)

    !  ****** contribution to SF due to initial state gluon-boson fusion diagrams,
    !  without the EW coupling--hence this is for one particular quark-antiquark pair
    !              (i.e. the result is flavor-independent);
    !  Also not included is the overall \alpha_s / 2pi factor.
     
    !          Formulas from Furmanski & Petronzio
    !           ---------------------------------------------- 
      Use mDisPac
      Implicit None
	real*8 :: tmp, tmp1, smpsnf, a5
C >>                                       ! dummy variable
      integer, intent(in) :: iHad
      Real*8, intent (in) :: XX, Q, vMu
      Real*8, intent(out) :: SfnG(5)
C                                          ! Global variables:
      Integer :: iHadn
      Real*8  :: x, aMu
      Common / CSfZm1G / x, aMu, iHadn
C --                                     ! Local variables    
      Integer IR
      Real*8 :: F2G1, ER, F2MF1G1, Qmu, Smll=1d-6, OmgLo(5)
C --                                     ! Functions
      Real*8 :: aC2G1, aC2MC1G1, aPgqFx, ADZINT_CT
      
      EXTERNAL aC2G1, aC2MC1G1, aPgqFx
C                            !    ------- fill integration common
      iHadn= iHad
      x    = XX
      aMu  = vMu
c htj
c      print*,ihadn,x,amu
c      stop
C                                                            ! F2 integral
c xx      F2G1 = ADZINT_CT(aC2G1,    x, 1d0, aErr, rErr, ER, IR, 2, 1)

	a5 = 1.d0 - 0.00000000001d0*(1.d0-x)
      F2G1 = 
     >   SMPSNF(aC2G1,x + 0.00000000001d0*(a5-x), 
     >	              x + 0.0000001d0*(a5-x),4, ER, IR)
     > + SMPSNF(aC2G1,x + 0.0000001d0*(a5-x),
     >                x + 0.000001d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2G1,x + 0.000001d0*(a5-x),
     >                x + 0.00001d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2G1,x + 0.00001d0*(a5-x),
     >                x + 0.0001d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2G1,x + 0.0001d0*(a5-x),
     >                x + 0.001d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2G1,x + 0.001d0*(a5-x),
     >                x + 0.01d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2G1,x + 0.01d0*(a5-x),
     >                x + 0.1d0*(a5-x),8, ER, IR) 
     > + SMPSNF(aC2G1,x + 0.1d0*(a5-x),
     >                x + 0.9d0*(a5-x),16, ER, IR) 
     > + SMPSNF(aC2G1,x + 0.9d0*(a5-x),
     >                x + 0.99d0*(a5-x),8, ER, IR) 
     > + SMPSNF(aC2G1,x + 0.99d0*(a5-x),
     >                x + 0.999d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2G1,x + 0.999d0*(a5-x),
     >                x + 0.9999d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2G1,x + 0.9999d0*(a5-x),
     >                x + 0.99999d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2G1,x + 0.99999d0*(a5-x),
     >                x + 0.999999d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2G1,x + 0.999999d0*(a5-x),a5,4, ER, IR)

C                                                            ! F2 - F1 integral
c xx      F2MF1G1 = ADZINT_CT(aC2MC1G1, x, 1d0, aErr, rErr, ER, IR, 2, 1)

	a5 = 1.d0 - 0.00000000001d0*(1.d0-x)
      F2MF1G1 = 
     >   SMPSNF(aC2MC1G1,x + 0.00000000001d0*(a5-x), 
     >	                 x + 0.0000001d0*(a5-x),4, ER, IR)
     > + SMPSNF(aC2MC1G1,x + 0.0000001d0*(a5-x),
     >                   x + 0.000001d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2MC1G1,x + 0.000001d0*(a5-x),
     >                   x + 0.00001d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2MC1G1,x + 0.00001d0*(a5-x),
     >                   x + 0.0001d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2MC1G1,x + 0.0001d0*(a5-x),
     >                   x + 0.001d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2MC1G1,x + 0.001d0*(a5-x),
     >                   x + 0.01d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2MC1G1,x + 0.01d0*(a5-x),
     >                   x + 0.1d0*(a5-x),8, ER, IR) 
     > + SMPSNF(aC2MC1G1,x + 0.1d0*(a5-x),
     >                   x + 0.9d0*(a5-x),16, ER, IR) 
     > + SMPSNF(aC2MC1G1,x + 0.9d0*(a5-x),
     >                   x + 0.99d0*(a5-x),8, ER, IR) 
     > + SMPSNF(aC2MC1G1,x + 0.99d0*(a5-x),
     >                   x + 0.999d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2MC1G1,x + 0.999d0*(a5-x),
     >                   x + 0.9999d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2MC1G1,x + 0.9999d0*(a5-x),
     >                   x + 0.99999d0*(a5-x),4, ER, IR) 
     > + SMPSNF(aC2MC1G1,x + 0.99999d0*(a5-x),a5,4, ER, IR)

C                                         ! F_transv = (F_R + F_L)/2= F1 /2 
      SfnG(1) = (F2G1 - F2MF1G1) /2d0
C                                         ! F_long = (F2 - F1) / 2
      SfnG(3) = F2MF1G1 /2d0
C                                      
      SfnG(2) = 0d0        ! Zero-mass case has no chirality-mixing terms
      SfnG(4) = 0d0        ! Zero-mass case has no chirality-mixing terms
      SfnG(5) = 0d0        ! Gluon term does not have parity violating term
      
                                     ! Log(Q/mu) term:
      If (Abs(Q/aMu - 1d0) > Smll) Then
c **      Qmu = ADZINT_CT(aPgqFx, X, 1d0, aErr, rErr, ER, IR, 1, 1)

	a5 = 1.d0 - 0.00000000001d0*(1.d0-x)

	Qmu = 
     >   SMPSNF(aPgqFx,x + 0.00000000001d0*(a5-x), 
     >	               x + 0.000001d0*(a5-x),4, Er, IR)
     > + SMPSNF(aPgqFx,x + 0.000001d0*(a5-x),
     >                 x + 0.00001d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aPgqFx,x + 0.00001d0*(a5-x),
     >                 x + 0.0001d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aPgqFx,x + 0.0001d0*(a5-x),
     >                 x + 0.001d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aPgqFx,x + 0.001d0*(a5-x),
     >                 x + 0.01d0*(a5-x),8, Er, IR) 
     > + SMPSNF(aPgqFx,x + 0.01d0*(a5-x),
     >                 x + 0.1d0*(a5-x),16, Er, IR) 
     > + SMPSNF(aPgqFx,x + 0.1d0*(a5-x),
     >                 x + 0.9d0*(a5-x),32, Er, IR) 
     > + SMPSNF(aPgqFx,x + 0.9d0*(a5-x),
     >                 x + 0.99d0*(a5-x),16, Er, IR) 
     > + SMPSNF(aPgqFx,x + 0.99d0*(a5-x),
     >                 x + 0.999d0*(a5-x),8, Er, IR) 
     > + SMPSNF(aPgqFx,x + 0.999d0*(a5-x),
     >                 x + 0.9999d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aPgqFx,x + 0.9999d0*(a5-x),
     >                 x + 0.99999d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aPgqFx,x + 0.99999d0*(a5-x),
     >                 x + 0.999999d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aPgqFx,x + 0.999999d0*(a5-x),a5,4, Er, IR)
      
      OmgLo = (/5d-1, 0d0, 0d0, 0d0, 5d-1/)
      SfnG = SfnG + 2d0 * OmgLo * Log(Q/aMu) * Qmu 
      End If

      SfnG = SfnG * TR     ! QCD color factor
      
      End Subroutine SfZm1G                ! ***************************
C==================                                        
      Function   aC2G1(z)   ! Integrand for the F2_G integral
      
      Implicit None
C                                          ! Global variables:
      Integer :: iHadn
      Real*8  :: x, aMu
      Common / CSfZm1G / x, aMu, iHadn
C                                          ! Local variables      
      Integer :: IR
      Real*8  :: z, Pdfnxz, aC2G1, C2G1A
C                                          ! External function
      Real*8  :: PDF_CT
      
C                                               executables   
      Pdfnxz = PDF_CT(iHadn, 0, x/z, aMu, IR) 
           
      C2G1A=((z**2+(1d0-z)**2)*log((1d0-z)/z)-1d0+8d0*z*(1d0-z))
      
      aC2G1 = C2G1A * Pdfnxz / z      
      
      End Function aC2G1                    ! ***************************
C==============                                        
      Function aC2MC1G1(z)   ! Integrand for the gluon F2 - F1 integral
      
      Implicit None
C                                          ! Global variables:
      Integer :: iHadn
      Real*8  :: x, aMu
      Common / CSfZm1G / x, aMu, iHadn
C                                          ! Local variables      
      Integer :: IR
      Real*8  :: z, Pdfnxz, aC2MC1G1, C2MC1G1A
C                                          ! External function
      Real*8  :: PDF_CT
      
      Pdfnxz = PDF_CT(iHadn, 0, x/z, aMu, IR)
      
      C2MC1G1A= 4d0 * z * (1d0-z)

      aC2MC1G1 = C2MC1G1A * Pdfnxz /z
     
      End Function aC2MC1G1                 ! ***************************
C======

      Function aPgqFx(z)   ! Integrand for the gluon Log(Q/mu) integral
      
      Implicit None
C                                          ! Global variables:
      Integer :: iHadn
      Real*8  :: x, aMu
      Common / CSfZm1G / x, aMu, iHadn
C                                          ! Local variables      
      Integer :: IR
      Real*8  :: z, Pdfnxz, SpltFn, aPgqFx
C                                          ! External function
      Real*8  :: PDF_CT
      
C                                           ! The gluon PDF
      Pdfnxz = PDF_CT(iHadn, 0, x/z, aMu, IR)
C                                           ! The gluon-->quark splitting fn
      SpltFn = ( z**2 + (1d0-z)**2 ) /2d0
C                                           ! The integrand      
      aPgqFx = SpltFn * Pdfnxz / z
      
      End function aPgqFx

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      Subroutine SfZm1Q  (iHad, iPrtn, XX, Q, vMu, SFnG)    
    
    !  ****** contribution to the general SfnG due to initial state quark parton iPrtn
    !  without the EW coupling--the result is flavor-independent;
    !  Also not included is the overall \alpha_s / 2pi factor.
    ! 
    !          Formulas from Furmanski & Petronzio
    !           ---------------------------------------------- 
      Use mDisPac
      Implicit None
	real*8 :: tmp, tmp1, smpsnf, a5
	integer :: Num
C >>                                       ! dummy variable
      integer, intent(in) :: iHad, iPrtn
      Real*8, intent (in) :: XX, Q, vMu
      Real*8, intent(out) :: SfnG(5)
C                                          ! Global variables:
      Integer :: iHadn, iPartn, Ir
      Real*8  :: x, aMu, Pdfnx
      Common / CSfZm1Q / x, aMu, Pdfnx, iHadn, iPartn
C --                                     ! Local variables    
      Real*8 :: F2Qa, F2Qb, F2M1Q, F1m3Q, ER, F1,f2,f3, OmgLo(5)
     >        , Qmu1, Qmu2, Smll = 1d-6, AdzIntRes
C --                                     ! External functions
      Real*8 :: aC2Q1, aF2Q1, G2M1Q1, G1M3Q1, aPqq, aPqqFx, ADZINT_CT
     >        , PDF_CT ,xLi

      EXTERNAL aC2Q1, aF2Q1, G2M1Q1, G1M3Q1, aPqq, aPqqFx, xLi
C                                !    ------- fill integration common
      iHadn= iHad
      iPartn = iPrtn
      x    = XX
      aMu  = vMu

      Pdfnx = PDF_CT(iHadn, iPartn, x, aMu, IR)
                                ! For F2Q, this is the \delta(1-z) term               
                        ! here we use the analytic result for the  
                        ! integral of aC2Q1 to speed up the computation  
cmarco11 Old implementation
c      F2Qa = ADZINT_CT(aC2Q1, 0d0, x, aErr, rErr, ER, IR, 1, 0) 
c     >     * Pdfnx
      AdzIntRes = (-2d0*Log(1d0 - x)**2 - Log(1d0 - x)*(-6d0 + 2d0*x + x
     >  **2 -4d0*Log(x)) + x*(7d0 + 2d0*x + (2d0 + x)*Log(x))+ 4d0*xLi(2
     >  ,x))/2d0
      F2Qa = AdzIntRes * Pdfnx
                                     ! The F2 regular convolution term
	a5 = 1.d0 - 0.00000000001d0*(1.d0-x)

	F2Qb = 
     >   SMPSNF(aF2Q1,x + 0.00000000001d0*(a5-x), 
     >	              x + 0.00000001d0*(a5-x),4, Er, IR)
     > + SMPSNF(aF2Q1,x + 0.00000001d0*(a5-x),
     >                x + 0.0000001d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.0000001d0*(a5-x),
     >                x + 0.000001d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.000001d0*(a5-x),
     >                x + 0.00001d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.00001d0*(a5-x),
     >                x + 0.0001d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.0001d0*(a5-x),
     >                x + 0.001d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.001d0*(a5-x),
     >                x + 0.01d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.01d0*(a5-x),
     >                x + 0.1d0*(a5-x),8, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.1d0*(a5-x),
     >                x + 0.9d0*(a5-x),16, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.9d0*(a5-x),
     >                x + 0.99d0*(a5-x),8, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.99d0*(a5-x),
     >                x + 0.999d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.999d0*(a5-x),
     >                x + 0.9999d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.9999d0*(a5-x),
     >                x + 0.99999d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.99999d0*(a5-x),
     >                x + 0.999999d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.999999d0*(a5-x),
     >                x + 0.9999999d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aF2Q1,x + 0.9999999d0*(a5-x),a5,4, Er, IR)



                                     ! The F2-F1 (2*Flong) term

c xx      F2m1Q = ADZINT_CT(G2M1Q1, X, 1d0, aErr, rErr, ER, IR, 0, 0)

      F2m1Q = SMPSNF(G2M1Q1,x,1.d0,32, Er, Ir)


c xx      F1m3Q = ADZINT_CT(G1M3Q1, X, 1d0, aErr, rErr, ER, IR, 0, 0)

                                     ! The F1-F3 (2*Fright) term
	a5 = 1.d0 - 0.00000000001d0*(1.d0-x)

      F1m3Q = 
     >   SMPSNF(G1M3Q1,x + 0.00000000001d0*(a5-x), 
     >	               x + 0.00000001d0*(a5-x),4, Er, Ir)
     > + SMPSNF(G1M3Q1,x + 0.00000001d0*(a5-x),
     >                 x + 0.0000001d0*(a5-x),4, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.0000001d0*(a5-x),
     >                 x + 0.000001d0*(a5-x),4, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.000001d0*(a5-x),
     >                 x + 0.00001d0*(a5-x),4, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.00001d0*(a5-x),
     >                 x + 0.0001d0*(a5-x),4, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.0001d0*(a5-x),
     >                 x + 0.001d0*(a5-x),4, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.001d0*(a5-x),
     >                 x + 0.01d0*(a5-x),4, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.01d0*(a5-x),
     >                 x + 0.1d0*(a5-x),8, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.1d0*(a5-x),
     >                 x + 0.9d0*(a5-x),16, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.9d0*(a5-x),
     >                 x + 0.99d0*(a5-x),8, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.99d0*(a5-x),
     >                 x + 0.999d0*(a5-x),4, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.999d0*(a5-x),
     >                 x + 0.9999d0*(a5-x),4, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.9999d0*(a5-x),
     >                 x + 0.99999d0*(a5-x),4, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.99999d0*(a5-x),
     >                 x + 0.999999d0*(a5-x),4, Er, Ir) 
     > + SMPSNF(G1M3Q1,x + 0.999999d0*(a5-x),a5,4, Er, Ir)


      F2 = F2Qb - F2Qa
      F1 = F2 - F2m1Q                ! F1 = F2 - (F2-F1)
      SfnG(1) = F1 / 2d0             ! SfnG_1 = F1/2
      SfnG(3) = F2m1Q / 2d0          ! SfnG_3 = Flong = (F2-F1)/2
      SfnG(2) = 0d0        ! Zero-mass case has no chirality-mixing terms
      SfnG(4) = 0d0        ! Zero-mass case has no chirality-mixing terms
      F3 =  F1 - F1m3Q                                   
      SfnG(5) = F3 /2d0              ! SfnG_5 = F3/2 
      
                                     ! Log(Q/mu) term:
      If (Abs(Q/aMu - 1d0) > Smll) Then
                         ! delta(1-z) term
                        ! here we use the analytic result for the  
                        ! integral of aPqq to speed up the computation  
cmarco11 Old implementation
c      Qmu1 =  ADZINT_CT(aPqq,  0d0, x, aErr, rErr, ER, IR, 1, 0)    
c     >     * Pdfnx                                          
      Qmu1 = (-0.5d0*x*(2d0 + x) - 2d0*Log(1d0 - x)) ! delta(1-z) term
     >     * Pdfnx                                          

c **      Qmu2 = ADZINT_CT(aPqqFx, X, 1d0, aErr, rErr, ER, IR, 1, 1)     ! regular term

	a5 = 1.d0 - 0.00000000001d0*(1.d0-x)

	Qmu2 = 
     >   SMPSNF(aPqqFx,x + 0.00000000001d0*(a5-x), 
     >	               x + 0.000001d0*(a5-x),4, Er, IR)
     > + SMPSNF(aPqqFx,x + 0.000001d0*(a5-x),
     >                 x + 0.00001d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aPqqFx,x + 0.00001d0*(a5-x),
     >                 x + 0.0001d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aPqqFx,x + 0.0001d0*(a5-x),
     >                 x + 0.001d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aPqqFx,x + 0.001d0*(a5-x),
     >                 x + 0.01d0*(a5-x),8, Er, IR) 
     > + SMPSNF(aPqqFx,x + 0.01d0*(a5-x),
     >                 x + 0.1d0*(a5-x),16, Er, IR) 
     > + SMPSNF(aPqqFx,x + 0.1d0*(a5-x),
     >                 x + 0.9d0*(a5-x),32, Er, IR) 
     > + SMPSNF(aPqqFx,x + 0.9d0*(a5-x),
     >                 x + 0.99d0*(a5-x),16, Er, IR) 
     > + SMPSNF(aPqqFx,x + 0.99d0*(a5-x),
     >                 x + 0.999d0*(a5-x),8, Er, IR) 
     > + SMPSNF(aPqqFx,x + 0.999d0*(a5-x),
     >                 x + 0.9999d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aPqqFx,x + 0.9999d0*(a5-x),
     >                 x + 0.99999d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aPqqFx,x + 0.99999d0*(a5-x),
     >                 x + 0.999999d0*(a5-x),4, Er, IR) 
     > + SMPSNF(aPqqFx,x + 0.999999d0*(a5-x),a5,4, Er, IR)
      
      OmgLo = (/5d-1, 0d0, 0d0, 0d0, 5d-1/)
      SfnG = SfnG + 2d0 * OmgLo * Log(Q/aMu) * (Qmu2 - Qmu1) 
      End If
              
      SfnG = SfnG * CF               ! QCD color factor

      End Subroutine SfZm1Q             !  *****************************
c=====                         ! ---- integrand for the \delta function term
      Function aC2Q1(z)
C                                     Local variable
      Real*8 z, aC2Q1
      
      aC2Q1 = ((1d0+ Z**2) *(LOG((1d0- Z)/Z) - 3d0/4d0) /(1d0 - Z))
     >           + (9d0+ 5d0*Z) / 4d0


      End Function aC2Q1                                    ! ----- *****
c=====                        ! ---- integrand for the F2 term
      Function aF2Q1 (z)
      
      Implicit None      
C                                          ! Global variables:
      Integer :: iHadn, iPartn
      Real*8  :: x, aMu, Pdfnx
      Common / CSfZm1Q / x, aMu, Pdfnx, iHadn, iPartn
C                                          ! Local variable
      Integer IR
      Real*8 :: z, Pdfnxz, aF2Q1
C                                          ! External functions
      Real*8 :: aC2Q1, PDF_CT
    
      Pdfnxz = PDF_CT(iHadn, iPartn, x/z, aMu, IR)

      aF2Q1 =   aC2Q1(z) * (Pdfnxz/z - Pdfnx)

      End Function aF2Q1                                    ! ----- *****
                       
c=====                 ! ---- integrand for the Flong term = (F2 - F1)/2
      Function G2M1Q1 (z)

      Implicit None      
C                                          ! Global variables:
      Integer :: iHadn, iPartn
      Real*8  :: x, aMu, Pdfnx
      Common / CSfZm1Q / x, aMu, Pdfnx, iHadn, iPartn
C                                          ! Local variable
      Integer IR
      Real*8 :: z, Pdfnxz, G2M1Q1 
C                                          ! External functions 
      Real*8 :: PDF_CT
      
      Pdfnxz = PDF_CT(iHadn, iPartn, x/z, aMu, IR)

      G2M1Q1 = 2d0 * Z * Pdfnxz / z
      
      End Function G2M1Q1                                  ! ----- *****
    
C=====            ! ---- integrand for the Fright term = (F1 - F3)/2
      Function G1M3Q1 (z)

      Implicit None
C                                          ! Global variables:
      Integer :: iHadn, iPartn
      Real*8  :: x, aMu, Pdfnx
      Common / CSfZm1Q / x, aMu, Pdfnx, iHadn, iPartn
C                                          ! Local variable
      Integer IR
      Real*8 z, Pdfnxz, G1M3Q1
C                                          ! External functions 
      Real*8 :: PDF_CT

      Pdfnxz = PDF_CT(iHadn, iPartn, x/z, aMu, IR)

      G1M3Q1 = (1d0 - Z) * Pdfnxz /z

      End Function G1M3Q1                                  ! ----- *****
C=====
      Function aPqq(z)

      Implicit None
      
      Real*8 :: z, aPqq
      
      aPqq = (1d0 + z**2) / (1d0 - z)
      
      End Function aPqq
C=====
      Function aPqqFx (z)
      
      Implicit None
C                                          ! Global variables:
      Integer :: iHadn, iPartn
      Real*8  :: x, aMu, Pdfnx
      Common / CSfZm1Q / x, aMu, Pdfnx, iHadn, iPartn
C                                          ! Local variable
      Integer IR
      Real*8 :: z, Pdfnxz, aPqqFx
C                                          ! External functions
      Real*8 :: aPqq, PDF_CT
    
      Pdfnxz = PDF_CT(iHadn, iPartn, x/z, aMu, IR)

      aPqqFx =   aPqq(z) * (Pdfnxz/z - Pdfnx)

      End Function aPqqFx
