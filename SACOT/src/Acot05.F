      SUBROUTINE Acot0Q (Ihad,Iprtn, F1M,F2M, X,Q, aMu, SfnG)
C-----------------------------------------------------------------------------
C      Computes LO Hadron Structure Functions in the general irreducible set based on
C        symmetry considerations (cf. DisPac05 notes, originated in ACOT-II paper).
C  SfnG(1): transverse,   parity conserving amp. ~ g_V^2 + g_A^2 (g_R^2 + g_L^2) (=
C  SfnG(2): transverse,   parity conserving amp. ~ g_V^2 - g_A^2 (   g_V*g_A   )
C  SfnG(3): longitudinal, parity conserving amp. ~ g_V^2 + g_A^2 (g_R^2 + g_L^2)
C  SfnG(4): longitudinal, parity conserving amp. ~ g_V^2 - g_A^2 (   g_V*g_A   )
C  SfnG(5): transverse,   parity violating  amp. ~    g_V*g_A    (g_R^2 - g_L^2)
C
C        (x, Q) are input kinematic variables;
C        aMu is the factorization scale, which does not have to be equal to Q.
C        F1m/F2m are the masses of the initial/final state quark masses.
C-----------------------------------------------------------------------------
       Implicit None
C                                              Dummy Arguments       
       Integer, Intent(In) :: Ihad, Iprtn
       Real*8,  Intent(In) :: F1M, F2M, X, Q, aMu
       Real*8,  Intent(Out):: SfnG(5)
C                                              Local variables
       Integer:: iRet, i
       Real*8 :: OmgLo(5), PdfTmp
C                                              External function
       Real*8 :: PDF_CT
C                                -------------

C                                              Compute the partonic irreducible amp's
       CALL Omg0 (Q, F1M, F2M, OmgLo)
C                                              PDF factor
       PdfTmp = PDF_CT(Ihad, Iprtn, X, aMu, iRet)
C                                              Result       
       SfnG = OmgLo * PdfTmp /2d0     ! 2d0 is the conversion factor due to def. of C_i

       do i = 1, 5
          If(isNaN(SfnG(i))) then
             print *,'fatal error #1: SfnG=',SfnG
             stop
          EndIf
       enddo
C                                              
       END Subroutine Acot0Q                ! *******************
C
      SUBROUTINE Omg0 (Q, F1M, F2M, OmgLo)
C-----------------------------------------------------------------------------
C      Irreducible LO quark-initiated amplitudes OmgLo from ACOT-I paper, 
C                       Simplified from Table 7.
C-----------------------------------------------------------------------------
      Implicit None
      
C                                                  Dummy variables      
      Real*8, Intent(In) :: F1M, F2M, Q
      Real*8, Intent(Out):: OmgLo(5)
C                                                  Local variables      
      Real*8 :: Delta, A, B, C, F1m2, F2m2, Q2, Del
      Real*8 :: delsq, tmp
      
c xx      DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2d0*(A*B+B*C+C*A))
C-----------------------------------------------------------------------------
      F1m2  = F1m**2
      F2m2  = F2m**2
      Q2    = Q**2
c xx      DEL   = DELTA(-Q2,F1M2,F2M2)
c replace delta by form that should be less sensitive to round-off (jcp 10/2013)
        A = -Q2
        B = F1M2
        C = F2M2
c first make A .ge. B .ge. C 
        If(B .gt. A) then
           tmp = A
           A = B
           B = tmp
        Endif
        If(C .gt. A) then
           tmp = A
           A = C
           C = tmp
        Endif
        If(C .gt. B) then
           tmp = B
           B = C
           C = tmp
        Endif
           
        delsq = (A-B-C)**2 - 4.d0*B*C
        if(delsq .ge. 0.d0) then
           del = sqrt(delsq)
        else
           print *,'Fatal delsq #1 =',delsq,' A,B,C=',A,B,C
           stop
        endif
      
C-----------------------------------------------------------------------------
      OmgLo(1) = (Q2+F1m2+F2m2)/DEL
      OmgLo(2) = -2d0 *F1m*F2m/DEL
      OmgLo(5) = 1d0              ! This corresponds to C_a=g2L^2-g2R^2 (opposite of ACOT-II def.)

      OmgLo(3) = ((F1m2+F2m2) + ((F1m2-F2m2)**2/Q2))/DEL
      OmgLo(4) = +2d0 *F1m*F2m/DEL

      Return
      End Subroutine Omg0                        ! *********************

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE Acot1G (iHdn, M2, M1, XX, QQ, vMu, SFnG)
C-----------------------------------------------------------------------------
C     Computes NLO Hadron general SF SfnG --- Gluon Fusion term
C        These results are independent of final-state quark flavors;
C
C        For notation, see comment lines in Subroutine Acot0Q
C  From ACOT-II paper, except for change of definitions noted in these comments.
C-----------------------------------------------------------------------------

      Use mDisPac
      Implicit None
      Real*8 :: tmp1, a5, ERR, smpsnf
      Integer :: IER
      Integer :: Num
C                                             Dummy Arguments       
      Integer, Intent(In) :: iHdn
      Real*8,  Intent(In) :: M1, M2, XX, QQ, vMu
      Real*8,  Intent(Out):: SFnG(5)
C                                             Global variables
      Integer :: iHad, Ired
      Real*8  :: x, Q, F1M, F2M, aMu
      Common / CxNloG / x, Q, F1M, F2M, aMu, iHad, Ired
C                                             Local variables
      Integer:: iRet, i
      Real*8 :: SFired(5), ErrEst, z0, z1
C                                             External function
      Real*8 :: PDF_CT, ADZINT_CT, Alpi, xNloG
      External xNloG
C                      ------- populate the common block for integrand


      x = XX
      Q = QQ
      F1m = M1
      F2m = M2
      aMu = vMu
      iHad = iHdn
C                                    Limits of convolution integration      
      z0 = x *( (F1M+F2M)**2 + Q**2 ) / Q**2
      z1 = 1d0
C                                            Test these limits
      IF (z0 >= z1) Then
        Print '(A / 2(1pE15.6))', 
     >        'z0 >= z1 in Acot1G convolution, z0, z1=', z0, z1
        Stop
      EndIf
C                   Calculate the convolution integrals for the Ireducible amps
      Do Ired = 1, 5

	a5 = 1.d0 - 0.00000000001d0*(1.d0-z0)

	SfIred(Ired) = 
     >   SMPSNF(xNloG,z0 + 0.00000000001d0*(a5-z0), 
     >	              z0 + 0.000001d0*(a5-z0),4, ERR, IER)
     > + SMPSNF(xNloG,z0 + 0.000001d0*(a5-z0),
     >                z0 + 0.00001d0*(a5-z0),4, ERR, IER) 
     > + SMPSNF(xNloG,z0 + 0.00001d0*(a5-z0),
     >                z0 + 0.0001d0*(a5-z0),4, ERR, IER) 
     > + SMPSNF(xNloG,z0 + 0.0001d0*(a5-z0),
     >                z0 + 0.001d0*(a5-z0),4, ERR, IER) 
     > + SMPSNF(xNloG,z0 + 0.001d0*(a5-z0),
     >                z0 + 0.01d0*(a5-z0),8, ERR, IER) 
     > + SMPSNF(xNloG,z0 + 0.01d0*(a5-z0),
     >                z0 + 0.1d0*(a5-z0),16, ERR, IER) 
     > + SMPSNF(xNloG,z0 + 0.1d0*(a5-z0),a5,48, ERR, IER) 

      EndDo
C                                        Put on the QCD color factor
      SfnG = SfIred * TR

      do i = 1, 5
         If(isNaN(SfnG(i))) then
            print *,'fatal error #2: SfnG=',SfnG
            stop
         EndIf
      enddo

      END Subroutine Acot1G                  ! *************************
C      
      Function xNloG(z)
C----------------------------------------------------------------------------
C     Integrand for the Acot1G convolution integral
C----------------------------------------------------------------------------
      Implicit None
      
      Real*8  :: xNloG
      Real*8, intent(in) :: z
C                                 Common variables
      Integer :: iHad, Ired
      Real*8  :: x, Q, F1M, F2M, aMu
C                                 Local variables
      Integer :: iRet
      Real*8  :: xHat, PdfTmp, Ftem, OmgNlo(5)
C                                 External function      
      Real*8  :: PDF_CT
            
      Common / CxNloG / x, Q, F1M, F2M, aMu, iHad, Ired

C                                                ! The gluon PDF
      PdfTmp  = PDF_CT(iHad, 0, z, aMu, iRet)
C                                                ! The Wilson Coeff.      
      xHat = x/z
      Call OmgG1 (xHat, Q, F1M, F2M, OmgNlo)
      Ftem = OmgNlo(Ired)                        ! Ired is the irreducible amp label
C                                                ! The integrand      
      xNloG = Ftem * PdfTmp / z
      
      END Function xNloG                        ! **************************
C
      SUBROUTINE OmgG1 (X, Q, F1M, F2M, OmgNlo)
C-----------------------------------------------------------------------------
C      Irreducible NLO gluon-initiated amplitudes Omga_G: Eqs.(18)-(19). ACOT-II
C        Needed in the NLO convolution integral
C-----------------------------------------------------------------------------
      Implicit None      
C                                                  Dummy variables      
      Real*8, Intent(In) :: F1M, F2M, x, Q
      Real*8, Intent(Out):: OmgNlo(5)
C                                                  Local variables      
      Real*8 :: Delta, A, B, C, F1m2, F2m2, Q2, Del, E1, E2, Eq,
     >  xMax, sMin, qMin, s, RS, tLog, uLog, k, p
      Real*8 :: delsq, tmp
      
      Delta(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2d0*(A*B+B*C+C*A))
C-----------------------------------------------------------------------------
      F1m2  = F1m**2
      F2m2  = F2m**2
      Q2    = Q**2
c xx      DEL   = DELTA(-Q2,F1M2,F2M2)
c replace delta by form that should be less sensitive to round-off (jcp 10/2013)
        A = -Q2
        B = F1M2
        C = F2M2
c first make |A| .ge |B| .ge. |C| 
        If(B .gt. A) then
           tmp = A
           A = B
           B = tmp
        Endif
        If(C .gt. A) then
           tmp = A
           A = C
           C = tmp
        Endif
        If(C .gt. B) then
           tmp = B
           B = C
           C = tmp
        Endif
           
        delsq = (A-B-C)**2 - 4.d0*B*C
        if(delsq .ge. 0.d0) then
           del = sqrt(delsq)
        else
           print *,'Fatal delsq #2 =',delsq,' A,B,C=',A,B,C
           stop
        endif
      
C-----------------------------------------------------------------------------

      xMax=Q2/((F1M+F2M)**2+Q2)
      IF((X.GE.1d0).OR.(X.GE.xMax)) RETURN

      sMin=(F1M+F2M)**2
      qMin=SQRT( sMin*X/(1d0-X) )
      IF(Q.LE.qMin) RETURN

C     Cm-Energy for the hard process for this SHAT=S
      S = Q2 * (1d0/X - 1d0)
      IF(S.LE.sMin) RETURN
      RS= SQRT(S)
 
c xx      Del  = Delta(S, F1m2, F2m2)
C-----------------------------------------------------------------------------
c replace delta by form that should be less sensitive to round-off (jcp 10/2013)
        A = S
        B = F1M2
        C = F2M2
c first make A .ge. B .ge. C 
        If(B .gt. A) then
           tmp = A
           A = B
           B = tmp
        Endif
        If(C .gt. A) then
           tmp = A
           A = C
           C = tmp
        Endif
        If(C .gt. B) then
           tmp = B
           B = C
           C = tmp
        Endif
           
        delsq = (A-B-C)**2 - 4.d0*B*C
        if(delsq .ge. 0.d0) then
           del = sqrt(delsq)
        else
           print *,'Fatal delsq #3 =',delsq,' A,B,C=',A,B,C
           stop
        endif
C-----------------------------------------------------------------------------

      tLog = - Log(4*F1m2*s/(s+F1m2-F2m2+Del)**2)   ! L_t of ACOT-II
      uLog = - Log(4*F2m2*s/(s-F1m2+F2m2+Del)**2)   ! L_u of ACOT-II

      k  = ( Q2 + S)          /(2d0*RS)             !
      p  = Del/(2d0*RS)                             !
      E1 = ( F1m2 - F2m2 + S) /(2d0*RS)             ! Eq.(16) & above
      E2 = (-F1m2 + F2m2 + S) /(2d0*RS)             !
      Eq = (-Q2 + S)          /(2d0*RS)             !

      OmgNlo(1) = +tLog * (1d0/2 + E1*(E1/k-1d0)/k)                       !
     >            +uLog * (1d0/2 + E2*(E2/k-1d0)/k)                       !
     >            - 2d0 *p *Eq**2 /(k**2 *RS)                             !
                                                                          !
      OmgNlo(2) = (2*F1M*F2M/(4*k**2*s))                                  !
     >          * ((tLog + uLog)*(-F1m2 - F2m2 + s) - 4*p*RS)             !   Eq.(18)
                                                                          !
      OmgNlo(5) = p*(F2m2 - F1m2)/(k**2*RS)                               !
     > + tLog * (1d0/2 + E1*(E1/k-1d0)/k + F1m2*(F2m2-F1m2)/(2*k**2*s))   !
     > - uLog * (1d0/2 + E2*(E2/k-1d0)/k + F2m2*(F1m2-F2m2)/(2*k**2*s))   !
                                                                          
      OmgNlo(3) = p*((F2m2-F1m2)**2 + Q2*(2*Q2-F1m2-F2m2))/(k**2*Q2*RS)   !
     > - (tLog - uLog)*Eq*(F2m2 - F1m2)/(k**2*RS)                         !
     > - (tLog + uLog)*(                                                  !
     >   (F1m2 + F2m2)*(Q2**2 -(F2m2-F1m2)**2 -2*s*k**2) /(4*s*Q2*k**2)   !
     >    - Eq *(F1m2 + F2m2 - (F2m2-F1m2)**2/Q2) /(2*k**2*RS)            !   Eq.(19)
     >    + F1m2*F2m2 /(k**2*s) )                                         !
                                                                          !
      OmgNlo(4) = F1M*F2M* ( 2 *p /(k**2*RS)                              !
     >          - (tLog + uLog) *(1/Q2 + (s -F1m2 -F2m2)/(2*k**2*s)) )    !

      End Subroutine OmgG1                    ! ************************
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE AcotSub (iHdn,iPartn, F1m,F2m, XX,Q, vMu, Subtr)
C-----------------------------------------------------------------------------
C     Computes singular part of NLO Gluon Fusion Structure Functions Acot05G
C        for one single quark flavor contribution;
C        If this quark flavor is "light", and it contributes to the summation
C        over the incoming partons in the factorization formula, then this contribution
C        must be subtracted.
C        For a given gauge-invariant pair of gluon-fusion diagrams, there are two of
C        singular contributions---one for each quark flavor (the same for NC case).
C        The calling program supplies g2chi, which contains all the flavor dependence;
C        summation over quark flavors is managed in the calling program.
C
C        For notation, see comment lines in Subroutine Acot0Q
C  From ACOT-II paper, except for change of definitions noted in these comments.
C-----------------------------------------------------------------------------

      Use mDisPac
      Implicit None
      Real*8 :: tmp1, smpsnf
	integer Num
C                                             Dummy Arguments
      Integer, Intent(In) :: iHdn, iPartn
      Real*8,  Intent(In) :: F1m, F2m, XX, Q, vMu
      Real*8,  Intent(Out):: Subtr(5)
C                                             Global variables
      Integer :: iHad
      Real*8  :: x, aMu
      Common / CxActSub / x, aMu, iHad
C                                             Local variables
      Integer:: iRet
      Real*8 :: OmgLo(5), ErrEst, z0, z1, PertPdf, aM1
C                                             External function
      Real*8 :: PDF_CT, ADZINT_CT, Alpi, xActSub
      External xActSub
C                      -----------------------
C                                              Set up the common block
      x = XX
      aMu = vMu
      iHad = iHdn
C                                               -----------------
C                                    Limits of convolution integration
      z0 = x
      z1 = 1d0
C                                            Test these limits
      IF (z0 >= z1) Then
        Print '(A / 2(1pE15.6))',
     >        'z0 >= z1 in Acot1Sub convolution, z0, z1=', z0, z1
        Stop
      EndIf
C                   Calculate the convolution integrals for the perturbative Pdf


c Simpson's rule is good enough for this one...
      PertPdf = SMPSNF(xActSub, z0, z1, 48, ErrEst, iRet)

C                                          ! Add the pertinent factors
      PertPdf = PertPdf * Log((aMu/F1m)**2)
C                                              ---------------------
C                                      Get the LO irreducible parton amp Omg
C                                      In the SACOT scheme, set quark-parton mass = 0
      aM1 = 0
      CALL Omg0 (Q, aM1, F2M, OmgLo)
C                                       ---- Put the two factors together, ----

      Subtr = OmgLo *PertPdf *TR

      END Subroutine AcotSub                  ! ************************

      Function xActSub(z)
C----------------------------------------------------------------------------
C     Integrand for the AcotSub convolution integral
C----------------------------------------------------------------------------
      Implicit None

      Real*8  :: xActSub
      Real*8, intent(in) :: z
C                                             Global variables
      Integer :: iHad
      Real*8  :: x, aMu
      Common / CxActSub / x, aMu, iHad
C                                             Local variables
      Integer :: iRet
      Real*8  :: GluPdf, SpltFn
C                                             External function
      Real*8  :: PDF_CT
C                             -------------
C                                           ! The gluon-->quark splitting fn
      SpltFn = ( z**2 + (1d0-z)**2 ) /2d0
C                                           ! The gluon PDF
      GluPdf  = PDF_CT(iHad, 0, x/z, aMu, iRet)
C                                           ! The integrand
      xActSub = SpltFn * GluPdf / z

      END Function xActSub                        ! ********************

