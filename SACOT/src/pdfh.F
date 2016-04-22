      FUNCTION PdfH (Ihadron, Iparton, X, Q)
C     -------------------------------------------------
C     Parton Distribution Functions inside Hadrons.

C     Revised 4/4/94 by HLL & WKT: 
C        PdfH retains its name and argument list for compatibility with all
C        existing programs which Call this function; PdfH switches between
C        target hadrons; 
C        PdfP is for proton target; it switches between different Iset's.
C
C     Ihadron = -1,  0,  1,  2,  3,  4,  5,  6 : 
C              pbar, n,  p,  D, Cu,  C, ~D, Fe 

C               5 is "isoscalar-corrected iron" hence = D

C         In all cases, adjust the Iparton label
C         to convert to the corresponding proton distribution which is
C         given in Function PdfP;

C     ---------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     Character Msg*80

      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx

      DATA
     1   HUGE, DUM,  D0m,   D1p 
     1 / 1D10, 0.D0, -1D-6, 1.000001D0 /
     1   IW1,  IW2    / 2*0  /
C                                                -------  Check x-range
      IF (X.LE.D0m .or. X.GT.D1p) Then
        Call WARNR(IW1,NWRT,'X out of range in PdfH.', 'X', X,
     >               D0, D1, 1)
        TEM = DUM
      EndIf
C                   --- Conversion of  Ihadron  to proton distributions,
C                                     if necessary ---- 
      If (Ihadron.gt.6 .or. Ihadron.lt.-1) then
         Call WARNI(IW, NWRT,
     >    'Only Ihardon=-1,0,1,2,...,6 (pbar,n,p,D,*,) are active',
     >    'Ihadron', Ihadron, -1,6,1)
       PdfH=0.D0
       Return
      Endif
C                                 
      Jp=Abs(Iparton)
      Neff = NFL(Q)
c nfl(q) returns the number of `light' flavors at scale Q - effective
      If ( Jp .gt. NEFF) then
c 		if Jp > Neff, then set PdfH=0 and return
c         Call WARNI(IW, NWRT,
c     >    'Iparton out of range',
c     >    'Iparton', Iparton, -Neff,Neff,1)
	 PdfH = 0D0
         Return
      Endif

      If (Jp.eq.1 .or. Jp.eq.2) Then
C                                   Use Isospin symmetry n<->p  == u<->d
         Ipartner=3-Jp
         If (Iparton.lt.0) Ipartner=-Ipartner

         If (Ihadron.eq.1) then
            Tem= PdfP(Iparton, X, Q)
         Elseif (Ihadron.eq.-1) then
            Tem= PdfP(-Iparton, X, Q)
         Elseif (Ihadron.eq.0) then
            Tem= PdfP(Ipartner, X, Q)
         Elseif (Ihadron.eq.2 .or.Ihadron.eq.4 .or.Ihadron.eq.5) then
            Tem=( PdfP(Iparton, X, Q)
     >           +PdfP(Ipartner, X, Q) )/2.D0
         Elseif (Ihadron.eq.3) then
            Tem=( 29.D0* PdfP(Iparton, X, Q)
     >           +35.D0* PdfP(Ipartner, X, Q) )/64.D0
         Elseif (Ihadron.eq.6) then
            Tem=( 26.D0* PdfP(Iparton, X, Q)
     >           +30.D0* PdfP(Ipartner, X, Q) )/56.D0
         Endif
      Else

         Tem= PdfP(Iparton, X, Q)

      Endif
C                                      --- Make sure PdfH >= 0 --------
C                                       (unless Iknl<0 - polarized pdf)
cpn11 Not true for heavy-quark PDFs at NNLO
c      IF (TEM .LT. D0 .and. Iknl .ge. 1) Then
c        IF (TEM .LT. D0m .AND. X .LE. 0.55D0) Then
c        Call WARNR(IW2,NWRT,'PdfH < 0; Set --> 0', 'PDF',TEM,D0,HUGE,1)
c        WRITE (NWRT, '(A, I5, 2(1PE15.3))')
c     >      ' Iparton, X, Q = ', Iparton, X, Q
c        EndIf
c        TEM = D0
c      EndIf
C                        -------- Return function value and error code
      PdfH = TEM

      RETURN
C                        ****************************
      END
