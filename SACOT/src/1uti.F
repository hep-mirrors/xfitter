      FUNCTION SMPSNF(FN, A, B, NX, ERR, IER)
c jcp: IER is not used -- should remove it here and in all calls!
C                                                  -=-=- smpsnf
C
C                       Does integral of FN(X)*dx from A TO B by SIMPSON'S METHOD
C
C                       Double precision version of SMPSN
C
C                       Input:          External function:      FN
C                                       Lower limit      :      A
C                                       Upper limit      :      B
C                                       Number of points :      Nx
C
C                       Uses (Nx-1) evenly spaced intervals.
C
C                       Output:         error estimate:         ERR
C                                       error code    :         IER
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      PARAMETER (MXPT = 1000)
C
      DIMENSION X(MXPT)
      external fn

C
      IF (NX .EQ. 1) Then
         SmpsnF = 0
         RETURN
      EndIf

      IF (NX .LT. 0 .OR. NX .GT. MXPT) GOTO 99
C
      DX = (B - A) / (NX-1)
      IF (DX .LE. 0) THEN
        WRITE (NOUT, *) 'DX .LE. 0 in SMPSNF, DX =', DX
        SMPSNF = 0
        RETURN
      ENDIF
C
      DO 10 I = 1, NX
      X(I) = (A*(NX-I) + B*(I-1)) / (NX-1)
   10 CONTINUE
C
      IF (NX .GT. 4) GOTO 50
C
      GOTO (20, 30, 40), NX-1
   20 SMPSNF = (FN(X(1)) + FN(X(2))) * DX / 2D0
      RETURN
   30 SMPSNF = (FN(X(1)) + 4D0 * FN(X(2)) + FN(X(3))) * DX / 3D0
      RETURN
   40 SMPSNF = (( FN(X(1)) + 4D0 * FN(X(2)) +     FN(X(3))) / 3D0
     > + (-FN(X(2)) + 8D0 * FN(X(3)) + 5D0 * FN(X(4))) / 12D0 ) * DX
      RETURN
C
   50 SE = FN(X(2))
      SO = 0
      NM1 = NX - 1
      DO 60 I = 4, NM1, 2
      IM1 = I - 1
      SE = SE + FN(X(I))
      SO = SO + FN(X(IM1))
   60 CONTINUE
      MS = MOD (NX, 2)
      IF (MS .EQ. 1) THEN
        SMPSNF = (FN(X(1)) + 4D0 * SE + 2D0 * SO + FN(X(NX))) * DX / 3D0
        TRPZ = (FN(X(1)) + 2D0 * (SE + SO) + FN(X(NX))) * DX / 2D0
      ELSE
        SMPSNF =(FN(X(1)) + 4D0 * SE + 2D0 * SO + FN(X(NM1))) * DX / 3D0
     > +(-FN(X(NM1-1)) + 8D0 * FN(X(NM1)) + 5D0 * FN(X(NX))) * DX / 12D0
        TRPZ = (FN(X(1)) + 2D0 * (SE + SO + FN(X(NM1))) + FN(X(NX)))
     >          * DX / 2D0
      ENDIF
C
      ERR = SMPSNF - TRPZ
C
      RETURN
C
   99 WRITE (NOUT, 999) NX
  999 FORMAT (/ 5X, 'NX = ', I6,
     >              'out of range in SIMPSON INTEGRATION ROUTINE')
      STOP
C                        ****************************
      END


      
      function xLi(n,x)
C                                                   -=-=- xli
CPN-------------------------------------------------------------------------
C     Polylogarithm function for real x,-1<x<1.For n=2 the precision is about
C     12 significant figures. For n>2 the precision is smaller, but still
C     better than 4 digits.
C     Dilogarithm part is the translation of Carl Schmidt's C code.

c jcp: Some double precision typing added.  But overall accuracy not checked.
c jcp: A faster version might be available from K.S. K"olbig, J.A. Mignaco and E. Remiddi, BIT 10 (1970) 38
c jcp: Or from CERNLIB dilog64.f, cgplg64.f


      implicit NONE

      integer NCUT, i,n,m3
      double precision xLi,Out,x,pi2by6,zeta3,c1,c2
      double precision r,xt,L,xln1m

      parameter (m3=8)
      dimension c1(2:m3),c2(2:m3)

      data NCUT/27/
      data c1/0.75d0,-0.5833333333333333d0,0.4548611111111111d0,
     >        -0.3680555555555555d0,0.3073611111111111d0,
     >        -0.2630555555555555d0,0.2294880243764172d0/

      data c2/-0.5d0,0.5d0,-0.4583333333333333d0,0.416666666666666d0,
     >        -0.3805555555555555d0,0.35d0,-0.3241071428571428d0/

      data zeta3,pi2by6 /1.202056903159594285d0,1.644934066848226436d0/
      L=0.0
      i=0
      r=1.0d0

C     Check if x lies in the correct range
      if (abs(x).gt.r) then
        print *,'xLi: x out of range (-1,1) , x=',x
        STOP
      endif

      if (n.lt.0) then
       print *,'The polylogarithm Li undefined for n=',n
       STOP
      elseif (n.eq.0) then
       Out=x/(1d0-x)
      elseif (n.eq.1) then
       Out=-dlog(1-x)
      elseif (n.eq.2) then
                                                !Calculate dilogarithm
                                                !separately for x<0.5 and x>0.5
      if (x.ge.(-0.5d0).and.x.le.0.5d0) then

         do while(i.le.NCUT)
          i=i+1
          r=r*x

          L=L+r/i/i
         enddo
         Out=L
       elseif (x.eq.0) then
         Out=0d0
       elseif(x.gt.0.5d0) then !n.eq.2,x>0.5
         xt = 1.0d0-x
         L = pi2by6 - dlog(x)*dlog(xt)

         do while(i.le.NCUT)
          i=i+1
          r=r*xt

          L=L-r/i/i
         enddo
         Out=L
       elseif (x.lt.(-0.5d0)) then
         xt=-x/(1d0-x)
         L=-0.5d0*dlog(1.d0-x)**2

         do while (i.le.NCUT)
          i=i+1
          r=r*xt
          L=L-r/i/i
         enddo
         Out=L
       endif
      elseif (n.eq.3.and.x.ge.0.8) then !use the expansion of Li3 near x=1
       L=zeta3+pi2by6*dlog(x)
       xt=(1d0-x)
       xln1m=dlog(xt)

       do i=2,m3
        L=L+(c1(i)+c2(i)*xln1m)*xt**i
       enddo
       Out=L
      else !n>3 or x=3,x<0.8

         do while(i.le.NCUT)
          i=i+1
          r=r*x

          L=L+r/dble(i)**dble(n)
         enddo
         Out=L

      endif

      xLi=Out
C                      **********************
      End ! xLi

      FUNCTION ADZINT_CT (F, A, B, AERR, RERR, ERREST, IER, IACTA
     >                   , IACTB)
C                                                   -=-=- adzint

C===========================================================================
C GroupName: Adzint
C Description: adaptive integration
C ListOfFiles: adzint adzcal adzspl intusz SGLINT_CT totalz
C===========================================================================
C                                  Authors: Wu-Ki Tung and John C. Collins
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Adzint.f,v 1.1 97/12/21 21:19:04 wkt Exp $
C #Log: Adzint.f,v $
c Revision 1.1  97/12/21  21:19:04  wkt
c Initial revision
c 

C     FUNCTION   ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZSPL (F, I, IER)
C     SUBROUTINE ADZCAL (F,I)
C     SUBROUTINE SGLINT_CT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOTALZ
C     FUNCTION   INTUSZ (X, FX)
C
C     COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), NUMINT,
C    > ICTA, ICTB, FA, FB, IB
C                        ------------------------

C     Adaptive integration routine which allows the integrand to be 
C     indeterminant at the lower and/or the upper ends of integration. 

C     Can self-adjust to any integrable singularity at the ends and compute 
C     the closest approximant, hence achieve the required accuracy efficiently
C     (provided the switch(s) IACTA (IACTB) are set to 2).
 
C     Input switches for end-treatment:
C        IACTA = 0 :   Use closed lower-end algorithm 
C                1 :   Open lower-end -- use open quadratic approximant
C                2 :   Open lower-end -- use adaptive singular approximant

C        IACTB = 0, 1, 2   (same as above, for the upper end)
 
C                Integral of F(X) from A to B, with error
C                less than ABS(AERR) + ABS(RERR*INTEGRAL)
C                Best estimate of error returned in ERREST.
CError code is IER: 0 :  o.k.
C                1 :  maximum calls to function reached before the 
C                     error criteria are met;
C                2 :  IACTA out of range, set to 1;
C                3 :  IACTB out of range, set to 1.
C                4 :  Error on Limits : B < A ; zero result returned.
C                5 :  Range of integration DX zero or close to roundoff
C                     returns DX * F(A+DX/2)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      DATA SMLL,Sml / 1d-20,1d-12 /
      data fac0 /2d0/ 

      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in ADZINT!'
        
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call', 
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB
 
      DDX = B - A
      If (DDX .Le. 0D0) Then
        ADZINT_CT = 0D0
        Ier = 4
        If (DDX .Lt. 0D0) 
     >     Print '(/A/)', 'B < A in AdzInt; check limits!!'
        Return
      ElseIf (DDX .Le. Sml) Then
        ADZINT_CT = F(A + DDX/2) * DDX
        Ier = 5
        Return
      EndIf

      NUMINT = 3
      DX = DDX/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A 
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1))
             ELSE 
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2.)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I))
             ELSE
               IB = I
               FB = F(B-DX/2.)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I))
          ENDIF
          CALL ADZCAL_CT(F,I)
   10     CONTINUE
       CALL TOTALZ_CT
C                                                   Adaptive procedure:
cLai revised to ensure the accuracy requested and improve the speed
       facMax= fac0*MAXINT
       Do
 20       TARGET = ABS(AERR) + ABS(RERR * RES)  !targeted total error
          if(ERS <= TARGET .or. Ier>0) goto 30  
          NUMOLD = NUMINT
          facNum = Sqrt(dble(NUMOLD)*fac0)
          Do
             DO I = 1, NUMOLD
                IF (ERR(I)*facNum > TARGET) CALL ADZSPL_CT(F,I,IER)
             Enddo              ! I
             if(NUMOLD .ne. NUMINT) then  !additional intervals were added
                goto 20
             else
                facNum=facNum*fac0
                if(facNum>facMax) goto 30
             endif
          Enddo !do loop if NUMOLD==NUMINT
       enddo !do loop if .not.(ERS<=Target .or. Ier>0)
 30    continue
      ADZINT_CT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE ADZCAL_CT (F,I)
C                                                   -=-=- adzcal
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB
 
      SAVE / ADZWRK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.
     
      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL SGLINT_CT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGLINT_CT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                       Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF
 
      RESULT(I) = TEM         
      ERR   (I) = ABS (ER)
 
      RETURN
C                        ****************************
      END

      SUBROUTINE SGLINT_CT (IACT, F1, F2, F3, DX, FINT, ESTER)
C                                                   -=-=- sglint

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
C                                                         Use quadratic formula
      TEM = DX * (4.*F1 + 3.*F2 + 2.*F3) / 9.
C                 Error est based on Diff between quadratic and linear integrals
      ER  = DX * (4.*F1 - 6.*F2 + 2.*F3) / 9.

C                          Invoke adaptive singular parametrization if IACT = 2
C                      Algorithm is based on the formula F(x) = AA + BB * x **CC
C                 where AA, BB & CC are determined from F(Dx/4), F(Dx/2) & F(Dx)

      IF (IACT .EQ. 2) THEN
          T1 = F2 - F1
          T2 = F3 - F2
          IF (T1*T2 .LE. 0.) GOTO 7
          T3  = T2 - T1
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -0.8D0)  GOTO 7
          BB  = T1**2 / T3
          AA  = (F1*F3 - F2**2) / T3
C                                          Estimated integral based on A+Bx**C
          TMP = DX * (AA + BB* 4.**CC / (CC + 1.))
C                                       Error estimate based on the difference
          ER = TEM - TMP
C                                              Use the improved integral value
          TEM= TMP 
      ENDIF

    7 FINT = TEM
      ESTER= ER
      RETURN
C                        ****************************
      END

      SUBROUTINE TOTALZ_CT
C                                                   -=-=- totalz
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

      SUBROUTINE ADZSPL_CT (F, I, IER)
C                                                   -=-=- adzspl
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      DATA TINY / 1.D-20 /
     
      IF (NUMINT .GE. MAXINT)  THEN
          IER = 1
          RETURN
          ENDIF
      NUMINT = NUMINT + 1
C                                                         New interval NUMINT
      IF (I .EQ. IB) IB = NUMINT
      U(NUMINT) = (U(I) + V(I)) / 2.
      V(NUMINT) = V(I)
 
      FU(NUMINT) = FW(I)
      FV(NUMINT) = FV(I)
C                                                             New interval I
       V(I) =  U(NUMINT)
      FV(I) = FU(NUMINT)
C                                                    Save old Result and Error
      OLDRES = RESULT(I)
      OLDERR = ERR(I)
     
      CALL ADZCAL_CT (F, I)
      CALL ADZCAL_CT (F, NUMINT)
C                                                               Update result
      DELRES = RESULT(I) + RESULT(NUMINT) - OLDRES
      RES = RES + DELRES
C                                  Good error estimate based on Simpson formula
      GODERR = ABS(DELRES) 
C                                                             Update new global 
      ERS = ERS + GODERR - OLDERR
C                                  Improve local error estimates proportionally
      SUMERR = ERR(I) + ERR(NUMINT)
      IF (SUMERR .GT. TINY) THEN
         FAC = GODERR / SUMERR 
      ELSE
         FAC = 1.
      ENDIF
      
      ERR(I)      = ERR(I) * FAC
      ERR(NUMINT) = ERR(NUMINT) * FAC
 
      RETURN
C                        ****************************
      END

      FUNCTION SMPSN2(FN, A, B, NX, ERR, IER)
c jcp: IER is not used -- should remove it here and in all calls!
C                                                  -=-=- smpsn2
C
C                       Does integral of FN(X)*dx from A TO B by SIMPSON'S METHOD
C
C                       Double precision version of SMPSN
C
C                       Input:          External function:      FN
C                                       Lower limit      :      A
C                                       Upper limit      :      B
C                                       Number of points :      Nx
C
C                       Uses (Nx-1) evenly spaced intervals.
C
C                       Output:         error estimate:         ERR
C                                       error code    :         IER
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      PARAMETER (MXPT = 1000)
C
      DIMENSION X(MXPT)
      external fn

C
      IF (NX .EQ. 1) Then
         Smpsn2 = 0
         RETURN
      EndIf

      IF (NX .LT. 0 .OR. NX .GT. MXPT) GOTO 99
C
      DX = (B - A) / (NX-1)
      IF (DX .LE. 0) THEN
        WRITE (NOUT, *) 'DX .LE. 0 in SMPSN2, DX =', DX
        SMPSN2 = 0
        RETURN
      ENDIF
C
      DO 10 I = 1, NX
      X(I) = (A*(NX-I) + B*(I-1)) / (NX-1)
   10 CONTINUE
C
      IF (NX .GT. 4) GOTO 50
C
      GOTO (20, 30, 40), NX-1
   20 SMPSN2 = (FN(X(1)) + FN(X(2))) * DX / 2D0
      RETURN
   30 SMPSN2 = (FN(X(1)) + 4D0 * FN(X(2)) + FN(X(3))) * DX / 3D0
      RETURN
   40 SMPSN2 = (( FN(X(1)) + 4D0 * FN(X(2)) +     FN(X(3))) / 3D0
     > + (-FN(X(2)) + 8D0 * FN(X(3)) + 5D0 * FN(X(4))) / 12D0 ) * DX
      RETURN
C
   50 SE = FN(X(2))
      SO = 0
      NM1 = NX - 1
      DO 60 I = 4, NM1, 2
      IM1 = I - 1
      SE = SE + FN(X(I))
      SO = SO + FN(X(IM1))
   60 CONTINUE
      MS = MOD (NX, 2)
      IF (MS .EQ. 1) THEN
        SMPSN2 = (FN(X(1)) + 4D0 * SE + 2D0 * SO + FN(X(NX))) * DX / 3D0
        TRPZ = (FN(X(1)) + 2D0 * (SE + SO) + FN(X(NX))) * DX / 2D0
      ELSE
        SMPSN2 =(FN(X(1)) + 4D0 * SE + 2D0 * SO + FN(X(NM1))) * DX / 3D0
     > +(-FN(X(NM1-1)) + 8D0 * FN(X(NM1)) + 5D0 * FN(X(NX))) * DX / 12D0
        TRPZ = (FN(X(1)) + 2D0 * (SE + SO + FN(X(NM1))) + FN(X(NX)))
     >          * DX / 2D0
      ENDIF
C
      ERR = SMPSN2 - TRPZ
C
      RETURN
C
   99 WRITE (NOUT, 999) NX
  999 FORMAT (/ 5X, 'NX = ', I6,
     >              'out of range in SIMPSON INTEGRATION ROUTINE')
      STOP
C                        ****************************
      END

      Real*8 Function snp12(x)

c      Use mDisPac

c Nielsen function
      Real*8, intent(in)::x

      Real*8:: GausInt, snp12int
      external snp12int
      Real*8, Parameter :: D0=0d0,D1=1d0
      Integer :: IER
      Real*8 :: ERREST, XX
c v05p original      Real*8  :: AERR = 1d-16, RERR = 1d-4 
      Real*8  :: AERR = 1d-16, RERR = 1d-10 
cmarco rErr must be set to 1d-5 or 1d-6 to avoid spikes in F2b
c jcp rErr=1d-10 is OK here -- not measurably slower.
      common /snp/ XX

      XX=X
      snp12= .5d0* GausInt
     >  (snp12int, D0,D1, AERR,RERR,ERREST,IER)
      end function snp12

      Real*8 function snp12int(z)
      Real*8, intent(in):: z
      Real*8 :: XX,dlxz
      common /snp/ XX
      dlxz=dlog(1d0-XX*z)
      snp12int=dlxz*dlxz/z
      end function snp12int

      FUNCTION GausInt(F,XL,XR,AERR,RERR,ERR,IRT)
C                                                   -=-=- gausint

c *** Note, the coefficients appear to be single precision ***
C                                           Adptive Gaussian integration
C     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
        DIMENSION XLIMS(100), R(93), W(93)
        INTEGER PTR(4),NORD(4),NIN,NOUT,NWRT

        external f

        COMMON/IOUNIT/NIN,NOUT,NWRT
        DATA PTR,NORD/4,10,22,46,  6,12,24,48/
        DATA R/.2386191860d0,.6612093865d0,.9324695142d0,
     1 .1252334085d0,.3678314990d0,.5873179543d0,.7699026742d0,
     1 .9041172563d0,
     1 .9815606342d0,.0640568929d0,.1911188675d0,.3150426797d0,
     1 .4337935076d0,
     1 .5454214714d0,.6480936519d0,.7401241916d0,.8200019860d0,
     1 .8864155270d0,
     1 .9382745520d0,.9747285560d0,.9951872200d0,.0323801710d0,
     1 .0970046992d0,
     1 .1612223561d0,.2247637903d0,.2873624873d0,.3487558863d0,
     1 .4086864820d0,
     1 .4669029048d0,.5231609747d0,.5772247261d0,.6288673968d0,
     1 .6778723796d0,
     1 .7240341309d0,.7671590325d0,.8070662040d0,.8435882616d0,
     1 .8765720203d0,
     1 .9058791367d0,.9313866907d0,.9529877032d0,.9705915925d0,
     1 .9841245837d0,
     1 .9935301723d0,.9987710073d0,.0162767488d0,.0488129851d0,
     1 .0812974955d0,
     1 .1136958501d0,.1459737146d0,.1780968824d0,.2100313105d0,
     1 .2417431561d0,
     1 .2731988126d0,.3043649444d0,.3352085229d0,.3656968614d0,
     1 .3957976498d0,
     1 .4254789884d0,.4547094222d0,.4834579739d0,.5116941772d0,
     1 .5393881083d0,
     1 .5665104186d0,.5930323648d0,.6189258401d0,.6441634037d0,
     1 .6687183100d0,
     1 .6925645366d0,.7156768123d0,.7380306437d0,.7596023411d0,
     1 .7803690438d0,
     1 .8003087441d0,.8194003107d0,.8376235112d0,.8549590334d0,
     1 .8713885059d0,
     1 .8868945174d0,.9014606353d0,.9150714231d0,.9277124567d0,
     1 .9393703398d0,
     1 .9500327178d0,.9596882914d0,.9683268285d0,.9759391746d0,
     1 .9825172636d0,
     1 .9880541263d0,.9925439003d0,.9959818430d0,.9983643759d0,
     1 .9996895039d0/
        DATA W/.4679139346d0,.3607615730d0,.1713244924d0,
     1 .2491470458d0,.2334925365d0,.2031674267d0,.1600783285d0,
     1 .1069393260d0,
     1 .0471753364d0,.1279381953d0,.1258374563d0,.1216704729d0,
     1 .1155056681d0,
     1 .1074442701d0,.0976186521d0,.0861901615d0,.0733464814d0,
     1 .0592985849d0,
     1 .0442774388d0,.0285313886d0,.0123412298d0,.0647376968d0,
     1 .0644661644d0,
     1 .0639242386d0,.0631141923d0,.0620394232d0,.0607044392d0,
     1 .0591148397d0,
     1 .0572772921d0,.0551995037d0,.0528901894d0,.0503590356d0,
     1 .0476166585d0,
     1 .0446745609d0,.0415450829d0,.0382413511d0,.0347772226d0,
     1 .0311672278d0,
     1 .0274265097d0,.0235707608d0,.0196161605d0,.0155793157d0,
     1 .0114772346d0,
     1 .0073275539d0,.0031533461d0,.0325506145d0,.0325161187d0,
     1 .0324471637d0,
     1 .0323438226d0,.0322062048d0,.0320344562d0,.0318287589d0,
     1 .0315893308d0,
     1 .0313164256d0,.0310103326d0,.0306713761d0,.0302999154d0,
     1 .0298963441d0,
     1 .0294610900d0,.0289946142d0,.0284974111d0,.0279700076d0,
     1 .0274129627d0,
     1 .0268268667d0,.0262123407d0,.0255700360d0,.0249006332d0,
     1 .0242048418d0,
     1 .0234833991d0,.0227370697d0,.0219666444d0,.0211729399d0,
     1 .0203567972d0,
     1 .0195190811d0,.0186606796d0,.0177825023d0,.0168854799d0,
     1 .0159705629d0,
     1 .0150387210d0,.0140909418d0,.0131282296d0,.0121516047d0,
     1 .0111621020d0,
     1 .0101607705d0,.0091486712d0,.0081268769d0,.0070964708d0,
     1 .0060585455d0,
     1 .0050142027d0,.0039645543d0,.0029107318d0,.0018539608d0,
     1 .0007967921d0/
c **    DATA TOLABS,TOLREL,NMAX/1.D-35,5.D-4,100/
c remove values that are reset below (jcp) DATA TOLABS,TOLREL,NMAX/1.D-99,1.D-6,100/
        DATA NMAX/100/
C
C
        TOLABS=AERR
        TOLREL=RERR
     
        GausInt=0.
        NLIMS=2
        XLIMS(1)=XL
        XLIMS(2)=XR
C
10      AA=(XLIMS(NLIMS)-XLIMS(NLIMS-1))/2D0
        BB=(XLIMS(NLIMS)+XLIMS(NLIMS-1))/2D0
        TVAL=0.
        DO 15 I=1,3
15      TVAL=TVAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
        TVAL=TVAL*AA
        DO 25 J=1,4
        VAL=0.
        DO 20 I=PTR(J),PTR(J)-1+NORD(J)
20      VAL=VAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
        VAL=VAL*AA
        TOL=MAX(TOLABS,TOLREL*ABS(VAL))
        IF (ABS(TVAL-VAL).LT.TOL) THEN
                GausInt=GausInt+VAL
                NLIMS=NLIMS-2
                IF (NLIMS.NE.0) GO TO 10
                RETURN
                END IF
25      TVAL=VAL
        IF (NMAX.EQ.2) THEN
                GausInt=VAL
                RETURN
                END IF
        IF (NLIMS.GT.(NMAX-2)) THEN
                WRITE(NOUT,50) GausInt,NMAX,BB-AA,BB+AA
                RETURN
                END IF
        XLIMS(NLIMS+1)=BB
        XLIMS(NLIMS+2)=BB+AA
        XLIMS(NLIMS)=BB
        NLIMS=NLIMS+2
        GO TO 10
C
50      FORMAT (' GausInt FAILS, GausInt,NMAX,XL,XR=',G15.7,I5,2G15.7)
C                        ****************************
        END
 
C                                                          =-=-= Fintrp

      FUNCTION ZBRNT(FUNC, X1, X2, TOLIN, IRT)
C                                                   -=-=- zbrnt
 
C                          Return code  IRT = 1 : limits do not bracket a root;
C                                             2 : function call exceeds maximum
C appears to be Wijngaarden-Dekker-Brent algorithm zbrent from Numerical Recipes.

c www.ulib.org/webRoot/Books/Numerical_Recipes/bookfpdf/f9-3.pdf 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ITMAX = 1000, EPS = 3.E-12)
      external func

      IRT = 0
      TOL = ABS(TOLIN)    !modified so it can't change input TOLIN (jcp 11/4/02)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.)  THEN
        PRINT *, 'Root must be bracketed for ZBRNT. Set = 0'
        IRT = 1
        ZBRNT=0.
        PRINT *,'ZBRNT fatal error'     !make it a fatal error
        stop

c **    RETURN

      ENDIF
      C=B     !Previously missing -- would matter if FC=0 -- (jcp 11/12/01)
      FC=FB
      DO 11 ITER=1,ITMAX

        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRNT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      PRINT *, 'Warning: ZBRNT exceeding maximum iterations.'
      IRT = 2
      ZBRNT=B
      RETURN
      END

      SUBROUTINE WARNR (IWRN, NWRT, MSG, NMVAR, VARIAB,
C                                                   -=-=- warnr
     >                  VMIN, VMAX, IACT)
 
C      Subroutine to handle warning messages.  Writes the (warning) message
C      and prints out the name and value of an offending variable to SYS$OUT
C      the first time, and to output file unit # NWRT in subsequent times.
C      
C      The switch IACT decides whether the limits (VMIN, VMAX) are active or
C      not.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
 
      CHARACTER*(*) MSG, NMVAR

      Save Iw

      Data Nmax / 100 /

      IW = IWRN
      VR = VARIAB
 
      IF  (IW .EQ. 0) THEN
         PRINT '(1X, A/1X,2A,1PD16.7/A,I4)', MSG, NMVAR, ' = ', VR,
     >         ' For all warning messages, check file unit #', NWRT
         IF (IACT .EQ. 1) THEN
         PRINT       '(A/2(1PE15.4))', ' The limits are: ', VMIN, VMAX
         WRITE (NWRT,'(A/2(1PE15.4))') ' The limits are: ', VMIN, VMAX
         ENDIF
      ENDIF
 
      If (Iw .LT. Nmax) Then
         WRITE (NWRT,'(I5, 2A/1X,2A,1PD16.7)') IW, '   ', MSG,
     >                  NMVAR, ' = ', VR
      Elseif (Iw .Eq. Nmax) Then
         Print '(/A/)', '!!! Severe Warning, Too many errors !!!'
         Print '(/A/)', '    !!! Check The Error File !!!'
         Write (Nwrt, '(//A//)')
     >     '!! Too many warnings, Message suppressed from now on !!'
      Endif

      IWRN = IW + 1
 
      RETURN
C                  ****************************
      END

C                                                          =-=-= Inputs
      SUBROUTINE WARNI (IWRN, NWRT, MSG, NMVAR, IVAB,
C                                                   -=-=- warni
     >                  IMIN, IMAX, IACT)

C     t++++++++++++++++++++++++++++++++     Routines to handle Warnings
C                                                  Integer version
      CHARACTER*(*) MSG, NMVAR

      Save Iw

      Data Nmax / 100 /

      IW = IWRN
      IV = IVAB
      
      IF  (IW .EQ. 0) THEN
         PRINT '(1X,A/1X, 2A,I10 /A,I4)', MSG, NMVAR, ' = ', IV,
     >         ' For all warning messages, check file unit #', NWRT
         IF (IACT .EQ. 1) THEN
         PRINT       '(A/2I10)', ' The limits are: ', IMIN, IMAX
         WRITE (NWRT,'(A/2I10)') ' The limits are: ', IMIN, IMAX
         ENDIF
      ENDIF

      If (Iw .LT. Nmax) Then
         WRITE (NWRT,'(1X,A/1X,2A, I10)') MSG, NMVAR, ' = ', IV
      Elseif (Iw .Eq. Nmax) Then
         Print '(/A/)', '!!! Severe Warning, Too many errors !!!'
         Print '(/A/)', '    !!! Check The Error File !!!'
         Write (Nwrt, '(//A//)')
     >     'Too many warnings, Message suppressed !!'
      Endif

      IWRN = IW + 1

      RETURN
C               *************************
      END

      