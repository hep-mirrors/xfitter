      Subroutine StEwCpl2
C                                                   -=-=- stewcpl2

C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: EwCpl2
C Description: Setup + Func. Ewk couplings^2 for DY: Scat. and Annih.
C ListOfFiles: stewcpl2 ewcplg2 
C ====================================================================
C Entry points (ewcpl2sc ewcpl2an)

C #Header: /Net/d2a/wkt/1hep/2ewk/RCS/EwCpl2.f,v 6.1 98/08/16 17:24:34 wkt Exp $
C #Log:	EwCpl2.f,v $
c Revision 6.1  98/08/16  17:24:34  wkt
c Re-organization; rationalization; initialization for DIS & DY
c 

C          Setup EW SQuaRed coupling constants for Vector Boson Production
C          in the Annihilation and Compton channels; put in Common /EwCplC/. 

C                      --  Unpolarized - (Gv**2 + Ga**2)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NFL = 6, NBN = 4)

      COMMON / EWCPLC / CPLANH(-NFL:NFL, NBN, -NFL:NFL), 
     >                  CPLSCT(-NFL:NFL, NBN)

      Data Lprt / 0 /

ccpy 03/14/2012
      ID_WZ=1
      Call SetEwk(ID_WZ)

      NFLT = NFLTOT()

      If (Lprt .ge. 1)   PRINT '(A/)', 
     >  'Squared EW Couplings for Annih. and Compton processes:'
      DO 4 IBN = 1, NBN
         If (Lprt .ge. 1)   PRINT '(A/9A9)', ' ANNIH:'
     >, 'cb','cb','db','ub','G','u','d','s','c'
      DO 5 IP1 = -NFLT, NFLT
         IF (IP1 .EQ. 0) GOTO 5
C                                                                Scattering
C                                                        Final state summed
         CPLSCT (IP1, IBN) = 0.
         DO 6 IP2 = 1, NFLT
           IF  (IP1 .GT. 0) THEN
             CPLV = GEWQT (IP1, IBN, 1, IP2)
             CPLA = GEWQT (IP1, IBN,-1, IP2)
           ELSE
             CPLV = GEWQT (IP2, IBN, 1,-IP1)
             CPLA = GEWQT (IP2, IBN,-1,-IP1)
           ENDIF
           CPL2 = CPLV**2 + CPLA**2
           CPLSCT (IP1, IBN) = CPLSCT (IP1, IBN) + CPL2
    6 CONTINUE
C                                                                Annihilation
         DO 7 IP2 = -NFLT, NFLT
         IF (IP2 .EQ. 0) GOTO 7
C                                  To annihilate, must have quark - antiquark
         IF (SIGN(1,IP2) .EQ. SIGN(1,IP1)) THEN
            CPLANH (IP1, IBN, IP2) = 0.0
            GOTO 7
         ENDIF

         IF (IP1 .GT. 0) THEN
           CPLV = GEWQT (IP1, IBN, 1,-IP2)
           CPLA = GEWQT (IP1, IBN,-1,-IP2)
         ELSE
           CPLV = GEWQT (IP2, IBN, 1,-IP1)
           CPLA = GEWQT (IP2, IBN,-1,-IP1)
         ENDIF

         CPL2 = CPLV**2 + CPLA**2
         CPLANH (IP1, IBN, IP2) = CPL2
    7 CONTINUE

      If (Lprt .ge. 1)
     $  PRINT '(9F8.3)', (CPLANH(IP1, IBN, IP2), IP2=-4,4)
    5 CONTINUE
      If (Lprt .ge. 1)
     $  PRINT '(/A,9F7.3/)', ' Compt.:', (CPLSCT(IP1, IBN), IP1=-4,4)
    4 CONTINUE

      Return
C                        ****************************
      End

