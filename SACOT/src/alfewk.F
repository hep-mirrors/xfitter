      Function ALFEWK (IBOSON)
C                                                   -=-=- alfewk

C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: EwCpl0
C Description: Function calls to extract the Ewk couplings + Running alfa_em
C ListOfFiles: alfewk alfaem
C ====================================================================
C Entry function not shown: (vbnmas swg2f gewlt gewlh gewqt gewqh) 

C #Header: /Net/d2a/wkt/1hep/2ewk/RCS/EwCpl0.f,v 6.1 98/08/16 17:21:35 wkt Exp $
C #Log:	EwCpl0.f,v $
c Revision 6.1  98/08/16  17:21:35  wkt
c Re-organization; rationalization; initialization for DIS & DY
c 

C     +++++++++++++++++++++++++++++  Functions to extract EWK coefficients

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (NSP = 2, NGN = 3, NBN = 4, NPOL = 3)
      PARAMETER (NFL = NSP * NGN)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / EWKPAR / WMS, ZMS, SWG2, ALFE, ALFEW(NBN)
      COMMON / KMATRX / VKM (NGN, NGN)

C            Left, Right, Vector, and Axial-vector couplings of Leptons & Quarks

      COMMON / EW1LCP / GLL(NSP,NBN,NSP), GLR(NSP,NBN,NSP), 
     >                  GLV(NSP,NBN,NSP), GLA(NSP,NBN,NSP)
      COMMON / EW1QCP / GQL(NSP,NBN,NSP), GQR(NSP,NBN,NSP), 
     >                  GQV(NSP,NBN,NSP), GQA(NSP,NBN,NSP)
      COMMON / EW2QCP / HQL(NFL,NBN,NFL), HQR(NFL,NBN,NFL), 
     >                  HQV(NFL,NBN,NFL), HQA(NFL,NBN,NFL)
 

C                                 "Alpha QED" for all Electo-Weak Vector Bosons
      ALFEWK = ALFEW(IBOSON)

      Return
C                        ****************************

      Entry VBNMAS (IBOSON)

      IF     (IBOSON .EQ. 1) THEN
        BM = 0.0
      ELSEIF (IBOSON .EQ. 2 .OR. IBOSON .EQ. 3) THEN
        BM = WMS
      ELSEIF (IBOSON .EQ. 4) THEN
        BM = ZMS
      ELSE
        WRITE (NOUT, '(2A, I5)') ' Vector Boson Index out of range',
     >  ' in VBNMAS; IBOSON =', IBOSON
      ENDIF

      VBNMAS = BM

      Return
C                        ****************************

      Entry SWG2F ()

      SWG2F = SWG2

      Return
C                        ****************************

      Entry GEWLT (IT1, IBS, IBT, IT2)

C           (G) Electro-Weak coupling for Leptons in the Tensor base.
C            -  -       -                 -              -

C               To be used in the calculation of general E-WK matrix elements.

C           IT1, IT2 = 1,2 :  Weak Isospin (T3) of the two leptons
C           IBS =  1 -4    :  Boson label (see SETEWK)
C           IBT =             Boson polarization (tensor) label
C                  1       :  vector
C                 -1       :  axial-vector
C                  other   :  illegal

C                   For now, this still applies to fields rather then particles;
C                       hence there is no distinction between particle & anti-p.

C                 To distinquish between them, either write another module, or
C                       use the convention, IT1 = -1, -2 for the anti-part;
C                            and test the signs upon Entry to determine the 
C                            channel (i.e. scattering, decay, or production).

      IF ((IT1.LT.1 .OR. IT1.GT.NSP) .OR. (IT2.LT.1 .OR. IT2.GT.NSP)) 
     >  THEN
        WRITE (NOUT,*)' Lepton Label out of range in GEWLT; IT1,IT2 =',
     >  IT1, IT2
        STOP
      ENDIF

      IF     (IBT .EQ. 1) THEN
         TEM = GLV (IT1, IBS, IT2)
      ELSEIF (IBT .EQ.-1) THEN
         TEM = GLA (IT1, IBS, IT2)
      ELSE
         WRITE (NOUT, *) 
     > ' IBoson Tensor index out of range in GEWLT; IBT = ', IBT
         STOP
      ENDIF

      GEWLT = TEM

      Return
C                        ****************************

      Entry GEWLH (IT1, IBS, IBH, IT2)

C           (G) Electro-Weak coupling for Leptons in the cHirality base.
C            -  -       -                 -               -

C               To be used in the calculation of general E-WK matrix elements.

C           IT1, IT2 = 1,2 :  Weak Isospin (T3) of the two leptons
C           IBS =  1 -4    :  Boson label (see SETEWK)
C           IBH =             Boson polarization (helicity) label
C                  1       :  right-handed
C                 -1       :  left -handed
C                  other   :  illegal

C                 For now, this still just apply to fiels rather then particles;
C                       hence there is no distinction between particle & anti-p.

C                 To distinquish between them, either write another module, or
C                       use the convention, IT1 = -1, -2 for the anti-part;
C                            and test the signs upon Entry to determine the 
C                            channel (i.e. scattering, decay, or production).

C         When make the field ---> particle transition, cHirality ---> Helicity

      IF ((IT1.LT.1 .OR. IT1.GT.NSP) .OR. (IT2.LT.1 .OR. IT2.GT.NSP)) 
     >  THEN
        WRITE (NOUT,*)' Lepton Label out of range in GEWLH; IT1,IT2 =',
     >  IT1, IT2
        STOP
      ENDIF

      IF     (IBH .EQ. 1) THEN
         TEM = GLR (IT1, IBS, IT2)
      ELSEIF (IBH .EQ.-1) THEN
         TEM = GLL (IT1, IBS, IT2)
      ELSE
         WRITE (NOUT, *) 
     > ' IBoson cHirality index out of range in GEWLH; IBH = ', IBH
         STOP
      ENDIF

      GEWLH = TEM

      Return
C                        ****************************

      Entry GEWQT (IQ1, IBS, IBT, IQ2)

C           (G) Electro-Weak coupling for Quarks in the Tensor base.
C            -  -       -                 -             -

C               To be used in the calculation of general E-WK matrix elements.

C           IBT =  1     :  Vector
C                 -1     :  Axial-vector
C                  other :  illegal

      IF ((IQ1.LT.1 .OR. IQ1.GT.NFL) .OR. (IQ2.LT.1 .OR. IQ2.GT.NFL)) 
     >  THEN
        WRITE (NOUT,*)' Quark Label out of range in GEWQT; IQ1,IQ2 =',
     >  IQ1, IQ2
        STOP
      ENDIF

      IF     (IBT .EQ. 1) THEN
         TEM = HQV (IQ1, IBS, IQ2)
      ELSEIF (IBT .EQ.-1) THEN
         TEM = HQA (IQ1, IBS, IQ2)
      ELSE
         WRITE (NOUT, *) 
     > ' IBoson Tensor index out of range in GEWQT; IBT = ', IBT
         STOP
      ENDIF

      GEWQT = TEM

      Return
C                        ****************************

      Entry GEWQH (IQ1, IBS, IBH, IQ2)

C           (G) Electro-Weak coupling for Quarks in the cHirality base.
C            -  -       -                 -             -

C               To be used in the calculation of general E-WK matrix elements.

C           IBH =  1     :  right-handed
C                 -1     :  left -handed
C                  other :  illegal

      IF ((IQ1.LT.1 .OR. IQ1.GT.NFL) .OR. (IQ2.LT.1 .OR. IQ2.GT.NFL)) 
     >  THEN
        WRITE (NOUT,*)' Quark Label out of range in GEWQH; IQ1,IQ2 =',
     >  IQ1, IQ2
        STOP
      ENDIF

      IF     (IBH .EQ. 1) THEN
         TEM = HQR (IQ1, IBS, IQ2)
      ELSEIF (IBH .EQ.-1) THEN
         TEM = HQL (IQ1, IBS, IQ2)
      ELSE
         WRITE (NOUT, *) 
     > ' IBoson cHirality index out of range in GEWQH; IBH = ', IBH
         STOP
      ENDIF

      GEWQH = TEM

      Return
C                        ****************************
      END

