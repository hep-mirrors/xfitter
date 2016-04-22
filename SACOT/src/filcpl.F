      Subroutine FILCPL (CH, GV, GA, GL, GR)
C                                                   -=-=- filcpl

C            Given the charge of the fermion, compute the ElectroWeak couplings 
C            for a given generation.  Used in SetEwk.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NSP = 2, NGN = 3, NBN = 4, NPOL = 3)

      COMMON / EWKPAR / WMS, ZMS, SWG2, ALFE, ALFEW(NBN)

      DIMENSION CH(NSP), GV(NSP,NBN,NSP), GA(NSP,NBN,NSP), 
     >          T3(NSP), GL(NSP,NBN,NSP), GR(NSP,NBN,NSP)

      DATA T3 / 0.5, -0.5 /

      DO 31 IS1 = 1, NSP
      DO 31 IS2 = 1, NSP
        IF     (IS1 .EQ. IS2) THEN
            GV (IS1, 1, IS2) = CH(IS1)
            GA (IS1, 1, IS2) = 0.0
            GV (IS1, 2, IS2) = 0.0
            GA (IS1, 2, IS2) = 0.0
            GV (IS1, 3, IS2) = 0.0
            GA (IS1, 3, IS2) = 0.0
            GV (IS1, 4, IS2) = T3(IS1) - 2.* CH(IS1) * SWG2
            GA (IS1, 4, IS2) =-T3(IS1)
        ELSEIF (IS1 .GT. IS2) THEN
            GV (IS1, 1, IS2) = 0.0
            GA (IS1, 1, IS2) = 0.0
            GV (IS1, 2, IS2) = 1.0
            GA (IS1, 2, IS2) =-1.0
            GV (IS1, 3, IS2) = 0.0
            GA (IS1, 3, IS2) = 0.0
            GV (IS1, 4, IS2) = 0.0
            GA (IS1, 4, IS2) = 0.0
        ELSEIF (IS1 .LT. IS2) THEN
            GV (IS1, 1, IS2) = 0.0
            GA (IS1, 1, IS2) = 0.0
            GV (IS1, 2, IS2) = 0.0
            GA (IS1, 2, IS2) = 0.0
            GV (IS1, 3, IS2) = 1.0
            GA (IS1, 3, IS2) =-1.0
            GV (IS1, 4, IS2) = 0.0
            GA (IS1, 4, IS2) = 0.0
        ENDIF

        DO 32 IBN = 1, NBN
            GL (IS1,IBN,IS2) = (GV(IS1,IBN,IS2) - GA(IS1,IBN,IS2)) / 2.
            GR (IS1,IBN,IS2) = (GV(IS1,IBN,IS2) + GA(IS1,IBN,IS2)) / 2.
   32   CONTINUE

   31 CONTINUE


      Return
C                        ****************************
      END

