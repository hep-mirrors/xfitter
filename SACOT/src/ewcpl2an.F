      Function EwCpl2An (JP1, JBN, JP2)
C                                                   -=-=- ewcplg2

C     Function call of the squared EW coupling constants for Drell-Yan
C     processes in the Annihilation channel.

C     EwCplC common block is set up by the SetEwCpl2 subroutine

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (NFL = 6, NBN = 4)

      COMMON / EWCPLC / CPLANH(-NFL:NFL, NBN, -NFL:NFL), 
     >                  CPLSCT(-NFL:NFL, NBN)

      EwCpl2An = CPLANH (JP1, JBN, JP2)

      Return
C                        ****************************

      Entry EwCpl2Cn (JP, JBN)

C     Function call of the squared EW coupling constants for Drell-Yan
C     processes in the Compton Sc. channel.

      EwCpl2Cn = CPLSCT (JP, JBN)

      Return
C                        ****************************
      END

