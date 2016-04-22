      Function EwCplSc (Iboson, Ihelic, Ipartn)
C                                                   -=-=- ewcplg

C                Electro-weak coupling for VecBos + Parton --> Parton2

C                Use helicity basis for the boson polarization index.

C                       Iboson (1:4) |  Ihelic (-1:1)
C                                      |
C                                       |______________ 
C                                       |         (Parton2, redundant)
C                                      |
C                      Ipartn (-6:6)  |

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      COMMON / WKCPLG / QCPLN(4, -1:1, -6:6)

      SAVE

      If (  (Iboson.Lt. 1 .or. Iboson.Gt.4) 
     > .or. (Ihelic.Lt.-1 .or. Ihelic.Gt.1) 
     > .or. (Ipartn.Lt.-6 .or. Ipartn.Gt.6) ) Then
         Print '(/ A / A, 3I4 /)'
     >,'Iboson, Ihelicity, Iparton out of range in EwCplSc; '
     >,'incorrect input values : ', Iboson, Ihelic, Ipartn
         Stop
      EndIf

      tem = QCPLN (Iboson, Ihelic, Ipartn)

      EwCplSc = tem

      RETURN
C                        ****************************

      Entry EwCplAn (Iboson, Ihelic, Ipartn)

C           Electro-weak coupling for Parton + Parton2bar ---> VecBoson

C                      Ipartn (-6:6)  | 
C                                      |
C                                       |______________Iboson (1:4) 
C                                       |              Ihelic (-1:1)
C                                      |
C                         Parton2bar  |

      Print '(/ A / A /)', 'EwCplAn has not yet been implemented.'
     > , 'Can be obtained by crossing from EwCplSc above.'

      EwCplAn = 0.

      END

C                                                          =-=-= EwCpl2
