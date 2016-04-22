      SUBROUTINE QCDOUT(NOUT)
C                                                   -=-=- qcdout
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Prints out the values of parameters to unit NOUT.
C
      logical AlphaSet, AmHatSet
      COMMON /QCDPAR_CT/ FSWITCH, AL, NF, NORDER, ALPHASET
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet

      Double precision POLEMASS, MSBMASS
      Integer iMass
      Logical MassesSet
      COMMON / COMQMS/ POLEMASS(9), MSBMASS(9), iMass, MassesSet

C
      WRITE (NOUT,110) AL, NF, NORDER,FSWITCH, iMass
 110   FORMAT(
     1 ' Lambda (MSBAR) =',G13.5,', NFL (total # of Flavors) =',I3,
     2 ', Order (loops) =', I2,', Fswitch = ',G13.5,', iMass =',I1)
cmarco here we inserted POLEMASS
       WRITE (NOUT,120) (I,POLEMASS(I),I=1,NF)
 120   FORMAT (3(' poleM', I1, '=', G13.5, :, ','))
cmarco here we inserted MSbar MASS
       WRITE (NOUT,140) (I,MSBMASS(I),I=1,NF)
 140   FORMAT (3(' MSbM', I1, '=', G13.5, :, ','))

      IF (NHQ .GT. 0)
     1   WRITE (NOUT,130) (I, ALAMF(I), I = NF-NHQ, NF)
 130   FORMAT (: ' ! Effective lambda given number of light quarks:'/
     >    (2(' ! ', I1, ' quarks => lambda = ', G13.5 : '; ')) )
      RETURN
      END
C
C***************************************************************
C
