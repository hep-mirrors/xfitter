      Subroutine SetQCD
C                                                   -=-=- setqcd
C These comments are included in the lead subprogram to survive forsplit.

C===========================================================================
C GroupName: Setqcd
C Description: Set up the qcdpac of programs, initiate common blocks
C ListOfFiles: setqcd
C===========================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Setqcd.f,v 1.1 97/12/21 20:35:00 wkt Exp $ 
C #Log:	Setqcd.f,v $
c Revision 1.1  97/12/21  20:35:00  wkt
c Initial revision
c 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      External DatQCD
      LOGICAL ALPHASET                      !pn10 explicit common seems to be
      COMMON / QCDPAR_CT / FSWITCH, AL, NF, NORDER, ALPHASET !asked for by the latest g95

      ALPHASET = .false.

      Return
C                        ****************************
      END

      BLOCK DATA DATQCD
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
c#include "couplings.inc"
      LOGICAL ALPHASET, AmHatSet
      DOUBLE PRECISION POLEMASS, MSBMASS
      INTEGER IMASS
      LOGICAL MassesSet
      COMMON / COMQMS/ POLEMASS(9), MSBMASS(9), IMass, MassesSet

      COMMON /COMQCH/ VALQCH(9) !fractional quark electric charges
      COMMON /QCDPAR_CT/ fswitch, AL, NF, NORDER, ALPHASET
      COMMON /COMALP/ ALPHA
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
C
      DATA Fswitch, AL, NF, NORDER, ALPHASET /0d0, .226, 5, 2, .FALSE. /
      DATA VALQCH/ 0.66666667, -0.33333333,
     >  -0.33333333, 0.66666667,
     >  -0.33333333, 0.66666667,
     >  3*0./

cpn10 To choose the correct maximal number of active flavors, NfMx, 
c all switching points AMHAT are very large initially, 
c until they are explicitly set to be in the physical range in SetAmhatf  
      data AMHAT, Amn, AmHatSet /10*1e+15,0.5d0, .false./

cmarco iMass specifies the type of quark masses provided by fitinp. 
c      iMass=0 or 1 for pole (=on-shell) or MSbar input quark masses.
c
c      The logical flag MassesSet specifies whether the values of the 
c      pole and and MSbar masses in COMMON/COMQMS are compatible. If
c      they are compatible, MassesSet = .true. However,
c      when the value of the pole or MSbar masses is changed 
c      by calling QCDSET, all quark masses need to be re-synchronized,
c      and the evolution of alpha_s and running masses must be performed
c      using the new mass values, before any of the new mass values
c      can be safely retrieved back. Technically, this is achieved
c      by setting MassesSet=.false. every time when the masses in 
c      COMMON/COMQMS are changed in QCDSET. After that,
c      the pole masses, quark masses, and alpha_s must be resynchronized
c      by calling the subroutine MassSync, where MassesSet is set to .true.. 

      data iMass, MassesSet /1, .FALSE./
cmarco Default pole and MSbar masses
      DATA POLEMASS/  2*0.001, 0.2, 1.3, 4.5, 174., 3*0./
      DATA MSBMASS/  0.00249, 0.00505, 0.101, 1.27, 4.19, 172., 3*0./

C                       ******************************
      END

C                                                          =-=-= Alphas
