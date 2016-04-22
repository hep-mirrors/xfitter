      subroutine SetAmhatf
C Sets switching points Amhatf between Nf and Nf+1 flavors, given the
c input quark masses QMS and rescaling factor Fswitch. The type of
c the input quark mass (on-shell or MSbar) is determined 
c by the flag iMass inside fitinp
      implicit none
      double precision fswitch, al, alam, amhat, amn
      integer nf, norder, nhq, n, iwrn
      logical set, AmHatSet

      COMMON / QCDPAR_CT / FSWITCH, AL, NF, NORDER, Set
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
      Double precision POLEMASS, MSBMASS, qms(9)
      Integer iMass
      Logical MassesSet
      COMMON / COMQMS/ POLEMASS(9), MSBMASS(9), iMass, MassesSet

      integer nin, nout, nwrt
      COMMON / IOUNIT / NIN, NOUT, NWRT

      AMHAT(0) = 0.; NHQ=0
cpn10 Do a preliminary separation between the light and heavy quarks, 
c    by assuming that quarks with  masses lighter than the
c    hard-wired scale AMN=0.5 GeV (in DatQCD) are automatically light, 
c    and the other quarks are heavy. If Lambda_QCD is calculated
c    explicitly in SetL1, the distinction will be refined in SetL1

      DO N = 1, NF
cpn11   Decide if the input mass is an on-shell or MSbar mass
       IF(iMass.EQ.0) THEN
          QMS(N)=POLEMASS(N)
        ELSE
          QMS(N)=MSBMASS(N)
        ENDIF

        if (QMS(N).le.Amn) then 
cpn10 Quarks with on-shell masses below Amn are automatically treated
c     as massless 
          AMHAT(N)=0d0
        else 
cpn10 Quarks with on-shell masses above Amn are treated as massless
c     only if the switching point is below Amn 
          AMHAT(N) = EXP(fswitch*dlog(10d0)/2d0)*QMS(N)
          if (AMHAT(N) .le. Amn) then 
            CALL WARNI (IWRN, NWRT, 
     >        'Nonpert. AMHAT(N) in SetAmhatf; setting AMHAT(N)=0',
     >        'N =', N, 0,9, 0)
            AMHAT(N)=0d0
          else                  !increment the number of heavy quarks
            NHQ=NHQ+1
          endif                 !amhat(N).le. Amn
        endif                   !qms(N).le. Amn
      enddo                     !n=1,nf

      AmHatSet=.true.
      return 
      End                       !SetAmhatf
