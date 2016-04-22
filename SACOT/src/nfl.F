      FUNCTION NFL(AMU)
C                                                   -=-=- nfl
      implicit none
C           NFL returns the number of active quark flavors.
cpn10 Revised and simplified to be used independently from CTEQ
c     alpha_s evolution
      double precision alam, amhat, amn, fswitch, AL, amu
      integer nhq, NF, NORDER, nfl, i
      LOGICAL AlphaSet, AmHatSet
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet
      COMMON / QCDPAR_CT / fswitch, AL, NF, NORDER, ALPHASET
C
      if (.not. AmHatSet) call SetAmhatf !set switching points

      NFL=NF-NHQ
cpn10 either no heavy quarks, or amu in the non-pert. region
      if (NHQ.eq.0.or.(AMU.le.AMN)) return
      
      do i=nf-nhq+1, nf
        if (amu.ge.amhat(i)) then 
          nfl=i
        else 
          return
        endif                   !amu
      enddo                     !i=nf-nhq+1
      
      return 
      END
C
C***************************************************************
C
