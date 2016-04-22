      SUBROUTINE ALFSET (QSin, AlfSin)

C   Given the value of Alpha(strong), ALFS, at the scale QS, 
C   this routine determines the QCD running coupling at any other scale.
C   If the input Qs=0, AlfSet will use the boundary condition that was
c   stored from the previous call of AlfSet.  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      double precision Qs, AlfS, alam1, blam1, err, AlfSin
      EXTERNAL RTALF
      COMMON / RTALFC / ALFST, JORD, NEFF

      integer iAlsEvl, iPDFEvl 
      Common / Evolpak / iAlsEvl, iPDFEvl

      logical AlphaSet, AmHatSet
      COMMON /QCDPAR_CT/ fswitch, AL, NF, NORDER, ALPHASET
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ, AmHatSet

      data qs, alfs, alam1, blam1,err /2*99d0, 0.01d0, 10d0, 1d-3/
      save Qs, AlfS

      if (.not. AmHatSet) call SetAmhatf !set switching points

cpn11   Initialize alpha_s using the boundary condition AlfSIn(QsIn) 
c       passed as the arguments of AlfSet. Use the initial condition
c       Alfs(Qs) from the previous call of AlfSet if QsIn=0
      if (abs(QsIn).gt.1d-7) then
        Qs=QsIn; Alfs=AlfsIn
      endif 

c htj, comment out for xfitter
c      If (iAlsEvl.eq.0) then !our own alpha_s -- truncated
c        ALFST = ALFS
c        CALL PARQCD (2, 'ORDR', ORDR, IR1)
c        JORD  = nint(ORDR)
c        
c        NEFF = NFL(QS)
c        
c        EFLLN  = QZBRNT (RTALF, ALAM1, BLAM1, ERR, IR2)
c        EFFLAM = QS / EXP (EFLLN)
c        
c        CALL SETL1 (NEFF, EFFLAM)
c htj, comment out for xfitter
c      elseif(iAlsEvl.eq.1) then ! QCDNUM
cpn11 temporary stop
c      stop 'Fix switching points and NFL before calling alfset'

c        QS2=QS*QS
c        Call setalf(ALFS,QS2)
c htj, comment out for xfitter
c      ElseIf (iAlsEvl.eq.2) then ! Hoppet
c        Call SetAlf_hoppet(Qs, Alfs)         
c      else
c        print *, 'Alfset: iAlsEvl=',iAlsEvl,' is not supported'
c        stop
c      Endif                     !iAlsEvl
 
      ALPHASET =.TRUE.

      END
 
C
C**************************************************************
C
