      Function xSecLh (iLptn1, iLptn2, iHadn, x, Q, y, iFlv, iXsc)

! 2006.12.06:
             ! Common block to store heavy flavor S.F.s deleted.
             ! Logic of the calculation of heavy flavor xSec drastically simplified by HLL.

! 2005.12.25: original version.

C     ----------------------------
C     Reduced cross section dSig / dxdQ^2 for general lepton-hadron scattering.
C
C         First written 2006.02.09 by wkt: a part of the DisPac05;
C         revised by HLL for the reduced cross section, 2006.03;
C         Re-written again by wkt for general cross sections, according to
C         the tech note, DisPac05a.tex, 2006.03.20.
C
C         Formulas straight from the ACOT-I paper, in particular Eq.(14).
C         Use Helicity Structure Functions for definiteness (no ambiguity on defs)
C         and for their direct association with the group-theoretical kinematics.
C     ----------------------------

C                            /  iLptn2
C              iLptn1 ______/
C                           \  iBsn / iCmptn
C                            \ ______  iFlv
C                             |------  X
C                            /
C                     iHadn /
C
C     The process will be specified by:

C     iLptn1,2 : external lepton : (1) nu,      (2) e-/mu-,
C                                , (-1) nubar, (-2) e+/mu+

C     iFlv :   0   : total inclusive
C            -6:6  : semi-inclusive of flavor iFlv
C           11:16  : semi-inclusive (quark + anti-quark) of flavor 1:6

C     iXsc :   0   : DESY Reduced Cross sections
C              1   : My Effective xSec. defined in the DisPac05 tech note
C              2   : d Sig / dx dy
C              3   : d Sig / dx dQ^2
C              4   : .... add as you please!

      Use mDisPac
      Implicit None
                                          ! These are physics variables
      integer, intent(in) :: iLptn1, iLptn2, iHadn, iFlv, iXsc
      Real*8, intent(in)  :: X, y, Q
      Real*8 :: xSecLh
                                          ! These are input theory parameters
      integer ::  iOrdr, iSchm
      Real*8  ::  Q0, ScalFc, achi
      Common / cDisCal / Q0, ScalFc, achi, iOrdr, iSchm

      Integer :: nPol(-2:2) = (/2, 1, 0, 1, 2/)  ! only needed in the absolute xSec section
     >         , nSgn(-2:2) = (/-1,-1,0, 1, 1/)

      Real*8  :: Q2, chPsi, chPsi2, shPsi2, zeta, wMas2
      Real*8  :: rXsec, EfXsc

      Real*8, PARAMETER:: hMass = 0.938, Ssml=1d-20

      Real*8 :: AlfEwk, vBnMas             !, EffXsc

      Real*8  :: Sfn123(3), SfnHel(-1:1)


C                        -----------------------------
      If (.Not. DisSet) Call mSetDis
                                          ! Kinematical variables
      Q2 = Q * Q

      chPsi = (2d0-y) /y /Sqrt(1d0 + (2d0*x*hMass/Q)**2)
      chPsi2 = chPsi**2
      shPsi2 = chPsi2 - 1d0

      Current: If(iLptn1 == iLptn2) then                       ! neutral current case

        EfXsc = EffXsc(1)                   ! This is the gamma-gamma term

        zeta = AlfEwk(4) / AlfEwk(1)          ! ratio of gauge couplings

        zeta = zeta * Q2 / (Q2 + vBnMas(4)**2)          ! ratio of propagators

cdump
c      print*,' In xSecLh, Q, zeta =', sqrt(Q2), zeta
c      print*,' EfXsc, EffXsc(5),EffXsc(4) =', EfXsc,EffXsc(5),EffXsc(4) 
      
CCPY 03/12/2012 A missing factor 2 in front of the photon-Z interference
c        is restored
C        If (zeta >= 5d-3) Then                ! add the gamma-Z term, if needed:
C           EfXsc = EfXsc + zeta * EffXsc(5)
C                                              ! add the ZZ term, if needed:
C           If (zeta >= 3d-2) EfXsc = EfXsc + zeta**2 * EffXsc(4)
C        End If
        EfXsc = EfXsc + 2d0*zeta * EffXsc(5) + zeta**2 * EffXsc(4)

c htj
c        print*,"==========1=========="
c        print*,efxsc, EffXsc(5), EffXsc(4), EffXsc(1)
c        stop

        rXsec = EfXsc * 4d0 * x / (1d0 + chPsi2)       ! DESY Reduced xSec.

        If (iXsc > 0)                            ! My Effective xSce. needed only
     >  EfXsc = EfXsc * (4d0*Pi *AlfEwk(1))**2   ! for iXsc > 0 calculations later

c htj
c        print*,EfXsc,AlfEwk(1)
c        stop

      Else If (Abs(iLptn1-iLptn2) ==  1) Then   Current     ! Charged Current
        If (iLptn1-iLptn2 == -1) Then           ! W+ exchange
          EfXsc = EffXsc(2)

        Else If (iLptn1-iLptn2 ==  1)  Then     ! W- exchange
          EfXsc = EffXsc(3)
        End If
        rXsec = EfXsc * x * y**2 / 2d0 / nPol(iLptn1) ! DESY-type Reduced xSec.

        If (iXsc > 0) Then                           ! My Effective xSce. needed only
        wMas2 = vBnMas(2)**2                         ! for iXsc > 0 calculations later
        EfXsc = EfXsc * (wMas2 *Q2 *Gf /(Q2 +wMas2))**2 /2d0
        End If
      Else                                     Current
        rXsec = 0
        EfXsc = 0
        Print *, "Warning: iLptn1,2 =", iLptn1,iLptn2,
     >           " illegal in xSecLh call!"
      End If Current

      Select Case (iXsc)   ! Return Cross sections

      Case (0)                              ! DESY Reduced xSec
        xSecLh = rXsec
      Case (1)                              ! My Effective xSce.
        xSecLh = EfXsc
      Case (2)                              ! dSig/dx dy
        xSecLh = EfXsc * y /(2d0*Pi*Q2*nPol(iLptn1)) * XsGv2pb
      Case (3)                              ! dSig/dx dQ^2
        xSecLh = EfXsc * y**2 /(2d0*Pi*Q2*Q2*nPol(iLptn1)) * XsGv2pb
      Case (4)                              ! dSig/dE dcos
cpn12        xSecLh = EfXsc * tem * XsGv2pb
        stop 'temporary stop in lhXsec: tem is not defined'
cgino marco
      Case (5)                              ! H1 combined data
        xSecLh = rXsec
      Case Default
        Print *, "iXsc =", iXsc, "  illegal in xSecLh call!"
        Stop
      End Select
                                      ! Avoid returning zero answer
      If(xSecLh < Ssml) xSecLh = Ssml

      CONTAINS
                              ! --------------
        Function EffXsc(iCC)     ! Calculates the "Effective cross-section",
                                 ! defined in the DisPac05 Tech Note.
        Integer, intent(in) :: iCC          ! iCompton : 1 - 5
        Real*8 :: EffXsc

                                     ! Calculate the structure functions
        Call SfnGm(iHadn, iCC, iFlv, X, Q, SFn123, SfnHel)

                                ! Calculate the "reduced cross-section"
        EffXsc = g2Lm(iLptn1,iCC) *(SfnHel(-1) -SfnHel(1)) *chPsi
     >         + g2Lp(iLptn1,iCC) *
     >  (SfnHel(0) *shPsi2 +(SfnHel(1)+SfnHel(-1)) *(1d0 +chPsi2) /2d0)

c htj
c        print*,"===== icc ====="
c        print*,icc,Effxsc
c        print*,SfnHel(-1),SfnHel(1),SfnHel(0)
c        print*,g2Lm(iLptn1,iCC),g2Lp(iLptn1,iCC),chPsi,chPsi2
c        if (icc==5) stop
c        stop

        End Function EffXsc
                               ! ***************
      End Function xSecLh
