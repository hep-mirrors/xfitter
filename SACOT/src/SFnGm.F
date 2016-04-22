! ----------------------------------------------------------------------------------
      Subroutine SfnGm(iHadn,iCmptn,iFlv, X,Q, Sfn123,SfnHel)
! ----------------------------------------------------------------------------------
Cmarco Returns arrays of structure functions Sfn123, SfnHel
! 2008.10.02:
             ! ZM calculations (mSwtch = 0, 1) re-written.
             !   mSwtch = 0 should now correspond to the "conventional ZM scheme":
             !              same active partons, determined by initial state fac. scheme only;
             !              ignore final-state phase space considerations all together.
             !          = 1,8 corresponds to a possible version of "improved ZM scheme":
             !              respects final-state phase-space restrictions, including using Acot-chi;
             !              select partons according to gluon-fusion final-state, add corresponding
             !              LO term by appropriate adjustment of scale (to be above the threshold).
! 2006.12.06:
             ! New input argument iFlv introduced; common block to store heavy flavor S.F.s deleted.
             ! Logic of the calculation of heavy flavor S.F.s drastically simplified by HLL.
             ! New option for choice of factorization scale (iFscle=3) for HQ production added.
! 2006.03.13:
             ! Transformed the HvFlvSf output to general "semi-inclusive" flavor SF's in /SiFlvSf/
! 2006.03.10:
             ! Expanded to output heavy quark structure functions in common / HvFlvSf /
! 2006.01.24:
!             Expanded to include iCmptn = -1/0 for (w+ -/+ W-)/2

                    ! ===== General Mass Hadron Structure Functions
! 2005.12.29:
             ! -------This is the main program module to calculation S.F.'s--------!

!    Forward Compton Sc.:

!         iBsn(1) + iHadn --> iBsn(2) + iHadn  , with iFlv + X (summed) in the intermediate state.

!                   iBsn1 ________ __ __iFlv__ __ _____iBsn1
!                                 |\/|        |\/|
!                   iHadn ________|/\|________|/\|_____iHadn
!                                        X

!                         iCmptn =   1,     2,     3,     4,        5
!                   for iBsn(1,2)= (1,1), (2,2), (3,3), (4,4), (1,4)+(4,1)
!                                  photon   W+     W-      Z   photon-Z int.

!         iFlv : 0   (total inclusive---sum over all active quark-parton flavors);
!               -6:6 (individual flavors);
!               11:16 (quark+antiquark).

                          ! -----------------------------------!

!         iOrdr: 0  :   LO : alpha_s^0
!                1  :  NLO : alpha_s^1
!                2  : NNLO : alpha_s^2

!         iSchm: 4-digit integer tkmn
!                n : factorization number scheme switch:
!                  0   : VFNS : # of flavors depends on the factorization scale;
!                  1-6 : FFNS : fixed flavor number = n;
!                m : mass-effect switch:  (Cf. ACOT-chi: Kretzer/Schmidt/Tung)
!                  0   : ZM (zero-mass) partons (kinematic limits neglected)
!                  1,8  : ZM - kinematic limits respected (klr.), for two possible
!                   definitions of "unsubtracted" one-loop gluon matrix element
!                   below the heavy-quark production threshold
!                  2   : GM1- klr. + Massive Wilson Coef. up to alphaS^1
!                             + ZM of alphaS^2 if required (SACOT + ACOT-chi)
!                  3   : GM2- klr. + Massive Wilson Coef. up to alphaS^2
!                             (SACOT + ACOT-chi)
!                k : factorization scale switch:
!                  0 : default : factorization scale aMu = Q ;
!                  1 : aMu = (Q**5 + Q0**5)**(1d0/5d0)  ;
!                  2 : aMu = ScalFc * Q    where ScalFc is supplied via Common/cDisCal/
!                  3 : aMu = sqrt (Q^2 + M1^2 + M2^2) for heavy flavor production.
!                  other options, not yet implemented;
!                t : switch for specific test/debug function: e.g.
!                  0 : default
!                  1 : print out comparison between SfZM05 and MVV-nlo
                         ! ----------------------------------------- !
!    Outputs F_(1,2,3) as the 3-dim array        Sfn123( 1/2/3),
!    and the helicity SFs F(Left,Long.,Right) as SfnHel(-1/0/1).

      Use mDisPac
!      Use GMalf2
!      Use ZMalf2fcharged


      Implicit None
!                                         Dummy (physical) variables
      Integer, Intent(in) :: iHadn, iCmptn, iFlv
      Real*8,  Intent(in) :: X, Q
      Real*8, Intent(Out) :: Sfn123(3), SfnHel(-1:1)

!           Global variables that depend on theoretical scheme used;
!     Should be set once and for all in the front-end (or interface) program (chidis.f)
      Integer :: iOrdr, iSchm
      Real*8  :: Q0, ScalFc, achi     ! Q0 is the Qmin below which PDFs are not defined;
                                ! With iFscle = 1, aMu is defined to be always > Q0.
                                ! ScalFc is the (dimensionless) factorization scale;
                                ! With iFscle = 2,3, aMu ~ ScalFc * the physical scale. 
                                !achi is the shape parameter in the rescaling variable chi
      Common / cDisCal / Q0, ScalFc, achi, iOrdr, iSchm

!                                                Local variables
      Logical :: CC, SubAdd, Qpdf(-6:6,-6:6)
      integer, parameter :: nPrtn2=6
      Integer :: iFschm =0, iFscle =0, mSwtch =0,     
     >  nPrtn1 =3, iFv, kFv,           
     >  iPrtn1, iPrtn2, aPrtn1, aPrtn2, NFAboveThr=0        
      Real*8  :: aMu, Zeta, aM1, aM2, W 
      Real*8  :: Q0a(5),Q0b(5),Q1a(5),Q1b(5),Qtem(5),Gtem(5), RG(5)
      Real*8  :: SfT2(3), SfH2(-1:1), SfT0(3), SfH0(-1:1), SfTs(3),  
     >  SfTq(3), SfHq(-1:1), SfTg(3), SfHg(-1:1), SfHs(-1:1),  
     >  SfT3(3), SfH3(-1:1)

      logical :: AboveThreshold(nPrtn2, nPrtn2)

      Integer :: nHeavy = 4     ! iParton >= Nheavy are "heavy" (wrt to Lambda_QCD)
      Real*8  :: Smll = 1d-8, sml = 1d-5 ! small numbers cut-offs
      Integer ::iprtn3, NFL , i !htj
      Real*8  :: Alpi, aMass, fchi, bMu1, bMu2, FLo3, F2o3,
     >  F2tmp3, FLtmp3, e2, SiF2l, SiFLl
      integer iw,  NIN, NOUT, NWRT
      data iw /0/
      COMMON / IOUNIT / NIN, NOUT, NWRT

cpn11 SemiInclusiveF2c=.TRUE. (default) 
c                       to return semi-inclusive F2c for heavy quarks
c     SemiInclusiveF2c=.FALSE. to return the heavy-quark part F_h 
c     of the inclusive F2
      logical :: SemiInclusiveF2c=.TRUE. !DO NOT CHANGE MINDLESSLY     

!     ------------------------------------- Initialization

      If (.Not. DisSet) Call mSetDis
      Qpdf = .false.            ! tag quark contribution channels (for subtraction purpose)

      Select Case (iCmptn)
      Case (-1,0,2,3)           ! Charged Current process
        CC = .true.
      Case Default              ! Neutral Current process
        CC = .false.
      End Select
  
      iFschm = Mod (iSchm, 10)  ! flavor # factorization scheme
      mSwtch = Mod (iSchm/10, 10) ! mass treatment switch
      iFscle = Mod (iSchm/100, 10) ! factorization scale switch

c htj
c      print*,iFschm,mSwtch,iFscle  
c      stop
c 2/23/09 Lai disable below warning for the moment  
c      !pn2009 Improved ZM makes sense only at NLO 
c      WrongMSwtch: If (iOrdr > 1 .and. 
c     >  (mSwtch /=0 .or. mSwtch /= 2 .or. mSwtch /=3)) then
c        print *, 
c     >    'STOP in SFnGm: improved ZM scheme is only implemented at NLO'
c        print *, 'iOrdr, mSwtch = ', iOrdr, mSwtch
c      End If WrongMSwtch


      kFv = iFlv / 10             ! kFv = 1 for iFlv 14:16 ; 0 otherwise
      iFv = Mod(iFlv,10)        ! iFv -6:6
      
                          ! --------- Set factorization scale for the process
      FacScle: select case (iFscle)
      case (0,3)
        aMu = Q
        If (iFscle==3 .and. abs(iFv) >= nHeavy) ! Redefine the fac scale for semi-inclusive HQ prod.
     >    aMu = Sqrt(Q**2 + ScalFc * aMass(1, abs(iFv))**2) !use the MSbar mass in 
                                             ! the fact. scale for no deep reason     
      case (1)
        aMu = (Q**5 + Q0**5)**(1d0/5d0) ! aMu is bounded below by Q0, for any Q
      case (2)
        aMu = ScalFc * Q        ! Simple example; could be any function of (x,Q), i.e. Sqrt(Q^2+M_c^2)
      case default
        print *, 'Stop: iFscle = ', iFscle, 
     >    'has not been implemented in SFnGm.f!'
        stop 
      end select FacScle

c htj
c      print*,iflv,amu
c      stop
      ! Define number of incoming quark-partons nPrtn1   ---  FFNS and VFNS switch
                                         ! This is done with the FFNS label iSchm
      If (iFschm >= 1 .and. iFschm <= 6) Then
        nPrtn1 = iFschm
      ElseIf (iFschm == 0) then
        nPrtn1 = NFL(aMu)       !  VFNS: number of partons at the scale Q
      EndIf

!      nPrtn2 = NflTot()    ! number of final state quark-flavors to be summed
                           ! (default = 6)

      W = Q *SQRT(1d0/x -1d0)   ! Used to determine available final state phase space
                               ! Physically, this is really SQRT(W^2-M_N^2)
      call ExamineMassThresholds() 

c htj
c      print*,iprtn1,iprtn2
c      stop
                                             ! -----------Start calculation
              ! Part I : initial state QUARK-parton
      SfT0 = 0
      SfH0 = 0
      SfTq = 0
      SfHq = 0
      Part1: DO iPrtn1 = -nPrtn1, nPrtn1 ! Sum over initial state PARTONs only
      If (iPrtn1 == 0) Cycle    ! skip gluon (done separately later)
     
      aPrtn1 = abs(iPrtn1)
      

                                ! Store default ZM (zero-mass) results
      Call SfZm0Q (iHadn, Iprtn1,  X, Q, aMu, Q0b) ! ZM default
      If (iOrdr > 0) 
     >  Call SfZm1Q (iHadn, Iprtn1,  X, Q, aMu, Q1b) ! ZM default

c htj
c      print*,iHadn, Iprtn1,  X, Q, aMu, Q1b
c      stop

      Part2: Do iPrtn2 = -nPrtn2, nPrtn2 ! Sum over final state FLAVORs (all)
      If (iPrtn2 == 0) Cycle    ! skip gluon
        
      If(iFlv /= 0) then        ! semi-inclusive case : skip channels not needed.
        If (kFv==0 .and. iPrtn2 /= iFv) cycle ! iFv is -6:6 in this case
        If (kFv==1 .and. abs(iPrtn2) /= iFv) cycle ! iFv is  4:6 in this case
      Endif
        
      If (g2q(0,iPrtn1,iCmptn,iprtn2) < Smll) Cycle ! skip uncoupled channels
      
      Qpdf(iPrtn1, iprtn2) = .true. ! tag channels with the LO quark included
      aPrtn2 = abs(iPrtn2)
      aM1 = aMass(0, aPrtn1)  !use pole masses in all DIS kinematics
      aM2 = aMass(0, aPrtn2)    
      if (mSwtch.eq.1) then  ! set a special factorization 
                             !scale for some realizations of improved ZM
        bMu1 = Max(aMu, aM1)
      else
        bMu1 = amu
      endif

                                 ! LO  --------------
      Q0a = 0e0                 !LO quark contribution = 0, unless explicitly computed
  ! Treat the naive ZM scheme differently from other schemes: neglect threshold effects, do not
  ! distinguish between the summation over initial- and final-state flavors, and use massless 
  ! expressions for any x, Q 
      NaiveZM: If (mSwtch == 0) then 
        if (abs(iprtn2) <= nPrtn1) Q0a = Q0b
      Else NaiveZM              !implement threshold considerations       
c htj 
c      if (iordr==2) then
c         print*,aPrtn1,aPrtn2,AboveThreshold(aPrtn1,aPrtn2),nHeavy
c         print*,Q0b
c         stop
c      end if

        If (AboveThreshold(aPrtn1,aPrtn2)) Then ! Above production threshold
          Heavy: If (Max(aPrtn1,aPrtn2)<nHeavy) Then ! for both quarks light;
            Q0a = Q0b           ! Use ZM defaults
          Else
            Zeta = fchi(x,Q,am1,am2) !X ! Rescaling variable (generalized chi)
            If (CC .and.  (mSwtch==2.or.mSwtch==3)) Then !GM scheme
              Call Acot0Q(iHadn, iPrtn1, 0d0,aM2, Zeta,Q,aMu, Q0a)
            Else                !improved ZM; massless LO with zeta=chi variable
              Call SfZm0Q (iHadn, iPrtn1, Zeta, Q, bMu1, Q0a)
            End If              !CC
          End If Heavy
        End If !AboveThreshold
      End If NaiveZM
                                    ! --------- Assemble irreducible amplitudes,
         ! ---- fold in EW coupling to produce contributing term to Structure Functions;
                                              ! then, accumulate to Structure Functions
      Call SFnConv (Q0a, g2Q(1, iPrtn1, iCmptn, iprtn2), SfT2, SfH2)
      SfT0 = SfT0 + SfT2
      SfH0 = SfH0 + SfH2

c htj
c      print*,SfT2

        ! NLO Quark term ---------------
      Q1a = 0e0                 !NLO quark contribution = 0, unless explicitly computed

      NloQ: If (iOrdr >= 1) Then
        NaiveZM1: If (mSwtch == 0) then !Naive ZM; use massless expressions for 
                                !active quark flavors regardless of kinematics
          if (abs(iprtn2) <= nPrtn1) Q1a = Q1b    
        Else NaiveZM1           !implement threshold considerations           
          If (AboveThreshold(aPrtn1,aPrtn2)) Then 
            Heavy1: If (Max(aPrtn1,aPrtn2)<nHeavy) Then ! for both quarks light;
              Q1a = Q1b         ! Use ZM defaults
            Else Heavy1         ! at least one heavy quark
              Zeta =fchi(x,Q,am1,am2)!X ! Rescaling variable (generalized chi)
              If (CC .and. (mSwtch==2.or.mSwtch==3)) Then !GM scheme
                Call GotSK1Q (iHadn, iPrtn1, aM2, Zeta,Q,aMu, Q1a)
              Else              !improved ZM; massless quark NLO with zeta=chi variable
                Call SfZm1Q (iHadn, iPrtn1, Zeta, Q, bMu1, Q1a)
              End If
            End If Heavy1
          End If !AboveThreshold
        End If NaiveZM1

        Q1a = Q1a *Alpi(aMu)/2d0 !multiply the alpi here
        Call SFnConv (Q1a, g2Q(1, iPrtn1, iCmptn, iprtn2), SfT2, SfH2)
        SfTq = SfTq + SfT2
        SfHq = SfHq + SfH2

c htj
c        print*,Q1a
c        print*,alpi(amu),amu
c        print*,SfTq
c        stop

      End If NloQ
        
      End do Part2
      End Do Part1

c htj test
c      open(1001,file="alpi_test.xfitter.txt",status="unknown")
c      do i=1,185
c        amu=1.3d0+exp(dble(i)*0.05d0)-1d0
c        write(1001,'(3ES16.8)') amu,alpi(amu)
c        write(*,'(3ES16.8)') amu,alpi(amu)
c      end do
c      close(1001)
c      stop

               ! Part II : GLUON initial state  ---------------
      SfTg = 0                  ! gluon-fusion terms
      SfHg = 0
      SfTs = 0                  ! subtraction terms
      SfHs = 0
      Gluon: If (iOrdr >=1 ) then ! (gluon does not contribute to L0)
                                ! NC/CC processes are treated the same way; default
                                ! Naive ZM result is independent of quark flavors, store as RG
        Call SfZm1G (iHadn,  X, Q, aMu, RG)
        RG = RG * 2D0
        ! Sum over distinct final state pairs
        ! Calculated here is the t- and u-channel diagram contributions to
        !           vector-boson + gluon --> iPrtn1 + iPrtn2-bar
        ! (The virtual exchanged quark lines are iPrtn2/iPrtn1-bar respectively.)
        ParG1: Do iPrtn1 = 1,  nPrtn2 ! quark in the final state     :  iPrtn1
        ParG2: Do iPrtn2 = 1,  nPrtn2 ! anti-quark in the final state: -iPrtn2
           
        If(iFlv /= 0) then      ! semi-inclusive case : skip channels not needed.
          If (kFv==0 .and. iPrtn1 /= iFv .and. -iPrtn2 /= iFv) Cycle
          If (kFv==1 .and. iPrtn1 /= iFv .and. -iPrtn2 /=-iFv) Cycle
        Endif
           
        If (g2q(0,iPrtn2,iCmptn,iprtn1) < Smll) Cycle ! skip uncoupled channels
           
        aM1 = aMass(0, iPrtn1)  !use pole masses in all DIS kinematics
        aM2 = aMass(0, iPrtn2)

        if (mSwtch.eq.1) then ! set a special factorization scale 
                              !for some realizations of improved ZM
          bMu1 = Max(aMu, aM1); bMu2=Max(aMu,aM2)
        else
          bMu1 = amu; bMu2=amu
        endif

        Q0a = 0e0;  Q0b = 0e0   ! (Use arrays Q0a, Q0b, Q1a,Q1b to store  )
        Q1a = 0e0;  Q1b = 0e0   ! (subtraction terms, if needed.)
        Gtem = 0e0; Qtem=0e0    !NLO gluon and subtraction = 0, unless explicitly computed
        SubAdd = .false.        ! tag subtraction calculation (for gluon fusion processes)

        NaiveZM2: If (mSwtch == 0) then !If naive ZM, use massless NLO expressions for 
                                            !diagrams with active quarks regardless of kinematics
          if (iPrtn1 <= nPrtn1 .and. iPrtn2 <= nPrtn1)  Gtem = RG  
        Else NaiveZM2           !implement threshold considerations           
          If (AboveThreshold(iPrtn1,iPrtn2)) Then
            mSwG: Select Case (mSwtch)
          Case (1,8)  mSwG      ! improved ZM: honor kinematic limits (ACOT-chi)
                                !Evaluate the ZM gluon contribution (with the g->q qbar collinear
                                !term subtracted)
            Heavy21: If (Max(iPrtn1,iPrtn2) < nHeavy) Then ! For light flavors
              Gtem = RG
            Else Heavy21
              Zeta =fchi(x,Q,am1,am2) !X !fchi(x,Q,am1,am2) ! Rescaling variable (generalized chi)
              Call SfZm1G (iHadn,  Zeta, Q, bMu1, Q0a) 
              Call SfZm1G (iHadn,  Zeta, Q, bMu2, Q0b)   
              Gtem = Q0a+Q0b
                                !Restore the subtracted g->q qbar collinear term when the corresponding
                                !LO quark contribution does not exist (e.g., if mu < below the switching point mu_0
                                !for the incoming PDF flavor)
              If (.not.Qpdf(-iPrtn1,-iprtn2)) then
                SubAdd = .true.
                If (mSwtch==1) then
                  Call SfZm0Q (iHadn,-Iprtn1,Zeta, Q, bMu1, Q1a)  
                Else If (mSwtch==8) then
                  Call AcotSub(iHadn,-iPrtn1, aM1, 0d0, Zeta, Q, 
     >              bMu1, Q1a)
                End If
                Q1a(5) = -Q1a(5)
              EndIf             !.not.Qpdf(-iPrtn1,-iprtn2)
              
              If (.not.Qpdf(iPrtn2,iprtn1)) then
                SubAdd = .true.
                If (mSwtch==1) then
                  Call SfZm0Q (iHadn, Iprtn2, Zeta, Q, bMu2, Q1b)
                Else If (mSwtch==8) then
                  Call AcotSub(iHadn, iPrtn2, aM2, 0d0, Zeta, Q, 
     >              bMu2, Q1b)
                EndIf
              Endif             ! .not.Qpdf(iPrtn2,iprtn1)
            
              Qtem = -(Q1a+Q1b) ! -sign to ADD rather than SUBTRACT
            End If Heavy21
          Case (2,3) mSwG       ! general-mass scheme
          Heavy22: If (Max(iPrtn1,iPrtn2) < nHeavy) Then ! For light flavors
              Gtem = RG         ! Use default (ZM) Wilson coef.                 
            Else Heavy22        ! Include full mass effects in Wilson Coefficient
                                !   if either final quark is heavy
              Call Acot1G (iHadn, aM1, aM2, X, Q, aMu, Gtem)
                                ! Need mass subtraction to avoid double-counting if line-
                                ! reversed quark (-iPrtn1) is an initial-state parton
                                ! Quark contribution exists if Qpdf(-iPrtn1,-iprtn2)=.true.
              Zeta = fchi(x,Q,am1,am2)!X !fchi(x,Q,am1,am2) ! Rescaling variable (generalized chi)
              If (iPrtn1 <= nPrtn1 .and. Qpdf(-iPrtn1,-iprtn2)) Then
                SubAdd = .true.
                If (.not.CC) Then ! Neutral Current case:
                  Call AcotSub  ! Set mass=0, except for the log (m/mu) factor
     >              (iHadn, -iPrtn1, aM1, 0d0, Zeta,Q,aMu,Q1a)
                Else            ! Charge Current case
                  Call AcotSub  ! Keep final-state mass to match LO term.
     >              (iHadn, -iPrtn1, aM1, aM2, Zeta,Q,aMu,Q1a)
                End If          !iPrtn1 <= nPrtn1
                                ! in order to use the EW coupling below, need to reverse the
                                ! sign of the parity violating term for anti-quark -iPrtn1;
                                ! (CP symmetry ==> R<-->L)
                Q1a(5) = -Q1a(5)
              EndIf
                                ! Likewise for the other (crossed) diagram
                                !  which corresponds to initial state parton iPrtn2
                                ! Quark contribution exists if Qpdf(iPrtn2, iprtn1)=.true.
              If (iPrtn2 <= nPrtn1 .and. Qpdf(iPrtn2, iprtn1)) Then
                SubAdd = .true.
                If (.not.CC) Then ! Neutral Current case:
                  Call AcotSub 
     >              (iHadn, iPrtn2, aM2, 0d0, Zeta, Q, aMu, Q1b)
                Else            ! Charge Current case
                  Call AcotSub  ! Keep final-state mass to match LO term.
     >              (iHadn, iPrtn2, aM2, aM1, Zeta, Q, aMu, Q1b)
                End If
              EndIf             !iPrtn2 <= nPrtn1
              Qtem = Q1a + Q1b  ! Keep subtraction term in parallel with Gtem
            EndIf Heavy22
          End Select mSwG
        End If !AboveThreshold 2
      End If NaiveZM2

                               ! Accumulate gluon fusion term (always)
      Gtem = Gtem *Alpi(aMu)/2d0 !multiply the alpi here
      Call SFnConv (Gtem, g2Q(1,iPrtn2,iCmptn,iprtn1), SfT2, SfH2)
      SfTg = SfTg + SfT2
      SfHg = SfHg + SfH2
            ! Accumulate subtraction term (when applicable)
      If (SubAdd) then
        Qtem = Qtem *Alpi(aMu)/2d0 !multiply the alpi here

        Call SFnConv (Qtem, g2Q(1,iPrtn2,iCmptn,iprtn1), SfT2, SfH2)
        SfTs = SfTs + SfT2                         
        SfHs = SfHs + SfH2
        SubAdd = .FALSE.
      EndIf

      End Do ParG2
      End Do ParG1
      End If Gluon

                                ! NNLO calculation
      SfT3 = 0; SfH3 = 0; FLo3=0d0; F2o3=0d0
 
      HighOrd: If(iOrdr>1.and.ICmptn==1) then !start NNLO, neutral-current DIS
                                             !only 
        PartNNLO: Do iPrtn2 = 1, nPrtn2

          If(iFlv /= 0) then !for the semi-inclusive case
            If (iprtn2 /= abs(iFv)) Cycle ! skip channels not needed
            If (kFv /= 1) then !only semi-inclusive q+qbar is implemented 
              Call WARNI(IW, NWRT, 
     >        'STOP: NNLO semi-incl. F2 is defined only for kFv=10',
     >        'kFv', kFv, 0,0,0)
             stop  
           End If 
         End If                 !iFlv /=0


C Skip uncoupled channels or contributions below the mass threshold
            aM1 = aMass(0, iprtn2) !use pole masses in all DIS kinematics
            aM2 = aMass(0, iprtn2)
            If (g2q(0,iprtn2,iCmptn,iprtn2) < Smll .or. 
     >        .not. AboveThreshold(iprtn2,iprtn2)) Cycle 

            F2tmp3 =0d0
            FLtmp3 =0d0

            call Fnxlo(iprtn2,F2tmp3,FLtmp3);
            
            e2 = (g2Q(1, iprtn2, iCmptn, iprtn2) + g2Q(2, iprtn2, iCmptn
     >        , iprtn2))/2d0

            F2tmp3 = e2*F2tmp3; FLtmp3 = e2*FLtmp3

c  For the semi-inclusive structure function, add the real non-singlet 
c  contribution from quarks other than the tagged quark flavor.
c  This non-singlet piece contributes with the same factor both 
c  for abs(iFlv)=1...6 and ifLv=11...16

            if (iFlv /= 0 .and. SemiInclusiveF2c ) then
              do iprtn3 =1, nprtn1
                if (iprtn3 .eq. iprtn2) Cycle

                e2 = (g2Q(1, iPrtn3, iCmptn, iprtn3) + g2Q(2, iPrtn3,
     >            iCmptn, iprtn3))/2d0

                F2tmp3 = F2tmp3+ e2*SiF2l(X, Q, aMu, aM1, Ihadn, iprtn3)
                FLtmp3 = FLtmp3+ e2*SiFLl(X, Q, aMu, aM1, Ihadn, iprtn3)
              end do
            end if 

            F2o3 = F2o3 + F2tmp3; FLo3 = FLo3 + FLtmp3

cpn11 debug          write(15,'(4(1X,g12.6))') x, Q, F2o3, FLo3

        End Do PartNNLO

        Call F2LConvNC(F2o3,FLo3,SfT3,SfH3)
       
      End if HighOrd            !iOrdr>1
 
                            ! Return cumulative structure functions
      Sfn123 = (SfT0 + SfTq + SfTg - SfTs) + SfT3 
      SfnHel = (SfH0 + SfHq + SfHg - SfHs) + SfH3 

c htj
c      print*,Sfn123
c      print*,SfnHel
c      stop


      CONTAINS

      subroutine Fnxlo(iprtn,F2out,FLout)
      implicit none
      integer, Intent(in) :: iprtn
      real*8, Intent(out):: F2out, FLout
   
      real*8 :: cNS, cPS, cG, Subtraction, aM
      real*8 :: cLNS, cLPS, cLG, SubtractionL
      real*8 :: FLZMNS2, F2ZMNS2, F2lNSSub, F2ZMPS2, FLZMPS2
      real*8 :: F2LNS2, F2hlPS2, F2ZMG2, FLLNS2, F2hg2, FLZMG2, 
     >  FLhlPS2, FLhg2, F2HQSub2, FLHQSub2, 
     >  F2ZMPSsave, F2ZMGsave, FLZMPSsave, FLZMGsave,
     >  xold=-1d0, Qold=-1d0, amuold=-1d0
      integer :: iprtn3, iHadnOld = -1000

c Pavel's fix (11/16/2012) to previous bug.
      cNS=0d0; cPS=0d0; cG=0d0; cLNS=0d0; cLPS=0d0; cLG=0d0

      if ( xold /= x .or. Qold /= Q .or. 
     >  amuold /= amu .or. IhadnOld /= Ihadn) then
cpn11   Recycle ZM singlet contributions, all from ZmNNLO.f
        F2ZMPSsave= F2ZMPS2(X,Q,aMu,iHadn) 
        F2ZMGsave = F2ZMG2(X,Q,aMu,iHadn) 
        FLZMPSsave= FLZMPS2(X,Q,aMu,iHadn) 
        FLZMGsave=  FLZMG2(X,Q,aMu,iHadn) 
        xold =x; Qold=Q; amuold=amu; IhadnOld = Ihadn
      endif                     !iFlv xold

      Subtraction=0d0; SubtractionL=0d0;

      mSw2: Select Case (mSwtch)  !choose the mass scheme

      Case (3) !GM scheme
        IprtnIf: if(iprtn.lt.Nheavy) then !coefficient functions for F_l

          cLNS = FLZMNS2(X,Q,amu,iHadn,iprtn,Nheavy-1) !fLNS_{l,l,light} with Nf=Nlight
          cNS = F2ZMNS2(X,Q,amu,iHadn,iprtn,Nheavy-1)  !fNS_{l,l,light} with
                                                       ! Nf=Nlight, in GmJS.f
          do iprtn3=Nheavy,NFAboveThr                       !Contributions of heavy states

            if (.not. AboveThreshold(iprtn3,iprtn3)) Cycle !pn11: precaution, if quark masses
                                                      !are not arranged in the ascending order
            aM = aMass(0, iprtn3) !use pole masses in all DIS kinematics
            Zeta = fchi(x,Q,am,am) 
            
            cNS = cNS + F2LNS2(x,Q,aMu,aM,iHadn,iprtn) ! in GmJS.f, Eq. (25)
     >        - F2lNSSub(X ,Q,aMu,aM,iHadn,iPrtn) ! conv. of A2NSll_heavy,in GmJS.f, GmJSa.f;
                                !F2lNSsub pn must take X as the input (rescaling done inside)!
            
            cLNS = cLNS + FLLNS2(x,Q,aMu,aM,iHadn,iprtn) ! in GmJS.f, Eq. (25)
                                                         ! no subtraction c^(0)_L=0
          end do                !iprtn3

          cPS = F2ZMPSsave;  cG  = F2ZMGsave
          cLPS = FLZMPSsave; cLG  = FLZMGsave
          
        Else                    !coefficient functions for F_h
          cPS = F2hlPS2(X,Q,aMu,aM2,iHadn) !in GmJS.f
          cG  = F2hg2  (X,Q,aMu,aM2,iHadn) !in GmJS.f
          
          cLPS = FLhlPS2(X,Q,aMu,aM2,iHadn) !in GmJS.f
          cLG  = FLhg2  (X,Q,aMu,aM2,iHadn) !in GmJS.f

          if (iprtn <= Nprtn1) then 
            Zeta = fchi(x,Q,am1,am1) !X ! Rescaling variable (generalized chi)   
            cNS = F2ZMNS2(Zeta,Q,amu, iHadn,iprtn,NFAboveThr) 
            cLNS = FLZMNS2(Zeta,Q,amu,iHadn,iprtn,NFAboveThr)
          endif !iprtn <= Nprtn1

          If (Qpdf(iPrtn, iPrtn)) then ! subtraction
                                !Zeta defined inside F2HQsub2
            Subtraction = F2HQSub2(X,Q,aMu,aM2,iHadn)
            SubtractionL = FLHQSub2(X,Q,aMu,aM2,iHadn)
          End If
        End If IprtnIf

      Case (0) !ZM scheme
          cPS = F2ZMPSsave;  cG  = F2ZMGsave
          cLPS = FLZMPSsave; cLG  = FLZMGsave

         if (iprtn <= Nprtn1) then 
           cNS = F2ZMNS2(X,Q,amu,iHadn,iprtn, NFAboveThr) 
           cLNS = FLZMNS2(X,Q,amu,iHadn,iprtn,NFAboveThr) 
         endif !iprtn <= Nprtn1

       Case (1,2,8) mSw2 !Improved ZM scheme
        IprtnIf2: if(iprtn.lt.Nheavy) then !coefficient functions for F_l

          cNS = F2ZMNS2(X,Q,amu,iHadn,iprtn,Nheavy-1) !fNS_{l,l,light} with
          cLNS= FLZMNS2(X,Q,amu,iHadn,iprtn,Nheavy-1) ! Nf=Nlight, now in GmJS.f
 
          do iprtn3=Nheavy,Nprtn1 !Contributions of heavy states
            if (.not. AboveThreshold(iprtn3,iprtn3)) Cycle
            aM = aMass(0, iprtn3) !use pole masses in all DIS kinematics
            Zeta = fchi(x,Q,am,am) 
            cNS = cNS + F2ZMNS2(Zeta,Q,amu,iHadn,iprtn,1)
            cLNS = cLNS + FLZMNS2(Zeta,Q,amu,iHadn,iprtn,1)

          end do                !iprtn3

          cPS = F2ZMPSsave;  cG  = F2ZMGsave
          cLPS = FLZMPSsave; cLG  = FLZMGsave

        Else                             !coefficient functions for F_h
          Zeta = fchi(x,Q,am1,am1) !X ! Rescaling variable (generalized chi)                                  

          cPS = F2ZMPS2(Zeta,Q,aMu,iHadn)
          cG  = F2ZMG2 (Zeta,Q,aMu,iHadn)

          cLPS = FLZMPS2(Zeta,Q,aMu,iHadn)
          cLG  = FLZMG2 (Zeta,Q,aMu,iHadn)

          if (iprtn <= Nprtn1) then 
            cNS = F2ZMNS2(Zeta,Q,amu,iHadn,iprtn,Nprtn1) !Nf=Nprtn1, now in GmJS.f
            cLNS = FLZMNS2(Zeta,Q,amu,iHadn,iprtn,Nprtn1) !Nf=Nprtn1, now in GmJS.f
          endif                 !iprtn <= Nprtn1

        End If IprtnIf2
  
      End Select mSw2

      F2out = cNS + cPS + cG - Subtraction
      FLout = cLNS + cLPS + cLG - SubtractionL

      End Subroutine Fnxlo

      subroutine ExamineMassThresholds()
cpn11 1. Sets AboveThreshold(ip21,ip22)=.TRUE. for each pair of quark flavors
c   ip1 and ip2 that can be kinematically produced at the given energy W, 
c   i.e., W > aM(ip21) + aM(ip22). Computes the total number NFAboveThr
c   of quark-antiquark pairs (with the same flavor) that can be
c   physically produced from gluon splittings at this energy (needed for
c   the non-singlet NNLO contribution)
      integer ip21, ip22
      
      NFAboveThr = 0
      do ip21=1, Nprtn2
        do ip22=1, Nprtn2
          AboveThreshold(ip21,ip22) = 
     >     ( W > (aMass(0,ip21) + aMass(0,ip22))*(1d0+sml))
        enddo                   !ip22
       
        if (AboveThreshold(ip21,ip21)) NFAboveThr=NFaboveThr+1
      enddo                     !ip21      

      End Subroutine ExamineMassThresholds
      
      End Subroutine SfnGm         ! =====  *****************

      Real*8 function fchi (x,Q,am1,am2) 
! For given x, Q, particle masses aM1 and M2, and shape parameter achi
! returns the rescaled light-cone variable chi. This function must be
! called only above the heavy-quark threshold
!
! Valid ranges of achi    Returned value 
! 0.0 :                   Original rescaling variable chi
! 0.0:1.0:       Generalized rescaling variable zeta(lambda) from 
!                arXiv:0903.2667 by Nadolsky and Tung
! 100:            zeta=x
! Other values:   abort with an error message 

      implicit none
      Real*8,  Intent(in) :: X, Q, am1,am2
      real*8, external:: zbrnt, zetadiff
      real*8 :: xtmp, chi0, rerr=1e-10, eps, ax, achi1, akappa, 
     > sml=1.0e-8
      Integer :: iOrdr, iSchm
      Real*8  :: Q0, ScalFc, achi 
      Common / cDisCal / Q0, ScalFc, achi, iOrdr, iSchm
      integer :: irt
      common /cchi /ax, achi1, akappa
      integer iw,  NIN, NOUT, NWRT
      data iw /0/
      COMMON / IOUNIT / NIN, NOUT, NWRT

      
      akappa = (aM1+aM2)**2/Q**2

      chi0=x*(1e0+akappa)       !chi variable from hep-ph/0110247
                                !by Tung, Kretzer, Schmidt
      if (chi0 >=1) stop 'STOP in fchi: chi0 > 1'

      ax = x; achi1 = achi      
c 	print *,'tmp: achi1=',achi1

      CheckAchi: if (abs(achi1) <= sml) then
        fchi = chi0
        return
      else if (abs(achi-100.0).lt.sml) then
        fchi = x
        return
      else if (achi1.gt.0d0.and.achi1.lt.1d0) then
cpn2009 Sanity check: the relation between chi and x is one-to-one for
c      chi<1 only if achi < 1 + 1/akappa
        if (achi >= 1e0 + 1e0/akappa) then
          print *, 'STOP in fchi: x is not a single-value function',
     >      ' of chi in the region chi<1' 
          print *, 'achi, akappa =', achi, akappa
	  stop		!missing in previous versions jcp
        endif
        eps = x*rerr
        xtmp = zbrnt(zetadiff, x, chi0, eps, irt)
        if (irt /= 0) then
          print *, 'Error in fchi: irt = ',irt
          print *, 'for x, Q, aM1, aM2 = ',x,Q,aM1,aM2
          stop
        End If
      else CheckAchi
       Call WARNR(IW, NWRT,
     >    'STOP: Wrong value of the rescaling parameter lambda',
     >    'achi', achi, 0,0,0)
       stop
      End if CheckAchi

      fchi = xtmp
      End function fchi


      Real*8 function zetadiff(zeta)
      implicit none
      real*8, intent(in) :: zeta
      real*8 :: ax,achi1, akappa
      common /cchi /ax,achi1, akappa

      zetadiff = zeta-ax*(1d0 + akappa* zeta**achi1)

      End function zetadiff


