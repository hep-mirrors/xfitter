      Module mDisPac
   
      Implicit None

      Logical :: DisSet = .false.

                          ! Structure Function level: Virtual Forward Compton Amplitude 
                          ! incoming Boson=iBsn(1) and outgoing boson=iBsn(2);
                          ! iBsn = 1/2/3/4/  for photon/W+/W-/Z;
      ! iCmptn = 1/2/3/4 for the diagonal processes; and 5 for photon-Z interference
      !          in addition, -1, 0 for W+ -/+ W- .
   
                                                             ! QCD constants; pi
      Real*8,PARAMETER:: CF =4d0/3d0, CG =3d0, TR =5d-1, 
     >  Pi=3.1415926536d0, Tf=0.5D0
      Real*8, parameter ::  XsGv2pb = 0.389379323D9          ! GeV^2 pb
     >                     ,Gf = 1.16637D-5                  ! G-Fermi
C                                                   ! Global: Integration parameters    
      Real*8  :: aErr = 1d-16, rErr = 1d-6 
cmarco rErr must be set to 1d-5 or 1d-6 to avoid spikes in F2b
c jcp that may be cured by using Simpson integration in routines of flib/3dis.
c jcp However, when the Simpson methods are used, rErr = 1d-6 here doesn't add much extra time.

              ! Electroweak couplings for Dis Structure Functions, 
              ! They are set in subroutine mSetDis

      Integer :: iBsn(2, -1:5)         ! vector boson pair for each of the iCmptn choices
      
      Real*8  :: g2Bn(-1:5)            ! Gauge coupling squared
      
      Real*8  :: g2Lp(-2:2, 5), g2Lm(-2:2, 5)  ! parity conserving/violating Lepton coupling combinations

      Real*8  g2p(-6:6, 5, -6:6)       ! parity conserving combination   g_V g_V+g_A g_A=g_L g_L + g_R g_R
     >      , g2m(-6:6, 5, -6:6)       ! parity violating  combination   g_V g_A+g_A g_V=g_L g_L - g_R g_R
     >      , g2x(-6:6, 5, -6:6)       ! chirality violating combination g_V g_V-g_A g_A=g_L g_R + g_R g_L
     >      , g2q(0:3, -6:6, -1:5, -6:6)       ! All above combined; 0th element = Sum Abs
     >      , g2Qpc (-6:6, 5, 6)    ! Accumulated parity conserving combination for quark partons
     >      , g2Qpv (-6:6, 5, 6)    ! Accumulated parity violating  combination for quark partons
     >      , g2Gpc    (6, 5, 6)   ! Accumulated parity conserving combination for gluon parton
     >      , g2Gpv    (6, 5, 6)   ! Accumulated parity violating  combination for gluon parton
  
                                             ! Some convenient arrays for initiation   
      Real*8 :: AA( 1:3) = (/ 1d0, 2d0, 3d0 /) 
     >        , AB(-1:1) = (/-1d0, 0d0, 1d0 /)
   
      !! The definition of structure functions calculated in this DisPac is the same as
      !! that of Furmanski and Petronzio. These are related to the "conventional" SFs by:
      !!        F_1 = 2 * F_1^conv ;
      !!        F_2 = F_2^conv / x ;
      !!        F_3 = F_3          .
      !! The relation between these "tensor" SFs and the helicity SFs are:
      !!        F_1 = F_L + F_R   ;   F_L = (F_1 + F_3) / 2   ;
      !!        F_2 = F_1 + F_0 *2;   F_0 = (F_2 - F_1) / 2   ;
      !!        F_3 = F_L - F_R   ;   F_R = (F_1 = F_3) / 2   .   

      End Module mDisPac             !  ***************************
  
      Subroutine mSetDis

                     ! Compute relevant EW coupling combinations for DIS processes.
                     ! Call the Ewk routine once, and put results in Module mDisPac
                                    ! Note: the helicity couplings in EwkPac are a factor of
                                    ! 2 different from those of the ACOT1 paper.
      Use mDisPac
      
      Implicit None
C                                                ! Local variables
      Integer :: jCmptn, iBn, nFlv, k, iPrtn1, iPrtn2, iLep1, iLep2
      Real*8  :: gL(2), gR(2)
      Real*8  :: g2Qpc_tmp1, g2Qpc_tmp2, g2Qpv_tmp1, g2Qpv_tmp2
     >         , g2Gpc_tmp, g2Gpv_tmp, s = 1d-3
     
      Integer :: iLptn2(2,5), j
      
C                                                ! External function 
      Integer:: NflTot                           ! (from QcdPac)
      Real*8 :: gEWqH, gEWLH, AlfEwk                            ! (from EwkPac)

      nFlv = NflTot()               ! in principle
      nFlv=6                        ! in practice, Use the full range in the SM
      
      Do iBn = 1, 4
        g2Bn(iBn) = AlfEwk(iBn) * 4d0 * Pi
      End Do
      g2Bn(-1:0) = g2Bn(2)   ! these are for W+/W- combos => the same gauge coupling
      
      Do jCmptn = 1, 5       ! The forward Compton Amplitudes have 5 indep. combinations
                             ! of incoming and outgoing vector bosons:
        If (jCmptn <= 4) then
            iBsn(1:2, jCmptn) = (/jCmptn, jCmptn/)
        ElseIf (jCmptn == 5) then
            iBsn(1:2, jCmptn) = (/1, 4/)
        EndIf
      End Do
                       ! Compute and store Lepton-vector-boson coupling combinations
      g2Lp = 0d0
      g2Lm = 0d0
      iLptn2 (1,1:5) = (/1, 2, 0, 1, 1/)     ! final-state lepton for iLptn1 = 1
      iLptn2 (2,1:5) = (/2, 0, 1, 2, 2/)     ! final-state lepton for iLptn1 = 2
      Leptn: Do iLep1 = 1, 2
        Do jCmptn = 1, 5
          iLep2 = iLptn2(iLep1,jCmptn)
                                     ! For the lepton vertex, the Boson is "outgoing".
                                     ! In the EwkPac def., the boson is "incoming".
                                     ! Therefore, we need to interchange W+/W- here:
          If (jCmptn == 2) Then
             j = 3
          Else If (jCmptn == 3) Then
             j = 2
          Else
             j = jCmptn
          End If                           
          
          If (iLep2 /= 0) Then
             gL(1) = gEWlH ( iLep1, iBsn(1,j), -1,  iLep2) 
             gR(1) = gEWlH ( iLep1, iBsn(1,j),  1,  iLep2)
             If(jCmptn == 5) then
                gL(2)= gEWlH ( iLep1, iBsn(2,j), -1,  iLep2)
                gR(2)= gEWlH ( iLep1, iBsn(2,j),  1,  iLep2)
             Else
                gL(2)=gL(1)
                gR(2)=gR(1)
             Endif
          Else
             gL = 0
             gR = 0
          End If
      
          g2Lp(iLep1, jCmptn) = !HLL 2d0*
     >      (gL(1)*gL(2) + gR(1)*gR(2))   ! parity conserving coupling combination
                                           ! g_V g_V + g_A g_A (1,2 labels omitted)
          g2Lm(iLep1, jCmptn) = !HLL 2d0*
     >      (gL(1)*gL(2) - gR(1)*gR(2))   ! parity violating  coupling combination
                                           ! g_V g_A + g_A g_V
          g2Lp(-iLep2, jCmptn) = g2Lp(iLep1, jCmptn) 
          g2Lm(-iLep2, jCmptn) =-g2Lm(iLep1, jCmptn) 
        End Do
      End Do Leptn
       
                ! Compute and store vector-boson quark scattering coupling combinations
      Part1: Do iPrtn1 = -nFlv, nFlv
      Part2: Do iPrtn2 = -nFlv, nFlv
      Bsns : Do jCmptn = 1, 5
  
        Do k = 1, 2
        If (iPrtn1 > 0 .and. iPrtn2 > 0) Then
          gL(k) = gEWqH ( iPrtn1, iBsn(k,jCmptn), -1,  iPrtn2) 
          gR(k) = gEWqH ( iPrtn1, iBsn(k,jCmptn),  1,  iPrtn2)
        ElseIf (iPrtn1 < 0 .and. iPrtn2 < 0) Then
          gL(k) = gEWqH (-iPrtn2, iBsn(k,jCmptn),  1, -iPrtn1) 
          gR(k) = gEWqH (-iPrtn2, iBsn(k,jCmptn), -1, -iPrtn1) 
        Else
          gL(k) = 0d0
          gR(k) = 0d0
        EndIf
        End Do
C                                            Compute coupling combinations and store in mDisPac                                      
        g2p(iPrtn1, jCmptn, iPrtn2) = 2d0*
     >   (gL(1)*gL(2) + gR(1)*gR(2))   ! parity conserving coupling combination
                                           ! g_V g_V + g_A g_A (1,2 labels omitted)
        g2m(iPrtn1, jCmptn, iPrtn2) = 2d0*
     >   (gL(1)*gL(2) - gR(1)*gR(2))   ! parity violating  coupling combination
                                           ! g_V g_A + g_A g_V
        g2x(iPrtn1, jCmptn, iPrtn2) = 2d0*
     >   (gL(1)*gR(2) + gR(1)*gL(2))   ! chirality violating coupling combination
                                           ! g_V g_V - g_A g_A (1,2 labels omitted)
        
                      ! ---- Combine them ------
        g2q(1,iPrtn1,jCmptn,iPrtn2) = g2p(iPrtn1,jCmptn,iPrtn2)
        g2q(2,iPrtn1,jCmptn,iPrtn2) = g2x(iPrtn1,jCmptn,iPrtn2)
        g2q(3,iPrtn1,jCmptn,iPrtn2) = g2m(iPrtn1,jCmptn,iPrtn2)
                      ! convenient measure of the "strength" of the coupling
                      ! (Useful as test for uncoupled channels.)
        g2q(0,iPrtn1,jCmptn,iPrtn2)     = g2p(iPrtn1,jCmptn,iPrtn2)**2
     >   + g2x(iPrtn1,jCmptn,iPrtn2)**2 + g2m(iPrtn1,jCmptn,iPrtn2)**2
      EndDo Bsns
                ! Form (W+ +/- W-) combinations, useful for some neutrino experiments
                ! Notice that the parity-violating F_3 changes sign in these combinations
      g2q(1:3,iPrtn1, 0,iPrtn2) = 0.5d0 *
     >  (g2q(1:3,iPrtn1,2,iPrtn2) + g2q(1:3,iPrtn1,3,iPrtn2)*(/1,1,-1/))
      g2q(1:3,iPrtn1,-1,iPrtn2) = 0.5d0 *
     >  (g2q(1:3,iPrtn1,2,iPrtn2) - g2q(1:3,iPrtn1,3,iPrtn2)*(/1,1,-1/))
     
      g2q(0,iPrtn1,-1:0,iPrtn2) = g2q(1,iPrtn1,-1:0,iPrtn2)**2 +
     >   g2q(2,iPrtn1,-1:0,iPrtn2)**2 + g2q(3,iPrtn1,-1:0,iPrtn2)**2  
      EndDo Part2
      EndDo Part1  
                    ! Now compute the combination of EW couplings that enter the 
                    ! DIS SF calculation for each contributing incoming parton
      g2Gpc = 0
      g2Gpv = 0
      g2Qpc = 0
      g2Qpv = 0
      
      Comptn: Do jCmptn=1,5
        
        g2Gpc_tmp=0
        g2Gpv_tmp=0
        
        Quark1: Do iPrtn1 = 1, nFlv             ! For each flavor, fill the quark and
                                                ! and the anti-quark entries in parallel
        g2Qpc_tmp1=0
        g2Qpv_tmp1=0
        g2Qpc_tmp2=0
        g2Qpv_tmp2=0
        Quark2: Do iPrtn2 = 1, nFlv
                                   ! make sure coefficients have the requisite symmetry
      if(abs(g2p(iPrtn1,jCmptn,iPrtn2)-g2p(-iPrtn2,jCmptn,-iPrtn1))>s
     >.or.abs(g2x(iPrtn1,jCmptn,iPrtn2)-g2x(-iPrtn2,jCmptn,-iPrtn1))>s
     >.or.abs(g2m(iPrtn1,jCmptn,iPrtn2)+g2m(-iPrtn2,jCmptn,-iPrtn1))>s
     >) then
      print '(A, 3I5, 3(/2F15.7))', 
     > 'Bomb: ', iPrtn1, jCmptn, iPrtn2
     > , g2p( iPrtn1, jCmptn, iPrtn2),g2p(-iPrtn2, jCmptn,-iPrtn1)
     > , g2x( iPrtn1, jCmptn, iPrtn2),g2x(-iPrtn2, jCmptn,-iPrtn1)
     > , g2m( iPrtn1, jCmptn, iPrtn2),g2m(-iPrtn2, jCmptn,-iPrtn1)
      End If
                                   ! Quarks: for each flavor, sum over the final state iPrtn2:
          g2Qpc_tmp1 = g2Qpc_tmp1 + g2p( iPrtn1, jCmptn, iPrtn2)
          g2Qpv_tmp1 = g2Qpv_tmp1 + g2m( iPrtn1, jCmptn, iPrtn2)
                                   !  - iPrtn1 case is obtained by CP symmetry, when needed       
          g2Qpc_tmp2 = g2Qpc_tmp2 + g2p(-iPrtn1, jCmptn,-iPrtn2)
          g2Qpv_tmp2 = g2Qpv_tmp2 + g2m(-iPrtn1, jCmptn,-iPrtn2)
       
          g2Qpc( iPrtn1, jCmptn, iPrtn2) = g2Qpc_tmp1
          g2Qpv( iPrtn1, jCmptn, iPrtn2) = g2Qpv_tmp1      
          g2Qpc(-iPrtn1, jCmptn, iPrtn2) = g2Qpc_tmp2          ! Here iPrtn2 is the flavor #
          g2Qpv(-iPrtn1, jCmptn, iPrtn2) = g2Qpv_tmp2          ! It is always positive

                                ! Gluon: for each flavor, sum over both quark and anti-quark
          g2Gpc_tmp = g2Gpc_tmp + g2p(  iPrtn1,  jCmptn,   iPrtn2)
          g2Gpv_tmp = g2Gpv_tmp + g2m(  iPrtn1,  jCmptn,   iPrtn2)
       
          g2Gpc_tmp = g2Gpc_tmp + g2p( -iPrtn1,  jCmptn,  -iPrtn2)
          g2Gpv_tmp = g2Gpv_tmp + g2m( -iPrtn1,  jCmptn,  -iPrtn2)
        
          g2Gpc( iPrtn1, jCmptn, iPrtn2)=g2Gpc_tmp
          g2Gpv( iPrtn1, jCmptn, iPrtn2)=g2Gpv_tmp
        
        EndDo Quark2
        
        EndDo Quark1
      Enddo Comptn
     
      DisSet = .true.
  
      End Subroutine mSetDis            ! ******************************
         
