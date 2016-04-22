      Function StrFnGm (iHadn,iCmptn,X,Q, iSfn, iFlv)
C     ! Invokes SFnGm.f to get the full set of tensor and helicity S.F.s F_{1,2,3} and F_{-1,0,1};
C     ! Returns the relevant S.F. to the calling program, according to the argument iSfn:

C     ! iSfn : 1, 2, 3    ==> F_{1,2,3}/Prefactor
C     !       -1, 0, -2   ==> F_{left, longitudinal, right}/Prefactor  ("1" is used for F_1)
C     !
C     ! Prefactor is defined in CalculateDIS (2*x for F_L, x otherwise)
C
C     ! iFlv : 0   (total inclusive---sum over all active quark-parton flavors);
C     !       -6:6 (individual flavors);
C     !       11:16 (quark+antiquark).                            ! ---------------------!


C 2006.12.06:
C             ! New input argument iFlv introduced; common block to store heavy flavor S.F.s deleted.
C             ! Logic of the calculation of heavy flavor S.F.s drastically simplified by HLL.
C             ! New option for choice of factorization scale (iFscle=3) for HQ production added.

C 2005.12.25: original version.

C                            ! ---------------------!

      Implicit None

      Integer, Intent(in) :: iHadn, iCmptn, iSfn, iFlv
      Real*8,  Intent(in) :: X, Q
      Real*8, PARAMETER:: Ssml=1d-20

                                          ! Total Inclusive StrFn's:
      Real*8 :: Sfn123(3), SfnHel(-1:1), StrFnGm
      Real*8 :: Tem
      
                          ! ---------------------------
                                                      ! Calculate the structure functions
      Call SfnGm(iHadn,iCmptn,iFlv, X,Q, Sfn123,SfnHel)
                                                ! Assign the relevant SF to the function call
      Strfn: Select case (iSfn)
      Case (0)                 ! Longitudinal Structure function
        Tem = SfnHel(0)
      Case (-1)                 ! F_Left
        Tem = SfnHel(-1)
      Case (-2)                 ! F_Right
        Tem = SfnHel(1)
      Case (4)                 ! Ratio = sig_long/sig_trans = F_L / (2xF_1) = F_0/((F_+ F_-)/2)
        Tem = SfnHel(0)/ ((SfnHel(1)+SfnHel(-1))/2d0)
      Case (1,2,3)             ! Conventional F_1,2,3
        Tem = Sfn123(iSfn)
      Case Default
        Print *, 'iSfn = ', iSfn, ' is not allowed in StrFnGm!'
      End Select Strfn

      If (Abs(Tem) < Ssml) Tem = Ssml     ! Avoid returning zero

      StrFnGm = Tem

      Return
                         ! **************************
      End Function StrFnGm
