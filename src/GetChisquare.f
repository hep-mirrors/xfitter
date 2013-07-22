*     ---------------------------------------------------------  	 
*     Calculate chisquare:
*      - first get error matrix
*      - invert error matric to get errors
*      - calculate chisquare
*     ---------------------------------------------------------

      subroutine GetNewChisquare(flag_in,n0_in,fchi2_in,rsys_in,ersys_in,pchi2_in,
     $     fcorchi2_in)
      implicit none

      include 'ntot.inc'
      include 'steering.inc'
      include 'systematics.inc'
      
      integer n0_in, flag_in
      double precision fchi2_in, ERSYS_in(NSYSMax), RSYS_in(NSYSMax)
      double precision pchi2_in(nset), fcorchi2_in

      double precision chi2_log
      integer i,j,jsys,k

      double precision ScaledGamma(NSysMax,Ntot) ! Scaled Gamma matrix
      double precision ScaledGammaSav(NSysMax,Ntot) ! Scaled Gamma matrix, saved

      double precision ScaledOmega(NSysMax,Ntot) ! Scaled Omega matrix
      
      double precision ScaledErrors(Ntot)  ! uncorrelated uncertainties, diagonal

      double precision ScaledErrorMatrix(NCovarMax,NCovarMax) ! stat+uncor error matrix
      double precision ScaledSystMatrix(NCovarMax,NCovarMax)  ! syst. covar matrix
      double precision ScaledTotMatrix(NCovarMax,NCovarMax)   ! stat+uncor+syst covar matrix

      integer NDiag, NCovar   ! Number of diagonal and full covariance input data points
      integer List_Diag(Ntot), List_Covar(Ntot), List_Covar_inv(Ntot)

      integer Iterate
      logical LFirst
      data LFirst /.true./

      integer omegaIteration 
      Logical doMatrix, doNuisance, doExternal, LStop
      
C----------------------------------------------------------------------------

         ! --> WS debug
         if(LDEBUG) print*,'GetNewChisquare flag_in=',flag_in
         
c Global initialisation 

      if (LFirst) then
         LFirst = .false.

C    !> Determine which mechanisms for syst. errors should be used:
         Call Init_Chi2_calc(doMatrix, doNuisance, doExternal) 

C    !> Determine which errors are diagonal and which are using covariance matrix
         Call init_chi2_stat(NDiag, NCovar, List_Diag, List_Covar,
     $        List_Covar_inv,n0_in)

         ! print*,'NDiag=',NDiag,'NCovar=',NCovar

         do i=1,NCovarMax
            do j=1,NCovarMax
               ScaledSystMatrix(i,j) = 0.
            enddo
         enddo
      endif  ! LFirst
      
C 
      do jsys=1,nsys
         rsys_in(jsys) = 0.d0
         ersys_in(jsys) = 0.d0
      enddo

      do i=1,nset
         pchi2_in(i)=0.d0
      enddo
      
      fchi2_in = 0.d0
      fcorchi2_in = 0.d0

C   Determine if we need to iterate for stat. errors:
      if (Chi2ExtraSystRescale) then
         Iterate = 1
      else
         Iterate = 0
      endif

      if ( Chi2FirstIterationRescale .and. flag_in.gt.1 ) then
  ! Reset iterations:
         Iterate = 0
      endif


C !> Read external (minuit) systematic sources if present:
      if (doExternal) then
         call Chi2_calc_readExternal( rsys_in, ersys_in, flag_in )
      endif

      if (.not. Chi2FirstIterationRescale  .or. flag_in.eq.1) then
C !> Calculated scaled syst. uncertainties:
         call Chi2_calc_GetGamma(ScaledGamma, ScaledOmega)
         
C !> Store rescaled gamma (important for asymmetric errors ):
         do k=1,nsys
            do i=1,n_syst_meas(k)
               j =  syst_meas_idx(i,k)
               ScaledGammaSav(k,j) = ScaledGamma(k,j)
            enddo
         enddo            


        ! print *,' --- ScaledGamma'
        ! do i=1,n0_in
          ! print *,(ScaledGamma(j,i),j=1,nsys)
        ! enddo


C !> Rebuild syst. covariance matrix 

         if ( doMatrix ) then
            Call Chi2_calc_covar(ScaledGamma
     $           ,ScaledSystMatrix
     $           ,List_Covar_Inv,n0_in)
         endif
          ! print *,' --- ScaledSystMatrix'
          ! do i=1,6
            ! print *,(ScaledSystMatrix(j,i),j=1,6)
          ! enddo
      else 
C !> Restore saved gamma:
         do k=1,nsys
            do i=1,n_syst_meas(k)
               j =  syst_meas_idx(i,k)
               ScaledGamma(k,j) = ScaledGammaSav(k,j)
            enddo
         enddo
      endif


C !> Get uncor errors/nuisance parameters.
      do while ( Iterate.ge.0 )

c !> First recalc. stat. and bin-to-bin uncorrelated uncertainties:
         if (.not. Chi2FirstIterationRescale .or. flag_in.eq.1 .or.
     $  Chi2OffsRecalc) then
            Call Chi2_calc_stat_uncor(ScaledErrors
     $           ,ScaledErrorMatrix
     $           ,rsys_in,n0_in, NCovar, List_Covar, Iterate)
     
          ! print *,' --- ScaledErrors'
            ! print *,(ScaledErrors(j),j=1,6)
          ! print *,' --- ScaledErrorMatrix'
          ! do i=1,6
            ! print *,(ScaledErrorMatrix(j,i),j=1,6)
          ! enddo

C  !> Sum covariance matricies and invert the total:

            if ( doMatrix .or. NCovar .gt. 0 ) then
               Call Chi2_calc_SumCovar(ScaledErrorMatrix, 
     $              ScaledSystMatrix, 
     $              ScaledTotMatrix, NCovar)
     
              ! print *,' --- ScaledTotMatrix Inv.'
              ! do i=1,6
                ! print *,(ScaledTotMatrix(j,i),j=1,6)
              ! enddo
            endif

C !> same for diagonal part:
            do i=1,n0_in
               ScaledErrors(i) = 1.D0 
     $              / (ScaledErrors(i)*ScaledErrors(i))
            enddo

         endif


C !> Next determine nuisance parameter shifts
         omegaIteration = 1
         do 
            Call Chi2_calc_syst_shifts(
     $           ScaledErrors
     $           ,ScaledTotMatrix
     $           ,ScaledGamma
     $           ,rsys_in,ersys_in,list_covar_inv, flag_in, n0_in)

C !> Asymmetric errors loop:
            Call UseOmegaScale(ScaledGamma
     $           ,ScaledGammaSav
     $           ,ScaledOmega
     $           ,rsys_in
     $           ,omegaIteration,
     $           LStop)         
            if (LStop) Exit 
            omegaIteration = omegaIteration + 1
         enddo 


C !> See if we want to use asymmetric errors

         Iterate = Iterate - 1
      enddo   ! while ( Iterate.ge.0 )

C !> Calculate chi2
      call chi2_Calc_chi2(
     $     ScaledErrors,
     $     ScaledGamma,
     $     ScaledTotMatrix,
     $     rsys_in,
     $     ndiag, list_diag, ncovar, list_covar,
     $     fchi2_in, pchi2_in, fcorchi2_in)


C !> Add log term
      if ( Chi2PoissonCorr ) then
         call chi2_calc_PoissonCorr(ScaledErrors, chi2_log, n0_in)
         fchi2_in = fchi2_in + chi2_log
      endif

       ! print*,'fchi2_in=',fchi2_in


C
C !> Store extra output for FCN = 3:
C
      if (Flag_In.eq.3) then
         Call Chi2_calc_FCN3(ScaledErrors,ScaledGamma,RSys_in,n0_in)
      endif


      return 
      end


C------------------------------------------------------------------
C
C> @brief Check for of systematic uncertainties, which methods should be used
C
C> @param doMatrix switch on covariance matrixw method
C> @param doNuisance switch on hessian method
C> @param doExternal switch on external (minuit) method
C
C------------------------------------------------------------------
      subroutine Init_Chi2_calc(doMatrix, doNuisance, doExternal) 

      implicit none
      logical doMatrix, doNuisance, doExternal
      ! logical doOffset
      include 'ntot.inc'
      include 'systematics.inc'
      integer k, n_m, n_n, n_e, n_o
      character*64 Msg
C----------------------
      doMatrix   = .false.
      doNuisance = .false.
      doExternal = .false.
      doOffset = .false.
      n_m = 0
      n_n = 0
      n_e = 0
      n_o = 0
      do k=1,nsys
         if ( SysForm(k) .eq. isMatrix) then
            doMatrix = .true.
            n_m = n_m + 1
            if (n_m .gt. NCovarMax ) then
               print *,'ERROR ERROR ERROR'
               print *,'Number of points used for covariance matrix'//
     $              ' exceeds NCovarMax = ', NCovarMax
               print *,'Increase NCovarMax in ntot.inc and recompile'
               print *,'STOP'
               call HF_stop
            endif
         endif
         if ( SysForm(k) .eq. isNuisance) then
            doNuisance = .true.
            n_n = n_n + 1
         endif
         if ( SysForm(k) .eq. isExternal) then
            doExternal = .true.
            n_e = n_e + 1
         endif
         if ( SysForm(k) .eq. isOffset) then
            doOffset = .true.
            n_o = n_o + 1
         endif

      enddo
C
C Add some info messages:
C
      if (doMatrix) then
         write (Msg,'(''I: Use matrix method for '',i4,'' sources'')')
     $        n_m
         call HF_errlog(271120122,Msg)
      endif
      if (doNuisance ) then
         write (Msg,
     $  '(''I: Use hessian method for'',i4,'' sources'')') n_n
        call HF_errlog(271120123,Msg)
      endif
      if (doOffset) then
         write (Msg,
     $  '(''I: Use offset method for'',i4,'' sources'')') n_o
        call HF_errlog(70120131,Msg)
      endif

      if (doExternal) then
         write (Msg,
     $  '(''I: Use external (minuit) method for'',i4,'' sources'')') n_e
        call HF_errlog(271120124,Msg)
      endif

      end


C-----------------------------------------------------------------------------------------
C
C> @brief Calculate number of diagonal and covariance input data and create corresponding lists
C
C> @param NDiag number of diagonal input data
C> @param NCovar number of covariance input data
C> @param List_Diag list of diagonal input data
C> @param List_Covar list of covariance input data
C> @param List_Covar_inv inverted list of covariance input data
C> @param n0_in total number of input data
C
C------------------------------------------------------------------------------------------
      subroutine Init_chi2_stat(NDiag, NCovar, List_Diag, List_Covar, 
     $     List_Covar_inv, n0_in)

      implicit none
      include 'ntot.inc'
      include 'indata.inc'
      include 'systematics.inc'
      include 'steering.inc'
      integer NDiag, NCovar, List_Diag(NTOT), List_Covar(NTOT), 
     $     List_Covar_inv(NTOT), n0_in
      integer i,k,l
      logical isCov
C-----------------------------------------------------------------------

      NDiag  = 0
      NCovar = 0
      do i=1,n0_in
         List_Covar_Inv(i) = 0   ! Reset inverted list
         isCov = .false.

C Check if already requested to be 
         if ( is_covariance(i) ) then
            isCov = .true.
         else
C Check systematic sources, if a matrix source point to point i
             do k=1,nsys
                if (SysForm(k) .eq. isMatrix) then
                   do l=1,n_syst_meas(k)
                      if (syst_meas_idx(l,k).eq.i) then
                         isCov = .true.
                      endif
                   enddo
                endif 
             enddo
         endif


         if (isCov) then
            NCovar = NCovar + 1

            List_Covar(NCovar) = i
            List_Covar_inv(i) = NCovar
            ! Add a check:
            if (NCovar .gt. NCovarMax) then
               print *,'ERROR ERROR ERROR'
               print *,'Number of points used for covariance matrix'//
     $              ' exceeds NCovarMax = ', NCovarMax
               print *,'Increase NCovarMax in ntot.inc and recompile'
               print *,'STOP'
               call HF_stop
            endif
         else
            NDiag = NDiag + 1
            List_Diag(NDiag) = i
         endif
         
      enddo

      if (lDebug) then
         print *,'DEBUG from Init_chi2_stat'
         print *,'Ncovar= ',Ncovar,' NDiag=',NDiag
      endif

      end

C----------------------------------------------------------
C
C> @brief Get external (minuit) parameters 
C
C> @param rsys_in
C> @param ersys_in
C> @param IFlag
C
C---------------------------------------------------------

      subroutine Chi2_calc_readExternal( rsys_in, ersys_in, IFlag )

      implicit none
      include 'ntot.inc'
      include 'systematics.inc'
      include 'endmini.inc'
      include 'extrapars.inc'
      double precision rsys_in(nsysmax), ersys_in(nsysmax)
      integer i,idx, iflag
      integer GetParameterIndex
      character*80 parname
      double precision val, err, xlo, xhi
      integer ipar
C-------------------------------------------------
      do i=1,NSys
         if (SysForm(i) .eq. isExternal ) then
            idx = GetParameterIndex(system(i))
            if (idx.eq.0) then
               print *,'ERROR ERROR ERROR'
               print *,'External systematics ',system(i),' not found'
               print *,'on the list of external parameters'
               print *,'Contact herafiter-help@desy.de with ! Stop.'
               print *,' ' 
               call HF_stop
            endif
            rsys_in(i) = parminuitsave( iExtraParamMinuit(idx) )
            if (IFlag.eq.3) then
               call mnpout(iExtraParamMinuit(idx)
     $              ,parname,val,err,xlo,xhi,ipar)
               ersys_in(i) = err
            endif
         endif
      enddo

      end

C---------------------------------------------------------
C
C> @brief Calculate re-scaled effect of systematic error sources
C
C> @param ScaledGamma
C> @param ScaledOmega
C
C---------------------------------------------------------

      subroutine chi2_calc_GetGamma(ScaledGamma, ScaledOmega)

      implicit none

      include 'ntot.inc'
      include 'systematics.inc'

      double precision ScaledGamma(NSysMax,Ntot) !> Scaled Gamma matrix
      double precision ScaledOmega(NSysMax,Ntot) !> Scaled Omega matrix

      include 'indata.inc'
      include 'theo.inc'
c      include 'steering.inc'

      integer i1,j1,i,j,k
      integer scaling_type
      double precision scale
C-----------------------------------------------------
      do k=1,NSYS
         scaling_type = SysScalingType(k) 

         do i1=1,n_syst_meas(k)
            i = syst_meas_idx(i1,k)

            if (scaling_type.eq. isNoRescale) then
               scale = daten(i)
            elseif (scaling_type.eq. isLinear) then
               scale = theo(i)
            elseif (scaling_type.eq. isPoisson) then
               scale = sqrt(theo(i)*daten(i))
            elseif (scaling_type.eq. isLogNorm) then
               scale = theo(i)
               call HF_errlog(271120121,
     $              'S: Not implemented LogNormal scaling requested.'//
     $              ' Using linear instead')
            endif

            ! The scaled syst. errors:
            ScaledGamma(k,i) = beta(k,i)*scale
            ScaledOmega(k,i) = omega(k,i)*scale
         enddo
      enddo
C-----------------------------------------------------
      end

C--------------------------------------------------------------------------------------------
C
C> @brief Calculate covariance matrix for systematic error sources which are treated using covariance matrix approach
C
C> @param ScaledGamma
C> @param ScaledSystMatrix
C> @param List_Covar_Inv
C> @param n0_in
C
C---------------------------------------------------------------------------------------------
      subroutine chi2_calc_covar(ScaledGamma,ScaledSystMatrix, List_Covar_Inv, n0_in)

      implicit none
C
      include 'ntot.inc'
      include 'systematics.inc'
      include 'indata.inc'
      include 'theo.inc'
      double precision ScaledGamma(NSysMax,Ntot) ! Scaled Gamma matrix
      double precision ScaledSystMatrix(NCovarMax,NCovarMax)  ! syst. covar matrix
      integer List_Covar_Inv(NTOT)
C--
      integer n0_in
      integer i,j,k, i1, j1, i2, j2
C----------------------------------------------------------------------
      do i=1,n0_in
         do j=i,n0_in
            ScaledSystMatrix(i,j) = 0
         enddo
      enddo

      do k=1,nsys
         if (SysForm(k).eq.isMatrix) then            

! The covariance matrix:
            do i2=1,n_syst_meas(k)
               i1 = syst_meas_idx(i2,k)  ! data point index
               i = list_covar_inv(i1)    ! cov. matrix index 

               do j2=i2,n_syst_meas(k)
                  j1 = syst_meas_idx(j2,k)  ! data point idx
                  j = list_covar_inv(j1)    ! cov. matrix idx

                  ScaledSystMatrix(i,j) =  
     $                 ScaledSystMatrix(i,j) 
     $                 + ScaledGamma(k,i1)*ScaledGamma(k,j1)

               enddo
            enddo
         endif
      enddo

C----------------------------------------------------------------------
      end

      
C===========================================================
C
C> @brief Get uncertainties
C
C> @param[in] Idx data point index
C> @param[out] Stat absolute errors
C> @param[out] StatConst absolute errors
C> @param[out] Uncor absolute errors
C
C-----------------------------------------------------------
      Subroutine GetPointErrors(Idx, Stat, StatConst, Uncor)

      implicit none
      integer Idx
      double precision Stat, StatConst, Uncor
      include 'ntot.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'steering.inc'
      double precision d,t,mix
C-------------------------------------------------------------
      d = daten(idx)
      t = theo(idx)

      if ( t.le.0 ) then
         t = d
         call HF_errlog(13011601,
     $ 'W: Negative or zero prediction.'//
     $ ' Reset to data for error scaling.')
      endif

      mix = sqrt(abs(d*t))


      stat       = e_stat_poisson(idx)*mix
      statconst  = e_stat_const(idx)*d

      Uncor      = sqrt((e_uncor_mult(idx)*t)**2+
     $     (e_uncor_const(idx)*d)**2+(e_uncor_poisson(idx)*mix)**2)

C-------------------------------------------------------------
      end

C----------------------------------------------------------------------------
C
C> @brief Sum covariance matrices and invert the result
C
C> @param ScaledErrorMatrix
C> @param ScaledSystMatrix
C> @param ScaledTotMatrix
C> @param NCovar
C
C----------------------------------------------------------------------------
      Subroutine Chi2_calc_SumCovar(ScaledErrorMatrix, ScaledSystMatrix, 
     $              ScaledTotMatrix, NCovar)

      implicit none
      include 'ntot.inc'

      double precision ScaledErrorMatrix(NCovarMax,NCovarMax) ! stat+uncor error matrix
      double precision ScaledSystMatrix(NCovarMax,NCovarMax)  ! syst. covar matrix
      double precision ScaledTotMatrix(NCovarMax,NCovarMax)   ! stat+uncor+syst covar matrix
      integer NCovar
      double precision Array(NCovarMax*2)
      integer IFail

      integer i,j
C-----------------------------
      do i=1,NCovar
         do j=i,NCovar
            ScaledTotMatrix(i,j) = ScaledErrorMatrix(i,j) 
     $           + ScaledSystMatrix(i,j)
            ScaledTotMatrix(j,i) = ScaledTotMatrix(i,j)
         enddo
      enddo
      
        ! print *,' --- ScaledTotMatrix'
        ! do i=1,6
          ! print *,(ScaledTotMatrix(j,i),j=1,6)
        ! enddo

C-----------------------------   
      Call DInv(NCovar,ScaledTotMatrix,NCovarMax,Array,IFail)
C      print *,IFail,NCovar
   
      end

C-----------------------------------------------------------------------
C
C> @brief Scale covariance matrix and/or diagonal uncertainties
C 
C> @param ScaledErrors
C> @param ScaledErrorMatrix
C> @param rsys_in
C> @param n0_in
C> @param NCovar
C> @param List_Covar
C> @param Iterate
C
C-----------------------------------------------------------------------
      subroutine chi2_calc_stat_uncor(ScaledErrors, ScaledErrorMatrix, 
     $     rsys_in,n0_in, NCovar, List_Covar, Iterate)

      implicit none
      include 'ntot.inc'
      include 'systematics.inc'
      include 'steering.inc'
      include 'theo.inc'
      include 'covar.inc'

      double precision ScaledErrors(NTot), ScaledErrorMatrix(NCovarMax
     $     ,NCovarMax)
      double precision ScaledErrorsStat(NTot), ScaledErrorsSyst(NTot)
      double precision rsys_in(NSYS)
      integer n0_in, NCovar, List_Covar(NTot), iterate

      integer i,j,i1,j1
      double precision Stat, StatConst, Unc, Sum
c       
      
      include 'indata.inc'
      double precision Offs
C-------------------------------------------------------

C
C Start with diagonal part
C
      do i=1,n0_in
         Call GetPointErrors(i, Stat, StatConst, Unc)
         sum = 1.
         if (Chi2ExtraSystRescale .and. Iterate.eq.0) then
C Re-scale for systematic shifts:            
            do j=1,NSYS
               if ( (SysForm(j) .eq. isNuisance)
     $              .and. (SysScalingType(j) .eq. isLinear ) ) then
                  Sum = Sum - beta(j,i)*rsys_in(j)
               endif
            enddo
         endif
         Offs = 0.
         if (Chi2OffsFinal) then
            do j=1,NSYS
               if (SysForm(j) .eq. isOffset) then
                  Offs = Offs + beta(j,i)**2
               endif
            enddo
         endif
         ScaledErrors(i) = sqrt((Stat*Sum)**2+StatConst**2+Unc**2+Offs*daten(i)**2)
         ScaledErrorsStat(i) = sqrt((Stat*Sum)**2+StatConst**2)
         ScaledErrorsSyst(i) = sqrt(Unc**2+Offs*daten(i)**2)
      enddo

C
C Do also covariance part:
C
      do i1=1,NCovar
         i = List_Covar(i1)
         do j1=i1,NCovar
            j = List_Covar(j1) 

            ScaledErrorMatrix(i1,j1) = 
     $          ScaledErrorsStat(i)*ScaledErrorsStat(j)*corr_stat(i,j) +
     $          ScaledErrorsSyst(i)*ScaledErrorsSyst(j)*corr_syst(i,j) +
     $          cov(i,j)

         enddo
      enddo
C--------------------------------------------------------
      end

C----------------------------------------------------------------------------------
C
C> @brief Determine shifts of nuisance parameters
C
C> @param ScaledErrors
C> @param ScaledTotMatrix
C> @param ScaledGamma
C> @param rsys_in
C> @param ersys_in
C> @param list_covar_inv
C> @param iflag
C> @param n0_in
C
C----------------------------------------------------------------------------------
      subroutine chi2_calc_syst_shifts(
     $     ScaledErrors
     $     ,ScaledTotMatrix
     $     ,ScaledGamma
     $     ,rsys_in,ersys_in,list_covar_inv,  iflag, n0_in)

      implicit none
      include 'ntot.inc'
      include 'systematics.inc'
      include 'theo.inc'
      include 'indata.inc'
      include 'steering.inc'
C
      double precision ScaledErrors(NTOT)
      double precision ScaledTotMatrix(NCovarMax,NCovarMax)   !> stat+uncor+syst covar matrix
      double precision ScaledGamma(NSysMax,Ntot) !> Scaled Gamma matrix

      double precision rsys_in(NSYSMax), ERSYS_in(NSYSMax)
      integer list_covar_inv(NTOT),  iflag, n0_in
      logical doExternal

      integer k,l, i1,j1,i,j, j2, i2
      double precision A(NSYSMax,NSYSMax), C(NSysMax)
      double precision d_minus_t
      double precision ShiftExternal(NTOT)
      
      integer com_list(NTot),n_com_list  !> List of affected data, common for two sources.
      integer IR(2*NSysMax), Ifail

      logical lfirst
      data lfirst /.true./
      save lfirst
C-
      logical HaveCommonData(NsysMax, NsysMax)
C--------------------------------------------------------

C Determine pairs of syst. uncertainties which share  data
      if (LFirst .or. ResetCommonSyst) then
         LFirst = .false.
         ResetCommonSyst = .false. 
         do l=1,nsys
            do k=l,nsys
               Call Sys_Data_List12(l,k,n_com_list,com_list)
               if (n_com_list.gt.0) then
                  HaveCommonData(k,l) = .true.
               else
                  HaveCommonData(k,l) = .false.
               endif
            enddo
         enddo
      endif

C Get extra piece, from external systematics:
      do i=1,n0_in
         ShiftExternal(i) = 0.0D0
      enddo

      do l=1,nsys
         if (SysForm(l) .eq. isExternal ) then
            do i1 = 1, n_syst_meas(l)
               i  = syst_meas_idx(i1,l)
               ShiftExternal(i) = ShiftExternal(i) 
     $              + ScaledGamma(l,i)*rsys_in(l)
            enddo
         endif
      enddo

  ! A system of  "number  of isNuisance systematics" equations, indexed using "l":
  !
  !    A * Shift = C
  !


C Reset the matricies:
      do i=1,nsys
         C(i) = 0.0D0
         do j=1, nsys
            A(i,j) = 0.0D0
         enddo
         A(i,i)  = 1.D0
      enddo

      do l=1,nsys
         if ( SysForm(l) .eq. isNuisance ) then
C Start with "C"
            do i1=1,n_syst_meas(l)
               i = syst_meas_idx(i1,l)
               if (FitSample(i) ) then

                  if ( list_covar_inv(i) .eq. 0) then
C Diagonal error:
                     d_minus_t = daten(i) - theo(i) + ShiftExternal(i)
                     C(l) = C(l) + ScaledErrors(i)
     $                    *ScaledGamma(l,i)*( d_minus_t )
                  else
C Covariance matrix, need more complex sum:
                     i2 = list_covar_inv(i)
                     do j1=1,n_syst_meas(l)
                        j = syst_meas_idx(j1,l)

                        if (FitSample(j)) then
                           d_minus_t = daten(j) - theo(j) 
     $                          + ShiftExternal(i)
                           j2 = list_covar_inv(j) 
                           if (j2 .gt. 0) then
                              C(l) = C(l) + ScaledTotMatrix(i2,j2)
     $                             *ScaledGamma(l,i)*( d_minus_t )
                           endif
                        endif
                     enddo
                  endif
               endif
            enddo

C Now A:

            do k=l,NSys
C
               if ( (sysform(k) .eq. isNuisance ) .and.
     $              HaveCommonData(k,l) ) then

                  do i1 = 1,n_syst_meas(k)
                     i = syst_meas_idx(i1,k)
                     if ( FitSample(i) ) then
                        if (  list_covar_inv(i) .eq. 0) then
C Diagonal error:
                           A(k,l) = A(k,l) +
     $                          ScaledErrors(i)
     $                          *ScaledGamma(l,i)
     $                          *ScaledGamma(k,i)
                        else
C Covariance matrix:
                           i2 = list_covar_inv(i)
                           do j1=1,n_syst_meas(l)
                              j = syst_meas_idx(j1,l)
C                           do j=1,n0_in
                              if ( FitSample(j) ) then
                                 j2 = list_covar_inv(j)
                                 if (j2 .gt. 0) then
                                    A(k,l) = A(k,l) +
     $                                   ScaledTotMatrix(i2,j2)
     $                                   *ScaledGamma(l,i)
     $                                   *ScaledGamma(k,j)
                                 endif
                              endif
                           enddo

                        endif
                     endif
                  enddo

               endif               
            enddo
         endif
      enddo

C
C Under diagonal:
C
      do l=1,nsys
         do k=1,l-1
            A(k,l) = A(l,k)
         enddo
      enddo

C Ready to invert
      if (nsys.gt.0) then

         if (LDebug) then
            print *,'DUMP of Syst. shifts matrix'
            do l=1,nsys
               print *,'l=',l,C(l)
               do k=1,nsys
                  print *,l,k,A(l,k)
               enddo
            enddo
         endif

         
         if (iflag.eq.3) then
            Call DEQInv(Nsys,A,NsysMax,IR, IFail, 1, C)
         else
            Call DEQN(Nsys,A,NsysMax,IR,IFail,1,C)
         endif

         do l=1,nsys
            if ( Sysform(l) .eq. isNuisance) then
               rsys_in(l) = - C(l)
               if (iflag.eq.3) then
                  ersys_in(l) = sqrt(A(l,l))
               endif
            endif
         enddo
      endif
C--------------------------------------------------------
      end


C------------------------------------------------------------------------------------
C
C> @brief For input sources isys1 and isys2 build a compressed list of affected data points.
C> @param isys1
C> @param isys2
C> @param n_list
C> @param i_list
C
C------------------------------------------------------------------------------------
      subroutine Sys_Data_list12(isys1,isys2,n_list,i_list)

      implicit none
      include 'ntot.inc'
      include 'systematics.inc'
      integer isys1, isys2, i_list(NTot), n_list
      integer i1,j1, i,j
C---------------------------------------------------------------
      n_list = 0
      do i1=1,n_syst_meas(isys1)
         i = syst_meas_idx(i1,isys1)
         do j1=1,n_syst_meas(isys2) 
            j = syst_meas_idx(j1,isys2)
            if ( i.eq.j ) then
               n_list = n_list + 1 
               i_list(n_list) = i
               goto 17
            endif
         enddo
 17      continue
      enddo
C---------------------------------------------------------------

      end
      

C----------------------------------------------------------------------
C
C> @brief Calculate chi2.
C
C> @param ScaledErrors
C> @param ScaledGamma
C> @param ScaledTotMatrix
C> @param rsys_in
C> @param NDiag
C> @param List_Diag
C> @param NCovar
C> @param List_Covar
C> @param fchi2_in
C> @param pchi2_in
C> @param fcorchi2_in
C
C----------------------------------------------------------------------
      subroutine chi2_calc_chi2(ScaledErrors,ScaledGamma,
     $     ScaledTotMatrix,rsys_in
     $     ,NDiag, List_Diag, NCovar, List_Covar
     $     ,fchi2_in, pchi2_in, fcorchi2_in)

      implicit none
      include 'ntot.inc'
      include 'systematics.inc'
      include 'theo.inc'
      include 'indata.inc'

      double precision ScaledGamma(NSysMax,Ntot) ! Scaled Gamma matrix
      double precision ScaledErrors(Ntot)  ! uncorrelated uncertainties, diagonal
      double precision ScaledTotMatrix(NCovarMax,NCovarMax)   ! stat+uncor+syst covar matrix
      double precision rsys_in(NSysMax)
      integer NDiag, list_covar(NTot), NCovar, list_diag(NTot)
      double precision fchi2_in, pchi2_in(nset), fcorchi2_in

      integer i,j, i1, j1, k
      double precision d,t, chi2, sum
      
      double precision SumCov(NCovarMax)

C---------------------------------------------------------------------------
      fchi2_in = 0.0D0
 ! Also zero fit/control sample chi2s
      chi2_fit = 0.
      chi2_cont = 0.

C Diagonal part:
      do i1=1,NDiag
         i = list_diag(i1)
         d = DATEN(i)
         t = THEO (i)
         sum = 0.0D0
         do k = 1,NSys
            Sum = Sum + ScaledGamma(k,i)*rsys_in(k)
         enddo
C Chi2 per point:
         chi2 = (d - t + Sum)**2 * ScaledErrors(i)

C Sums:
         if ( FitSample(i) ) then
            chi2_fit  = chi2_fit  + chi2
            fchi2_in  = fchi2_in  + chi2
            pchi2_in(JSET(i)) = pchi2_in(JSET(i)) + chi2
         else
            chi2_cont = chi2_cont + chi2
         endif
      enddo

       ! print*,'chi2_calc1: ',fchi2_in
      
C Covariance matrix part

C 1) Pre-compute sums of systematic shifts:
      do i1=1,NCovar
         i = list_covar(i1) 
         sumcov(i1) = Daten(i) - Theo(i)
         do k = 1,NSys
            SumCov(i1) = SumCov(i1) + ScaledGamma(k,i)*rsys_in(k)
         enddo        
      enddo


C 2) Actual chi2 calculation:
      do i1=1,NCovar
         i = list_covar(i1)
         Chi2 = 0

         do j1 = 1, NCovar
            j = list_covar(j1)
            Chi2 = Chi2 + SumCov(i1)*SumCov(j1)*ScaledTotMatrix(i1,j1)
         enddo
C Sums:
         if ( FitSample(i) ) then
            chi2_fit  = chi2_fit  + chi2
            fchi2_in  = fchi2_in  + chi2
            pchi2_in(JSET(i)) = pchi2_in(JSET(i)) + chi2
         else
            chi2_cont = chi2_cont + chi2
         endif         

      enddo

       ! print*,'chi2_calc2: ',fchi2_in

C Correlated chi2 part:
      fcorchi2_in = 0.d0
      do k=1,NSys
         fcorchi2_in = fcorchi2_in + rsys_in(k)**2
      enddo
      fchi2_in = fchi2_in + fcorchi2_in

       ! print*,'chi2_calc3: ',fchi2_in

C---------------------------------------------------------------------------
      end

C--------------------------------------------------------------------------
C
C> @brief Calculate additional log correction factor
C> @param ScaledErrors uncertainties
C> @param chi2_log log correction
C> @param n0_in number of points
C
C--------------------------------------------------------------------------
      subroutine chi2_calc_PoissonCorr(ScaledErrors, chi2_log, n0_in)

      implicit none
      include 'ntot.inc'
      double precision ScaledErrors(Ntot)
      double precision chi2_log
      integer n0_in

      include 'indata.inc'
      include 'systematics.inc'
      integer i
C-------------------------------------------------------------------------
      chi2_log = 0.D0
      do i=1,n0_in
         if (FitSample(i)) then
            chi2_log = chi2_log - log( alpha(i)*alpha(i) 
     $           * ScaledErrors(i))
         endif
      enddo

      end


C------------------------------------------------------------------------
C
C> @brief Extra output for iteration #3
C> @param ScaledErrors uncorrelated uncertainties, diagonal
C> @param ScaledGamma scaled Gamma matrix
C> @param RSys_in
C> @param n0_in
C
C------------------------------------------------------------------------
      Subroutine Chi2_calc_FCN3(ScaledErrors,ScaledGamma,RSys_in, n0_in)

      implicit none
      include 'ntot.inc'
      include 'systematics.inc'
      include 'theo.inc'
      double precision ScaledErrors(Ntot)  ! uncorrelated uncertainties, diagonal
      double precision ScaledGamma(NSysMax,Ntot) ! Scaled Gamma matrix
      double precision rsys_in(NSysMax)
      integer n0_in

      integer i,k
C------------------------------------------------------------------
      do i=1,n0_in
         ALPHA_MOD(i) =  1.D0/sqrt(ScaledErrors(i))
         THEO_MOD(i)  = THEO(i)
         do k=1,NSYS
            THEO_MOD(i) = THEO_MOD(i) - ScaledGamma(k,i)*RSys_in(k)
         enddo
      enddo
      
      ! CALL cvfillgamma(nsys,n0_in,ScaledGamma,NSYSMAX)
      ! CALL cvfillgamma(ScaledGamma,nsys,n0_in,NTOT)

C------------------------------------------------------------------
      end

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      subroutine GetChisquare(flag_in,n0_in,fchi2_in,rsys_in,ersys_in,pchi2_in,fcorchi2_in)
      implicit none

      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'for_debug.inc'
      include 'systematics.inc'
      INCLUDE 'theo.inc'


      integer ir(nsysmax),n0_in 
      integer isys,ipoint,jpoint,ifail,flag_in
      double precision chisq,fchi2_in,chi2error,fchi2_error
      double precision d,t,error,errorunc,errorconst
      double precision errorsta, fac, facd, fcorchi2_in
      double precision d_i, t_i, error_i, d_j, error_j, t_j         
      integer i,j,jsys,h1iset
      double precision pchi2_in(nset)
      double precision EBSYS_in(NSYSMax),ERSYS_in(NSYSMax)
      double precision BSYS_in(NSYSMax),RSYS_in(NSYSMax)
      double precision sub

      double precision dNEvt, tNEvt

      double precision factor2, factor_1, factor_2

      integer npoisson, ngauss

      integer k,l
C      integer getcachesize

*     ----------------------------------------------------------
*     Initialise
*     ----------------------------------------------------------
      npoisson = 0
      ngauss   = 0

      chisq=0.d0

      fchi2_error = 0.d0

      fchi2_in = 0.d0

 ! Also zero fit/control sample chi2s
      chi2_fit = 0.
      chi2_cont = 0.


      sub = 0.d0


      do jsys=1,nsys
         ir(nsys)=0.d0
         bsys_in(jsys) = 0.d0
         rsys_in(jsys) = 0.d0
         ebsys_in(jsys) = 0.d0
         ersys_in(jsys) = 0.d0
      enddo

      do i=1,nset
         pchi2_in(i)=0.d0
      enddo
*     -- for ICHI2=2, the matrix sysa is filled already in Systematics

      if (ICHI2.ne.2) then
         do isys=1,nsys
            do jsys=1,nsys
               sysa(isys,jsys) = 0.d0
               if (jsys.eq.isys) sysa(isys,jsys) = 1.d0
            enddo
         enddo
      endif


*     ----------------------------------------------------------
*     Pascaud-like chi2 plus error scaled modifications
*     ----------------------------------------------------------
      if (mod(ICHI2, 10).eq.1) then



*     ---------------------------------------------------------
*     now calculate the bsys(nsys) and the sysa matrix
*     ---------------------------------------------------------

            
         do ipoint=1,n0_in
               
            d = DATEN(ipoint)
            t = THEO(ipoint)
            error = ALPHA(ipoint)
                  

***   scale errors for chi2 calculation
***   in principle:  unc*(t/d),  sta*dsqrt(t/d)
            if (ICHI2.eq.11 .or. ICHI2.eq.41) then
               if ( (alpha(ipoint)/d.gt.Chi2MaxError)
     $              .or. (ControlFitSplit.and. .not. FitSample(ipoint))
     $              ) then
C     Turn off the point for the syst. errors shift estimation:
                  error = 1.D010
               else
***   mixed scaling - decompose - scale - recombine
                  errorunc = E_UNC(ipoint)*d/100.
                  errorconst = E_STA_CONST(ipoint)*d/100.

                  if (errorunc.gt.error) then
                     errorsta = 0.
                  else
                     errorsta = error**2-errorunc**2-errorconst**2
                     if (errorsta.gt.0) then
                        errorsta = sqrt(errorsta)
                     else
                        errorsta = 0.
                     endif
                  endif
                  if (t.gt.0) then
                     errorsta = errorsta*dsqrt(abs(t/d))
                     errorunc = errorunc*(abs(t/d))
                  endif
                  error = dsqrt(errorsta**2+errorunc**2+errorconst**2)
               endif

            else if (ICHI2.eq.21) then
***   linear scaling
               error = error*(abs(t/d))
            else if (ICHI2.eq.31) then
***   sqrt scaling
               error = error*dsqrt(abs(t/d))
            endif

            

            factor2 = 1.D0/error**2

            do isys = 1,nsys                                 

               if (beta(isys,ipoint).ne.0) then
                  if (SysScalingType(isys) .eq. isNoRescale) then
                     factor_1 = -d
                  else
                     factor_1 = t
                  endif

                  bsys_in(isys) = bsys_in(isys) 
     +                 + factor_1*(d-t)*BETA(isys,ipoint)*factor2
               

c                  ebsys_in(isys) = ebsys_in(isys)     !> ????
c     +                 + t * BETA(isys,ipoint)/error
               
                  do  jsys=isys,nsys
                     if (beta(jsys,ipoint).ne.0) then
                        if (SysScalingType(jsys) .eq. isNoRescale ) then
                           factor_2 = -d
                        else
                           factor_2 = t
                        endif

                        sysa(isys,jsys) = sysa(isys,jsys)
     +                       + beta(isys,ipoint)*beta(jsys,ipoint)
     $                       *factor_1*factor_2*factor2

                     endif
                  enddo
               endif
            enddo

         enddo


         do isys=1,nsys
C            print '(5F10.2)',(sysa(isys,jsys),jsys=1,5)
            do jsys=isys+1,nsys
               sysa(jsys,isys) = sysa(isys,jsys)
            enddo
         enddo
*     ---------------------------------------------------------
*     inverse sysa and find the shifts
*     ---------------------------------------------------------


         
         if (nsys.gt.0) then

            if (flag_in.eq.3) then
               Call DEQINV(NSys,sysa,NSYSMAX,IR,IFAIL,1,bsys_in)
            else
c               CALL DPOSV('Upper',Nsys,1,sysa,NSYSMax,Bsys_in,NSYSmax
c     $              ,Ifail)
               Call DEQN(NSys,sysa,NSYSMAX,IR,IFAIL,1,bsys_in)
            endif
         endif

         do isys=1,nsys
            rsys_in(isys) = -bsys_in(isys)
         enddo

         if (DEBUG.and.flag_in.eq.1) then
            write(78,*)
            do isys=1,nsys
               write(78,*) 'isys rsys ',isys,rsys_in(isys)
            enddo
            write(78,*)
         endif




*     ---------------------------------------------------------
*     now calculate the chi2
*     ---------------------------------------------------------

         do ipoint=1,n0_in

            h1iset = JSET(ipoint)
            d = daten(ipoint)
            t = theo(ipoint)

*** Factor to scale theory/data:
            fac  = 1.d0
            facd = 1.d0
            do isys=1,nsys
               if ( SysScalingType(isys).eq. isLinear ) then
                  fac  = fac  - rsys_in(isys)*beta(isys,ipoint)                
               else
                  facd = facd + rsys_in(isys)*beta(isys,ipoint)                
               endif
            enddo 


            error = alpha(ipoint)

***   scale errors for chi2 calculation - as above!
            if (ICHI2.eq.11 .or. ICHI2.eq.41) then
***   mixed scaling - decompose - scale - recombine
               errorunc = E_UNC(ipoint)*d/100.
               errorconst = E_STA_CONST(ipoint)*d/100.
               if (errorunc.gt.error) then
                  errorsta = 0.
               else
                  errorsta = error**2-errorunc**2-errorconst**2
                  if (errorsta.gt.0) then
                     errorsta = sqrt(errorsta)
                  else
                     errorsta = 0.
                  endif
               endif
               if (t.gt.0) then
                  if (iDH_MOD.ne.0) then
                     errorsta = errorsta*dsqrt(abs(t*fac/d))
                  else
                     errorsta = errorsta*dsqrt(abs(t/d))
                  endif
                  errorunc = errorunc*(abs(t/d))
               endif
               error = dsqrt(errorsta**2+errorunc**2+errorconst**2)
               ! Extra contribution due to 2xlog sigma term:
               chi2error =  2.*log( error/alpha(ipoint)) ! subtract un-modified error such that delta chi2=0 if errors are not modified.


            else if (ICHI2.eq.21) then
***   linear scaling
               error = error*(abs(t/d))
            else if (ICHI2.eq.31) then
***   sqrt scaling
               error = error*dsqrt(abs(t/d))
            endif

                                  
            if (DEBUG.and.flag_in.eq.1) then
               write(78,*) 'ipoint fac ',ipoint,fac
            endif


            t = t*fac
            d = d*facd

            THEO_MOD(ipoint)=t
            ALPHA_MOD(ipoint)=error

C
C Add Poisson log likelihood
C
            if (alpha(ipoint)/d.gt.Chi2MaxError) then
C use Poisson formula:
               dNEvt = (d*d)/(alpha(ipoint)*alpha(ipoint))
               tNEvt = (t*d)/(alpha(ipoint)*alpha(ipoint))
               chisq = 2* ( tNEvt - dNEvt - dNEvt*log(tNevt/dNEvt))
               
C               print *,dnevt,tnevt,chisq
               chi2error = 0.
               npoisson = npoisson+1
            else
               chisq = (d-t)**2/error**2
               ngauss = ngauss+1
            endif

            if (ControlFitSplit) then
               if (FitSample(ipoint)) then
                  chi2_fit  = chi2_fit + chisq
                  fchi2_in  = fchi2_in + chisq
               else
                  chi2_cont = chi2_cont + chisq
               endif
      
            else               
               fchi2_in = fchi2_in + chisq
            endif

            if ( ICHI2.eq.41) then
               fchi2_error = fchi2_error + chi2error
               fchi2_in = fchi2_in + chi2error
            endif

            if (flag_in.eq.3) then
               pchi2_in(h1iset) = pchi2_in(h1iset)+chisq
            endif

            if (flag_in.eq.1) then
               write(87,880) h1iset, 
     $              AbstractBins(1,ipoint),
     $              AbstractBins(2,ipoint),AbstractBins(3,ipoint),
     $              t,d,(d-t)/error
            endif
            if (flag_in.eq.3) then
               write(88,880) h1iset,
     $              AbstractBins(1,ipoint),
     $              AbstractBins(2,ipoint),AbstractBins(3,ipoint),
     $              t,d,(d-t)/error
 880           format(1x, i2, 2x, G12.6, 2x, G12.4, 2x, G12.6, 3(2x, G12.4))
            endif

*     -- errors on the shifts 

            if (flag_in.eq.3) then
               do isys=1,nsys
                  ersys_in(isys) = sqrt(sysa(isys,isys))
               enddo
            endif

         enddo

         do isys=1,nsys
            fchi2_in = fchi2_in + rsys_in(isys)**2
         enddo

         

c....print out the correlated chi2
         if (flag_in.eq.3) then
            fcorchi2_in=0.0d0
            do isys=1,nsys
               fcorchi2_in= fcorchi2_in +rsys_in(isys)**2
            enddo
         endif


c...........................


      elseif (ICHI2.eq.2) then  !  CTEQ-like chi2

         do ipoint=1,n0_in

            d_i = daten(ipoint)
            t_i = theo(ipoint)
            error_i = alpha(ipoint)

            chisq = (d_i-t_i)**2 / error_i**2
            fchi2_in = fchi2_in + chisq

            do jsys=1,nsys
               bsys_in(jsys) = bsys_in(jsys)
     +              + beta(jsys,ipoint)*(d_i-t_i)/error_i**2
            enddo

            if (flag_in.eq.3) then
               h1iset = JSET(ipoint)
               pchi2_in(h1iset) = pchi2_in(h1iset)+chisq
            endif


         enddo

*     
*     - Now calculate the term to subtract from chi2 (cf CTEQ)

         do i=1,nsys
            do j=1,nsys
               sub = sub + bsys_in(i) * sysa(i,j) * bsys_in(j)
            enddo
         enddo


c        print *,'haha',fchi2_in,sub
         fchi2_In = fchi2_in - sub


*     -- and get the systematic shifts :

         do i=1,nsys
            RSYS_in(i) = 0.
            do j=1,nsys
               rsys_in(i) = rsys_in(i) + sysa(i,j) * bsys_in(j)
            enddo
         enddo

c....print out the correlated chi2
         if (flag_in.eq.3) then
            fcorchi2_in=0.0d0
            do isys=1,nsys
               fcorchi2_in= fcorchi2_in +rsys_in(isys)**2
            enddo
         endif
c...........................

      elseif (mod(ICHI2, 10).eq.3) then ! Offset
*     ---------------------------------------------------------
*     now calculate the chi2
*     ---------------------------------------------------------

         fcorchi2_in=0.0d0
         
         do ipoint=1,n0_in

            h1iset = JSET(ipoint)
            d = daten(ipoint)
            t = theo(ipoint)
            error = alpha(ipoint)

            if (DEBUG.and.flag_in.eq.1) then
               write(78,*) 'ipoint fac ',ipoint,fac
            endif


            THEO_MOD(ipoint)=t
            ALPHA_MOD(ipoint)=error

C
C Add Poisson log likelihood
C
            if (alpha(ipoint)/d.gt.Chi2MaxError) then
C use Poisson formula:
               dNEvt = (d*d)/(alpha(ipoint)*alpha(ipoint))
               tNEvt = (t*d)/(alpha(ipoint)*alpha(ipoint))
               chisq = 2* ( tNEvt - dNEvt - dNEvt*log(tNevt/dNEvt))
               
C               print *,dnevt,tnevt,chisq
               chi2error = 0.
               npoisson = npoisson+1
            else
               chisq = (d-t)**2/error**2
               ngauss = ngauss+1
            endif
            fchi2_in = fchi2_in + chisq
            if (flag_in.eq.3) then
              fcorchi2_in = fcorchi2_in + ((1-t/d)/E_TOT(ipoint)*1.d2)**2
            endif

            if (flag_in.eq.3) then
               pchi2_in(h1iset) = pchi2_in(h1iset)+chisq
            endif

            if (flag_in.eq.1) then
               write(87,880) h1iset, 
     $              AbstractBins(1,ipoint),
     $              AbstractBins(2,ipoint),AbstractBins(3,ipoint),
     $              t,d,(d-t)/error
            endif
            if (flag_in.eq.3) then
               write(88,880) h1iset,
     $              AbstractBins(1,ipoint),
     $              AbstractBins(2,ipoint),AbstractBins(3,ipoint),
     $              t,d,(d-t)/error
cws     880           format(1x, i2, 2x, G12.6, 2x, G12.4, 2x, G12.6, 3(2x, G12.4))
            endif

         enddo

      endif

      if (ichi2.eq.41) then
         print '(''Chi2 due to 2xlog sigma term'',F6.1)', fchi2_error
      endif

      if (npoisson.gt.0) then
         print '(''Use Poisson stats for '',I6,'', use Gauss stats for''
     $        ,I6,'' events'')',npoisson,ngauss
      endif

      return
      end

      
C---------------------------------------------------------------------
C
C> @brief 
C> @param ScaledGamma scaled Gamma matrix
C> @param ScaledGammaSav scaled Gamma matrix
C> @param ScaledOmega scaled Omega matrix
C> @param rsys_in
C> @param Iteration iteration index
C> @param LStop
C
C---------------------------------------------------------------------
      subroutine UseOmegaScale(ScaledGamma,ScaledGammaSav,ScaledOmega,
     $     rsys_in,Iteration,LStop) 
      
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'systematics.inc'
      double precision ScaledGamma(NSysMax,Ntot) ! Scaled Gamma matrix
      double precision ScaledGammaSav(NSysMax,Ntot) ! Scaled Gamma matrix
      double precision ScaledOmega(NSysMax,Ntot) ! Scaled Omega matrix
      double precision  RSYS_in(NSYSMax)
      double precision rsys_save(NSYSMax)
      integer Iteration
      logical LStop
      integer i,j,k, iter 
      double precision shift

      integer iterMax
      parameter (iterMax = 10)

C----------------------------------------------------------
      if ( AsymErrorsIterations.eq.0) then
         LStop = .true.
         Return
      endif

      if (Iteration.eq.1) then
         do i=1,nsys
            rsys_save(i) = 0.
         enddo
      endif

C calculate shift in rsys:
      shift = 0.
      do i=1,nsys
         shift = max(shift, abs(rsys_in(i)-rsys_save(i)))
      enddo

      do i=1,nsys
         rsys_save(i) = rsys_in(i)
      enddo

C recalulate
      do k=1,nsys
         do i=1,n_syst_meas(k)
            j =  syst_meas_idx(i,k)
            ScaledGamma(k,j) = ScaledGammaSav(k,j) 
     $           + ScaledOmega(k,j)*rsys_in(k)
         enddo
      enddo            
      
      if (LDebug) then
         print *,'shift',iteration, shift,nsys
      endif

      LStop = .false.
      if (iteration.ge. AsymErrorsIterations ) then

C ! Check if max. shift is small:
         if ( shift.gt. 0.05) then
            print *,'ERROR: Large nuisance parameter shift =',shift
            print *,'CONSIDER INCREASING AsymErrorsIterations to'
     $           ,AsymErrorsIterations+5
            call HF_errlog(13053001,
     $ 'W:UseOmegaScale: Max shift exeeds 5%. Consider increasing'
     $           //' AsymErrorsIterations')
         endif

         LStop = .true.
      endif

      end
