C--------------------------------------------------------
C> @Brief Calculate chisquare
C
C>      - first get error matrix
C>      - invert error matric to get errors
C>      - calculate chisquare
C
C
C---------------------------------------------------------

      subroutine GetNewChisquare(flag_in,n0_in,fchi2_in,rsys_in,ersys_in,pchi2_in,
     $     fcorchi2_in)
      implicit none

#include "ntot.inc"
#include "steering.inc"
#include "systematics.inc"
#include "indata.inc"
#include "datasets.inc"

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
c Global initialisation

      if (LFirst) then
         LFirst = .false.

C    !> Determine which mechanisms for syst. errors should be used:
         Call Init_Chi2_calc(doMatrix, doNuisance, doExternal)

C    !> Determine which errors are diagonal and which are using covariance matrix
         Call init_chi2_stat(NDiag, NCovar, List_Diag, List_Covar,
     $        List_Covar_inv,n0_in)

         print*,'NDiag=',NDiag,'NCovar=',NCovar


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
C !> Rebuild syst. covariance matrix

         if ( doMatrix ) then
            Call Chi2_calc_covar(ScaledGamma
     $           ,ScaledSystMatrix
     $           ,List_Covar_Inv,n0_in)
         endif
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

C  !> Sum covariance matricies and invert the total:

            if ( doMatrix .or. NCovar .gt. 0 ) then
               Call Chi2_calc_SumCovar(ScaledErrorMatrix,
     $              ScaledSystMatrix,
     $              ScaledTotMatrix, NCovar)
            endif

C !> same for diagonal part:
            do i=1,n0_in
               if(ScaledErrors(i)/=0d0)then
                 ScaledErrors(i)=ScaledErrors(i)**(-2)
               else
                 ScaledErrors(i)=1d0 !When is this necessary? --Ivan
                 if(NCovar==0)then
                   !no cov matrix and no ScaledErrors errors, break
                   print*,'GetNewChisquare: no stat and unc errors in
     $ data! (possibly cov matrix forgot to be included?)'
                   call hf_stop
                 endif
               endif
            enddo
         endif

C !> Next determine nuisance parameter shifts
         omegaIteration = 1
         do
            if ( LConvertCovToNui .and. do_reduce
     $           .and. flag_in .ne. 3 ) then
                  ! use simplified (slightly) faster version of the code
               call chi2_calc_syst_shifts_simple(
     $              ScaledErrors
     $              ,ScaledGamma
     $              ,rsys_in
     $              ,n0_in
     $              )
c               stop
            else
               Call Chi2_calc_syst_shifts(
     $              ScaledErrors
     $              ,ScaledTotMatrix
     $              ,ScaledGamma
     $              ,rsys_in,ersys_in,list_covar_inv, flag_in, n0_in
     $              ,scaledOmega)
            endif

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

C !> For asymmetric erros and exteral systematic sources we need to modify ScaledGamma:
      if (doExternal.and.AsymErrorsIterations.gt.0) then
         call chi2_calc_asymError_external(ScaledGamma, ScaledOmega, rsys_in)
      endif

C !> Calculate chi2
      call chi2_calc_chi2(
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
         if (lDebug) print '(''Log term contribution='',F6.2)',chi2_log
      else
         do i=1,NDATASETS
            chi2_poi(i) = 0.D0
         enddo
         chi2_poi_tot = 0.D0
      endif
       ! print*,'fchi2_in=',fchi2_in


C
C !> Store extra output for FCN = 3:
C
      if (Flag_In.eq.3) then
         Call Chi2_calc_FCN3(ScaledErrors,ScaledGamma,RSys_in,n0_in)
      endif

c     export systematic uncertainties and shifts
      do k=1,nsys
         do j=1,n_syst_meas(k)
            i = syst_meas_idx(j,k)

            scgamma(k,i) = ScaledGammaSav(k,i)
            if (AsymErrorsIterations.eq.0) then
               scomega(k,i) = 0d0
            else
               scomega(k,i) = ScaledOmega(k,i)
            endif
            sysshift(k) = rsys_in(k)
         enddo
      enddo

c     export uncorrelated errors
      do i=1,n0_in
         scerrors(i) = 1./sqrt(ScaledErrors(i))
c         print *,i,1./sqrt(ScaledErrors(i))
      enddo


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
#include "ntot.inc"
#include "systematics.inc"
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
#include "ntot.inc"
#include "indata.inc"
#include "systematics.inc"
#include "steering.inc"
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
#include "ntot.inc"
#include "systematics.inc"
#include "endmini.inc"
#include "extrapars.inc"
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

#include "ntot.inc"
#include "systematics.inc"

      double precision ScaledGamma(NSysMax,Ntot) !> Scaled Gamma matrix
      double precision ScaledOmega(NSysMax,Ntot) !> Scaled Omega matrix

#include "indata.inc"
#include "theo.inc"
c#include "steering.inc"

      integer i1,j1,i,j,k
      integer scaling_type
      double precision scale
      logical lfirstPass/.true./
C-----------------------------------------------------
      do k=1,NSYS
         scaling_type = SysScalingType(k)

         do i1=1,n_syst_meas(k)
            i = syst_meas_idx(i1,k)

            if ( (scaling_type.eq. isNoRescale)
     $           .or. LForceAdditiveData(i) ! force additive scaling for all syst. sources for this datapoint
     $           ) then
               scale = daten(i)
               if ( LForceAdditiveData(i) .and. k.eq.1
     $              .and. lfirstPass) then
                  call hf_errlog(2018121101,
     $                 'I: Force additive systematics for a datapoint')
               endif
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
      lfirstPass = .false.
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
#include "ntot.inc"
#include "systematics.inc"
#include "indata.inc"
#include "theo.inc"
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
C> @brief      Get uncertainties
C> @param[in]  Idx data point index
C> @param[out] Stat absolute errors
C> @param[out] StatConst absolute errors
C> @param[out] Uncor absolute errors
C
C-----------------------------------------------------------
      Subroutine GetPointErrors(Idx, Stat, StatConst, Uncor)

      implicit none
      integer Idx
      double precision Stat, StatConst, Uncor
#include "ntot.inc"
#include "indata.inc"
#include "theo.inc"
#include "steering.inc"
      double precision d,t,mix
C-------------------------------------------------------------
      d = daten(idx)
      t = theo(idx)

      if ( t.le.0 ) then
         t = d
         if (e_stat_poisson(idx).ne.0) then ! warning only if scaling is needed/
            call HF_errlog(13011601,
     $           'W: Negative or zero prediction.'//
     $           ' Reset to data for error scaling.')
         endif
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
#include "ntot.inc"

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
C> @param[out] ScaledErrors
C> @param[out] ScaledErrorMatrix
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
#include "ntot.inc"
#include "systematics.inc"
#include "steering.inc"
#include "theo.inc"
#include "covar.inc"

      double precision ScaledErrors(NTot), ScaledErrorMatrix(NCovarMax
     $     ,NCovarMax)
      double precision ScaledErrorsStat(NTot), ScaledErrorsSyst(NTot)
      double precision rsys_in(NSYS)
      integer n0_in, NCovar, List_Covar(NTot), iterate

      integer i,j,i1,j1
      double precision Stat, StatConst, Unc, Sum
c

#include "indata.inc"
      double precision Offs
C-------------------------------------------------------

C
C Start with diagonal part
C
      do i=1,n0_in
         Call GetPointErrors(i, Stat, StatConst, Unc)
         sum=0.
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

         sum=exp(sum)
         ScaledErrorsStat(i)=sqrt(Stat**2*sum+StatConst**2)
         ScaledErrorsSyst(i)=sqrt(Unc**2+Offs*daten(i)**2)
         ScaledErrors(i)=sqrt(ScaledErrorsStat(i)**2+ScaledErrorsSyst(i)**2)
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
     $          ScaledErrors(i)*ScaledErrors(j)*corr(i,j) +
     $          cov(i,j)
         enddo
      enddo
C--------------------------------------------------------
      end

C-----
C
C> @brief extend lists of data-syst sources for data points connected via cov. matrix
C
      subroutine expand_syst_lists(tot_matrix,list_covar_inv,n0_in)
      implicit none
#include "ntot.inc"
#include "systematics.inc"
      double precision tot_matrix(NCovarMax,NCovarMax)   !> stat+uncor+syst covar matrix
      integer list_covar_inv(NTOT)
      integer n0_in
      integer l,j,j1,i,ic,n,jc,k
      logical flag
C-------------------------------------------------
      do l=1,nsys
         n = n_syst_meas(l)
         j1 = 1

         do while ( j1 .le. n)
            flag = .false.
C loop over all data, find non-zero correlations
            j = syst_meas_idx(j1,l)     ! j -> index of the data
            jc = list_covar_inv(j)
            if (jc.gt.0) then
C Check all data which may have correlations with this data point
               do i=1,n0_in
                  ic = list_covar_inv(i)
                  if (ic .gt.0 ) then
                     if (tot_matrix(jc,ic).ne.0.0) then
C Check if "i" point is on the list already
                        do k=1,n
                           if ( i.eq.syst_meas_idx(k,l)) then
                              goto 17
                           endif
                        enddo
C New point, add to the lists:
                        n = n + 1
                        syst_meas_idx(n,l) = i
c                        print *,'EXPAND LIST',l,i
                        flag = .true.
 17                     Continue
                     endif
                  endif
               enddo
            endif
            if (.not.flag) then
               j1 = j1 + 1
            endif
         enddo

         n_syst_meas(l) = n

      enddo
C-------------------------------------------------
      end



      subroutine chi2_calc_syst_shifts_simple(
     $     ScaledErrors
     $     ,ScaledGamma
     $     ,rsys_in,   n0_in)

      implicit none
#include "ntot.inc"
#include "systematics.inc"
#include "theo.inc"
#include "indata.inc"
#include "steering.inc"
C
      double precision ScaledErrors(NTOT)
      double precision ScaledGamma(NSysMax,Ntot) !> Scaled Gamma matrix
      double precision rsys_in(NSYSMax)
      double precision A(NSYSMax,NSYSMax), C(NSysMax)

      double precision AS(n0_in,NSysMax)  ! automatic, scaled sys.

      double precision ASp(n0_in*(NsysMax+1)/2),
     $     SGp(n0_in*(NsysMax+1)/2)

      double precision d_minus_t1
      integer   n0_in
      integer i,j,l,i1,k
      integer ifail
      integer IR(2*NSysMax)

      double precision time1, time2

C--------------------------------------------------------------------------------------------
C Reset the matricies:
      do i=1,nsys
         C(i) = 0.0D0
         do j=1, nsys
            A(i,j) = 0.0D0
         enddo
      enddo

      do i=1,n0_in
         do j=1,nsys
            AS(i,j) =
     $           ScaledErrors(i)
     $           *ScaledGamma(j,i)
         enddo
      enddo


      do l=1,nsys
C Start with "C"

         do i1=1,n_syst_meas(l) ! loop over all data affected by this source
            i = syst_meas_idx(i1,l) ! i -> index of the data
c         do i=1,n0_in

            d_minus_t1 = daten(i) - theo(i)

C  Diagonal error:
            C(l) = C(l) +  AS(i,l)
     $           *( d_minus_t1 )

         enddo
      enddo

      call cpu_time(time1)


      if ( .not. UseBlas ) then

!$OMP PARALLEL DO

         do i=1,n0_in
            do l=1,nsys
               do k=l,NSys
c Diagonal error:
                  A(k,l) = A(k,l) +
     $                 AS(i,l)
     $                 *ScaledGamma(k,i)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO


      else
C symmetric matrix does not work
c      print *,A(1,1),A(nsys,nsys)
Cuse BLAS: L L; R L; R U; L U
c      call cublas_dsymm('L','U',nsys,n0_in, 1.0D0, ScaledGamma, n0_in, AS
c      call dsymm('L','U',nsys,n0_in, 1.0D0, ScaledGamma, n0_in, AS
C use BLAS: L L; R L; R U; L U
c      call dsymm('R','U',nsys,nsys, 1.0D0, AS, n0_in, ScaledGamma
c     $     , nsysmax
c     $     , 0.D0, A, nsysmax)

         call dgemm('N','N',nsys,nsys, n0_in, 1.0D0
C         call cublas_dgemm('N','N',nsys,nsys, n0_in, 1.0D0
     $     , ScaledGamma
     $        , nsysmax
     $     , AS
     $     , n0_in
     $     , 0.D0, A, nsysmax)

      endif

C Penalty term, unity by default
      do i=1,nsys
         A(i,i) = A(i,i) + SysPriorScale(i)
      enddo

c      print *,A(1,1),A(nsys,nsys)

      call cpu_time(time2)
      print *,'CPU LOOP=',time2-time1
      call flush
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

         Call DEQN(Nsys,A,NsysMax,IR,IFail,1,C)

         do l=1,nsys
            rsys_in(l) = - C(l)
         enddo
      endif

      call cpu_time(time1)
C--------------------------------------------------------------------------------------------
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
C> @param ScaledOmega
C
C----------------------------------------------------------------------------------
      subroutine chi2_calc_syst_shifts(
     $     ScaledErrors
     $     ,ScaledTotMatrix
     $     ,ScaledGamma
     $     ,rsys_in,ersys_in,list_covar_inv,  iflag, n0_in, ScaledOmega)

      implicit none
#include "ntot.inc"
#include "systematics.inc"
#include "theo.inc"
#include "indata.inc"
#include "steering.inc"
C
      double precision ScaledErrors(NTOT)
      double precision ScaledTotMatrix(NCovarMax,NCovarMax)   !> stat+uncor+syst covar matrix
      double precision ScaledGamma(NSysMax,Ntot) !> Scaled Gamma matrix
      double precision ScaledOmega(NSysMax,Ntot) ! Scaled Omega matrix

      double precision rsys_in(NSYSMax), ERSYS_in(NSYSMax)
      integer list_covar_inv(NTOT),  iflag, n0_in
      logical doExternal

      integer k,l, i1,j1,i,j, j2, i2
      double precision A(NSYSMax,NSYSMax), C(NSysMax)

      double precision, allocatable :: AA(:,:)
      double precision, allocatable :: AA2(:,:)
      double precision, allocatable :: RR(:,:)

      double precision d_minus_t1, d_minus_t2,add
      double precision ShiftExternal(NTOT)

      integer com_list(NTot),n_com_list  !> List of affected data, common for two sources.
      integer IR(2*NSysMax), Ifail,  Npdf

      integer nsystheo, itheoisys(NSysMax)
      integer nsys_sav, n0_in_sav

      logical lfirst
      data lfirst /.true./
      data nsys_sav,n0_in_sav/0,0/
      save lfirst,nsys_sav,n0_in_sav
C-
      logical HaveCommonData(NsysMax, NsysMax)
C--------------------------------------------------------
C Check if number of sources/data points change:
      ResetCommonSyst = (nsys.ne.nsys_sav) .or. (n0_in.ne.n0_in_sav)
      nsys_sav = nsys
      n0_in_sav = n0_in 
C Determine pairs of syst. uncertainties which share  data


      if (LFirst .or. ResetCommonSyst) then
         LFirst = .false.
         ResetCommonSyst = .false.


         call expand_syst_lists(scaledtotmatrix,list_covar_inv,n0_in)

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
  ! Consider asymmetric uncertainties:
               if (AsymErrorsIterations.eq.0) then
                  ShiftExternal(i) = ShiftExternal(i)
     $                 + ScaledGamma(l,i)*rsys_in(l)
               else
                  ShiftExternal(i) = ShiftExternal(i)
     $                 + ScaledGamma(l,i)*rsys_in(l)
     $                 + ScaledOmega(l,i)*rsys_in(l)*rsys_in(l)
               endif
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
C Penalty term, unity by default
         A(i,i)  =  SysPriorScale(i)
      enddo

!$OMP PARALLEL DO

      do l=1,nsys
         if ( SysForm(l) .eq. isNuisance ) then
C Start with "C"

            do i1=1,n_syst_meas(l)         ! loop over all data affected by this source
               i = syst_meas_idx(i1,l)     ! i -> index of the data
c            do i=1,n0_in
               if (FitSample(i) ) then

                  d_minus_t1 = daten(i) - theo(i) + ShiftExternal(i)

                  if ( list_covar_inv(i) .eq. 0) then
C Diagonal error:
                     C(l) = C(l) +  ScaledErrors(i)
     $                    *ScaledGamma(l,i)*( d_minus_t1 )
                  else
C Covariance matrix, need more complex sum:
                     i2 = list_covar_inv(i)  ! i2 -> covar. matrix index for i.
                     do j1=1,n_syst_meas(l)
                        j = syst_meas_idx(j1,l) ! j -> index of the data
c                     do j = 1, n0_in
                        if (j.ge.i) then
                           if (FitSample(j)) then
                              d_minus_t2 = daten(j) - theo(j)
     $                             + ShiftExternal(j)
                              j2 = list_covar_inv(j)
                              if (j2 .gt. 0) then
                                 add =  ScaledTotMatrix(i2,j2)
     $                                *( ScaledGamma(l,i)*d_minus_t2
     $                                + ScaledGamma(l,j)*d_minus_t1 )
                                 if (i.ne.j) then
                                    C(l) = C(l) + add
                                 else
                                    C(l) = C(l) + 0.5*add
                                 endif
                              endif
                           endif
                        endif
                     enddo
                  endif
               endif
            enddo
C Now A:

            do i=1,n0_in
               do k=l,NSys
C
               if ( (sysform(k) .eq. isNuisance ) ! ) then
     $              .and.HaveCommonData(k,l) ) then

c                  do i1 = 1,n_syst_meas(k)
c                     i = syst_meas_idx(i1,k)
c
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
C                            do j=i,n0_in
                              if ( j.ge.i .and. FitSample(j) ) then
                                 j2 = list_covar_inv(j)
                                 if (j2 .gt. 0) then
                                    add =
     $                                   ScaledTotMatrix(i2,j2)
     $                             *( ScaledGamma(l,i)*ScaledGamma(k,j)
     $                               +ScaledGamma(l,j)*ScaledGamma(k,i))
                                    if ( i.ne.j) then
                                       A(k,l) = A(k,l) + add
                                    else
                                       A(k,l) = A(k,l) + 0.5*add
                                    endif
                                 endif
                              endif
                           enddo

                        endif
                     endif
c                  enddo

                  endif
               enddo
            enddo
         endif
      enddo

!$OMP END PARALLEL DO

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

C Also dump correlation matrix for PDF eigenvectors, if present
         if (iflag.eq.3) then

C Loop over all sources, find theory sources, count them.
            nsystheo = 0
            do l=1,nsys
               if ( ISystType(l) .eq. iTheorySyst) then
                  nsystheo = nsystheo + 1
                  itheoisys(nsystheo) = l  ! reference from "theory" index to "sys" index
               endif
            enddo


            if (nsystheo.gt.0) then

               npdf = nsystheo

               Allocate(AA(npdf,npdf))

               open (52,file=trim(OutDirName)//'/pdf_shifts.dat',
     $              status='unknown')
               write (52,'(''LHAPDF set='',A32)')
     $              trim(adjustl(LHAPDFSET))
               write (52,'(i3)') npdf

               do l=1,npdf
                  write (52,'(i3,2F8.4)') l,
     $                 rsys_in(itheoisys(l)),
     $                 ersys_in(itheoisys(l))
               enddo
               close (52)


               open (52,file=trim(OutDirName)//'/pdf_vector_cor.dat'
     $              ,status='unknown')
               write (52,'(i3)') nsys-nsysdata
               do l=1,npdf
                  write (52,'(i3,200F8.4)')  l, (
     $                 A(itheoisys(k),itheoisys(l))
     $                 /ersys_in(itheoisys(k))
     $                 /ersys_in(itheoisys(l)),
     $                 k=1,npdf)
                  do k=1,npdf
                     AA(k,l) = A(itheoisys(k),itheoisys(l))
                  enddo
               enddo
               close (52)

C Also rotation matrix:
               Call MyDSYEVD(Npdf,AA,Npdf,C,ifail)

C scale to take into account error reduction
               do i=1,Npdf
                  do j=1,Npdf
                     AA(j,i) = AA(j,i)*sqrt(C(i))
                  enddo
               enddo

C We want to preserve original directions as much as possible

               Allocate(RR(Npdf,Npdf))
               Allocate(AA2(Npdf,Npdf))

               do k=npdf,1,-1

                  do i=1,npdf
                     do j=1,npdf
                        RR(i,j) = 0.
                        if (i.eq.j) then
                           RR(i,j) = 1.
                        endif

                        RR(i,j) = RR(i,j) + AA(k,i)*AA(k,j)

                     enddo
                  enddo


                  Call MyDSYEVD(k,RR,Npdf,C,ifail)
C rotate rotation matrix:
                  do i=1,k
                     do j=1,k
                        AA2(i,j) = 0.
                        do l=1,k
                           AA2(i,j) = AA2(i,j) + AA(i,l)*RR(l,j)
                        enddo
                     enddo
                  enddo


                  do i=1,k
                     do j=1,k
                        AA(i,j) = AA2(i,j)
                     enddo
                  enddo
               enddo ! loop over k.

C Last loop to keep the direction of the original vectors
               do i=1,npdf
                  if (AA(i,i).lt.0) then
                     do j=1,npdf
                        AA(j,i) = -AA(j,i)
                     enddo
                  endif
               enddo


               open (52,file=trim(OutDirName)//'/pdf_rotate.dat'
     $              ,status='unknown')
               write (52,'(''LHAPDF set='',A32)')
     $              trim(adjustl(LHAPDFSET))
               write (52,'(i4)') Npdf
               do i=1,Npdf
C                  print *,'haha',i,C(i),ifail
                  write (52,'(i5,200F10.6)') i,
     $                 (AA(j,i),j=1,Npdf)
               enddo
               close (52)

               DeAllocate(AA,AA2,RR)

            endif
         endif
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
#include "ntot.inc"
#include "systematics.inc"
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
#include "ntot.inc"
#include "systematics.inc"
#include "theo.inc"
#include "indata.inc"
#include "steering.inc"

      double precision ScaledGamma(NSysMax,Ntot) ! Scaled Gamma matrix
      double precision ScaledErrors(Ntot)  !1/d^2, where d is scaled uncorrelated uncertainty, for each datapoint
      double precision ScaledTotMatrix(NCovarMax,NCovarMax)   ! stat+uncor+syst covar matrix
      double precision rsys_in(NSysMax)
      integer NDiag, list_covar(NTot), NCovar, list_diag(NTot)
      double precision fchi2_in, pchi2_in(nset), fcorchi2_in

      integer i,j, i1, j1, k
      double precision d,t, chi2, sum
      integer offdiag

      double precision SumCov(NCovarMax)

C---------------------------------------------------------------------------
      fchi2_in = 0.0D0
 ! Also zero fit/control sample chi2s
      chi2_fit = 0.
      chi2_cont = 0.
      offdiag = 0

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
         residuals(i)=(d-t+Sum)*sqrt(ScaledErrors(i))
         chi2=residuals(i)**2
C     Sums:
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
         Chi2 = 0d0
         do j1 = 1, NCovar
            j = list_covar(j1)
            Chi2 = Chi2 + SumCov(i1)*SumCov(j1)*ScaledTotMatrix(i1,j1)
            if ( ( JSET(i) .ne. JSET(j) )
     $           .and. (ScaledTotMatrix(i1,j1) .ne. 0d0 ) ) then
               if ( offdiag .eq. 0 ) then
                  call hf_errlog(15090916,
     $                 'I: Offdiag. elements in
     $ inverse covariance. Partial chisq values are set to zero')
               endif
               offdiag = offdiag+1
            endif
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

c reset partial chisq to 0, if there are offdiagonal elements
c partial chisq are not reasonably defined
      if ( offdiag .ne. 0 ) then
         do i1=1,NCovar
            i = list_covar(i1)
            pchi2_in(JSET(i)) = 0d0
         enddo
      endif

       ! print*,'chi2_calc2: ',fchi2_in

C Correlated chi2 part:
      fcorchi2_in = 0.d0
      do k=1,NSys
         fcorchi2_in = fcorchi2_in
     $        + rsys_in(k)**2 * SysPriorScale(k)
C Also, store as residuals:
         residuals(ndiag+k) = rsys_in(k)*sqrt(SysPriorScale(k))
      enddo
      fchi2_in = fchi2_in + fcorchi2_in

       ! print*,'chi2_calc3: ',fchi2_in

C---------------------------------------------------------------------------
      end

C--------------------------------------------------------------------------
C
C> @brief Calculate additional log correction factor
C> @param ScaledErrors uncertainties
C> @param[out] chi2_log log correction
C> @param n0_in number of points
C
C Fills vector chi2_poi -- log corrections for each dataset
C--------------------------------------------------------------------------
      subroutine chi2_calc_PoissonCorr(ScaledErrors, chi2_log, n0_in)

      implicit none
#include "ntot.inc"

      double precision ScaledErrors(Ntot)
      double precision chi2_log
      integer n0_in !FIXME: variable npoints should be used to get number of points


#include "indata.inc"
#include "systematics.inc"
#include "datasets.inc"
      integer i
      double precision dchi2
C-------------------------------------------------------------------------
      chi2_log = 0.D0
      do i=1,NDATASETS
         chi2_poi(i) = 0.D0
      enddo

      do i=1,n0_in
         if (FitSample(i)) then
            if ( alpha(i).gt.0 ) then
               dchi2=-log(alpha(i)*alpha(i)*ScaledErrors(i))
               chi2_poi_data(i) = dchi2
               chi2_log = chi2_log + dchi2
               chi2_poi(JSET(i)) = chi2_poi(JSET(i)) + dchi2
            endif
         endif
      enddo
      chi2_poi_tot = chi2_log

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
#include "ntot.inc"
#include "systematics.inc"
#include "theo.inc"
      double precision ScaledErrors(Ntot)  ! uncorrelated uncertainties, diagonal
      double precision ScaledGamma(NSysMax,Ntot) ! Scaled Gamma matrix
      double precision rsys_in(NSysMax)
      integer n0_in

      integer i,k
C------------------------------------------------------------------
      do i=1,n0_in
         if(ScaledErrors(i).ne.1.0D0) then
           ALPHA_MOD(i) =  1.D0/sqrt(ScaledErrors(i))
        else
c special case if no scaled errors given, i.e. given total cov matrix
            ALPHA_MOD(i) =  0.0D0
        endif
         THEO_MOD(i)  = THEO(i)
         do k=1,NSYS
            THEO_MOD(i) = THEO_MOD(i) - ScaledGamma(k,i)*RSys_in(k)
         enddo
      enddo
C------------------------------------------------------------------
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
#include "ntot.inc"
#include "steering.inc"
#include "systematics.inc"
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

      character*64 largershiftsyst
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
         if (abs(rsys_in(i)-rsys_save(i)).gt.shift) then
            shift = max(shift, abs(rsys_in(i)-rsys_save(i)))
            largershiftsyst = system(i)
         endif
c         print *, i, system(i), rsys_in(i),rsys_save(i),
c     .        abs(rsys_in(i)-rsys_save(i))
      enddo

      do i=1,nsys
         rsys_save(i) = rsys_in(i)
      enddo

C recalulate
      do k=1,nsys
         if (SysForm(k) .ne. isExternal) then
            do i=1,n_syst_meas(k)
               j =  syst_meas_idx(i,k)
               ScaledGamma(k,j) = ScaledGammaSav(k,j)
     $              + ScaledOmega(k,j)*rsys_in(k)
            enddo
         endif
      enddo

      if (LDebug) then
         print *,'shift',iteration, shift,nsys
      endif

      LStop = .false.
      if (iteration.ge. AsymErrorsIterations ) then

C ! Check if max. shift is small:
         if ( shift.gt. 0.05) then
            print *,'ERROR: Large nuisance parameter change in shift'
            print *,'from last iteration for ',largershiftsyst
            print *,'DeltaShift=',shift
            print *,'CONSIDER INCREASING AsymErrorsIterations to'
     $           ,AsymErrorsIterations+5
            call HF_errlog(13053001,
     $ 'W:UseOmegaScale: Max shift exeeds 5%. Consider increasing'
     $           //' AsymErrorsIterations')
         endif

         LStop = .true.
      endif

      end

C> @brief For asymmetric external systematic unertainty sources, update ScaledGamma matrix.
C> @param ScaledGamma scaled Gamma matrix (input-output)
C> @param ScaledOmega scaled Omega matrix (input)
C> @param rsys_in     shifts of systematic uncertaintis (input)
      Subroutine chi2_calc_asymerror_external(ScaledGamma, ScaledOmega
     $     , rsys_in)

      implicit none
C------------
#include "ntot.inc"

#include "systematics.inc"
      double precision ScaledGamma(NSysMax,Ntot) ! Scaled Gamma matrix
      double precision ScaledOmega(NSysMax,Ntot) !
      double precision rsys_in(NSysMax)
      integer i, l, j
C----
      do l=1,Nsys
         if ( SysForm(l) .eq. isExternal) then
            do i=1,n_syst_meas(l)
               j = syst_meas_idx(i,l)
               ScaledGamma(l,j) =  ScaledGamma(l,j)
     $              +  ScaledOmega(l,j)*rsys_in(l)
            enddo
         endif
      enddo
      end

C------------------------------------------------------------------------------
C
C> @brief Conversion from covariance matrix to nuisance parameter representation.
C
C> @param NDimCovar   -- Dimension of the covariance matrix
C> @param NDimSyst    -- Dimension of systematics matrix
C> @param NCovar -- Actual number of elements in the covariance matrix
C> @param Covar  -- Input covariance matrix. Output: nuisance parameters.
C> @param ANuisance -- Output nuisance parameter representation
C> @param Tolerance -- fractional sum of eigenvalues for the sourced treated as uncorrelated uncertainty. 0: NCorrelated = NCovar, 1: NCorrelated = 0.
C> @param Ncorrelated -- Output number of correlated nuisance parameters
C> @param Uncor      -- Output uncorrelated uncertainty
C  @param LSepDiag   -- Separate diagonal part
C
C--------------------------------------------------------------------------------
      subroutine GetNuisanceFromCovar( NDimCovar, NDimSyst, NCovar,
     $     Covar, ANuisance, Tolerance,
     $     Ncorrelated, Uncor, LSepDiag)
      implicit none
C--------------------------------------------------------------------------------
      integer NDimCovar, NDimSyst, NCovar
      double precision Covar(NDimCovar, NDimCovar)
      double precision ANuisance(NDimSyst, NDimCovar)
      double precision Tolerance
      integer Ncorrelated
      double precision Uncor(NDimCovar)
      logical LSepDiag

      double precision Eigenvalues(NDimCovar)
      integer ifail

      double precision factor, facMax, facMin

      double precision, allocatable :: testm(:,:),diag(:)
      double precision Sum,Run
      integer i,j,k

C--------------------------------------------------------------------------------

C Try to remove diagonal term first:

      if ( LSepDiag ) then

         allocate(testm(NCovar,NCovar))
         allocate(diag(NCovar))

         facMax = 1.0D0
         facMin = 0.0D0

C First check if the matrix positive definite

         do i=1,NCovar
            do j=1,NCovar
               testm(i,j) = Covar(i,j)
            enddo
         enddo
         Call MyDSYEVD(NCovar,testm,NCovar, EigenValues,IFail)

         if (EigenValues(1).lt.0) then
            print
     $   '(''Negative eigenvalue for the covariance matrix '',G12.3)',
     $           Eigenvalues(1)
            print *,'List of eigenvalues'
            do i=1,NCovar
               print *,i,Eigenvalues(i)
            enddo
            Call hf_errlog(2015050701,
     $           'S: Covariance matrix is not positive definite')

         endif


         do i=1,NCovar
            do j=1,NCovar
               testm(i,j) = Covar(i,j)
            enddo
         enddo
         do while ((facMax-facMin.gt.0.01).or.(Eigenvalues(1).lt.0))
            factor = 0.5*(facMax + facMin)
            do i=1,NCovar
               do j=1,NCovar
                  testm(i,j) = Covar(i,j)
               enddo
            enddo
            do j=1,NCovar
               diag(j) = sqrt(factor*covar(j,j))
               testm(j,j) = Covar(j,j) - diag(j)*diag(j)
            enddo

            Call MyDSYEVD(NCovar,testm,NCovar, EigenValues,IFail)


            if (EigenValues(1).lt.0) then
               facMax = factor
            else
               facMin = factor
            endif
         enddo
         ! Ok, subtract diagonal:
         do j=1,NCovar
            Covar(j,j) = Covar(j,j) - diag(j)*diag(j)
         enddo
         DeAllocate(testm)

      endif

      Call MyDSYEVD(NCovar,Covar,NDimCovar, EigenValues,IFail)

      Sum = 0
      do i=1,NCovar
c         print *,'Eig',i,EigenValues(i)
         Sum = Sum + EigenValues(i)
      enddo

      Ncorrelated = NCovar
      Run = 0.0D0
      do i=NCovar,1,-1
         Run = Run + EigenValues(i)/Sum
C         print *,i,EigenValues(i)/Sum,Run, 1.D0-Tolerance
         if (Run  .gt. 1.D0 - Tolerance) then
            Ncorrelated = NCovar - i + 1
            goto 17
         endif
      enddo
 17   continue
c      print *,NCorrelated

      do i=1,NCovar
         do j=1,NCovar
            Covar(j,i) = Covar(j,i)*sqrt(max(0.0D0,EigenValues(i)))
            ANuisance(i,j) = Covar(j,i)
C            print *,j,i, ANuisance(i,j),Covar(i,j)
         enddo
      enddo

      do j=1, NCovar
         do i=1,NCorrelated
            ANuisance(i,j) = Covar(j,NCovar-i+1)
         enddo
         do i=NCorrelated+1,NCovar
            Uncor(j) = Uncor(j) + Covar(j,NCovar-i+1)**2
         enddo
         Uncor(j) = sqrt(Uncor(j))
      enddo

      if (LSepDiag) then
       ! Add diag back to uncor:
         do j=1, NCovar
            Uncor(j) = sqrt( Uncor(j)**2 + diag(j)**2 )
         enddo
         DeAllocate(diag)
      endif

      end

      subroutine CovMatrixConverter(fileName)
      implicit none
#include "ntot.inc"
#include "systematics.inc"
#include "covar.inc"
      character*(*) fileName
      integer ncovar
      double precision tolerance
      namelist/Covar/NCovar,Tolerance
      integer i,j,k, ncorr
      double precision test,devmax,dev
C----------------------------------------

C      call hf_errlog(1,'I:Read covariance matrix from file')
      open (51,file=trim(FileName),status='old',err=91)
      read (51,nml=Covar,err=92,end=93)


      do i=1,NCovar
         read (51,*,err=99) ( Cov(i,j),j=1,NCovar)
      enddo
      do i=1,NCovar
         do j=1,NCovar
            corr_syst(i,j) = cov(i,j)
         enddo
      enddo

      close(51)

      Call GetNuisanceFromCovar(NTot,NSysMax,NCovar,Cov,Beta,
     $     Tolerance, NCorr, alpha, .false.)

      print *,'Nuisance paramters (point, Uncor, Corr1, ... CorrN):'
      print *,'NCorr=',NCorr
      do i=1,NCovar
         print '(I4,300E12.4)',i,alpha(i),(beta(j,i),j=1,NCorr)
      enddo

      dev = 0.
      print *,'  '
      print *,'Test Covariance matrix:'
      print '(2A4,3A12)','i','j',' nui ',' orig ','corr.diff%'
      do i=1,NCovar
         do j=1,NCovar
            test = 0.
            do k=1,NCorr
               test = test + beta(k,i)*beta(k,j)
            enddo
            if (i.eq.j) then
               test = test + alpha(i)*alpha(i)
            endif
            dev = (test-corr_syst(i,j))/sqrt(corr_syst(i,i)*corr_syst(j,j))*100.
            print '(2i4,3F12.2)',i,j,test,corr_syst(i,j)
     $           ,dev
            devmax = max(abs(dev),devmax)
         enddo
      enddo
      print '(''Maximum deviation from the original correlation = ''
     $     ,F10.2,''%'')',devmax

      return
 91   call hf_errlog(1,'F:Can not open covar.in file')
 92   call hf_errlog(2,'F:Can not read Covar namelist')
 93   call hf_errlog(3,'F:Can not find Covar namelist')
 99   call hf_errlog(3,'F:Error reading cov. matrix')
      end

C---------------------------------------------------------
C
C @brief redunce number of nuisance parameters by first constructing the covariance matrix and keeping only impprtant vectors
C
      subroutine reduce_nui(UncorNew,UncorConstNew
     $     , UncorPoissonNew)
      implicit none

#include "ntot.inc"
#include "systematics.inc"
#include "indata.inc"
      integer iCovarType

      double precision UncorNew(NTot),UncorConstNew(NTot),
     $     StatNew(NTot), StatConstNew(NTot), UncorPoissonNew(Ntot)

      integer isys_scaling, isys, nExternal

      double precision, allocatable :: C(:,:)    ! covariance matrix
      double precision, allocatable :: S(:,:,:)  ! nuisance param representation of it.

      double precision uncor_loc(NTOT, 0:n_sys_scaling_max-1)
      integer nui_cor(0:n_sys_scaling_max-1)
      logical l_present(0:n_sys_scaling_max-1)

      integer i,j

      character*3 name_t(0:n_sys_scaling_max-1), name_n
      data name_t/':A',':M',':P'/

      character*80 name_s

      double precision tolerance
      logical lfirst
      data lfirst/.true./
      namelist/ReduceSyst/do_reduce,tolerance, useBlas, nExternal

C------------------------------------------------
      useBlas = .false.
      if (lfirst) then
         lfirst = .false.
         Tolerance = 0.
         do_reduce = .false.
         nExternal = 0
         open (51,file='steering.txt',status='old')
         read (51,NML=ReduceSyst,end=19,err=17)
 19      continue
         close (51)
      endif

      if (.not. do_reduce) return

      if (nExternal .eq. 0) then
         LConvertCovToNui = .true.
      endif

      ! Allocate covariance matrix for the data
      Allocate(C(npoints,npoints))
      ! Allocate space to save syst. vecctors
      Allocate(S(npoints,npoints,n_sys_scaling_max))

      do isys_scaling=0,n_sys_scaling_max-1  ! loop over scaling type
         ! Clean the covariance matrix
         C = 0.0
         l_present(isys_scaling) = .false.
         do isys=1,nSys
            if ( SysScalingType(isys).eq.isys_scaling ) then
               l_present(isys_scaling) = .true.

               ! add to covariance matrix
               do i=1,NPoints
                  do j=1,NPoints
                     C(i,j) = C(i,j) + beta(isys,i)*daten(i)
     $                                *beta(isys,j)*daten(j)
                  enddo
               enddo
            endif
         enddo
         ! translate to
         if (l_present(isys_scaling)) then
            call GetNuisanceFromCovar(NPoints,NPoints,NPoints,
     $           C, S(1,1,isys_scaling+1), tolerance,
     $           Nui_cor(isys_scaling),
     $           Uncor_loc(1, isys_scaling), .false.)
         else
            Nui_cor(isys_scaling) = 0
         endif

      enddo


      ! now re-set all systematic sources and set them new.
      NSys = 0
      N_Syst_Meas = 0

      UncorNew = sqrt( UncorNew**2 + Uncor_loc(:,isLinear)**2)
      UncorPoissonNew = sqrt( UncorPoissonNew**2
     $     + Uncor_loc(:,isPoisson)**2)
      UncorConstNew = sqrt( UncorConstNew**2
     $     + Uncor_loc(:,isNoRescale)**2)


      do isys_scaling=0,n_sys_scaling_max-1
         ! use type to define new name
c         if (isys_scaling .eq. isNoRescale ) then
c            name_t = ':A'
c         endif
         do isys=1,Nui_Cor(isys_scaling)
            if (isys.lt.10) then
               write (name_n,'(''00'',i1)') isys
            elseif (isys.lt.100) then
               write (name_n,'(''0'',i2)') isys
            elseif (isys.lt.1000) then
               write (name_n,'(i3)') isys
            endif
            if (isys .le. nExternal) then
               name_s = 'reduced_'//name_n//':E'
            else
               name_s = 'reduced_'//name_n
            endif
            call AddSystematics(trim(name_s)//name_t(isys_scaling))
            do i=1,NPoints
               n_syst_meas(NSYS) = n_syst_meas(NSYS) + 1
               syst_meas_idx(n_syst_meas(NSYS),NSYS) = i
               beta(NSYS,i) =  S(isys,i,isys_scaling+1) /daten(i)
            enddo
         enddo
      enddo

      Deallocate(C)
      Deallocate(S)

      goto 18
 17   continue
      call hf_errlog(1,
     $     'F:Error reading ReduceSyst Namelist ! Stop')
 18   continue

      end

C---------------------------------------------------------------
C> @brief Convert covariance matricies to nuisance param.
C>
C---------------------------------------------------------------

      subroutine covar_to_nui(UncorNew,UncorConstNew,StatNew
     $     ,StatConstNew, UncorPoissonNew)
      implicit none

#include "ntot.inc"
#include "covar.inc"
#include "indata.inc"
#include "theo.inc"
#include "systematics.inc"
#include "datasets.inc"

      integer NCovar
      integer List_Covar(NCovarMax)

      double precision UncorNew(NTot),UncorConstNew(NTot),
     $     StatNew(NTot), StatConstNew(NTot), UncorPoissonNew(Ntot)

      double precision cov_loc(NCovarMax,NCovarMax)
      double precision anui_loc(NCovarMax,NCovarMax)
      double precision uncor(NCovarMax), Diag(NCovarMax)
      double precision uncorSt(NCovarMax)  ! from stat. subtraction

      double precision unc(NTot) ! input uncor. errors
      double precision sta(NTot) ! input stat. errors

      integer nui_cor
      integer i,j,i1,j1
      integer iCovType, iCovBit

      double precision cor
      double precision Stat,StatConst,tot1,tot2

      logical LFirst
      data LFirst /.true./

      character *80 DataName(NSET)
      character *3  DataSystType(NSET)

      double precision Tolerance
      logical LSubtractStat
      namelist/CovarToNuisance/Tolerance, LConvertCovToNui,
     $     DataName,DataSystType, LSubtractStat

      character*8 sys_prefix(NCovTypeMax)
      data sys_prefix/'CTot_','CTot_','CSyst_','CStat_','CSyst_'/

      character*80 name_s
      character*3  name_n, name_t
      integer imaskSta, imaskUnc
      character*80 message

C--------------------------------------------------------
      if (LFirst) then
         LFirst = .false.
         Tolerance   = 0.
         LConvertCovToNui = .false.
         LSubtractStat = .false.
         do i=1,NSET
            DataName(i) = ''
            DataSystType(i) = ':M'
         enddo

         open (51,file='steering.txt',status='old')
         read (51,NML=CovarToNuisance,end=19,err=17)
 19      continue
         close (51)


      endif


      do i=1,Npoints
         theo(i) = daten(i) ! Set theory = data for error scaling
         call GetPointErrors(i, Stat, StatConst, Unc(i))
         Sta(i) = sqrt(Stat**2+StatConst**2)

         theo(i) = 0.0D0

C Set "new" uncorrelated uncertainties to the "old" values:
         UncorNew(i) =      e_uncor_mult(i)
         StatNew(i)  =      e_stat_poisson(i)
         UncorConstNew(i) = e_uncor_const(i)
         StatConstNew(i)  = e_stat_const(i)
         UncorPoissonNew(i) = e_uncor_poisson(i)

      enddo

      do iCovType = 1, NCovTypeMax
         iCovBit = IShft(1,ICovType-1)

C Create list first:
         NCovar = 0
         do i=1,Npoints
            if ( IAND(icov_type(i),iCovBit).eq.iCovBit ) then
               NCovar = NCovar + 1
               List_Covar(NCovar) = i
            endif
         enddo

c         print *,iCovType,Icovbit,NCovar

         if (NCovar.eq.0) then
            cycle   ! nothing to be done
         endif

C Generate compact matrix:
         do i1=1,NCovar
            i = List_Covar(i1)
            do j1=1,NCovar
               j = List_Covar(j1)
C Use proper source:
               if (iCovBit .eq. iCovSyst) then
                  cov_loc(i1,j1) = cov(i,j)
               elseif (iCovBit .eq. iCovSystCorr) then
                  cov_loc(i1,j1) = corr_syst(i,j)*Unc(i)*Unc(j)
               elseif (iCovBit .eq. iCovStatCorr) then
                  cov_loc(i1,j1) = corr_stat(i,j)*Sta(i)*Sta(j)

               elseif (iCovBit .eq. iCovTotal) then
                  cov_loc(i1,j1) = cov(i,j)
               elseif (iCovBit .eq. iCovTotalCorr) then
                  tot1 = sqrt(sta(i)**2+unc(i)**2)
                  tot2 = sqrt(sta(j)**2+unc(j)**2)
                  cov_loc(i1,j1) = corr(i,j)*tot1*tot2

               endif
            enddo
         enddo
c      endif

         do j=1,NCovar
            Uncor(j) = 0.D0
            UncorSt(j) = 0.D0
         enddo

         if ( (iCovBit.eq.iCovSyst) .or. (iCovBit.eq.iCovSystCorr) )
     $        then
C Direct diagonalisation:
            Call GetNuisanceFromCovar(NCovarMax, NCovarMax,NCovar,
     $           cov_loc, anui_loc, Tolerance,Nui_cor,Uncor,.false.)
         elseif ( (iCovBit.eq.iCovStatCorr)
     $           .or.(iCovBit.eq.iCovTotal)
     $           .or.(iCovBit.eq.iCovTotalCorr) ) then
C Subtract diagonal as much as possible:

            if (LSubtractStat) then
               Call SubtractStat(cov_loc,sta,NCovarMax,NCovar, UncorSt)
            endif

            Call GetNuisanceFromCovar(NCovarMax, NCovarMax,NCovar,
     $           cov_loc, anui_loc, Tolerance,Nui_cor,Uncor,.true.)

            if (LSubtractStat) then
               do j=1,NCovar
                  Uncor(j) = sqrt(Uncor(j)**2 + UncorSt(j)**2)
               enddo
            endif

         else
            Call hf_errlog(142107,
     $   'W:Nuisance rep. code not ready for this  error type yet')
         endif

C         print *,'hihi',Nui_cor,tolerance, ncovar


C Define the scaling property based on the first point:
         imaskSta = iStatTypesBitMask(JSet(List_Covar(1)))
         imaskUnc = iUncorTypesBitMask(JSet(List_Covar(1)))


         if ((iCovBit.eq.iCovSyst).or.(iCovBit.eq.iCovSystCorr))   then
                  ! Multiplicative is default for syst.
            name_t = ':M'

            ! Check bits
            if (iAnd(imaskUnc,ibLinear).eq.imaskUnc) then
               name_t = ':M'
            elseif (iAnd(imaskUnc,ibConst).eq.imaskUnc) then
               name_t = ':A'
            else
               Call hf_errlog(14090401,
     $ 'W: inconsistent use of uncor and '//
     $ 'uncor const for dataset: "'//
     $              Trim(DataSetLabel(JSet(List_Covar(1))))
     $              //'".  Use multiplicaiive errors')
            endif

         elseif ( (iCovBit.eq.iCovStatCorr) ) then
                  ! Poisson is default for stat.
            name_t = ':P'

            ! Check bits
            if (iAnd(imaskSta,ibPoisson).eq.imaskSta) then
               name_t = ':P'
            elseif (iAnd(imaskSta,ibConst).eq.imaskSta) then
               name_t = ':A'
            else
               Call hf_errlog(14090402,
     $ 'W: inconsistent use of stat and '//
     $ 'stat const for dataset "'//
     $              Trim(DataSetLabel(JSet(List_Covar(1))))
     $              //'".  Use Poisson errors')
            endif

         else                   ! Additive is defualt for full
            name_t = ':A'
         endif

         do i=1,NSET
            if ( DataName(i).eq.DataSetLabel(JSet(List_Covar(1))))
     $           then
               name_t = DataSystType(i)
            endif
         enddo

         do j1=1,Nui_cor
            if (j1.lt.10) then
               write (name_n,'(''00'',i1)') j1
            elseif (j1.lt.100) then
               write (name_n,'(''0'',i2)') j1
            else
               write (name_n,'(i3)') j1
            endif
            name_s = trim(Sys_prefix(icovtype))
     $           // name_n // '_'
     $           // DataSetLabel( JSet(List_Covar(1)) )


            Call AddSystematics(Trim(name_s)//name_t)
            do i1=1,NCovar
               i = List_Covar(i1)  ! point to the data
               n_syst_meas(NSYS) = n_syst_meas(NSYS) + 1
               syst_meas_idx(n_syst_meas(NSYS),NSYS) = i
               beta(NSYS,i) =  anui_loc(j1,i1)/daten(i)
            enddo
         enddo

         do i1=1,NCovar
            i = List_Covar(i1)
            if ( (iCovBit.eq.iCovSyst) .or. (iCovBit.eq.iCovSystCorr) )
     $           then
C Re-set uncorrelated systematics:
               if (name_t .eq. ':M') then
                  UncorNew(i) = Uncor(i1)/daten(i)
                  UncorPoissonNew(i) = 0.0D0
                  UncorConstNew(i) = 0.0D0
               elseif (name_t .eq. ':A') then
                  UncorNew(i) = 0.0D0
                  UncorPoissonNew(i) = 0.0D0
                  UncorConstNew(i) = Uncor(i1)/daten(i)
               elseif (name_t .eq. ':P') then
                  UncorNew(i) = 0.0D0
                  UncorPoissonNew(i) = Uncor(i1)/daten(i)
                  UncorConstNew(i) = 0.0D0
               else
                  call hf_errlog(22071401,
     $                 'S: Unknowns scale '//name_t//' STOP')
               endif
            elseif (iCovBit.eq.iCovStatCorr) then
               if (name_t .eq. ':A') then
                  StatNew(i) = 0.0D0
                  StatConstNew(i) = Uncor(i1)/daten(i)
               elseif (name_t .eq. ':P') then
                  StatNew(i) = Uncor(i1)/daten(i)
                  StatConstNew(i) = 0.0D0
               else
                  call hf_errlog(22071401,
     $                 'S: Unknowns scale '//name_t//' STOP')
               endif
            elseif ( (iCovBit.eq. iCovTotalCorr)
     $              .or.(iCovBit.eq. iCovTotal) ) then
               if (name_t .eq. ':A') then
                  StatNew(i) = 0.0D0
                  StatConstNew(i) = Uncor(i1)/daten(i)
                  UncorNew(i) = 0
                  UncorPoissonNew(i) = 0
                  UncorConstNew(i) = 0
               elseif (name_t .eq. ':P') then
                  StatNew(i) = Uncor(i1)/daten(i)
                  StatConstNew(i) = 0
                  UncorNew(i) = 0
                  UncorPoissonNew(i) = 0
                  UncorConstNew(i) = 0
               else
                  call hf_errlog(22071401,
     $                 'S: Unknowns scale '//name_t//' STOP')
               endif
            endif
         enddo

      enddo
      goto 18
 17   continue
      call hf_errlog(1,
     $     'F:Error reading CovarToNuisance Namelist ! Stop')
 18   continue
      end


      subroutine SubtractStat(cov,sta,ncovMax,ncov,unc)
      implicit none
      integer ncovMax, ncov
      double precision cov(ncovMax,ncovMax),sta(*),unc(*)
      integer i,j, IFail
      double precision, allocatable :: c_loc(:,:), eig(:)
C----------------------------------------------------------

C First check if this is possible:
      allocate( c_loc(ncov,ncov) )
      allocate( eig(ncov) )
      do i=1,ncov
         do j=1,ncov
            c_loc(i,j) = cov(i,j)
         enddo
         c_loc(i,i) = c_loc(i,i) - sta(i)*sta(i)
      enddo

      call MyDSYEVD(ncov, c_loc, ncov, eig, IFail)

      if (eig(1).gt.0) then
         do i=1,ncov
            unc(i) = sta(i)
            cov(i,i) = cov(i,i) - sta(i)*sta(i)
         enddo
      else
         do i=1,ncov
            unc(i) = 0.0
         enddo
         call hf_errlog(18081601,
     $'W:Can not subtract stat. component from covariance matrix:'//
     $' perhaps it is not diagonal?')
      endif

      deallocate (c_loc)
      deallocate ( eig )
      end


C--------------------------------------------------------
C> @Brief Interface to lapack, to dynamically allocate work arrays
      subroutine MyDSYEVD(NCovar,Covar,NDimCovar, EigenValues,ifail)
      implicit none
      integer NCovar, NDimCovar
      double precision Covar(NDimCovar,NDimCovar), EigenValues(NCovar)
      integer IFail
      double precision Work
      integer IWork
      Character*80 msg
C---------------------------------------------------------------
C Determine optimal size of the work array:
      Call DSYEVD('V','U',NCovar,Covar,NDimCovar, EigenValues, Work,
     $     -1, IWork, -1, IFail)


      write (msg,
     $ '(''I:MyDSYEVD: optimal dimensions for work arrays:'',2i6)')
     $     int(work)+1, iwork
      call HF_ERRLOG(14121701,msg)
      call MyDSYEVD2(NCovar,Covar,NDimCovar, EigenValues,
     $     int(work)+1,iwork,ifail)

      end

      subroutine MyDSYEVD2(NCovar,Covar,NDimCovar, EigenValues, nwork,
     $     nlwork,ifail)
      implicit none
      integer NCovar, NDimCovar
      double precision Covar(NDimCovar,NDimCovar), EigenValues(NCovar)
      integer nwork, nlwork
      double precision Work(nwork)  ! Dynamic array
      integer IWork(nlwork)         ! Dynamic array
      integer IFail
C---------------------------------------------------------------------
      Call DSYEVD('V','U',NCovar,Covar,NDimCovar, EigenValues, Work,
     $     nwork, IWork, nlwork, IFail)


      end
