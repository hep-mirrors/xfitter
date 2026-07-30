C--------------------------------------------------------
C> @file GetChisquare_new.f
C> @Brief Calculate chisquare using BLAS-3 calls wherever possible.
C
C> This is a copy of GetChisquare.f that is intended to be edited.
C> Implemented with _new subroutines so the old method
C> can still be used whenever needed. 
C> In the future we could possibly just use this mehtod,
C> or at least set it as the default. But for now, it will 
C> be the experimental path where optimisations are tested
C> and added.
C
C
C---------------------------------------------------------

C---------------------------------------------------------------------------
C> @brief Runtime switch: true if OpenMP.useNewChi2 is set (cached).
C---------------------------------------------------------------------------
      logical function UseNewChi2Flag()
      implicit none
      integer GetParamI
      external GetParamI
      logical lfirst, flag
      data lfirst /.true./, flag /.false./
      save lfirst, flag
      if (lfirst) then
         lfirst = .false.
         flag = (GetParamI('UseNewChi2') .ne. 0)
         if (flag) then
            print '(A)',
     $       ' Chi2: NEW calculation path (GetChisquare_new.f)'
         endif
      endif
      UseNewChi2Flag = flag
      end

   

      subroutine GetNewChisquare_new(flag_in,n0_in,fchi2_in,rsys_in,ersys_in,pchi2_in,
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

C Allocated on first call rather than declared as fixed-size locals.
C NSysMax x Ntot arrays push the libxfitter.dylib image past the 4 GiB that
C dyld can load on arm64 (throws a "Symbol not found").
      double precision, allocatable :: ScaledGamma(:,:)    ! Scaled Gamma matrix
      double precision, allocatable :: ScaledGammaSav(:,:) ! Scaled Gamma matrix, saved

      double precision, allocatable :: ScaledOmega(:,:)    ! Scaled Omega matrix
      save ScaledGamma, ScaledGammaSav, ScaledOmega

      double precision ScaledErrors(Ntot)  ! uncorrelated uncertainties, diagonal

C Covariance matrices, allocated on first call with the actual number of
C covariance points (NCovarDim), so large stat. correlation matrices
C (e.g. 2000x2000) fit without a compile-time NCovarMax limit.
      double precision, allocatable :: ScaledErrorMatrix(:,:) ! stat+uncor error matrix
      double precision, allocatable :: ScaledSystMatrix(:,:)  ! syst. covar matrix
      double precision, allocatable :: ScaledTotMatrix(:,:)   ! stat+uncor+syst covar matrix
      integer NCovarDim       ! Leading dimension of the covariance matrices

      integer NDiag, NCovar   ! Number of diagonal and full covariance input data points
      integer List_Diag(Ntot), List_Covar(Ntot), List_Covar_inv(Ntot)
      save ScaledErrorMatrix, ScaledSystMatrix, ScaledTotMatrix
      save NCovarDim, NDiag, NCovar, List_Diag, List_Covar,
     $     List_Covar_inv

      integer Iterate
      logical LFirst
      data LFirst /.true./

      integer omegaIteration
      Logical doMatrix, doNuisance, doExternal, LStop
C These are set on the first call only, so they must survive between calls
C (not relying on -fno-automatic for this):
      save doMatrix, doNuisance, doExternal

C Timing variables for profiling (xf_wtime = wall-clock if OpenMP, cpu_time otherwise)
      double precision t_start, t_gamma, t_covar, t_stat, t_sumcovar
      double precision t_syst_shifts, t_chi2, t_total
      double precision xf_wtime
      external xf_wtime
      logical lPrintTiming
      data lPrintTiming/.true./  ! Set to .false. to disable timing output

c Global initialisation

      if (lPrintTiming) t_start = xf_wtime()

      if (LFirst) then
         LFirst = .false.

C    !> Initialize LAPACK solver settings
         call init_lapack_solver()

C    !> Determine which mechanisms for syst. errors should be used:
         Call Init_Chi2_calc(doMatrix, doNuisance, doExternal)

C    !> Determine which errors are diagonal and which are using covariance matrix
         Call init_chi2_stat(NDiag, NCovar, List_Diag, List_Covar,
     $        List_Covar_inv,n0_in)

         print*,'NDiag=',NDiag,'NCovar=',NCovar

C    !> Allocate the scaled syst. matrices (see declarations above)
         allocate(ScaledGamma(NSysMax,Ntot))
         allocate(ScaledGammaSav(NSysMax,Ntot))
         allocate(ScaledOmega(NSysMax,Ntot))

C    !> Allocate covariance matrices with the actual size
         NCovarDim = max(NCovar,1)
         allocate(ScaledErrorMatrix(NCovarDim,NCovarDim))
         allocate(ScaledSystMatrix(NCovarDim,NCovarDim))
         allocate(ScaledTotMatrix(NCovarDim,NCovarDim))

         do i=1,NCovarDim
            do j=1,NCovarDim
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
         if (lPrintTiming) then
            t_gamma = xf_wtime()
            print '(A,F8.3,A)', ' TIMING: GetGamma      = ',
     $           t_gamma-t_start,' s'
         endif

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
     $           ,ScaledSystMatrix, NCovarDim
     $           ,List_Covar_Inv,n0_in)
            if (lPrintTiming) then
               t_covar = xf_wtime()
               print '(A,F8.3,A)', ' TIMING: Chi2_covar    = ',
     $              t_covar-t_gamma,' s'
            endif
         endif
      else
C !> Restore saved gamma:
         do k=1,nsys
            do i=1,n_syst_meas(k)
               j =  syst_meas_idx(i,k)
               ScaledGamma(k,j) = ScaledGammaSav(k,j)
            enddo
         enddo
         if (lPrintTiming) then
            t_covar = xf_wtime()
         endif
      endif


C !> Get uncor errors/nuisance parameters.
      do while ( Iterate.ge.0 )

c !> First recalc. stat. and bin-to-bin uncorrelated uncertainties:
         if (.not. Chi2FirstIterationRescale .or. flag_in.eq.1 .or.
     $  Chi2OffsRecalc) then
            if (lPrintTiming) t_covar = xf_wtime()  ! Reset timing reference
            Call Chi2_calc_stat_uncor(ScaledErrors
     $           ,ScaledErrorMatrix, NCovarDim
     $           ,rsys_in,n0_in, NCovar, List_Covar, Iterate)
            if (lPrintTiming) then
               t_stat = xf_wtime()
               print '(A,F8.3,A)', ' TIMING: stat_uncor    = ',
     $              t_stat-t_covar,' s'
            endif

C  !> Sum covariance matricies and invert the total:

            if ( doMatrix .or. NCovar .gt. 0 ) then
               Call Chi2_calc_SumCovar_new(ScaledErrorMatrix,
     $              ScaledSystMatrix,
     $              ScaledTotMatrix, NCovarDim, NCovar)
               if (lPrintTiming) then
                  t_sumcovar = xf_wtime()
                  print '(A,F8.3,A)', ' TIMING: SumCovar+Inv  = ',
     $                 t_sumcovar-t_stat,' s'
               endif
            else
               if (lPrintTiming) then
                  t_sumcovar = xf_wtime()
               endif
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
     $              ,ScaledTotMatrix, NCovarDim
     $              ,ScaledGamma
     $              ,rsys_in,ersys_in,list_covar_inv, flag_in, n0_in
     $              ,scaledOmega)
            endif
            if (lPrintTiming) then
               t_syst_shifts = xf_wtime()
               print '(A,F8.3,A)', ' TIMING: syst_shifts   = ',
     $              t_syst_shifts-t_sumcovar,' s'
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
     $     ScaledTotMatrix, NCovarDim,
     $     rsys_in,
     $     ndiag, list_diag, ncovar, list_covar,
     $     fchi2_in, pchi2_in, fcorchi2_in)
      if (lPrintTiming) then
         t_chi2 = xf_wtime()
         print '(A,F8.3,A)', ' TIMING: chi2_calc     = ',
     $        t_chi2-t_syst_shifts,' s'
         t_total = xf_wtime()
         print '(A,F8.3,A)', ' TIMING: TOTAL chi2    = ',
     $        t_total-t_start,' s'
         print *, '----------------------------------------'
      endif


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


C------------------------------------------------------------------------------------
C> Now add new functions in one by one, if they need to be changed.
C> Unchanged functions can remain in GetChisquare.f 
C> The leftover functions are:
C>   - subroutine Init_Chi2_calc(doMatrix, doNuisance, doExternal)
C>   - subroutine Init_chi2_stat(NDiag, NCovar, List_Diag, List_Covar,
C>     $     List_Covar_inv, n0_in)
C>   - subroutine Chi2_calc_readExternal( rsys_in, ersys_in, IFlag )
C>   - subroutine chi2_calc_GetGamma(ScaledGamma, ScaledOmega)
C>   - subroutine chi2_calc_covar(ScaledGamma,ScaledSystMatrix,
C>     $     NCovarDim, List_Covar_Inv, n0_in)
C>   - Subroutine GetPointErrors(Idx, Stat, StatConst, Uncor)
C>   - Subroutine Chi2_calc_SumCovar(ScaledErrorMatrix, ScaledSystMatrix,
C>     $              ScaledTotMatrix, NCovarDim, NCovar)
C>   - subroutine chi2_calc_stat_uncor(ScaledErrors, ScaledErrorMatrix,
C>     $     NCovarDim, rsys_in,n0_in, NCovar, List_Covar, Iterate)
C>   - subroutine expand_syst_lists(tot_matrix,NCovarDim,
C>     $     list_covar_inv,n0_in)
C>   - subroutine chi2_calc_syst_shifts_simple(
C>     $     ScaledErrors
C>     $     ,ScaledGamma
C>     $     ,rsys_in,   n0_in)
C>   - subroutine chi2_calc_syst_shifts(
C>     $     ScaledErrors
C>     $     ,ScaledTotMatrix, NCovarDim
C>     $     ,ScaledGamma
C>     $     ,rsys_in,ersys_in,list_covar_inv,  iflag, n0_in, ScaledOmega)
C>   - subroutine Sys_Data_list12(isys1,isys2,n_list,i_list)
C>   - subroutine chi2_calc_chi2(ScaledErrors,ScaledGamma,
C>     $     ScaledTotMatrix,NCovarDim,rsys_in
C>     $     ,NDiag, List_Diag, NCovar, List_Covar
C>     $     ,fchi2_in, pchi2_in, fcorchi2_in)
C>   - subroutine chi2_calc_PoissonCorr(ScaledErrors, chi2_log, n0_in)
C>   - Subroutine Chi2_calc_FCN3(ScaledErrors,ScaledGamma,RSys_in, n0_in)
C>   - subroutine UseOmegaScale(ScaledGamma,ScaledGammaSav,ScaledOmega,
C>     $     rsys_in,Iteration,LStop)
C>   - Subroutine chi2_calc_asymerror_external(ScaledGamma, ScaledOmega
C>     $     , rsys_in)
C>   - subroutine GetNuisanceFromCovar( NDimCovar, NDimSyst, NCovar,
C>     $     Covar, ANuisance, Tolerance,
C>     $     Ncorrelated, Uncor, LSepDiag)
C>   - subroutine reduce_nui(UncorNew,UncorConstNew
C>     $     , UncorPoissonNew)
C>   - subroutine covar_to_nui(UncorNew,UncorConstNew,StatNew
C>     $     ,StatConstNew, UncorPoissonNew)
C>   - subroutine MyDSYEVD(NCovar,Covar,NDimCovar, EigenValues,ifail)
C>   I think that's it!
C------------------------------------------------------------------------------------

C Start with faster inversion in Chi2_calc_SumCovar
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
      Subroutine Chi2_calc_SumCovar_new(ScaledErrorMatrix, ScaledSystMatrix,
     $              ScaledTotMatrix, NCovarDim, NCovar)

      implicit none
#include "ntot.inc"

      integer NCovarDim
      double precision ScaledErrorMatrix(NCovarDim,NCovarDim) ! stat+uncor error matrix
      double precision ScaledSystMatrix(NCovarDim,NCovarDim)  ! syst. covar matrix
      double precision ScaledTotMatrix(NCovarDim,NCovarDim)   ! stat+uncor+syst covar matrix
      integer NCovar
      double precision Array(NCovarDim*2)
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
C This could be a cholesky factorisation here instead. 
C And just pass the triangles to the solvers downstream. 
      Call DINV_AUTO(NCovar,ScaledTotMatrix,NCovarDim,Array,IFail)
C      print *,IFail,NCovar

      end
