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
               Call Chi2_calc_syst_shifts_new(
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
C      Call DINV_AUTO(NCovar,ScaledTotMatrix,NCovarDim,Array,IFail)
C Do the cholesky factorisation, where L can be the lower triangle. 
C We will pass that to everything later, but the upper triangle will still 
C keep Vcov as is. Which I think we need for later anyway. 
      Call DPOTRF('L', NCovar, ScaledTotMatrix, NCovarDim, IFail)
      if (IFail .ne. 0) then
         Call HF_ERRLOG(26073101,
     $ 'S: Chi2_calc_SumCovar_new: covariance matrix not pos-def.'//
     $ 'If persistent, use the old chi2 method, useNewChi2=false')
      endif
      
C      print *,IFail,NCovar

      end
C-----------------------------------------------------------------


C Now we want to modify the chi2 shifts calculation to deal with the
C both covariance and nuisance paths better. Initially, the intention 
C was to just deal with covariance - but nuisance can also be done similarly. 
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
      subroutine chi2_calc_syst_shifts_new(
     $     ScaledErrors
     $     ,ScaledTotMatrix, NCovarDim
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
      integer NCovarDim
      double precision ScaledTotMatrix(NCovarDim,NCovarDim)   !> stat+uncor+syst covar matrix
      double precision ScaledGamma(NSysMax,Ntot) !> Scaled Gamma matrix
      double precision ScaledOmega(NSysMax,Ntot) ! Scaled Omega matrix

      double precision rsys_in(NSYSMax), ERSYS_in(NSYSMax)
      integer list_covar_inv(NTOT),  iflag, n0_in
      logical doExternal

      integer k,l, i1,j1,i,j, j2, i2
C Make A allocatable.
      double precision, allocatable :: A(:,:)
      double precision C(NSysMax)
      save A

      double precision, allocatable :: AA(:,:)
      double precision, allocatable :: AA2(:,:)
      double precision, allocatable :: RR(:,:)

C Some quantities for the new system
C diag - for the diagonal parts
C cov - for the covariance parts
C For the diagonal parts 
C Zdiag(i,l) = sqrt(w_i) * Gamma(l,i) where w_i = ScaledErrors(i)
C which is just the inverse variance
C So then A is just Zdiag^T Zdiag. 
C And C = Zdiag^T (sqrt(w) rdiag)
C From Cdiag = Gammadiag sqrt(w) sqrt(w) rdiag
C In some of my working out, I have referred to sqrt(w)*rdiag as zdiag
C I have done the same in the code now.
C-------------------------------
C for the covariance path, the idea is much the same
C instead of sqrt(w) sqrt(w), we can have Lcov Lcov^T
C in place of inverse covariance Vcov^-1 where Lcov is cholesky factor of Vcov
C so Zcov = L^-1 Gammacov^T (I have mixed Zcov and Zcov^T in some workings - but equivalent)
C if x=(rcov + Gammacov^T * b)
C Using the chi2 form x^T Vcov^-1 x, and the definition above,
C C = Zcov^T Lcov^-1 rcov
C A = Zcov^T Zcov
C these are I think equivalent to what came before, but now
C we can use BLAS-3 calls to speed the whole thing up.
C In this Z formalism , rows are data points, columns are systematic sources. 

      double precision, allocatable :: Zcov(:,:)
      double precision, allocatable :: Zdiag(:,:)
      double precision, allocatable :: rcov(:,:)
      double precision, allocatable :: zdiag(:,:)

      integer nd, ndd
      double precision wgt

      double precision d_minus_t1, d_minus_t2,add
      double precision ShiftExternal(NTOT)

      integer com_list(NTot),n_com_list  !> List of affected data, common for two sources.
      integer IR(2*NSysMax), Ifail,  Npdf

      integer nsystheo, itheoisys(NSysMax)
      integer nsys_sav, n0_in_sav

C Timing variables (xf_wtime = wall-clock if OpenMP, cpu_time otherwise)
      double precision t1, t2, t3, t4
      double precision xf_wtime
      external xf_wtime

      logical lfirst
      data lfirst /.true./
      data nsys_sav,n0_in_sav/0,0/
      save lfirst,nsys_sav,n0_in_sav
C-
      logical HaveCommonData(NsysMax, NsysMax)
C--------------------------------------------------------
      t1 = xf_wtime()
      if (.not. allocated(A)) then
         allocate(A(NSysMax,NSysMax))
      endif

C Check if number of sources/data points change:
      ResetCommonSyst = (nsys.ne.nsys_sav) .or. (n0_in.ne.n0_in_sav)
      nsys_sav = nsys
      n0_in_sav = n0_in 
C Determine pairs of syst. uncertainties which share  data


      if (LFirst .or. ResetCommonSyst) then
         LFirst = .false.
         ResetCommonSyst = .false.


         call expand_syst_lists(scaledtotmatrix,NCovarDim,
     $        list_covar_inv,n0_in)

C Parallelize HaveCommonData computation (O(nsys^2) calls)
!$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(k,n_com_list,com_list)
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
!$OMP END PARALLEL DO
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


C Reset the matrices:
      do i=1,nsys
         C(i) = 0.0D0
         do j=1, nsys
            A(i,j) = 0.0D0
         enddo
C Penalty term, unity by default
         A(i,i)  =  SysPriorScale(i)
      enddo

      t2 = xf_wtime()
      print '(A,F8.3,A,I6,A,I6)', '   syst_shifts init:   ',
     $     t2-t1,' s  (nsys=',nsys,' n0_in=',n0_in,')'

C Diagonal contributions first
C as stated above
C for Zdiag, we use sqrt of ScaledErrors, and we use Gamma(l,i)
C Once Zdiag is constructed, construct A with rank-k update
C Construct C with GEMV
      ndd = max(n0_in - NCovar, 1)
      nd = 0
      if ( n0_in - NCovar .gt. 0 ) then
         allocate(Zdiag(ndd,nsys))
         allocate(zdiag(ndd))

         do i=1,n0_in
            if ( list_covar_inv(i) .eq. 0 ) then
               nd = nd +1
               if ( FitSample(i) ) then
                  wgt = sqrt(ScaledErrors(i))
C zdiag = sqrt(w_i) * rdiag
                  zdiag(nd) = wgt * ( daten(i) - theo(i) + ShiftExternal(i) )
                  do l=1,nsys
                     if (SysForm(l) .eq. isNuisance ) then
                        Zdiag(nd,l) = wgt * ScaledGamma(l,i)
                     else
                        Zdiag(nd,l) = 0.0D0
                     endif
                  enddo
               else
                  zdiag(nd) = 0.0D0
                  do l=1,nsys
                     Zdiag(nd,l) = 0.0D0
                  enddo
               endif
            endif
         enddo

C Now construct A via symmetric rank-k update
C DSYRK performs alpha*A*A^T + beta*C
C or the 'T' version which is A^T*A ordered.
C Lower or Upper
C 'T'ranspose or 'N'ot? idk
C Order of C, row of A (if 'T'), alpha, A, LDA as declared in calling program
C beta, C, LDC
         call DSYRK('L', 'T', nsys, nd, 1.0D0, Zdiag, ndd,
     $    1.0D0, A, NSysMax)
C Construct C = Zdiag^T * zdiag (+C) i.e. gamma * sqrt(w) * sqrt(w) * r
C which just gives us the old gamma * w * r
C DGEMV does alpha*A*x + beta*y or alpha*A**T*x + beta*y
C 'T'ranspose, rows in A, columns in A, alpha, A, LDA
C x, increment for elements of x, beta, y, increment for elements of y.
         call DGEMV('T', nd, nsys, 1.0D0, Zdiag, ndd,
     $    zdiag, 1, 1.0D0, C, 1)

         deallocate(Zdiag)
         deallocate(zdiag)
      endif

C Covariance bits next

      t3 = xf_wtime()
      print '(A,F8.3,A)', '   syst_shifts OMP loop:',t3-t2,' s'

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
            Call DEQINV_AUTO(Nsys,A,NsysMax,IR, IFail, 1, C)
         else
            Call SPDSOLVE_AUTO(Nsys,A,NsysMax,IR,IFail,1,C)
         endif

         t4 = xf_wtime()
         print '(A,F8.3,A)', '   syst_shifts matrix solve:',t4-t3,' s'

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
