      subroutine Error_Bands_Pumplin

      implicit none

#include "steering.inc"
#include "endmini.inc"
#include "alphas.inc"
#include "thresholds.inc"
#include "ntot.inc"
#include "systematics.inc"
#include "g_offset.inc"
#include "fcn.inc"
#include "theo.inc"

      integer shift_dir
      double precision a
      dimension a(MNE)

      integer i,j,npar,idx,idx2,kflag,ii
      character*48 name,name2
      character*300 base,base2
      character tag(40)*3
      data (tag(i),i=1,40) /'s01','s02','s03','s04','s05',
     +     's06','s07','s08','s09','s10',
     +     's11','s12','s13','s14','s15',
     +     's16','s17','s18','s19','s20',
     +     's21','s22','s23','s24','s25',
     +     's26','s27','s28','s29','s30',
     +     's31','s32','s33','s34','s35',
     +     's36','s37','s38','s39','s40'/


#include "couplings.inc"

C SG: x-dependent fs:
      double precision fs0
      double precision fshermes


      character*20 parname
      integer  ind,ind2,jext,iint

      integer nparmax
      parameter (nparmax=1000)

      integer  iunint(nparmax)  ! internal param. number
      integer  iexint(nparmax)  ! external param. number
      double precision parval,parerr,parlolim,parhilim

      integer mpar
      double precision chichi
      double precision chi2data_theory !function

      double precision shift
      double precision GetUmat !function
      ! double precision DecorVarShift

C for uncertainty bands with increased tolerance
      double precision fmin, fedm, errdef
      integer npari, nparx, istat

C for theory errors:
      double precision, allocatable :: TheoVars(:,:,:)

C---------------------------------------------------------------

C
C  Fix relation between internal and external params.
C
      mpar = 0

      do ind=1,MNE
         call mnpout(ind,parname,parval,parerr,parlolim,
     $        parhilim,iunint(ind))
         if (iunint(ind).gt.0) then
            write (6,*) 'Parameter',ind,' name=',parname
            write (6,*) 'Internal index=',iunint(ind)
            write (6,*) ' '
            mpar = mpar + 1
         endif
      enddo

      print *,mpar,' variable parameters'


      do ind=1,mpar
         do ind2=1,MNE
            if (iunint(ind2).eq.ind) then
               iexint(ind) = ind2
            endif
         enddo
      enddo

      print *,'Relation internal->external index:'

      do ind=1,mpar
         write (6,*) 'internal=',ind,' external=',iexint(ind)
      enddo


      npar = MNE !> npar runs over external parameters.

      allocate(TheoVars(NTOT,2,mpar))
C
C Loop over de-correlated (diagonalised) errors:
C
      do j=1,mpar

         jext = iexint(j)

         base = TRIM(OutDirName)//'/pdfs_q2val_'//tag(j)
         idx = index(base,' ')-1

         base2 = TRIM(OutDirName)//'/pdfs_'//tag(j)
         idx2  = index(base2,' ')-1

         do shift_dir=-1,1,2
            if (shift_dir.eq.-1) then
               if (idx.gt.0) then
                  name  = base(1:idx)//'m_'
                  name2 = base2(1:idx2)//'m.lhgrid'
               else
                  name = base//'m_'
                  name2 = base2//'m.lhgrid'
               endif
            else
               if (idx.gt.0) then
                  name = base(1:idx)//'p_'
                  name2 = base2(1:idx2)//'p.lhgrid'
               else
                  name = base//'p_'
                  name2 = base2//'m.lhgrid'
               endif
            endif


C
C Shift parameters by the j-th de-correlated error:
C
            do i=1,npar
               a(i) = pkeep(i)
               iint = iunint(i)
               if (iint.gt.0) then
                  if(doOffset) then
                    shift = shift_dir * DecorVarShift(iint, j)
                  else
                    call MNSTAT(fmin, fedm, errdef, npari, nparx, istat)   !> MW&FG for scaling with DeltaChi2>1.0                                           
                    shift = shift_dir * GetUmat(iint,j)*SQRT(errdef)       !> MW&FG for scaling with DeltaChi2>1.0                                      
                  endif
                  a(i) = a(i) + shift
               endif
            enddo  ! i

            call copy_minuit_extrapars(a) !Set shifted parameters

            ifcncount = ifcncount+1
            chichi = chi2data_theory(2)

            TheoVars(:,(shift_dir+1)/2+1,j) = THEO

C
C Write results out:
C
            open (76,file=name2,status='unknown')
            call store_pdfs(name)
            close (76)
            if(shift_dir.eq.-1) then
              call save_data_lhapdf6(j*2-1)
              call save_data_tmd(j*2-1)
              call error_band_action(j*2-1)
            else
              call save_data_lhapdf6(j*2)
              call save_data_tmd(j*2)
              call error_band_action(j*2)
            endif


         enddo  ! shift_dir

      enddo  ! j

C write out once more (with theory errors filled)

      Theo = TheoFCN3           ! restore
      Theo_mod = TheoModFCN3
      ALPHA_Mod = ALphaModFCN3

      call GetTheoErrorsAsym(TheoVars,ntot,mpar)
      call writefittedpoints

      deallocate(TheoVars)

      return
      end

C-----------------------------------------------------------
C> @brief Compute asymmetric uncertainties for theory predictions based on eigenvector variations, asymmetric hessian
C
C  @param TheoVars 3-D array of theory variations shaped as ndata, 2, nvector
C  @param nd number of data points, ndata
C  @param nv number of eigenvectors, nvector
C-----------------------------------------------------------
      subroutine GetTheoErrorsAsym(TheoVars,nd,nv)
      implicit none
      integer nd,nv
      double precision TheoVars(nd,2,nv)
#include "ntot.inc"
#include "theo.inc"
      integer i,j
      double precision up,dn,eu,ed
C-----------------------------------------
      do i=1,nd
         up = 0.
         dn = 0.
         do j=1,nv
            eu = (max(max(0.,  TheoVars(i,1,j)-Theo(i) )
     $           , TheoVars(i,2,j)-Theo(i)))**2
            ed = (max(max(0., -TheoVars(i,1,j)+Theo(i) )
     $           ,-TheoVars(i,2,j)+Theo(i)))**2

            up = up + eu
            dn = dn + ed
         enddo
         theo_tot_up(i) = sqrt(up)
         theo_tot_down(i) = sqrt(dn)
      enddo

      end

C-----------------------------------------------------------
C> @brief Compute asymmetric uncertainties for theory predictions based on eigenvector variations, symmetric hessian
C
C  @param TheoVars 3-D array of theory variations shaped as ndata, nvector
C  @param nd number of data points, ndata
C  @param nv number of eigenvectors, nvector
C-----------------------------------------------------------
      subroutine GetTheoErrorsSym(TheoVars,nd,nv)
      implicit none
      integer nd,nv
      double precision TheoVars(nd,nv)
#include "ntot.inc"
#include "theo.inc"
      integer i
      double precision e(NTOT)
C-----------------------------------------
      e = 0.0  ! set to 0
      do i=1,nv
         e = e + (TheoVars(:,i)-THEO)**2
      enddo

      theo_tot_up = sqrt(e)
      theo_tot_down = sqrt(e)

      end


!> =================================================
!> Generate error bands for symmetrian hessian case
!> =================================================
      subroutine ErrBandsSym
      implicit none
#include "endmini.inc"
#include "fcn.inc"
#include "steering.inc"
#include "ntot.inc"
#include "systematics.inc"
#include "theo.inc"
      external fcn
      integer icond
      double precision fmin, fedm, errdef
      integer npari, nparx, istat, ifail
      integer i,j, k, ind, ind2, mpar, jext
      double precision, allocatable :: Amat(:,:)
      double precision, allocatable :: eigenvalues(:)
C
      integer  iunint(MNE)  ! internal param. number
      integer  iexint(MNE)  ! external param. number
      double precision
     $     parval, parerr,parlolim,parhilim

      double precision a(MNE)
      integer idx,idx2,iint,kflag

      character *64 parname

      character*48 name,name2
      character*48 base,base2
      character tag(40)*3
      data (tag(i),i=1,40) /'s01','s02','s03','s04','s05',
     +     's06','s07','s08','s09','s10',
     +     's11','s12','s13','s14','s15',
     +     's16','s17','s18','s19','s20',
     +     's21','s22','s23','s24','s25',
     +     's26','s27','s28','s29','s30',
     +     's31','s32','s33','s34','s35',
     +     's36','s37','s38','s39','s40'/

C Function
      double precision GetUmat

      double precision chichi
      double precision chi2data_theory ! function
C for theory errors:
      double precision, allocatable :: TheoVars(:,:)


C------------------------------------------------------------------------

      if (ReadParsFromFile) then
         call ReadPars(ParsFileName,pkeep)
         call MNSTAT(fmin, fedm, errdef, npari, nparx, istat)
      else
C         call MNCOMD(fcn,'SET ERRDEF 9',icond,0)
         call MNCOMD(fcn,'HESSE',icond,0)
C         call MNCOMD(fcn,'ITERATE 10',icond,0)
C     Check the covariance matrix:
         call MNSTAT(fmin, fedm, errdef, npari, nparx, istat)
         print *,'Covariance matrix status =',istat,npari

         if (istat .ne. 3) then
            call hf_errlog(16042702,
     $           'S:Problems with error matrix, can not produce bands')
         endif
      endif

      mpar = 0
      do ind=1,nparx
         call mnpout(ind,parname,parval,parerr,parlolim,
     $        parhilim,iunint(ind))
         if (iunint(ind).gt.0) then
            write (6,*) 'Parameter',ind,' name=',parname
            write (6,*) 'Internal index=',iunint(ind)
            pkeep(ind) = parval
            write (6,*) ' '
            mpar = mpar + 1
         endif
      enddo

      do ind=1,mpar
         do ind2=1,nparx
            if (iunint(ind2).eq.ind) then
               iexint(ind) = ind2
            endif
         enddo
      enddo


      Allocate(Amat(Npari, Npari))
      Allocate(Eigenvalues(Npari))


      if (ReadParsFromFile) then
         call ReadParCovMatrix(CovFileName, Amat, Npari)
      else
         call MNEMAT( Amat, Npari)
      endif

C Diagonalize:
      call MyDSYEVD( Npari, Amat, Npari, Eigenvalues, ifail)

C scale the matirx
      do i=1,npari
         do j=1,npari
            Amat(j,i) = Amat(j,i) * sqrt(Eigenvalues(i))
         enddo
      enddo

      allocate(TheoVars(NTOT,Npari))

C
C Loop over de-correlated errors:
C
      do j=1,Npari
         jext = iexint(j)
         base = TRIM(OutDirName)//'/pdfs_q2val_'//tag(j)
         idx = index(base,' ')-1

         base2 = TRIM(OutDirName)//'/pdfs_'//tag(j)
         idx2  = index(base2,' ')-1

         if (idx.gt.0) then
            name  = base(1:idx)//'s_'
            name2 = base2(1:idx2)//'s.lhgrid'
         else
            name = base//'s_'
            name2 = base2//'s.lhgrid'
         endif


C
C Shift parameters by the j-th de-correlated error:
C
         do i=1,MNE
            a(i) = pkeep(i)
            iint = iunint(i)
            if (iint.gt.0) then
               a(i) = a(i) + Amat(iint,j)
C               a(i) = a(i) + GetUmat(iint,j)
            endif
         enddo

         call copy_minuit_extrapars(a) !Set shifted parameters

         ifcncount = ifcncount+1
         chichi = chi2data_theory(2)  ! sum-rules and evolution are inside

         TheoVars(:,j) = THEO        ! save for error calc.

C
C Write results out:
C
         open (76,file=name2,status='unknown')
         call store_pdfs(name)
         close (76)

         call save_data_lhapdf6(j)
         call save_data_tmd(j)

         call error_band_action(j)
      enddo                     ! j


C write out once more (with theory errors filled)

      Theo = TheoFCN3           ! restore
      Theo_mod = TheoModFCN3
      ALPHA_Mod = ALphaModFCN3

      call GetTheoErrorsSym(TheoVars,ntot,Npari)
      call writefittedpoints

      deallocate(TheoVars)
      end

!> read parameter values from the pars out file
! Probably borken since 2.2.0 --Ivan
      subroutine ReadPars(FileName, pvals)
      implicit none
      character*(*) FileName
      double precision pvals(*)
      integer IStatus

      character*120 buff
      double precision fmin, fedm, errdef
      integer npari, nparx, istat, ind,ii
      character*100 parname
      double precision parval, parerr, parlolim, parhilim
C------------------------------------------------
      print *,'Reading parameter values from '//trim(FileName)
      open (51,file=FileName, status='old',err=3)
 1    read (51,'(A120)',end=2, err=4) buff

      call MNPARS(buff,IStatus)
      goto 1
C Decode

 2    close (51)

      call MNSTAT(fmin, fedm, errdef, npari, nparx, istat)
      do ind=1,nparx
         call mnpout(ind,parname,parval,parerr,parlolim,
     $        parhilim,ii)
         pvals(ind) = parval
      enddo



      return
 3    call hf_errlog(16042810,'F: Can not find parameters file = '
     $     //trim(FileName)//', STOP')
 4    call hf_errlog(16042811,'F: Can not read parameters file = '
     $     //trim(FileName)//', STOP')
      end

      subroutine ReadParCovMatrix( FileName, Cov, Npars)
      implicit none
      character *(*) FileName
      integer NPars
      double precision Cov(Npars,Npars)
      integer i,j
C---------------------------------------------------
      open (51, file=FileName, status='old', err=1)

      print *,npars
      do i=1,NPars
         read (51,*,err=2,end=3) ( Cov(j,i),j=1,NPars )
         print '(20E10.2)' ,( Cov(j,i),j=1,NPars )
      enddo
      print *,'Read covariance matrix from '//trim(FileName)
      close (51)
      return
 1    call hf_errlog(16042820,
     $    'F: Can not open file with parameters covariance matrix '
     $     //trim(FileName)//', STOP')
 2    call hf_errlog(16042821,
     $  'F: Error while reading file with parameters covariance matrix '
     $     //trim(FileName)//', STOP')
 3    call hf_errlog(16042822,
     $    'F: Unexpected end of file with parameters covariance matrix '
     $     //trim(FileName)//', STOP')
      end
