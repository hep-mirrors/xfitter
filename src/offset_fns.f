c ===============================================
      Subroutine Offset_SaveStatCov
      implicit none
#include "fcn.inc"
#include "steering.inc"
      integer j,k
c      integer nparFCN ! number of fit parameters
      ! character*16 OffsLabel
      character*300 OutFile
      double precision Cov(nparFCN,nparFCN)
      do j=1,nparFCN
        do k=1,nparFCN
          Cov(j,k) = 0.d0
        enddo
      enddo
      call MNEMAT(Cov,nparFCN)
      ! write(*,*) '------- Stat. covariance matrix -------'
      ! write(*,*) Cov
      OutFile = TRIM(OutDirName)//'/statcov_0.txt'
      print *,' '
      print *,'==> Saving covariance matrix to ',Trim(OutFile)
      call flush(6)
      OPEN(86,file=OutFile,form='formatted',status='replace')
      write(86,*) nparFCN
      write(86,*) Cov
      CLOSE(86)
      
      return
      end

c ===============================================
      Subroutine Offset_SaveParams(mu)
      implicit none
#include "fcn.inc"
#include "steering.inc"
      integer mu
      double precision VAL,ERR,XLOLIM,XUPLIM
      integer IU,IUINT
      character*16 chnam
      character*16 OffsLabel
      character*64 OutFile
      double precision vp(nparFCN)
      do iu=1,nparFCN
        call MNPOUT(-iu,CHNAM,VAL,ERR,XLOLIM,XUPLIM,IUINT)
        vp(iu) = val
      enddo
      ! write(*,*) '------- Fitted parameters ',mu, ' -------'
      ! write(*,*) vp
      OutFile = TRIM(OutDirName)//'/params_' // OffsLabel(mu,'.txt')
      OPEN(86,file=OutFile,form='formatted',status='replace')
      write(86,*) nparFCN
      write(86,*) vp
      CLOSE(86)

      return
      end


c ===============================================
      Subroutine Offset_Finalize(iErr)
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "systematics.inc"
      integer iErr
      integer OffsetCollect
      integer nOffset
      
       ! Collect results from all Offset fits
       ! ------------------------------------
c       CALL GetNOffset(nOffset)
       iErr = OffsetCollect(TRIM(OutDirName)//CHAR(0))
       if(iErr.eq.1) then
         ! print *,'WARNING, RC=',iErr,' collecting offset results.'
         print *,'WARNING, still not all Offset method results available.'
         call hf_errlog(12110210,
     +    'W: not enough data to finalize Offset calculation')
         ! call HF_stop
         return
       elseif(iErr.eq.2) then
         print *,'WARNING, bad data in some Offset method results.'
         call hf_errlog(12110211,
     +    'W: bad data in some Offset method results')
         return
       endif
       call hf_errlog(12102801,
     +    'I: full Offset method results saved to offset.save.txt')

      return
      end

c ===============================================
      Subroutine RecovCentrPars
      implicit none
      ! include 'ntot.inc'
#include "steering.inc"
#include "iofnames.inc"
#include "for_debug.inc"
      ! include 'systematics.inc'
      ! integer iErr
      character*32 Suffix
      external fcn
      
        Suffix = '.txt'
      
        ! Recover last central fit results
        ! ------------------------------------
        CorSysIndex = 0
        MinuitIn='minuit.temp.in.txt'
        ! --- prepare minuit input for the final run
        Call RecoverParams(TRIM(OutDirName), MinuitIn)
        ResultsFile = TRIM(OutDirName)//'/Results'//Suffix
        MinuitOut = TRIM(OutDirName)//'/minuit.out'//Suffix
        MinuitSave = TRIM(OutDirName)//'/minuit.save'//Suffix
        call minuit_ini  ! opens Minuit i/o files
        lprint = .true.
        ! lprint = .false.
        call minuit(fcn,0)
        close(85)
        close(24)
        close(25)
        close(7)

      return
      end

