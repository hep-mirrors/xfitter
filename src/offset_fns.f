c ===============================================
      Subroutine Offset_SaveStatCov
      implicit none
      include 'fcn.inc'
c      integer nparFCN ! number of fit parameters
      character*16 OffsLabel
      character*64 OutFile
      double precision Cov(nparFCN,nparFCN)
      call MNEMAT (Cov,nparFCN)
      write(*,*) '------- Stat. covariance matrix -------'
      write(*,*) Cov
      OutFile = 'output/statcov_0.txt'
      OPEN(86,file=OutFile,form='formatted',status='replace')
      write(86,*) nparFCN
      write(86,*) Cov
      CLOSE(86)
      
      return
      end

c ===============================================
      Subroutine Offset_SaveParams(mu)
      implicit none
      include 'fcn.inc'
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
      ! call OffsetNewParams(vp, mu)
      write(*,*) '------- Fitted parameters ',mu, ' -------'
      write(*,*) vp
      OutFile = 'output/params_' // OffsLabel(mu,'.txt')
      OPEN(86,file=OutFile,form='formatted',status='replace')
      write(86,*) nparFCN
      write(86,*) vp
      CLOSE(86)

      return
      end


c ===============================================
      Subroutine Offset_Finalize
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'systematics.inc'
      integer icond
      integer OffsetCollect
      
        ! Collect results from all Offset fits
        ! ------------------------------------
        ! --- nSys is fixed by read_data
        ! --- and stored in common/systema/ (systematics.inc)
        icond = OffsetCollect(nSys, 'output'//CHAR(0))
        if(icond.gt.0) then
          print *,'ERROR ',icond,' collecting offset results.'
          call HF_stop
        endif
        call hf_errlog(12102801,
     +     'I: full Offset method results saved to offset.save.txt')
        if (DOBANDS) then
          ! fill params and umat
          call hf_errlog(12102802,
     +     'W: error bands not yet implemented for Offset method')
          ! call Error_Bands_Pumplin   
        endif

      return
      end
