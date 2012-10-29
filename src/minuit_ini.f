      subroutine minuit_ini
C======================================================================
C
C 2 Oct 2011: add extra parameters before starting minuit
C WS: modified to have different file names for the offset method
C
C----------------------------------------------------------------------
      implicit none
      include 'steering.inc'

      character*72 MinuitIn
      character*72 MinuitOut
      character*72 MinuitSave      
      character*32 OffsLabel	      
      character*32 Suffix	      
      
      if(CorrSystByOffset) then
        Suffix='_'//OffsLabel(CorSysIndex,'.txt')
      else
        Suffix = '.txt'
      endif
      
      MinuitOut = 'output/minuit.out'//Suffix
      open ( 25, file=MinuitOut )

      MinuitIn='minuit.in.txt' 
      write(6,*) ' read minuit input params from file ',MinuitIn
      call HF_errlog(12020504,
     +     'I: read minuit input params from file '//MinuitIn) 
      open ( 24, file=MinuitIn )

      MinuitSave = 'output/minuit.save'//Suffix
      open (  7, file=MinuitSave)
      
      call mintio(24,25,7)
      
      return
      end

      
      subroutine ExtraParam
C======================================================================
C
C MINUIT module 'minuit.F' is modified
C to call ExtraParam after reading parameters
C
C----------------------------------------------------------------------
      implicit none
      include 'extrapars.inc'
      integer i, ierrf

C Add extra parameter:

      do i = 1,nExtraParam
         call mnparm(100+i,ExtraParamNames(i)
     $        ,ExtraParamValue(i)
     $        ,ExtraParamStep(i)
     $        ,ExtraParamMin(i)
     $        ,ExtraParamMax(i)
     $        ,ierrf)
         if (ierrf.ne.0) then
            print *,'Error adding extra parameter',i
            print *,'name, value, step, min, max are:',
     $           ExtraParamNames(i)
     $        ,ExtraParamValue(i)
     $        ,ExtraParamStep(i)
     $        ,ExtraParamMin(i)
     $        ,ExtraParamMax(i)
            print *,'Error code=',ierrf
            call HF_errlog(12020505,'F: Error in ExtraParam')
         else
            iExtraParamMinuit(i) = 100+i
         endif
      enddo
      end

      
C ==================================================
      CHARACTER*(*) FUNCTION OffsLabel(mu, tail)
C
C 2012-06-25 ws
C Generate trailing part of an output file name
C
C---------------------------------------------------
      implicit none
      integer mu
      character*(*) tail
      integer amu,smu,ndig
      parameter(ndig=3) ! as NSYSMAX = 300
      character str*(ndig)
      character*8 fmt
      character*2 mp
      parameter(mp = 'mp')
      
      if(mu.eq.0) then
        OffsLabel = '0' // tail
        return
      end if
      
      smu = (ISIGN(1,mu)+1)/2 +1
      amu = IABS(mu)
      write(fmt,'(a,i1,a,i1,a)') '(i',ndig,'.',ndig,')'
      write(str,fmt) amu
      OffsLabel = str // mp(smu:smu) // tail
      return
      end

c ===============================================
      Subroutine Do_Fit

      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'
      include 'systematics.inc'
      include 'indata.inc'
      include 'for_debug.inc'
      ! include 'theo.inc'
      include 'fcn.inc'
      ! include 'endmini.inc'
     
      external fcn

      integer amu
      double precision daten_notshifted(NTOT)
      double precision musign

      character*72 ResultsFile
      character*16 OffsLabel

      integer j
      logical file_exists
      ! .......................................

       ResultsFile = 'output/Results'
       if(CorrSystByOffset) then
         ResultsFile = TRIM(ResultsFile)//'_' // OffsLabel(CorSysIndex,'.txt')
       else
         ResultsFile = TRIM(ResultsFile)//'.txt'
       endif
       
       if(UsePrevFit .eq. 2) then
         INQUIRE(FILE=ResultsFile, EXIST = file_exists)
         if(file_exists) return
       endif
       
       if(CorrSystByOffset .and. CorSysIndex .ne. 0) then
         ! --- store original data
         do j=1,npoints
            daten_notshifted(j) = daten(j)
         enddo
         ! --- shift by CorSysIndex
         amu = IABS(CorSysIndex)
         ! write(*,*)'--> OFFSET CorrSrc = ',CorSysIndex
         musign = ISIGN(1,CorSysIndex)
         do j=1,npoints
               ! write(*,*)'npoints ',j,DATEN (j), beta(CorSysIndex,j)
               daten (j) = 
     &         daten_notshifted(j)*(1 + musign*beta(amu,j))
         enddo
       endif

*     ------------------------------------------------
*     Do the fit
*     ------------------------------------------------
       OPEN(85,file=ResultsFile,form='formatted',status='replace')
       write(*,*) 'ResultsFile: ',ResultsFile
       call minuit_ini  ! opens Minuit i/o files
cws       lprint = .true.
       lprint = .false.
       call minuit(fcn,0)       
       close(85)
*     ------------------------------------------------
       
       if(CorrSystByOffset) then
         ! --- these units are used by MNCOMD for DOBANDS
         close(24)
         close(25)
         close(7)
         call Offset_SaveParams(CorSysIndex)
         if(CorSysIndex .eq. 0) then
           call Offset_SaveStatCov
           call write_pars(0)        ! Write out central parameters
         else
           ! --- restore original (unshifted) data
           do j=1,npoints
              daten(j) = daten_notshifted(j)
           enddo
         endif
       endif
         
      return
      end
