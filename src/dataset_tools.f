      integer Function GetBinIndex(IDataSet,CName)
C-----------------------------------------------------------------
C
C  Created 24/05/2011 by SG
C  Return bin index corresponding to name CName
C
C-----------------------------------------------------------------
      implicit none
#include "ntot.inc"
c#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
      integer IDataSet
      character *(*) CName
      integer i
C----------------------------------------------------------------

      do i=1,DATASETBinningDimension(IDataSet)
         if (CName.eq. DATASETBinNames(i,IDataSet)) then
            GetBinIndex = i
            Return
         endif
      enddo
      
      GetBinIndex = 0
c      print *,'ERROR: Could not find index for variable',cname

      end

      integer Function GetInfoIndex(IDataSet,CName)
C-----------------------------------------------------------------
C
C  Created 24/05/2011 by SG
C  Return info index corresponding to name CName
C
C-----------------------------------------------------------------
      implicit none
#include "ntot.inc"
c#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
      integer IDataSet
      character *(*) CName
      integer i
C----------------------------------------------------------------

      do i=1,DATASETInfoDimension(IDataSet)
         if (CName.eq. DATASETInfoNames(i,IDataSet)) then
            GetInfoIndex = i
            Return
         endif
      enddo
      
      GetInfoIndex = 0
C      print *,'ERROR: Could not find index for information ',cname

      end

      integer Function GetParameterIndex(CName)
!>
!> Created 24 Sept 2011
!>
!>
      implicit none
#include "extrapars.inc"
      character*(*) CName
      integer i
C--------------------------------------------------
      do i=1,nExtraParam
         if (ExtraParamNames(i).eq.cName) then
            GetParameterIndex = i
            Return
         endif
      enddo

C Not found
      GetParameterIndex = 0
C--------------------------------------------------
      end

! [--- WS 2015-10-04
! =========================================================
!> \brief Get current value of an extra-parameter by name.
!> \details During the fit reads current value stored in \c parminuitsave.
      double precision function XParValueByName(Pname)
        implicit none
        character*(*) Pname
#include "endmini.inc"
        ! for parminuitsave
#include "extrapars.inc"
#include "fcn.inc"
        integer idx
        ! functions:
        integer GetParameterIndex

        ! 1-based index of Pname in ExtraParams
        idx = GetParameterIndex(Pname)
        if (idx.eq.0) then
           print *,'Minuit Extra parameter "',Pname,'" not defined.'
           call HF_stop
        endif
        if(NparFCN.gt.0) then
          ! iExtraParamMinuit gives global Minuit index
          XParValueByName = parminuitsave(iExtraParamMinuit(idx))
        else
          XParValueByName = ExtraParamValue(idx)
        endif
      end
! ---]


      integer Function GetKFactIndex(IDataSet,CName)
C-----------------------------------------------------------------
C
C  Created 26/05/2011 by SG
C  Return info index corresponding to kfactor name CName
C
C-----------------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "for_debug.inc"
#include "datasets.inc"
      integer IDataSet
      character *(*) CName
      integer i
C----------------------------------------------------------------

      do i=1,DATASETNKfactors(IDataSet)
         if (CName.eq. DATASETKFactorNames(i,IDataSet)) then
            GetKFactIndex = i
            Return
         endif
      enddo
      
      GetKFactIndex = 0
      print *,'ERROR: Could not find index for kfactor ',cname

      end


      Logical Function FailSelectionCuts(reaction,nbin,bins,binnames,IndexDataset)
C----------------------------------------------------------------------
C  Created 24/05/11 by SG.
C
C  Apply selection cuts for reaction
C
C----------------------------------------------------------------------
      integer nbin                                  ! number of variables
      double precision bins(nbin)                   ! values of variables
      character *(*) reaction                       ! reaction type
      character *(*) binnames(nbin)                 ! variable names 
      integer IndexDataset                          ! dataset index

C Extra DIS selection
      logical FailDISSelection

      logical LFirst /.true./

      
C Namelist variables:
      integer NRulesMax
      parameter (NRulesMax = 100)               !> Maximum number of cut-rules
      integer NRules                            !> actual number of rules per process
      character*80 ProcessName(NRulesMax)       !> names of the process
      character *80 Variable(NRulesMax)         !> names of variables to apply cuts
      double precision CutValueMin(NRulesMax)   !> Min. value of the cut
      double precision CutValueMax(NRulesMax)   !> Max. value of the cut
      integer NDatasetMax
      parameter (NDatasetMax = 10)
      integer NDataset(NRulesMax)               !> actual number of provided dataset indices
      save NDataset                             
      integer Dataset(NDatasetMax,NRulesMax)    !> dataset indices

      namelist  /Cuts/Variable,CutValueMin,CutValueMax,ProcessName,Dataset

      integer i,j,k
      logical CheckThisDataSet
C------------------------------------------------

      FailSelectionCuts = .false. 

      if (LFirst) then
         LFirst = .false.
         print 
     $'('' First time in FailSelectionCuts. Read the cut definitions'')'

         do i=1,NRulesMax
            ProcessName(i) = ' '
            Variable(i)    = ' '
         enddo

         open (52,file='steering.txt',status='old')
         read (52,NML=Cuts,END=71,ERR=72)
         close (52)
C Count rules
         NRules = 0
         do i=1,NRulesMax
            if (ProcessName(i).ne.' ') then
               NRules = NRules + 1
               
C Count number of dataset indices   
               NDataset(i)  = 0
               do j=1,NDatasetMax
                 if(Dataset(j,i).eq.0) then
                   exit
                 else
                   NDataset(i) = NDataset(i) + 1            
                 endif
               enddo
               
               print ' ("ProcessName  ",A80)',ProcessName(i)
               print ' ("Variable:    ",A80)',Variable(i)
               print ' ("CutValueMin: ",E15.5)',CutValueMin(i)
               print ' ("CutValueMax: ",E15.5)',CutValueMax(i)
               do j=1,NDataset(i)
                 print ' ("Dataset(",I0,"): ",I6)',j,Dataset(j,i)
               enddo
               
            endif
         enddo

         print '('' Read '',i3,'' cut rules'')',NRules
      endif

C-- Run over all rules, check for appropriate process/variable
      do j=1,NRules
         if (reaction .eq. processname(j)) then
            if (Variable(j).eq.'Whad2') then
               FailSelectionCuts = FailDISSelection(nbin,bins,
     $              BinNames,CutValueMin(j))
               if (FailSelectionCuts.EQV..true.) then
                  Return
               endif
               goto 17          ! next rule               
            else

C Check dataset indices, if they were provided     
               CheckThisDataSet = .true.  ! by default check all datasets
               if(NDataset(j).gt.0) then
                 CheckThisDataSet = .false.  ! explicit list was provided, now default is false
                 do k=1,NDataset(j)
                   if(DataSet(k,j).eq.IndexDataset) then
                     CheckThisDataSet = .true.  ! item matches current datset, stop loop
                     exit
                   endif
                 enddo
               endif
               if(.not.CheckThisDataSet) then
                 goto 17  ! next rule
               endif

               do k=1,nbin
                  if (Variable(j).eq.BinNames(k)) then
                     if (
     $                    bins(k).lt.CutValueMin(j)
     $                    .or.bins(k).gt.CutValueMax(j)) then
                        FailSelectionCuts = .true. ! Fail Cut
                        Return
                     endif
                     goto 17    ! next rule
                  endif
               enddo
            endif
            call HF_errlog(12032201,
     $           'I: FailSelectionCuts can not cut on requested variable')
c            print 
c     $'(''Warning in FailSelectionCuts: variable '',A8,'' not found'')'
c     $           ,Variable(j)
c            print '(''Check reaction '',a16)', reaction
c     call HF_stop
 17            continue
         endif
      enddo
      
      Return
 71   continue
      print 
     $  '(''Error in FailSelectionCuts: EOD reading namelist Cuts'')'
      call HF_stop
 72   continue
      print 
     $  '(''Error in FailSelectionCuts: Error reading namelist Cuts'')'
      call HF_stop

      end

      logical Function FailDISSelection(NBin,Bins,BinNames,CutMin)
C---------------------------------------------------------------
C
C Created 26/05/11, Apply DIS selection
C
C--------------------------------------------------------------
      implicit none
      integer NBin
      Double Precision Bins(NBin)
      double precision CutMin
      character *(*) BinNames(NBin)
#include "ntot.inc"
#include "steering.inc"

      integer idxQ2,idxX,idxY,i
      real*4 q2,x,y, cut

      logical FailCuts
C---------------------------------------------------------------
      idxQ2 = 0
      idxX  = 0
      idxY  = 0

      do i=1,NBin
         if (BinNames(i).eq.'Q2') then
            idxQ2 = i
         endif
         if (BinNames(i).eq.'x') then
            idxX = i
         endif
         if (BinNames(i).eq.'y') then
            idxY = i
         endif
      enddo
      if (idxQ2.eq.0 .or. idxX.eq.0 .or. idxY.eq.0) then
         print 
     $    '(''ERROR in FAIL DIS SELECTION: missing, q2, x or Y'',3I4)'    
     $        ,idxQ2,idxX,idxY
         call HF_stop
      endif
      q2 = bins(idxQ2)
      X  = bins(idxX)
      Y  = bins(idxY)
      cut = Real (CutMin)

      FailDISSelection = FailCuts(q2,x,y,cut)

C--------------------------------------------------------------
      end


      logical function FailCuts(q2,x,y,cut)
C
C 26/07/2010: added pq2max cut
C
#include "steering.inc"
      real*4 q2,x,y, whad2, cut
      logical fail

      fail = .false.


C 24 Aug 2010: Add saturation inspired cut
      if (q2 .lt. asatur * x**lsatur ) then
         fail = .true.
         print *,'Fail saturation cut',x,q2
      endif

C 28 Oct 2010: fixed target data get out of higher twist...etc
      whad2=q2/x-q2+0.938
      if (whad2<cut) then
         fail=.true.
         print *, 'Failed Whad2 cut', x, q2, whad2
      endif

      FailCuts = fail

      return
      end

