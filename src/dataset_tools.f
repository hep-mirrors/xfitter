      integer Function GetBinIndex(IDataSet,CName)
C-----------------------------------------------------------------
C
C  Created 24/05/2011 by SG
C  Return bin index corresponding to name CName
C
C-----------------------------------------------------------------
      implicit none
      include 'ntot.inc'
c      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
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
      include 'ntot.inc'
c      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
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
      print *,'ERROR: Could not find index for information ',cname

      end


      integer Function GetKFactIndex(IDataSet,CName)
C-----------------------------------------------------------------
C
C  Created 26/05/2011 by SG
C  Return info index corresponding to kfactor name CName
C
C-----------------------------------------------------------------
      implicit none
      include 'ntot.inc'
c      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
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


      Logical Function FailSelectionCuts(reaction,nbin,bins,binnames)
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

C Extra DIS selection
      logical FailDISSelection

      logical LFirst /.true./

      
C Namelist variables:
      integer NProcMax,NRulesMax
      parameter (NProcMax  = 20)  ! Cuts for 20 process types at most
      parameter (NRulesMax = 20)  ! For each process, at most 20 rules
      integer NProcesses          ! actual number of processes
      character*80 ProcessName(NProcMax)         ! names of the processes
      integer NRules(NProcMax)    ! actual number of rules per process
      character *80 Variable(NRulesMax,NProcMax)  ! names of variables to apply cuts
      double precision CutValueMin(NRulesMax,NProcMax)  ! Min. value of the cut
      double precision CutValueMax(NRulesMax,NProcMax)  ! Max. value of the cut

      namelist  /Cuts/NProcesses, NRules,Variable,CutValueMin,CutValueMax,ProcessName

      integer i,j,k
C------------------------------------------------

      FailSelectionCuts = .false. 

      if (LFirst) then
         LFirst = .false.
         print 
     $ '(''First time in FailSelectionCuts. Read the cut definitions'')'

         open (52,file='steering.txt',status='old')
         read (52,NML=Cuts,END=71,ERR=72)
         print '(''Read  cut rules for '',i5,'' processes'')',NProcesses
         close (52)

      endif

C-- Check if the reaction is on the process list
      do i=1,NProcesses
         if (reaction .eq. processname(i)) then
C-- Run over all rules, check for appropriate variable
            do j=1,NRules(i)
               do k=1,nbin
                  if (Variable(j,i).eq.BinNames(k)) then
                     if (
     $                        bins(k).lt.CutValueMin(j,i)
     $                    .or.bins(k).ge.CutValueMax(j,i)) then
                        FailSelectionCuts = .true.  ! Fail Cut
                        Return
                     endif
                     goto 17  ! next rule
                  endif
               enddo
               print 
     $ '(''Error in FailSelectionCuts: variable '',A8,'' not found'')'
     $              ,Variable(j,i)
               print '(''Check reaction '',a16)', reaction
               stop
 17            continue
            enddo
         endif
      enddo

C Extra DIS cuts 
      if (Reaction .eq. 'NC e+-p') then
         FailSelectionCuts = FailDISSelection(nbin,bins,BinNames)
      else 
      endif
      
      Return
 71   continue
      print 
     $  '(''Error in FailSelectionCuts: EOD reading namelist Cuts'')'
      stop
 72   continue
      print 
     $  '(''Error in FailSelectionCuts: Error reading namelist Cuts'')'
      stop

      end

      logical Function FailDISSelection(NBin,Bins,BinNames)
C---------------------------------------------------------------
C
C Created 26/05/11, Apply DIS selection
C
C--------------------------------------------------------------
      implicit none
      integer NBin
      Double Precision Bins(NBin)
      character *(*) BinNames(NBin)
      include 'ntot.inc'
      include 'steering.inc'

      integer idxQ2,idxX,idxY,i
      real*4 q2,x,y

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
         stop
      endif
      q2 = bins(idxQ2)
      X  = bins(idxX)
      Y  = bins(idxY)


      FailDISSelection = FailCuts(0,q2,x,y)

C--------------------------------------------------------------
      end
