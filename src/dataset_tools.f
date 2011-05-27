      integer Function GetBinIndex(IDataSet,CName)
C-----------------------------------------------------------------
C
C  Created 24/05/2011 by SG
C  Return bin index corresponding to name CName
C
C-----------------------------------------------------------------
      implicit none
      include 'steering.inc'
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
      print *,'ERROR: Could not find index for variable',cname

      end

      integer Function GetInfoIndex(IDataSet,CName)
C-----------------------------------------------------------------
C
C  Created 24/05/2011 by SG
C  Return info index corresponding to name CName
C
C-----------------------------------------------------------------
      implicit none
      include 'steering.inc'
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
      include 'steering.inc'
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
      integer nbin
      double precision bins(nbin)
      character *(*) reaction,binnames(nbin)

      logical FailDISSelection
C------------------------------------------------

      FailSelectionCuts = .false. 

      if (Reaction .eq. 'NC e+-p') then
         FailSelectionCuts = FailDISSelection(nbin,bins,BinNames)
      else 
      endif
      
      
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
