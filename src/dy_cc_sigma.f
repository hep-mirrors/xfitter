      Subroutine GetDYCCXsection(IDataSet)
C-----------------------------------------------------------
C
C Created by SG, 26/05/2011
C Calculate DY W+, W- and asymmetry cross sections
C
C------------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'theo.inc'
      integer IDataSet

      integer NPmax
      parameter(NPmax=100)
      
      logical LFirstTime
      data LFirstTime/.true./
C
C TO BE IMPROVED
C
      double precision dy_mass(2)
      double precision dy_y(2)
      double precision pt_cut,eb(NPmax+1)
      double precision wm_bsigs(NPmax), wp_bsigs(NPmax)


      integer i,idx,idxEta1,idxEta2
      integer idxKfactWplus, idxKfactWminus
      logical LAsymmetry
C Functions:
      integer GetBinIndex
      integer GetInfoIndex
      integer GetKFactIndex
C------------------------------------------------------------

      if (NDATAPOINTS(IDataSet).gt.NPmax) then
         print *,'ERROR IN GetDYCCXsection'
         print *,'INCREASE NPMax to ',NDATAPOINTS(IDataSet)
         stop
      endif

      if (LFirstTime) then
         LFirstTime = .false.

         pt_cut=25.
         dy_mass(1) = 1.
         dy_mass(2) = 7000.

         dy_y(1) = -10.
         dy_y(2) =  10.

C Get indicies:
         idxEta1 = GetBinIndex(IDataSet,'eta1')
         idxEta2 = GetBinIndex(IDataSet,'eta2')

         if (idxEta1.eq.0 .or. idxEta2.eq.0) then
            print 
     $     '(''ERROR in GetDYCCXsection, can not find bin index for Eta1, Eta2'',2i6)'
     $           ,idxEta1,idxEta2
            stop
         endif

C Define bins:
         do i=1,NDATAPOINTS(IDataSet) 
            idx =  DATASETIDX(IDataSet,i)
            eb(i)   =  AbstractBins(idxEta1,idx)
            eb(i+1) =  AbstractBins(idxEta2,idx)
         enddo

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C It seems we will be able to run only one W asymmetry at a time this way!
C
         call w_set_etabins(dy_mass, dy_y, pt_cut, NDATAPOINTS(IDataSet)
     $        , eb)
      endif

C Calculate prediction:
      call w_get_etabins_xs(wm_bsigs, wp_bsigs)

      idxKfactWplus  = GetKFactIndex(IDataSet,'Wplus')
      idxKfactWminus = GetKFactIndex(IDataSet,'Wminus') 

C Apply k-factors:
      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)
         if (idxKfactWminus.gt.0) then
            wm_bsigs(i) = wm_bsigs(i)*Kfactors(idxKfactWminus,idx)
         endif
         if (idxKfactWplus.gt.0) then
            wp_bsigs(i) = wp_bsigs(i)*Kfactors(idxKfactWplus ,idx)
         endif
      enddo

C Check type of the data
      LAsymmetry =  
     $     DATASETInfo( GetInfoIndex(IDataSet,'asymmetry'), IDataSet).gt.0

      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)
         if (LAsymmetry) then
            THEO(idx) = (wp_bsigs(i)-wm_bsigs(i) )/(wp_bsigs(i)+wm_bsigs(i))
         endif
      enddo
      end
