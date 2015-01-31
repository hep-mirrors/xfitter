	subroutine GetGridkt
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "fcn.inc"
#include "qcdnumhelper.inc"
      

      integer IDataSet
      integer idxQ2, idxX, idxY, i,  idx
      
      double precision X(NPMaxDIS),Y(NPMaxDIS),Q2(NPMaxDIS),XSec(NPMaxDIS)
C Functions:
      integer GetBinIndex
      integer GetInfoIndex

      logical IsReduced
      
      Real RVQ2(ntot),xx,q2x

      Real sq2,xmymin
      Integer nq2point
      Common /myq2grid/ nq2point,sq2(200),xmymin
      
      Integer Output(ntot),j
      Real q2old
      Logical first
      Data first/.true./
      
      if(first) then 
C
C Get indexes for Q2, x and y bins:
C
      xmymin = 1.
      j=0
      do Idataset = 1,NDATASETS
         idxQ2 = GetBinIndex(IDataSet,'Q2')
         idxX  = GetBinIndex(IDataSet,'x')
         idxY = GetBinIndex(IDataSet,'y')
c      IsReduced = DATASETInfo( GetInfoIndex(IDataSet,'reduced'), IDataSet).gt.0
          if (idxQ2.eq.0 .or. idxX.eq.0) then
            Return
         endif
         do i=1,NDATAPOINTS(IDataSet)
      
C
C Reference from the dataset to a global data index:
C
         idx =  DATASETIDX(IDataSet,i)
C
C Local X,Y,Q2 arrays, used for QCDNUM SF caclulations:
C
         X(i)   = AbstractBins(idxX,idx)
         if(idxY.ne.0) Y(i)   = AbstractBins(idxY,idx)
         Q2(i)  = AbstractBins(idxQ2,idx)
         xx = x(i)
         q2x = q2(i)
c         write(6,*) ' getgridkt x,q2 ',xx,q2x,i,idx
         j=j+1
         RVQ2(j) = real(q2(i))
         if(x(i).LE.xmymin) xmymin=x(i)

         End do
      End do
      
      
      call SORTRX(npoints,RVQ2,Output)

      write(6,*) ' Getgridkt: q2 ordered array for datasets',NDATASETS
      q2old = 0
      j=0
      do i=1,npoints

      if(rvq2(output(i)).ne.q2old) then
        q2old = rvq2(output(i))
        j=j+1
        if(j.gt.200) then
          write(6,*) ' FATAL in getnxxskt: j>200 '
          write(6,*) ' program stop '
          stop
        endif
        sq2(j)= real(rvq2(output(i)))
      endif
c      write(6,*) ' i = ',i,' q2 = ', rvq2(output(i)),output(i)
      end do
      nq2point = j
c      write(6,*) ' getgridkt nq2point = ',nq2point

      first=.false.
      endif 
      
      Return
      End
