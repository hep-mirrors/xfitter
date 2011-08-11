      SUBROUTINE WRITEFITTEDPOINTS

      implicit none
      
      include 'steering.inc'
      include 'ntot.inc'
      include 'datasets.inc'
      INCLUDE 'indata.inc'
      include 'systematics.inc'
      INCLUDE 'theo.inc'
      
      integer i,j,index

      open(90,file='output/fittedresults.txt')
      write(90,*)ndatasets


      do i=1,ndatasets
         write(90,*)DATASETNUMBER(i)
         write(90,*)DATASETLABEL(i)
c         write(90,*) '     q2          x        y    data     +- uncorr.err'//
c     &        '   +-toterr      theory      pull     dataset'

         write (90,17) (DATASETBinNames(j,i),j=1,3),'data    '
     $        ,' +- uncor  ',' +- tot   ',' theory   ', ' pull   ', 'iset'
 17      format(1X,8(A11,1X),A4)

         do j=1,NDATAPOINTS(i)
            index = DATASETIDX(i,j)
            write(90,'(1X,8(e11.5,1X),i4)') 
     $              AbstractBins(1,index),
     $              AbstractBins(2,index),AbstractBins(3,index),
     &           DATEN(index),ALPHA(index),
     &           E_TOT(index)/100.*DATEN(index),THEO(index),
     &           (DATEN(index)-THEO(index))/ALPHA(index),
     &           DATASETNUMBER(i)
cv
c            write(44,111) VQ2(index),VX(index), f2sh(index),flsh(index),
c     &           xf3sh(index)
         enddo
cv         write(34,*), index,i,DATASETNUMBER(i)
      enddo
  111  format(1X, F10.3, 2X, F12.6, 2X, 3(F12.6,2X))
      close(90)
   
      
      RETURN
      END
