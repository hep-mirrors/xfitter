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
         write(90,*)NQ2BINS(i)
         index = 0
         do j=1,NQ2BINS(i)
            index = index + NXBINS(i,j)
            write(90,*)VAL_Q2(i,index),NXBINS(i,j)
         enddo
         write(90,*) '     q2          x        y    data     +- uncorr.err'//
     &        '   +-toterr      theory      pull     dataset'
         do j=1,NDATAPOINTS(i)
            index = DATASETIDX(i,j)
            write(90,'(1X,8(e11.5,1X),i4)') VQ2(index),VX(index), VY(index),
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
