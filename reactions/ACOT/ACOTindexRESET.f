      subroutine ACOTindexRESET()
c     
c     ACOT uses a k-factor to speed up the calculations.
c     Each call to the ACOT code increments index
c     which is used to identify each data point,
c     and this is used to store a k-factor in a local array.
c     
c     index is reset by this function when parameters are updated
c     and we start a new chi2 calculation.
c     We assume data are always called in the same order
c      
      implicit none
      integer index
      common /acotIndex/  index

      index=0
      return
      end
