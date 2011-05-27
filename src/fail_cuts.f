
      logical function FailCuts(iset,q2,x,y)
C
C 26/07/2010: added pq2max cut
C
      include 'steering.inc'
      real*4 q2,x,y, whad2
      logical fail

      fail = .false.

      if (pxmax.gt.0.0.and.x.gt.pxmax) fail = .true.
      if (pxmin.gt.0.0.and.x.lt.pxmin) fail = .true.
      if ( (pq2min.gt.0.0.and.q2.lt.pq2min).or.
     &     (pq2max.gt.0.0.and.q2.gt.pq2max) )
     & fail = .true.

C 24 Aug 2010: Add saturation inspired cut
      if (q2 .lt. asatur * x**lsatur ) then
         fail = .true.
         print *,'Fail saturation cut',x,q2
      endif

C 28 Oct 2010: fixed target data get out of higher twist...etc
      whad2=q2/x-q2+0.938
      if (whad2<15) then
         fail=.true.
         print *, 'Failed Whad2 cut', x, q2, whad2
      endif

      FailCuts = fail

      return
      end

