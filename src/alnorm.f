      function alnorm(am,as)
      implicit none
cv      Program voica
C---------------------------------------------------
C Created by SG, 23 Apr 2008 following
C
C http://www.brighton-webs.co.uk/distributions/lognormal.asp
C
C Log-normal random number generator
C
C Input:  am -- mean value 
C         as -- RMS
C----------------------------------------------------
      real pi
      parameter (pi=3.141592)

      real am,as
      integer ntime, ndate, isrnd,is
      real normrnd1,normrnd2
      real stdlog, amlog,r1,r2,rr, alnorm
      COMMON/SLATE/IS(40)
cv     am=1
cv      as=1
 
C SG: Comment out initialization of the seed, already done in read_data !
Csg      call datime(ndate,ntime)
Csg      ntime = ntime*100+is(6)
Csg      isrnd = ntime
      
Csg      call rmarin(isrnd,0,0)
cv      call rnorml(normrnd1,1)
cv      call rnorml(normrnd2,1)
      call ranmar(normrnd1,1)
      call ranmar(normrnd2,1)	

cv      r1 = rand()
cv      r2 = rand()

      r1 = normrnd1
      r2 = normrnd2
      
      rr = sqrt(-2*log(r1))*sin(2*pi*r2)

      stdlog = sqrt(log(1+(as/am)**2 ) )
      amlog  = log(am) - 0.5 * log(1+(as/am)**2)
 

cv      stdlog=0.548662
cv      amlog =-0.150515
      rr = amlog + rr*stdlog

      alnorm = dble(exp(rr))

      
cv      print*,'voica gets the lognorml distribution....',alnorm
      end
