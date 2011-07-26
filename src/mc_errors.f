
      Subroutine MC_method
C------------------------------------------------------------
C
C MC method for propagating of the data uncertainties. Creat a replica
C of the data, which fluctuates accoding to their uncertainteis.
C
C------------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'datasets.inc'      
      include 'ntot.inc'
      include 'systematics.inc'
      include 'indata.inc'
      include 'for_debug.inc'

C To be used as a seed:
      integer icount
      integer isys,ntime,ndate
C Common from CERNLIB datime:
      integer IS
      common/SLATE/ IS(40)

      integer n0
      double precision s,voica,dummy_st

C Single precision here:
      real rndsh,ranflat
C 
      double precision rand_shift(NSYS)
      double precision r_sh_fl(NSYS)
      double precision f_un
      parameter (f_un = 2.0)   ! translate 0.:-1 to -1.:1. 
C For log normal random shifts:
      real lsig, lmu,lrunif
C functions:
      real logshift
      double precision alnorm
C------------------------------------------------------------

      if (iseedmc.ne.0) then
C Seed from the steering:
         icount = iseedmc         
      else
C Seed from current time
         call datime(ndate,ntime)
         ntime = ntime*1000000+is(6)
         icount=ntime
      endif
      

cv initialise the random shifts
      do isys=1,nsys
         rand_shift(isys)=0.
         r_sh_fl(isys)=0.
      enddo

cv initialise the random seed gener

      call rmarin(icount,0,0)
      call rluxgo(3,icount,0,0)
      print*,'initialize smeering with a seed isdrn = ',icount

C
C Loop over systematic sources:
C         
      do isys=1,nsys
         call rnorml(rndsh,1)   ! gauss random number
         call ranlux(ranflat,1) ! uniform random number

         rand_shift(isys) = rndsh
         r_sh_fl(isys) = (ranflat-0.5)*f_un

         print '(''random numbers: sys, gauss, flat '',2i6,2F8.2)',
     $        isys,isys, rand_shift(isys),
     $        r_sh_fl(isys)
      enddo

C
C Loop over the data:
C
      do n0=1,npoints
         call rnorml(rndsh,1)   
         call ranlux(ranflat,1)

         s = DATEN(n0)

         do isys=1,nsys

cv  test different distributions
cv  first for systematic uncert, then for stat.                    

            if (systype.eq.1) then ! gauss
               s = s*(1.+ beta(isys,n0) * rand_shift(isys))
               
            elseif (systype.eq.2) then ! uniform
               s = s*(1. + beta(isys,n0) * r_sh_fl(isys))
               
            elseif (systype.eq.3) then ! lognormal
               if (beta(isys,n0).ne.0) then
                  lsig=beta(isys,n0) 
                  lmu=1.
                  lrunif=r_sh_fl(isys)
                  s=s*logshift(lmu,lsig,lrunif)
                  print*,'log...', n0,isys,
     $                 lrunif, beta(isys,n0), 
     $                 s,logshift(lmu,lsig,lrunif)
               endif
            endif               ! endif (sys for systematic shifts)
         enddo                  ! end loop over the systematic shifts
               
         voica=s                ! save cross section before the stat shift

CV now choose sta (advised gauss)  
              
         if (statype.eq.1) then ! gauss
            s = s + rndsh * alpha(n0)
         elseif (statype.eq.3.) then ! lognormal
            dummy_st=alpha(n0)
            s= s*alnorm(1.,dummy_st)
         elseif (statype.eq.2) then ! uniform
            s = s + alpha(n0)*f_un*(ranflat-0.5)
         endif
         
 
         print 
     $ '(''Original, systematics and stat. shifted data:'',i4,3E12.4)'
     $        , n0,DATEN(n0), voica,s
         DATEN(n0) = s
      enddo   

C------------------------------------------------------------
      end



*     ---------------------------------------------
      function alnorm(am,as)
*     ---------------------------------------------
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



*     -------------------------------------------
        real function logshift(mmu,ssig,rrunif) 
*     -------------------------------------------
c        real function logshift(mu,sig,runif) 
C-----------------------------------------------------------------------
C-
C-   Purpose and Methods: 
C-
C-   Inputs  :
C-   Outputs :
C-   Controls:
C-
C-   Created  12-JUN-2008   Voica Radescu
C-
C-----------------------------------------------------------------------
      IMPLICIT NONE

      
      real zeroth, ANS,ex2,runif
      real mu, sig,x2, mu2, sig2,z1,z2
      external zeroth
      real mmu,rrunif, ssig
c      integer iseed,isdrn,i


      COMMON/PARAM/mu,sig,runif

      mu=mmu
      sig=ssig
      runif=rrunif
c      print*,'logshift', mmu, ssig,rrunif

c      call system_clock(iseed)
c      isdrn=iseed
c      call rmarin(isdrn,0,0)
c      call rluxgo(3,isdrn,0,0)
c      print*,'initialize smeering with a seed isdrn',isdrn

c      do i=1,1000
c         call ranlux(runif,1)
c         print*,'uniform random numbers   ', runif

c         z1 = zeroth(0.0001)
c         z2 = zeroth(10000.)

c      print *,'zeroth(0)=',z1,  ' zeroth(10000.)',z2


c   CALL tZERO(A_in,B_in,X0_out,R_out,EPS_in,MAXF_in,F_in)
      call rzero(0.0001,10000.,x2,ex2,0.0001,10000,zeroth)
c         print *,'result:',x2,ex2
c         z1 = zeroth(x2)
      logshift=x2


 



C-----------------------------------------------------------------------
 999  CONTINUE
      END


*     -------------------------------------------
      REAL FUNCTION ZEROTH(XX,i)
C-----------------------------------------------------------------------
C-
C-   Purpose and Methods: 
C-
C-   Inputs  :
C-   Outputs :
C-   Controls:
C-
C-   Created  12-JUN-2008   Voica Radescu
C-
C-----------------------------------------------------------------------
      IMPLICIT NONE

      real runif,erf
      integer i, isdrn, iseed
      real mmu, ssig
      real xx, mu, sig,stdlog,amlog
      real zeroth1

      COMMON/PARAM/mu,sig,runif

cv transform the formula from mean, std of x to log(x) 
      stdlog = sqrt(log(1+(sig/mu)**2 ) )
      amlog  = log(mu) - 0.5 * log(1+(sig/mu)**2)

      zeroth1=0.5+0.5*erf((log(xx)-amlog)/stdlog/1.4142)
c      if (sig.gt.0) then
      zeroth=zeroth1-runif
c     else
c         zeroth = (1-zeroth1)-runif
c      endif
C-----------------------------------------------------------------------
 999  RETURN
      END


*     -------------------------------------------




