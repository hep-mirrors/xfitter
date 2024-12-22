C------------------------------------------------------------
C
!> MC method for propagating of the data uncertainties.
!> Creat a replica of the data, which fluctuates accoding to their uncertainteis.
C
C------------------------------------------------------------
      Subroutine MC_method()

      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "systematics.inc"
#include "indata.inc"
#include "for_debug.inc"
#include "theo.inc"

      integer istage,systypeMod

C To be used as a seed:
      integer isys

      integer n0
      double precision s,voica,dummy_st,sorig

C Single precision here:
      real rndsh(3)  ! additive, poisson, linear
     $     ,ranflat

      double precision rand_shift(NSYS)
      double precision r_sh_fl(NSYS)
      double precision f_un
      parameter (f_un = 2.0)   ! translate 0.:-1 to -1.:1.

      real amu
      integer npoi, ierr
C For log normal random shifts:
      real lsig, lmu,lrunif

      double precision epsilon    ! estimated acceptance/lumi correction
      double precision data_in
      double precision estat_in, ecor_in, euncor_in, etot_in !> Input uncertainites
      double precision estat_out,euncor_out  ! recalculated stat. error

      double precision scaleF
      integer scaling_type

C functions:
      real logshift
      double precision alnorm
C------------------------------------------------------------



cv initialise the random shifts
      do isys=1,nsys
         rand_shift(isys)=0.
         r_sh_fl(isys)=0.
      enddo


C
C Loop over systematic sources:
C
      open(90,file=TRIM(OutDirName)//'/mc_rnd.txt')
      do isys=1,nsys
         call rnorml(rndsh,1)   ! gauss random number
         call ranlux(ranflat,1) ! uniform random number

         rand_shift(isys) = rndsh(1)
         r_sh_fl(isys) = (ranflat-0.5)*f_un

         print '(''random numbers: sys, gauss, flat '',2i6,2F8.2)',
     $        isys,isys, rand_shift(isys),
     $        r_sh_fl(isys)
         write(90,'(''random numbers: sys, gauss, flat '',2i6,2F8.2)')
     $        isys,isys, rand_shift(isys),
     $        r_sh_fl(isys)
      enddo
      close(90)

C
C Loop over the data:
C
      do n0=1,npoints
         call rnorml(rndsh,3)
         call ranlux(ranflat,1)

         if (lrandData) then
            s = DATEN(n0)
         else
            s = THEO(n0)
         endif
         sorig = s

         do isys=1,nsys

cv  test different distributions
cv  first for systematic uncert, then for stat.

            if (systype.eq.1) then ! gauss syst
C ! Introduce asymmetric errors, for Gaussian case only:
               if ( .not. LAsymSyst(isys) ) then
                  s = s*(1.+ beta(isys,n0) * rand_shift(isys))
               else
                  if ( rand_shift(isys).gt. 0) then
                     s = s*(1.+ BetaAsym(isys,1,n0) * rand_shift(isys))
                  else
                     s = s*(1.+ BetaAsym(isys,2,n0) * rand_shift(isys))
                  endif
               endif

            elseif (systype.eq.2) then ! uniform
               s = s*(1. + beta(isys,n0) * r_sh_fl(isys))

            elseif (systype.eq.3) then ! lognormal
               if (beta(isys,n0).ne.0) then
                  lsig=beta(isys,n0)
                  lmu=1.
                  lrunif=r_sh_fl(isys)/f_un + 0.5  ! Expect random number between 0 and 1.
                  s=s*logshift(lmu,lsig,lrunif)
c                  print*,'log...', n0,isys,
c     $                 lrunif, beta(isys,n0),
c     $                 s,logshift(lmu,lsig,lrunif)
               endif
            endif               ! endif (sys for systematic shifts)
         enddo                  ! end loop over the systematic shifts

         voica=s                ! save cross section before the stat shift

CV now choose sta (advised gauss OR poisson)

         if (statype.eq.1) then ! gauss
C do separate fluctuations for stat-const, stat-poisson and stat-linear pieces
            s = s
     $         + sqrt( e_uncor_const(n0)**2 + e_stat_const(n0)**2)
     $              * daten(n0)*rndsh(1)
     $         + sqrt( e_uncor_poisson(n0)**2 + e_stat_poisson(n0)**2)
     $              * sqrt(abs(daten(n0)*sorig))*rndsh(2)
     $         + e_uncor_mult(n0)*sorig*rndsh(3)
         elseif (statype.eq.3.) then ! lognormal
            lsig = alpha(n0)
            lmu=1.
            call ranlux(lrunif,1)
            s=s*logshift(lmu,lsig,lrunif)
         elseif (statype.eq.2) then ! uniform
            s = s + alpha(n0)*f_un*(ranflat-0.5)
         elseif (statype.eq.4 .or. statype.eq.14
     $           .or. statype.eq.34) then         ! poisson & mixed
C Poisson require theory:
            if (lranddata) then
               print *,'Poisson errors require theory predictions!'
               call HF_stop
            endif

            data_in    = daten(n0)
            estat_in   = alpha(n0)
            euncor_in  = e_unc(n0)/100.D0*data_in
            etot_in    = e_tot(n0)/100.D0*data_in
            ecor_in    = sqrt(etot_in**2-estat_in**2)

            if (statype.eq.14 .or. statype.eq.34) then
C Calculate stat:
               if (estat_in.gt.euncor_in) then
                  estat_in = sqrt(estat_in**2-euncor_in**2)
               else
                  call HF_errlog(12020520,
     +            'W: mc_errors: uncor. error larger than stat')
                  estat_in = 0.
               endif
            else
C Reset uncor:
               euncor_in = 0.
            endif

C Get acceptance/lumi correction, called "epsilon"
            epsilon = data_in/estat_in**2

C Expected number of events:
            amu = epsilon*theo(n0)
            call RNPSSN(amu, Npoi, Ierr)

            s = (s/THEO(n0)) * Npoi/epsilon

C Also apply fluctuations due to uncorrelated systematics:

C New absolute uncor:
            euncor_out = euncor_in / data_in * s ! rescale to new value

            if (statype.eq.14) then
               s = s + rndsh(1)*euncor_in
            else if (statype.eq.34) then
               lsig = euncor_in
               lmu=1.
               call ranlux(lrunif,1)
               s=s*logshift(lmu,lsig,lrunif)
            endif

C  Redefine stat errors:

            if (Npoi.le.0) then
C Set small value and large uncertainty, effectively excluding the measurement:
               s = 1.0D-4
               estat_out = 10.0D0
            else
               estat_out   = sqrt(Dble(Npoi))/epsilon
            endif


C Store uncor in %:
            e_unc(n0) = euncor_out/s*100.0


            alpha(n0) = sqrt(euncor_out**2+estat_out**2)
            e_tot(n0) = sqrt(euncor_out**2+estat_out**2+ecor_in**2)
     $           /s*100.0
         endif


         print
     $ '(''Original, systematics and stat. shifted data:'',i4,5E12.4)'
     $        , n0,sorig, voica,s,alpha(n0),e_unc(n0)/100.*s


         if (s.eq.0) then
            call hf_errlog(1302201902,
     $          'S: ToyMC cross section with exact ZERO value, stopOB')
         endif

C     Scale relative error sources, depending on scaling rule defined in chi2-related section of steering or in data files
C         For :
C     - additive        ("NoRescale") keep absolute errors unmodified
C     - multiplicaitive ("Linear")    keep relative errors unmodified
C     - poisson         ("Poisson")   keep (relative error)*sqrt(value) unmodified

         scaleF = DATEN(n0)/s !=oldValue/newValue

         if (s .lt. 0) then
            call hf_errlog(1302201901,
     $          'W: ToyMC cross section with negative value')
            e_uncor_const(n0) = sqrt(e_uncor_const(n0)**2
     $           +e_uncor_poisson(n0)**2 )
            e_stat_const(n0)  = sqrt(e_stat_const(n0)**2+
     $           e_stat_poisson(n0)**2)
            e_stat_poisson(n0)  = 0.
            e_uncor_poisson(n0) = 0.
         else
            e_stat_poisson(n0) = e_stat_poisson(n0)  * sqrt(scaleF)
            e_uncor_poisson(n0) = e_uncor_poisson(n0)  * sqrt(scaleF)
         endif
         e_uncor_mult(n0) = e_uncor_mult(n0) !  NO scale, keep relative
         e_uncor_const(n0) = e_uncor_const(n0) * scaleF
         e_stat_const(n0)  = e_stat_const(n0) * scaleF
         e_tot(n0) = e_tot(n0) * scaleF

C     Also correlated systematicss:
         do isys=1,nsys

            scaling_type = SysScalingType(isys)

            if (
     $           (scaling_type .eq. isNoRescale)
     $           .or. (LForceAdditiveData(n0) )
     $           ) then         ! additive, keep absolute
               beta(isys,n0) = beta(isys,n0) * scaleF
               omega(isys,n0) = omega(isys,n0) * scaleF
            elseif (scaling_type.eq. isLinear) then  ! mult, do nothing
               beta(isys,n0) = beta(isys,n0)
               omega(isys,n0) = omega(isys,n0)
            elseif (scaling_type.eq. isPoisson) then
               beta(isys,n0) = beta(isys,n0) * sqrt(scaleF)
               omega(isys,n0) = omega(isys,n0) * sqrt(scaleF)
            endif
         enddo

         DATEN(n0) = s
C update alpha:
         alpha(n0) =  sqrt(e_uncor_mult(n0)**2
     $        +e_stat_poisson(n0)**2
     $        +e_uncor_const(n0)**2
     $        +e_stat_const(n0)**2
     $        +e_uncor_poisson(n0)**2)
     $        *daten(n0)

      enddo
      end
C---------------------------------------------------
C Created by SG, 23 Apr 2008 following
C
C http://www.brighton-webs.co.uk/distributions/lognormal.asp
C
!> Log-normal random number generator
!> @param[in] am mean value
!> @param[in] as RMS
C
C Input:  am -- mean value
C         as -- RMS
C----------------------------------------------------
      function alnorm(am,as)
*     ---------------------------------------------
      implicit none
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
*     -------------------------------------------
        real function logshift(mmu,ssig,rrunif)
*     -------------------------------------------

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

      call rzero(0.0001,10000.,x2,ex2,0.0001,10000,zeroth)
      logshift=x2

C-----------------------------------------------------------------------
 999  CONTINUE
      END


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
!>
!> @param XX
!> @param i
C----------------------------------------------------------------------
      REAL FUNCTION ZEROTH(XX,i)

      IMPLICIT NONE

      real runif,erf
      integer i, iseed
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


C------------------------------------------------------------------------
C
!>  Separate random number seed initialisation.
C
C------------------------------------------------------------------------
      subroutine init_rnd_seeds()

      implicit none
#include "steering.inc"
      integer icount,ntime,ndate
C Common from CERNLIB datime:
      integer IS
      common/SLATE/ IS(40)

C-------------------------------------------
      if (iseedmc.ne.0) then
C Seed from the steering:
         icount = iseedmc
      else
C Seed from current time
C Does not actually work, since datime is a dummy placeholder --Ivan
         call datime(ndate,ntime)
         ntime = ntime*1000000+is(6)
         icount=ntime
      endif

cv initialise the random seed gener

      call rmarin(icount,0,0)
      call rluxgo(3,icount,0,0)
      print*,'initialize smearing with a seed isdrn = ',icount


C--------------------------------------------------------------


      end
