C---------------------------------------------------
C 
!> Read steering file steering.txt
C
C---------------------------------------------------
      subroutine read_steer_py

      call read_steer

      end

      subroutine read_data_py

      call read_data

      end


      subroutine init_theory_py
      call init_theory_modules
      end
     
      subroutine Fit()
      call do_fit
      call close_theor_eval
      call HF_errsum(6)
      end
      
      subroutine Input_Data_Table( dataIn, statIn, systIn,
     $     nData, nSyst)
      implicit none 
      integer nData
      integer nSyst
      real*8 dataIn(nData)
      real*8 statIn(nData)
      real*8 systIn(nSyst,ndata)
Cf2py intent(in) dataIn 
Cf2py intent(in) statIn 
Cf2py intent(in) systIn 
cf2py intent(in) nData 
cf2py intent(in) nSyst

#include "ntot.inc"
#include "steering.inc"
#include "systematics.inc"
#include "indata.inc"
      integer i,j

c-----------------------------------------------------------------------

      Npoints = nData
      do i = 1, nData
         DATEN(i) = dataIn(i)
         E_STAT_CONST(i) = StatIn(i)
         do j=1,nsyst
            beta(j,i) = systIn(j,i)
            n_syst_meas(j) = nData
            syst_meas_idx(i,j) = i
         enddo
      enddo
      
      NSYS = Nsyst

C
      end


      subroutine Input_Theory_Table( theoIn, theoEig,
     $     nData, nTheoEig)
      implicit none 
      integer nData
      integer nTheoEig
      real*8 theoIn(nData)
      real*8 theoEig(nTheoEig,ndata)
Cf2py intent(in) theoIn
Cf2py intent(in) theoEig
cf2py intent(hide) nData  
cf2py intent(hide) nTheoEig

#include "ntot.inc"
#include "steering.inc"
#include "theo.inc"
      integer i,j
      do i=1,ndata
         THEO(i)  = TheoIn(i)
      enddo
      end


      
      subroutine Profile( Chi2 )

      implicit none 

      integer nData
      integer nTheoEig

      real*8 Chi2

cf2py intent(out) Chi2

#include "ntot.inc"
#include "indata.inc"
#include "systematics.inc"
C #include "f2pyblock.inc"
C

      real*8 sys_shifts(2000),sys_shifts_err(2000)
      real*8 dataset_chi2(150), corr_chi2
      
      common/profile_out/sys_shifts, sys_shifts_err, 
     $     dataset_chi2, corr_chi2


      call GetNewChisquare( 3, NPoints, Chi2, sys_shifts, 
     $     sys_shifts_err, dataset_chi2, corr_chi2)
      end


