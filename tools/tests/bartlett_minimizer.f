C Regression test for non-MINUIT Bartlett preparation.
C After building/installing xFitter and sourcing setup.sh, compile with:
C gfortran -cpp -ffixed-line-length-none -Iinclude -Wl,--export-dynamic
C   tools/tests/bartlett_minimizer.f -Llib -lxfitter -o /tmp/bartlett-test
C Run with the installed lib directory on LD_LIBRARY_PATH, and create
C temp/bartlett-minimizer-output first.
C The mock minimizer exposes two POIs and one free external nuisance;
C another external nuisance is fixed and absent from its parameter array.
      program test_bartlett_minimizer
      implicit none
#include "ntot.inc"
#include "endmini.inc"
#include "systematics.inc"
#include "bartlett_fd.inc"
#include "theo.inc"
#include "indata.inc"
#include "steering.inc"
      integer index_out
      character*80 name_out
      double precision value_out, error_out
      logical finite_unc
      common /mock_uncertainty/ finite_unc
      double precision center(MNE), current(MNE)
      common /mock_parameters/ current
      finite_unc=.false.
      nsys=3
      system(1)='external_free'
      system(2)='external_fixed'
      system(3)='profiled'
      SysForm(1)=isExternal
      SysForm(2)=isExternal
      SysForm(3)=isNuisance
      call Bartlett_CountFreePOI
      if (.not.BartlettHaveNPOI) stop 1
      if (BartlettNPOI.ne.2) stop 2
      if (any(BartlettPOIMinuit(1:2).ne.(/1,3/))) stop 3
      if (SysExtFixed(1).or..not.SysExtFixed(2)) stop 4
      if (SysExtFixed(3)) stop 5
      npoints=2
      center=0d0
      center(1)=2d0
      center(2)=7d0
      center(3)=3d0
      current=center
      BartlettEnabled=.true.
      EoEEnabled=.true.
      call Bartlett_ComputeD(center)
      if (.not.BartlettHaveD) stop 6
      if (maxval(abs(BartlettD(1:2,1)-(/2d0,4d0/)))
     $    .gt.1d-9) stop 7
      if (maxval(abs(BartlettD(1:2,2)-(/3d0,-1d0/)))
     $    .gt.1d-9) stop 8
      if (any(current.ne.center)) stop 9
      if (maxval(abs(theo(1:2)-(/13d0,5d0/))).gt.1d-9) stop 10
      finite_unc=.true.
      OutDirName='temp/bartlett-minimizer-output'
      BartlettLRFactor=0.2d0
      BartlettExtErr=0d0
      BartlettExtErr(1)=0.7d0
      call write_pars(0)
      open(90,file=trim(OutDirName)//'/parsout_0',status='old')
      read(90,*) index_out,name_out,value_out,error_out
      if (index_out.ne.0.or.name_out.ne.'poi_a') stop 11
      if (value_out.ne.2d0) stop 12
      if (abs(error_out-2d0*sqrt(1.1d0)).gt.1d-6) stop 13
      read(90,*) index_out,name_out,value_out,error_out
      if (index_out.ne.1.or.name_out.ne.'external_free') stop 14
      if (value_out.ne.7d0.or.error_out.ne.0.7d0) stop 15
      close(90)
C With EoE disabled, preserve the minimizer's own parameter output.
      open(90,file=trim(OutDirName)//'/parsout_0',status='replace')
      write(90,*) 'raw-minimizer-output'
      close(90)
      EoEEnabled=.false.
      call write_pars(0)
      open(90,file=trim(OutDirName)//'/parsout_0',status='old')
      read(90,*) name_out
      close(90)
      if(name_out.ne.'raw-minimizer-output') stop 16
      print *, 'PASS: POIs, derivatives, restoration and parameter output'
      end

      integer function minimizerusesminuit()
      minimizerusesminuit=0
      end

      integer function getminimizernpars()
      getminimizernpars=3
      end

      subroutine getminimizerparname(i,name)
      integer i
      character*(*) name
      if (i.eq.1) name='poi_a'
      if (i.eq.2) name='external_free'
      if (i.eq.3) name='poi_b'
      end

      double precision function getparamunc(name)
      use, intrinsic :: ieee_arithmetic
      character*(*) name
      logical finite_unc
      common /mock_uncertainty/ finite_unc
C Exercise the fallback step when covariance is unavailable.
      getparamunc=ieee_value(0d0,ieee_quiet_nan)
      if (finite_unc) getparamunc=2d0
      end

      subroutine setfittedparamsfromarray(p)
      implicit none
#include "endmini.inc"
      double precision p(*),current(MNE)
      common /mock_parameters/ current
      current=p(1:MNE)
      end

      subroutine set_scan_parameters(p)
      double precision p(*)
C Deliberately does not update the minimizer's live parameter map.
      end

      subroutine update_theory_iteration
      implicit none
#include "ntot.inc"
#include "endmini.inc"
#include "theo.inc"
      double precision current(MNE)
      common /mock_parameters/ current
      theo(1)=2d0*current(1)+3d0*current(3)
      theo(2)=4d0*current(1)-current(3)
      end

      double precision function getfittedparamd(name)
      implicit none
#include "endmini.inc"
      character*(*) name
      double precision current(MNE)
      common /mock_parameters/ current
      getfittedparamd=0d0
      if(name.eq.'poi_a') getfittedparamd=current(1)
      if(name.eq.'external_free') getfittedparamd=current(2)
      if(name.eq.'poi_b') getfittedparamd=current(3)
      end
