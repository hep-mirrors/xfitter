      program FitPDF

      implicit none
      external fcn

      include 'steering.inc'
      include 'thresholds.inc'
      include 'couplings.inc'
      include 'lprint.inc'

      integer icond

*     ------------------------------------------------
*     Read the steering file steering.txt
*     ------------------------------------------------
      call read_steer

*     ------------------------------------------------
*     Read the measured data points
*     ------------------------------------------------
      call read_data

*
* Initialise EW parameters
*
      call Init_EW_parameters


*     ------------------------------------------------
*     Initialise theory module
*     ------------------------------------------------
C-1- 29/07/2010 -----------------------------------
C---  added '.or.itheory.eq.3' into the condition
C---  to enable the Dipole model
      if(itheory.eq.0.or.itheory.eq.3) then
C-2- 29/07/2010: end of the addition --------------
         call qcdnum_ini
      elseif(itheory.eq.1) then       
c          here goes a call to a CASCADE ini subroutine, if needed         
      endif

*
* Initialise calculations for each dataset:
*
      call Init_theory_datasets

*     ------------------------------------------------
*     Do the fit
*     ------------------------------------------------
      open(85,file='output/Results.txt')
      call minuit_ini
      lprint = .true.
      call minuit(fcn,0)
      close(85)
      if (DOBANDS) then
         write(6,*) ' --- Calculate error bands ...'
         lprint = .false.    
         call MNCOMD(fcn,'ITERATE 10',icond,0)
         call MNCOMD(fcn,'MYSTUFF 1000',icond,0)
         call MNCOMD(fcn,'MYSTUFF 2000',icond,0)
         call Error_Bands_Pumplin
      endif
      close (24)
      close (25)

*     ------------------------------------------------
*     Done
*     ------------------------------------------------
      stop
      end

