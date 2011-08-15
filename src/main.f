      program H1Fitter 
C--------------------------------------------------------
C
C> H1 PDF Fit Program
C
C-------------------------------------------------------

      implicit none
      external fcn

      include 'steering.inc'
      include 'thresholds.inc'
      include 'couplings.inc'
      include 'for_debug.inc'

      integer icond
C-----------------------------------------------------
*     ------------------------------------------------
*     Print info message
*     ------------------------------------------------
      call h1fitterInfo

*     ------------------------------------------------
*     Read the steering file steering.txt
*     ------------------------------------------------
      call read_steer

*     ------------------------------------------------
*     Read the measured data points
*     ------------------------------------------------
      call read_data

*     ------------------------------------------------
*     Initialise theory modules
*     ------------------------------------------------
      call init_theory_modules

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

C-----------------------------------------------------



      Subroutine H1fitterinfo
*     ------------------------------------------------
      print *,' '
      print *,' '
      print *,'----------------------------------------------------------------------'
      print *,'                                                                      '
      print *,'  HH   HH    11     FFFFFFF  II  TTTTTTTT  TTTTTTTT  EEEEEEE  RRRRR   '
      print *,'  HH   HH   111     FF       II     TT        TT     EE       RR   RR '
      print *,'  HHHHHHH 11 11     FFFFF    II     TT        TT     EEEEEE   RR   RR '
      print *,'  HH   HH    11     FF       II     TT        TT     EE       RRRRR   '
      print *,'  HH   HH    11     FF       II     TT        TT     EE       RR  RR  '
      print *,'  HH   HH    11     FF       II     TT        TT     EEEEEEE  RR   RR '
      print *,'                                                                      '
      print *,'  Version 0.1.0                                                     '
      print *,'  https://svnsrv.desy.de/desy/h1fitter            http://h1.desy.de   '
      print *,'----------------------------------------------------------------------'
      print *,' '
      print *,' '
c      stop
      end

*     ------------------------------------------------
