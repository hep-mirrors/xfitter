
      subroutine init_theory_modules
*     ------------------------------------------------

      implicit none 
      include 'steering.inc'
      include 'thresholds.inc'
      include 'couplings.inc'


*     ------------------------------------------------
*     Initialise EW parameters
*     ------------------------------------------------

      call Init_EW_parameters

*     ------------------------------------------------
*     Initialise qcdnum
*     ------------------------------------------------

c itheory=3 for dipole -- do we want this here?
      if(itheory.eq.0.or.itheory.eq.3) then
         call qcdnum_ini
      elseif(itheory.eq.1) then       
c          here goes a call to a CASCADE ini subroutine, if needed         
      endif

*     ------------------------------------------------
*     Initialise calculations for each dataset:
*     ------------------------------------------------

      call Init_theory_datasets


      return
      end
