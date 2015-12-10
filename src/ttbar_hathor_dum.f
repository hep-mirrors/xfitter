c----------------------------------------------------------
c     These routines are called in case the xFitter was
c     not compiled with --enable-hathor option but the
c     ttbar reaction type is chosen in a config file
c----------------------------------------------------------


      Subroutine GetHathorXsection(IDataSet)
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"

      integer IDataSet

      print *, '--------------------------------------------------'
      print *, 'You have chosen to use Hathor but xFitter is not'
      print *, 'compiled with --enable-hathor option.'
      
      call hf_stop ! Fatal error
      
      return

      end


      Subroutine hathorinit()
      print *, '--------------------------------------------------'
      print *, 'You have chosen to use Hathor but xFitter is not'
      print *, 'compiled with --enable-hathor option.'

      call hf_stop ! Fatal error

      return
      end
