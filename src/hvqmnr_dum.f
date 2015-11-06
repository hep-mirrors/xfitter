c----------------------------------------------------------
c     This routine is called in case the program was
c     not compiled with --enable-hvqmnr option but the
c     corresponding reaction type is chosen in a config file
c----------------------------------------------------------
      subroutine GetHVQMNRXsection(IDataSet)
      implicit none
      integer IDataSet

      print *, '--------------------------------------------------'
      print *, 'For ''HVQMNR pp QQbar'' reaction please configure with --enable-hvqmnr'
      call hf_stop ! Fatal error
      return
      end
