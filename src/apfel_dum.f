c----------------------------------------------------------
c     this file contains dummy APFEL interface routines
c     which are called in case HERAfitter was not compiled
c     with --enable-apfel option
c----------------------------------------------------------
      subroutine apfel_ini
      call print_apfel_error_message
      return 
      end

      subroutine print_apfel_error_message
      print *, '--------------------------------------------------'
      print *, 'APFEL: You have chosen to use APFEL but HERAfitter'
      print *, 'is not compiled with --enable-apfel option.'
      call exit(-10)
      return 
      end

      subroutine alphaqcd
      end


      subroutine setpdfset
      end

      subroutine xpdfall
      end

      subroutine setalphaqcdref
      end


      subroutine evolveapfel
      end

