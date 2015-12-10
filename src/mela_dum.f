c----------------------------------------------------------
c     this file contains dummy MELA interface routines
c     which are called in case xFitter was not compiled
c     with --enable-mela option
c----------------------------------------------------------
      subroutine print_mela_error_message
      print *, '--------------------------------------------------'
      print *, 'MELA: You have chosen to use MELA but xFitter'
      print *, 'is not compiled with --enable-mela option.'
      call exit(-10)
      return 
      end

      subroutine mela_init
      call print_mela_error_message
      return 
      end

      subroutine readparameters
      call print_mela_error_message
      return 
      end

      subroutine SetxFitterParametersMELA
      call print_mela_error_message
      return 
      end

      subroutine xstructurefunctions
      call print_mela_error_message
      return 
      end
