c----------------------------------------------------------
c     this file contains dummy lhapdf interface routines
c     which are called in case h1fitter was not compiled
c     with --enable-lhapdf option
c----------------------------------------------------------

      subroutine print_lhapdf_messg
      print *,'calling lhapdf function but h1fitter compiled without --enable-lhapdf switch'
      return
      end
      
      subroutine initpdfsetbyname(a)
      call print_lhapdf_messg
      stop
      return 
      end

      subroutine initpdf(a)
      call print_lhapdf_messg
      stop
      return 
      end

      subroutine initpdfset(a)
      call print_lhapdf_messg
      stop
      return 
      end

      subroutine evolvepdf(a)
      call print_lhapdf_messg
      stop
      return 
      end


