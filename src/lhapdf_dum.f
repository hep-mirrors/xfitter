c----------------------------------------------------------
c     this file contains dummy lhapdf interface routines
c     which are called in case HERAfitter was not compiled
c     with --enable-lhapdf option
c----------------------------------------------------------

      subroutine print_lhapdf_messg
      print *,'calling lhapdf function but HERAfitter compiled without --enable-lhapdf switch'
      return
      end
      
      subroutine initpdfsetbyname(a)
      call print_lhapdf_messg
      call HF_stop
      return 
      end

      subroutine initpdf(a)
      call print_lhapdf_messg
      call HF_stop
      return 
      end

      subroutine initpdfset(a)
      call print_lhapdf_messg
      call HF_stop
      return 
      end

      subroutine evolvepdf(a)
      call print_lhapdf_messg
      call HF_stop
      return 
      end


      double precision function alphasPDF(a)
      call print_lhapdf_messg
      call HF_stop
      return 
      end

      Subroutine numberpdf()
      call print_lhapdf_messg
      call HF_stop
      return 
      end
