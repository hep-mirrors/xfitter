c----------------------------------------------------------
c     this file contains dummy lhapdf interface routines
c     which are called in case xFitter was not compiled
c     with --enable-lhapdf option
c----------------------------------------------------------

      subroutine print_lhapdf_messg
      call hf_errlog(14060201, 'S: Call to lhapdf function but xFitter compiled without --enable-lhapdf switch')
      return
      end
      
      subroutine initpdfsetbyname(a)
      call print_lhapdf_messg
      return 
      end

      subroutine initpdf(a)
      call print_lhapdf_messg
      return 
      end

      subroutine initpdfset(a)
      call print_lhapdf_messg
      return 
      end

      subroutine evolvepdf(a)
      call print_lhapdf_messg
      return 
      end

      subroutine evolvepdfphoton(a)
      call print_lhapdf_messg
      return 
      end


      double precision function alphasPDF(a)
      call print_lhapdf_messg
      return 
      end

      Subroutine numberpdf()
      call print_lhapdf_messg
      return 
      end

      logical function has_photon()
      call print_lhapdf_messg
      has_photon = .false.
      end

