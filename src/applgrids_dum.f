c----------------------------------------------------------
c     this file contains dummy applgrid interface routines
c     which are called in case h1fitter was not compiled
c     with --enable-applgrid option
c----------------------------------------------------------
      subroutine getAPPLgrids(ng)
      integer ng

      call print_ag_messg
      return 
      end



      subroutine ag_ngrids( ng )
      integer ng

      call print_ag_messg
      return 
      end


      subroutine ag_gridids(igrids)

      integer igrids(100)

      call print_ag_messg
      return 
      end

      subroutine ag_getnbins(igrid)

      integer igrids

      call print_ag_messg
      return 
      end

      subroutine ag_convolute(igrid, xsec)

      double precision xsec
      integer igrid

      call print_ag_messg
      return 
      end

      subroutine ag_releasegrids
      call print_ag_messg
      return 
      end

      subroutine print_ag_messg
      print *, '--------------------------------------------------'
      print *, 'APPLGRID: You have chosen to use applgrid but h1fitter is not'
      print *, 'compiled with --enable-applgrid option.'
      return 
      end

