c----------------------------------------------------------
c     this file contains dummy applgrid interface routines
c     which are called in case xFitter was not compiled
c     with --enable-applgrid option
c----------------------------------------------------------
      subroutine getAPPLgrids(ng)
      integer ng

      call print_ag_messg
      return 
      end

      subroutine update_theor_ckm
c      call print_ag_messg
      return
      end


      subroutine appl_readfastnlogrids
      end


      subroutine ag_ngrids( ng )
      integer ng

      call print_ag_messg
      return 
      end


      Subroutine Calc_pdf_applgrid_fast
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

      subroutine ag_convolute(igrid, iorder, muR, muF, xsec)

      double precision muR,muF,xsec
      integer igrid,iorder

      call print_ag_messg
      return 
      end

      subroutine ag_releasegrids
      call print_ag_messg
      return 
      end

      subroutine print_ag_messg
      print *, '--------------------------------------------------'
      print *, 'APPLGRID: You have chosen to use applgrid but xFitter is not'
      print *, 'compiled with --enable-applgrid option.'
      return 
      end


      subroutine appl_readgrid
      print *,'CALL APPLGRID with DUMMY APPLGRID interface!!'
      call HF_stop
      end

      subroutine appl_getbinnumber
      print *,'CALL APPLGRID with DUMMY APPLGRID interface!!'
      call HF_stop
      end

      subroutine appl_getbinlowedge
      print *,'CALL APPLGRID with DUMMY APPLGRID interface!!'
      call HF_stop
      end

      subroutine appl_getbinwidth
      print *,'CALL APPLGRID with DUMMY APPLGRID interface!!'
      call HF_stop
      end

      subroutine appl_ngrids(n)
      integer n
C
      n = 0
      print *,'NO APPLGRIDs INITIALIZED'
      print *,'NO APPLGRID interface comiled'
      print *,'DUMMY version applgrids_dum has been called!'
      call HF_errlog(12020601,'W: DUMMY version applgrids_dum called!')
C
      end



