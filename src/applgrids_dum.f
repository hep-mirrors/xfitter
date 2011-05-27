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

c----------------------------------------------------------
c     routines that call  the pdf and alphas routines
c     from QCDNUM
c----------------------------------------------------------
      double precision function fnalphas(Q)
      double precision Q, Q2      ! fact scale
      double precision mur2       ! renorm scale
      double precision alphaspdf
      integer nfl, ialerr
      Q2 = Q**2
      mur2 = RFROMF(Q2) ! get renorm scale from factscale in QCDNUM

      ialerr = 0
c     get alpha_s interpolation
      fnalphas = ASFUNC(mur2, nfl, ialerr)
      if ( ialerr .eq. 1) then
        print *, 'alpha_s convolution error. Renorm scale too low:'
	print *, 'rens^2 = ',mur2, ' < 0.1GeV^2'
	print *, 'setting alpha_s at rens^2 = 0.1GeV^2'
	mur2 = 0.1
	fnalphas = ASFUNC(mur2, nfl, ialerr)
      endif

      return
      end

      subroutine fnpdf(x, Q, xf)
      double precision x, Q, Q2
      double precision xf(13)
      integer iqnset, iqnchk

      Q2 = Q**2
      iqnset = 1
      iqnchk = 0
      call fpdfxq(iqnset, x, Q2, xf, iqnchk)
      return
      end
c----------------------------------------------------------


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

