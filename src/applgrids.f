c----------------------------------------------------------
c     read ngrids applgrids from "applgrids_list.txt"
c----------------------------------------------------------
      subroutine getAPPLgrids(ng)
      integer ng, n
      integer fid
      integer dummy
      parameter (fid = 44)
      character(256) grid_file

      open(unit=fid,FILE='applgrids_list.txt',STATUS='unknown')
      do ig=1,ng
        read(fid, '(A)') grid_file
	call appl_readgrid( dummy, grid_file//char(0) )
      enddo

      call appl_ngrids(n)
      print *, n, ' applgrid grids have been read'
      close(fid)

      return 
      end

c----------------------------------------------------------
c     routines that call  the pdf and alphas routines
c     from QCDNUM
c----------------------------------------------------------
      double precision function appl_fnalphas(Q)
      double precision Q, Q2      ! fact scale
      double precision mur2       ! renorm scale
      double precision alphaspdf
      integer nfl, ialerr
      Q2 = Q**2
      mur2 = RFROMF(Q2) ! get renorm scale from factscale in QCDNUM

      ialerr = 0
c     get alpha_s interpolation
      appl_fnalphas = ASFUNC(mur2, nfl, ialerr)
      if ( ialerr .eq. 1) then
        print *, 'alpha_s convolution error. Renorm scale too low:'
	print *, 'rens^2 = ',mur2, ' < 0.1GeV^2'
	print *, 'setting alpha_s at rens^2 = 0.1GeV^2'
	mur2 = 0.1
	appl_fnalphas = ASFUNC(mur2, nfl, ialerr)
      endif

      return
      end

      subroutine appl_fnpdf(x, Q, xf)
      double precision x, Q, Q2
      double precision xf(-6:6)
      integer iqnset, iqnchk

      iqnset = 1
      iqnchk = 0
      do ifl=-6,6
        xf(ifl)=0.d0
      enddo
      if ( x .lt. 1.d-7 .or. x .gt. 1d0-1d-7 ) return

      Q2 = Q*Q

      call fpdfxq(iqnset, x, Q2, xf, iqnchk)
      return
      end
c----------------------------------------------------------

      subroutine ag_ngrids( ng )
c      integer ng
      call appl_ngrids(ng)

      return 
      end


      subroutine ag_gridids(igrids)

c      integer igrids(100)
      call appl_gridids(igrids)

      return 
      end

      integer function ag_getnbins(igrid)

      integer appl_getnbins
      ag_getnbins = appl_getnbins(igrid)
      return 
      end

      subroutine ag_convolute(igrid, xsec)

      double precision xsec(100)
c      integer igrid

      call appl_convolute(igrid,xsec)
      return 
      end

      subroutine ag_releasegrids
      call appl_releasegrids
      return 
      end
