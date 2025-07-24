


      program lhapdf_example
      use, intrinsic :: iso_c_binding
      implicit none

! LHAPDF interface
      
! Parameters
      character(len=40) :: setname
      real(8) :: x, Q
      integer :: pid
      integer iset

      double precision xf(13)

      integer gid

      double precision data(200)
      double precision datas(200)
      double precision dataserr(200)
      
      double precision ref(200)
      double precision referr(200)

      integer getnbins

      integer nbins

      integer i
      integer ig
      
! LHAPDF init
!     setname = "CT10nlo"
      setname = "NNPDF31_nnlo_as_0118"
!     setname = "NNPDF40_nnlo_as_01180"

      iset = 0
      
      call InitPDFset(setname)
      call InitPDF(iset)
      
      gid=0

      call readgrid( gid, "/Users/sutt/nlo2/pineappl-grids/Z_7TeV_NNLO/cc_etay_bin_0.appl" )
      call readgrid( gid, "/Users/sutt/nlo2/pineappl-grids/Z_7TeV_NNLO/cc_etay_bin_1.appl" )
      call readgrid( gid, "/Users/sutt/nlo2/pineappl-grids/Z_7TeV_NNLO/cc_etay_bin_2.appl" )
      call readgrid( gid, "/Users/sutt/nlo2/pineappl-grids/Z_7TeV_NNLO/cf_etay_bin_0.appl" )
      call readgrid( gid, "/Users/sutt/nlo2/pineappl-grids/Z_7TeV_NNLO/cf_etay_bin_1.appl" )
      
!      call readgrid( gid, "/Users/sutt/nlo2/pineappl-grids/W_Z_7TeV_NNLO/cc_etay_bin_0.appl" )
!      call readgrid( gid, "/Users/sutt/nlo2/pineappl-grids/W_Z_7TeV_NNLO/cc_etay_bin_1.appl" )
!      call readgrid( gid, "/Users/sutt/nlo2/pineappl-grids/W_Z_7TeV_NNLO/cc_etay_bin_2.appl" )
!      call readgrid( gid, "/Users/sutt/nlo2/pineappl-grids/W_Z_7TeV_NNLO/cf_etay_bin_0.appl" )
!      call readgrid( gid, "/Users/sutt/nlo2/pineappl-grids/W_Z_7TeV_NNLO/cf_etay_bin_1.appl" )

!      call readgrid( gid, "/Users/sutt/nlo2/pineappl-grids/W_Z_7TeV_NNLO/Wm_etay.appl" )
!      call readgrid( gid, "/Users/sutt/nlo2/pineappl-grids/W_Z_7TeV_NNLO/Wp_etay.appl" )

      ! these two are full LO, NLO, NNLO reference grids with stored smoothing matrix
      call readgrid( gid, "/Users/sutt/nlo2/incjets/total-fixed-refs/NNLO.y8_ptj.appl-smooth-nlo.appl")
      call readgrid( gid, "/Users/sutt/nlo2/incjets/total-fixed-refs/NNLO.y0_ptj.appl-smooth-nlo.appl")

      
      do ig=0,gid
      
         nbins = getnbins(ig)

         print *, "nbins: ", nbins

         call getreference( ig, ref, referr )
                  
         call convolute( ig, data )
         
         call sconvolute( ig, datas, dataserr )

         
         do i=1,nbins
            print *, i, datas(i)/ref(i),  " +- ", dataserr(i)/ref(i), &
                 " ::  ", data(i)/ref(i), " +- ", referr(i)/ref(i), &
                 " ::  ", ref(i), " +- ", referr(i)
                 
                  
         end do

      end do

         
      end program lhapdf_example



      double precision function fnalphas(Q)
        implicit none 
        double precision, intent(in) :: Q
        double precision alphaspdf       
        
        fnalphas = alphaspdf(Q)
      end function fnalphas


      subroutine fnpdf(x,Q,xf)
        implicit none 
        double precision, intent(in) :: x
        double precision, intent(in) :: Q
        double precision xf(13)

        call evolvepdf(x,Q,xf)
      end 

