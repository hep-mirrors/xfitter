C
C Buffer for fast PDF calls for APPLGRID
C
      integer NAPPLPDFMAX  !> max. number of entries
      parameter (NAPPLPDFMAX=100000)   ! 17 Feb 2017: reduce default size (rarely used)
      integer NAPPLPDF     !> actual number of entries
      integer NAPPLPDFINT  !> number of entries in internal storage

!> @brief  reference for applgrid buffer:
      integer iRefApp(2*NAPPLPDFMAX)

!> @brief x,Q2 values for applgrid buffer:
      double precision XAPPLPDF(NAPPLPDFMAX), QAPPLPDF(NAPPLPDFMAX)      

C PDF store:
      double precision APPLPDF(NAPPLPDFMAX,-6:6)
      common/cAPPLPDF/APPLPDF,XAPPLPDF,QAPPLPDF,NAPPLPDF,IRefApp
     $     ,NAPPLPDFINT
