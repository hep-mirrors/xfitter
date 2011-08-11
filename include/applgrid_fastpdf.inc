C
C Buffer for fast PDF calls for APPLGRID
C
      integer NAPPLPDFMAX  !> max. number of entries
      parameter (NAPPLPDFMAX=500000) 
      integer NAPPLPDF     !> actual number of entries


C x,Q2 values:
      double precision XAPPLPDF(NAPPLPDFMAX), QAPPLPDF(NAPPLPDFMAX)      

C PDF store:
      double precision APPLPDF(NAPPLPDFMAX,-6:6)
      common/cAPPLPDF/APPLPDF,XAPPLPDF,QAPPLPDF,NAPPLPDF
