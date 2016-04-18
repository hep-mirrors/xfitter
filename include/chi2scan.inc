#include "dimensions.h"
C
C Steering for chi2 scan
C

      integer Nchi2points
      parameter (Nchi2points = NCHI2POINTS_C)

c label for the parameter
      character*64 label
c central value (to be used for reference theory)
      double precision central
c List of parameter points
      double precision values(Nchi2points)

c In case the parameter scan is done with applgrids or tables, need to specify a list of grids or tables per dataset
c Generalise to n datasets
      integer dataid(NSET)
c Generalise to m terms
      character*8 term(NTermsMax,NSET)
c Generalise to n datasets and m terms
      character*1000 TheorySources(Nchi2points,NTermsMax,NSET)

      logical SCAN
      logical pdferrors
      logical pdfprofile
      logical scaleprofile

      character*128  chi2lhapdfref        !> lhapdf reference set (to be used for reference theory)
      character*128  chi2lhapdfset        !> lhapdf test set
      character*128  chi2lhapdfvarset     !> lhapdf additional variation set
      integer	     chi2nparvar		  !> For LHAPDF VAR PDF, number of parameter variations

      integer	     chi2parpoint	         !> Set this parameter point and do a chi2 fit

      common/chi2scan/           
     $     label,central,values,
     $     dataid,term,TheorySources,SCAN,pdferrors,
     $	   pdfprofile,scaleprofile,
     $	   chi2lhapdfref,chi2lhapdfset,chi2lhapdfvarset,chi2nparvar,
     $	   chi2parpoint
