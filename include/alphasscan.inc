#include "dimensions.h"
C
C Steering for alphas scan
C

      logical ASSCAN
      logical aspdfprofile
      logical asscaleprofile

      character*128  alphaslhapdf         !> lhapdf alphas scan set
      character*128  aslhapdfset        !> lhapdf set containing all PDF variations
      character*128  aslhapdfvarset     !> lhapdf additional variation set
      integer	     asnparvar		  !> For LHAPDF VAR PDF, number of parameter variations

      character*128  aslhapdfref        !> lhapdf reference set (to be used for reference theory)

      common/alphasscan/           
     $     ASSCAN,
     $	   aspdfprofile,asscaleprofile,
     $	   alphaslhapdf,
     $ 	   aslhapdfset,aslhapdfvarset,asnparvar,
     $	   aslhapdfref
     