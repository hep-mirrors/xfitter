C
C Steering for PDF reweighting (switched on or off)
C

	integer RWMETHOD
	integer RWOUTREPLICAS
	character*128  RWPDFSET
	character*128  RWDATA

      common/reweighting/            !> Should be initialized via steering
     $     RWMETHOD,
     $     RWPDFSET,RWDATA,RWREPLICAS,RWOUTREPLICAS

