C
C Steering for NNPDF reweighting (switched on or off)
C

	integer NNPDFREWEIGHTMETHOD
	integer NNPDFOUTREPLICAS
	character*128  NNPDFSET
	character*128  NNPDFRWDATA

      common/nnpdf/            !> Should be initialized via steering
     $     NNPDFREWEIGHTMETHOD,
     $     NNPDFSET,NNPDFRWDATA,NNPDFOUTREPLICAS

