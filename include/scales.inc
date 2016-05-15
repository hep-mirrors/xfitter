C 
C !> Define data-file dependent scales
C	
      double precision DataSetMuR(nset)         !> renormalisation scale
      double precision DataSetMuF(nset)         !> factorisation scale
      double precision DataSetMuRes(nset)       !> resummation scale
      double precision DataSetMuC3(nset)        !> soft resummation scale
      integer          DataSetIOrder(nset)      !> Can change calculation order. 1-LO, 2-NLO, 3-NNLO.

      common/CScales/DataSetMuR,DataSetMuF,DataSetMuRes,
     $ DataSetMuC3,DataSetIOrder
