C
C Store useful arrays for QCDNUM calculations and maximum number of points in a grid
C

C QCDNUM:
      double precision cnep2f(-6:6)
      double precision cnem2f(-6:6)
      double precision cnep3f(-6:6)
      double precision cnem3f(-6:6)
      DATA CNEP2F/0.,0.,1.,0.,1.,0.,0.,0.,1.,0.,1.,0.,0./
      DATA CNEM2F/0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0./
      DATA CNEP3F/0.,0.,-1.,0.,-1.,0.,0.,0.,1.,0.,1.,0.,0./
      DATA CNEM3F/0.,-1.,0.,-1.,0.,-1.,0.,1.,0.,1.,0.,1.,0./

      double precision ccep2f(-6:6)
      double precision ccem2f(-6:6)
      double precision ccep3f(-6:6)
      double precision ccem3f(-6:6)
      DATA CCEP2F/0.,0.,1.,0.,1.,0.,0.,1.,0.,1.,0.,0.,0./
      DATA CCEM2F/0.,0.,0.,1.,0.,1.,0.,0.,1.,0.,1.,0.,0./
      DATA CCEP3F/ 0.,0.,-1.,0.,-1.,0.,0.,1.,0.,1.,0.,0.,0./
      DATA CCEM3F/0.,0. ,0.,-1.,0.,-1.,0.,0.,1.,0.,1.,0.,0./

      integer NPmaxDIS
      parameter (NPmaxDIS =  2000)
