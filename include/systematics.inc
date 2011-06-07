      integer NSYS

      parameter (NSYS = 200)

      character*16 SYSTEM
      dimension SYSTEM(NSYS)
      double precision BETA(NSYS,NTOT)
      double precision ALPHA(NTOT)
      double precision sysa(nsys,nsys)
      integer npoints,ndis
      common/systema/ALPHA,BETA,sysa,npoints,ndis,SYSTEM

      
