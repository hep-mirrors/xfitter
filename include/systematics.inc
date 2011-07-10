      integer NSYSMAX

      parameter (NSYSMAX = 300)

      integer nSys

      character*16 SYSTEM
      dimension SYSTEM(NSYSMAX)
      double precision BETA(NSYSMAX,NTOT)
      double precision ALPHA(NTOT)
      double precision sysa(nsysMAX,nsysMAX)
      integer npoints,ndis

      common/systema/ALPHA,BETA,sysa,npoints,ndis,SYSTEM,NSys

      
