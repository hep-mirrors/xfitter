      integer NSYSMAX

      parameter (NSYSMAX = 200)

      integer nSys

      character*16 SYSTEM
      dimension SYSTEM(NSYSMAX)
      double precision BETA(NSYSMAX,NTOT)
      double precision ALPHA(NTOT)
      double precision sysa(nsysMAX,nsysMAX)
      integer npoints,ndis

      integer CompressIdx(NSYSMAX)
      common/systema/ALPHA,BETA,sysa,npoints,ndis,SYSTEM,NSys,
     $     CompressIdx

      
