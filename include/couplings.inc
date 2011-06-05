C alpha electromagnetic, pi and conversion factor:
      double precision alphaem,pi,convfac
      parameter (alphaem = 1./137.0D0)
      parameter (pi = 3.141592653589793D0)
      parameter (convfac = 3.893796D8)

C Fermi constant:
      Double Precision Gf
      parameter (Gf = 1.16637d-5)

C Electric charges:
      double precision euq,edq, e2u,e2d 
      parameter (euq =  2./3.)
      parameter (edq = -1./3.)
      parameter (e2u = euq*euq)
      parameter (e2d = edq*edq)      

      double precision sin2thw,cos2thw,Mz,Mw,cau,cvu,cad,cvd
      common/couplings/sin2thw,cos2thw,Mz,Mw,
     +     cau,cvu,cad,cvd
