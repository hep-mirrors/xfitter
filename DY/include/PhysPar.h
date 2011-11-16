#ifndef namespace_PhysPar_h
#define namespace_PhysPar_h

/*!
  @file PhysPar.h
  @author Andrey Sapronov <Andrey.Sapronov@cern.ch>
  @date Sun May 17 2009
 
  Namespace for useful physical parameters.  
 */

namespace PhysPar {
  // constants
  extern double conhc; // conversion constant (hc)^2, GeV*pb (?)
  extern double alpha;
  extern double alphai;
  extern double gfermi;
  
  // widths
  extern double wz;
  extern double ww;
  extern double wh;
  extern double wtp;
  
  // boson masses
  extern double mh;
  extern double m2h;
  extern double mw;
  extern double m2w;
  extern double mz;
  extern double m2z;
  
  // fermion masses
  extern double men;
  extern double mel;
  extern double mmn;
  extern double mmo;
  extern double mtn;
  extern double mta;
  extern double mup;
  extern double mdn;
  extern double mch;
  extern double mst;
  extern double mtp;
  extern double mbt;

  // fermion charges
  extern double chq[2];
  extern double chl;

  // isospin
  extern double I3q[2];
  extern double I3l;

  // couplings
  extern double sw2;
  extern double cw2;
  extern double vq[2];
  extern double vl;

  // ckm
  extern double V2[6];
  				
  // current run parameters
  extern int idlept;

  extern void setPhysPar();
};

#endif
