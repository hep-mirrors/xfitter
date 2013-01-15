/*!
  @file PhysPar.cc
  @author Andrey Sapronov <Andrey.Sapronov@cern.ch>
  @date Sun May 17 2009
 
  Function to set DY electroweak parameters from herafitter EW common blocks
 
 */

#include <cmath>

#include "couplings.icc"

namespace PhysPar {
  // constants
  double conhc		= 0.389379323e9; // conversion constant (hc)^2, GeV*pb (?)
  double alphai		= 128.89;
  double alpha		= 1./alphai;
  double alphas		= 0.1176;
  double gfermi		= 1.16637e-5;
  
  // widths
  double wz		= 2.4935;
  double ww		= 2.1054;
  double wh		= 1e-3;
  double wtp		= 1.551;
  
  // boson masses
  double mh		= 115;
  double m2h		= pow(mh,2);
  double mw		= 80.398;
  double m2w		= pow(mw,2);
  double mz		= 91.188;
  double m2z		= pow(mz,2);
  double ma		= 0.;
  double mv		= 91.188;
  
  // fermion masses
  double men		= 1e-10;
  double mel		= 0.51099892e-3;
  double mmn		= 1e-10;
  double mmo		= 0.105658369;
  double mtn		= 1e-10;
  double mta		= 1.77699;
  double mup		= 0.06983;
  double mdn		= 0.06983;
  double mch		= 1.2;
  double mst		= 0.150;
  double mtp		= 174.;
  double mbt		= 4.6;

  // fermion charges
  double chq[2]		= {-1./3.,2./3.}; // {d,u}
  double chl		= -1;

  // isospin
  double I3q[2]		= {-0.5, 0.5};
  double I3l		= -0.5;

  // weinberg angle
  double cw2		= pow(mw/mz,2);
  double sw2		= 1 - cw2;
  // couplings
  double vq[2]		= { -1. + 4./3.*sw2, 1.-8./3.*sw2 };
  double vl		= I3l-chl*sw2;

  // ckm
  double V2[6]		= { pow(0.97419,2), pow(0.97334,2), 
  			    pow(0.22570,2), pow(0.22560,2), 0.,0. }; // MCFM
  				
  // current run parameters
  int idlept		= 11;

  // Function to set parameters according to HERAFITTER electroweak
  // common block
  void setPhysPar();
};

void PhysPar::setPhysPar()
{
  // constants
  conhc		= constants_.convfac;
  gfermi	= constants_.gf;
  
  // widths
  wz		= widths_.wz;
  ww		= widths_.ww;
  wh		= widths_.wh;
  wtp		= widths_.wtp;
  
  // boson masses
  mh		= boson_masses_.mh;
  m2h		= pow(mh,2);
  mw		= boson_masses_.mw;
  m2w		= pow(mw,2);
  mz		= boson_masses_.mz;
  m2z		= pow(mz,2);
  
  // fermion masses
  men		= fermion_masses_.men;
  mel		= fermion_masses_.mel;
  mmn		= fermion_masses_.mmn;
  mmo		= fermion_masses_.mmo;
  mtn		= fermion_masses_.mtn;
  mta		= fermion_masses_.mta;
  mup		= fermion_masses_.mup;
  mdn		= fermion_masses_.mdn;
  mch		= fermion_masses_.mch;
  mst		= fermion_masses_.mst;
  mtp		= fermion_masses_.mtp;
  mbt		= fermion_masses_.mbt;

  // fermion charges
  chq[0]	= -1./3.; // d
  chq[1]	= 2./3.;  // u
  chl		= -1;

  // isospin
  I3q[0]	= -0.5;
  I3q[1]	= 0.5;
  I3l		= -0.5;

  // currently DY uses electroweak G-fermi scheme as in MCFM 
  // weinberg angle
  cw2		= pow(mw/mz,2);
  sw2		= 1 - cw2;
  // alpha_em derived from Gf
  alpha		= 1./M_PI*gfermi*sqrt(2.)*sw2*cw2*m2z;
  alphai	= 1./alpha;

  // couplings derived from sin2thw
  vq[0]		=  -1. + 4./3.*sw2;
  vq[1]		= 1.-8./3.*sw2;
  vl		= I3l-2.*chl*sw2;

  // ckm
  V2[0]		= pow(ckm_matrix_.Vud,2);
  V2[1]		= pow(ckm_matrix_.Vcs,2);
  V2[2]		= pow(ckm_matrix_.Vus,2);
  V2[3]		= pow(ckm_matrix_.Vcd,2);
  V2[4]		= pow(ckm_matrix_.Vub,2);
  V2[5]		= pow(ckm_matrix_.Vcb,2);
  				
  // current run parameters
  idlept		= 11;
};

