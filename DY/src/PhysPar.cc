/*!
  @file PhysPar.cc
  @author Andrey Sapronov <Andrey.Sapronov@cern.ch>
  @date Sun May 17 2009
 
  Copyright (c) 2008-2009 Andrey Sapronov
 
  The values are assigned here.
 
  Adjusted to comply with MCFM
 
 */
 /*****************************************************************************/

#include <cmath>

namespace PhysPar {
  int ih1		= 1;
  int ih2		= 1;
  double ebeam		= 7000.;
  // scales
  double omega		= 1e-5;

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
  //double V2[4]		= { pow(0.9742,2), pow(0.9733,2), pow(0.226,2), pow(0.226,2) };
  double V2[4]		= { pow(0.97419,2), pow(0.97334,2), pow(0.22570,2), pow(0.22560,2) }; // MCFM
  				
  // current run parameters
  int idlept		= 11;
  double ml		= mel;
  double m2l		= pow(mel,2);
  double mq		= mup;
  double m2q		= pow(mup,2);
};

