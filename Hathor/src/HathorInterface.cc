#include <iostream>

#include "Hathor.h"
#include "../interface/HERAFitterPdf.h"

extern "C" {
  int hathorinit_(const double& sqrtS, const bool& ppbar, const double& mt,
		  const unsigned int& pertubOrder, const unsigned int& precisionLevel);
  int hathorcalc_(const int *idataset, double *xsec);
}

extern "C" {
  int rlxd_size(void);
  void rlxd_get(int state[]);
  void rlxd_reset(int state[]);
  void rlxd_init(int level,int seed);
}

Hathor* hathor;   // FIXME: delete pointers at the end! (in some hathordestroy_ or so)
HERAFitterPdf* pdf; // FIXME: delete pointers at the end! (in some hathordestroy_ or so)
double mtop;

int *RndStore; // FIXME: delete at the end

int hathorinit_(const double& sqrtS, const bool& ppbar, const double& mt,
		const unsigned int& pertubOrder, const unsigned int& precisionLevel) {
  pdf = new HERAFitterPdf();
  hathor = new Hathor(*pdf);

  if(ppbar)
    hathor->setColliderType(Hathor::PPBAR);
  else
    hathor->setColliderType(Hathor::PP);

  hathor->setSqrtShad(sqrtS);

  unsigned int scheme = Hathor::LO;
  if(pertubOrder>1)
    scheme = scheme | Hathor::NLO;
  if(pertubOrder>2)
    scheme = scheme | Hathor::NNLO;
  hathor->setScheme(scheme);

  hathor->PrintOptions();

  hathor->setPrecision(pow(10,2+precisionLevel));

  mtop = mt;
  std::cout << " Top mass and renorm./fact. scale used for Hathor [GeV]: " << mtop << std::endl;

  // Random number setup:
  rlxd_init(1,1);
  int nRnd = rlxd_size();

  std::cout << " Size of random number array = " << nRnd << "\n";
  RndStore = new int [nRnd];
  rlxd_get(RndStore);

  return 0;
}

int hathorcalc_(const int *idataset, double *xsec) {

  // Reset random numbers
  rlxd_reset(RndStore);

  hathor->getXsection(mtop, mtop, mtop);

  double val,err;
  hathor->getResult(0,val,err);

  xsec[0] = val;
  
  return 0;
}
