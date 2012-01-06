#include <iostream>

#include "Hathor.h"
#include "../interface/H1FitterPdf.h"

extern "C" {
  int hathorinit_(const double& sqrtS, const bool& ppbar, const double& mt,
		  const unsigned int& pertubOrder, const unsigned int& precisionLevel);
  int hathorcalc_(const int *idataset, double *xsec);
}

Hathor* hathor;   // FIXME: delete pointers at the end! (in some hathordestroy_ or so)
H1FitterPdf* pdf; // FIXME: delete pointers at the end! (in some hathordestroy_ or so)
double mtop;

int hathorinit_(const double& sqrtS, const bool& ppbar, const double& mt,
		const unsigned int& pertubOrder, const unsigned int& precisionLevel) {
  pdf = new H1FitterPdf();
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

  return 0;
}

int hathorcalc_(const int *idataset, double *xsec) {
  hathor->getXsection(mtop, mtop, mtop);

  double val,err;
  hathor->getResult(0,val,err);

  //  std::cout << val << " +/- " << err << std::endl;

  // rounding precision to be automized!!!
  val *= 10;
  val = floor(val+0.5);
  val /= 10;
  xsec[0] = val;

  //  std::cout << val << std::endl;

  //  xsec[0] = floor(val+0.5);
  
  return 0;
}
