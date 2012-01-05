#include <iostream>

#include "Hathor.h"
#include "../interface/H1FitterPdf.h"

extern "C" {
  int hathorinit_(const double& sqrtS, const bool& ppbar);
  int hathorcalc_(const int *idataset, double *xsec);
}

Hathor* hathor;
H1FitterPdf* pdf;
// FIXME: delete pointers at the end! (in some hathordestroy_ or so)

int hathorinit_(const double& sqrtS, const bool& ppbar) {
  pdf = new H1FitterPdf();
  hathor = new Hathor(*pdf);

  if(ppbar)
    hathor->setColliderType(Hathor::PPBAR);
  else
    hathor->setColliderType(Hathor::PP);

  hathor->setSqrtShad(sqrtS);

  hathor->setScheme(Hathor::LO | Hathor::NLO | Hathor::NNLO // make me configurable from steering file!
		    | Hathor::POLE_MASS);                   // make me configurable from steering file!
  hathor->setPrecision(Hathor::MEDIUM);                     // make me configurable from steering file!

  return 0;
}

int hathorcalc_(const int *idataset, double *xsec) {
  const double mt = 173.;  // make me configurable!
  hathor->getXsection(mt,mt,mt);

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
