/**
  @file ext_pdfs.cc

  @brief Generic external PDFs and alphaS calls for ACOT codes

  Provides an interface to call PDFs.

  @version 0.1
  @date 2017/04/16
 */

#include "xfitter_cpp_base.h"
#include "ext_pdfs.h"


// Store global pointers here

pXFXlike     gACOT_XFX = nullptr;
pOneParFunc  gACOT_alphaS = nullptr;

void acot_get_pdfs_(const double& x, const double& q, double* pdfs) {
  (*gACOT_XFX) (x,q,pdfs);
}

double acot_get_alphas_(const double& q) {
  return (*gACOT_alphaS) (q);
}

void acot_set_pdfs_alphaS( pXFXlike xfx, pOneParFunc aS) {
  gACOT_XFX = xfx;
  gACOT_alphaS = aS;
}
