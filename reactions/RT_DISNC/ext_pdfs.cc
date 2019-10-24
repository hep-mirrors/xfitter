/**
  @file ext_pdfs.cc

  @brief Generic external PDFs and alphaS calls for RT codes

  Provides an interface to call PDFs.

  @version 0.1
  @date 2017/04/16
 */

#include "xfitter_cpp_base.h"
#include "ext_pdfs.h"


// Store global pointers here

pXFXlike     gRT_XFX = nullptr;
pOneParFunc  gRT_alphaS = nullptr;

void rt_get_pdfs_(const double& x, const double& q, double* pdfs) {
  (*gRT_XFX) (x,q,pdfs);
}

double rt_get_alphas_(const double& q) {
  return (*gRT_alphaS) (q);
}

void rt_set_pdfs_alphaS( pXFXlike xfx, pOneParFunc aS) {
  gRT_XFX = xfx;
  gRT_alphaS = aS;
}
