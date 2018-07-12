#include "xfitter_cpp_base.h"
#include "ext_pdfs.h"
#include "iostream"

// Store global pointers here

pXFXlike gAPFEL_XFX = nullptr;

void externalsetapfel1_(const double& x, const double& q, double* pdfs) {
  (*gAPFEL_XFX) (x,q,pdfs);
}

void APFEL_set_pdfs(pXFXlike xfx) {
  gAPFEL_XFX = xfx;
}
