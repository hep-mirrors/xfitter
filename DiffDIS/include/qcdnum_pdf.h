#ifndef QCDNUM_PDF_H_
#define QCDNUM_PDF_H_

#include "pdf_base.h"

extern "C" {
  // --- Fortran routines
  void qcdnumgetall_(double *x, double *QQ, double *f);
  double qcdnumget_(double *x, double *QQ, int *iparton);
}

class qcdnum_pdf_t : public pdf_base_t {
  // ================================================================
  void Vals(double x, double QQ, double f[], int xpow_out=0) {
    memset(f, 0, (MAX_FLAVORS+1)*sizeof(f[0]));
    if(fabs(1-x) < 1e-8) return;
    qcdnumgetall_(&x, &QQ, f);
    for(int j=0; j <= 5; j++) f[j] *= pow(x, xpow_out-1);
  }
  
  double Val(double x, double QQ, int pn, int xpow_out=0) {
    if(fabs(1-x) < 1e-8) return 0;
    return qcdnumget_(&x, &QQ, &pn) * pow(x, xpow_out-1);
  }
};

#endif
