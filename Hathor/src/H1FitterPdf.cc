#include "../interface/H1FitterPdf.h"

#include <iostream>
#include <vector>

extern "C" {
  double asfunc_(const double* r2, int *nf, int *ierr);
  void fpdfxq_(int *iset, const double *x, const double *q2, double *pdfs, int *ichk);
}

H1FitterPdf::H1FitterPdf(const std::string str) {
  PDFname = str;
  imember = 0;
}

void
H1FitterPdf::GetPdf(double x, double muf, double h[13]){
  int iqnset = 1;
  int iqnchk = 0;

  const double muf2 = muf*muf;

  std::vector<double> pdf;
  pdf.resize(13);
  // ordering is:
  // 0   1   2   3   4   5    6   7   8   9   10  11  12
  // tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t

  fpdfxq_(&iqnset, &x, &muf2, &pdf[0], &iqnchk); 

  h[Hathor::ABOTTOM]  = pdf[1];
  h[Hathor::ACHARM]   = pdf[2];
  h[Hathor::ASTRANGE] = pdf[3];
  h[Hathor::AUP]      = pdf[4];
  h[Hathor::ADOWN]    = pdf[5];

  h[Hathor::GLUON]    = pdf[6];

  h[Hathor::DOWN]     = pdf[7];
  h[Hathor::UP]       = pdf[8];
  h[Hathor::STRANGE]  = pdf[9];
  h[Hathor::CHARM]    = pdf[10];
  h[Hathor::BOTTOM]   = pdf[11];

  for (int i=0; i<13; i++){
    h[i] /= x;
  }
}

double
H1FitterPdf::GetAlphas(double mu){
  double mu2 = mu*mu;
  int nf, ierr;
  return asfunc_(&mu2,&nf,&ierr);
}

void
H1FitterPdf::InitMember(int i){imember=i;}

int
H1FitterPdf::NumberPdf(void){return 1;}

std::string H1FitterPdf::GetName(void){return(PDFname);}
