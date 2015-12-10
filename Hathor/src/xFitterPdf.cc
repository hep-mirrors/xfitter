#include "../interface/xFitterPdf.h"

#include "get_pdfs.h"
#include <iostream>
#include <vector>

xFitterPdf::xFitterPdf(const std::string str) {
  PDFname = str;
  imember = 0;
}

void
xFitterPdf::GetPdf(double x, double muf, double h[13]){

  const double muf2 = muf*muf;

  std::vector<double> pdf;
  pdf.resize(13);
  // ordering is:
  // 0   1   2   3   4   5    6   7   8   9   10  11  12
  // tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t

  HF_GET_PDFS_UNCACHED_WRAP(&x, &muf2, &pdf[0]); 

  h[AbstractHathor::ATOP]     = pdf[0];
  h[AbstractHathor::ABOTTOM]  = pdf[1];
  h[AbstractHathor::ACHARM]   = pdf[2];
  h[AbstractHathor::ASTRANGE] = pdf[3];
  h[AbstractHathor::AUP]      = pdf[4];
  h[AbstractHathor::ADOWN]    = pdf[5];

  h[AbstractHathor::GLUON]    = pdf[6];

  h[AbstractHathor::DOWN]     = pdf[7];
  h[AbstractHathor::UP]       = pdf[8];
  h[AbstractHathor::STRANGE]  = pdf[9];
  h[AbstractHathor::CHARM]    = pdf[10];
  h[AbstractHathor::BOTTOM]   = pdf[11];
  h[AbstractHathor::TOP ]     = pdf[12];

  for (int i=0; i<13; i++){
    h[i] /= x;
  }
}

double
xFitterPdf::GetAlphas(double mu){
  double mu2 = mu*mu;
  return HF_GET_ALPHAS_WRAP(&mu2); 
}

void
xFitterPdf::InitMember(int i){imember=i;}

int
xFitterPdf::NumberPdf(void){return 1;}

std::string xFitterPdf::GetName(void){return(PDFname);}
