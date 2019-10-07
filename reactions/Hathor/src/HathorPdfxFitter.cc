#include "HathorPdfxFitter.h"
#include "ReactionHathor.h"
#include <iostream>
#include <valarray>

HathorPdfxFitter::HathorPdfxFitter(ReactionHathor *ptrReactionTheory)
{
  PDFname = "xFitterPdf";
  imember = 0;
  _reactionTheory = ptrReactionTheory;
  IsValid = false;
}

void HathorPdfxFitter::GetPdf(double x, double muf, double h[13])
{
  /*static double xmin = 1.0;
  static double xmax = 0.0;
  if(x < xmin)
    xmin = x;
  if(x > xmax)
    xmax = x;
  printf("HATHORGetPdf x = %f  [%12.8f%12.8f]\n", x, xmin, xmax);*/

  if(!IsValid)
    return;

  std::valarray<double> pdf(0.0, 13);
  //_reactionTheory->xfx(x, muf, &pdf[0]);
  pdf_xfxq_wrapper_(x, muf, &pdf[0]);
  for(auto& val : pdf)
    val /= x;

  // ordering is:
  // 0   1   2   3   4   5    6   7   8   9   10  11  12
  // tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t

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
}

double HathorPdfxFitter::GetAlphas(double mu)
{
  if(!IsValid)
    return 0.0;

  //return _reactionTheory->alphaS(mu);
  return alphas_wrapper_(mu);
}

void HathorPdfxFitter::InitMember(int i){ imember=i; }

int HathorPdfxFitter::NumberPdf(void){ return 1; }

std::string HathorPdfxFitter::GetName(void){ return(PDFname); }
