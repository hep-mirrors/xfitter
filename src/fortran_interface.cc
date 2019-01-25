#include"xfitter_steer.h"
#include"BaseEvolution.h"
#include<cmath>
//Functions to access various things from fortran
//PDFs and alpha_s are taken from default evolution
using namespace std;
//PDFs
extern "C" void hf_get_pdfsq_(double const&x,double const&Q,double*pdfs){
  xfitter::defaultEvolution->xfxQArray()(x,Q,pdfs);//returns through *pdfs
}
extern "C" void hf_get_pdfs_(double const&x,double const&Q2,double*pdfs){
  xfitter::defaultEvolution->xfxQArray()(x,sqrt(Q2),pdfs);//returns through *pdfs
}
//alpha_s
extern "C" double hf_get_alphasq_(double const&Q){
  return xfitter::defaultEvolution->AlphaQCD()(Q);
}
extern "C" double hf_get_alphas_(double const&Q2){
  return xfitter::defaultEvolution->AlphaQCD()(sqrt(Q2));
}
