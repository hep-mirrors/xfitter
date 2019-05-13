#include"xfitter_steer.h"
#include"BaseEvolution.h"
#include<cmath>
//Functions to access various things from fortran
//PDFs and alpha_s are taken from default evolution
using namespace std;
//PDFs
extern "C" void hf_get_pdfsq_(double const&x,double const&Q,double*pdfs){
  //printf("hf_get_pdfsq_ should not be used: use pdf_xfxq_wrapper_\n");
  xfitter::defaultEvolution->xfxQarray(x,Q,pdfs);//returns through *pdfs
}
extern "C" void hf_get_pdfs_(double const&x,double const&Q2,double*pdfs){
  //printf("hf_get_pdfs_ should not be used: use pdf_xfxq_wrapper_\n");
  xfitter::defaultEvolution->xfxQarray(x,sqrt(Q2),pdfs);//returns through *pdfs
}
//alpha_s
extern "C" double hf_get_alphasq_(double const&Q){
  //printf("hf_get_alphasq_ should not be used: use alphas_wrapper_\n");
  return xfitter::defaultEvolution->getAlphaS(Q);
}
extern "C" double hf_get_alphas_(double const&Q2){
  //printf("hf_get_alphas_ should not be used: use alphas_wrapper_\n");
  return xfitter::defaultEvolution->getAlphaS(sqrt(Q2));
}
