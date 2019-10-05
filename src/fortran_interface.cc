#include"xfitter_steer.h"
#include"BaseEvolution.h"
#include<cmath>
//Functions to access various things from fortran
//PDFs and alpha_s are taken from default evolution
using namespace std;
//PDFs
extern "C" void hf_get_pdfsq_(double const&x,double const&Q,double*pdfs){
  xfitter::defaultEvolutionInstance()->xfxQarray(x,Q,pdfs);//returns through *pdfs
}
extern "C" void hf_get_pdfs_(double const&x,double const&Q2,double*pdfs){
  xfitter::defaultEvolutionInstance()->xfxQarray(x,sqrt(Q2),pdfs);//returns through *pdfs
}
//alpha_s
extern "C" double hf_get_alphasq_(double const&Q){
  return xfitter::defaultEvolutionInstance()->getAlphaS(Q);
}
extern "C" double hf_get_alphas_(double const&Q2){
  return xfitter::defaultEvolutionInstance()->getAlphaS(sqrt(Q2));
}
