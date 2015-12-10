#include "FastNLOxFitter.h"
#include "get_pdfs.h"

using namespace std;

extern "C"{
  //void fpdfxq_(int *iset, const double *x, const double *q2, double *pdfs,int *ichk);
  //void evolution_();
  //double asfunc_( double* r2, int* nf  , int* ierr);
  
}


FastNLOxFitter::FastNLOxFitter(string name) : fastNLOReader(name) {
   // --- fastNLO user: if you have your own alpha_s routing in FastNLOxFitter::EvolveAlphas(double,double)
   //     it is convenient to automatically interface it here.
   //     Otherwise the FastNLO alpha_s evolution code Alphas.cc is used, or
   //     you have to call SetAlphasEvolution.
   //     It might be also convenient to make your scale choices here!
   //FillAlphasCache();
}



//______________________________________________________________________________


double FastNLOxFitter::EvolveAlphas(double Q) const {
  // --- fastNLO user: 
  // Implementation of Alpha_s evolution as function of the
  // factorization scale [and alphas(Mz)].
  // here we access the getter method from xfitter
   double mu2 = Q*Q;
   return HF_GET_ALPHAS_WRAP( &mu2 );
}


//______________________________________________________________________________


bool FastNLOxFitter::InitPDF(){
   // --- fastNLO user: 
   //  Initalize PDF parameters if necessary
   //
  return true;
   // It might be necessary that the PDF grid is recalculated/generated.
   //evolution_();
}


//______________________________________________________________________________



vector<double> FastNLOxFitter::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  x-proton and factorisation scale.
   //
   double muf2 = muf*muf;
   vector < double > xfx(14);// 13+1 for photon
   HF_GET_PDFS_WRAP(&xp, &muf2, &xfx[0]);
   return xfx;
}

