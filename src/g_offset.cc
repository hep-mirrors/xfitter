/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include "offset.h"
#include "FTNFitPars.h"
#include "decor.h"

// static OffsetCalc_t *gOffs;
static Decor_t DeCor;

// ============================================================
static int CollectOffsetFits(int nCSS, const string& outdir) {
  string fpath(outdir);
  if(*(fpath.end()-1) != '/') fpath.append("/");
  ifstream dat;
  int n;
  
  #ifdef DEBUG_OFFSET
    cout << "FillFromFiles opening: " << (fpath+"statcov_0.txt") << endl;
  #endif
  dat.open((fpath+"statcov_0.txt").c_str());
  if(!dat) return 1;
  n = -1;
  dat >> n;
  dat.close();
  
  FTNFitPars_t mfp;
  mfp.Read((fpath+"MI_saved_0.txt").c_str(), "tp");
  // vector<int> VI = mfp.GetVarIndices();
  
  
  #ifdef DEBUG_OFFSET
    cout << "FillFromFiles nVarPar = " << n << endl;
  #endif
  if(n <= 0) {return 2;}
  OffsetCalc_t *gOffs;
  gOffs = new OffsetCalc_t(nCSS, n);
  int rc = gOffs->FillFromFiles(outdir);
  if(!rc) {
    double par[n], err[n];
    gOffs->GetParams(par);
    gOffs->GetErrors(err);
    ofstream covdat((fpath+"offset.save.txt").c_str());
    covdat << "Parameters" << endl;
    for(int j=0; j < n; j++) covdat << setw(3) << (j+1) << "  "  << mfp.UID(mfp.VarIndex(j)) << "  " << scientific << par[j] << "  " << err[j] << endl;
    covdat << endl;
    covdat << "StatCov:\n" << gOffs->StatCov << endl;
    covdat << "CosyCov:\n" << gOffs->CosyCov << endl;
    covdat << "Cov:\n" << gOffs->Cov << endl;
    covdat.close();
  }
  DeCor.SetCov(gOffs->Cov);
  delete gOffs;
  // gOffs = NULL;
  return rc;
}

// --- FORTRAN stub

extern "C" {
  // ===================================
  int offsetcollect_(int *ncs, const char* outdir) {
    return CollectOffsetFits(*ncs, outdir);
  }
  
  // ===================================
  void decordiag_(double *dchi2) {
    DeCor.SetDeltaChiSq(*dchi2);
    DeCor.Diag();
    // --- sort acc. to J. Pumplin convention
    DeCor.Sort();
    cout << "\n==> De-correlated (diagonal) shifts for ErrDef = " << *dchi2 << endl;
    for(int k=0; k < DeCor.GetN(); k++) cout << (k+1) <<": "<< DeCor.GetEigShift(k) << endl;
  }
  
  // ===================================
  double decorvarshift_(int *iV, int *iE) {
    return DeCor.VarError(*iV-1, *iE-1);
  }
  
}