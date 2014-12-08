/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2012--2013
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include "offset.h"
#include "FTNFitPars.h"
// #include "decor.h"

// static OffsetCalc_t *gOffs;
static CovMatrix_t CovTot;
static const int isOffset   = 3;
static vector<int> OffsetIndex;

// ============================================================
int SelectOffsetSources(int sysform[], int n) {
  OffsetIndex.clear();
  for(int j=0; j < n; j++) if(sysform[j] == isOffset) OffsetIndex.push_back(j);
  return OffsetIndex.size();
}

// ============================================================
static int CollectOffsetFits(const string& outdir) {
  START_XCODE
  string fpath(outdir);
  if(!fpath.empty() && *(fpath.end()-1) != '/') fpath.append("/"); 
  CovMatrix_t Cstat, Csyst;
  // string fn(fpath+"statcov_0.txt");
  
  // --- read C_stat from the central fit
  
  if(Cstat.Read(fpath+"statcov_0.txt")) return 1;
  int nVarPar = Cstat.num_rows();
  #ifdef DEBUG_OFFSET
    cout << "CollectOffsetFits nVarPar = " << nVarPar << endl;
  #endif
  if(nVarPar <= 0) {return 2;}
  
  // --- read parameters from the central fit
  
  FTNFitPars_t mfp;
  mfp.Read((fpath+"MI_saved_0.txt").c_str(), "tp");  
  
  OffsetCalc_t gOffs(OffsetIndex.size(), nVarPar);
  int rc = gOffs.FillFromFiles(outdir);
  if(rc) return rc;
  Csyst = gOffs.GetCov();
  CovTot = Cstat+Csyst;

  // double par[nVarPar], err[nVarPar];
  // gOffs.GetParams(par);
  // gOffs.GetErrors(err);
  ofstream covdat((fpath+"offset.save.txt").c_str());
  covdat << "Parameters  " << nVarPar << endl;
  covdat << scientific;
  for(int j=0; j < nVarPar; j++) {
    int par_index = mfp.VarIndex(j);
    double tot_err = sqrt(CovTot[j][j]);
    mfp.SetErr(par_index, tot_err); // --- sloppy
    covdat << setw(3) << (j+1)
           << "  "  << mfp.UID(par_index)
           << "  "  << gOffs.Param(j)
           << "  " << tot_err
           << endl;
  }
  covdat << endl;
  Cstat.ShowL(covdat, "Stat. Covariance");
  Csyst.ShowL(covdat, "Syst. Covariance");
  CovTot.ShowL(covdat, "Covariance");
  covdat << fixed;
  covdat.precision(3);
  // covdat.width(6);
  CovTot.ShowL(covdat, "Correlation", true);
  covdat.close();
  
  mfp.Write((fpath+"MI_saved_final.txt").c_str(), "p");
  
  // CovTot.deCorrelate(5);
  // cout.unsetf(ios::floatfield);
  // cout << "Decorrelated shifts  " << CovTot.EigShifts << endl;

  // DeCor.SetCov(gOffs.Cov);
  // delete gOffs;
  // gOffs = NULL;
  return rc;
  END_XCODE
}

// --- FORTRAN stub

extern "C" {
  // ===================================
  int probeoffset_(int *sys_form, int *n) {
    return SelectOffsetSources(sys_form, *n);
  }
  
  // ===================================
  void getnoffset_(int *n) {
    *n = OffsetIndex.size();
  }
  
  // ===================================
  int offsetindex_(int *j) {
    return OffsetIndex.at(*j-1)+1;
  }
  
  // ===================================
  int offsetcollect_(const char* outdir) {
    return CollectOffsetFits(outdir);
  }
  
  // ===================================
  void decordiag_(double *dchi2) {
    CovTot.deCorrelate(*dchi2);
    cout << "\n==> Decorrelated (diagonal) shifts for ErrDef = " << *dchi2 << endl;
    cout << CovTot.EigShifts << endl;
    // for(int k=0; k < DeCor.GetN(); k++) cout << (k+1) <<": "<< DeCor.GetEigShift(k) << endl;
  }
  
  // ===================================
  double decorvarshift_(int *iV, int *iE) {
    return CovTot.VarError(*iV-1, *iE-1);
  }
  
  // ===================================
  void cvfillgamma_(int *nSys, int *nData, double* fortranArray, int *nr_) {
    for(vector<int>::iterator it = OffsetIndex.begin() ; it != OffsetIndex.end(); ++it) cout << ' ' << *it;
    cout << endl;
    CovTot.FillGamma(fortranArray, *nSys, *nData, *nr_);
  }
  
}
