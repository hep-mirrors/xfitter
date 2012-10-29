/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

/*
N = number of variable parameters
M = number of correlated sources

1. 
Run the 'central' fit, i.e. with correlated sytematics set to zero.
2.
call OffsetInit(M, N)
call OffsetStatCovSet(Cov)
where Cov = covariance matrix for the 'central' fit.
3.
For each mu=-M,...,-1,  1,...,M run the fit and
call OffsetNewParams(p, mu)
where p = fitted variable parameters; double precision p(N)
4.
call OffsetCovGet(Cov)
which fills Cov with the full covariance matrix.
5.
call OffsetQuit()
this discards the Offset calculator object
*/

#include "../include/offset.h"

static OffsetCalc_t *gOffs;


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
    #ifdef DEBUG_OFFSET
      cout << "FillFromFiles nVarPar = " << n << endl;
    #endif
    if(n <= 0) {return 2;}
    gOffs = new OffsetCalc_t(nCSS, n);
    int rc = gOffs->FillFromFiles(outdir);
    if(!rc) {
      double par[n], err[n];
      gOffs->GetParams(par);
      gOffs->GetErrors(err);
      ofstream covdat((fpath+"offset.save.txt").c_str());
      covdat << "Parameters" << endl;
      for(int j=0; j < n; j++) covdat << setw(3) << (j+1) << "  " << scientific << par[j] << "  " << err[j] << endl;
      covdat << endl;
      covdat << "StatCov:\n" << gOffs->StatCov << endl;
      covdat << "CosyCov:\n" << gOffs->CosyCov << endl;
      covdat << "Cov:\n" << gOffs->Cov << endl;
      covdat.close();
    }
    delete gOffs;
    gOffs = NULL;
    return rc;
  }
  
extern "C" {

  // ===================================
  void offsetinit_(int *ncs, int *nvp) {
    gOffs = new OffsetCalc_t(*ncs, *nvp);
  }

  // ===================================
  void offsetquit_() {
    delete gOffs;
  }

  // ===================================
  void offsetstatcovset_(double* cov) {
    gOffs->SetCovStat(cov);
  }

  // ===================================
  void offsetnewparams_(double* p, int *mu) {
    gOffs->NewParams(p, *mu);
  }

  // ===================================
  void offsetcovget_(double* p) {
    gOffs->GetCov(p);
  }

  // ===================================
  int offsetcollect_(int *ncs, const char* outdir) {
    // return gOffs->FillFromFiles(outdir);
    return CollectOffsetFits(*ncs, outdir);
  }

}
