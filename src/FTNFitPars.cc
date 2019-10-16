#include "FTNFitPars.h"
// #include <cstring>

extern "C" {
  // -------------------------------------------------------------------
        // SUBROUTINE MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,ISTAT)
  // CC       Provides the user with information concerning the current status
  // CC          of the current minimization. Namely, it returns:
  // CC        FMIN: the best function value found so far
  // CC        FEDM: the estimated vertical distance remaining to minimum
  // CC        ERRDEF: the value of UP defining parameter uncertainties
  // CC        NPARI: the number of currently variable parameters
  // CC        NPARX: the highest (external) parameter number defined by user
  // CC        ISTAT: a status integer indicating how good is the covariance
  // CC           matrix:  0= not calculated at all
  // CC                    1= approximation only, not accurate
  // CC                    2= full matrix, but forced positive-definite
  // CC                    3= full accurate covariance matrix
  void mnstat_(double *fmin, double *fedm, double *errdef, int *npari, int *nparx, int *istat);
  
  
  // -------------------------------------------------------------------
        // SUBROUTINE MNPOUT(IUEXT,CHNAM,VAL,ERR,XLOLIM,XUPLIM,IUINT)
  // #include "./d506dp.inc"
  // CC     User-called
  // CC   Provides the user with information concerning the current status
  // CC          of parameter number IUEXT. Namely, it returns:
  // CC        CHNAM: the name of the parameter
  // CC        VAL: the current (external) value of the parameter
  // CC        ERR: the current estimate of the parameter uncertainty
  // CC        XLOLIM: the lower bound (or zero if no limits)
  // CC        XUPLIM: the upper bound (or zero if no limits)
  // CC        IUINT: the internal parameter number (or zero if not variable,
  // CC           or negative if undefined).
  // CC  Note also:  If IUEXT is negative, then it is -internal parameter
  // CC           number, and IUINT is returned as the EXTERNAL number.
  // CC     Except for IUINT, this is exactly the inverse of MNPARM
  void mnpout_(int *iuext, char *chnam, double *val, double *err, double *xlolim, double *xuplim, int *iuint, int len);
  
  
  // -------------------------------------------------------------------
  // integer nExtraParamMax   !> Maximum number of parameters
  // parameter (nExtraParamMax=50)

  // integer nExtraParam      !> Actual number of parameters
  // character*32 ExtraParamNames(nExtraParamMax)     !> Names of extra pars
  // double precision ExtraParamValue(nExtraParamMax) !> Initial values
  // double precision ExtraParamStep (nExtraParamMax) !> Initial step size
  // double precision ExtraParamMin  (nExtraParamMax) !> Min value 
  // double precision ExtraParamMax  (nExtraParamMax) !> Max value
  // integer iExtraParamMinuit       (nExtraParamMax) !> Minuit param. index

  // common/ExtraPars/ExtraParamNames,ExtraParamValue,ExtraParamStep,
  // $     ExtraParamMin,ExtraParamMax,iExtraParamMinuit,nExtraParam
  
  const int nExtraParamMax=50;
  struct COMMON_ExtraPars_t {
    char Names[nExtraParamMax][32];
    double Value[nExtraParamMax], Step[nExtraParamMax], Min[nExtraParamMax], Max[nExtraParamMax];
    int iExtraParamMinuit[nExtraParamMax], nExtraParam;
  };
  // --- read_steer.f DOES NOT fill iExtraParamMinuit
    
}

extern COMMON_ExtraPars_t extrapars_;
const int XPbase=101;

static FTNFitPars_t MInput;

// ==========================================================
/**
  @brief Replace all parameters by those defined in Minuit.
  @details Do nothing if the Minuit parameters undefined.
*/
void FTNFitPars_t::GetMinuitParams() {
  double fmin, fedm, errdef;
  int npari, nparx, istat;
  mnstat_(&fmin, &fedm, &errdef, &npari, &nparx, &istat);
  // Provides the user with information concerning the current status
  //    of the current minimization. Namely, it returns:
  //  FMIN: the best function value found so far
  //  FEDM: the estimated vertical distance remaining to minimum
  //  ERRDEF: the value of UP defining parameter uncertainties
  //  NPARI: the number of currently variable parameters
  //  NPARX: the highest (external) parameter number defined by user
  //  ISTAT: a status integer indicating how good is the covariance
  //     matrix:  0= not calculated at all
  //              1= approximation only, not accurate
  //              2= full matrix, but forced positive-definite
  //              3= full accurate covariance matrix
  
  if(!nparx) return;
  
  const int namlen=10;
  int iuext, iuint;
  double val, err, xlolim, xuplim;
  Reset();
  for(int i=1; i <= nparx; i++) {
    char chnam[namlen+8];
    iuext=i;
    mnpout_(&iuext, chnam, &val, &err, &xlolim, &xuplim, &iuint, namlen);
    if(iuint < 0) continue;
    // chnam[namlen] = '\0';
    Xstring name(chnam, namlen);
    name.TrimRight();
    if(!xlolim && !xuplim) AddParam(i, name, val, err);
    else AddParam(i, name, val, err, xlolim, xuplim);
  }
  fix_v2g();
}

// ==========================================================
void FTNFitPars_t::SetExtraParams() {
  for(int i=0; i < extrapars_.nExtraParam; i++) {
    int ind = Index(XPbase+i);
    if(ind < 0) continue;
    extrapars_.Value[i] = Value(ind);
    extrapars_.Step[i] = Error(ind);
    if(M_pars[ind].HasLimits()) {
      extrapars_.Min[i] = M_pars[ind].LowerLimit();
      extrapars_.Max[i] = M_pars[ind].UpperLimit();
    }
  }
}

// --------------------------------------------------------------------------------

// ============================================================
static void RecoverCentralParams(const char* outdir, const char* fn) {
  START_XCODE
  FTNFitPars_t mp;
  string fpath(outdir);
  if(*(fpath.end()-1) != '/') fpath.append("/");  
  mp.ReadMnSaved(fpath+"MI_saved_0.txt");
  mp.SetTitle("Final fit parameters");
  mp.SetCommands("call fcn 3;return"); // make central pdfs
  // mp.Write((fpath+fn).c_str());
  mp.Write(fn);
  mp.SetExtraParams();
  return;
  END_XCODE
}



// ==========================================================
Xstring FTNstring(const char* fn, int fnLEN) {
  // cout << "FTNstring: " << fnLEN <<"\n'" << string(fn, fnLEN) <<"'"<<endl;
  return Xstring(fn,fnLEN).TrimRight();
}



extern "C" {
  void mntinpread_(const char* fn, const char* what, int* verb, int fnLEN, int whatLEN) {
    MInput.SetVerbose(*verb);
    MInput.Read(FTNstring(fn, fnLEN).c_str(), FTNstring(what, whatLEN));
  }
  // void mntinpreadpar_(const char* fn, int fnLEN) {MInput.Read(FTNstring(fn, fnLEN).c_str(), "p");}
  void mntinphascmd_(const char* cmd, int* ans, int cmdLEN) {*ans = MInput.HasCommand(FTNstring(cmd, cmdLEN).c_str());}
  void mntinpsetcmd_(const char* cmd, int cmdLEN) {MInput.SetCommands(FTNstring(cmd, cmdLEN).c_str());}
  void mntinpwrite_  (const char* fn, int fnLEN) {MInput.Write(FTNstring(fn, fnLEN).c_str());}
  void mntinpwritepar_  (const char* fn, int fnLEN) {MInput.Write(FTNstring(fn, fnLEN).c_str(), "p");}
  void mntinpfixxtra_() {MInput.SetExtraParams();}
  void mntinpgetparams_() {MInput.GetMinuitParams();}
  void mntinpsetparams_(const double* vals) {MInput.SetVarValues(vals);}
  void recoverparams_(const char* outdir, const char* fn, int oLEN, int fnLEN) {
    RecoverCentralParams(FTNstring(outdir, oLEN).c_str(), FTNstring(fn, fnLEN).c_str());
  }
  double mntparamvalue_(const char* name, int* idx, int nLEN) {
    *idx = MInput.Index(FTNstring(name, nLEN).c_str());
    return *idx < 0 ? 0 : MInput.Value(*idx);
  }
  void mntchitrace_(const char* fn, int fnLEN) {MInput.OpenTrace(FTNstring(fn, fnLEN));}
  void mnttrclose_() {MInput.CloseTrace();}
  void mntshowvvalues_(double* chi2) {MInput.ShowVVals(*chi2);}
  void mntshowvnames_(int* ndf) {MInput.ShowVNames(*ndf);}
}
