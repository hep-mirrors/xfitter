/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2012-2015
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef FTN_FITPARS_H_
#define FTN_FITPARS_H_

#include <iostream>
#include "FitPars_base.h"

// #define TST_SHOW(a) cout << #a" = '" << (a) <<"'"<< endl;


// -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
class FTNFitPars_t : public FitPars_base_t {
  string M_tit;
  StringList M_cmd;
  ofstream M_trfil;
  int Nchi_tr;
  
  public:
  // ================================================
  ~FTNFitPars_t() {
    if(M_trfil.is_open()) M_trfil.close();
  }
  
  // ================================================
  string ReadTitle(const char * fnsave) {
    M_tit.erase();
    TextReader_t tx(fnsave);
    tx.Open();
    tx.Read("");
    // tx.Close();
    int iLine = tx.Find("SET TIT*",0,false);
    if(iLine >= 0) M_tit = tx[iLine+1];
    return M_tit;
  }

  // ================================================
  string GetTitle() {return M_tit;}

  // ================================================
  void SetTitle(const string& tit) {M_tit = tit;}

  // ================================================
  void ReadCommands(const char * fnsave) {
    M_cmd.clear();
    TextReader_t tx(fnsave);
    tx.Open();
    tx.Read();
    int N = tx.NLines();
    int iLine = tx.Find("PARAM*",0,false);
    if(iLine < 0) return;
    while(iLine < N && !(tx[iLine].IsBlank())) iLine++;
    // cout << tx[iLine] << endl;}
    while(iLine < N && (tx[iLine].IsBlank())) iLine++;
    if(iLine == N) return;
    // return tx.Get(iLine);
    for(; iLine < N; iLine++) {
      if(tx[iLine].IsBlank()) continue;
      Xstring c = tx[iLine];
      c.Trim();
      if(c.BeginsWith("*")) continue;
      M_cmd.push_back(c);
      if(c.Match("return",false)) break;
    }
  }

  // ================================================
  void Read(const char* fn, Xstring what="tpc") {
    // START_XCODE
    if(M_verbose) cout << "===> Reading ["<< what << "] from '" << fn <<"'"<< endl;
    if(what.ToLower().Contains("t")) ReadTitle(fn);
    if(what.ToLower().Contains("p")) {
      ReadMnSaved(fn);
      if(M_verbose) cout << M_pars.size() << " parameters read."<< endl;
    }
    if(what.ToLower().Contains("c")) {
      ReadCommands(fn);
      if(M_verbose) cout << M_cmd.size() << " commands read."<< endl;
    }
    // END_XCODE
  }

  // ================================================
  StringList GetCommands() {
    return M_cmd;
  }

  // ================================================
  void SetCommands(const Xstring& cmds) {
    M_cmd = cmds.Split(";\n");
  }

  // ================================================
  bool HasCommand(const Xstring& cmd) {
    return M_cmd.Index(cmd,false) != string::npos;
  }

  // ================================================
  void Write(ostream &ostr=cout, Xstring what="tpc") {
    // cout << "  what: " << what << endl;
    if(what.ToLower().Contains("t")) ostr <<"SET TITLE\n"<< M_tit << "\n";
    if(what.ToLower().Contains("p")) ShowF(ostr);
    if(what.ToLower().Contains("c")) ostr << M_cmd.Join("\n") << endl;
    ostr << endl;
  }

  // ================================================
  void Write(const char* fn, const Xstring& what="tpc") {
    if(M_verbose) cout << "===> Writing ["<< what << "] to '" << fn <<"'"<< endl;
    ofstream out(fn);
    Write(out, what);
    out.close();
  }
  
  // ================================================
  void OpenTrace(const string& fn) {
    if(!fn.empty()) M_trfil.open(fn.c_str());
    // if(M_trfil.is_open()) ShowVarNames(M_trfil);
    // else ShowVarNames(cout);
    Nchi_tr = 0;
  }
  
  // ================================================
  void CloseTrace() {
    M_trfil.close();
  }
  
  // ================================================
  void ShowVNames(int ndf) {
    int nv = GetNVar();
    if(M_trfil.is_open()) {
      M_trfil << "#$NDF = "<< ndf << endl;
      M_trfil << "#---\n";
      M_trfil << "#     "<< setw(15) << "chi2";
      ShowVarNames(M_trfil);
      M_trfil << "#  1";
      for(int j=0; j <= nv; j++) M_trfil << setw(15) << j+2;
      M_trfil << endl;
    }
    else {
      cout << "#     "<< setw(15) << "chi2";
      ShowVarNames(cout);
    }
  }
  
  // ================================================
  void ShowVVals(double chi2) {
    if(M_trfil.is_open()) {
      M_trfil << setw(6) << ++Nchi_tr << setw(15) << chi2;
      ShowVarVals(M_trfil);
      // if(!(Nchi_tr % 50)) M_trfil.flush();
      if(!(Nchi_tr % 50)) Write("minuit.all.txt", "p");
    }
    else {
      cout << setw(6) << ++Nchi_tr << setw(15) << chi2;
      ShowVarVals(cout);
    }
  }
  
  // ==========================================================
  void GetMinuitParams();
  void SetExtraParams();

};

#endif
