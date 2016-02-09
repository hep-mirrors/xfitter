/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2011--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/


// #include <map>
// #include <cmath>
#include "FitPars_base.h"


// ================================================
void FitPars_base_t::Show(ostream &ostr, bool full_err) const {
  int j, N = M_pars.size();
  int wName = 0;
  for(j=0; j < N; j++) wName = max(static_cast<int>(M_pars[j].M_name.size()), wName);
  ostr << "\nPARAMETERS C" << endl;
  for(j=0; j < N; j++) {
    ostr << setw(2) << j
      <<" "<< setw(wName+2) << M_pars[j].M_name 
      <<" "<< M_pars[j].FlagString()
      <<"  "<< setw(14) << M_pars[j].M_value;
    if(M_pars[j].IsVar()) {
      ostr <<"  "<< (full_err ? GetError(j) : M_pars[j].Error());
      if(M_pars[j].HasLimits()) {
        ostr <<"  [ ";
        if(M_pars[j].HasLowerLimit()) ostr << M_pars[j].M_lo;
        ostr <<" : ";
        if(M_pars[j].HasUpperLimit()) ostr << M_pars[j].M_up;
        ostr <<" ]";
      }
    }
    ostr << endl;
  }
  ostr << endl;
}

// ================================================
void FitPars_base_t::ShowF(ostream &ostr, bool full_err) const {
  int j, N = M_pars.size();
  int wName = 10;
  // for(j=0; j < N; j++) wName = max(static_cast<int>(M_pars[j].M_name.size()), wName);
  ostr << "PARAMETERS" << endl;
  for(j=0; j < N; j++) {
    ostr << setw(6)<<right << M_pars[j].M_id
      <<" '"<< setw(wName) << left << M_pars[j].M_name.substr(0,wName)
      // <<"'"<< M_pars[j].FlagString()
      <<"'  "<< setw(14)
      // << scientific
      << M_pars[j].M_value << "  ";
    if(M_pars[j].IsVar()) {
      ostr << (full_err ? GetError(j) : M_pars[j].Error());
      if(M_pars[j].HasLimits()) {
        // ostr <<"  [ ";
        ostr <<"  ";
        if(M_pars[j].HasLowerLimit()) ostr << M_pars[j].M_lo;
        // ostr <<" : ";
        ostr <<"  ";
        if(M_pars[j].HasUpperLimit()) ostr << M_pars[j].M_up;
        // ostr <<" ]";
      }
    }
    else ostr << 0;
    ostr << endl;
  }
  ostr << endl;
}

// // ================================================
// void FitPars_base_t::ShowVarVals(ostream &ostr) {
  // int N = M_pars.size();
  // for(int j=0; j < N; j++) {
    // if(!M_pars[j].IsVar()) continue
    // ostr <<" "<< M_pars[j].M_value;
  // }
  // ostr << endl;
// }

// ================================================
// Uses \c M_var2glb.
void FitPars_base_t::ShowVarVals(ostream &ostr) const {
  int N = M_var2glb.size();
  for(int j=0; j < N; j++) {
    ostr << setw(15)
    // <<" "
    << M_pars[M_var2glb[j]].M_value;
  }
  ostr << endl;
}

// ================================================
// Uses \c M_var2glb.
void FitPars_base_t::ShowVarNames(ostream &ostr) const {
  int N = M_var2glb.size();
  // ostr << "# ";
  for(int j=0; j < N; j++) {
    ostr << setw(15)
    // <<" "
    << M_pars[M_var2glb[j]].M_name;
  }
  ostr << endl;
}

// ================================================
int FitPars_base_t::read_FTN(TextReader_t& tx, int ioffs) {
  // TextReader_t tx(fnsave);
  // tx.Open();
  // tx.SkipUntil("PARAMETERS");
  // tx.Read("");
  int N = tx.NLines();
  int id, npar=0;
  // 7'q2        '  0.70472E+00  0.10058E-01
  // 8'q3        '  0.70251E-11  0.10957E-02  0.00000E+00  0.10000E+02
  // 9'q4        '  0.00000E+00  0.00000E+00   

  for(int j=ioffs; j < N; j++) {
    if(tx[j].IsBlank()) break;
    npar++;
    StringList toks = Xstring(tx[j]).Split("' \t");
    // cout << tx[j] <<"\n"
        // << toks.size()
        // <<"  "<< toks[0]
        // << endl;
    bool HasEnum = toks[0].IsInt();
    if(HasEnum) {
      id = toks[0].GetInt();
      toks.erase(toks.begin());
    }
    else id = M_pars.size();
    Xstring name(NewName(toks[0]));
    switch(toks.size()) {
      case 2:
        M_pars.push_back(FitParam_t(id, name, toks[1].GetReal()));
        break;
      case 3:
        M_pars.push_back(FitParam_t(id, name, toks[1].GetReal(), toks[2].GetReal()));
        break;
      case 5:
        M_pars.push_back(FitParam_t(id, name, toks[1].GetReal(), toks[2].GetReal(), toks[3].GetReal(), toks[4].GetReal()));
        break;
      default: throw Fiasco("Illegal PARAMETER specification:\n\"%s\"", tx[j].c_str());
    }
    M_InNames.push_back(name);
  }
  
  return npar;
}

// ================================================
int FitPars_base_t::read_Cws(TextReader_t& tx) {
  int N = tx.NLines();
  // 0  G1       ---L  0.300541  0.0243123  [ 0 :  ] 
  // 4  G5       C---  0
  // 5  q1       --UL  0.150967  0.00510392  [ 0 : 1.98403 ]
  // 6  q2       ----  1.234  0.028723 
  // or
  // 0  G1         0.300541  0.0243123  [ 0 :  ] 
  // 4  G5         0
  // 5  q1         0.150967  0.00510392  [ 0 : 1.98403 ]
  // 6  q2         1.234  0.028723 

  for(int j=0; j < N; j++) {
    Xstring pdat(tx[j]);
    int fLoUp=0;
    double lo=0,up=0;
    if(pdat.Contains("[")) {
      if(!pdat.Match("*[*:*]*")) throw Fiasco("Illegal PARAMETER specification:\n\"%s\"", tx[j].c_str());
      pdat.Replace(":"," : ");
      StringList xdat = pdat.Split("[:]");
      if(xdat.size() != 3) throw Fiasco("Illegal PARAMETER specification:\n\"%s\"", tx[j].c_str());
      pdat = xdat[0].Trim();
      if(!xdat[1].IsBlank()) {fLoUp |= 1; lo = xdat[1].Trim().GetReal();}
      if(!xdat[2].IsBlank()) {fLoUp |= 2; up = xdat[2].Trim().GetReal();}
    }
    StringList toks = pdat.Split("' \t");
    bool HasEnum = toks[0].IsInt();
    if(HasEnum) toks.erase(toks.begin());
    Xstring name(NewName(toks[0]));
    bool HasFlags = !toks[1].IsReal();
    if(HasFlags) toks.erase(toks.begin()+1);
    double val = toks[1].GetReal();
    double err;
    if(toks.size() == 3) err = toks[2].GetReal();
    FitParam_t fp(name, toks[1].GetReal());
    if(toks.size() == 3) {
      fp = FitParam_t(name, toks[1].GetReal(), toks[2].GetReal(), lo, up);
      fp.M_flags.set(FitParam_t::f_LimLo, fLoUp & 1);
      fp.M_flags.set(FitParam_t::f_LimUp, fLoUp & 2);
      // fp.M_flags.set(FitParam_t::f_Fixed, mpi->IsFixed());
    }
    else if(toks.size() != 2) throw Fiasco("Illegal PARAMETER specification:\n\"%s\"", tx[j].c_str());
    M_pars.push_back(fp);
    M_InNames.push_back(name);
  }
  
  return N;
}

// ================================================
int FitPars_base_t::ReadMnSaved(const string& fnsave, const string& old_new) {
  if(!old_new.empty()) SetMapping(old_new);
  if(M_verbose > 1) cout << "\n==> Reading fit parameters from '" << fnsave <<"'..." << endl;
  int N;
  // switch(fmt) {
    // case 'F': N = ReadMinuitFTN(fnsave); break;
    // case 'C': N = ReadMinuitCws(fnsave); break;
    // default: throw Fiasco("ReadMnSaved: unknown format '%c'", fmt);
  // }
  N = read_Minuit_saved(fnsave.c_str());
  // cout << "Format = '" << M_fmt_in <<"'"<< endl;
  if(M_verbose > 1) cout << N << " parameters read." << endl;
  if(!N) throw "No parameters read.";
  
  // --- M_nameList defines the order of parameters
  if(M_nameList.empty()) M_nameList = M_InNames;
  if(N < M_nameList.size()) throw "Too little parameters read.";
  N = M_nameList.size();
  int pind;
  for(int j=0; j < N; j++) {
    pind = M_InNames.Index(M_nameList[j]);
    if(pind == string::npos) throw Fiasco("Parameter \"%s\" not found in \"%s\"", M_nameList[j].c_str(), fnsave.c_str());
    
    // if(M_pars[pind].M_flags[FitParam_t::f_Const]) Add(M_pars[pind].M_name, M_pars[pind].M_value);
    // else {
      // Add(M_pars[pind].M_name, M_pars[pind].M_value, M_pars[pind].M_err);
      // // if(M_pars[pind].M_flags[FitParam_t::f_LimLo]) SetLowerLimit(j, M_pars[pind].M_lo); //--- clears Upper Limit
      // // if(M_pars[pind].M_flags[FitParam_t::f_LimUp]) SetUpperLimit(j, M_pars[pind].M_up); //--- clears Lower Limit
      // switch(M_pars[pind].M_flags[FitParam_t::f_LimLo] + 2*M_pars[pind].M_flags[FitParam_t::f_LimUp]) {
        // case 1: SetLowerLimit(j, M_pars[pind].M_lo); break; //--- clears Upper Limit
        // case 2: SetUpperLimit(j, M_pars[pind].M_up); break; //--- clears Lower Limit
        // case 3: SetLimits(j, M_pars[pind].M_lo, M_pars[pind].M_up); break;
      // }
    // }
  }
  fix_v2g();
  return N;
}

// =============================================================
void FitPars_base_t::Import(const FitPars_base_t& fp, int id_min) {
  int id0=id_min, N = fp.GetN();
  for(int j=0; j < N; j++) {
    M_pars.push_back(fp.M_pars[j]);
    // --- find first free UID
    while(Index(id0) >= 0) id0++;
    M_pars.back().M_id = id0;
    M_InNames.push_back(fp.M_pars[j].M_name);
  }
}

