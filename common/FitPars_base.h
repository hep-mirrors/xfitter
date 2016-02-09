/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2011--2015
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef FIT_PARS_BASE_H_
#define FIT_PARS_BASE_H_

#include <map>
#include <cmath>
#include "TextReader.h"

// -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
class FitParam_t {
  // --- FORTRAN Minuit parameter is identified by id
  int M_id;
  double M_value, M_err, M_syserr, M_lo, M_up;
  string M_name;
  enum {f_LimLo, f_LimUp, f_Fixed, f_Const};
  // --- Flags: Const  Fixed  Upper  Lower
  bitset<4> M_flags; //--- initialized with 0 by bitset ctor
  
// protected:
public:
  // ========================================
  bool IsConst() const {return M_flags.test(f_Const);}
  bool IsFixed() const {return M_flags.test(f_Fixed);}
  bool HasLowerLimit() const {return M_flags.test(f_LimLo);}
  bool HasUpperLimit() const {return M_flags.test(f_LimUp);}
  bool HasLimits() const {return M_flags.test(f_LimLo) || M_flags.test(f_LimUp);}
  bool IsVar() const {return !M_flags.test(f_Const) && !M_flags.test(f_Fixed);}
  double LowerLimit() const {return M_lo;}
  double UpperLimit() const {return M_up;}
  double Value() const {return M_value;}
  double Error() const {return M_err;}
  
  // ========================================
  /// @brief Return 4-letter string of active flags.
  string FlagString() const {
    string fs("CFUL");
    string::reverse_iterator rit;
    int k=0;
    for (rit = fs.rbegin(); rit < fs.rend(); rit++, k++) if(!M_flags.test(k)) *rit = '-';
    return fs;
  }
  
  // // ========================================
  // FitParam_t(const string& name, double val) {
    // M_name = name;
    // M_value = val;
    // M_flags.set(f_Const);
  // }
  
  // ========================================
  FitParam_t(const string& name, double val, double step = 0.0) {
    M_id = 0;
    M_name = name;
    M_value = val;
    M_err = step;
    M_syserr = 0;
    if(M_err == 0.0) M_flags.set(f_Const);
  }
  
  // ========================================
  FitParam_t(const string& name, double val, double step, double Lo, double Hi) {
    M_id = 0;
    M_name = name;
    M_value = val;
    M_err = step;
    M_syserr = 0;
    M_lo = Lo;
    M_up = Hi;
    // if(M_err == 0.0) throw
    M_flags.set(f_LimLo);
    M_flags.set(f_LimUp);
  }
  
  // ========================================
  FitParam_t(int id, const string& name, double val, double step = 0.0) {
    M_id = id;
    M_name = name;
    M_value = val;
    M_err = step;
    M_syserr = 0;
    if(M_err == 0.0) M_flags.set(f_Const);
  }
  
  // ========================================
  FitParam_t(int id, const string& name, double val, double step, double Lo, double Hi) {
    M_id = id;
    M_name = name;
    M_value = val;
    M_err = step;
    M_syserr = 0;
    M_lo = Lo;
    M_up = Hi;
    // if(M_err == 0.0) throw
    M_flags.set(f_LimLo);
    M_flags.set(f_LimUp);
  }
  
  friend class FitPars_base_t;
  
};


// -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o
class FitPars_base_t {
protected:
  vector<FitParam_t> M_pars;
  StringList M_nameList;
  StringList M_InNames;
  map<string,string> M_NamesMap;
  Xstring M_fmt_in;
  vector<int> M_var2glb;
  int M_verbose;
  
public:
  // ================================================
  FitPars_base_t(const string& names = "") {
    if(!names.empty()) SetNames(names);
    M_verbose = 1;
  }
  
  // ================================================
  void SetVerbose(int v) {M_verbose = v;}
    
  // ================================================
  /**
    @brief Define order of parameters by name, before their creation.
  */
  void SetNames(const string& names) {
    M_nameList = Xstring(names).Split(" \t,");
  }
  
  // ================================================
  void Reset() {
    M_pars.clear();
    M_InNames.clear();
  }
  
  // ================================================
  Xstring GetFmtIn() const {
    return M_fmt_in;
  }
  
  // // ================================================
  // string FlagString(int j) {
    // ROOT::Minuit2::MinuitParameter p = Parameter(j);
    // string fs("CFUL");
    // int k=0;
    // if(!p.IsConst()) fs[0] = '-';
    // if(!p.IsFixed()) fs[1] = '-';
    // if(!p.HasUpperLimit()) fs[2] = '-';
    // if(!p.HasLowerLimit()) fs[3] = '-';
    // return fs;
  // }
  
  // ================================================
  void Show(ostream &ostr=cout, bool full_err=false) const;
  
  // ================================================
  void ShowF(ostream &ostr=cout, bool full_err=false) const;
  
  // ================================================
  void ShowVarVals(ostream &ostr=cout) const;
  void ShowVarNames(ostream &ostr=cout) const;
  
  // ================================================
  /// @brief Rename some parameters.
  void SetMapping(const string& old_new) {
    StringList renmap = Xstring(old_new).Split(" \t,=");
    if(renmap.size() & 1) throw Fiasco("Uneven # items in the name map list");
    StringList::iterator si;
    for(si= renmap.begin(); si != renmap.end(); si+=2) {
      M_NamesMap[*si] = si[1];
    }
  }
  
  // =============================================================
  void AddParam(int id, const string& name, double val, double step = 0.0) {
    M_pars.push_back(FitParam_t(id, name, val, step));
    M_InNames.push_back(name);
  }
  
  // =============================================================
  void AddParam(int id, const string& name, double val, double step, double Lo, double Hi) {
    M_pars.push_back(FitParam_t(id, name, val, step, Lo, Hi));
    M_InNames.push_back(name);
  }
  
  // =============================================================
  void Import(const FitPars_base_t& fp, int id_min = 0);
  
  // =============================================================
  void Import(const string& fn, int id_min = 0) {
    FitPars_base_t mip;
    mip.ReadMnSaved(fn);
    Import(mip, id_min);
  }
  
private:
  // ================================================
  string NewName(const string& old) const {
    // if(M_NamesMap.count(old)) return M_NamesMap[old];
    if(M_NamesMap.count(old)) return M_NamesMap.at(old);
    return old;
  }
  
  // ================================================
  int read_Minuit_saved(const char * fnsave) {
    // *this = ROOT::Minuit2::MnUserParameters();
    M_pars.clear();
    M_InNames.clear();
    TextReader_t tx(fnsave);
    tx.Open();
    // Xstring hdr(tx.SkipUntil("PARAMETERS*"));
    tx.Read("");
    int iLine = tx.Find("PARAMETERS*",0,false);
    if(iLine < 0) throw Fiasco("Cannot find 'PARAMETERS' in %s", fnsave);
    Xstring hdr = tx[iLine].substr(10);
    M_fmt_in = hdr.Trim();
    if(M_fmt_in.empty() || M_fmt_in[0]=='F') return read_FTN(tx, iLine+1);
    // if(hdr == "C") return read_Cws(tx);
    if(M_fmt_in[0] == 'C') return read_Cws(tx);
    return 0;
  }
  
  // ===>
  int read_FTN(TextReader_t& tx, int ioffs=0);
  int read_Cws(TextReader_t& tx);

public:
  // ===>
  int ReadMnSaved(const string& fnsave, const string& old_new="");

  // // ================================================
  // void Scale(const string& name, double s) {
    // int ind = Index(name);
    // ROOT::Minuit2::MinuitParameter p = Parameter(ind);
    // SetValue(ind, s*p.Value());
    // if(p.IsConst()) return;
    // SetError(ind, s*p.Error());
    // if(p.HasLowerLimit() && p.HasUpperLimit()) SetLimits(ind, s*p.LowerLimit(), s*p.UpperLimit());
    // else if(p.HasLowerLimit()) SetLowerLimit(ind, s*p.LowerLimit());
    // else if(p.HasUpperLimit()) SetUpperLimit(ind, s*p.UpperLimit());
  // }
  
  // ================================================
  int Index(int Id) const {
    int j, N = M_pars.size();
    for(j=0; j < N; j++) if(M_pars[j].M_id == Id) return j;
    return -1;
  }
  
  // ================================================
  /**
  @brief Case sensitive search for a first match.
  
  Same as \c M_InNames.Index(const string& snam) if \c M_NamesMap empty.
  */
  int Index(const string& snam) const {
    int j, N = M_pars.size();
    for(j=0; j < N; j++) if(M_pars[j].M_name == snam) return j;
    return -1;
  }
  
  // ================================================
  void fix_v2g() {
    M_var2glb.clear();
    // vector<ROOT::Minuit2::MinuitParameter> MnPars = Parameters();
    int j, N = M_pars.size();
    for(j=0; j < N; j++) {
      if(M_pars[j].IsVar()) M_var2glb.push_back(j);
    }
    // if(M_var2glb.size() != VariableParameters()) throw Fiasco("Internal ERROR in fix_v2g");
  }
  
  // ================================================
  int GetNVar() {
    fix_v2g();
    return M_var2glb.size();
  }
  
  // ================================================
  int GetN() const {return M_pars.size();}
  
  // ================================================
  vector<int> GetVarIndices() {
    fix_v2g();
    return M_var2glb;
  }
  
  // ================================================
  int VarIndex(int ivar) {
    fix_v2g();
    return M_var2glb.at(ivar);
  }
  
  // ================================================
  string Name(int ig) const {return M_pars[ig].M_name; }
  int UID(int ig) const {return M_pars[ig].M_id; }
  double Value(int ig) const {return M_pars[ig].M_value; }
  double Error(int ig) const {return M_pars[ig].M_err; }
  double SysError(int ig) const {return M_pars[ig].M_syserr; }
  
  // ================================================
  double GetError(int ig) const { return hypot(M_pars[ig].M_err, M_pars[ig].M_syserr); }
  
  // ================================================
  void SetVal(int ig, double v) { M_pars[ig].M_value = v; }
  void SetErr(int ig, double e) { M_pars[ig].M_err = e; }
  void SetSysErr(int ig, double se) { M_pars[ig].M_syserr = se; }
  
  // ================================================
  void SetVarValues(const double* vals) {
    // fix_v2g();
    int N = M_var2glb.size();
    for(int ig=0; ig < N; ig++) M_pars[ig].M_value = *vals++;
  }
  
};


#endif
