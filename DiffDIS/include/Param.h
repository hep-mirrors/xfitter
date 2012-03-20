/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
  \version 2.1.4
  
  Changes:
    2012-03-02
      added SetVal(name, double/int)
_____________________________________________________________*/

#ifndef CLS_PARAM_H_
#define CLS_PARAM_H_

/** @addtogroup gr_coca
 *  @{
 */

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <bitset>

#include "genut.h"
#include "Xstring.h"

#include "dbgtools.h"

using namespace std;

// class ParamVec_t;

/**
  \brief Typed parameter; INT, REAL or STRING.
*/
// oooooooooooooooooooooooooooooooooooooooooooooooooooooo
class Param_t {
  static char QuoteChar[2];
	// static const string NullString;
  #define NullString Xstring::null_string
	string M_name;
	vector<string> M_alias;
  
  bitset<4> M_flags;
  #define f_SET 0
  #define f_UPD 1
  
	int M_i;
	double M_r;
	string M_s;
	enum ValType_t {t_NONE, t_INT, t_REAL, t_STR} M_t;

// protected:
	string comment;
  // static int outwid;

// protected:
public:
  // static bool QuotedStrings;

  /// @cond NO_DOX_
  //================================================================
  bool GetFlag(size_t i) const {
    return M_flags[i];
  }
  /// @endcond

  // =======================================================
	Param_t(const string& name=NullString) {
		M_name.assign(name);
		M_t = t_STR;
		// M_s = "--no-value--";
  }
	
  // =======================================================
	Param_t(const string& name, int val, const string& c=NullString) {
		M_name.assign(name);
		M_t = t_INT;
		M_i = val;
		comment = c;
    M_flags.set(f_SET);
	}
	
  // =======================================================
	Param_t(const string& name, double val, const string& c=NullString) {
		M_name.assign(name);
		M_t = t_REAL;
		M_r = val;
		comment = c;
    M_flags.set(f_SET);
	}

  // =======================================================
	Param_t(const string& name, const string& val, const string& c=NullString) {
		M_name.assign(name);
		M_t = t_STR;
		M_s = val;
		comment = c;
    M_flags.set(f_SET);
	}

  /// Add alternative name.
  // =======================================================
	void AddAlias(const string& a) {
		M_alias.push_back(a);
	}

  // // -------------------------------
  // friend class Error;
  // class Error : public Fiasco {
  // public:
    // // Error(const Param_t* p, const char* s) : Fiasco("Param_t: '%s' is not %s.", p->M_name.c_str(), s) {}
    // Error(const char* s) : Fiasco("Param_t: '%s' is not %s.", this->M_name.c_str(), s) {}
  // };
  
  //================================================================
  int getI() const {
    if(M_t != t_INT) throw Fiasco("Param_t: '%M_s' is not INT.", M_name.c_str());
    // if(M_t != t_INT) throw Error("INT");
    return M_i;
  }

  //================================================================
  double getR() const {
    if(M_t != t_REAL) throw Fiasco("Param_t: '%s' is not REAL.", M_name.c_str());
    return M_r;
  }

  //================================================================
  string getS() const {
    if(M_t != t_STR) throw Fiasco("Param_t: '%s' is not STR.", M_name.c_str());
    return M_s;
  }

  /// Quoted strings in output
  // =======================================================
	static void SetOutQuote(char qc) {QuoteChar[0] = qc;}
	static char GetOutQuote() {return QuoteChar[0];}

  // =======================================================
  friend ostream &operator<<(ostream &ostr, const Param_t& a) {
    // if(!a.comment.empty()) ostr << "; "<< a.comment << ":\n";
    // ostr << setw(outwid)<<left << a.M_name <<" ";
    ostr << a.M_name << " ";
    // ostr.setf(ios::showpoint);
    switch(a.M_t) {
      case t_INT: ostr << a.M_i; break;
      case t_REAL:
        {
        // stringstream buf(stringstream::in);
        stringstream buf;
        buf << a.M_r;
        if(buf.str().find_first_of(".eE") == string::npos) buf << ".0";
        ostr << buf.str();
        }
        // ostr << a.M_r;
        break;
      case t_STR: 
        ostr << Param_t::QuoteChar << a.M_s << Param_t::QuoteChar;
        // if(QuotedStrings) ostr <<"'"<< a.M_s <<"'";
        // else ostr << a.M_s;
        break;
    }
    #ifdef _DEBUG_
    ostr << " @" << a.M_t;
    #endif
    return ostr;
    // return ostr << endl;
  }
  
  // ..........................
  friend class ParamVec_t;
  // class ParamVec_t;
  // friend ostream& ParamVec_t::operator<<(ostream &ostr, const ParamVec_t& a);

};

/**
  \brief A list of typed parameters.
*/
// oooooooooooooooooooooooooooooooooooooooooooooooooooooo
class ParamVec_t : public vector<Param_t> {
  static string BadChars;
  
	//================================
	void CheckLegalName(const string& name) {
    if(!BadChars.empty() && name.find_first_of(BadChars) != string::npos)
      throw Fiasco("Illegal name: '%s'", name.c_str());
  }

protected:
  string M_QuoteChars;
  string M_OutLeader;
  bool M_OutAligned;

  public:
  
    enum mode_t {m_NONE, m_DEF, m_NEW, m_SET};
    /*
      create or overwrite
      must not exist
      must exist
    */
    
  private:
    mode_t M_mode;
    
    // =======================================================
    int NameWidth() const {
      ParamVec_t::const_iterator iter;
      // vector<Param_t>::iterator iter;
      int w,wid=0;
      for(iter = begin(); iter != end(); iter++) if((w = iter->M_name.length()) > wid) wid = w;
      // for(int i=0; i < size(); i++) if((w = at(i).M_name.length()) > wid) wid = w;
      return wid;
    }
  
  public:
  
    // -------------------------------
    // friend class Error;
    // class Error : public Fiasco {
    // public:
      // Error(const char* t, const char* s) : Fiasco("ParamVec_t: %s '%s'.", t, s) {}
    // };
  
    //================================================
    int FindParam(const string& name, bool strict=true) const {
      int ind, last=size();
      vector<string>::const_iterator iter;
      for(ind = 0; ind < last; ind++) {
        // if(AllParams[ind].M_name.compare(name)==0) return ind;
        if(at(ind).M_name == name) return ind;
        for(iter = at(ind).M_alias.begin(); iter != at(ind).M_alias.end(); iter++)
          if(*iter == name) return ind;
      }
      if(strict) throw Fiasco("Unknown parameter: '%s'", name.c_str());
      return -1;
    }
    
    // ===================================================
    ParamVec_t() {
      M_mode = m_DEF;
      M_OutAligned = false;
    }
    
    //=====================================
    void SetMode(mode_t m) {M_mode = m;}
    
    //=====================================
    void SetQuoteChars(const string& qc=Xstring::quote_chars) {
      M_QuoteChars = qc;
      if(M_QuoteChars.empty()) Param_t::QuoteChar[0] = '\0';
      else Param_t::QuoteChar[0] = M_QuoteChars[0];
      // if(!M_QuoteChars.empty()) Param_t::SetOutQuote(M_QuoteChars[0]);
    }
    
    //=====================================
    void SetOutLeader(const string& qc) {M_OutLeader = qc;}
    
    //=====================================
    void SetOutAligned(bool oa=true) {M_OutAligned = oa;}
    
    //=====================================
    // bool QuotedStrings() const {return !M_QuoteChars.empty();}
    
    //=============================================
    template <typename T>
    Param_t& Define(const string& name, T val, mode_t m1 = m_NONE, bool set_val=true) {
      CheckLegalName(name);
      // mode_t morg = M_mode;
      // if(m1 != m_NONE) SetMode(m1);
      int ind = FindParam(name, m1 == m_SET);
      if(ind < 0) {
        ind = size();
        push_back(Param_t(name, val));
      }
      else {
        if(m1 == m_NEW) throw Fiasco("'%s' already defined.", name.c_str());
        at(ind) = Param_t(name, val);
      }
      if(!set_val) at(ind).M_flags.reset();
      return at(ind);
    }
      
    //=============================================
    template <typename T>
    Param_t& Create(const string& name, T val, bool set_val=true) {return Define(name, val, m_NEW, set_val);}
      
    //=============================================
    void Set(const ParamVec_t& pvec) {
      ParamVec_t::const_iterator iter;
      int w,wid=0;
      for(iter = pvec.begin(); iter != pvec.end(); iter++) {
        int ind = FindParam(iter->M_name, false);
        if(ind < 0) {
          ind = size();
          push_back(*iter);
        }
        else {
          at(ind) = *iter;
        }
        // at(ind).M_flags.reset();
      }
    }
      
    //=============================================
    void Delete(const string& name) {
      // ignore non-existing
      int ind = FindParam(name, false);
      if(ind >= 0) erase(begin()+ind);
    }
      
    //=============================================
    bool Exists(const string& name, bool strict=false) {
      return(FindParam(name, strict) >= 0);
    }
      
    //================================================================
    void AutoDefine(const string& name, const string& tx_in, mode_t m1 = m_NONE) {
      // --- type set acc. to the value of tx_in
      Xstring tx(tx_in);
      // if(!M_QuoteChars.empty()) tx.Unquote(M_QuoteChars);
      if(tx.IsInt()) Define(name, tx.GetInt(), m1);
      else if(tx.IsReal()) Define(name, tx.GetReal(), m1);
      else Define(name, tx.Unquote(M_QuoteChars), m1);
    }

    //================================================================
    void SetVal(int ind, const string& tx_in, bool quoted=false) {
      // --- the param type is checked first and cannot be changed
      Xstring tx(tx_in);
      switch(at(ind).M_t) {
        case Param_t::t_INT:
          at(ind).M_i = tx.GetInt();
          break;
        case Param_t::t_REAL:
          at(ind).M_r = tx.GetReal();
          break;
        case Param_t::t_STR:
          // if(!M_QuoteChars.empty()) tx.Unquote(M_QuoteChars);
          // at(ind).M_s = quoted? tx.Token(tx.substr(0,1)) : tx;
          at(ind).M_s = tx.Unquote(M_QuoteChars);
          break;
      }
      at(ind).M_flags.set(f_SET);
    }

    //================================================================
    void SetVal(const string& name, const string& tx, bool quoted=false) {
      // --- must exist
      SetVal(FindParam(name), tx, quoted);
    }

    //================================================================
    void SetVal(int ind, double val) {
      stringstream vs(ios_base::out);
      switch(at(ind).M_t) {
        case Param_t::t_INT:
          at(ind).M_i = val;
          break;
        case Param_t::t_REAL:
          at(ind).M_r = val;
          break;
        case Param_t::t_STR:
          vs << val;
          at(ind).M_s = vs.str();
          break;
      }
      at(ind).M_flags.set(f_SET);
    }

    //================================================================
    void SetVal(const string& name, double val) {
      // --- must exist
      SetVal(FindParam(name), val);
    }

    //================================================================
    void SetVal(int ind, int val) {
      stringstream vs(ios_base::out);
      switch(at(ind).M_t) {
        case Param_t::t_INT:
          at(ind).M_i = val;
          break;
        case Param_t::t_REAL:
          at(ind).M_r = val;
          break;
        case Param_t::t_STR:
          vs << val;
          at(ind).M_s = vs.str();
          break;
      }
      at(ind).M_flags.set(f_SET);
    }

    //================================================================
    void SetVal(const string& name, int val) {
      // --- must exist
      SetVal(FindParam(name), val);
    }

    //================================================================
    int GetInt(const string& name) const {
      return at(FindParam(name)).getI();
    }

    //================================================================
    double GetReal(const string& name) const {
      return at(FindParam(name)).getR();
    }

    //================================================================
    string GetString(const string& name) const {
      return at(FindParam(name)).getS();
    }

    //================================================================
    Param_t& operator[] (const string& name) {return at(FindParam(name));}
    const Param_t& operator[] (const string& name) const {return at(FindParam(name));}

    //================================================================
    Param_t& operator[] (int ind) {return at(ind);}
    const Param_t&  operator[] (int ind) const {return at(ind);}

    // =======================================================
    friend ostream &operator<<(ostream &ostr, const ParamVec_t& pv) {
      // ios_base::fmtflags orgflags = ostr.setf(ios::showpoint);
      int wid = pv.M_OutAligned ? pv.NameWidth() : 0;
      for(int i=0; i < pv.size(); i++) if(pv[i].GetFlag(f_SET)) ostr << pv.M_OutLeader << setw(wid) << left << pv[i] << endl;
      // ostr.flags(orgflags);
      return ostr;
    }
  
};

// #undef real_type

///@}

#endif
