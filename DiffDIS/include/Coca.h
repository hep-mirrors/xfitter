/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
  \version 2.1.4
  
  Changes:
    2012-03-02
      replace ReadSaved
      by ReadNamed
_____________________________________________________________*/

#ifndef CLS_COCA_H_
#define CLS_COCA_H_

/** @defgroup gr_coca Coca
  Control-cards reader/writer
  @{
*/

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
// #include <string.h>

// #include "emsg.h"
#include "genut.h"
#include "Param.h"
#include "TextReader.h"

// #define WHSP Xstring::white_space

// #define _DEBUG_COCA_
// #undef SHOW
// #ifdef _DEBUG_COCA_
// #define SHOW(a) cout << #a" = '" << (a) <<"'"<< endl;
// #else
// #define SHOW(a)
// #endif

using namespace std;

/// Control-cards file io
// oooooooooooooooooooooooooooooooooooooooo
class Coca_t : public ParamVec_t {
  public:
    enum mode_t {Strict, Auto, Ignore};

  // --- DATA
  // ____________
   
  private:
    mode_t mode;
    // char M_CmntChar;
    string M_CmntStr;
    string M_Leader;
    int M_verbose;
    
    // --- output formatting
    string M_Name, M_Head, M_Tail;
    string M_afterLeader;
    // int M_MaxNameLen;
    bool show_time;

  // --- METHODS
  // ______________
   
  private:
    void F_Init(const string& name="", const string& ldr="") {
      M_Name = name;
      SetLeader(ldr);
      SetVerbose(1);
      SetMode(Strict);
      SetTimeStamp();
      SetHead("-->  test1 " + M_Name);
      SetTail("--- end-of-params");
    }
    
    
  public:
    void ParseList(const TextReader_t& TR, size_t istart=0, size_t count=string::npos);
    
    // ===================================================
    Coca_t(const string& name, char cc, const string& ldr="") : ParamVec_t() {
      SetCmntChar(cc);
      F_Init(name, ldr);
    }
    
    Coca_t(const string& name="", const string& cc="#", const string& ldr="") : ParamVec_t() {
      SetCmntStr(cc);
      F_Init(name, ldr);
      // M_Name = name;
      // SetLeader(ldr);
      // SetVerbose(1);
      // SetMode(Strict);
      // SetTimeStamp();
      // SetHead("--> " + M_Name);
      // SetTail("--- end-of-params");
    }
    
    //================================
    //void IgnoreEmptyLines(bool b=true) {StopOnBlank = !b;}

    //================================
    void SetVerbose(int v) {M_verbose = v;}
    void SetMode(mode_t m) {mode = m;}
    void SetCmntChar(char c) {M_CmntStr = c;}
    void SetCmntStr(const string& s) {M_CmntStr = s;}
    void SetName(const string& name) {
      M_Name = name;
      SetHead("--> " + M_Name);
    }
    
    //================================
    void SetLeader(const string& s, const string& whsp="") {
      M_Leader = s; M_afterLeader = whsp;
      SetOutLeader(M_Leader);
      if(!M_Leader.empty()) SetOutLeader(M_Leader+M_afterLeader);
    }
    
    // bool IgnoreEmptyLines(bool yes=true) {ign_empty_line = yes;}
    void SetHead(const string& t) {M_Head = t;}
    void SetTail(const string& t) {M_Tail = t;}
    void SetTimeStamp(bool y=true) {show_time = y;}
    
    //=========================================
    string GetName() const {return M_Name;}
    string GetCmntStr() const {return M_CmntStr;}
    string GetLeader() const {return M_Leader;}
    
    //================================================================
    // void Show(ostream &ostr=cout, int cnt=-1, int start=0);
    void Show(ostream &ostr=cout);

    //================================================================
    void Save(const char *fname=NULL) {
      const char* fn = fname? fname : M_Name.c_str();
      // string fn(fname? fname : M_Name.c_str());
      if(!*fn) throw Fiasco("Cannot save unnamed Coca_t to unnamed file.");
      // if(!*fn) fn = tmpnam(NULL);
      // if(fn.empty()) fn = tmpnam(NULL);
      // if(*fn == '\\' || *fn == '/') fn++;
      // cout << "fn '"<< fn <<"'"<< endl;
      ofstream f(fn);
      Show(f);
      f.close();
    }

    //================================================================
    void Read(const char *fname, const char *BODpat=NULL, const char *EODpat=NULL);

    /// Read using current Head and Tail for BODpat and EODpat
    //================================================================
    void ReadSaved(const char *fname) {ReadNamed(fname);}
    
    void ReadNamed(const char *fname) {
      Read(fname, (M_CmntStr+M_Head).c_str(), (M_CmntStr+M_Tail).c_str());
    }
};

///@}

#endif
