/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2013
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef CLS_TEXTREADER_H_
#define CLS_TEXTREADER_H_

/** @addtogroup gr_coca
 *  @{
 */

#include <vector>
#include <bitset>
#include <iostream>
#include <iomanip>
#include <fstream>
// #include "genut.h"
#include "Xstring.h"

#include "dbgtools.h"

using namespace std;

/**
  \brief Read a [part of] text file into a \c StringList of lines
*/
// ooooooooooooooooooooooooooo
class TextReader_t {

  // --- DATA
  // ____________

  protected:
    string M_fname;
    ifstream M_ifil;
    StringList M_text;
    int M_LineNo;
    int LineNumBase; ///< line enumeration base; default = 0
    bitset<4> M_flags;
    bool FailOnError; ///< if true then throw exception else set ErrorCode
    int M_err_code;
    #define f_TRIM_LEFT 0
    #define f_TRIM_RIGHT 1   
  
  // --- METHODS
  // ______________
  
  public:
    // ===================================================
    TextReader_t(const string& fname="") : LineNumBase(0), FailOnError(true), M_err_code(0) {
      DBG_SHOW(M_ifil.is_open())
      // M_flags[f_TRIM_LEFT] = M_flags[f_TRIM_RIGHT] = 1;
      SetTrim();
      M_fname = fname;
      // M_verbose = 0;
    }
    
    // ===================================================
    ~TextReader_t() {
      Close();
    }
    
    /// Set trimming of leading and/or trailing whitespace.
    // ===================================================
    void SetTrim(bool left=1, bool right=1) {
      M_flags[f_TRIM_LEFT]  = left;
      M_flags[f_TRIM_RIGHT] = right;
    }
    
    //================================
    // void SetVerbose(int v) {M_verbose = v;}
    void SetLineNumBase(int lb) {LineNumBase = lb;}

    //================================
    void SetFailOnError(bool foe=true) {FailOnError = foe;}

    //================================
    int ErrorCode() {return M_err_code;}

    /**
      Clear stored text lines
    */
    //================================================================
    void Clear() {M_text.clear();}

    //================================================================
    int NLines() const {return M_text.size();}

    //================================================================
    const Xstring& operator[] (int ind) const {return M_text.at(ind-LineNumBase);} 
    // Xstring& operator[] (int ind) {return M_text.at(ind-LineNumBase);} 

    //================================================================
    int Find (const string& gpat, int start=0, bool case_sensitive = true) {
      for(int j=start; j < M_text.size(); j++) if(M_text[j].Match(gpat, case_sensitive)) return j;
      return -1;
    } 

    //================================================================
    void Open(const char *fname = NULL) {
      Close();
      if(fname) M_fname = fname;
      M_ifil.open(M_fname.c_str());
      if(!M_ifil) {
        if(FailOnError) throw Fiasco("Cannot open file '%s'", M_fname.c_str());
        else {M_err_code = 1; return;}
      }
      M_LineNo = 0;
    }
    
    /**
      Safe close. No error if not opened.
    */
    //================================================================
    void Close() {
      M_err_code = 0;
      DBG_SHOW(M_ifil.is_open())
      if(M_ifil.is_open()) M_ifil.close();
      DBG_SHOW(M_ifil.is_open())
    }
    
    //================================================================
    ifstream& Istream() {
      return M_ifil;
    }
    
    /**
      Skip lines until the last read line matches \c gpat (glob-style pattern).
      An empty \c gpat matches any blank (whitespace) line.
    */
    //================================================================
    Xstring SkipUntil(const char* gpat) {
      if(!gpat) return "";
      Xstring buf;
      // good() returns true if none of the stream's error flags (eofbit, failbit and badbit) are set.
      // fail() checks if either failbit or badbit is set <==> !M_ifil

      while(getline(M_ifil, buf).good()) {
        M_LineNo++;
        if(!*gpat) {if(buf.IsBlank()) break;}
        else if(buf.Match(gpat)) break;
      }
      // if(!M_ifil.good()) throw Fiasco("No lines matched gpat: '%s'", gpat);
      if(!M_ifil.good()) {
        if(FailOnError) throw Fiasco("EOF hit in SkipUntil");
        else {M_err_code = 2;}
      }
      return buf;
    }
    
    /**
      Read and store text lines as long as they do not match \c EODpat (glob-style pattern) or EOF is reached.
      For NULL \c EODpat the file is read until EOF.
      An empty \c EODpat ("") matches any blank (empty or whitespace-only) line.
      
      Return count of read lines.
      
      The leading and/or trailing spaces are trimmed \em after \c EODpat test,
      according to the current setting (default is: trim both) — see \c SetTrim().
    */
    //================================================================
    int Read(const char* EODpat=NULL) {
      // emsg(0,"Trying to read parameters from '%s' ...", fname);
      Xstring buf;
      int n=0;
      while(getline(M_ifil, buf).good()) {
        M_LineNo++;
        n++;
        if(EODpat) {
          if(!*EODpat) {if(buf.IsBlank()) break;}
          else if(buf.Match(EODpat)) break;
        }
        if(M_flags[f_TRIM_LEFT]) buf.TrimLeft();
        if(M_flags[f_TRIM_RIGHT]) buf.TrimRight();
        M_text.push_back(buf);
      }
      // if(M_verbose) cout << n << " line(s) read." << endl;
      return n;
    }
    
    /**
      Read and store \c nL text lines.
      
      Return count of read lines — can be < \c nL if EOF is reached.
      
      The leading and/or trailing spaces are trimmed,
      according to the current setting (default is: trim both) — see \c SetTrim().
    */
    //================================================================
    int Read(int nL) {
      Xstring buf;
      int n=0;
      while(nL-- && getline(M_ifil, buf).good()) {
        M_LineNo++;
        n++;
        if(M_flags[f_TRIM_LEFT]) buf.TrimLeft();
        if(M_flags[f_TRIM_RIGHT]) buf.TrimRight();
        M_text.push_back(buf);
      }
      // if(M_verbose) cout << n << " line(s) read." << endl;
      return n;
    }
    
    /**
      Write stored lines to \c ostream (\c cout by default).
    */
    //================================================================
    void Dump(ostream &ostr=cout) {
      StringList::iterator iter;
      for(iter = M_text.begin(); iter != M_text.end(); iter++) ostr << *iter << endl; 
    }
    
    /**
      Get a subrange of stored lines.
    */
    //================================================================
    StringList Get(size_t start=0, int count=-1) const {
      return M_text.Range(start, count);
    }
    
    // friend class Coca_t;
    friend class DataTable_t;
    
    #undef f_TRIM_LEFT
    #undef f_TRIM_RIGHT   
};

///@}

#endif
