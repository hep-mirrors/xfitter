/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include <time.h>
#include <cstdlib>
#include "Coca.h"

// ======================================================
// void Coca_t::Show(ostream &ostr, int cnt, int start) {
void Coca_t::Show(ostream &ostr) {
  // int ind, n=cnt;
  // if(n < 0) n = size() - start;
  // SetHead("--> " + M_Name);
  ostr << M_CmntStr << M_Head << "\n";
  
  // --- time stamp
  if(show_time) {
    time_t t;
    time(&t);
    char tbuf[20];
    strftime(tbuf, sizeof(tbuf), "%Y-%m-%d %H:%M:%S", localtime(&t));
    ostr << M_CmntStr << " " << tbuf << "\n";
  }

  // int MaxNameLen = NameWidth();
  // string sep(M_Leader);
  // if(!sep.empty()) sep.append(M_afterLeader);
  // for(ind = 0; ind < n; ind++) ostr << sep << setw(MaxNameLen) << left << at(start+ind) << endl;
  ostr << *this;
  
  if(!M_Tail.empty()) ostr << M_CmntStr << M_Tail <<endl;
}

// ======================================================
void Coca_t::ParseList(const TextReader_t& TR, size_t istart, size_t count) {
  string sep(Xstring::white_space + "=&");
  int Np=0;
  size_t j1=istart+count;
  if(count == string::npos || j1 > TR.NLines()) j1 = TR.NLines();
  for(size_t j=istart; j < j1; j++) {
    Xstring tx(TR[j]);
    if(tx.empty()) continue;
    if(M_Leader.empty()) {
      if(tx.BeginsWith(M_CmntStr)) continue;
    }
    else {
      if(tx.BeginsWith(M_Leader)) tx.erase(0, M_Leader.size());
      else if(tx.BeginsWith(M_CmntStr)) continue; //--- BeginsWith("") is always false
      else throw Fiasco("Syntax error in '%s'", TR[j].data());
    }

    string nam = tx.Token(sep);
    // tx = tx.Remainder(sep);
    // tx = tx.Remainder(); //--- only whsp trimmed
    tx.TrimLeft();

    if(tx.empty()) throw Fiasco("No value in '%s'", TR[j].data());
    bool cat = tx[0] == '&';
    tx.TrimLeft(sep);
    if(tx.empty()) throw Fiasco("No value in '%s'", TR[j].data());
    
    cout << flush;
    #define COSTR cout
    if(cat) {
      // if(!M_QuoteChars.empty()) 
      tx.Unquote(M_QuoteChars);
      SetVal(nam, GetString(nam)+tx);
    }
    else {
      if(Exists(nam, mode == Strict)) SetVal(nam, tx, true);
      else {
        if(mode == Ignore) {
          if(M_verbose) COSTR << "IGNORED: " << TR[j] << endl;
          continue;
        }
        else {
          AutoDefine(nam, tx);
          if(M_verbose) COSTR << "ADDED: " << TR[j] << endl;
        }
      }
    }
    
    Np++;
  }
  // emsg(0,"%d parameter(s) read.\n", Np);
  COSTR << flush;
  cout << flush;
  if(M_verbose > 1) cout << Np << " parameter(s) read." << endl;
}

/**
  A file with 'control cards' is organized by lines.
  A line with the first non-blank characters = CommentString is a comment; 
  see the constructor and \c SetCmntStr.
  
  A data line must have the format:
  
  \<Leader\>\<name\>{\<whitespace\> | =}\<value\>
  
  Any whitespace before, between and after the above tokens will be ignored.
  
  Optional args BODpat (Begining-Of-Data pattern) and EODpat (End-Of-Data pattern) are glob-style patterns.
  An empty pattern ("") matches blank line.
  All initial lines not matching BODpat will be skipped before reading the data and the first line matching EODpat will end reading.
*/
// ======================================================
void Coca_t::Read(const char *fname, const char *BODpat, const char *EODpat) {
  if(M_verbose > 1) cout << "\n==> Reading parameters from '" << fname <<"' ..." << endl;
  TextReader_t InTx(fname);
  // InTx.SetLineNumBase(1);
  InTx.Open();
  if(BODpat) {
    Xstring bod = InTx.SkipUntil(BODpat);
    if(bod.Match("?--> *")) {
      bod = bod.substr(5);
      // cout <<"=>"<< bod << endl;
      SetName(bod.Trim());
    }
  }
  InTx.Read(EODpat);
  // InTx.Close();
  // cout << InTx.NLines() << endl;
  ParseList(InTx);
}
