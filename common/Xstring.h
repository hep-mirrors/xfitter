/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
  \version 2.2.1
  
  Changes:
  2012-02-15
    Token returns first token and chops it off *this.
    Remainder removed.
  2012-03-31
    Match improved to accept escaped \? and \*
_____________________________________________________________*/

#ifndef X_STRING_H_
#define X_STRING_H_

#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <cstdio>
#include "fiasco.h"

using namespace std;

class Xstring;
// typedef vector<Xstring> StringList;

/// List of \c Xstring objects.
// oooooooooooooooooooooooooooooooooooooooooooooooooo
class StringList : public vector<Xstring> {
public:
  Xstring Join(const string& sep=" ") const;
  void Append(const StringList& slist);
  
  // ==============================================
  StringList Range(size_t start=0, int count=-1) const {
    StringList A(*this);
    A.erase(A.begin(),A.begin()+start);
    if(count >= 0) A.erase(A.begin()+count,A.end());
    return A;
  }
  
  // =============================================
  size_t Index(const string& glob_patt, bool case_sensitive = true);
  
  // =======================================================
  friend ostream& operator<<(ostream &ostr, const StringList& a);
  // {
    // int N=a.size();
    // for(int j=0; j < N; j++) ostr << j <<": "<< string(a.at(j)) << endl;
    // return ostr;
  // }
  
};

/// Extended string class
// ooooooooooooooooooooooooooooooooooooooo
class Xstring : public string {
  size_t lastpos;
public:
  static const string white_space;//(" \t\n\r\f");
	static const string empty_line;
	static const string null_string;
	static const string quote_chars;
  Xstring(const char* s=NULL) : string(s? s : "") {lastpos = string::npos;}
  Xstring(const char* s, size_t n) : string(s,n) {lastpos = string::npos;}
  Xstring(const string& s) : string(s) {lastpos = string::npos;}

  //==============================================
  bool IsBlank() const {
    return find_first_not_of(white_space) == string::npos;
  }

  //==============================================
  bool IsInt() const {
    // if(empty()) return true;
    if(empty()) return false;
    int v;
    char c;
    return 1 == sscanf(c_str(), "%d%c", &v, &c);
  }

  //==============================================
  int GetInt() const {
    // try {
      // if(empty()) throw "empty";// return 0;//
      if(empty()) throw Fiasco("*** ERROR in Xstring::GetInt: empty string.");
      //oz();
      int v;
      char c;
      if(1 != sscanf(c_str(), "%d%c", &v, &c)) throw Fiasco("*** ERROR in Xstring::GetInt: %s is not int.", c_str());
      return v;
    // }
    // catch(const char* msg) {
      // cout <<"*** ERROR in GetInt: '"<< *this <<"' is "<< msg << endl;
      // exit(999);
    // }
  }

  //==============================================
  bool IsReal() const {
    // if(empty()) return true;
    if(empty()) return false;
    double v;
    char c;
    return 1 == sscanf(c_str(), "%lf%c", &v, &c);
  }

  //==============================================
  double GetReal() const {
    // try {
      if(empty()) throw Fiasco("*** ERROR in Xstring::GetReal: empty string.");
      //oz();
      double v;
      char c;
      if(1 != sscanf(c_str(), "%lf%c", &v, &c)) throw Fiasco("*** ERROR in Xstring::GetReal: %s is not real.", c_str());
      return v;
    // }
    // catch(const char* msg) {
      // cout <<"*** ERROR in GetReal: '"<< *this <<"' is "<< msg << endl;
      // exit(999);
    // }
  }

  //==============================================
  Xstring& TrimLeft(const string& chset=Xstring::white_space) {
    if(!empty()) {
      int first = find_first_not_of(chset);
      erase(0, first);
    }
    return *this;
  }

  //==============================================
  Xstring& TrimRight(const string& chset=white_space) {
    if(!empty()) {
      int last = find_last_not_of(chset);
      erase(last+1);
    }
    return *this;
  }

  //==============================================
  Xstring& Trim(const string& chset=white_space) {
    return TrimRight(chset).TrimLeft(chset);
  }

private:
  //==============================================
  void Unquote1(char q) {
    // if(size() < 2) return;
    if(empty()) return;
    if(at(0) != q) return;
    if(size() == 1) {
      clear();
    } else {
      size_type pos = find(q, 1);
      if(pos != npos && (pos+1) != size()) return;
      erase(begin());
      erase(end()-1);
    }
  }

public:
  //==============================================
  Xstring& UnquoteFirst() {
    if(empty()) return *this;
    Unquote1(at(0));
    return *this;
  }

  //==============================================
  Xstring& Unquote(const string& chset=quote_chars) {
    if(empty() || chset.empty()) return *this;
    // for(int i=0; i < chset.size(); i++) if(at(0) == chset[i]) {Unquote1(chset[i]); break;}
    char qc;
    if((qc = BeginsWithChr(chset))) Unquote1(qc);
    return *this;
  }

  /// returns first char if contained in \c chset, 0 otherwise
  //==============================================
  char BeginsWithChr(const string& chset=quote_chars) const {
    if(empty()) return '\0';
    size_t p = chset.find(at(0));
    return  p == npos ? '\0' : chset[p];
  }

  //==============================================
  bool BeginsWith(const string& str) const {
    return !(str.empty() || compare(0, str.size(), str));
    // find(str) == 0;
  }

  //==============================================
  bool EndsWith(const string& str) const {
    size_t len = str.size();
    return size() >= len && compare(size() - len, len, str) == 0;
  }

  //==============================================
  bool Contains(const string& str) const {
    return find(str) != string::npos;
  }

private:
  //==============================================
  bool Match_qm(const string& qs, size_type start) const {
    const char QUOT='\002';
    if(qs.find(QUOT) == string::npos) return !compare(start, qs.size(), qs);
    string S = substr(start);
    if(S.size() < qs.size()) return false;
    char *qp;
    const char *sp = S.c_str();
    string Q(qs);
    int i,j;
    for(i=j=0; i < qs.size(); i++, j++) if(qs[i] != QUOT && qs[i] != S[i]) return false;
    return true;
  }

public:

  /**
    \brief Does \c *this match \c glob_patt?
    
    \c glob_patt is a glob-style pattern.
    \c *this cannot contain \001 or \002
  */
  //==============================================
  bool Match(const string& glob_patt, bool case_sensitive = true) const {
    Xstring Ms(glob_patt);
    // cout << "Match: " << Ms << endl;
    if(!case_sensitive) {
      return Xstring(*this).ToLower().Match(Ms.ToLower());
    }
    if(glob_patt.find_first_of("*?") == string::npos) return *this == glob_patt;
    // const string STAR("\001");
    // const string QUOT("\002");
    const char STAR='\001';
    const char QUOT='\002';
    Ms.Replace("*",STAR);
    Ms.Replace(string("\\")+STAR,'*');
    Ms.Replace("?",QUOT);
    Ms.Replace(string("\\")+QUOT,'?');
    // cout << "Match: " << Ms << endl;
    StringList Ls = Ms.Split(string(1,STAR));
    size_type cp=0;
    int cStar=0;
    int st, maxstart;
    for(cStar=0; cStar < Ls.size(); cStar++) {
      maxstart = size() -cp - Ls[cStar].size();
      if(!cStar && Ms[0] != STAR) maxstart = 0;
      for(st=0; st <= maxstart; st++) if(Match_qm(Ls[cStar],cp+st)) break;
      if(st > maxstart) return false;
      cp += st+Ls[cStar].size();
    }
    return (Ms.end()[-1] == STAR) ? cp <= size() : cp == size();
    // return true;
    
    if(Ms[0] != STAR) {
      size_type p1 = find(Ls[cStar]);
      if(p1 == string::npos) return false;
      cp = p1 + Ls[cStar].length();
    }
    // int starPos = find("*",
    // return find(str) != string::npos;
  }

  //==============================================
  // Xstring Token(const string& sep=white_space, size_type start=0);
  Xstring Token(const string& sep=white_space);
  Xstring Token(char ch) {return Token(string(1,ch));}

  // //==============================================
  // Xstring Remainder(const string& sep=white_space) const {
    // if(lastpos == string::npos || lastpos >= length()) return "";
    // return Xstring(substr(lastpos, string::npos)).TrimLeft(sep);
  // }

  //==============================================
  StringList Split(const string& sep=white_space) const;
  
  /// Split quoted strings into a list
  StringList SplitQuoted(const string& sep = white_space, const string& quot = quote_chars) const;
  vector<int> SplitInt(const string& sep = white_space) const;
  vector<double> SplitReal(const string& sep = white_space) const;

  //==============================================
  // Xstring& Join(const StringList& sl, const string& sep=" ");

  Xstring& ToUpper();
  Xstring& ToLower();
  Xstring& Replace(const string& sold, const string& snew);
  Xstring& Replace(const string& sold, char c) {return Replace(sold, string(1,c));}
  // enum {};
  Xstring Justified(int typ, size_type wid) const;

};

#endif
