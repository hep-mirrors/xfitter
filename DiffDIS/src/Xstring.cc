/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
  \version 2.2.1
_____________________________________________________________*/

#include "Xstring.h"

const string Xstring::white_space(" \t\n\r\f");
const string Xstring::quote_chars("'`\"");
const string Xstring::empty_line("\n");
const string Xstring::null_string;

/**
  \brief Chop off the first token acc. to separators \c sep.
  
  Multiple separators act as one. Leading separators are ignored.
  Returns the token.
*/
//==============================================
Xstring Xstring::Token(const string& sep) {
  Xstring tok;
  // static int lastpos;
  size_type first, lastpos;
  //cout << start << endl;
  // if(start == string::npos) start = lastpos;
  // --- skip leading sep
  first = find_first_not_of(sep);
  //cout <<dec<< "first: "<< first << endl;
  // lastpos = string::npos;
  erase(0, first);
  if(first == string::npos) return tok;
  lastpos = find_first_of(sep);
  tok = substr(0, lastpos);
  erase(0, lastpos);
  return tok;
}
#if 0
//==============================================
Xstring Xstring::Token(const string& sep, string::size_type start) {
  // static int lastpos;
  int first;
  //cout << start << endl;
  if(start == string::npos) start = lastpos;
  // --- skip leading sep
  first = find_first_not_of(sep, start);
  //cout <<dec<< "first: "<< first << endl;
  lastpos = string::npos;
  if(first == string::npos) return "";
  lastpos = find_first_of(sep, first);
  return substr(first, lastpos == string::npos ? lastpos : lastpos - first);
}
#endif

/// Split string acc. to separators \c sep. Multiple separators act as one.
//==============================================
StringList Xstring::Split(const string& sep) const {
  StringList lst;
  Xstring s(*this), t;
  t = s.Token(sep);
  if(t.empty()) return lst;
  do {
    lst.push_back(t);
    t = s.Token(sep);
  } while(!t.empty());
  return lst;
}

// #define SHOW(a) cout << #a" = (" << (a) <<")"<< endl;
// #define SHOW(a)

//==============================================
StringList Xstring::SplitQuoted(const string& sep, const string& quot) const {
  Xstring A(*this), tok;
  StringList res;
  // string quot("'\"");
  // string quot(Xstring::quote_chars);
  // string sep(" \t,");
  char c;
  
  A.Trim();
  A.Trim(sep);
  if(A.empty()) return res;
  for(;;) {
    // SHOW(A)
    c = A.BeginsWithChr(quot);
    // SHOW(A.BeginsWithChr(quot))
    if(c) {
      if(A[1] == c) {
        tok = "";
        A.erase(0,2);
      } else {
        tok = A.Token(c);
        A.erase(0,1);
      }
    } else {
      tok = A.Token(sep);
    }
    // SHOW(tok)
    res.push_back(tok);
    if(A.empty()) return res;
    // SHOW(A)
    c = A.BeginsWithChr(sep);
    if(!c) throw Fiasco("Syntax error in SplitQuoted: %s", c_str());
    A.TrimLeft(sep);
  }
}

// ==========================================
vector<int> Xstring::SplitInt(const string& sep) const {
  Xstring A(*this);
  vector<int> res;
  // StringList sl;
  // string sep(" \t,");
  
  A.Trim();
  if(A.empty()) return res;
  StringList sl = A.Split(sep);
  for(int j=0; j < sl.size(); j++) res.push_back(sl[j].GetInt());
  return res;
}

// ==========================================
vector<double> Xstring::SplitReal(const string& sep) const {
  Xstring A(*this);
  vector<double> res;
  // StringList sl;
  // string sep(" \t,");
  
  A.Trim();
  if(A.empty()) return res;
  StringList sl = A.Split(sep);
  for(int j=0; j < sl.size(); j++) res.push_back(sl[j].GetReal());
  return res;
}



// //====================================
// Xstring& Xstring::Join(const StringList& sl, const string& sep) {
  // clear();
  // StringList::const_iterator iter;
  // bool nxt=false;
  // for(iter = sl.begin(); iter != sl.end(); iter++ ) {
    // if(nxt) append(sep);
    // append(*iter);
    // nxt = true;
  // }
  // return *this;
// }

//====================================
Xstring& Xstring::ToUpper() {
  string::iterator iter;
  for(iter = begin(); iter != end(); iter++ )
    *iter = toupper(*iter);
  return *this;
}

//====================================
Xstring& Xstring::ToLower() {
  string::iterator iter;
  for(iter = begin(); iter != end(); iter++ )
    *iter = tolower(*iter);
  return *this;
}

//============================================================
Xstring& Xstring::Replace(const string& sold, const string& snew) {
  size_type pos, start=0;
  while((pos = find(sold, start)) != string::npos) {
    size_type olen = sold.length();
    size_type nlen = snew.length();
    replace(pos, olen, snew, 0, nlen);
    start = pos+nlen;
  }
  return *this;
}

/**
  \c typ  neg.,0, pos. == left,center, right
*/
//============================================================
Xstring Xstring::Justified(int typ, size_type wid) const {
  if(wid <= length()) return *this;
  // int nsp = typ ? (wid - length()) : (wid - length())/2;
  int nsp = wid - length();
  if(typ > 0) return string(nsp,' ') + *this;
  else if(typ < 0) return *this + string(nsp,' ');
  else return string(nsp/2,' ') + *this + string((nsp+1)/2,' ');
}



//====================================
Xstring StringList::Join(const string& sep) const {
  Xstring xs;
  StringList::const_iterator iter;
  bool nxt=false;
  for(iter = begin(); iter != end(); iter++ ) {
    if(nxt) xs.append(sep);
    xs.append(*iter);
    nxt = true;
  }
  return xs;
}

//====================================
void StringList::Append(const StringList& slist) {
  StringList::const_iterator iter;
  for(iter = slist.begin(); iter != slist.end(); iter++ ) {
    push_back(*iter);
  }
}

// =============================================
size_t StringList::Index(const string& glob_patt, bool case_sensitive /* = true */) {
  size_t j, N = size();
  for(j=0; j < N; j++) if(at(j).Match(glob_patt, case_sensitive)) return j;
  return string::npos;
}

  // =======================================================
  ostream& operator<<(ostream &ostr, const StringList& a) {
    int N=a.size();
    // for(int j=0; j < N; j++) ostr << j <<": "<< string(a.at(j)) << endl;
    for(int j=0; j < N; j++) ostr << j <<": '"<< a[j] <<"'"<< endl;
    return ostr;
  }
  
