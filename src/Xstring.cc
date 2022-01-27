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
  size_type first, lastpos;
  // --- skip leading sep
  first = find_first_not_of(sep);
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
  int first;
  if(start == string::npos) start = lastpos;
  // --- skip leading sep
  first = find_first_not_of(sep, start);
  lastpos = string::npos;
  if(first == string::npos) return "";
  lastpos = find_first_of(sep, first);
  return substr(first, lastpos == string::npos ? lastpos : lastpos - first);
}
#endif

//==============================================
bool Xstring::IsBlank() const {
  return find_first_not_of(white_space) == string::npos;
}

//==============================================
bool Xstring::IsInt() const {
  if(empty()) return false;
  int v;
  char c;
  return 1 == sscanf(c_str(), "%d%c", &v, &c);
}

//==============================================
int Xstring::GetInt() const {
  if(empty()) throw Fiasco("*** ERROR in Xstring::GetInt: empty string.");
  int v;
  char c;
  if(1 != sscanf(c_str(), "%d%c", &v, &c)) throw Fiasco("*** ERROR in Xstring::GetInt: %s is not int.", c_str());
  return v;
}

  //==============================================
bool Xstring::IsReal() const {
  if(empty()) return false;
  double v;
  char c;
  return 1 == sscanf(c_str(), "%lf%c", &v, &c);
}

//==============================================
double Xstring::GetReal() const {
  if(empty()) throw Fiasco("*** ERROR in Xstring::GetReal: empty string.");
  double v;
  char c;
  if(1 != sscanf(c_str(), "%lf%c", &v, &c)) throw Fiasco("*** ERROR in Xstring::GetReal: %s is not real.", c_str());
  return v;
}

//==============================================
Xstring& Xstring::TrimLeft(const string& chset) {
  if(!empty()) {
    int first = find_first_not_of(chset);
    erase(0, first);
  }
  return *this;
}

//==============================================
Xstring& Xstring::TrimRight(const string& chset) {
  if(!empty()) {
    int last = find_last_not_of(chset);
    erase(last+1);
  }
  return *this;
}

//==============================================
Xstring& Xstring::Trim(const string& chset) {
  return TrimRight(chset).TrimLeft(chset);
}

//==============================================
Xstring& Xstring::UnquoteFirst() {
  if(empty()) return *this;
  Unquote1(at(0));
  return *this;
}

//==============================================
Xstring& Xstring::Unquote(const string& chset) {
  if(empty() || chset.empty()) return *this;
  char qc;
  if((qc = BeginsWithChr(chset))) Unquote1(qc);
  return *this;
}

/// returns first char if contained in \c chset, 0 otherwise
//==============================================
char Xstring::BeginsWithChr(const string& chset) const {
  if(empty()) return '\0';
  size_t p = chset.find(at(0));
  return  p == npos ? '\0' : chset[p];
}

//==============================================
bool Xstring::BeginsWith(const string& str) const {
  return !(str.empty() || compare(0, str.size(), str));
}

//==============================================
bool Xstring::EndsWith(const string& str) const {
  size_t len = str.size();
  return size() >= len && compare(size() - len, len, str) == 0;
}

//==============================================
bool Xstring::Contains(const string& str) const {
  return find(str) != string::npos;
}


/**
   \brief Does \c *this match \c glob_patt?
   
   \c glob_patt is a glob-style pattern.
   \c *this cannot contain \001 or \002
*/
//==============================================
bool Xstring::Match(const string& glob_patt, bool case_sensitive) const {
  Xstring Ms(glob_patt);
  if(!case_sensitive) {
    return Xstring(*this).ToLower().Match(Ms.ToLower());
  }
  if(glob_patt.find_first_of("*?") == string::npos) return *this == glob_patt;
  
  const char STAR='\001';
  const char QUOT='\002';
  Ms.Replace("*",STAR);
  Ms.Replace(string("\\")+STAR,'*');
  Ms.Replace("?",QUOT);
  Ms.Replace(string("\\")+QUOT,'?');
  
  StringList Ls = Ms.Split(string(1,STAR));
  size_type cp=0;
  size_t cStar=0;
  int st, maxstart;
  for(cStar=0; cStar < Ls.size(); cStar++) {
    maxstart = size() -cp - Ls[cStar].size();
    if(!cStar && Ms[0] != STAR) maxstart = 0;
    for(st=0; st <= maxstart; st++) if(Match_qm(Ls[cStar],cp+st)) break;
    if(st > maxstart) return false;
    cp += st+Ls[cStar].size();
  }
  return (Ms.end()[-1] == STAR) ? cp <= size() : cp == size();
  
  if(Ms[0] != STAR) {
    size_type p1 = find(Ls[cStar]);
    if(p1 == string::npos) return false;
    cp = p1 + Ls[cStar].length();
  }
}



//==============================================
void Xstring::Unquote1(char q) {
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

//==============================================
bool Xstring::Match_qm(const string& qs, size_type start) const {
  const char QUOT='\002';
  if(qs.find(QUOT) == string::npos) return !compare(start, qs.size(), qs);
  string S = substr(start);
  if(S.size() < qs.size()) return false;
  
  string Q(qs);
  size_t i,j;
  for(i=j=0; i < qs.size(); i++, j++) if(qs[i] != QUOT && qs[i] != S[i]) return false;
  return true;
}


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




// ==============================================                                                                                                                                
StringList StringList::Range(size_t start=0, int count=-1) const {
  StringList A(*this);
  A.erase(A.begin(),A.begin()+start);
  if(count >= 0) A.erase(A.begin()+count,A.end());
  return A;
}
