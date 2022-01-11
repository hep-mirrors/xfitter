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

/// List of \c Xstring objects.
// oooooooooooooooooooooooooooooooooooooooooooooooooo
class StringList : public vector<Xstring> {
public:
  Xstring Join(const string& sep=" ") const;
  void Append(const StringList& slist);
  
  // ==============================================
  StringList Range(size_t start, int count) const;
  
  // =============================================
  size_t Index(const string& glob_patt, bool case_sensitive = true);
  
  // =======================================================
  friend ostream& operator<<(ostream &ostr, const StringList& a);
  
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

  bool IsBlank() const;
  bool IsInt() const;
  int GetInt() const;
  bool IsReal() const;
  double GetReal() const;
  Xstring& Trim(const string& chset=Xstring::white_space);
  Xstring& TrimLeft(const string& chset=Xstring::white_space);
  Xstring& TrimRight(const string& chset=Xstring::white_space);
  bool Contains(const string& str) const;
  bool EndsWith(const string& str) const;
  bool BeginsWith(const string& str) const;
  char BeginsWithChr(const string& chset=quote_chars) const;
  Xstring& Unquote(const string& chset=quote_chars);
  Xstring& UnquoteFirst();
  Xstring Token(const string& sep=white_space);
  Xstring Token(char ch) {return Token(string(1,ch));}
  StringList Split(const string& sep=white_space) const;
  /// Split quoted strings into a list
  StringList SplitQuoted(const string& sep = white_space, const string& quot = quote_chars) const;
  vector<int> SplitInt(const string& sep = white_space) const;
  vector<double> SplitReal(const string& sep = white_space) const;
  Xstring& ToUpper();
  Xstring& ToLower();
  Xstring& Replace(const string& sold, const string& snew);
  Xstring& Replace(const string& sold, char c) {return Replace(sold, string(1,c));}
  Xstring Justified(int typ, size_type wid) const;
  bool Match(const string& glob_patt, bool case_sensitive = true) const;
private:
  void Unquote1(char q);
  bool Match_qm(const string& qs, size_type start) const;
};

#endif
