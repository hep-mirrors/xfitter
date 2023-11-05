#include"expression_utils.h"
#include<algorithm>
#include <string>
#include <vector>

namespace xfitter {
  // Return a list of all parameter names present in expression, given as string s, excluding duplicates
  // Parameter is defined as string of alphanumeric characters and '_', not beginning with a digit, and not a builtin or "x"
  void extractParameterNames(const string& s, std::vector<std::string>& ret) {
  // fills ret
  static const std::vector<std::string> ignored = {
    "abs", "acos", "asin", "atan", "atan2", "ceil", "cos", "cosh", "e", "exp", "fac",
    "floor", "ln", "log", "log10", "ncr", "npr", "pi", "pow", "sin", "sinh", "sqrt",
    "tan", "tanh", "x"
  };

  const char* p = s.c_str();
  // p is pointer to next character to be read
  while (true) {
    while (true) {
      char c = *p;
      if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c == '_')) break;
      if (c == 0) return;
      ++p;
    }
    const char *b = p;//beginning of substring
    ++p;
    while ((*p >= 'a' && *p <= 'z') || (*p >= 'A' && *p <= 'Z') || (*p == '_') || (*p >= '0' && *p <= '9')) ++p;
    std::string name = std::string(b, p - b);

    //check if this name should be skipped
    if (std::binary_search(ignored.begin(), ignored.end(), name)) goto skip_append;    
    
    //check if this name is already in list
    if (std::find(ret.begin(), ret.end(), name) != ret.end()) goto skip_append;

    //else append
    ret.push_back(name);

  skip_append:
    if(*p == 0) return;
    ++p;
  }
  //unreachable
 }
}
