/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include <string>
#include <fstream>

// ==========================================
bool FileReadable(const std::string& fn) {
  std::ifstream ifile(fn.c_str());
  ifile.close();
  return !ifile.fail();
}
