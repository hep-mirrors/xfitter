#pragma once
#include"hf_errlog.h"

// C wrapper for base fortran functions

extern"C" {
  //Error logging function
  void hf_errlog_(const int &id, const char text[], int);
}

// Some basic other functions
int OrderMap(std::string ord);
// Some utility functions for working with strings
std::string stringFromFortran(char*fortran_string,size_t size);
bool beginsWith(const std::string&str,const std::string&prefix);
bool beginsWith(const std::string&str,const char*       prefix);
bool beginsWith(const char*       str,const char*       prefix);
void stripString(std::string&s);//Remove leading and trailing whitespace in string

namespace xfitter{
/*!
\brief Returns path to directory in which all xfitter output should be stored
\details The output directory is determined by setting OutDirName in steering.txt
*/
std::string getOutDirName();
}
