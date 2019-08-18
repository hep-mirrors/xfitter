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
void stringToFortran(char* destination, size_t destination_size, const std::string& s);
bool beginsWith(const std::string&str,const std::string&prefix);
bool beginsWith(const std::string&str,const char*       prefix);
bool beginsWith(const char*       str,const char*       prefix);
void stripString(std::string&s);//Remove leading and trailing whitespace in string

bool fileExists(const std::string& fileName);

namespace xfitter{
/*!
\brief Returns path to directory in which all xfitter output should be stored
\details The output directory is determined by "OutputDirectory:" in parameters.yaml
*/
std::string getOutDirName();
}
