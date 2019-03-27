#pragma once

#include <string>

// C wrapper for base fortran functions

extern"C" {
  //Error logging function
  void hf_errlog_(const int &id, const char text[], int);
}


// Some basic other functions
int OrderMap(std::string ord);
void hf_errlog(int id, const std::string& message);
std::string stringFromFortran(char*fortran_string,size_t size);
