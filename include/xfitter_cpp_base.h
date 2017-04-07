#pragma once

// C wrapper for base fortran functions

extern"C" {
  //Error logging function
  void hf_errlog_(const int &id, const char text[], int);
}
