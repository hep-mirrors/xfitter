#ifndef DYcalc_h
#define DYcalc_h 1

#include <string>

#include "IntSteps.h"
#include "BinMatrix.h"
#include "PDFconv.h"

class DYcalc : public IntSteps
{
 public:
  DYcalc(){};
  ~DYcalc();

  DYcalc( BinMatrix *, PDFconv*, const IntSteps* ist );

  int Integrate();
 private:
  int (DYcalc::*intY)(const int, double *);
  int intY_W(const int , double *);
  int intY_Z(const int , double *);
  int intYbins_Z(const int , double *);

 private:
  BinMatrix *_bm;
  PDFconv *_pc;

  int _nbins;
  double *_bin_int;
  
 public:
  BinMatrix* getBM() { return _bm; }
  PDFconv* getPC() { return _pc; }
  void getCalcRes(double *);
};

#endif
