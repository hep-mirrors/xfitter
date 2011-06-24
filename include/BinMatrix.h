#ifndef BinMatrix_h
#define BinMatrix_h 1

#include <string>

#include "IntSteps.h"

typedef double**** double4d;

class BinMatrix : public IntSteps
{
 public:
  BinMatrix(){};
  ~BinMatrix();

  BinMatrix(const double *be, const IntSteps* ist):IntSteps(*ist)
  { _beam_en = *be; }

 public:
  double4d BM;

 protected:
  void BuildBM_Z();
  void BuildBM_W();
  double CosthAnIntW(const double &, const double &, const int );

  double _beam_en;

  std::string _var_name;
  int _nbins;
  double *_bins;

 public:
  double getBeamEn() { return _beam_en;}
  void setBins(const std::string&, const int, const double*);
  int getNbins() { return _nbins; }

 private:
  double costhLT(const double &, const double &);
};



#endif
