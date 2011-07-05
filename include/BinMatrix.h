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

  BinMatrix(const double *be, const IntSteps* ist);

 public:
  double4d BM;

 protected:
  void BuildBM_W_eta();
  void BuildBM_Z_eta();
  void BuildBM_Z_y();
  double CosthAnIntW(const double &, const double &, 
    const double &, const double&);
  double CosthAnIntZ(const double &, const double &, const int,
    const double &, const double&);

  double _beam_en;

 public:
  double getBeamEn() { return _beam_en;}
  int getNbins() { return _nbins; }

 private:
  double costhLT(const double &, const double &);
};



#endif
