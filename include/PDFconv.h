#ifndef PDFconv_h
#define PDFconv_h 1

#include <string>

#include "IntSteps.h"

class PDFconv:public IntSteps
{
 public:
  PDFconv(){};
  ~PDFconv();

  PDFconv(const int , const double*, const IntSteps*);

  int interpPDF();

  int (PDFconv::*getPDFconv)(const int, const int, double, const double&,
      double &, double &);
  int getPDFconvZ(const int, const int, double, const double&,
      double &, double &);
  int getPDFconvW(const int, const int, double, const double&,
      double &, double &);

 private:
  int _NINT;
  int **_IDX;

  int _chg_prod;
  double _beam_en;

  static const int _nfl;
  double *_fdef;

 protected:
  double *_XINT;
  double *_Q2INT;
  double **_XFINT;

 private:
  void init();

 public:
  bool isSameBeam(const int chg_prod, const double *beam_en){
    if ( chg_prod == _chg_prod && *beam_en == _beam_en ) return true;
    else return false;
  }
};

#endif
