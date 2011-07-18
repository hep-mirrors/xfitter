#ifndef DYcalc_h
#define DYcalc_h 1

#include <string>

#include "IntSteps.h"
#include "BinMatrix.h"
#include "PDFconv.h"

//! Performs the Simpson calculation of DY cross section
/*!
 Inherits form IntSteps, aggregates to BinMatrix and PDFconv
 */

class DYcalc : public IntSteps
{
 public:
  DYcalc(){};
  ~DYcalc();

  //! Proper constructor
  /*!
   Assigns corresponding pointers to y integration function for W and Z
   \param bm  Bin matrix
   \param pc  PDF convolution
   \param ist Integration steps
   */
  DYcalc(BinMatrix *bm, PDFconv* pc, const IntSteps* ist);

  //! Simpson integration in m range
  int Integrate();
 private:
  //! Pointer to rapidity integration
  int (DYcalc::*intY)(const int, double *);
  //! Rapidity integration for W process
  int intY_W(const int , double *);
  //! Rapidity integration for Z process
  int intY_Z(const int , double *);
  //! Rapidity integration for Z process in case of y bins
  int intYbins_Z(const int , double *);

 private:
  //! Bin matrix pointer
  BinMatrix *_bm;
  //! PDF convolution pointer
  PDFconv *_pc;

  //! Array with integration results for bins
  double *_bin_int;
  
 public:
  //! Returns pointer to bin matrix
  BinMatrix* getBM() { return _bm; }
  //! Returns pointer to pdf convolution
  PDFconv* getPC() { return _pc; }
  //! Access to integration results
  /*!
   Standar bin arrangment for Z. For W - even elements Wmins, odd - Wplus
   */
  void getCalcRes(double *);
};

#endif
