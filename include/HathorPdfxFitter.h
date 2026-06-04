#ifndef Hathor_xFitterPdf_
#define Hathor_xFitterPdf_

#include "HathorPdf.h"

class HathorPdfxFitter : public Pdf
{
protected:

  std::string PDFname;
  int imember;

public:

  HathorPdfxFitter();
  void GetPdf(double x, double muf, double h[13]);
  double GetAlphas(double mu);
  void InitMember(int i);
  int NumberPdf(void);
  std::string GetName(void);
  void Foo();

  // flag which determines if GetPdf() and GetAlphas() methods should call proper routines or return dummy results
  // (Hathor uses these at the very early initialisation step, need protection against QCDNUM not initialized)
  bool IsValid;
};

#endif
