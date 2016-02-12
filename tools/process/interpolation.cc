#include "xfitter_cpp.h"
#define LHAPDF6 6

#if LHAPDF_FAMILY == LHAPDF6
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/PDF.h"
#include "LHAPDF/Factories.h"
#endif
#include <iostream>
#include <cstring>
#include <libgen.h>
#include <vector>
#include <cstdlib>
#include <map>

using namespace std;
#if LHAPDF_FAMILY == LHAPDF6
using namespace LHAPDF;
#endif

extern "C" {
  void interpolation(double x, double Q, char* pdfset_path, int n_flavours, int* pdf_flavours, int pdf_number, double* values);
}

#if LHAPDF_FAMILY == LHAPDF6
void interpolation(double x, double Q, char* pdfset_path,int n_flavours, int* pdf_flavours, int pdf_number, double* values){
  
  string pdfset_path_str(pdfset_path);
  string base_name(basename(pdfset_path));

  static PDF* regridded_pdf(NULL);
  
  if ((!regridded_pdf) || (regridded_pdf->memberID() != pdf_number)) {
    regridded_pdf = mkPDF(base_name, pdf_number);
  }
  
  map<int,double> flavours;
  regridded_pdf->xfxQ(x, Q, flavours);

  int i, fl;
  for (i=0;i<n_flavours; i++) {
	fl = pdf_flavours[i];
	values[i]=flavours[fl];
  }

}
#else 
void interpolation(double x, double Q, char* pdfset_path,int n_flavours, int* pdf_flavours, int pdf_number, double* values){
  cerr<< "S: the grid combination is not applicable to LHAPDFv < 6.x" << endl;
  exit(1);
}
#endif
