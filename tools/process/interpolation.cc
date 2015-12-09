#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/PDF.h"
#include "LHAPDF/Factories.h"
#include <iostream>
#include <cstring>
#include <libgen.h>
#include <vector>
#include <map>

using namespace std;
using namespace LHAPDF;

extern "C" {
  void interpolation(double x, double Q, char* pdfset_path, int n_flavours, int* pdf_flavours, int pdf_number, double* values);
}

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
