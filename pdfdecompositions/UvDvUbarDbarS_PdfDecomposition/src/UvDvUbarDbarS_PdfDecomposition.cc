 
/*
   @file UvDvUbarDbarS_PdfDecomposition.cc
   @date 2018-07-11
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on 2018-07-11
*/

#include "UvDvUbarDbarS_PdfDecomposition.h"

// Constructor
UvDvUbarDbarS_PdfDecomposition::UvDvUbarDbarS_PdfDecomposition(const std::string& inName) : BasePdfDecomposition(inName) {
    
}

// Init at start:
void UvDvUbarDbarS_PdfDecomposition::initAtStart(const std::string & pars) {

  // Get PDFs for uv, dv, ubar, dbar and s:
  
  
  return;
}

// Compute PDF in a physical base in LHAPDF format for given x and Q
std::map<int,double>  UvDvUbarDbarS_PdfDecomposition::compute ( const double& x, const double& Q) const
{
  std::map<int,double> out;
  return out;
}



