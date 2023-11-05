#include "CheckForPDF.h"
#include "xfitter_cpp.h"
#include "LHAPDF/Paths.h"
#include <iostream>
#include <algorithm>
#include <string>

void CheckForPDF(char const* pdfname) {
  bool found = false;
  std::string spdfname(pdfname);
  spdfname.erase(std::remove_if(spdfname.begin(), spdfname.end(), ::isspace), spdfname.end());

  for (const std::string& s : LHAPDF::availablePDFSets()) {
    if (s == spdfname) {
      found = true;
      break;
	  }
  }

  if (!found){
    std::cout << "PDF: " << spdfname << " not found" << std::endl;
    std::cout << "List of available PDFs in this installation:" << std::endl;
    for (const string& s : LHAPDF::availablePDFSets()) {
      std::cout << " " << s << std::endl;
    }

    std::cout << "You can try installing your chosen pdf with: ./tools/download-lhapdf.sh " << spdfname  << std::endl;

    std::exit(0);

  } else {
    std::string msg = "I: PDF: " + spdfname + " found in this installation";
    hf_errlog_(19092501, msg.c_str(), msg.size());
  }
}


// Fortran wrapper

extern "C" {
  void checkforpdf_(char *pdfname, long int length){
    char tmp[length];
    memcpy(tmp,pdfname,length);
    tmp[length-1] = '\0';
    CheckForPDF(tmp);
  }
}
