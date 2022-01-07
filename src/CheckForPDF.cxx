#include "CheckForPDF.h"
#include "xfitter_cpp.h"
#include "LHAPDF/Paths.h"
#include <iostream>
#include <algorithm>
#include <cstring>

using namespace std;

void CheckForPDF(char const*pdfname){
  bool found= false;
  string spdfname=string(pdfname);
  spdfname.erase(std::remove_if(spdfname.begin(), spdfname.end(), ::isspace),spdfname.end());
  for (const string& s : LHAPDF::availablePDFSets()){
    if (s==spdfname) found=true;
	    }
  if (!found){
    cout << "PDF: " << spdfname << " not found"<< endl;
    cout << "List of available PDFs in this installation:" << endl;
    for (const string& s : LHAPDF::availablePDFSets())
      cout << " " << s << endl;

    cout << "You can try installing your chosen pdf with: ./tools/download-lhapdf.sh " << spdfname  << endl;

    exit(0);
  } else {
    string msg = "I: PDF: " + spdfname + " found in this installation";
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
