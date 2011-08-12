#include <stdlib.h>
#include <iostream>
#include <TSystem.h>
#include <TString.h>
using std::cout;
using std::cerr;
using std::endl;

#include "H1FitterPainter.h"

int main(int argc, char **argv) {

  TString OutputPath("output/");
  TString OutputPathRef("output/");
  
  if(argc == 2) {
    OutputPath.Form(argv[1]);
    OutputPathRef.Form(argv[1]);
  }
  if(argc == 3) {
    OutputPath.Form(argv[1]);
    OutputPathRef.Form(argv[2]);
  }


  H1FitterPainter* painter = new H1FitterPainter;
  painter->SetPath(OutputPath);
  painter->SetPathRef(OutputPathRef);
  painter->Draw();
  delete painter;
  return 0;

}
