#include <stdlib.h>
#include <iostream>
#include <TSystem.h>
#include <TString.h>
#include <getopt.h>
#include <TError.h>

using std::cout;
using std::cerr;
using std::endl;

#include "Painter.h"

#include "PdfTable.h"


int main(int argc, char **argv) {

  TString OutputPath("output/");
  TString OutputPathRef("output/");
  
  gErrorIgnoreLevel=1001;

//  int c; 
  bool DrawBands = false;
  bool DrawExp = false;
  bool DrawModel = false;
  bool DrawParam = false;
  bool DrawBase = false;

  int ioptband = -1;
  int ioptexp = -1;
  int ioptmodel = -1;
  int ioptparam = -1;

  int ibase = 0;
  int ibase2 = 0;

  if ( argc == 1 ) {
    printf("program usage:\n DrawResults [--bands] [dir1] [dir2]\n      OR \n DrawResults [--exp] [exp_dir] [--model] [model_dir] [--param] [param_dir] [basedir]\n");
    return -1;
  }

  if ( argc > 2 ) {
    for (int iar=1; iar<argc; iar++) {      
      if ( ! TString(argv[iar]).CompareTo("--bands") ) {
	if (DrawBands) {
	  printf("Do not use --bands option more then once.\n");
	  return -1;
	}
	DrawBands = true;
	ioptband = iar;
      }
      if ( ! TString(argv[iar]).CompareTo("--exp") ) {
	if (DrawExp) {
	  printf("Do not use --exp option more then once.\n");
	  return -1;
	}
	DrawExp = true;
	ioptexp = iar;
      }
      if ( ! TString(argv[iar]).CompareTo("--model") ) {
	if (DrawModel) {
	  printf("Do not use --model option more then once.\n");
	  return -1;
	}
	DrawModel = true;
	ioptmodel = iar;
      }
      if ( ! TString(argv[iar]).CompareTo("--param") ) {
	if (DrawParam) {
	  printf("Do not use --param option more then once.\n");
	  return -1;
	}
	DrawParam = true;
	ioptparam = iar;
      }
    }
  }

  if (DrawBands && (DrawExp||DrawModel||DrawParam)) {
    printf("Do not use --bands option together with --exp, --model or --param option.\n");
    return -1;
  }

  if ( argc > 1 ) {
    for (int iar=1; iar<argc; iar++) {
      if (iar!=ioptband && iar!=ioptband+1 && iar!=ioptexp && iar!=ioptexp+1 && iar!=ioptmodel && iar!=ioptmodel+1 && iar!=ioptparam && iar!=ioptparam+1 ) {
	if (!DrawBase) {
	  ibase = iar;
	  DrawBase = true;
	} else {
	  if ( ibase2==0 && !(DrawBands||DrawExp||DrawModel||DrawParam) ) {
	    ibase2 = iar;
	  } else {
	    printf("program usage:\n DrawResults [--bands] [dir1] [dir2]\n      OR \n DrawResults [--exp] [exp_dir] [--model] [model_dir] [--param] [param_dir] [basedir]\n");
	    return -1;
	  }
	}
      }
    }
  }

  if ( !(DrawExp||DrawModel||DrawParam) ) {
    
    if (DrawBase) {
      OutputPathRef.Form(argv[ibase]);
    } else if (ioptband > 0) {
      OutputPathRef.Form(argv[ioptband+1]);
    } else {
      printf("Error: no base or bands directory.\n");
      return -1;
    }
    
    if (ibase2!=0) {
      OutputPath.Form(argv[ibase2]);
    } else if (ioptband > 0) {
      OutputPath.Form(argv[ioptband+1]);
    } else if (DrawBase) {
      OutputPath.Form(argv[ibase]);
    }
    
    Painter* painter = new Painter(DrawBands);
    painter->SetPath(OutputPath);
    painter->SetPathRef(OutputPathRef);
    painter->Draw();
    delete painter;
    
  } else {
    
    Painter* painter = new Painter(false, DrawBase, DrawExp, DrawModel, DrawParam);
    if (DrawBase) painter->SetPathBase(argv[ibase]);
    if (DrawExp) painter->SetPathExp(argv[ioptexp+1]);
    if (DrawModel) painter->SetPathModel(argv[ioptmodel+1]);
    if (DrawParam) painter->SetPathParam(argv[ioptparam+1]);
    painter->Draw();
    delete painter;
    
  }
  //call here new drawing functions, need to change options parsing
  return 0;

}
