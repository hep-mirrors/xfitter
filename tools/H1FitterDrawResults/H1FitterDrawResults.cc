#include <stdlib.h>
#include <iostream>
#include <TSystem.h>
#include <TString.h>
#include <getopt.h>


using std::cout;
using std::cerr;
using std::endl;

#include "H1FitterPainter.h"

#include "PdfTable.h"


int main(int argc, char **argv) {

  TString OutputPath("output/");
  TString OutputPathRef("output/");
  

  int c; 
  Int_t nBands = 0;

  while (1)
    {
      static struct option long_options[] =
	{
	  {"bands",  required_argument, 0, 'b'},
	};
      /* getopt_long stores the option index here. */
      int option_index = 0;
      
      c = getopt_long (argc, argv, "bands:",
                            long_options, &option_index);
      /* Detect the end of the options. */
      if (c == -1)
	break;
      switch (c) 
	{
	case 'b':
	  nBands = atoi(optarg); 
	  break;
	case '?':
	  printf("program usage: DrawResults --bands=n [dir1] [dir2]\n");
	  exit(1);
	   
	}

    }
  
  printf ("NBands=%i\n",nBands);

  if ( optind < argc) {

    printf ("%i %i %s\n",optind, argc, argv[optind]);
    if(argc-optind == 1 ) {
      OutputPath.Form(argv[optind]);
      OutputPathRef.Form(argv[optind]);
    }
    if(argc-optind == 2) {
      OutputPath.Form(argv[optind]);
      OutputPathRef.Form(argv[optind+1]);
    }
  }


  H1FitterPainter* painter = new H1FitterPainter(nBands);
  painter->SetPath(OutputPath);
  painter->SetPathRef(OutputPathRef);
  painter->Draw();
  delete painter;
  return 0;

}
