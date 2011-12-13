// Author: Krzysztof Nowak
// DESY, 01/08/2011

//  Version 0.1, 

////////////////////////////////////////////////////////////////////////
//
//   FastNLOInterface
//                                                                      //
//  The interface through which fortran based h1fitter interacts with   //
//  c++ version of FastNLOReader.                                       // 
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <map>
#include <FastNLOReader.h>
#include <cmath>

using namespace std;


extern "C" {
   int fastnloinit_(const char *s, const int *idataset, const char *thfile );
   int fastnlocalc_(const int *idataset, double *xsec);
   int getalf_( double* alfs, double* r2 );
}

map<int, FastNLOReader*> gFastNLO_array;

int fastnloinit_(const char *s, const int *idataset, const char *thfile  ) {

   map<int, FastNLOReader*>::const_iterator FastNLOIterator = gFastNLO_array.find(*idataset);
   if(FastNLOIterator != gFastNLO_array.end( )) 
      return 1;
  
   FastNLOReader* fnloreader = NULL;
   if(string(s).compare( "H1 inclusive jet 99-00 data"))
      fnloreader = new FastNLOReader( thfile );  
   else {
      cerr << "fastnloinit did not recognize data set <" << s << ">" <<  endl;
      return 1;
   }
   
   fnloreader->SetPDFInterface(FastNLOReader::kH1FITTER);
   //fnloreader->SetAlphasEvolution(FastNLOReader::kQCDNUMInternal); // fully consistend alpha_s evolution has to be implemented.
   //fnloreader->SetScaleVariation(0);

   gFastNLO_array.insert(pair<int, FastNLOReader*>(*idataset, fnloreader) );
   return 0;
}


int fastnlocalc_(const int *idataset, double *xsec) {


   map<int, FastNLOReader*>::const_iterator FastNLOIterator = gFastNLO_array.find(*idataset);
   if(FastNLOIterator == gFastNLO_array.end( )) 
      return 1;

   FastNLOReader* fnloreader = FastNLOIterator->second;

   double alfsMz= 0;
   double Mz2= 0.;
   getalf_(&alfsMz,&Mz2);
   fnloreader->SetAlphasMz( alfsMz );
   fnloreader->FillPDFCache();			// pdf is 'external'! you always have to call FillPDFCache();

   fnloreader->CalcCrossSection();

   vector < double > xs = fnloreader->GetXSection();
 
   for ( unsigned i=0;i<xs.size();i++){
      xsec[i] = xs[i];
   }
  
   return 0;
}
