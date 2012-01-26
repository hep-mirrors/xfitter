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

   //cout << "FastNLOINterface::fastnloinit_ idataset = "<<*idataset<<endl;

   map<int, FastNLOReader*>::const_iterator FastNLOIterator = gFastNLO_array.find(*idataset);
   if(FastNLOIterator != gFastNLO_array.end( )) 
      return 1;
  
   FastNLOReader* fnloreader = NULL;
   fnloreader = new FastNLOReader( thfile );  
   //    if(string(s).compare( "H1 inclusive jet 99-00 data"))
   //       fnloreader = new FastNLOReader( thfile );  
   //    else {
   //       cerr << "fastnloinit did not recognize data set <" << s << ">" <<  endl;
   //       return 1;
   //    }
   
   
   fnloreader->SetUnits(FastNLOReader::kAbsoluteUnits);
   if ( string(s).find("D0 pp jets") == 0 ){;
      cout << "D0 -> Setting to publication units."<<endl;
      fnloreader->SetUnits(FastNLOReader::kPublicationUnits);
   }
   else if ( string(s).find("CMS inclusive jets") == 0 ){;
      cout << "D0 -> Setting to publication units."<<endl;
      fnloreader->SetUnits(FastNLOReader::kPublicationUnits);
   }
   else{
      cout << "could not recognize D0 jets. Using absolute units."<<endl;
   }

   fnloreader->SetPDFInterface(FastNLOReader::kH1FITTER);
   fnloreader->SetAlphasMz( 0.1180 ); // just for initalisation
   fnloreader->SetAlphasEvolution(FastNLOReader::kQCDNUMInternal); // fully consistend alpha_s evolution has to be implemented.

   // looking for scale factor = 1!
   int nscale = fnloreader->GetNScaleVariations();
   vector<double> scalefactors = fnloreader->GetScaleFactors();
   for (int iscale = 0;iscale < nscale ;iscale++) {
      if ( fabs(scalefactors[iscale] - 1.) < 0.0001 ){
	 fnloreader->SetScaleVariation(iscale);
      }
   }

   // switching non-pert corr off
   fnloreader->SetContributionON(FastNLOReader::kNonPerturbativeCorrection,0,false);
   fnloreader->SetContributionON(FastNLOReader::kNonPerturbativeCorrection,1,false);

   // no threshold corrections
   //fnloreader->SetContributionON(FastNLOReader::kThresholdCorrection,0,false);


   fnloreader->FillAlphasCache();
   fnloreader->FillPDFCache();			// pdf is 'external'! you always have to call FillPDFCache();
   fnloreader->CalcCrossSection();
   //fnloreader->PrintCrossSectionsLikeFreader();



   gFastNLO_array.insert(pair<int, FastNLOReader*>(*idataset, fnloreader) );
   return 0;
}


int fastnlocalc_(const int *idataset, double *xsec) {

   cout << "FastNLOINterface::fastnlocalc_ idataset = " <<*idataset<<endl;

   // call QCDNUM::_evaluate here!
   
   map<int, FastNLOReader*>::const_iterator FastNLOIterator = gFastNLO_array.find(*idataset);
   if(FastNLOIterator == gFastNLO_array.end( )) 
      return 1;

   FastNLOReader* fnloreader = FastNLOIterator->second;

   double alfsMz= 0;
   double Mz2= 0.;
   getalf_(&alfsMz,&Mz2);
   fnloreader->SetAlphasMz( alfsMz );
   fnloreader->FillAlphasCache();
   fnloreader->FillPDFCache();			// pdf is 'external'! you always have to call FillPDFCache();

   fnloreader->CalcCrossSection();
   //fnloreader->PrintCrossSectionsLikeFreader();

   vector < double > xs = fnloreader->GetCrossSection();
 
   for ( unsigned i=0;i<xs.size();i++){
      xsec[i] = xs[i];
   }
  
   return 0;
}
