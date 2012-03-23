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
typedef vector<bool> BoolArray; 


extern "C" {
  int fastnloinit_(const char *s, const int *idataset, const char *thfile, bool *PublicationUnits , double* murdef, double *mufdef);
  int fastnlocalc_(const int *idataset, double *xsec);
  int getalf_( double* alfs, double* r2 );
  int fastnlopointskip_(const int *idataset, int *point, int *npoints);
  int hf_errlog_(const int* ID, const char* TEXT, long length);
  int hf_stop_();
}

map<int, FastNLOReader*> gFastNLO_array;
map<int, BoolArray*>     gUsedPoints_array;
int CreateUsedPointsArray(int idataset, int npoints);

int fastnloinit_(const char *s, const int *idataset, const char *thfile, bool *PublicationUnits , double* murdef, double *mufdef ) {

  
  //cout << "FastNLOINterface::fastnloinit_ idataset = "<<*idataset<< ", PublicationUnits "<< *PublicationUnits << endl;

   map<int, FastNLOReader*>::const_iterator FastNLOIterator = gFastNLO_array.find(*idataset);
   if(FastNLOIterator != gFastNLO_array.end( )) {
     int id = 12032301;
     char* text = "I: Double initialization of the same fastnlo data set!";
     hf_errlog_(&id, text, (long)strlen(text));
     //hf_stop_();
     return 1;
   }
  
   FastNLOReader* fnloreader = NULL;
   fnloreader = new FastNLOReader( thfile );  
   
   if(*PublicationUnits)
     fnloreader->SetUnits(FastNLOReader::kPublicationUnits);
   else 
     fnloreader->SetUnits(FastNLOReader::kAbsoluteUnits);

   if(*murdef>=0.)
     fnloreader->SetMuRFunctionalForm((FastNLOReader::EScaleFunctionalForm) ((int) (*murdef)) , false);
   if(*mufdef>=0.)
     fnloreader->SetMuFFunctionalForm((FastNLOReader::EScaleFunctionalForm) ((int) (*mufdef)) , false);

   fnloreader->SetPDFInterface(FastNLOReader::kH1FITTER);
   fnloreader->SetAlphasMz( 0.1180 ); // just for initalisation
   fnloreader->SetAlphasEvolution(FastNLOReader::kQCDNUMInternal); // fully consistend alpha_s evolution has to be implemented.

   // looking for scale factor = 1!
   int nscale = fnloreader->GetNScaleVariations();
   vector<double> scalefactors = fnloreader->GetScaleFactors();
   for (int iscale = 0;iscale < nscale ;iscale++) {
      if ( fabs(scalefactors[iscale] - 1.) < 0.0001 ){
	fnloreader->SetScaleVariation(iscale, false);
      }
   }

   // switching non-pert corr off - Done by default now
   //fnloreader->SetContributionON(FastNLOReader::kNonPerturbativeCorrection,0,false);
   //fnloreader->SetContributionON(FastNLOReader::kNonPerturbativeCorrection,1,false);

   // no threshold corrections
   //fnloreader->SetContributionON(FastNLOReader::kThresholdCorrection,0,false);


   //fnloreader->FillAlphasCache();
   //fnloreader->FillPDFCache();			// pdf is 'external'! you always have to call FillPDFCache();
   //fnloreader->CalcCrossSection();
   //fnloreader->PrintCrossSectionsLikeFreader();

   gFastNLO_array.insert(pair<int, FastNLOReader*>(*idataset, fnloreader) );
   return 0;
}


int fastnlocalc_(const int *idataset, double *xsec) {

  //cout << "FastNLOINterface::fastnlocalc_ idataset = " <<*idataset<<endl;

   map<int, FastNLOReader*>::const_iterator FastNLOIterator = gFastNLO_array.find(*idataset);
   map<int, BoolArray*>::const_iterator UsedPointsIterator = gUsedPoints_array.find(*idataset);
   if(FastNLOIterator == gFastNLO_array.end( )) {
     int id = 12032302;
     char text[256];
     sprintf(text, "S: Can not find FastnloReader for DataSet: %d", *idataset);
     hf_errlog_(&id, text, (long)strlen(text)); // this terminates the program by default
   }
   
   FastNLOReader* fnloreader = FastNLOIterator->second;
   
   if(UsedPointsIterator == gUsedPoints_array.end( )) 
     CreateUsedPointsArray(*idataset, fnloreader->GetNObsBins());
   UsedPointsIterator = gUsedPoints_array.find(*idataset);
   
   if(UsedPointsIterator == gUsedPoints_array.end( )) {
     int id = 12032303;
     char text[256];
     sprintf(text, "S: Can not find proper UsedPointsIterator for DataSet: %d", *idataset);
     hf_errlog_(&id, text, (long)strlen(text)); // this terminates the program by default
   }

   BoolArray*     usedpoints = UsedPointsIterator->second;

   double alfsMz= 0;
   double Mz2= 0.;
   getalf_(&alfsMz,&Mz2);
   fnloreader->SetAlphasMz( alfsMz );
   fnloreader->FillAlphasCache();
   fnloreader->FillPDFCache();			// pdf is 'external'! you always have to call FillPDFCache();

   fnloreader->CalcCrossSection();
   //fnloreader->PrintCrossSectionsLikeFreader();

   vector < double > xs = fnloreader->GetCrossSection();
 
   int outputidx = 0;
   for ( unsigned i=0;i<xs.size();i++){
     if(usedpoints->at(i)) {
       xsec[outputidx] = xs[i];
       outputidx++;
     }
   }
 
   return 0;
}


int fastnlopointskip_(const int *idataset, int *point, int *npoints) {
  map<int, BoolArray*>::const_iterator UsedPointsIterator = gUsedPoints_array.find(*idataset);
  if(UsedPointsIterator == gUsedPoints_array.end( )) 
    CreateUsedPointsArray(*idataset, *npoints);

  UsedPointsIterator = gUsedPoints_array.find(*idataset);
  if(UsedPointsIterator == gUsedPoints_array.end( )) {
    int id = 12032304;
    char text[256];
    sprintf(text, "S: fastnlopointskip: Can not find proper UsedPointsIterator for DataSet: %d",*idataset);
    hf_errlog_(&id, text, (long)strlen(text)); // this terminates the program by default
  }
  
  BoolArray*     usedpoints = UsedPointsIterator->second;
  usedpoints->at(*point-1) = false;
}

int CreateUsedPointsArray(int idataset, int npoints) {
  cout << "creating new table..."<<endl;
  BoolArray* usedpoints = new BoolArray;
  for (int i=0; i<npoints; i++)
    usedpoints->push_back(true);
  gUsedPoints_array.insert(pair<int, BoolArray*>(idataset, usedpoints) );
}
