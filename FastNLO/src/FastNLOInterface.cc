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

using namespace std;


extern "C" {
  int fastnloinit_(const char *s, const int *idataset);
  int fastnlocalc_(const int *idataset, double *xsec);
}

map<int, FastNLOReader*> gFastNLO_array;

int fastnloinit_(const char *s, const int *idataset) {

  map<int, FastNLOReader*>::const_iterator FastNLOIterator = gFastNLO_array.find(*idataset);
  if(FastNLOIterator != gFastNLO_array.end( )) 
    return 1;
  
  FastNLOReader* fnloreader = NULL;
  if(string(s).compare( "H1 inclusive jet 99-00 data"))
    fnloreader = new FastNLOReader( "FastNLO/tables/fnhdesy07073incjets_10G.tab" );  
  else {
    cerr << "fastnloinit did not recognize data set <" << s << ">" <<  endl;
    return 1;
  }
    
  //fnloreader->Print();
  fnloreader->SetPDFInterface(FastNLOReader::kH1FITTER);
  fnloreader->SetAlphasEvolution(FastNLOReader::kNLOJET);
  fnloreader->SetScaleVariation(0);
  fnloreader->SetAlphasMz( 0.1178 );

  fnloreader->FillPDFCache();   // pdf is 'external'! you always have to call FillPDFCache();
  
  gFastNLO_array.insert(pair<int, FastNLOReader*>(*idataset, fnloreader) );
  return 0;
}


int fastnlocalc_(const int *idataset, double *xsec) {

  map<int, FastNLOReader*>::const_iterator FastNLOIterator = gFastNLO_array.find(*idataset);
  if(FastNLOIterator == gFastNLO_array.end( )) 
    return 1;

  FastNLOReader* fnloreader = FastNLOIterator->second;
  fnloreader->CalcCrossSection();
  
  vector < double > xs = fnloreader->GetXSection();
  vector < double > xsref = fnloreader->GetReferenceXSection();

  for ( unsigned i=0;i<xs.size();i++){
    xsec[i] = xs[i];
    //cout << "i="<<i<<"\txs="<<xs[i]<<"\tref="<<xsref[i]<<endl;
  }
  
  return 0;
}
