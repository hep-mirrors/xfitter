// DB 08/2017 
/*
   @file ReactionfastNLO.cc
   @date 2016-12-06
   @author  AddReaction.py
   Created by  AddReaction.py on 2016-12-06
*/

#include "ReactionfastNLO.h"


//______________________________________________________________________________
// the class factories
extern "C" ReactionfastNLO* create() {
  return new ReactionfastNLO();
}


//______________________________________________________________________________
// Initialize at the start of the computation
int ReactionfastNLO::initAtStart(const string &s) {
   // Nothing todo.
   return 0;
}

//______________________________________________________________________________
// Initialise for a given dataset:
void ReactionfastNLO::setDatasetParamters(int dataSetID, map<string,string> pars, map<string, double> parsDataset) {
   if ( !pars.count("Filename") )  
      hf_errlog(17082503,"F:No fastNLO file specified. Please provide 'Filename' to reaction fastNLO.");

   // --- Read file and instantiate fastNLO
   const std::string filename = pars["Filename"];
   hf_errlog(17082501,"I: Reading fastNLO file "+filename);
   fastNLOTable::ffilename = filename;
   this->ReadTable();
   this->SetFilename(filename);
   
   // --- Set order of calculation
   /*
   // for ( auto ipars : pars  ) cout<<ipars.first<<"\t\t"<<ipars.second<<endl;
   // for ( auto ipars : parsDataset   cout<<ipars.first<<"\t\t"<<ipars.second<<endl;

   int order = OrderMap( GetParamS("Order"));  // Global order
   if (pars.find("Order") != pars.end() ) { // Local order 
      int localOrder = OrderMap( pars["Order"] );
      order = localOrder>order ? order : localOrder;
   }
   */

   // --- Set scales
   // todo
   hf_errlog(17082506,"W: ReactionfastNLO. scale settings, order, etc. not yet implemented");

   // --- Set futher parameters
   // todo

}
   

//______________________________________________________________________________
// Main function to compute results at an iteration
int ReactionfastNLO::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
   // this->FillAlphasCache();
   // this->FillPDFCache();
   this->CalcCrossSection();
   std::vector<double> cs = this->GetCrossSection();
   for ( std::size_t i=0; i<cs.size(); i++)  val[i]=cs[i];
   return 0;
}


//______________________________________________________________________________
double ReactionfastNLO::EvolveAlphas(double Q ) const {
   // Implementation of Alpha_s evolution as function of the
   // factorization scale [and alphas(Mz)].
   // here we access the getter method from xfitter
   return this->alphaS(Q);
}


//______________________________________________________________________________
bool ReactionfastNLO::InitPDF(){
   //  Initalize PDF parameters if necessary. Not needed for xFitter
   return true;
}


//______________________________________________________________________________
vector<double> ReactionfastNLO::GetXFX(double xp, double muf) const {
   //  GetXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  x-proton and factorisation scale.
   return this->xfx(xp,muf);
}
