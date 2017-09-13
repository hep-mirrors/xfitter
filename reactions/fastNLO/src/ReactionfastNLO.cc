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
   if ( pars.count("Order") ) { // Local order 
      hf_errlog(17090501,"I: Setting fastNLO order: "+pars["Order"]);
      bool success=true;
      // fastNLO default is 'NLO'
      if (pars["Order"]=="NNLO" ) success&=this->SetContributionON(fastNLO::kFixedOrder, 2, true);
      else if (pars["Order"]=="NLO" ) {;}
      else if (pars["Order"]=="LO" )  success&=this->SetContributionON(fastNLO::kFixedOrder, 1, false);
      else if (pars["Order"]=="Thr" ) success&=this->SetContributionON(fastNLO::kThresholdCorrection, 1, true);
      else {
	 hf_errlog(17090502,"E: fastNLO. Unrecognized order: "+pars["Order"]);
      }
      if (!success) {
	 hf_errlog(17090503,"W: fastNLO. Requested order cannot be set.");
      }
   }

   // --- Set scale factors
   double cmur=1, cmuf=1;
   if ( pars.count("ScaleFacMuR") ) { // Local order 
      hf_errlog(17090504,"I: Setting fastNLO scale factor mu_R: "+pars["ScaleFacMuR"]);
      cmur=GetParam("ScaleFacMuR");
   }
   if ( pars.count("ScaleFacMuF") ) { // Local order 
      hf_errlog(17090505,"I: Setting fastNLO scale factor mu_F: "+pars["ScaleFacMuF"]);
      cmuf=GetParam("ScaleFacMuF");
   }
   if ( cmur!=1 || cmuf!=1 ) 
      this->SetScaleFactorsMuRMuF(cmur,cmuf);

   // --- Set scale choice
   if ( !this->GetIsFlexibleScaleTable() && 
	(pars.count("ScaleChoiceMuR") || pars.count("ScaleChoiceMuF") ) ) {
      hf_errlog(17090508,"W: fastNLO. Scale choice requested, but this is not a flexible scale table.");
   }
   else {
      const std::map<std::string,fastNLO::EScaleFunctionalForm> sclmap{
	 {"kScale1"              ,fastNLO::kScale1},
	 {"kScale2"              ,fastNLO::kScale2},              
	 {"kQuadraticSum"        ,fastNLO::kQuadraticSum},        
	 {"kQuadraticMean"       ,fastNLO::kQuadraticMean},       
	 {"kQuadraticSumOver4"   ,fastNLO::kQuadraticSumOver4},   
	 {"kLinearMean"          ,fastNLO::kLinearMean},          
	 {"kLinearSum"           ,fastNLO::kLinearSum},           
	 {"kScaleMax"            ,fastNLO::kScaleMax},            
	 {"kScaleMin"            ,fastNLO::kScaleMin},            
	 {"kProd"                ,fastNLO::kProd},                
	    // {"kS2plusS1half"        ,fastNLO::kS2plusS1half},        
	    // {"kPow4Sum"             ,fastNLO::kPow4Sum},             
	    // {"kWgtAvg"              ,fastNLO::kWgtAvg},              
	    // {"kS2plusS1fourth"      ,fastNLO::kS2plusS1fourth},      
	    // {"kExpProd2"            ,fastNLO::kExpProd2},            
	    // {"kExtern"              ,fastNLO::kExtern},              
	    // {"kConst"               ,fastNLO::kConst},
	    };              
      // set mu_r
      if ( pars.count("ScaleChoiceMuR") ) { // 
	 hf_errlog(17090504,"I: Setting fastNLO scale choice mu_R: "+pars["ScaleChoiceMuR"]);
	 if ( sclmap.count(GetParamS("ScaleFacMuR"))==0 )
	      hf_errlog(17090504,"W: fastNLO. Scale choice for mu_R not available: "+GetParamS("ScaleFacMuR"));
	 else
	    this->SetMuRFunctionalForm(sclmap.at(GetParamS("ScaleFacMuR")));
      }
      // set mu_f
      if ( pars.count("ScaleChoiceMuF") ) { // Local order 
	 hf_errlog(17090506,"I: Setting fastNLO scale choice mu_F: "+pars["ScaleChoiceMuF"]);
	 if ( sclmap.count(GetParamS("ScaleFacMuF"))==0 )
	    hf_errlog(17090507,"W: fastNLO. Scale choice for mu_F not available: "+GetParamS("ScaleFacMuF"));
	 else 
	    this->SetMuFFunctionalForm(sclmap.at(GetParamS("ScaleFacMuF")));
      }
   }
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



