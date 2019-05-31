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
void ReactionfastNLO::setDatasetParameters(int ID, map<string,string> pars, map<string, double> parsDataset) {
   // ID=dataSetID
   if ( pars.count("Filename") ) {
      // --- Read file and instantiate fastNLO
      const std::string filename = pars["Filename"];
      hf_errlog(17082501,"I: Reading fastNLO file "+filename);
      //ffnlo.insert(std::make_pair(ID,vector<fastNLOReaction>{fastNLOReaction(filename,this)} ));
      ffnlo.insert(std::make_pair(ID,vector<fastNLOReaction*>{new fastNLOReaction(filename,this)} ));
   }
   else if ( pars.count("Filenames") ) {
      const std::string filenames = pars["Filenames"];
      std::stringstream ss(filenames);
      std::string filename;
      while (std::getline(ss, filename, ',')) {
	 hf_errlog(17082501,"I: Reading fastNLO file "+filename);
	 ffnlo[ID].push_back(new fastNLOReaction(filename,this));
      }
   }
   else {
      hf_errlog(17082503,"F:No fastNLO file specified. Please provide 'Filename' or 'Filenames' to reaction fastNLO.");
      return;
   }

   // fastNLOTable::ffilename = filename;
   // this->ReadTable();
   // this->SetFilename(filename);
   
   for ( fastNLOReaction* fnlo : ffnlo[ID] ) {
      // --- Set order of calculation 
      if ( pars.count("Order") ) { // Local order 
	 hf_errlog(17090510,"W: Ignoring key 'Order' in .dat file. Only global parameter 'Order' is used.");
      }
      string order = GetParamS("Order");  // Global order
      //hf_errlog(17090501,"I: Setting fastNLO order: "+order);
      //hf_errlog(17090501,"I: Setting fastNLO order: "+order);
      bool success=true;
      // fastNLO default is 'NLO'
      if (order=="NNLO" ) success &= fnlo->SetContributionON(fastNLO::kFixedOrder, 2, true); // swith NNLO ON
      else if (order=="NLO" ) success &= fnlo->SetContributionON(fastNLO::kFixedOrder,2,false); // switch 'off' NNLO ;
      else if (order=="LO" ) {
        success &= fnlo->SetContributionON(fastNLO::kFixedOrder, 1, false); // switch NLO OFF
        success &= fnlo->SetContributionON(fastNLO::kFixedOrder, 2, false); // switch NNLO OFF
      }
      else    hf_errlog(17090502,"E: fastNLO. Unrecognized order: "+order);
      // --- threshold corrections
      if ( pars.count("ThresholdCorrection")) {
	 int iThr = std::stoi(pars["ThresholdCorrection"]);
	 hf_errlog(17090511,"I: fastNLO. Activate threshold corrections.");
	 success &= fnlo->SetContributionON(fastNLO::kThresholdCorrection, iThr, true);
      }
      //if (!success)  hf_errlog(17090503,"W: fastNLO. Requested order "+order+" cannot be set.");
   
      // --- Set Units
      if ( pars.count("Units") ) { // Local order 
	 string units = pars["Units"] ;
	 if ( units=="absolute" ) 
	    fnlo->SetUnits(fastNLO::kAbsoluteUnits);
	 else if ( units=="publication" ) 
	    fnlo->SetUnits(fastNLO::kPublicationUnits);
	 else
	    hf_errlog(17090514,"E: fastNLO. Unrecognized parameter for key Units");
      }

      // --- Set scale factors
      double cmur=1, cmuf=1;
      if ( pars.count("ScaleFacMuR") ) {
	 hf_errlog(17090504,"I: Setting fastNLO scale factor mu_R: "+pars["ScaleFacMuR"]);
	 cmur=std::stod(pars["ScaleFacMuR"]);//GetParam("ScaleFacMuR");
      }
      if ( pars.count("ScaleFacMuF") ) {
	 hf_errlog(17090505,"I: Setting fastNLO scale factor mu_F: "+pars["ScaleFacMuF"]);
	 cmuf=std::stod(pars["ScaleFacMuF"]);//GetParam("ScaleFacMuF");
      }
      if ( cmur!=1 || cmuf!=1 ) 
	 fnlo->SetScaleFactorsMuRMuF(cmur,cmuf);

      // --- Set scale choice
      if ( !fnlo->GetIsFlexibleScaleTable() && 
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
	    if ( sclmap.count(pars["ScaleChoiceMuR"])==0 )
	       hf_errlog(17090522,"F: fastNLO. Scale choice for mu_R not available: "+pars["ScaleChoiceMuR"]);
	    else 
	       fnlo->SetMuRFunctionalForm(sclmap.at(pars["ScaleChoiceMuR"]));
	 }
	 // set mu_f
	 if ( pars.count("ScaleChoiceMuF") ) { // Local order 
	    hf_errlog(17090506,"I: Setting fastNLO scale choice mu_F: "+pars["ScaleChoiceMuF"]);
	    if ( sclmap.count(pars["ScaleChoiceMuF"])==0 ) 
	       hf_errlog(17090523,"F: fastNLO. Scale choice for mu_F not available: "+pars["ScaleChoiceMuF"]);
	    else 
	       fnlo->SetMuFFunctionalForm(sclmap.at(pars["ScaleChoiceMuF"]));
	 }
      }
   }
}

//______________________________________________________________________________
// Main function to compute results at an iteration
int ReactionfastNLO::compute(int ID, valarray<double> &val, map<string, valarray<double> > &err)
{
   // ffnlo[ID][0]->FillAlphasCache();
   // ffnlo[ID][0]->FillPDFCache();
   vector<double> valfnlo;
   valfnlo.reserve(val.size());
   for ( fastNLOReaction* fnlo : ffnlo[ID] ) {
      fnlo->CalcCrossSection();
      //ffnlo[ID][0]->PrintCrossSections();
      std::vector<double> cs = fnlo->GetCrossSection();
      for ( std::size_t i=0; i<cs.size(); i++) valfnlo.push_back(cs[i]);
   }

   if ( val.size()!=valfnlo.size() )
      hf_errlog(17112501,"F: Size of fastNLO cross section array does not match data.");
   for ( std::size_t i=0; i<val.size(); i++) val[i]=valfnlo[i];
   return 0;
}



