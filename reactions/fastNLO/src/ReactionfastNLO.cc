// DB 08/2017
/*
   @file ReactionfastNLO.cc
   @date 2016-12-06
   @author  AddReaction.py
   Created by  AddReaction.py on 2016-12-06
*/

#include "ReactionfastNLO.h"
#include "xfitter_cpp.h"


//______________________________________________________________________________
// the class factories
extern "C" ReactionfastNLO* create() {
  return new ReactionfastNLO();
}

void ReactionfastNLO::initTerm(TermData *td)
{
  unsigned ID = td->id;

  if ( td->hasParam("Filename") ) {
    // --- Read file and instantiate fastNLO
    const std::string filename = td->getParamS("Filename");
    hf_errlog(17082501,"I: Reading fastNLO file "+filename);
    //ffnlo.insert(std::make_pair(ID,vector<fastNLOReaction>{fastNLOReaction(filename,this)} ));
    ffnlo.insert(std::make_pair(ID,vector<fastNLOReaction*>{new fastNLOReaction(filename,this)} ));
  }
  else if ( td->hasParam("Filenames") ) {
    const std::string filenames = td->getParamS("Filenames");
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
    if ( td->hasParam("Order") ) { // Local order
      hf_errlog(17090510,"W: Ignoring key 'Order' in .dat file. Only global parameter 'Order' is used.");
    }
    string order = td->getParamS("Order");  // Global order
    //hf_errlog(17090501,"I: Setting fastNLO order: "+order);
    //hf_errlog(17090501,"I: Setting fastNLO order: "+order);
    bool success=true;
    // fastNLO default is 'NLO'
    if (order=="NNLO" ) {
      if(fnlo->ContrId(fastNLO::kFixedOrder,fastNLO::kNextToNextToLeading) < 0)
        hf_errlog(19053100,"F: fastNLO. Requested order "+order+" cannot be set.");
      else
        success &= fnlo->SetContributionON(fastNLO::kFixedOrder, 2, true); // swith NNLO ON
    }
    else if (order=="NLO" ) {
      if(fnlo->ContrId(fastNLO::kFixedOrder,fastNLO::kNextToLeading) < 0)
        hf_errlog(19053100,"F: fastNLO. Requested order "+order+" cannot be set.");
      if(fnlo->ContrId(fastNLO::kFixedOrder,fastNLO::kNextToNextToLeading) >= 0)
        success &= fnlo->SetContributionON(fastNLO::kFixedOrder,2,false); // switch NNLO OFF ;
    }
    else if (order=="LO" ) {
      if(fnlo->ContrId(fastNLO::kFixedOrder,fastNLO::kLeading) < 0)
        hf_errlog(19053100,"F: fastNLO. Requested order "+order+" cannot be set.");
      if(fnlo->ContrId(fastNLO::kFixedOrder,fastNLO::kNextToLeading) >= 0)
        success &= fnlo->SetContributionON(fastNLO::kFixedOrder, 1, false); // switch NLO OFF
      if(fnlo->ContrId(fastNLO::kFixedOrder,fastNLO::kNextToNextToLeading) >= 0)
        success &= fnlo->SetContributionON(fastNLO::kFixedOrder, 2, false); // switch NNLO OFF
    }
    else    hf_errlog(17090502,"E: fastNLO. Unrecognized order: "+order);
    // --- threshold corrections
    if ( td->hasParam("ThresholdCorrection")) {
      int iThr = std::stoi(td->getParamS("ThresholdCorrection"));
      hf_errlog(17090511,"I: fastNLO. Activate threshold corrections.");
      success &= fnlo->SetContributionON(fastNLO::kThresholdCorrection, iThr, true);
    }
    if (!success)  hf_errlog(17090503,"W: fastNLO. Requested order "+order+" cannot be set.");

    // --- Set Units
    if (td->hasParam("Units") ) { // Local order
      string units = td->getParamS("Units");
      if ( units=="absolute" )
        fnlo->SetUnits(fastNLO::kAbsoluteUnits);
      else if ( units=="publication" )
        fnlo->SetUnits(fastNLO::kPublicationUnits);
      else
        hf_errlog(17090514,"E: fastNLO. Unrecognized parameter for key Units");
    }

    // --- Set scale factors
    double cmur=1, cmuf=1;
    if ( td->hasParam("ScaleFacMuR") ) {
      hf_errlog(17090504,"I: Setting fastNLO scale factor mu_R: "+td->getParamS("ScaleFacMuR"));
      cmur=std::stod(td->getParamS("ScaleFacMuR"));//GetParam("ScaleFacMuR");
    }
    if (td->hasParam("ScaleFacMuF") ) {
      hf_errlog(17090505,"I: Setting fastNLO scale factor mu_F: "+td->getParamS("ScaleFacMuF"));
      cmuf=std::stod(td->getParamS("ScaleFacMuF"));//GetParam("ScaleFacMuF");
    }
    if ( cmur!=1 || cmuf!=1 )
      fnlo->SetScaleFactorsMuRMuF(cmur,cmuf);

    // --- Set scale choice
    if ( !fnlo->GetIsFlexibleScaleTable() &&
         (td->hasParam("ScaleChoiceMuR") || td->hasParam("ScaleChoiceMuF") ) ) {
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
      if ( td->hasParam("ScaleChoiceMuR") ) { //
        hf_errlog(17090504,"I: Setting fastNLO scale choice mu_R: "+td->getParamS("ScaleChoiceMuR"));
        if ( sclmap.count(td->getParamS("ScaleChoiceMuR"))==0 )
          hf_errlog(17090522,"F: fastNLO. Scale choice for mu_R not available: "+td->getParamS("ScaleChoiceMuR"));
        else
          fnlo->SetMuRFunctionalForm(sclmap.at(td->getParamS("ScaleChoiceMuR")));
      }
      // set mu_f
      if ( td->hasParam("ScaleChoiceMuF") ) { // Local order
        hf_errlog(17090506,"I: Setting fastNLO scale choice mu_F: "+td->getParamS("ScaleChoiceMuF"));
        if ( sclmap.count(td->getParamS("ScaleChoiceMuF"))==0 )
          hf_errlog(17090523,"F: fastNLO. Scale choice for mu_F not available: "+td->getParamS("ScaleChoiceMuF"));
        else
          fnlo->SetMuFFunctionalForm(sclmap.at(td->getParamS("ScaleChoiceMuF")));
      }
    }
  }
}

//______________________________________________________________________________
// Main function to compute results at an iteration
void ReactionfastNLO::compute(TermData* td, valarray<double> &val, map<string, valarray<double> > &err)
{
  td->actualizeWrappers();
  unsigned ID = td->id;
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
}



