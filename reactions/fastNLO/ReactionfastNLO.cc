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

struct DatasetData{
  vector<fastNLOReaction*>ffnlo;
};

void ReactionfastNLO::initTerm(TermData *td)
{
  DatasetData* data = new DatasetData;
  td->reactionData = (void*)data;
  vector<fastNLOReaction*>& ffnlo=data->ffnlo;

  BaseEvolution* pdf = td->getPDF();
  if(pdf!=td->getPDF(1)){
    hf_errlog(2020022500,"F: fastNLO reaction does not support two different PDFs");
  }

  if ( td->hasParam("Filename") ) {
    // --- Read file and instantiate fastNLO
    const std::string filename = td->getParamS("Filename");
    hf_errlog(17082501,"I: Reading fastNLO file "+filename);
    ffnlo.push_back(new fastNLOReaction(filename,pdf));
  }
  else if ( td->hasParam("Filenames") ) {
    const std::string filenames = td->getParamS("Filenames");
    std::stringstream ss(filenames);
    std::string filename;
    while (std::getline(ss, filename, ',')) {
      hf_errlog(17082501,"I: Reading fastNLO file "+filename);
      ffnlo.push_back(new fastNLOReaction(filename,pdf));
    }
  }
  else {
    hf_errlog(17082503,"F: No fastNLO file specified. Please provide 'Filename' or 'Filenames' to reaction fastNLO.");
    return;
  }

  for ( fastNLOReaction* fnlo : ffnlo ) {
    // --- Set order of calculation
    string order = td->getParamS("Order");
    bool success=true;
    // fastNLO default is 'NLO'
    if (order=="NNLO" ) {
      if(fnlo->ContrId(fastNLO::kFixedOrder,fastNLO::kNextToNextToLeading) < 0)
        hf_errlog(19073100,"W: fastNLO. Requested order "+order+" cannot be set. Using NLO only!");
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
      int iThr = td->getParamI("ThresholdCorrection");
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
      cmur=*td->getParamD("ScaleFacMuR");
    }
    if (td->hasParam("ScaleFacMuF") ) {
      hf_errlog(17090505,"I: Setting fastNLO scale factor mu_F: "+td->getParamS("ScaleFacMuF"));
      cmuf=*td->getParamD("ScaleFacMuF");
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

void ReactionfastNLO::freeTerm(TermData* td)
{
  //Free allocated resources
  DatasetData* data = (DatasetData*)td->reactionData;
  for(const fastNLOReaction* fastnlo:data->ffnlo){
    delete fastnlo;
  }
  delete data;
  td->reactionData = nullptr;
}

//______________________________________________________________________________
// Main function to compute results at an iteration
void ReactionfastNLO::compute(TermData* td, valarray<double> &val, map<string, valarray<double> > &err)
{
  size_t N=val.size();
  vector<fastNLOReaction*>& ffnlo=((DatasetData*)(td->reactionData))->ffnlo;
  vector<double> valfnlo;
  valfnlo.reserve(N);
  for ( fastNLOReaction* fnlo : ffnlo ) {
    fnlo->CalcCrossSection();
    vector<double> cs = fnlo->GetCrossSection();
    //append to valfnlo
    valfnlo.insert(valfnlo.end(),cs.begin(),cs.end());
  }

  bool symmetrise=td->hasParam("Symmetrise") && 2*N==valfnlo.size();
  if ( !(N == valfnlo.size()) && !symmetrise)
    hf_errlog(17112501,"F: Size of fastNLO cross section array does not match data.");
  
  if (symmetrise){
    for ( std::size_t i=0; i<N; i++) val[i]=valfnlo[i+N]+valfnlo[N-1-i];
  }else{
    for ( std::size_t i=0; i<N; i++) val[i]=valfnlo[i];
  }
}
