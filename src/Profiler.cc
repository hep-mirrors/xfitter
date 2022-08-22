#include "Profiler.h"
#include "xfitter_cpp_base.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include "xfitter_cpp.h"
#include "xfitter_steer.h"
#include "BaseEvolution.h"
#include "xfitter_cpp.h"
#include <algorithm>
#include <string.h>

extern "C" {
  void update_theory_iteration_();
  void addsystematics_(char const* name, int len);
  void store_pdfs_(char const* base, int len);
  void writefittedpoints_();
}
  
namespace xfitter
{
  std::valarray<double> Profiler::evaluatePredictions() {
    if (_getChi2) {
      int save_nsys = systema_.nsys;
      systema_.nsys = _nSourcesOrig; // reset to the original sources
      double chi2 = chi2data_theory_(2);
      systema_.nsys = save_nsys;
      std::cout << std::fixed << std::setprecision(2) << "Chi2 = " << chi2 << std::endl;
    }else{
      update_theory_iteration_();
    }

    //Return theory predictions
    return valarray<double>(c_theo_.theo, cndatapoints_.npoints);
  }

  void Profiler::storePdfFiles(int imember, int iPDF, std::string const& type) {
    string filename = _outputDir + "/pdfs_q2val_";

    char tag[20];

    if (imember>0) {
      if (type == "hessian") {
	sprintf (tag, "s%02d", (imember+1) / 2);
	filename +=  tag;
	filename +=  imember%2 == 1 ? "p_" : "m_";
      }
      else if ( type == "symmhessian") {
	sprintf (tag, "s%02d", imember );
	filename +=  tag;
	filename +=  "s_";
      }
      else if ( type == "replicas") {
	sprintf (tag, "mc%03d", imember );
	filename +=  tag;
	filename +=  "s_";
      }
    }
    writefittedpoints_();
    store_pdfs_(filename.c_str(),filename.size());
    sprintf (tag, "_%04d", imember);
    bool cp = system(((string)"cp " + _outputDir + "/fittedresults.txt "
		 + _outputDir + "/fittedresults.txt_set" + tag).c_str());

  }
  
  void Profiler::addSystematics( std::string const& name, std::valarray<double> uncertainties ) {

    // Fortran inteface
    size_t npt = uncertainties.size();
    int nsys = systema_.nsys;
    addsystematics_(name.c_str(), name.size());
    sysmeas_.n_syst_meas[nsys] = npt;
    for (int j = 0; j < npt; j++) {
      sysmeas_.syst_meas_idx[nsys][j] = j + 1;
      systema_.beta[j][nsys] =  -uncertainties[j];
      // also store asymmetric errors:
      systasym_.betaasym[j][0][nsys] = -uncertainties[j];
      systasym_.betaasym[j][1][nsys] =  uncertainties[j];
      systasym_.omega[j][nsys] = 0.0;
    }
    systasym_.lasymsyst[nsys] = false;
    systscal_.sysscalingtype[nsys] = 1;  //Apply linear scaling
    return;
  }

  void Profiler::addSystematics( std::string const& name, std::valarray<double> uncertaintiesP,   std::valarray<double> uncertaintiesM) {

    // Fortran inteface
    size_t npt = uncertaintiesP.size();
    int nsys = systema_.nsys;
    addsystematics_(name.c_str(), name.size());
    sysmeas_.n_syst_meas[nsys] = npt;
    for (int j = 0; j < npt; j++) {
      sysmeas_.syst_meas_idx[nsys][j] = j + 1;
      
      systasym_.betaasym[j][0][nsys] = - uncertaintiesP[j];
      systasym_.betaasym[j][1][nsys] = - uncertaintiesM[j];

      systema_.beta[j][nsys] = 0.5*(systasym_.betaasym[j][0][nsys] - systasym_.betaasym[j][1][nsys]);
      systasym_.omega[j][nsys] = 0.5*(systasym_.betaasym[j][0][nsys] + systasym_.betaasym[j][1][nsys]);
    };
    systasym_.lasymsyst[nsys] = true;
    systscal_.sysscalingtype[nsys] = 1;  //Apply linear scaling
    return;
  }

  void Profiler::doProfiling(){
    using namespace std;
    const YAML::Node node=XFITTER_PARS::rootNode["Profiler"];
    if(!node)return;
    if(!node.IsMap())hf_errlog(2018101201,"F: Cannot do profiling: Profiler node is not a YAML map");
    if (node["Status"]) {
      if (node["Status"].as<string>() == "Off") return;
    }
    if (node["getChi2"]) {
      _getChi2 = node["getChi2"].as<string>() == "On";
    }


    //rescaling  PDF eigenvectors
    if (node["scalePdfs"]) {
      _scalePdfs = node["scalePdfs"].as<float>();
      cerr << "[INFO] Rescaling PDF eigenvectors by factor: " <<node["scalePdfs"].as<float>()<<endl;
    }
    
    // extra info for xfitter-draw and xfitter-profile:
    if (node["enableExternalProfiler"]) {
      if (node["enableExternalProfiler"].as<string>() != "Off") {
	_storePdfs = true;
	_getChi2 = true;
      }
    }

    //    _outputDir = "output";    //default
    if(XFITTER_PARS::rootNode["OutputDirectory"]){
      _outputDir = XFITTER_PARS::rootNode["OutputDirectory"].as<string>();
    }

    int nsysloc = systema_.nsys;
    _nSourcesOrig = nsysloc;

    if (_getChi2) { // central prediction
      auto chi2tot = chi2data_theory_(1);
      auto fname = _outputDir + "/Results.txt";
      bool cp = system(( (string)"rm -f " + fname  ).c_str());
      fopen_(85, fname.c_str(), fname.size());
      double chi2Tot = chi2data_theory_(3);
      fclose_(85);
      cp = system(((string)"cp " + _outputDir + "/Results.txt " + _outputDir + "/Results_00.txt").c_str());
    }
    
    for(auto const&term:node){
      string name=term.first.as<string>();
      if(name=="Evolutions"){
        for(auto const&evol:term.second){
          profilePDF(evol.first.as<string>(),evol.second);
        }
      }
      else if (name=="Parameters") {
        for(auto const&param:term.second){
          string name = param.first.as<string>();
          if(XFITTER_PARS::gParameters.count(name)==0){
            cerr<<"[ERROR] Failed to profile parameter \""<<name<<"\": no such parameter"<<endl;
            hf_errlog(2018101202,"F: Failed to profile some parameter, see stderr");
          }
          profileParameter(name,param.second);
        }
      }
    }

    // Store theo file:
    if (node["WriteTheo"]) {
      if (node["WriteTheo"].as<string>() != "Off") {
        auto cent = evaluatePredictions();
        int ntot = systema_.nsys;
        systema_.nsys = nsysloc;
        writetheoryfiles_(ntot-nsysloc, &cent[0], node["WriteTheo"].as<string>() != "Asymmetric");
        systema_.nsys = ntot;
      }
    }

    if (_storePdfs) {
      auto fname = _outputDir + "/Results.txt";
      bool cp = system(( (string)"rm -f " + fname  ).c_str());
      systematicsflags_.resetcommonsyst = true;
      fopen_(85, fname.c_str(), fname.size());
      double chi2Tot = chi2data_theory_(3);
      fclose_(85);
      writefittedpoints_();
    }
  }

  void Profiler::profileParameter( std::string const& name, YAML::Node const& node) {
    using namespace std;
    if(!node.IsSequence()){
      cerr<<"[ERROR] Failed to profile parameter \""<<name<<"\": corresponding Profiler entry is not a sequence"<<endl;
      hf_errlog(2018101220,"F: Failed to profile some parameter, see stderr");
    }
    size_t len = node.size();
    if ( (len != 2) && (len != 3)  ) {
      hf_errlog( 2018082301,"S: Expected 2 or 3 parameters for a profiled variable. Check variable "+name);
    }

    // Evaluate varied predictions:
    std::vector< std::valarray<double> > preds;
    double *ppar = XFITTER_PARS::gParameters.at(name);

    double save = *ppar;
    for ( size_t i=0; i<len; i++) {
      *ppar = node[i].as<double>();
      updateAtConfigurationChange();
      preds.push_back( evaluatePredictions() );
    }
    *ppar = save;

    // Add to systematics list:
    if ( len == 2) {
      addSystematics(name+":T",(preds[1]-preds[0])/preds[0]);
    }
    else {
      addSystematics(name+":T",(preds[1]-preds[0])/preds[0], (preds[2]-preds[0])/preds[0]);
    }
  }

  void Profiler::profilePDF( std::string const& evolName, YAML::Node const& node) {
    // get evolution
    auto evol=get_evolution(evolName);
    // get corresponding yaml node:
    YAML::Node gNode=XFITTER_PARS::getEvolutionNode(evolName);
    YAML::Node const sets = node["sets"];
    YAML::Node const members = node["members"];
    YAML::Node const error_type_override = node["error_type_override"];
    //Sanity checks
    if ( !sets  ) {
      hf_errlog(2018082401,"S: Profiler: missing set parameters for evolution "+evolName);  // XXXXXXXXXXXXXXXX
    }
    if(!sets.IsSequence()||(error_type_override&&!error_type_override.IsSequence())){
      hf_errlog(2018082402,"S: Profiler: sets and (optional) error_type_override must be sequence");  // XXXXXXXXXXXXXXXX
    }
    if(error_type_override&&sets.size()!=error_type_override.size()){
      hf_errlog(2018082405,"S: Profiler: sets and error_type_override must be the same length");  // XXXXXXXXXXXXXXXX
    }


    size_t endi = sets.size();
    for (size_t i=0; i< endi; i++) {
      std::string pName = sets[i].as<string>();

      int central = 0;
      int first = 1;
      int last = 0;
      if ( members && members[i].IsSequence() ) {
        int msize = members[i].size();
        if(!members.IsSequence()){
          hf_errlog(2018082402,"S: Profiler: members must be sequence");  // XXXXXXXXXXXXXXXX
        }
        if(sets.size()!=members.size()){
          hf_errlog(2018082403,"S: Profiler: sets and members must be the same length");  // XXXXXXXXXXXXXXXX
        }
        if ( msize != 3) {
          hf_errlog(2018082404,"S: Profiler: sets must be sequence of length 3");  // XXXXXXXXXXXXXXXX
        }
        central = members[i][0].as<int>();
        first   = members[i][1].as<int>();
        try {
          last =  members[i][2].as<int>();
        }
        catch (...) {
          last = 0; /// auto-decodez
        }
      }
          
      // save original
      auto oSet   =Clone(gNode["set"]);
      auto oMember=Clone(gNode["member"]);

      if ( ! oSet  || ! oMember ) {
        hf_errlog(2018082410,"W: No central set or member variables for evolution : "+evolName);
      }

	if (oSet.as<string>() != pName) {
	  hf_errlog(2021011301,"W: Mismatch of the PDF set in the evolution and profiler: \033[1;31m" + oSet.as<string>() + " vs " + pName +  "\033[0m Could be ok, but beware.");
	}

        // Set central PDF and init 
        gNode["set"]   =pName;
        gNode["member"]=central;
        evol->atConfigurationChange();

	// Compatibility with fortran:
	auto c_str = pName.c_str();
	strcpy(clhapdf_.lhapdfset,c_str);
	std::fill(clhapdf_.lhapdfset+strlen(c_str),clhapdf_.lhapdfset+128,' ');

        // now we can get set properties: is it hessian asymmetric, MC or symmetric hessian
        std::string errorType = evol->getPropertyS("ErrorType");
        std::cout << "errorType: " << errorType << std::endl;

        if(error_type_override)
        {
            auto str = error_type_override[i].as<string>();
            if(str != "None")
            {
                errorType = str;
                std::cout << "errorType overwritten with: " << errorType << std::endl;
            }
        }
        if ( errorType != "symmhessian" && errorType != "hessian" && errorType != "replicas") {
           hf_errlog(2018082441,"S: Profiler Unsupported PDF error type : "+errorType);
        }

        // all predictions
        std::vector< std::valarray<double> > preds;

        preds.push_back(evaluatePredictions() );

	if (_storePdfs) {
	  storePdfFiles(0,i);
	}

        if ( last == 0) {
          // auto determine XXXXXXXXXXXXXXXX
          last = evol->getPropertyI("NumMembers")-1;
        }
        
        if ( last > evol->getPropertyI("NumMembers") ) {
          hf_errlog(2018082431,"W: Profiler: too many members requested, use max instead");
        };

        if (errorType == "hessian") {
          if ( (last - first + 1) % 2 != 0 ) {
            hf_errlog(2018082431,"S: Profiler: hessian sets need even number of members. Check your inputs");
          }
          if ( first % 2 == 0 ) {
            hf_errlog(2018082432,"S: Profiler: hessian error members should start from odd number. Check your inputs");
          }
        }
                
        // loop over all
        for (int imember = first; imember<=last; imember++) {
          gNode["member"] = imember;
          evol->atConfigurationChange();
          preds.push_back( evaluatePredictions() );
          //              for ( double th : preds[imember] ) {
          //std::cout << th << std::endl;
          //}
          //std::cout << imember << std::endl;
	  if (_storePdfs) {
	    storePdfFiles(imember,i,errorType);
	  }
        }

        // Restore original

	if ( oSet && oMember ) {
	  gNode["set"]   =oSet;
	  gNode["member"]=oMember;
	  evol->atConfigurationChange();
	}



        // Depending on error type, do nuisance parameters addition
        if ( errorType == "symmhessian" ) {
          for (int imember = first; imember<=last; imember++) {
            addSystematics("PDF_nuisance_param_"+std::to_string( ++_ipdf )+":T",(preds[imember-first+1]-preds[0])/preds[0]/_scalePdfs);
          }
        }
        else if ( errorType == "hessian") {
          for (int imember = first; imember<=last; imember += 2) {
            addSystematics("PDF_nuisance_param_"+std::to_string( ++_ipdf )+":T"
                           ,(preds[imember-first+1]-preds[0])/preds[0]/_scalePdfs
                           ,(preds[imember+1-first+1]-preds[0])/preds[0]/_scalePdfs);
          }
        }
        else if ( errorType == "replicas") {
          // construct average 
          for (size_t i=0; i<preds[0].size(); i++) {
            preds[0][i] = 0;
          }
          for ( int i=first; i<=last; i++) {
            preds[0] += preds[i-first+1];
          }
          preds[0] /= (last-first+1);
          
          // convert replicas to deviations from average:
          for ( int i=first; i<=last; i++) {
            preds[i-first+1] -= preds[0];
          }

          // convert to eigenvectors, add to list of systematics
          addReplicas(pName,preds);
        }
        else {
          hf_errlog(2018082441,"S: Profiler Unsupported PDF error type : "+errorType);
        }
        
    }
  }
  
  void Profiler::addReplicas(std::string const& pdfName,  std::vector< std::valarray<double> > const& uncertainties) {
    /// start with building the covariance matrix
    int ndata = uncertainties[0].size();
    int nrep  = uncertainties.size()-1;

    double *covar = new double [ndata*ndata];

    // potentially could be faster with BLAS, but this is not important, done once only;
    
    for ( int i=0; i<ndata; i++) {
      for ( int j=i; j<ndata; j++) {
        int id = i*ndata + j;
        covar[id] = 0;

        for ( int k=1; k<=nrep; k++) {
          covar[id] += uncertainties[k][i]*uncertainties[k][j] ;
        }

        covar[id] /= nrep;
        
        if ( i != j ) {
          int id2 = j*ndata + i;
          covar[id2] = covar[id];
        }       
      }
    }
    double *beta = new double[ndata*ndata];
    double alpha[ndata];
    int ncorr = 0;
    getnuisancefromcovar_(ndata,ndata,ndata,
                            covar,beta,0,
                            ncorr,alpha,0);
    //    std::cout << "NCorr = " << beta[0] << " " << ncorr << std::endl;
    // ready to add systematic sources

    for ( int i=0; i<ncorr; i++) {
      valarray<double> unc(ndata);
      for ( int j=0; j<ndata; j++) {
        unc[j] = beta[ndata*j + i];
        //      std::cout << unc[j] << std::endl;
      }
      addSystematics("PDF_nuisance_param_"+std::to_string(++_ipdf)+":T", unc/uncertainties[0] );
    }
    delete[] beta;
    delete[] covar;
  }
  
} //namespace xfitter

