#include "Profiler.h"
#include "xfitter_cpp_base.h"
#include <iostream>
#include <vector>
#include "xfitter_cpp.h"
#include "xfitter_steer.h"
#include "BaseEvolution.h"

extern "C" {
  void update_theory_iteration_();
  void addsystematics_(char const* name, int len);
}
  
namespace xfitter
{
  std::valarray<double> Profiler::evaluatePredictions() {
    update_theory_iteration_();

    int ndata =  cndatapoints_.npoints;
    
    std::valarray<double> out(ndata);

    for (int i =0; i<ndata; i++) {
      out[i] = c_theo_.theo[i];
    }
    return out;
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
    int nsysloc = systema_.nsys;

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
    };

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
      double val = node[i].as<double>();
      *ppar = val;
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
    //Sanity checks
    if ( !sets || !members  ) {
      hf_errlog(2018082401,"S: Profiler: missing set or member parameters for evolution "+evolName);  // XXXXXXXXXXXXXXXX
    }
    if(!sets.IsSequence()||!members.IsSequence()){
      hf_errlog(2018082402,"S: Profiler: sets and members must be sequence");  // XXXXXXXXXXXXXXXX
    }
    if(sets.size()!=members.size()){
      hf_errlog(2018082403,"S: Profiler: sets and members must be the same length");  // XXXXXXXXXXXXXXXX
    }


    size_t endi = sets.size();
    for (size_t i=0; i< endi; i++) {
      std::string pName = sets[i].as<string>();

      if ( members[i].IsSequence() ) {
        int msize = members[i].size();

        if ( msize != 3) {
          hf_errlog(2018082404,"S: Profiler: sets must be sequence of length 3");  // XXXXXXXXXXXXXXXX
        }

        int central = members[i][0].as<int>();
        int first   = members[i][1].as<int>();
        int last = 0;

        try {
          last =  members[i][2].as<int>();
        }
        catch (...) {
          last = 0; /// auto-decodez
        }
          
        // save original
        auto oSet   =Clone(gNode["set"]);
        auto oMember=Clone(gNode["member"]);

        if ( ! oSet  || ! oMember ) {
          hf_errlog(2018082410,"S: No set or member variables for evolution : "+evolName);
        }//This should be evolution's problem, not Profiler's --Ivan

        // all predictions
        std::vector< std::valarray<double> > preds;

        // Set central PDF and init 
        gNode["set"]   =pName;
        gNode["member"]=central;
        evol->atConfigurationChange();
        preds.push_back(evaluatePredictions() );


        // now we can get set properties: is it hessian asymmetric, MC or symmetric hessian
        std::string errorType = evol->getPropertyS("ErrorType");
        std::cout << errorType << std::endl;

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
        }

        // Restore original

        gNode["set"]   =oSet;
        gNode["member"]=oMember;
        evol->atConfigurationChange();


        // Depending on error type, do nuisance parameters addition
        if ( errorType == "symmhessian" ) {
          for (int imember = first; imember<=last; imember++) {
            addSystematics("PDF_nuisance_param_"+std::to_string( ++_ipdf )+":T",(preds[imember]-preds[0])/preds[0]);
          }
        }
        else if ( errorType == "hessian") {
          for (int imember = first; imember<=last; imember += 2) {
            addSystematics("PDF_nuisance_param_"+std::to_string( ++_ipdf )+":T"
                           ,(preds[imember]-preds[0])/preds[0]
                           ,(preds[imember+1]-preds[0])/preds[0]);
          }
        }
        else if ( errorType == "replicas") {
          // construct average 
          for (size_t i=0; i<preds[0].size(); i++) {
            preds[0][i] = 0;
          }
          for ( int i=first; i<=last; i++) {
            preds[0] += preds[i];
          }
          preds[0] /= (last-first+1);
          
          // convert replicas to deviations from average:
          for ( int i=first; i<=last; i++) {
            preds[i] -= preds[0];
          }

          // convert to eigenvectors, add to list of systematics
          addReplicas(pName,preds);
        }
        else {
          hf_errlog(2018082441,"S: Profiler Unsupported PDF error type : "+errorType);
        }
        
      }
      else {
        hf_errlog(2018082404,"S: Profiler: sets must be sequence of length 3");  // XXXXXXXXXXXXXXXX
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
    delete beta;
    delete covar;
  }
  
} //namespace xfitter

