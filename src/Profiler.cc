#include "Profiler.h"
#include "xfitter_cpp_base.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "xfitter_cpp.h"
#include "xfitter_steer.h"
#include "BaseEvolution.h"
#include "xfitter_cpp.h"
#include <algorithm>
#include <string.h>
#include <numeric>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>

extern "C" {
  void update_theory_iteration_();
  void addsystematics_(char const* name, int len);
  void store_pdfs_(char const* base, int len);
  void writefittedpoints_();
}
  
namespace xfitter
{
  vector <double> mcweights(vector<double> const& chi2, int ndata, bool GK_method)
  {
    vector<double> w;
    const int nrep = chi2.size();
  
    vector<double> logw;  // Calculate Ln(wnn) (nn - not normalised) -> use repnums as index for unweighted PDFs
    for (vector <double>::const_iterator c = chi2.begin(); c != chi2.end(); c++) {
      if ( (GK_method) == false) {
	logw.push_back( - (*c)/(2.0) +( (((double) ndata)-1.0)/2.)*log(*c)); //Bayesian
      }
      if ( (GK_method) == true)  logw.push_back(- (*c)/(2.0)); //Giele-Keller
    }
    // Get maximum value of log(w)
    double exp_avg = *(max_element(logw.begin(), logw.end()));
    //exp_avg =0;//only needed to normalise
    // Calculate weights
    for (size_t i=0 ;i<nrep; i++) {
      w.push_back(exp(logw[i] - exp_avg ));
    }
    //Drop any weights smaller than 1e-12
    double wtot = std::accumulate(w.begin(),w.end(),0.0); 
    for (size_t i=0; i<w.size(); i++) {
      if ((w[i]*(nrep/wtot)) < 1e-12)
	w[i]=0;
    }
    wtot=std::accumulate(w.begin(),w.end(),0.0);  // Normalise weights so Sum(weights)=N
    for (size_t i=0;i<w.size();i++){
      w[i]*=(nrep/wtot); 
    }
    return w;
  }
  
  std::pair < std::valarray<double>, double> Profiler::evaluatePredictions() {
    double chi2 = NAN;
    if (_getChi2) {
      int save_nsys = systema_.nsys;
      systema_.nsys = _nSourcesOrig; // reset to the original sources
      chi2 = chi2data_theory_(2);
      systema_.nsys = save_nsys;
      std::cout << std::fixed << std::setprecision(2) << "Chi2 = " << chi2 << std::endl;
    }else{
      update_theory_iteration_();
    }

    //Return theory predictions
    return std::pair<valarray<double>,double> ( valarray<double>(c_theo_.theo, cndatapoints_.npoints), chi2);
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

    if (node["threads"]) {
      _ncpu =  node["threads"].as<int>();
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
        auto cent = evaluatePredictions().first;
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
      auto pred = evaluatePredictions();
      preds.push_back( pred.first );
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


  void Profiler::compute_parallel(int NALL, int NPRED, int first, int iPdfSet,
			std::vector< std::valarray<double> >& preds,
			std::vector< double >& chi2vals,
				  YAML::Node gNode,
				  BaseEvolution* evol, const std::string& errorType)
  {
    
    // Shared memory for predictions
    int shmid;
    double* sharedArray;
    int ARRAY_SIZE = preds.size() * NPRED;
	
    // Create shared memory segment
    shmid = shmget(IPC_PRIVATE, sizeof(double) * ARRAY_SIZE, IPC_CREAT | 0666);
    if (shmid < 0) {
      hf_errlog(2023060200,"F: Failed to create shared memory segment");
    }
	
    // Attach shared memory segment
    sharedArray = static_cast<double*>(shmat(shmid, nullptr, 0));
    if (sharedArray == reinterpret_cast<double*>(-1)) {
      hf_errlog(2023060201,"F: Failed to attach shared memory segment");
    }

    // Shared memory for chi2s
    int shmid2;
    double* sharedArray2;
    int ARRAY_SIZE2 = chi2vals.size();
	
    // Create shared memory segment
    shmid2 = shmget(IPC_PRIVATE, sizeof(double) * ARRAY_SIZE2, IPC_CREAT | 0666);
    if (shmid < 0) {
      hf_errlog(2023060202,"F: Failed to create shared memory segment");
    }
	
    // Attach shared memory segment
    sharedArray2 = static_cast<double*>(shmat(shmid2, nullptr, 0));
    if (sharedArray == reinterpret_cast<double*>(-1)) {
      hf_errlog(2023060203,"F: Failed to attach shared memory segment");
    }


    // define Chunks

    std::cout << "N CPU: " << _ncpu << std::endl;

    int NCPU = _ncpu;
    int chunkSize = NALL / NCPU;
    int reminder  = NALL % NCPU; 
    int startIndex = 0;
    int endIndex = 0;
    
    // loop over all
    for (int icpu = 0; icpu<min(NCPU,NALL); icpu++) {
      startIndex = endIndex;
      endIndex   = startIndex + chunkSize;
      if (icpu < reminder) {
	endIndex += 1;
      }
      pid_t pid = fork();
      if ( pid == 0) {       
	for (int imember = first+startIndex; imember < first+endIndex; imember++) {
	  
	  gNode["member"] = imember;
	  evol->atConfigurationChange();
	  auto pred = evaluatePredictions();
	  
	  // store in shared memory
	  int idxOff = (imember-first)*NPRED;
	  for (size_t idx = 0; idx<NPRED; idx++) {
	    sharedArray[idxOff+idx] = pred.first[idx];
	  }
	      
	  sharedArray2[imember-first] = pred.second;

	  if (_storePdfs) {
	    storePdfFiles(imember,iPdfSet,errorType);
	  }
	      
	}
	exit(0);	    
      }
      else if (pid<0) {
	hf_errlog(2023060204,"F: Failed to create a fork process");	
      }
    }
	
    // Wait ...
    int status;
    while (wait(&status) > 0);
    
    // Store in vectors
    for (size_t imember = 0; imember<NALL; imember++) {
      chi2vals[imember] = sharedArray2[imember];
      
      std::valarray<double> temp;
      temp.resize(NPRED);
      for (size_t ipred = 0; ipred < NPRED; ipred++) {	    
	temp[ipred] = sharedArray[imember*NPRED+ipred];
      }
      preds[imember+1] = temp;
    }
    
    // Detach and remove shared memory segments
    shmdt(sharedArray);
    shmctl(shmid, IPC_RMID, NULL);
    shmdt(sharedArray2);
    shmctl(shmid2, IPC_RMID, NULL);
    
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
	std::vector< double > chi2vals;

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

	// Central+all PDFs
	preds.resize(last-first+2);
	// All PDFs only
	chi2vals.resize(last-first+1);
	
        auto pred = evaluatePredictions();
        preds[0] = pred.first;

	int NPRED = pred.first.size();
	int NALL = last-first+1;
	
	/// DO NOT store the central chi2:
	/// chi2vals.push_back(pred.second);
	
	if (_storePdfs) {
	  storePdfFiles(0,i);
	}

	if (_ncpu>0) {
	  // Use multiprocessing
	  compute_parallel(NALL, NPRED, first, i, preds, chi2vals, gNode, evol, errorType);
	}
	else {
	  // loop over all
	  for (int imember = first; imember<=last; imember++) {
	    gNode["member"] = imember;
	    evol->atConfigurationChange();
	    auto pred = evaluatePredictions();
	    preds[imember-first+1] = pred.first;
	    chi2vals[imember-first] = pred.second;
	    if (_storePdfs) {
	      storePdfFiles(imember,i,errorType);
	    }
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

	  // Write out aux files for reweighting method.
	  addMCweightsFiles(pName,chi2vals,preds[0].size(),preds.size()-1);
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

  void Profiler::addMCweightsFiles(std::string const& pdfName, std::vector<double>& chi2vals, int ndata, int nrep) {
    vector <double> weights = mcweights(chi2vals, ndata, false);
    std::ofstream chi2wf((_outputDir + "/pdf_BAYweights.dat").c_str());
    vector <double>::iterator ic = chi2vals.begin();
    vector <double>::iterator iw = weights.begin();
    string lhapdfsetname=pdfName;
    
    if (lhapdfsetname.find(".LHgrid")!=string::npos)   lhapdfsetname.erase(lhapdfsetname.find(".LHgrid"));
    chi2wf << "LHAPDF set=   " << lhapdfsetname<<endl;
    chi2wf << "Reweight method=   BAYESIAN"<<endl;
    chi2wf << "ndata=   " << ndata <<endl;
    chi2wf << nrep << endl;
    for (; ic != chi2vals.end(); ic++, iw++)
      chi2wf << (ic - chi2vals.begin()) << "\t" << *ic << "\t" << *iw << endl;
    chi2wf.close();
    
    vector <double> weights_GK = mcweights(chi2vals, ndata, true);
    std::ofstream chi2wf2((_outputDir + "/pdf_GKweights.dat").c_str());
    ic = chi2vals.begin();
    iw = weights_GK.begin();
    chi2wf2 << "LHAPDF set=   " << lhapdfsetname<<endl;
    chi2wf2 << "Reweight method=   GIELE-KELLER"<<endl;
    chi2wf2 << "ndata=   " << ndata <<endl;
    chi2wf2 << nrep << endl;
    for (; ic != chi2vals.end(); ic++, iw++)
      chi2wf2 << (ic - chi2vals.begin()) << "\t" << *ic << "\t" << *iw << endl;
    chi2wf2.close();
  }
  
} //namespace xfitter

