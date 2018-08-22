 
/*
   @file MINUITMinimizer.cc
   @date 2018-08-17
   @author  AddMinimizer.py
   Created by  AddMinimizer.py on 2018-08-17
*/

#include "MINUITMinimizer.h"
#include "FTNFitPars.h"          // tools to handle fortran minuit
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"

/// Fortran interfaces:
extern "C" {
  void generate_io_filenames_();
  void minuit_ini_();
  void minuit_(void fcn(const int&, const double&, double&, const double*, const int& , const double& ), int const& i);
  // FCN
  void fcn_(const int& npar, const double& dummy, double& chi2out, const double* pars, const int& iflag, const double& dummy2);  

  void mncomd_(void fcn(const int&, const double&, double&, const double*, const int& , const double& ),
	       const char command[], int& icond, int const& iin, int);
  
  /// Interface to minuit parameters
  void addexternalparam_(const char name[],  const double &val, 
                         const double  &step,
                         const double &min, const double &max, 
                         const double &prior, const double &priorUnc,
                         const int &add, 
                         map<std::string,double*> *map,
                         int len);

  ///
  void extraparam_();

  ///
  void  write_pars_(int const& i);
  void  error_bands_pumplin_();

  void  errbandssym_();
}

namespace xfitter {
  
/// the class factories, for dynamic loading
extern "C" MINUITMinimizer* create() {
    return new MINUITMinimizer();
}


// Constructor
MINUITMinimizer::MINUITMinimizer() : BaseMinimizer("MINUIT") {  
}

// Constructor
MINUITMinimizer::MINUITMinimizer(const std::string& inName) : BaseMinimizer(inName) 
{  
}

// Init at start:
void MINUITMinimizer::initAtStart() {
  // Call fortran interface
  generate_io_filenames_();
  minuit_ini_();
  return;
}

/// Miniminzation loop
void MINUITMinimizer::doMimimization() 
{
  extraparam_();
  auto text = XFITTER_PARS::gParametersY["MINUIT"]["Commands"];
  std::string cmds = text.as<string>();
  std::istringstream lineStream(cmds);
  std::string line;
  //  minuit_(fcn_,0);  // let fortran run ...
  while (std::getline(lineStream, line,'\n'))
    {
      int ires; 
      mncomd_(fcn_,line.c_str(),ires,0,line.size());
    }
  return;
}

/// Action at last iteration 
void MINUITMinimizer::actionAtFCN3() 
{
  //  double**p = getPars();
  //for (size_t i=0; i<getNpars(); i++) {
  //  std::cout << i << " " << *p[i] << "\n"; 
  //}
  //  exit(0);
  return;
}

/// Error analysis
void MINUITMinimizer::errorAnalysis() 
{
  auto errNode = XFITTER_PARS::gParametersY["MINUIT"]["doErrors"];
  if ( errNode ) {
    std::string bandType = errNode.as<std::string>();
    if ( bandType == "Pumplin" ) {
      hf_errlog(12020506, "I: Calculation of error bands required");
      std::string cmd = "ITERATE 10";
      int ires;
      mncomd_(fcn_, cmd.c_str(), ires, 0, cmd.size());
      cmd = "MYSTUFF 1000";
      mncomd_(fcn_, cmd.c_str(), ires, 0, cmd.size());
      cmd = "MYSTUFF 2000";
      mncomd_(fcn_, cmd.c_str(), ires, 0, cmd.size());

      write_pars_(0);
      error_bands_pumplin_();
    }
    else if ( bandType == "Hesse" ) {
      hf_errlog(12020506, "I: Calculation of symmetric error bands required");
      errbandssym_();
    }
    else if ( bandType == "None" ) {
      return;
    }
    else {
      hf_errlog(2018092201,"W: Unknown to MINUITminimizer error type requested: " + bandType);
    }
  }
  return;
}

/// parameters  
void MINUITMinimizer::addParameter(double par, std::string const &name, double step, double const* bounds , double  const* priors  )
{
  BaseMinimizer::addParameter(par,name,step,bounds,priors);
  double minv    =0;
  double maxv    =0;
  double priorVal=0;
  double priorUnc=0;
  int add = true;

  if ( bounds != nullptr ) {
    minv = bounds[0];
    maxv = bounds[1];
  }
  if ( priors != nullptr ) {
    priorVal = priors[0];
    priorUnc = priors[1];
  }
  addexternalparam_(name.c_str(),par,step,minv,maxv,priorVal,priorUnc,add,&XFITTER_PARS::gParameters,name.size());
  
  return;
}

  
} //namespace xfitter

