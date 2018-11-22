
/*
   @file CERESMinimizer.cc
   @date 2018-08-17
   @author  AddMinimizer.py
   Created by  AddMinimizer.py on 2018-08-17
*/

#include "CERESMinimizer.h"
#include "xfitter_cpp.h"
#include "xfitter_pars.h"
#include "xfitter_steer.h"

#include "ceres/ceres.h"
#include "glog/logging.h"
#include "ceres/dynamic_numeric_diff_cost_function.h"

// Fortran interface
extern "C" {
  // FCN
  void fcn_(const int& npar, const double& dummy, double& chi2out, const double* pars, const int& iflag, const double& dummy2);
  void iofilenamesmini_();
}

namespace xfitter {

/// the class factories, for dynamic loading
extern "C" CERESMinimizer* create() {
    return new CERESMinimizer();
}

/// connect fcn and minimizer
  void myFCN(double &chi2, double const* par, int iflag=2) {

    static int counter = 0;
    counter++;
    BaseMinimizer* mini = xfitter::get_minimizer();
    mini->setPars(par);

    int npar = mini->getNpars();

    double pp[200];  // 200 is needed for fcn ...
    for (int i=0; i<npar; i++) {
      pp[i] = par[i];
    }


    fcn_(npar, 0, chi2, pp, iflag, 0);
    std::cout << " Call again " << counter << " chi2=" << chi2;
    return;
}

//
// Connect fcn and xfitter pars
//
struct CostFunctiorData
{
   bool operator()(double const* const* parameters, double* residuals) const
   {
     double chi2;
     myFCN(chi2,parameters[0],2);
     chi2 = 0;
     for (int i = 0; i< cndatapoints_.npoints  + systema_.nsys; i++) {
       residuals[i] = c_resid_.residuals_[i];
       chi2 += residuals[i]*residuals[i];
     }
     std::cout << " residuals squared: " <<  chi2 << std::endl;
     return true;
   }
};



// Constructor
CERESMinimizer::CERESMinimizer() : BaseMinimizer("CERES")
{
}

// Constructor
CERESMinimizer::CERESMinimizer(const std::string& inName) : BaseMinimizer(inName)
{
}

// Init at start:
void CERESMinimizer::atStart() {
  //There's a suppressed warning here
  //warning: deprecated conversion from string constant to ‘char*’ [-Wwrite-strings]
  google::InitGoogleLogging((char*)"");
  // also init xfitter logging:
  iofilenamesmini_();
  return;
}

/// Minimization loop
void CERESMinimizer::doMinimization()
{
  int nData = cndatapoints_.npoints;
  int nSyst = systema_.nsys;
  int npars =  getNpars() ;
  std::cout << nData << " " << nSyst << " " << npars << "\n";

   // Cost function:
  ceres::DynamicNumericDiffCostFunction<CostFunctiorData>* dynamic_cost_function =
    new ceres::DynamicNumericDiffCostFunction<CostFunctiorData>(new CostFunctiorData);

  dynamic_cost_function->AddParameterBlock(npars);

  double parVals[npars];
  double**pars = getPars();
  for (int i =0; i<npars; i++) {
    parVals[i] = *pars[i];
    std::cout << " par " <<  parVals[i] << std::endl;
  }

  dynamic_cost_function->SetNumResiduals(nData+nSyst);

  ceres::Problem myProblem;
  myProblem.AddResidualBlock(dynamic_cost_function, NULL, parVals);

  ceres::Solver::Options myOptions;
  myOptions.minimizer_progress_to_stdout = true;

  myOptions.function_tolerance = 1.e-5; // typical chi2 is ~1000
  ceres::Solver::Summary mySummary;

  ceres::Solve(myOptions, &myProblem, &mySummary);

  // after mini actions
  double chi2;
  myFCN(chi2, parVals, 3);

  std::cout << mySummary.FullReport() << "\n";
  //remember convergence status to report later
  switch(mySummary.termination_type){
    case ceres::CONVERGENCE:
      convergence_status=ConvergenceStatus::SUCCESS;
      break;
    case ceres::NO_CONVERGENCE:
      convergence_status=ConvergenceStatus::NO_CONVERGENCE;
      break;
    case ceres::USER_SUCCESS:
      convergence_status=ConvergenceStatus::SUCCESS;
      break;
    default:
      convergence_status=ConvergenceStatus::ERROR;
  }
  return;
}

/// Action at last iteration
void CERESMinimizer::actionAtFCN3()
{
    return;
}

/// Error analysis
void CERESMinimizer::errorAnalysis()
{
    return;
}

ConvergenceStatus CERESMinimizer::convergenceStatus(){
  return convergence_status;
}

/// Parameter transfer
void CERESMinimizer::addParameter(double par, std::string const &name, double step, double const* bounds , double  const* priors  )
{
  BaseMinimizer::addParameter(par,name,step,bounds,priors);
  return;
}


}

