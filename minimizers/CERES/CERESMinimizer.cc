
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
#include "ceres/covariance.h"

#include <iomanip>
#include <fstream>

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
    std::cout << " Iteration " << counter << " FCN chi2=" << chi2;
    return;
}

//
// Connect fcn and xfitter pars
//
struct CostFunctiorData
{
   bool operator()(double const* const* parameters, double* residuals) const
   {
     cout << std::endl;
     double chi2;
     myFCN(chi2,parameters[0],2);
     chi2 = 0;
     for (int i = 0; i< cndatapoints_.npoints  + systema_.nsys; i++) {
       //residuals[i] = c_resid_.residuals_[i];
       if (i < cndatapoints_.npoints)
	 {
	   //std::cout << " res " << i << " " << c_resid_.residuals_[i] << " logchi2 " << cdatapoi_.chi2_poi_data_[i] << " sum2 " << pow(c_resid_.residuals_[i],2)+cdatapoi_.chi2_poi_data_[i] << std::endl;
	   residuals[i] = sqrt(max(0.,pow(c_resid_.residuals_[i],2)+cdatapoi_.chi2_poi_data_[i]));
	 }
       else
	 residuals[i] = c_resid_.residuals_[i];
       chi2 += residuals[i]*residuals[i];
     }
     //residuals[cndatapoints_.npoints  + systema_.nsys] = cdatapoi_.chi2_poi_tot_;
     //chi2 += residuals[cndatapoints_.npoints  + systema_.nsys]*residuals[cndatapoints_.npoints  + systema_.nsys];
     std::cout << " CERES residuals squared: " <<  chi2 << std::endl;
     return true;
   }
};


struct GUM {
  bool operator()(const double* parameters, double* cost) const {
     double chi2;
     myFCN(chi2,parameters,2);
     cost[0] = chi2;
     std::cout << " cost function: " <<  chi2 << std::endl;
     return true;
  }
  static ceres::FirstOrderFunction* Create(int npars) {
    const int kNumParameters = 14;
    return new ceres::NumericDiffFirstOrderFunction<GUM,
                                                    ceres::CENTRAL,
                                                    kNumParameters>(
        new GUM);
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

  double**pars = getPars();
  
  //double parVals[npars];
  double parVals[14];

  for (int i =0; i< npars; i++)
    parVals[i] = *pars[i];

  double covmat[npars * npars];
  fill(covmat,covmat+npars*npars, 0.);
  
  int strategy = 0;
  //Least squares
  if (strategy == 0)
    {
      ceres::Solver::Options myOptions;
      ceres::Solver::Summary mySummary;
  
      // Cost function:
      ceres::DynamicNumericDiffCostFunction<CostFunctiorData>* dynamic_cost_function =
	new ceres::DynamicNumericDiffCostFunction<CostFunctiorData>(new CostFunctiorData);

      dynamic_cost_function->AddParameterBlock(npars);

      dynamic_cost_function->SetNumResiduals(nData+nSyst);

      ceres::Problem myProblem;
      myProblem.AddResidualBlock(dynamic_cost_function, NULL, parVals);

      myOptions.minimizer_progress_to_stdout = true;

      myOptions.function_tolerance = 1.e-5; // typical chi2 is ~1000

      ceres::Solve(myOptions, &myProblem, &mySummary);

      //covariance
      ceres::Covariance::Options options;
      //options.algorithm_type = ceres::SPARSE_QR;
      //options.algorithm_type = ceres::DENSE_SVD;
      ceres::Covariance covariance(options);

      vector<pair<const double*, const double*> > covariance_blocks;
      covariance_blocks.push_back(make_pair(parVals, parVals));

      CHECK(covariance.Compute(covariance_blocks, &myProblem));

      covmat[npars * npars];
      covariance.GetCovarianceBlock(parVals, parVals, covmat);

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

      writePars(covmat);
      writeOutput(mySummary, covmat);
    }

  //General Unconstrained Minimization (to be used with PoissonCorr option)
  else if (strategy == 1)
    {
      ceres::GradientProblem problem(GUM::Create(npars));

      ceres::GradientProblemSolver::Options options;
      options.minimizer_progress_to_stdout = true;
      options.function_tolerance = 1.e-5;
      ceres::GradientProblemSolver::Summary Summary;
      ceres::Solve(options, problem, parVals, &Summary);

      std::cout << Summary.FullReport() << "\n";
    }
  
  // after mini actions
  double chi2;
  myFCN(chi2, parVals, 3);

  cout << endl;
  cout << "Fitted parameters" << endl;
  for (int i = 0; i < npars; i++)
    std::cout << setw(5) << i << setw(15) << _allParameterNames[i] << setw(15) <<  parVals[i] << " +/- " << sqrt(covmat[i*npars+i]) << std::endl;
  
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

void CERESMinimizer::writePars(const double* covmat)
{
  int npars =  getNpars() ;
  double parVals[npars];
  double**pars = getPars();
  for (int i =0; i<npars; i++)
    parVals[i] = *pars[i];

  std::ofstream f;
  f.open(stringFromFortran(coutdirname_.outdirname,sizeof(coutdirname_.outdirname))+"/parsout_0");
  if(!f.is_open()){
    hf_errlog(16042807,"W: Failed to open parsout_0 for writing");
    return;
  }

  for (int i =0; i<npars; i++)
    f << setw(5) << i << "   '" << _allParameterNames[i] << "'    " <<  parVals[i] << "    " << sqrt(covmat[i*npars+i]) << std::endl;
  
}

void CERESMinimizer::writeOutput(ceres::Solver::Summary mySummary, const double* covmat)
{
  std::ofstream f;
  f.open(stringFromFortran(coutdirname_.outdirname,sizeof(coutdirname_.outdirname))+"/ceres.out.txt");
  if(!f.is_open()){
    hf_errlog(16042807,"W: Failed to open ceres.out.txt for writing");
    return;
  }
    
  f << mySummary.FullReport() << "\n";

  int npars =  getNpars() ;

  //Write Covariance matrix
  f << std::endl << "COVARIANCE MATRIX " << std::endl;
  f << std::setw(14) << " ";
  for (int i = 0; i < npars; i++)
    f << std::setw(14) << _allParameterNames[i];
  f << std::endl;
  f << setprecision(8);
  for (int i = 0; i < npars; i++)
    {
      f << std::setw(14) << _allParameterNames[i];;
      for (int j =0; j<npars; j++)
	f << std::setw(14) << covmat[i*npars+j];
      f << std::endl;
    }

  //Write Correlation matrix
  f << std::endl << "CORRELATION MATRIX " << std::endl;
  f << std::setw(14) << " ";
  for (int i = 0; i < npars; i++)
    f << std::setw(14) << _allParameterNames[i];
  f << std::endl;
  f << setprecision(4);
  for (int i = 0; i < npars; i++)
    {
      f << std::setw(14) << _allParameterNames[i];;
      for (int j =0; j<npars; j++)
	f << std::setw(14) << covmat[i*npars+j]/sqrt(covmat[i*npars+i])/sqrt(covmat[j*npars+j]);
      f << std::endl;
    }
    
}

}

