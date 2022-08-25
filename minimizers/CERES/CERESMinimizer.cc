
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

#include <Eigen/Dense>

#include <iomanip>
#include <fstream>

//make this a yaml setting
const double glboff = 2.;

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

    double pp[200];  // 200 is needed for fcn ... //--> should use NEXTRAPARAMMAX_C from dimensions.h, and synchronize with MNE from endmini.inc
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
     double totoff = 0.;
     int nres = (chi2options_.chi2poissoncorr) ? cndatapoints_.npoints + systema_.nsys + cndatapoints_.npoints : cndatapoints_.npoints + systema_.nsys;
     for (int i = 0; i < nres; i++)
       {
	 if (i < cndatapoints_.npoints + systema_.nsys)
	   residuals[i] = sqrt(2.)*c_resid_.residuals_[i];
	 else
	   //Log penalty terms
	   {
	     //Add a positive offset to each squared residual to enforce a positive cost function also when the log penalty term is included
	     totoff += glboff;
	     int j = i-systema_.nsys-cndatapoints_.npoints;
	     residuals[i] = sqrt(2.)*sqrt(max(0.,cdatapoi_.chi2_poi_data[j]+glboff));
	     if (cdatapoi_.chi2_poi_data[j]+glboff < 0)
	       {
		 cout << "Warning: negative squared residual for i " << j << " offset " << glboff << " res^2 " << cdatapoi_.chi2_poi_data[j]+glboff << endl;
		 string message = "W: Negative squared residual in CERES minimisation, consider raising the chi-square offset";
		 hf_errlog_(22082101, message.c_str(), message.size());
	       }
	   }
	 chi2 += residuals[i]*residuals[i]/2.;
       }
     std::cout << " CERES residuals squared: " <<  chi2 - totoff << std::endl;
     return true;
   }
};

class PenaltyLog final: public ceres::LossFunction
{
public:
  void Evaluate(double s, double out[3]) const override
  {
    out[0] = s;
    if (chi2options_.chi2poissoncorr)
      out[0] += -2.*glboff*cndatapoints_.npoints;
    out[1] = 1.;
    out[2] = 0.;
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
  cout << endl;
  cout << "CERES minimisation" << endl;
  cout << "Data points: " << nData << "; Systematic uncertainties nuisance parameters: " << nSyst << "; Free parameters: " << npars;
  if (chi2options_.chi2poissoncorr)
    cout << "; Penalty log terms: " << nData;
  cout << endl;

  double**pars = getPars();
  
  double parVals[200]; // 200 --> should use NEXTRAPARAMMAX_C from dimensions.h, and synchronize NEXTRAPARAMMAX_C with MNE from endmini.inc

  for (int i = 0; i < npars; i++)
    parVals[i] = *pars[i];

  double covmat[npars * npars];
  fill(covmat,covmat+npars*npars, 0.);

  //First call to FCN for initialisation
  double chi2;
  myFCN(chi2,parVals,3);
  
  // Least squares minimisation
  ceres::Solver::Options soloptions;
  ceres::Solver::Summary summary;
  
  // Cost function:
  ceres::DynamicNumericDiffCostFunction<CostFunctiorData>* dynamic_cost_function =
    new ceres::DynamicNumericDiffCostFunction<CostFunctiorData>(new CostFunctiorData);

  dynamic_cost_function->AddParameterBlock(npars);

  int nres = (chi2options_.chi2poissoncorr) ? cndatapoints_.npoints + systema_.nsys + cndatapoints_.npoints : cndatapoints_.npoints + systema_.nsys;

  dynamic_cost_function->SetNumResiduals(nres);

  // Loss function
  ceres::LossFunction* loss_function(new PenaltyLog);

  ceres::Problem problem;
  problem.AddResidualBlock(dynamic_cost_function, loss_function, parVals);

  soloptions.minimizer_progress_to_stdout = true;

  soloptions.function_tolerance = 1.e-5; // typical chi2 is ~1000

  // --> Allow setting options from yaml
  /*
  soloptions.logging_type = ceres::SILENT;
  
  //multithreading
  //soloptions.num_threads = 4;

  //Ceres options
  soloptions.max_num_iterations = 100000;

  soloptions.minimizer_type = ceres::TRUST_REGION; //ceres::LINE_SEARCH;
  soloptions.linear_solver_type = ceres::DENSE_QR; //ceres::SPARSE_NORMAL_CHOLESKY;
  soloptions.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT; //ceres::DOGLEG;
  //soloptions.dogleg_type = ceres::SUBSPACE_DOGLEG;
  //soloptions.use_nonmonotonic_steps = true;
  
  //soloptions.max_num_consecutive_invalid_steps = 10;

  soloptions.dense_linear_algebra_library_type = ceres::EIGEN; //ceres::LAPACK
  */
  
  ceres::Solve(soloptions, &problem, &summary);

  cout << std::endl;
  cout << "CERES minimisation has converged" << std::endl;
  //if (covariance matrix is required)
  cout << std::endl;
  cout << "CERES Start calculation of covariance matrix" << std::endl;
      
  //covariance
  ceres::Covariance::Options covoptions;

  //covoptions.algorithm_type = ceres::SPARSE_QR;
  //covoptions.algorithm_type = ceres::DENSE_SVD;

  //multithreading
  //covoptions.num_threads = 4;
  ceres::Covariance covariance(covoptions);

  vector<pair<const double*, const double*> > covariance_blocks;
  covariance_blocks.push_back(make_pair(parVals, parVals));

  CHECK(covariance.Compute(covariance_blocks, &problem));

  covariance.GetCovarianceBlock(parVals, parVals, covmat);

  std::cout << summary.FullReport() << "\n";

  //remember convergence status to report later
  switch(summary.termination_type){
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
  writeOutput(summary, covmat);

  // after mini actions
  myFCN(chi2, parVals, 3);

  cout << endl;
  cout << "Fitted parameters" << endl;
  for (int i = 0; i < npars; i++)
    std::cout << setw(5) << i << setw(15) << _allParameterNames[i] << setw(15) <<  parVals[i] << " +/- " << sqrt(covmat[i*npars+i]) << std::endl;

  cout << endl;
  cout << "----- Parameters in YAML format (can copy-paste into parameters.yaml):" << endl;
  cout << "  Parameters:" << endl;
  for (int i = 0; i < npars; i++)
    {
      char val[15];
      char err[15];
      sprintf(val, "%.4f", parVals[i]);
      sprintf(err, "%.4f", sqrt(covmat[i*npars+i]));
      std::cout << "  " << _allParameterNames[i] << ": [ " << val << ", " << err << " ]" << std::endl;
    }
  cout << " ----- End of parameters in YAML format" << endl; 
  cout << endl;
  cout << std::endl << "Correlation matrix " << std::endl;
  cout << std::setw(14) << " ";
  for (int i = 0; i < npars; i++)
    cout << std::setw(14) << _allParameterNames[i];
  cout << std::endl;
  cout << setprecision(4);
  for (int i = 0; i < npars; i++)
    {
      cout << std::setw(14) << _allParameterNames[i];;
      for (int j =0; j<npars; j++)
	cout << std::setw(14) << covmat[i*npars+j]/sqrt(covmat[i*npars+i])/sqrt(covmat[j*npars+j]);
      cout << std::endl;
    }
  
  Eigen::MatrixXd cov(npars,npars);
  for (int i = 0; i < npars; i++)
    for (int j = 0; j < npars; j++)
      cov(i,j) = covmat[i+npars*j];
  
  Eigen::MatrixXd inv = cov.inverse();
  double rhok[npars];
  for (int i = 0; i < npars; i++)
    rhok[i] = sqrt(1. - 1./(cov(i,i)*inv(i,i)));

  cout << endl;
  cout << "Global correlation coefficient" << endl;
  for (int i = 0; i< npars; i++)
    cout << setw(5) << i << setw(15) << _allParameterNames[i] << setw(15) <<  rhok[i] << endl;
  
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

void CERESMinimizer::writeOutput(ceres::Solver::Summary summary, const double* covmat)
{
  std::ofstream f;
  f.open(stringFromFortran(coutdirname_.outdirname,sizeof(coutdirname_.outdirname))+"/ceres.out.txt");
  if(!f.is_open()){
    hf_errlog(16042807,"W: Failed to open ceres.out.txt for writing");
    return;
  }
    
  //f << summary.FullReport() << "\n";

  int npars =  getNpars() ;

  //Write Covariance matrix
  f << std::endl << "COVARIANCE MATRIX " << std::endl;
  f << std::setw(14) << " ";
  for (int i = 0; i < npars; i++)
    f << std::setw(14) << _allParameterNames[i];
  f << std::endl;
  f << setprecision(4);
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

  Eigen::MatrixXd cov(npars,npars);
  for (int i = 0; i < npars; i++)
    for (int j = 0; j < npars; j++)
      cov(i,j) = covmat[i+npars*j];
  
  Eigen::MatrixXd inv = cov.inverse();
  double rhok[npars];
  for (int i = 0; i < npars; i++)
    rhok[i] = sqrt(1. - 1./(cov(i,i)*inv(i,i)));

  f << "GLOBAL CORRELATION COEFFICIENT" << endl;
  for (int i = 0; i< npars; i++)
    f << setw(5) << i << setw(15) << _allParameterNames[i] << setw(15) <<  rhok[i] << endl;
}

}

