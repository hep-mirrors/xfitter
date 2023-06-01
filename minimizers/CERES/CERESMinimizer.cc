
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
#include <unistd.h>
#include <sys/wait.h>



// Fortran interface
extern "C" {
  // FCN
  void fcn_(const int& npar, const double& dummy, double& chi2out, const double* pars, const int& iflag, const double& dummy2);
  void iofilenamesmini_();
}

namespace xfitter {

double CERESMinimizer::glboff;
double *CERESMinimizer::offset;
double CERESMinimizer::totoffset;

/// the class factories, for dynamic loading
extern "C" CERESMinimizer* create() {
    return new CERESMinimizer();
}

/// connect fcn and minimizer
  void myFCN(double &chi2, double const* par, int iflag=2) {

    static int counter = 0;
    counter++;
    const BaseMinimizer* mini = xfitter::get_minimizer();
    mini->setPars(par);

    const int npar = mini->getNpars();

    double pp[200];  // 200 is needed for fcn ... //--> should use NEXTRAPARAMMAX_C from dimensions.h, and synchronize with MNE from endmini.inc
    for (int i=0; i<npar; i++) {
      pp[i] = par[i];
    }


    fcn_(npar, 0, chi2, pp, iflag, 0);
    std::cout << " Iteration " << counter << " FCN chi2=" << chi2 << std::endl;
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
     int nres = (chi2options_.chi2poissoncorr) ? cndatapoints_.npoints + systema_.nsys + cndatapoints_.npoints : cndatapoints_.npoints + systema_.nsys;
     bool calctotoff = true;
     if (CERESMinimizer::totoffset > 0.)
       calctotoff = false;
       
     for (int i = 0; i < nres; i++)
       {
	 if (i < cndatapoints_.npoints + systema_.nsys)
	   residuals[i] = sqrt(2.)*c_resid_.residuals_[i];
	 else
	   //Log penalty terms
	   {
	     int j = i-systema_.nsys-cndatapoints_.npoints;

	     //Add a positive offset to each squared residual to enforce a positive cost function also when the log penalty term is included
	     if (calctotoff)
	       {
		 if (CERESMinimizer::glboff > 0.)
		   CERESMinimizer::offset[j] = CERESMinimizer::glboff;
		 else
		   CERESMinimizer::offset[j] = 2.*max(0.,-cdatapoi_.chi2_poi_data[j]);
		 CERESMinimizer::totoffset += CERESMinimizer::offset[j];
	       }
	     residuals[i] = sqrt(2.)*sqrt(max(0.,cdatapoi_.chi2_poi_data[j]+CERESMinimizer::offset[j]));
	     if (cdatapoi_.chi2_poi_data[j]+CERESMinimizer::offset[j] < 0)
	       {
		 cout << "Warning: negative squared residual for i " << j << " offset " << CERESMinimizer::offset[j] << " res^2 " << cdatapoi_.chi2_poi_data[j]+CERESMinimizer::offset[j] << endl;
		 string message = "W: Negative squared residual in CERES minimisation, consider raising the chi-square offset";
		 hf_errlog_(22082101, message.c_str(), message.size());
	       }
    	   }
	 chi2 += residuals[i]*residuals[i]/2.;
       }
     //std::cout << " CERES total offset in Cost function: " <<  CERESMinimizer::totoffset << std::endl;
     std::cout << " CERES residuals squared: " <<  chi2 - CERESMinimizer::totoffset << std::endl;
     return true;
   }
};


  //
  // Connect fcn and xfitter pars
  //

  bool derivative(double const* const* parameters, double const* centralResiduals, int iPar, int nPar, int nRes, double* derivatives) {
    //
    double *pars = new double[nPar];
    double *resid = new double[nRes];
    for (int i=0; i<nPar; i+=1) {
      pars[i] = parameters[0][i];
    }
    
    
    double delta = abs(pars[iPar])<1.e-30 ? 1.e-30 :  abs(pars[iPar])* 1e-6;
    CostFunctiorData evaluate;
    pars[iPar] += delta;
    auto res = evaluate(&pars,resid);

    for (int i=0; i<nRes; i+=1) {
      derivatives[i] = (resid[i]-centralResiduals[i])/delta;
    }
    
    delete[] pars;
    delete[] resid;
    return true;
  }


  class CostFuntionrDataAndDerivative : public ceres::CostFunction
  {
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jakobian) const
    {     
      CostFunctiorData evaluate;
      // number of parameters:
      const int npar = parameter_block_sizes()[0];
      // number of residuals:
      const int nres = num_residuals();
     
      // central value:
      auto res= evaluate(parameters,residuals);

      for (int ipar=0; ipar<npar; ipar+=1) {
	pid_t pid = fork();
	if ( pid == 0) {
	  auto res = derivative(parameters,residuals, ipar, npar, nres, *(jakobian+ipar*nres));
	  exit(0);
	}
      }
      int status;
      while (wait(&status) > 0);
     
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
      //out[0] += -2.*glboff*cndatapoints_.npoints;
      out[0] += -2.*CERESMinimizer::totoffset;
    out[1] = 1.;
    out[2] = 0.;
    //std::cout << " CERES total offset in Loss function: " <<  CERESMinimizer::totoffset << std::endl;
    std::cout << " CERES Cost function: " <<  out[0]/2. << std::endl;
    CERESMinimizer::totoffset = 0.;
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

  YAML::Node ceresNode = XFITTER_PARS::rootNode["CERES"];
  glboff = ceresNode["offset"].as<double>();

  std::cout << " CERES global offset per bin: " <<  glboff << std::endl;

  if (chi2options_.chi2poissoncorr)
    offset = new double[nData];
  
  double**pars = getPars();
  
  double parVals[200]; // 200 --> should use NEXTRAPARAMMAX_C from dimensions.h, and synchronize NEXTRAPARAMMAX_C with MNE from endmini.inc

  for (int i = 0; i < npars; i++)
    parVals[i] = *pars[i];

  double covmat[npars * npars];
  fill(covmat,covmat+npars*npars, 0.);

  //First call to FCN for initialisation
  double chi2;
  myFCN(chi2,parVals,1);
  
  // Least squares minimisation
  ceres::Solver::Options soloptions;
  ceres::Solver::Summary summary;
  
  // Cost function:
  ceres::DynamicNumericDiffCostFunction<CostFunctiorData>* dynamic_cost_function =
    new ceres::DynamicNumericDiffCostFunction<CostFunctiorData>(new CostFunctiorData);
    
  dynamic_cost_function->AddParameterBlock(npars);

  int nres = (chi2options_.chi2poissoncorr) ? cndatapoints_.npoints + systema_.nsys + cndatapoints_.npoints : cndatapoints_.npoints + systema_.nsys;

  dynamic_cost_function->SetNumResiduals(nres);

  auto diffCostFunction = new CostFuntionrDataAndDerivative();
  //  diffCostFunction->AddParameterBlock(npars);
  // diffCostFunction->SetNumResiduals(nres);
  
  // Loss function to compensate offset for the log penalty terms
  ceres::LossFunction* loss_function(new PenaltyLog);

  ceres::Problem problem;
  problem.AddResidualBlock(dynamic_cost_function, loss_function, parVals);

  // Set CERES options
  soloptions.minimizer_progress_to_stdout = true;

  soloptions.function_tolerance = ceresNode["tolerance"].as<double>(); // typical chi2 is ~1000

  int strategy_type = ceresNode["strategy"].as<int>();
  switch (strategy_type)
    {
    case 0:
      soloptions.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
      break;
    case 1:
      soloptions.trust_region_strategy_type = ceres::DOGLEG;
      soloptions.dogleg_type = ceres::SUBSPACE_DOGLEG;
      break;
    default:
      soloptions.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    }

  // Additional options (to be checked if any of this is actually usable/usefull)
  //soloptions.logging_type = ceres::SILENT;
  
  //multithreading
  //soloptions.num_threads = 4;

  //Ceres options
  //soloptions.max_num_iterations = 100000;

  //soloptions.minimizer_type = ceres::TRUST_REGION; //ceres::LINE_SEARCH;
  //soloptions.linear_solver_type = ceres::DENSE_QR; //ceres::SPARSE_NORMAL_CHOLESKY;
  //soloptions.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT; //ceres::DOGLEG;
  //soloptions.dogleg_type = ceres::SUBSPACE_DOGLEG;
  //soloptions.use_inner_iterations = true;
  //soloptions.use_nonmonotonic_steps = true;
  
  //soloptions.max_num_consecutive_invalid_steps = 10;

  //soloptions.dense_linear_algebra_library_type = ceres::EIGEN; //ceres::LAPACK
  
  ceres::Solve(soloptions, &problem, &summary);

  cout << std::endl;
  cout << "CERES minimisation has converged" << std::endl;

  int docov = ceresNode["covariance"].as<int>();
  if (docov)
    {
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
    }

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
      sprintf(val, "%.6f", parVals[i]);
      sprintf(err, "%.6f", sqrt(covmat[i*npars+i]));
      std::cout << "  " << _allParameterNames[i] << ": [ " << val << ", " << err << " ]" << std::endl;
    }
  cout << "----- End of parameters in YAML format" << endl; 
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

  if (chi2options_.chi2poissoncorr)
    delete[] offset;

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

  f.close();

  f.open(stringFromFortran(coutdirname_.outdirname,sizeof(coutdirname_.outdirname))+"/pars.yaml");
  f << "Parameters:" << endl;
  for (int i = 0; i < npars; i++)
    {
      char val[15];
      char err[15];
      sprintf(val, "%.6f", parVals[i]);
      sprintf(err, "%.6f", sqrt(covmat[i*npars+i]));
      f << "  " << _allParameterNames[i] << ": [ " << val << ", " << err << " ]" << endl;
    }
  f.close();
  
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

  f.close();
}

}

