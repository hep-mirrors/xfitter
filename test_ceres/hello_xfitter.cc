// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: keir@google.com (Keir Mierle)
//
// A simple example of using the Ceres minimizer.
//
// Minimize 0.5 (10 - x)^2 using jacobian matrix computed using
// automatic differentiation.

#include "ceres/ceres.h"
#include "glog/logging.h"
#include "ceres/dynamic_numeric_diff_cost_function.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::DynamicNumericDiffCostFunction;
using ceres::EvaluationCallback;


extern "C" {
  void hfbanner_();
  void xfitter_();
  void read_steer_();
  void read_data_();
  void init_theory_modules_();
  void do_fit_();
  void close_theor_eval_();
  void hf_errsum_(const int& io);
  void init_func_map_();
  int read_reactions_();
  void parse_params_();
  void init_rnd_seeds_();
  void init_func_map_();
  // FCN
  void fcn_(const int& npar, const double& dummy, double& chi2out, const double* pars, const int& iflag, const double& dummy2);
}

// some more connections

void myFCN(const int& npar,double& chi2out, const double* pars, const int& iflag) {
  // move parameters around:
  double parsXF[200];
  for (int i=0; i<200; i++) {
    parsXF[i] = 0;
  }
  parsXF[2-1] = pars[0];
  parsXF[3-1] = pars[1];
  parsXF[7-1] = pars[2];
  parsXF[8-1] = pars[3];
  parsXF[9-1] = 25.0;
  parsXF[12-1] = pars[4];
  parsXF[13-1] = pars[5];
  parsXF[15-1] = pars[6];
  parsXF[22-1] = pars[7];
  parsXF[23-1] = pars[8];
  parsXF[33-1] = pars[9];
  parsXF[34-1] = pars[10];
  parsXF[41-1] = pars[11];
  parsXF[42-1] = pars[12];
  parsXF[43-1] = pars[13];
  parsXF[101] = 0.118;


  fcn_(npar, 0, chi2out, parsXF, iflag, 0);
};

void init_pars() {
  init_rnd_seeds_();
  parse_params_();
  int i = read_reactions_();
}

void init_theo() {
  init_func_map_();
  init_theory_modules_();
}

// Full connection
void startXF() {
  hfbanner_();
  read_steer_();
  init_pars();
  read_data_();
  init_theory_modules_();
}


#ifndef NSYSMAX_C
#define NSYSMAX_C 750
#endif
#ifndef NTOT_C
#define NTOT_C 2700
#endif
#ifndef NCHI2POINTS_C
#define NCHI2POINTS_C 100
#endif
#ifndef NEXTRAPARAMMAX
#define NEXTRAPARAMMAX_C 150
#endif
#ifndef NSET_C
#define NSET_C 150
#endif

extern struct {
  double theo_[NTOT_C];          // Theory predictions, filled for each iteration
  double theo_mod_[NTOT_C];      // Theory predictions, filled for each iteration
  double theo_fix_[NTOT_C];      // Fixed theory prediction (if given by &InTheory namelist)
  double theo_unc_[NTOT_C];      // Uncorrelated uncertainty on theory predictions 
  double theo_tot_up_[NTOT_C];         // Total up uncertainty on theory predictions 
  double theo_tot_down_[NTOT_C]; // Total down uncertainty on theory predictions
} c_theo_;


extern struct {
  double alpha[NTOT_C];       // Total uncorrelated errors
  double alpha_mod[NTOT_C];   // Total uncorrelated errors modified
  double beta[NTOT_C][NSYSMAX_C];   // Influence of systematic errors on measurements
  double sysa[NSYSMAX_C][NSYSMAX_C];    // Correlation matrix of systematics
  char system[NSYSMAX_C][64];     // Names of correlated systematic errors
  int nsys;                 // Actual number of correlated systematic sources
} systema_;

extern struct {
  int npoints_;       // Actual number of data points
} cndatapoints_;

extern struct {
  double residuals_[NTOT_C];
} c_resid_;



// get data
const int NPars = 14;

static int counter = 0;

struct CostFunctiorData
{
  bool operator()(double const* const* parameters, double* residuals) const
  {
    counter++;
    double chi2 = 0;
    myFCN(NPars, chi2, parameters[0], 2);
    std::cout << " Call again " << counter << " chi2=" << chi2;
    chi2 = 0;
    for (int i = 0; i< cndatapoints_.npoints_  + systema_.nsys; i++) {
      residuals[i] = c_resid_.residuals_[i];
      chi2 += residuals[i]*residuals[i];
    }
    std::cout << " residuals squared: " <<  chi2 << std::endl;
    return true;
  }
};

struct CostFunctorSyst
{
  bool operator()(double const* const* parameters, double* residuals) const
  {
    return true;
  }
};

struct MyDynamicNumericDiffCostFunction : DynamicNumericDiffCostFunction<CostFunctiorData>, EvaluationCallback {
  virtual void PrepareForEvaluation(bool evaluate_jacobians, bool new_evaluation_point) {
  }
};


int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);


  startXF();
  int nData = cndatapoints_.npoints_;
  int nSyst = systema_.nsys;
  
  // Cost function:
  DynamicNumericDiffCostFunction<CostFunctiorData>* dynamic_cost_function =
    new DynamicNumericDiffCostFunction<CostFunctiorData>(new CostFunctiorData);
  
  dynamic_cost_function->AddParameterBlock(NPars);  // For PDFs ... Can add another block for extra pars
  double pars[] = {-0.061953,5.562367,0.166118,-0.383100, 0.810476, 4.823512, 9.921366, 1.029995, 4.846279,
  		   7.059694,1.548098, 0.268798, -0.127297, 9.586246, 0., 0., 0.,  0., 0.};
  //  double pars[] = {-0.53061953,5.562367,0.,-0.383100, 0.810476, 4.823512, 9.921366, 1.029995, 4.846279,
  //		   17.059694,1.548098, 0.268798, 0.127297, 0.586246};

  // try to call once:
  double chi2 = 0;
  myFCN(NPars, chi2, pars, 1);
  std::cout << "FCN1: " << chi2 << std::endl;
  
  dynamic_cost_function->SetNumResiduals(nData+nSyst);  // For data
 
  
  Problem myProblem;
  myProblem.AddResidualBlock(dynamic_cost_function, NULL, pars);

  Solver::Options myOptions;
  myOptions.minimizer_progress_to_stdout = true;
  // myOptions.trust_region_strategy_type = ceres::DOGLEG;
  // myOptions.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
  
  Solver::Summary mySummary;
  Solve(myOptions, &myProblem, &mySummary); 

  std::cout << mySummary.FullReport() << "\n";
  for (int i =0; i<NPars; i++){
    std::cout << "Parameter " << i << " value " << pars[i] << std::endl;
  }

  /*
  std::cout << " \n \n \n Get Covariance \n \n \n";

  ceres::Covariance::Options covOptions;

  covOptions.algorithm_type = ceres::DENSE_SVD;
  
  ceres::Covariance covar(covOptions);
  std::vector<std::pair<const double*, const double*> > covar_blocks;
  covar_blocks.push_back(std::make_pair(pars,pars));

  CHECK(covar.Compute(covar_blocks,&myProblem));

  
  double cov_result[14*14];
  covar.GetCovarianceBlock(pars, pars, cov_result);
  for (int i=0; i<14; i++) {
    for (int j=0; j<14; j++) {
      std::cout << " cov i j " << i << " " <<j << " "<< cov_result[i+j*14] << std::endl;
    }
  }

  std::cout << "\n\n\n";
  for (int i =0; i<NPars; i++){
    std::cout << "Parameter " << i << " value " << pars[i] << " +- " << sqrt(cov_result[i+i*14]) << std::endl;
  }
  */
  return 0;
}
