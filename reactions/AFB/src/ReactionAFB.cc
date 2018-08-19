 
/*
   @file ReactionAFB.cc
   @date 2018-07-16
   @author  AddReaction.py
   Created by  AddReaction.py on 2018-07-16
*/

#include "ReactionAFB.h"
#include "iostream"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
//#include "/home/jf6g13/LHAPDF-6.2.1/include/LHAPDF_link/LHAPDF.h"

using namespace std;

// Constants
#define PI 3.14159265
#define GeVtofb 0.38937966e+12

// SM parameters
#define MZ 91.19
#define GammaZ 2.5
#define alphaEM 0.0078125
#define sthetaW 0.4817

// set collider energy and luminosity
#define energy 13000
#define lum 300

// define acceptance cuts
#define eta_cut 2.5
#define pT_cut 20

// kinematic cuts
#define ycut 1

// // PDF set and grid
// #define setname "CT14nnlo"
// #define PDF_set 0

// Setting of the integration
const int dim_integration = 2; // Only integration on yreduced variable
// Integration extremes
double yreducedmin = 0;
double yreducedmax = ycut;

// Integration number of calls 
size_t calls = 10000;

int   ichk = 0;
int   iset = 1;


double *ReactionAFB::propagators (double Minv)
{
    const double e = sqrt(4*PI*alphaEM);
    const double gsm = (e/(sthetaW*sqrt(1-pow(sthetaW,2))))*sqrt(1+pow(sthetaW,4));
    
    const double smangle = atan(-pow(sthetaW,2));
    
    // SM couplings
    static double *couplings_foton = new double[8];
    couplings_foton[0] = e*(2.0/3.0);
    couplings_foton[1] = 0;
    couplings_foton[2] = e*(-1.0/3.0);
    couplings_foton[3] = 0;
    couplings_foton[4] = e*(-1.0);
    couplings_foton[5] = 0;
    couplings_foton[6] = 0;
    couplings_foton[7] = 0;
    
    static double *couplings_Z = new double[8];
    couplings_Z[0] = (1.0/2.0)*gsm*(1.0/6.0)*(3*cos(smangle)+8*sin(smangle));
    couplings_Z[1] = (1.0/2.0)*gsm*(cos(smangle)/2.0);
    couplings_Z[2] = (1.0/2.0)*gsm*(1.0/6.0)*(-3*cos(smangle)-4*sin(smangle));
    couplings_Z[3] = (1.0/2.0)*gsm*(-cos(smangle)/2.0);
    couplings_Z[4] = (1.0/2.0)*gsm*((-cos(smangle)/2.0)+(-2*sin(smangle)));
    couplings_Z[5] = (1.0/2.0)*gsm*(-cos(smangle)/2.0);
    couplings_Z[6] = (1.0/2.0)*gsm*(cos(smangle)/2.0);
    couplings_Z[7] = (1.0/2.0)*gsm*(cos(smangle)/2.0);
        
    // Even combination of couplings
    double even_foton_up = (pow(couplings_foton[0],2)+pow(couplings_foton[1],2))*(pow(couplings_foton[4],2)+pow(couplings_foton[5],2));
    double even_foton_down = (pow(couplings_foton[2],2)+pow(couplings_foton[3],2))*(pow(couplings_foton[4],2)+pow(couplings_foton[5],2));
    double even_interf_up = ((couplings_foton[0]*couplings_Z[0])+(couplings_foton[1]*couplings_Z[1]))*((couplings_foton[4]*couplings_Z[4])+(couplings_foton[5]*couplings_Z[5]));
    double even_interf_down = ((couplings_foton[2]*couplings_Z[2])+(couplings_foton[3]*couplings_Z[3]))*((couplings_foton[4]*couplings_Z[4])+(couplings_foton[5]*couplings_Z[5]));
    double even_Z_up = (pow(couplings_Z[0],2)+pow(couplings_Z[1],2))*(pow(couplings_Z[4],2)+pow(couplings_Z[5],2));
    double even_Z_down = (pow(couplings_Z[2],2)+pow(couplings_Z[3],2))*(pow(couplings_Z[4],2)+pow(couplings_Z[5],2));
    
    // Odd combination of couplings
    double odd_foton_up = 4*couplings_foton[0]*couplings_foton[1]*couplings_foton[4]*couplings_foton[5];
    double odd_foton_down = 4*couplings_foton[2]*couplings_foton[3]*couplings_foton[4]*couplings_foton[5];
    double odd_interf_up = (couplings_foton[0]*couplings_Z[1]+couplings_foton[1]*couplings_Z[0])*(couplings_foton[4]*couplings_Z[5]+couplings_foton[5]*couplings_Z[4]);
    double odd_interf_down = (couplings_foton[2]*couplings_Z[3]+couplings_foton[3]*couplings_Z[2])*(couplings_foton[4]*couplings_Z[5]+couplings_foton[5]*couplings_Z[4]);
    double odd_Z_up = 4*couplings_Z[0]*couplings_Z[1]*couplings_Z[4]*couplings_Z[5];
    double odd_Z_down = 4*couplings_Z[2]*couplings_Z[3]*couplings_Z[4]*couplings_Z[5];
    
    // Propagators squared and interference
    double foton_squared = 1.0/pow(Minv,4);
    double interference = 2.0*(-pow(Minv,2)*(pow(MZ,2)-pow(Minv,2)))/(pow(Minv,4)*((pow(pow(MZ,2)-pow(Minv,2),2))+pow(MZ,2)*pow(GammaZ,2)));
    double Z_squared = 1.0/(pow(pow(MZ,2)-pow(Minv,2),2)+pow(MZ,2)*pow(GammaZ,2));
    
    static double *propagators = new double[4];
    
    propagators[0] = (even_foton_up * foton_squared)+(even_interf_up * interference) + (even_Z_up * Z_squared);
    propagators[1] = (odd_foton_up * foton_squared)+(odd_interf_up * interference) + (odd_Z_up * Z_squared);
    propagators[2] = (even_foton_down * foton_squared)+(even_interf_down * interference) + (even_Z_down * Z_squared);
    propagators[3] = (odd_foton_down * foton_squared)+(odd_interf_down * interference) + (odd_Z_down * Z_squared);
        
    return propagators;
}

// extern "C"{
// void allfxq_(int iset, double x, double q2, double *pdfs, int n, int ichk);
// }

// extern "C"{
// void allfxq_(int iset, double x, double q2, double *pdfs, int n, int ichk);
// }

extern "C"{
double fvalxq_(int iset, int id, double x, double q2, int ichk);
}

// extern "C"{
// void fpdfxq_(int *iset, const double *x, const double *q2, double *pdfs, int *ichk);
// }

////UUBAR EVEN FORWARD Matrix element
double ReactionAFB::uubarEF_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
//     double f1u = (pdf[PDF_set]->xfxQ(2, x1, Q))/x1;
//     double f1c = (pdf[PDF_set]->xfxQ(4, x1, Q))/x1;
//     double f2ubar = (pdf[PDF_set]->xfxQ(-2, x2, Q))/x2;
//     double f2cbar = (pdf[PDF_set]->xfxQ(-4, x2, Q))/x2;
    

    static double pdfs1[13];
    static double pdfs2[13];
    double q2 = pow(Q,2);
    
//     allfxq_(&iset,&x1,&q2,pdfs1,&ichk);
//     allfxq_(&iset,&x2,&q2,pdfs2,&ichk);
    
//     allfxq_(iset,x1,q2,pdfs1,0,ichk);
//     allfxq_(iset,x2,q2,pdfs2,0,ichk);
//     
// //     fpdfxq_(&iset,&x1,&q2,pdfs1,&ichk);
// //     fpdfxq_(&iset,&x2,&q2,pdfs2,&ichk);
//     
//     double f1u = pdfs1[7]/x1;
//     double f1c = pdfs1[9]/x1;
//     double f2ubar = pdfs2[4]/x2;
//     double f2cbar = pdfs2[2]/x2;
    
    double f1u = fvalxq_(iset,2,x1,q2,ichk)/x1;
    double f1c = fvalxq_(iset,4,x1,q2,ichk)/x1;
    double f2ubar = fvalxq_(iset,-2,x2,q2,ichk)/x2;
    double f2cbar = fvalxq_(iset,-4,x2,q2,ichk)/x2;
    
    
//     double f1u = xfx(x1,Q,2)/x1;
//     double f1c = xfx(x1,Q,4)/x1;
//     double f2ubar = xfx(x2,Q,-2)/x2;
//     double f2cbar = xfx(x2,Q,-4)/x2;

    // PDF combinations    
    double uubar_PDF = f1u*f2ubar + f1c*f2cbar;
    
    // Angular integration limits
    double qqbar_cos_theta_max = min(cos(2*atan(exp(-eta_cut-y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
    double qqbar_cos_theta_min = 0;
 
    double angular_integration_EF = (qqbar_cos_theta_max-qqbar_cos_theta_min)+(1.0/3.0)*(pow(qqbar_cos_theta_max,3)-pow(qqbar_cos_theta_min,3));
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_EF = dsigma*angular_integration_EF;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // U-UBAR
    double uubarEF = uubar_PDF*dsigma_EF;
    
    
  
    double *propagator = propagators (Minv);
    
    return uubarEF * propagator[0];
}

// double (ReactionAFB::*pointer_to_funct)(double *entries, size_t dim, void *params) = &ReactionAFB::uubarEF_funct;


////UUBAR EVEN FORWARD Integration
double ReactionAFB::integration_uubarEF (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};
    double uubarEF, error_uubarEF;
    
    // Initialization of the integration (quite a black box)
    
    gsl_monte_function Integrate_uubarEF = { &(ReactionAFB::uubarEF_funct), dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);   
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_uubarEF, integration_inf, integration_sup, dim_integration, calls, r, s, &uubarEF, &error_uubarEF);
    do
    {
        gsl_monte_vegas_integrate (&Integrate_uubarEF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &uubarEF, &error_uubarEF);
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);
    
    return 2*uubarEF;
}



double *ReactionAFB::observables (double Minv_inf, double Minv_sup)
{
    double uubarEF = integration_uubarEF (Minv_inf, Minv_sup);
//     double uubarEB = integration_uubarEB (Minv_inf, Minv_sup);
//     double uubarOF = integration_uubarOF (Minv_inf, Minv_sup);
//     double uubarOB = integration_uubarOB (Minv_inf, Minv_sup);
//     
//     double ubaruEF = integration_ubaruEF (Minv_inf, Minv_sup);
//     double ubaruEB = integration_ubaruEB (Minv_inf, Minv_sup);
//     double ubaruOF = integration_ubaruOF (Minv_inf, Minv_sup);
//     double ubaruOB = integration_ubaruOB (Minv_inf, Minv_sup);
//     
//     double ddbarEF = integration_ddbarEF (Minv_inf, Minv_sup);
//     double ddbarEB = integration_ddbarEB (Minv_inf, Minv_sup);
//     double ddbarOF = integration_ddbarOF (Minv_inf, Minv_sup);
//     double ddbarOB = integration_ddbarOB (Minv_inf, Minv_sup);
//         
//     double dbardEF = integration_dbardEF (Minv_inf, Minv_sup); 
//     double dbardEB = integration_dbardEB (Minv_inf, Minv_sup);
//     double dbardOF = integration_dbardOF (Minv_inf, Minv_sup);
//     double dbardOB = integration_dbardOB (Minv_inf, Minv_sup);
//     
//     // Reconstructed Forward and Backward
//     double Forward = uubarEF+ubaruEF+ddbarEF+dbardEF+uubarOF+ubaruOF+ddbarOF+dbardOF;
//     double Backward = uubarEB+ubaruEB+ddbarEB+dbardEB+uubarOB+ubaruOB+ddbarOB+dbardOB;
// 
//     // Cross section and statistical relative uncertainty
//     double XS = Forward + Backward;
//     double Epsilon_XS = 100.0 / sqrt(XS*lum);
//     
//     
//     // Reconstructed AFB and statistical absolute uncertainty
//     double AFB = (Forward - Backward) / (Forward + Backward);
//     double Delta_AFB = sqrt((1.0-pow(AFB,2)) / (XS*lum));
//     
//     
//     double *results = new double[4];
//     
//     results[0] = XS;
//     results[1] = Epsilon_XS;
//     results[2] = AFB;
//     results[3] = Delta_AFB;
    
    double *results = new double[4];
    results[0] = uubarEF;
    results[1] = uubarEF;
    results[2] = uubarEF;
    results[3] = uubarEF;
    
    
    return results;
}



// the class factories
extern "C" ReactionAFB* create() {
  return new ReactionAFB();
}

// Initialize at the start of the computation
int ReactionAFB::initAtStart(const string &s)
{
  return 0;
}

// Main function to compute results at an iteration
int ReactionAFB::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
    auto *Minv_min  = GetBinValues(dataSetID,"Minv_min"), *Minv_max  = GetBinValues(dataSetID,"Minv_max");  
    if (Minv_min == nullptr || Minv_max == nullptr) {
    std::cout << "\n\nFATAL ERROR" << std::endl;
    std::cout << "CHECK THE DATAFILE !!!" << std::endl;
    return 1;
  }

    auto min = *Minv_min, max = *Minv_max;
    
    printf("%1f \n", min[0]);
    printf("%1f \n", max[0]);
    
    printf("Calling AFB \n");
    
//     val[0] = *propagators (min[0]);
    
//     val[0] = integration_uubarEF(min[0], max[0]);
    
    double *results = observables (min[0], max[0]);
    
    val[0] = results[0];

    
  return 0;
}
