 
/*
   @file ReactionAFB.cc
   @date 2018-07-16
   @author  AddReaction.py
   Created by  AddReaction.py on 2018-07-16
*/

#include "ReactionAFB.h"
#include "iostream"
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

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

// PDF set and check
const int   iset = 1;
const int   ichk = 0;

// Setting of the integration
const int dim_integration = 2; // Integration on yreduced and Minv
// Integration extremes
double yreducedmin = 0;
double yreducedmax = ycut;

// Integration number of calls 
size_t calls = 10000;




//// Function retirning the combination of propagators
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
    
    
    
    // Call the PDFs
    double f1u = (ReactionTheory::xfx(x1,Q,2)) / x1;
    double f1c = (ReactionTheory::xfx(x1,Q,4)) / x1;
    double f2ubar = (ReactionTheory::xfx(x2,Q,-2)) / x2;
    double f2cbar = (ReactionTheory::xfx(x2,Q,-4)) / x2;

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
    
    return uubarEF * propagator[0]; // Multiply the PDFs combination with the correct propagator.
}

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
    
    return 2*uubarEF; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
}


// ////UUBAR EVEN BACKWARD Matrix element
// double ReactionAFB::uubarEB_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1u = (ReactionTheory::xfx(x1,Q,2)) / x1;
//     double f1c = (ReactionTheory::xfx(x1,Q,4)) / x1;
//     double f2ubar = (ReactionTheory::xfx(x2,Q,-2)) / x2;
//     double f2cbar = (ReactionTheory::xfx(x2,Q,-4)) / x2;
// 
//     // PDF combinations    
//     double uubar_PDF = f1u*f2ubar + f1c*f2cbar;
//     
//     // Angular integration limits
//     double qbarq_cos_theta_max = 0;
//     double qbarq_cos_theta_min = max(cos(PI - 2*atan(exp(-eta_cut-y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));   
//  
//     double angular_integration_EB = (qbarq_cos_theta_max-qbarq_cos_theta_min)+(1.0/3.0)*(pow(qbarq_cos_theta_max,3)-pow(qbarq_cos_theta_min,3));
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_EB = dsigma*angular_integration_EB;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // U-UBAR
//     double uubarEB = uubar_PDF*dsigma_EB;
//   
//     double *propagator = propagators (Minv);
//     
//     return uubarEB * propagator[0]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////UUBAR EVEN BACKWARD Integration
// double ReactionAFB::integration_uubarEB (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//     
//     double uubarEB, error_uubarEB;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_uubarEB = { &(ReactionAFB::uubarEB_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_uubarEB, integration_inf, integration_sup, dim_integration, calls, r, s, &uubarEB, &error_uubarEB);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_uubarEB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &uubarEB, &error_uubarEB);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*uubarEB; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////UUBAR ODD FORWARD Matrix element
// double ReactionAFB::uubarOF_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1u = (ReactionTheory::xfx(x1,Q,2)) / x1;
//     double f1c = (ReactionTheory::xfx(x1,Q,4)) / x1;
//     double f2ubar = (ReactionTheory::xfx(x2,Q,-2)) / x2;
//     double f2cbar = (ReactionTheory::xfx(x2,Q,-4)) / x2;
// 
//     // PDF combinations    
//     double uubar_PDF = f1u*f2ubar + f1c*f2cbar;
//     
//     // Angular integration limits
//     double qqbar_cos_theta_max = min(cos(2*atan(exp(-eta_cut-y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
//     double qqbar_cos_theta_min = 0;  
//  
//     double angular_integration_OF = pow(qqbar_cos_theta_max,2) - pow(qqbar_cos_theta_min,2);
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_OF = dsigma*angular_integration_OF;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // U-UBAR
//     double uubarOF = uubar_PDF*dsigma_OF;
//   
//     double *propagator = propagators (Minv);
//     
//     return uubarOF * propagator[1]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////UUBAR ODD FORWARD Integration
// double ReactionAFB::integration_uubarOF (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//        
//     double uubarOF, error_uubarOF;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_uubarOF = { &(ReactionAFB::uubarOF_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_uubarOF, integration_inf, integration_sup, dim_integration, calls, r, s, &uubarOF, &error_uubarOF);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_uubarOF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &uubarOF, &error_uubarOF);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*uubarOF; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////UUBAR ODD BACKWARD Matrix element
// double ReactionAFB::uubarOB_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1u = (ReactionTheory::xfx(x1,Q,2)) / x1;
//     double f1c = (ReactionTheory::xfx(x1,Q,4)) / x1;
//     double f2ubar = (ReactionTheory::xfx(x2,Q,-2)) / x2;
//     double f2cbar = (ReactionTheory::xfx(x2,Q,-4)) / x2;
// 
//     // PDF combinations    
//     double uubar_PDF = f1u*f2ubar + f1c*f2cbar;
//     
//     // Angular integration limits
//     double qbarq_cos_theta_max = 0;
//     double qbarq_cos_theta_min = max(cos(PI - 2*atan(exp(-eta_cut-y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2)))); 
//  
//     double angular_integration_OB = pow(qbarq_cos_theta_max,2) - pow(qbarq_cos_theta_min,2);
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_OB = dsigma*angular_integration_OB;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // U-UBAR
//     double uubarOB = uubar_PDF*dsigma_OB;
//   
//     double *propagator = propagators (Minv);
//     
//     return uubarOB * propagator[1]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////UUBAR ODD BACKWARD Integration
// double ReactionAFB::integration_uubarOB (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//     
//     double uubarOB, error_uubarOB;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_uubarOB = { &(ReactionAFB::uubarOB_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_uubarOB, integration_inf, integration_sup, dim_integration, calls, r, s, &uubarOB, &error_uubarOB);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_uubarOB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &uubarOB, &error_uubarOB);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*uubarOB; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////UBARU EVEN FORWARD Matrix element
// double ReactionAFB::ubaruEF_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1ubar = (ReactionTheory::xfx(x1,Q,-2)) / x1;
//     double f1cbar = (ReactionTheory::xfx(x1,Q,-4)) / x1;
//     double f2u = (ReactionTheory::xfx(x2,Q,2)) / x2;
//     double f2c = (ReactionTheory::xfx(x2,Q,4)) / x2;
// 
//     // PDF combinations    
//     double ubaru_PDF = f1ubar*f2u + f1cbar*f2c;
//     
//     // Angular integration limits
//     double qbarq_cos_theta_max = 0;
//     double qbarq_cos_theta_min = max(cos(PI - 2*atan(exp(-eta_cut-y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));   
//  
//     double angular_integration_EB = (qbarq_cos_theta_max-qbarq_cos_theta_min)+(1.0/3.0)*(pow(qbarq_cos_theta_max,3)-pow(qbarq_cos_theta_min,3));
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_EB = dsigma*angular_integration_EB;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // UBAR-U
//     double ubaruEF = ubaru_PDF*dsigma_EB;
//   
//     double *propagator = propagators (Minv);
//     
//     return ubaruEF * propagator[0]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////UBARU EVEN FORWARD Integration
// double ReactionAFB::integration_ubaruEF (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//     
//     double ubaruEF, error_ubaruEF;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_ubaruEF = { &(ReactionAFB::ubaruEF_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_ubaruEF, integration_inf, integration_sup, dim_integration, calls, r, s, &ubaruEF, &error_ubaruEF);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_ubaruEF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ubaruEF, &error_ubaruEF);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*ubaruEF; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////UBARU EVEN BACKWARD Matrix element
// double ReactionAFB::ubaruEB_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1ubar = (ReactionTheory::xfx(x1,Q,-2)) / x1;
//     double f1cbar = (ReactionTheory::xfx(x1,Q,-4)) / x1;
//     double f2u = (ReactionTheory::xfx(x2,Q,2)) / x2;
//     double f2c = (ReactionTheory::xfx(x2,Q,4)) / x2;
// 
//     // PDF combinations    
//     double ubaru_PDF = f1ubar*f2u + f1cbar*f2c;
//     
//     // Angular integration limits
//     double qqbar_cos_theta_max = min(cos(2*atan(exp(-eta_cut-y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
//     double qqbar_cos_theta_min = 0;  
//  
//     double angular_integration_EF = (qqbar_cos_theta_max-qqbar_cos_theta_min)+(1.0/3.0)*(pow(qqbar_cos_theta_max,3)-pow(qqbar_cos_theta_min,3));
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_EF = dsigma*angular_integration_EF;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // UBAR-U
//     double ubaruEB = ubaru_PDF*dsigma_EF;
//   
//     double *propagator = propagators (Minv);
//     
//     return ubaruEB * propagator[0]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////UUBAR EVEN BACKWARD Integration
// double ReactionAFB::integration_ubaruEB (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//     
//     double ubaruEB, error_ubaruEB;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_ubaruEB = { &(ReactionAFB::ubaruEB_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_ubaruEB, integration_inf, integration_sup, dim_integration, calls, r, s, &ubaruEB, &error_ubaruEB);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_ubaruEB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ubaruEB, &error_ubaruEB);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*ubaruEB; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////UBARU ODD FORWARD Matrix element
// double ReactionAFB::ubaruOF_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1ubar = (ReactionTheory::xfx(x1,Q,-2)) / x1;
//     double f1cbar = (ReactionTheory::xfx(x1,Q,-4)) / x1;
//     double f2u = (ReactionTheory::xfx(x2,Q,2)) / x2;
//     double f2c = (ReactionTheory::xfx(x2,Q,4)) / x2;
// 
//     // PDF combinations    
//     double ubaru_PDF = f1ubar*f2u + f1cbar*f2c;
//     
//     // Angular integration limits
//     double qbarq_cos_theta_max = 0;
//     double qbarq_cos_theta_min = max(cos(PI - 2*atan(exp(-eta_cut-y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2)))); 
//  
//     double angular_integration_OB = pow(qbarq_cos_theta_max,2) - pow(qbarq_cos_theta_min,2);
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_OB = dsigma*angular_integration_OB;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // UBAR-U
//     double ubaruOF = ubaru_PDF*dsigma_OB;
//   
//     double *propagator = propagators (Minv);
//     
//     return ubaruOF * propagator[1]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////UBARU ODD FORWARD Integration
// double ReactionAFB::integration_ubaruOF (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//        
//     double ubaruOF, error_ubaruOF;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_ubaruOF = { &(ReactionAFB::ubaruOF_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_ubaruOF, integration_inf, integration_sup, dim_integration, calls, r, s, &ubaruOF, &error_ubaruOF);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_ubaruOF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ubaruOF, &error_ubaruOF);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*ubaruOF; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////UBARU ODD BACKWARD Matrix element
// double ReactionAFB::ubaruOB_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1ubar = (ReactionTheory::xfx(x1,Q,-2)) / x1;
//     double f1cbar = (ReactionTheory::xfx(x1,Q,-4)) / x1;
//     double f2u = (ReactionTheory::xfx(x2,Q,2)) / x2;
//     double f2c = (ReactionTheory::xfx(x2,Q,4)) / x2;
// 
//     // PDF combinations    
//     double ubaru_PDF = f1ubar*f2u + f1cbar*f2c;
//     
//     // Angular integration limits
//     double qqbar_cos_theta_max = min(cos(2*atan(exp(-eta_cut-y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
//     double qqbar_cos_theta_min = 0; 
//  
//     double angular_integration_OF = pow(qqbar_cos_theta_max,2) - pow(qqbar_cos_theta_min,2);
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_OF = dsigma*angular_integration_OF;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // UBAR-U
//     double ubaruOB = ubaru_PDF*dsigma_OF;
//   
//     double *propagator = propagators (Minv);
//     
//     return ubaruOB * propagator[1]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////UBARU ODD BACKWARD Integration
// double ReactionAFB::integration_ubaruOB (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//     
//     double ubaruOB, error_ubaruOB;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_ubaruOB = { &(ReactionAFB::ubaruOB_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_ubaruOB, integration_inf, integration_sup, dim_integration, calls, r, s, &ubaruOB, &error_ubaruOB);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_ubaruOB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ubaruOB, &error_ubaruOB);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*ubaruOB; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////DDBAR EVEN FORWARD Matrix element
// double ReactionAFB::ddbarEF_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1d = (ReactionTheory::xfx(x1,Q,1)) / x1;
//     double f1s = (ReactionTheory::xfx(x1,Q,3)) / x1;
//     double f1b = (ReactionTheory::xfx(x1,Q,5)) / x1;
//     double f2dbar = (ReactionTheory::xfx(x2,Q,-1)) / x2;
//     double f2sbar = (ReactionTheory::xfx(x2,Q,-3)) / x2;
//     double f2bbar = (ReactionTheory::xfx(x2,Q,-5)) / x2;
// 
//     // PDF combinations    
//     double ddbar_PDF = f1d*f2dbar + f1s*f2sbar + f1b*f2bbar;
//     
//     // Angular integration limits
//     double qqbar_cos_theta_max = min(cos(2*atan(exp(-eta_cut-y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
//     double qqbar_cos_theta_min = 0;
//  
//     double angular_integration_EF = (qqbar_cos_theta_max-qqbar_cos_theta_min)+(1.0/3.0)*(pow(qqbar_cos_theta_max,3)-pow(qqbar_cos_theta_min,3));
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_EF = dsigma*angular_integration_EF;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // D-DBAR
//     double ddbarEF = ddbar_PDF*dsigma_EF;
//   
//     double *propagator = propagators (Minv);
//     
//     return ddbarEF * propagator[2]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////DDBAR EVEN FORWARD Integration
// double ReactionAFB::integration_ddbarEF (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//     
//     double ddbarEF, error_ddbarEF;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_ddbarEF = { &(ReactionAFB::ddbarEF_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_ddbarEF, integration_inf, integration_sup, dim_integration, calls, r, s, &ddbarEF, &error_ddbarEF);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_ddbarEF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ddbarEF, &error_ddbarEF);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*ddbarEF; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////DDBAR EVEN BACKWARD Matrix element
// double ReactionAFB::ddbarEB_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1d = (ReactionTheory::xfx(x1,Q,1)) / x1;
//     double f1s = (ReactionTheory::xfx(x1,Q,3)) / x1;
//     double f1b = (ReactionTheory::xfx(x1,Q,5)) / x1;
//     double f2dbar = (ReactionTheory::xfx(x2,Q,-1)) / x2;
//     double f2sbar = (ReactionTheory::xfx(x2,Q,-3)) / x2;
//     double f2bbar = (ReactionTheory::xfx(x2,Q,-5)) / x2;
// 
//     // PDF combinations    
//     double ddbar_PDF = f1d*f2dbar + f1s*f2sbar + f1b*f2bbar;
//     
//     // Angular integration limits
//     double qbarq_cos_theta_max = 0;
//     double qbarq_cos_theta_min = max(cos(PI - 2*atan(exp(-eta_cut-y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));   
//  
//     double angular_integration_EB = (qbarq_cos_theta_max-qbarq_cos_theta_min)+(1.0/3.0)*(pow(qbarq_cos_theta_max,3)-pow(qbarq_cos_theta_min,3));
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_EB = dsigma*angular_integration_EB;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // D-DBAR
//     double ddbarEB = ddbar_PDF*dsigma_EB;
//   
//     double *propagator = propagators (Minv);
//     
//     return ddbarEB * propagator[2]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////DDBAR EVEN BACKWARD Integration
// double ReactionAFB::integration_ddbarEB (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//     
//     double ddbarEB, error_ddbarEB;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_ddbarEB = { &(ReactionAFB::ddbarEB_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_ddbarEB, integration_inf, integration_sup, dim_integration, calls, r, s, &ddbarEB, &error_ddbarEB);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_ddbarEB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ddbarEB, &error_ddbarEB);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*ddbarEB; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////DDBAR ODD FORWARD Matrix element
// double ReactionAFB::ddbarOF_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1d = (ReactionTheory::xfx(x1,Q,1)) / x1;
//     double f1s = (ReactionTheory::xfx(x1,Q,3)) / x1;
//     double f1b = (ReactionTheory::xfx(x1,Q,5)) / x1;
//     double f2dbar = (ReactionTheory::xfx(x2,Q,-1)) / x2;
//     double f2sbar = (ReactionTheory::xfx(x2,Q,-3)) / x2;
//     double f2bbar = (ReactionTheory::xfx(x2,Q,-5)) / x2;
// 
//     // PDF combinations    
//     double ddbar_PDF = f1d*f2dbar + f1s*f2sbar + f1b*f2bbar;
//     
//     // Angular integration limits
//     double qqbar_cos_theta_max = min(cos(2*atan(exp(-eta_cut-y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
//     double qqbar_cos_theta_min = 0;  
//  
//     double angular_integration_OF = pow(qqbar_cos_theta_max,2) - pow(qqbar_cos_theta_min,2);
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_OF = dsigma*angular_integration_OF;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // D-DBAR
//     double ddbarOF = ddbar_PDF*dsigma_OF;
//   
//     double *propagator = propagators (Minv);
//     
//     return ddbarOF * propagator[3]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////DDBAR ODD FORWARD Integration
// double ReactionAFB::integration_ddbarOF (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//        
//     double ddbarOF, error_ddbarOF;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_ddbarOF = { &(ReactionAFB::ddbarOF_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_ddbarOF, integration_inf, integration_sup, dim_integration, calls, r, s, &ddbarOF, &error_ddbarOF);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_ddbarOF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ddbarOF, &error_ddbarOF);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*ddbarOF; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////DDBAR ODD BACKWARD Matrix element
// double ReactionAFB::ddbarOB_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1d = (ReactionTheory::xfx(x1,Q,1)) / x1;
//     double f1s = (ReactionTheory::xfx(x1,Q,3)) / x1;
//     double f1b = (ReactionTheory::xfx(x1,Q,5)) / x1;
//     double f2dbar = (ReactionTheory::xfx(x2,Q,-1)) / x2;
//     double f2sbar = (ReactionTheory::xfx(x2,Q,-3)) / x2;
//     double f2bbar = (ReactionTheory::xfx(x2,Q,-5)) / x2;
// 
//     // PDF combinations    
//     double ddbar_PDF = f1d*f2dbar + f1s*f2sbar + f1b*f2bbar;
//     
//     // Angular integration limits
//     double qbarq_cos_theta_max = 0;
//     double qbarq_cos_theta_min = max(cos(PI - 2*atan(exp(-eta_cut-y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2)))); 
//  
//     double angular_integration_OB = pow(qbarq_cos_theta_max,2) - pow(qbarq_cos_theta_min,2);
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_OB = dsigma*angular_integration_OB;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // D-DBAR
//     double ddbarOB = ddbar_PDF*dsigma_OB;
//   
//     double *propagator = propagators (Minv);
//     
//     return ddbarOB * propagator[3]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////DDBAR ODD BACKWARD Integration
// double ReactionAFB::integration_ddbarOB (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//     
//     double ddbarOB, error_ddbarOB;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_ddbarOB = { &(ReactionAFB::ddbarOB_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_ddbarOB, integration_inf, integration_sup, dim_integration, calls, r, s, &ddbarOB, &error_ddbarOB);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_ddbarOB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ddbarOB, &error_ddbarOB);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*ddbarOB; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////DBARD EVEN FORWARD Matrix element
// double ReactionAFB::dbardEF_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1dbar = (ReactionTheory::xfx(x1,Q,-1)) / x1;
//     double f1sbar = (ReactionTheory::xfx(x1,Q,-3)) / x1;
//     double f1bbar = (ReactionTheory::xfx(x1,Q,-5)) / x1;
//     double f2d = (ReactionTheory::xfx(x2,Q,1)) / x2;
//     double f2s = (ReactionTheory::xfx(x2,Q,3)) / x2;
//     double f2b = (ReactionTheory::xfx(x2,Q,5)) / x2;
// 
//     // PDF combinations    
//     double dbard_PDF = f1dbar*f2d + f1sbar*f2s + f1bbar*f2b;
//     
//     // Angular integration limits
//     double qbarq_cos_theta_max = 0;
//     double qbarq_cos_theta_min = max(cos(PI - 2*atan(exp(-eta_cut-y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));   
//  
//     double angular_integration_EB = (qbarq_cos_theta_max-qbarq_cos_theta_min)+(1.0/3.0)*(pow(qbarq_cos_theta_max,3)-pow(qbarq_cos_theta_min,3));
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_EB = dsigma*angular_integration_EB;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // DBAR-D
//     double dbardEF = dbard_PDF*dsigma_EB;
//   
//     double *propagator = propagators (Minv);
//     
//     return dbardEF * propagator[2]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////DBARD EVEN FORWARD Integration
// double ReactionAFB::integration_dbardEF (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//     
//     double dbardEF, error_dbardEF;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_dbardEF = { &(ReactionAFB::dbardEF_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_dbardEF, integration_inf, integration_sup, dim_integration, calls, r, s, &dbardEF, &error_dbardEF);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_dbardEF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &dbardEF, &error_dbardEF);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*dbardEF; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////DBARD EVEN BACKWARD Matrix element
// double ReactionAFB::dbardEB_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1dbar = (ReactionTheory::xfx(x1,Q,-1)) / x1;
//     double f1sbar = (ReactionTheory::xfx(x1,Q,-3)) / x1;
//     double f1bbar = (ReactionTheory::xfx(x1,Q,-5)) / x1;
//     double f2d = (ReactionTheory::xfx(x2,Q,1)) / x2;
//     double f2s = (ReactionTheory::xfx(x2,Q,3)) / x2;
//     double f2b = (ReactionTheory::xfx(x2,Q,5)) / x2;
// 
//     // PDF combinations    
//     double dbard_PDF = f1dbar*f2d + f1sbar*f2s + f1bbar*f2b;
//     
//     // Angular integration limits
//     double qqbar_cos_theta_max = min(cos(2*atan(exp(-eta_cut-y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
//     double qqbar_cos_theta_min = 0;  
//  
//     double angular_integration_EF = (qqbar_cos_theta_max-qqbar_cos_theta_min)+(1.0/3.0)*(pow(qqbar_cos_theta_max,3)-pow(qqbar_cos_theta_min,3));
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_EF = dsigma*angular_integration_EF;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // DBAR-D
//     double dbardEB = dbard_PDF*dsigma_EF;
//   
//     double *propagator = propagators (Minv);
//     
//     return dbardEB * propagator[2]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////DBARD EVEN BACKWARD Integration
// double ReactionAFB::integration_dbardEB (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//     
//     double dbardEB, error_dbardEB;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_dbardEB = { &(ReactionAFB::dbardEB_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_dbardEB, integration_inf, integration_sup, dim_integration, calls, r, s, &dbardEB, &error_dbardEB);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_dbardEB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &dbardEB, &error_dbardEB);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*dbardEB; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////DBARD ODD FORWARD Matrix element
// double ReactionAFB::dbardOF_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1dbar = (ReactionTheory::xfx(x1,Q,-1)) / x1;
//     double f1sbar = (ReactionTheory::xfx(x1,Q,-3)) / x1;
//     double f1bbar = (ReactionTheory::xfx(x1,Q,-5)) / x1;
//     double f2d = (ReactionTheory::xfx(x2,Q,1)) / x2;
//     double f2s = (ReactionTheory::xfx(x2,Q,3)) / x2;
//     double f2b = (ReactionTheory::xfx(x2,Q,5)) / x2;
// 
//     // PDF combinations    
//     double dbard_PDF = f1dbar*f2d + f1sbar*f2s + f1bbar*f2b;
//     
//     // Angular integration limits
//     double qbarq_cos_theta_max = 0;
//     double qbarq_cos_theta_min = max(cos(PI - 2*atan(exp(-eta_cut-y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2)))); 
//  
//     double angular_integration_OB = pow(qbarq_cos_theta_max,2) - pow(qbarq_cos_theta_min,2);
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_OB = dsigma*angular_integration_OB;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // DBAR-D
//     double dbardOF = dbard_PDF*dsigma_OB;
//   
//     double *propagator = propagators (Minv);
//     
//     return dbardOF * propagator[3]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////DBARD ODD FORWARD Integration
// double ReactionAFB::integration_dbardOF (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//        
//     double dbardOF, error_dbardOF;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_dbardOF = { &(ReactionAFB::dbardOF_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_dbardOF, integration_inf, integration_sup, dim_integration, calls, r, s, &dbardOF, &error_dbardOF);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_dbardOF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &dbardOF, &error_dbardOF);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*dbardOF; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }
// 
// ////DBARD ODD BACKWARD Matrix element
// double ReactionAFB::dbardOB_funct (double *entries, size_t dim, void *params)
// {
//     (void)(dim); /* avoid unused parameter warnings */
//     double yreduced = entries[0];
//     double Minv = entries[1];
//         
//     // Partonic cross section parameters
//     double Q = Minv;
//     double z = pow(Minv,2)/pow(energy,2);
//     double y = -(1.0/2.0)*log(z)*(yreduced);
//     double x1 = sqrt(z)*exp(y);
//     double x2 = sqrt(z)*exp(-y);
//     double dsigma_temp = pow(Minv,2)/(96*PI);
//     double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
//     
//     // Call the PDFs
//     double f1dbar = (ReactionTheory::xfx(x1,Q,-1)) / x1;
//     double f1sbar = (ReactionTheory::xfx(x1,Q,-3)) / x1;
//     double f1bbar = (ReactionTheory::xfx(x1,Q,-5)) / x1;
//     double f2d = (ReactionTheory::xfx(x2,Q,1)) / x2;
//     double f2s = (ReactionTheory::xfx(x2,Q,3)) / x2;
//     double f2b = (ReactionTheory::xfx(x2,Q,5)) / x2;
// 
//     // PDF combinations    
//     double dbard_PDF = f1dbar*f2d + f1sbar*f2s + f1bbar*f2b;
//     
//     // Angular integration limits
//     double qqbar_cos_theta_max = min(cos(2*atan(exp(-eta_cut-y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
//     double qqbar_cos_theta_min = 0; 
//  
//     double angular_integration_OF = pow(qqbar_cos_theta_max,2) - pow(qqbar_cos_theta_min,2);
//     
//     // Combination with angular integration (Forward - Backward for q-qbar)
//     double dsigma_OF = dsigma*angular_integration_OF;
// 
//     // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
//     // DBAR-D
//     double dbardOB = dbard_PDF*dsigma_OF;
//   
//     double *propagator = propagators (Minv);
//     
//     return dbardOB * propagator[3]; // Multiply the PDFs combination with the correct propagator.
// }
// 
// ////DBARD ODD BACKWARD Integration
// double ReactionAFB::integration_dbardOB (double Minv_inf, double Minv_sup)
// {
//     double integration_inf[2] = {yreducedmin, Minv_inf};
//     double integration_sup[2] = {yreducedmax, Minv_sup};
//     
//     double dbardOB, error_dbardOB;
//     // Initialization of the integration (quite a black box)
//     gsl_monte_function Integrate_dbardOB = { &(ReactionAFB::dbardOB_funct), dim_integration, 0 };
//     const gsl_rng_type *T;
//     gsl_rng *r;
//     gsl_rng_env_setup ();
//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);   
//     // Integration
//     {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
//     gsl_monte_vegas_integrate (&Integrate_dbardOB, integration_inf, integration_sup, dim_integration, calls, r, s, &dbardOB, &error_dbardOB);
//     do
//     {
//         gsl_monte_vegas_integrate (&Integrate_dbardOB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &dbardOB, &error_dbardOB);
//     }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//     gsl_monte_vegas_free (s);
//     }
// 
//     gsl_rng_free (r);
//     
//     return 2*dbardOB; // Factor 2 for the symmetric integration over yreduced: int_-1^1 dyreduced = 2 int_0^1 dyreduced.
// }


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
//     // Reconstructed AFB and statistical absolute uncertainty
//     double AFB = (Forward - Backward) / (Forward + Backward);
//     double Delta_AFB = sqrt((1.0-pow(AFB,2)) / (XS*lum));
    
//     double *results = new double[4];
//     
//     results[0] = XS;
//     results[1] = Epsilon_XS;
//     results[2] = AFB;
//     results[3] = Delta_AFB;
    
    
    
//     // Call the PDFs
//     double f1u = (ReactionTheory::xfx(0.01,100,2)) / 0.01;
//     printf("f1u(x1 = 0.01, 100, 2)/x1 = %1f \n", f1u);
//     
//     printf("uubarEF = %1f \n", uubarEF);
//     printf("uubarEB = %1f \n", uubarEB);
//     printf("uubarOF = %1f \n", uubarOF);
//     printf("uubarOB = %1f \n", uubarOB);
// 
//     printf("ubaruEF = %1f \n", ubaruEF);
//     printf("ubaruEB = %1f \n", ubaruEB);
//     printf("ubaruOF = %1f \n", ubaruOF);
//     printf("ubaruOB = %1f \n", ubaruOB);    
//     
//     printf("ddbarEF = %1f \n", ddbarEF);
//     printf("ddbarEB = %1f \n", ddbarEB);
//     printf("ddbarOF = %1f \n", ddbarOF);
//     printf("ddbarOB = %1f \n", ddbarOB);
// 
//     printf("dbardEF = %1f \n", dbardEF);
//     printf("dbardEB = %1f \n", dbardEB);
//     printf("dbardOF = %1f \n", dbardOF);
//     printf("dbardOB = %1f \n", dbardOB);
//     
//     printf("AFB = %1f \n", AFB);
    
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
    
    double *results = observables (min[0], max[0]);
    
    val[0] = results[2];

    
  return 0;
}
