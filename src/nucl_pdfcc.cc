// 09-02-2017 created by Marina Walt, University of Tuebingen
// module to summarize all modification to xFitter in c++ meant to implement the treatment of nuclear PDFs

// in order to be able to call the new routines from the old subprograms of xFitter
// nucl_pdfcc.cc needs to be added to the Makefile in /src folder 
// (and commands make clean, make, install to be executed after the modification of the Makefile).

#include <map>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

using std::cout;
using std::endl;
using namespace std;

#include "../include/xfitter_cpp.h"
//#include "gsl/gsl_integration.h"
//#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

// Define struct (for fortran common blocks)
extern "C" {
  extern struct{
    double ctglue[27];
    double ctuval[27];
    double ctdval[27];
    double ctubar[27];
    double ctdbar[27];
    double ctother[27];
    double ctstr[27];
    double ctphoton[27];
    int idxADataSet;
    int ratiostep;
  } ctpar_;
 extern struct{
     bool DEBUG;
  } for_debug_;  
}

// Fortran interface
extern "C" {
  void integratecncteq_(int &n, double ctpara[27], double &a, double &b, double &result, double &abserr);
  void integratecntuju_(int &n, double ctpara[27], double &a, double &b, double &result, double &abserr);  
}


void integrateab(double (* f) (double x, void * p), void * data, double &a, double &b, double &result, double &abserr)
{
  int integ_space = 10000; 
  double rel_error = 1e-3;
  double abs_error = 0;
  int key = 5;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(integ_space);
  gsl_function F;
  F.function = f;
  F.params = data;
  
  
  
// 11.01.2018 Error reporting
  
// 22.01.2018: error handler deactivated in order to prevent system from abort.
  gsl_set_error_handler_off(); 
  
  int status = gsl_integration_qag(&F, a, b, abs_error, rel_error, integ_space, 5, w, &result, &abserr);
  if (status) { /* an error occurred */
  cout << "GSL integration error occured:" << endl;
  cout << gsl_strerror(status) << endl; 
// Error handling  
  if (status == GSL_ERANGE) {
//      gsl_set_error_handler_off(); 
      cout << "error: GSL_ERANGE: result not representable, overflow/underflow" << endl;
      result = 0;
      return; 
  }
  else if (status == GSL_ENOMEM) {
//      gsl_set_error_handler_off();
      gsl_integration_workspace_free (w);
      key = 1;
      integ_space = 1000;
      gsl_integration_workspace * w = gsl_integration_workspace_alloc(integ_space);
      return;
  }
  else {
      cout << "Unknown error status" << endl;
  };
  };
//
  
  gsl_integration_qag(&F, a, b, abs_error, rel_error, integ_space, key, w, &result, &abserr);
  gsl_integration_workspace_free (w);
  return;
}

///////////////////////////

struct nCTEQ_params {double nparam; double ctpara[27];};

double nCTEQfunction (double x, void * p) //not used in this code, for reference purpose only
{
 nCTEQ_params * params = (nCTEQ_params*)p;
  // parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
 double nC = params->nparam; 
 double ctparaC[27];
 int i = 0;
 do {
   ctparaC[i] = params->ctpara[i];
   i=i+1;
  }while(i<26);
  
 //debugging
bool debugC = for_debug_.DEBUG;
 if(debugC) {
   cout << "ctparaC Test in nCTEQfunction" << endl;
   cout << ctparaC[1]<< endl; 
 };
 
 // function (without normalization constant)
 // parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
 return pow(x,(ctparaC[1]+nC))*pow((1-x),ctparaC[2])*exp(ctparaC[3]*x)*pow((1+exp(ctparaC[4])*x),ctparaC[5]);
}

////////////////////////////////////////////////////////////////////////////////////////////////
double nCTEQfunction0 (double x, void * p)
{
 nCTEQ_params * params = (nCTEQ_params*)p;
  // parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
 int nC = 0; 
 double ctparaC[27];
 int i = 0;
 do {
   ctparaC[i] = params->ctpara[i];
   i=i+1;
  }while(i<26);
 
  cout << "ctparaC[0] Test in nCTEQfunction0" << endl;
  cout << ctparaC[0]<< endl; 
  
  cout << "ctparaC[1] Test in nCTEQfunction0" << endl;
  cout << ctparaC[1]<< endl; 
  
  cout << "ctparaC[1+9] Test in nCTEQfunction0" << endl;
  cout << ctparaC[10]<< endl; 
  
 //debugging
bool debugC = for_debug_.DEBUG;
 if(debugC) {
  cout << "ctparaC Test in nCTEQfunction0" << endl;
  cout << ctparaC[1]<< endl; 
 };
 
 // function (without normalization constant)
 // parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
 return pow(x,(ctparaC[1]+nC))*pow((1-x),ctparaC[2])*exp(ctparaC[3]*x)*pow((1+exp(ctparaC[4])*x),ctparaC[5]);
}

////////////////////////////////////////////////////////////////////////////////////////////////
double nCTEQfunction1 (double x, void * p)
{
 nCTEQ_params * params = (nCTEQ_params*)p;
  // parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
 int nC = -1; 
 double ctparaC[27];
 int i = 0;
 do {
   ctparaC[i] = params->ctpara[i];
   i=i+1;
  }while(i<26);
  
  cout << "ctparaC[0] Test in nCTEQfunction1" << endl;
  cout << ctparaC[0]<< endl; 
  
  cout << "ctparaC[1] Test in nCTEQfunction1" << endl;
  cout << ctparaC[1]<< endl; 
  
  cout << "ctparaC[1+9] Test in nCTEQfunction1" << endl;
  cout << ctparaC[10]<< endl; 

 //debugging
 bool debugC = for_debug_.DEBUG;
 if(debugC) {
   cout << "ctparaC Test in nCTEQfunction1" << endl;
   cout << ctparaC[1]<< endl; 
 };
 
 // function (without normalization constant)
 // parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
 return pow(x,(ctparaC[1]+nC))*pow((1-x),ctparaC[2])*exp(ctparaC[3]*x)*pow((1+exp(ctparaC[4])*x),ctparaC[5]);
}
////////////////////////////////////////////////////////////////////////////////////////////////

void integratecncteq_(int &n, double ctpara[27], double &a, double &b, double &result, double &abserr)
{ 
  bool debugC = for_debug_.DEBUG;  
  double nparam = n;
  
  double ctparaC[27];
  int i = 0;
  do {
     ctparaC[i] = ctpara[i];
     i=i+1;
     }while(i<26);
    
//debugging
  if(debugC) {    
    cout << "ctparaC Test before nCTEQfunction" << endl;
    cout << ctparaC[1]<< endl;
  };
  
  //numerical integration
  //integrateab(nCTEQfunction, (&nparam, &ctpara), a, b, result, abserr);
  
  if (nparam ==0)
  {integrateab(nCTEQfunction0, &ctpara, a, b, result, abserr);
  }
  else if (nparam ==-1)
  {integrateab(nCTEQfunction1, &ctpara, a, b, result, abserr);
  }
  else
  {cout << "n-value not defined" <<endl;
  };
  
  
  //debugging
  if(debugC==1) {
    cout << "result in C:" << endl;
   // cout << a << endl;
   // cout << b << endl;
    cout << result << endl;
   // cout << abserr << endl;
  };
  return;
}

//////////////////////////////////////////
/////////////////////////////////////////
////////////////////////////////////////
// 16/04/2018

////////////////////////////////////////////////////////////////////////////////////////////////
double nTUJUfunction0 (double x, void * p)
{
 nCTEQ_params * params = (nCTEQ_params*)p;
  // parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
 int nC = 0; 
 double ctparaC[27];
 int i = 0;
 do {
   ctparaC[i] = params->ctpara[i];
   i=i+1;
  }while(i<26);
 
  cout << "ctparaC[0] Test in nTUJUfunction0" << endl;
  cout << ctparaC[0]<< endl; 
  
  cout << "ctparaC[1] Test in nTUJUfunction0" << endl;
  cout << ctparaC[1]<< endl; 
  
  cout << "ctparaC[1+9] Test in nTUJUfunction0" << endl;
  cout << ctparaC[10]<< endl; 
  
 //debugging
bool debugC = for_debug_.DEBUG;
 if(debugC) {
  cout << "ctparaC Test in nTUJUfunction0" << endl;
  cout << ctparaC[1]<< endl; 
 };
 
 // function (without normalization constant)
 // parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
 return pow(x,(ctparaC[1]+nC))*pow((1-x),ctparaC[2])*(1+ctparaC[3]*x+ctparaC[4]*pow(x,2));
}

////////////////////////////////////////////////////////////////////////////////////////////////
double nTUJUfunction1 (double x, void * p)
{
 nCTEQ_params * params = (nCTEQ_params*)p;
  // parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
 int nC = -1; 
 double ctparaC[27];
 int i = 0;
 do {
   ctparaC[i] = params->ctpara[i];
   i=i+1;
  }while(i<26);
  
  cout << "ctparaC[0] Test in nTUJUfunction1" << endl;
  cout << ctparaC[0]<< endl; 
  
  cout << "ctparaC[1] Test in nTUJUfunction1" << endl;
  cout << ctparaC[1]<< endl; 
  
  cout << "ctparaC[1+9] Test in nTUJUfunction1" << endl;
  cout << ctparaC[10]<< endl; 

 //debugging
 bool debugC = for_debug_.DEBUG;
 if(debugC) {
   cout << "ctparaC Test in nTUJUfunction1" << endl;
   cout << ctparaC[1]<< endl; 
 };
 
 // function (without normalization constant)
 // parameterisation where n = 0 (momentum sum-rule) or -1 ( counting sum rule)
 return pow(x,(ctparaC[1]+nC))*pow((1-x),ctparaC[2])*(1+ctparaC[3]*x+ctparaC[4]*pow(x,2));
}
////////////////////////////////////////////////////////////////////////////////////////////////





void integratecntuju_(int &n, double ctpara[27], double &a, double &b, double &result, double &abserr)
{ 
  bool debugC = for_debug_.DEBUG;  
  double nparam = n;
  
  double ctparaC[27];
  int i = 0;
  do {
     ctparaC[i] = ctpara[i];
     i=i+1;
     }while(i<26);
    
//debugging
  if(debugC) {    
    cout << "ctparaC Test before nTUJUfunction" << endl;
    cout << ctparaC[1]<< endl;
  };
  
  //numerical integration
  //integrateab(nCTEQfunction, (&nparam, &ctpara), a, b, result, abserr);
  
  if (nparam ==0)
  {integrateab(nTUJUfunction0, &ctpara, a, b, result, abserr);
  }
  else if (nparam ==-1)
  {integrateab(nTUJUfunction1, &ctpara, a, b, result, abserr);
  }
  else
  {cout << "n-value not defined" <<endl;
  };
  
  
  //debugging
  if(debugC==1) {
    cout << "result in C:" << endl;
   // cout << a << endl;
   // cout << b << endl;
    cout << result << endl;
   // cout << abserr << endl;
  };
  return;
}

