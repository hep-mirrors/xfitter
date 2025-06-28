/*
   @file Fantomas_PdfParam.cc
   @date 2022-04-07
   @author  Lucas Kotz
   Created by  Lucas Kotz on 2022-04-07
*/

#include "Fantomas_PdfParam.h"
#include "xfitter_cpp_base.h"
#include <cmath>
#include <iostream>
//lk25 changed from fantomas.cc to metamorphCollection.h after pavel's revisions
//#include "metamorphCollection.h"
#include "CToWrapper.h"
bool xFitterCollectionSet = false;
bool xFitterModulatorSet = false;

//lk25 removed fantomas.cc. Moved metacol object to Fantomas_PdfParam.cc
//metamorphCollection metacol = metamorphCollection();

namespace xfitter{
//for dynamic loading
extern"C" Fantomas_PdfParam*create(const char*name){
  return new Fantomas_PdfParam(name);
}
void Fantomas_PdfParam::atStart(){
  using namespace std;
  BasePdfParam::atStart();
  const size_t n=getNPar();
  int ifl = int (*pars[n-1]);
  if(n<3){
    cerr<<"[ERROR] Too few parameters given to parameterisation \""<<_name<<"\", expected at least 3, got "<<n<<endl;
    hf_errlog(18120700,"F: Wrong number of parameters for a parameterisation, see stderr");
  }
  if (xFitterCollectionSet == false)
  {  
    //metacol.ReadCard();
    readfantosteer_();
    xFitterCollectionSet = true;
  }
  if (xFitterCollectionSet == true)
  {
  }
  atIteration();
  xFitterModulatorSet = true;
}
  
// Update Fantomas parameters each time minuit varies them
void Fantomas_PdfParam::atIteration(){
  const unsigned int npar=getNPar();
  updateParameters();
}

void Fantomas_PdfParam::updateParameters(){
  if (xFitterCollectionSet == false)
  {
    std::cout << "Metamorph Collection not set. Call readfantosteer() before updating modulators." << std::endl;
  }
  const unsigned int n=getNPar();
  int ifl = int (*pars[n-1]);
  double parstmp[n-1]={0};
  for (int i = 0; i < n-1; i++)
    parstmp[i] = *pars[i];
  //std::cout << "[DEBUG] ifl = " << ifl << ", expected deltas size = " << (n-1) << std::endl;
  //metacol.UpdateParameters(ifl,parstmp);
  updatefantopars_(ifl,parstmp);
}

// Main function to compute PDF
double Fantomas_PdfParam::operator()(double x)const{
  if (xFitterModulatorSet == false)
  {
    std::cout << "Modulator functions have not been set. Make sure Fantomas_PdfParam::atStart() is being called." << std::endl;
  }
  const unsigned int npar = getNPar();

if (!pars || pars[npar - 1] == nullptr) {
    std::cerr << "[ERROR] Null pointer encountered in pars at index " << npar - 1 << std::endl;
    std::terminate();
}

  //std::cout << "npar: " << npar <<std::endl;
  int ifl = *pars[npar-1];
  //std::cout << "ifl: " << ifl << std::endl;
  //std::cout << "x: " << x << std::endl;
  // lk22 removed pars[0] from f since normalization is now added to metamorph function.
  //double f = metacol.f(ifl,x);
  double f = fantopara_(ifl,x);
  return f;
}
double Fantomas_PdfParam::moment(int n)const{
  //Integral of HERAPDF-style function is expressed in terms of Euler beta function:
  //beta(x,y)=int_0^1 t^(x-1)*(1-t)^(y-1) dx=gamma(x)*gamma(y)/gamma(x+y)
  //moment(n)=int_0^1 P[0]*x^(P[1]+n)*(1-x)^P[2]*(1+sum_{i=3}^N{P[i]*x^(i-2)})
  //Let A:=P[0], B:=P[1]+n+1, C:=P[2]+1, then
  //moment=int_0^1 A*x^(B-1)*(1-x)^(C-1)*(1+sum_{i=3}^N{P[i]*x^(i-2)})
  //=A*(beta(B,C)+sum_{i=3}^N{P[i]*beta(B+i-2,C)})
  //beta(B+1,C)=B/(B+C)*beta(B,C)
  //=> beta(B+n,C)=beta(B,C)*product_{k=0}^{n-1}{(B+k)/(B+C+k)}
  //=> beta(B+i-2,C)=beta(B,C)*product_{k=0}^{i-3}{(B+k)/(B+C+k)}
  //=> moment=A*beta(B,C)*(1+sum_{i=3}^N{P[i]*product_{k=0}^{i-3}{(B+k)/(B+C+k)}})=
  //=> moment=A*beta(B,C)*(1+P[3]*B/(B+C)+P[4]*B/(B+C)*(B+1)/(B+C+1)+...)
  //beta(B,C)=exp(lgamma(B)+lgamma(C)-lgamma(B+C))
  /*const double B=(*pars[1])+(n+1),C=(*pars[2])+1;
  if(B<=0.||C<=0.)return NAN;// integral does not converge
  const size_t N=getNPar();
  double sum=1;
  double prod=1;
  double a=B;
  double b=B+C;
  for(uint i=3;i<N;++i){
    prod=prod*a/b;
    sum+=(*pars[i])*prod;
    a++;
    b++;
  }
  return (*pars[0])*exp(lgamma(B)+lgamma(C)-lgamma(B+C))*sum;
  */
  const unsigned int npar=getNPar();
  int ifl = *pars[npar-1];
  // int npts = 5000; // Uncomment line to change number of integration points in adxmoment integration. Default is 10000.
  // lk22 removed pars[0] from moment since normalization is now added to metamorph function.
  //double moment = metacol.MellinMoment(ifl,n/*,npts*/);
  double moment = fantomellinmoment_(ifl,n);
  return moment;
}

void Fantomas_PdfParam::setMoment(int n,double val){
  BasePdfParam::setMoment(n,val);
  updateParameters();
}

}
