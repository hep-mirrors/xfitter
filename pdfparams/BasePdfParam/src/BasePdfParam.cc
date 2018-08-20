#include"BasePdfParam.h"
#include"BaseMinimizer.h"
#include"xfitter_pars.h"
#include"xfitter_steer.h"
#include<cmath>
#include<memory>
#include<iostream>

/// TEMPORARY XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
extern "C" {
  /// Interface to minuit parameters
  void addexternalparam_(const char name[],  const double &val, 
                         const double  &step,
                         const double &min, const double &max, 
                         const double &prior, const double &priorUnc,
                         const int &add, 
                         map<std::string,double*> *map,
                         int len);
}
/// Implement numeric integration
double BasePdfParam::moment(int iMoment)const{
  /// Simple rule, split log/lin spacing at xsplit=0.1

  const double xsplit = 0.1;
  const double xmin = 1e-6;
  const double xminlog = log10(xmin);
  const double xmaxlog = log10(xsplit);
  const int nlog = 100;
  const int nlin = 100;

  const double xsteplog = (xmaxlog-xminlog)/nlog;
  const double xsteplin = (1.0 - xsplit)/nlin;
  
  double sum = 0.;
  double x = xmin;

  // log x part:
  for (int i=0; i<nlog; i++) {
    double dx = pow(10,xminlog+(i+1)*xsteplog) - pow(10,xminlog+i*xsteplog);
    double val = (*this)(x+dx/2.)*pow(x+dx/2.,iMoment);
    x += dx;
    sum += dx*val;
  }
  // lin x part:
  for (int i=0; i<nlin; i++) {
    double dx = xsteplin;
    double val = (*this)(x+dx/2.)*pow(x+dx/2.,iMoment);
    x += dx;
    sum += dx*val;
  }
  
  return sum;
}
void BasePdfParam::setMoment(int nMoment,double value){
  *pars[0]=1;
  *pars[0]=value/moment(nMoment);
}
void BasePdfParam::initFromYaml(YAML::Node value) {
  // HARDWIRE A BIT FOR NOW XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  //TODO rewrite this for a different, new YAML format
  using namespace std;
  using uint=unsigned int;
  cout<<"DEBUG["<<_name<<"]: initFromYaml: value="<<value<<endl;
  if(value.IsSequence()){
    Npars=value.size();
    cout<<Npars<<endl;
    pars=new double*[Npars];
    // HARDWIRE old-way for now:  XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    for(uint i=0;i<Npars;++i){
      // get a good name
      const std::string pnam=_name+"_p"+std::to_string(i);
      double val     =value[i].as<double>();
      double step    =fabs(val)/100.;       /// if 0, parameter is fixed !!! 
      //double minv    =0;
      //double maxv    =0;
      //double priorVal=0;
      //double priorUnc=0;
      //int add = true;
      //      addexternalparam_(pnam.c_str(),val,step,minv,maxv,priorVal,priorUnc,add,&XFITTER_PARS::gParameters,pnam.size());

      xfitter::BaseMinimizer* minimizer = xfitter::get_minimizer();
      minimizer->addParameter(val,pnam,step,nullptr,nullptr);
      
      pars[i]=XFITTER_PARS::gParameters.at(pnam);
      cout<<pnam<<"="<<(*pars[i])<<endl;
    }
  }else{
    cout<<"ERROR["<<_name<<"]: initFromYaml: parameter is not a sequence!"<<endl;
  }
  //TODO: Handle possible errors
}

