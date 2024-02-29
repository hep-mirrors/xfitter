#include"BasePdfParam.h"
#include"BaseMinimizer.h"
#include"xfitter_pars.h"
#include"xfitter_steer.h"
#include"xfitter_cpp_base.h"
#include<cmath>
#include<memory>
#include<iostream>

namespace xfitter{
BasePdfParam::~BasePdfParam(){if(pars)delete[]pars;}
double BasePdfParam::moment(int iMoment)const{
	/// Numeric integration
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
void BasePdfParam::atStart(){
  using namespace std;
  YAML::Node node=XFITTER_PARS::getParameterisationNode(_name);
  YAML::Node parsNode=node["parameters"];
  if(!parsNode.IsSequence()){
    cerr<<"[ERROR] Bad \"parameters\" for parameterisation \""<<_name<<"\", expected sequence"<<endl;
    hf_errlog(18092420,"F: Bad parameters in parameterisation, see stderr");
  }
  Npars=parsNode.size();
  if (pars == nullptr)
    pars=new double*[Npars];
  for(unsigned int i=0;i<Npars;++i){
    try{
      pars[i]=XFITTER_PARS::gParameters.at(parsNode[i].as<string>());
    }catch(const YAML::BadConversion&ex){
      cerr<<"[ERROR] Bad name of parameter "<<i<<" for parameterisation \""<<_name<<"\": "<<ex.what()<<endl;
      hf_errlog(18092421,"F: Bad parameter name in parameterisation definition, see stderr");
    }catch(const out_of_range&ex){
      string parname=parsNode[i].as<string>();
      if(XFITTER_PARS::gParameters.count(parname)==0){
        cerr<<"[ERROR] Unknown parameter \""<<parname<<"\" in parameterisation \""<<_name<<"\""<<endl;
        hf_errlog(18092422,"F: Unknown parameter in parameterisation definition, see stderr");
      }
    }
  }
}
void BasePdfParam::atIteration(){}
}
