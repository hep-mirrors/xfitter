#include"BasePdfParam.h"
#include"BaseMinimizer.h"
#include"xfitter_pars.h"
#include"xfitter_steer.h"
#include"xfitter_cpp_base.h"
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
namespace xfitter{
/// Implement numeric integration
BasePdfParam::~BasePdfParam(){if(pars)delete[]pars;}
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
void BasePdfParam::atStart(){
  using namespace std;
  YAML::Node node=XFITTER_PARS::getParameterisationNode(_name);
  YAML::Node parsNode=node["parameters"];
  if(!parsNode.IsSequence()){
    cerr<<"[ERROR] Bad \"parameters\" for parameterisation \""<<_name<<"\", expected sequence"<<endl;
    hf_errlog(18092420,"F: Bad parameters in parameterisation, see stderr");
  }
  Npars=parsNode.size();
  pars=new double*[Npars];
  //TODO: destructor
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
  /*
  using uint=unsigned int;
  //cout<<"DEBUG["<<_name<<"]: initFromYaml: value="<<value<<endl;
  if(value.IsSequence()){
    Npars=value.size();
    //cout<<Npars<<endl;
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

    }
  }else{
    cout<<"ERROR["<<_name<<"]: initFromYaml: parameter is not a sequence!"<<endl;
  }
  */
}
}
