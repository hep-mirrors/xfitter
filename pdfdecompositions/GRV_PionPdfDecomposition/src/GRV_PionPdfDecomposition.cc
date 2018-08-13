/*
   @file GRV_PionPdfDecomposition.cc
   @date 2018-08-07
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on GRV_Pion
*/
#include"GRV_PionPdfDecomposition.h"
using uint=unsigned int;
using namespace std;
namespace xfitter{
class ParameterisationWrapper{
public:
	ParameterisationWrapper(BasePdfParam*base,const string decompositionName):base{base}{
		const uint N=base->getNPar();
		pars=new double[N];
		const string prefix=decompositionName+"_"+base->getName()+"_p";
		for(uint i=0;i<N;++i){
			pars[i]=XFITTER_PARS::gParameters[prefix+to_string(i)];
		}
	}
	~ParameterisationWrapper(){delete[]pars;}
	double operator(double x){
		std::unique_ptr<double[]>p(new double[N]);
		for(uint i=0;i<N;++i)p[i]=*pars[i];
		return base->compute(x,p);
	}
	double moment(int i){//return \int_0^1 x^i*f(x) dx
		const uint N=base->getNPar();
		std::unique_ptr<double[]>p(new double[N]);
		for(uint i=0;i<N;++i)p[i]=*pars[i];
		return base->moment(p,i);
	}
	void setMoment(int i,double S){
		*pars[0]=1;
		*pars[0]=S/moment(i);
	}
private:
	BasePdfParam*base;
	double**pars;
};
//For dynamic loading:
extern "C" GRV_PionPdfDecomposition*create(){
	return new GRV_PionPdfDecomposition();
}
GRV_PionPdfDecomposition::GRV_PionPdfDecomposition():BasePdfDecomposition("GRV_Pion"){
	par_v   =nullptr;
	par_qbar=nullptr;
	par_g   =nullptr;
}
GRV_PionPdfDecomposition::GRV_PionPdfDecomposition(const std::string& inName):BasePdfDecomposition(inName){}
GRV_PionPdfDecomposition::~GRV_PionPdfDecomposition(){
	if(par_v   )delete par_v;
	if(par_qbar)delete par_qbar;
	if(par_g   )delete par_g;
}
// Init at start:
void GRV_PionPdfDecomposition::initAtStart(const std::string & pars){
	//HARDCODE copied from UvDvUbarDbarS, then modified
	//TO BE REPLACED LATER
	for(const auto&node:XFITTER_PARS::gParametersY.at("GRV_PionPdfDecomposition")["HERAPDF_pdfparam"]){
		const string prmzName=node.first.as<string>();//Name of parameterisation
		BasePdfParam*pParam=new HERAPDF_PdfParam(prmzName);
		double*parValues=pParam->initFromYaml(node.second);//Parameters of the parameterisation
		//Register parameters of this parameterisation in the global parameter map
		for(int i=0;i<pParam->getNPar();i++){
			double val     =parValues[i];
			double step    =std::fabs(val)/100.;       /// if 0, parameter is fixed !!! 
			double minv    =0;//0 means no bounds for minuit
			double maxv    =0;
			double priorVal=0;
			double priorUnc=0;
			int    add     =true;
			std::ostringstream ss;
			ss<<getName()<<'_'<<prmzName<<"_p"<<i;
			const sring pnam=ss.str();
			std::cout<<"INFO[GRV_PionPdfDecomposition]: Registering parameter "<<ss<<"="<<val<<std::endl;
			/// Here it goes to minuit:
			addexternalparam_(pnam.c_str(),val,step,minv,maxv,priorVal,priorUnc,add,&XFITTER_PARS::gParameters,pnam.size() );
		}
    addParameterisation(pdfName,pParam);
  }
	par_v   =new ParameterisationWrapper(getPdfParam("v   "),getName());
	par_qbar=new ParameterisationWrapper(getPdfParam("qbar"),getName());
	par_g   =new ParameterisationWrapper(getPdfParam("g   "),getName());
}
void GRV_PionPdfDecomposition:initAtIteration() {
	//Enforce sum rules
	//Valence sum
	par_v->setMoment(0,2);
	//Momentum sum
	par_g->setMoment(1,1-4*par_qbar->moment(1));
}
// Returns a LHAPDF-style function, that returns PDFs in a physical basis for given x
std::function<std::map<int,double>(const double& x)>GRV_PionPdfDecomposition::f0()const{
	return [=](double const& x)->std::map<int, double>{
		double v   =(*par_v)(x);
		double qbar=(*par_qbar)(x);
		double g   =(*par_g)(x);
		double u=qbar-v/4;
		double d=qbar+v/4;
		std::map<int,double>res_={
			{-6,0},
			{-5,0},
			{-4,0},
			{-3,0},
			{-2,d},//ubar
			{-1,u},//dbar
			{ 1,d},
			{ 2,u},
			{ 3,0},
			{ 4,0},
			{ 5,0},
			{ 6,0},
			{21,g}
		};
		return res_;
	};
}
}
