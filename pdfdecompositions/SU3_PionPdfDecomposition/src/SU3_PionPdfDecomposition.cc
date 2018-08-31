/*
   @file SU3_PionPdfDecomposition.cc
   @date 2018-08-07
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on SU3_Pion
*/
#include"SU3_PionPdfDecomposition.h"
#include"PolySqrtPdfParam.h" //This is for hacks, remove later
#include"xfitter_pars.h"
using uint=unsigned int;
using namespace std;
namespace xfitter{
//For dynamic loading:
extern "C" SU3_PionPdfDecomposition*create(){
	return new SU3_PionPdfDecomposition();
}
SU3_PionPdfDecomposition::SU3_PionPdfDecomposition():BasePdfDecomposition("SU3_Pion"){
	par_v=nullptr;
	par_S=nullptr;
	par_g=nullptr;
}
SU3_PionPdfDecomposition::SU3_PionPdfDecomposition(const std::string& inName):BasePdfDecomposition(inName){}
SU3_PionPdfDecomposition::~SU3_PionPdfDecomposition(){
	if(par_v)delete par_v;
	if(par_S)delete par_S;
	if(par_g)delete par_g;
}
// Init at start:
void SU3_PionPdfDecomposition::initAtStart(const std::string & pars){
	//HARDCODE copied from UvDvUbarDbarS, then modified
	//The following is not very nice: Decomposition should not create or initialize parameterisations
	//Create and initialize paramterisations
	for(const auto&node:XFITTER_PARS::gParametersY.at("SU3_PionPdfDecomposition")["HARDWIRED_PolySqrt"]){
		const string prmzName=node.first.as<string>();//Name of parameterisation
		BasePdfParam*pParam=new PolySqrtPdfParam(prmzName);
		pParam->initFromYaml(node.second);
		addParameterisation(prmzName,pParam);
	}
	par_v=getPdfParam("v");
	par_S=getPdfParam("S");
	par_g=getPdfParam("g");
}
void SU3_PionPdfDecomposition::initAtIteration() {
	//Enforce sum rules
	//Valence sum
	par_v->setMoment(-1,1);
	//Momentum sum
	par_g->setMoment(0,1-6*par_S->moment(0)-2*par_v->moment(0));
}
// Returns a LHAPDF-style function, that returns PDFs in a physical basis for given x
std::function<std::map<int,double>(const double& x)>SU3_PionPdfDecomposition::f0()const{
	return [=](double const& x)->std::map<int, double>{
		double v=(*par_v)(x);
		double S=(*par_S)(x);
		double g=(*par_g)(x);
		double d=S+v;
		std::map<int,double>res_={
			{-6,0},
			{-5,0},
			{-4,0},
			{-3,S},//sbar
			{-2,d},//ubar
			{-1,S},//dbar
			{ 1,d},//d
			{ 2,S},//u
			{ 3,S},//s
			{ 4,0},
			{ 5,0},
			{ 6,0},
			{21,g}
		};
		return res_;
	};
}
}
