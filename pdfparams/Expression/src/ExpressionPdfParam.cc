#include"ExpressionPdfParam.h"
#include"xfitter_pars.h"
#include"xfitter_cpp_base.h"
#include"expression_utils.h"
#include<iostream>
using namespace std;
namespace xfitter{
//returns true if expression is not a sum
bool isProduct(const string&s){
  size_t d=0;//current bracket depth
  for(const char*p=s.c_str();*p;++p){
    if(*p=='(')d++;
    else if(*p==')')d--;
    else if(d==0&&(*p=='+'||*p=='-'))return false;
  }
  return true;
}
//if expression has form parname*other_stuff, return parname
//otherwise return ""
string getNormalizationParameter(const string&s){
  const char*p=s.c_str();
  while(*p==' '){//skip leading whitespace
    if(!*p)return"";
    ++p;
  }
  const char*b=p;//start of found substring
  while(isalnum(*p)||*p=='_')++p;//find end of parameter name
  size_t n=p-b;//substring size
  if(n==1&&*b=='x')return "";
  //make sure normalization factor is followed by '*' or '/'
  while(*p==' ')++p;//skip whitespace
  if(!(*p=='*'||*p=='/'))return "";
  return string(b,n);
}
//for dynamic loading
extern"C" ExpressionPdfParam*create(const char*s){return new ExpressionPdfParam(s);}
double ExpressionPdfParam::operator()(double _x)const{
  x=_x;
  return te_eval(expr);
}
void ExpressionPdfParam::atStart(){
  YAML::Node node=XFITTER_PARS::getParameterisationNode(_name)["expression"];
  string expression;
  try{
    expression=node.as<string>();
  //}catch(const YAML::InvalidNode&ex){
  }catch(const YAML::BadConversion&ex){
    if(!node){
      cerr<<"[ERROR] No expression given for parameterisation \""<<_name<<"\""<<endl;
      hf_errlog(19022600,"F: No expression given for parameterisation, see stderr");
    }else{
      cerr<<"[ERROR] Failed to convert expression given for decomposition \""<<_name<<"\" to string"<<endl;
      hf_errlog(19022601,"F: Invalid expression parameter for parameterisation, see stderr");
    }
  }
  vector<string>pars;
  extractParameterNames(expression,pars);
  size_t Npars=pars.size();
  te_variable*vars=new te_variable[Npars+1];//+1 because of x
  for(size_t i=0;i<Npars;++i){
    const string&parname=pars[i];
    vars[i].name=parname.c_str();
    try{
      vars[i].address=XFITTER_PARS::gParameters.at(parname);
    }catch(const std::out_of_range&e){
      cerr<<"[ERROR] Unknown parameter \""<<parname<<"\" in expression \""<<expression<<"\" for parameterisation \""<<_name<<"\""<<endl;
      hf_errlog(19022610,"F: Unknown parameter in parameterisation expression, see stderr");
    }
    vars[i].type=TE_VARIABLE;
    vars[i].context=nullptr;
  }
  vars[Npars].name="x";
  vars[Npars].address=&x;
  vars[Npars].type=TE_VARIABLE;
  vars[Npars].context=nullptr;
  int err;
  expr=te_compile(expression.c_str(),vars,Npars+1,&err);
  if(!expr){
    cerr<<"[ERROR] TinyExpr error while parsing expression \""<<expression<<"\" for parameterisation \""<<_name<<"\"; error at position "<<err<<endl;
    hf_errlog(19022611,"F: Failed to parse expression for parameterisation, see stderr");
  }
  delete[]vars;
  if(isProduct(expression)){
    //Try to guess which parameter should be used to set moment
    string normParName=getNormalizationParameter(expression);
    if(normParName=="")goto cant_moment;
    normPar=XFITTER_PARS::gParameters.at(normParName);
    return;
  }//else
  cant_moment:
  normPar=nullptr;
  hf_errlog(19022612,"W: Expression parameterization: cannot set moment");
}
void ExpressionPdfParam::setMoment(int n,double val){
  if(!normPar){
    cerr<<"[ERROR] Do not know which parameter to scale to set moment of parameterisation \""<<_name<<"\""<<endl;
    hf_errlog(19022613,"F: Expression parameterisation cannot set moment, see stderr");
  }
  *normPar=1;
  *normPar=val/moment(n);
}
}
