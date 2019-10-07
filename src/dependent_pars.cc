#include"dependent_pars.h"
#include"tinyexpr.h"
#include"xfitter_pars.h"
#include"xfitter_cpp_base.h"
#include"expression_utils.h"
#include<iostream>
#include<map>
using namespace std;
namespace xfitter{
//The Constraint class is used for dependent parameters
class Constraint{
public:
  Constraint(const string&parameterName,const string&expression,const vector<string>&dependencies);
  ~Constraint();
  void atIteration();
  double*parameter;//Constraint writes to this pointer at each iteration
private:
  te_expr*expr;//expression object
};
Constraint::Constraint(const string&parameterName,const string&expression,const vector<string>&deps){
  parameter=XFITTER_PARS::gParameters.at(parameterName);
  size_t Ndeps=deps.size();
  te_variable*vars=new te_variable[Ndeps];
  for(size_t i=0;i<Ndeps;++i){
    const string&parname=deps[i];
    vars[i].name=parname.c_str();
    try{
      vars[i].address=XFITTER_PARS::gParameters.at(parname);
    }catch(const std::out_of_range&e){
      cerr<<"[ERROR] Unknown parameter \""<<parname<<"\" in expression \""<<expression<<"\" for dependent parameter \""<<parameterName<<"\""<<endl;
      hf_errlog(18112000,"F: Unknown parameter in expression for dependent parameter, see stderr");
    }
    vars[i].type=TE_VARIABLE;
    vars[i].context=nullptr;
  }
  int err;
  expr=te_compile(expression.c_str(),vars,Ndeps,&err);
  if(!expr){
    cerr<<"[ERROR] TinyExpr error while parsing expression \""<<expression<<"\" for dependent parameter \""<<parameterName<<"\"; error at character "<<err<<endl;
    hf_errlog(18112001,"F: Failed to parse expression for dependent parameter, see stderr");
  }
  delete[]vars;
}
Constraint::~Constraint(){
  te_free(expr);
}
void Constraint::atIteration(){
  *parameter=te_eval(expr);
}
//array of all constraints, sorted in correct order
vector<Constraint*>constraints;
/*
  Dependent parameters may depend on other dependent parameters, which means that they need to be updated in a correct order.
  In order to determine that correct order, we build a directed acyclic graph (DAG) (dependency graph)
  Correct order of evaluation is then found by performing a topological sort of the graph.
*/
struct DAGnode{
  string name;//name of parameter
  string expression;
  vector<string>deps;//names of parameters this one depends on
  vector<DAGnode*>parents,children;//parents depend on this node, this node depends on children
};
//Remove leading whitespace and '='
string trimString(const string&s){
  const char*p=s.c_str();
  while(*p==' '||*p=='=')++p;
  return s.substr(size_t(p-s.c_str()));
}
void registerDependentParameters(const vector<DependentParameter>&dependentParameters){
  map<string,DAGnode*>nodemap;
  for(const DependentParameter&depPar:dependentParameters){
    DAGnode*node=new DAGnode();
    node->name=depPar.name;
    node->expression=trimString(depPar.expression);
    extractParameterNames(node->expression,node->deps);
    if(nodemap.count(node->name)>0){
      cerr<<"[ERROR] Redefinition of dependent parameter \""<<node->name<<'\"'<<endl;
      hf_errlog(18112002,"F: Redefinition of dependent parameter, see stderr");
    }
    nodemap[node->name]=node;
  }
  //connect nodes
  for(map<string,DAGnode*>::const_iterator it=nodemap.begin();it!=nodemap.end();++it){
    DAGnode*node=it->second;
    for(const string&deppar:node->deps){
      map<string,DAGnode*>::const_iterator it2=nodemap.find(deppar);
      if(it2!=nodemap.end()){
        DAGnode*child=it2->second;
        node->children.push_back(child);
        child->parents.push_back(node);
      }
    }
  }
  //topologically sort DAG, creating Constraint-s
  vector<DAGnode*>childless;
  for(map<string,DAGnode*>::const_iterator it=nodemap.begin();it!=nodemap.end();++it){
    DAGnode*node=it->second;
    if(node->children.empty())childless.push_back(node);
  }
  while(!childless.empty()){
    //remove this node
    DAGnode*node=childless.back();
    childless.pop_back();
    for(DAGnode*parent:node->parents){
      //remove pointer to this node from parent
      for(auto it=parent->children.begin();;++it)if(*it==node){
        *it=*(parent->children.end()-1);
        parent->children.pop_back();
        break;
      }
      //as parent now has one child less, it could become childless
      if(parent->children.empty())childless.push_back(parent);
    }
    //create constraint, delete this node
    constraints.push_back(new Constraint(node->name,node->expression,node->deps));
    //erase it from map
    nodemap.erase(node->name);
    delete node;
  }
  //by this point no node in DAG has children
  //If DAG is not empty, then there must be a cyclic dependency
  if(!nodemap.empty()){
    bool printcomma=false;
    cerr<<"[ERROR] Circular dependency between some of the following dependent parameters:";
    for(map<string,DAGnode*>::const_iterator it=nodemap.begin();it!=nodemap.end();++it){
      if(printcomma)cerr<<',';
      else printcomma=true;
      cerr<<' '<<it->first;
    }
    cerr<<';'<<endl;
    hf_errlog(18112003,"F: circular dependency of dependent parameters, see stderr");
  }
}
void updateDependentParameters(){
  for(Constraint*c:constraints)c->atIteration();
}
void resetDependentParameters(){
  for(Constraint*c:constraints)delete c;
  constraints.clear();
}
}
