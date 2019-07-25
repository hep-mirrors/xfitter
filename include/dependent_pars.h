#pragma once
#include<string>
#include<vector>
//We call a parameter dependent when it is expressed in terms of some other expression, e.g.
//A=(2*B+1)^2
//Function updateDependentParameters evaluates such expressions and sets values of all dependent parameters
//Expression is provided as a string
//We rely on "TinyExpr" library for parsing and evaluating expressions
namespace xfitter{
  struct DependentParameter{
    std::string name;//e.g. "Av"
    std::string expression;//e.g. "= B-4.0*C^2 +.5*E"
  };
  void registerDependentParameters(const std::vector<DependentParameter>&dependentParameters);//Call only once
  void updateDependentParameters();//Shall be called at each iteration, before parameterisations
  void resetDependentParameters();//Call at end to release memory (nobody calls it currently)
}
