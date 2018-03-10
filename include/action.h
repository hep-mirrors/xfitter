#pragma once
/*!
 @file action.h
 @date 10 Mar 2018
 @author SG

  Loadable "action" to study xFitter information during different stages of execution. 

*/

#include <iostream>

class Action
{
public:
  Action(){}
  virtual ~Action() {}
  virtual void Initialize() { }
  virtual void AtIteration() {std::cout << "At Iter " << std::endl;}
  virtual void Finalize() { }
};

