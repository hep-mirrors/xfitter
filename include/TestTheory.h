#pragma once

// Test simple class to use new theory interface

#include "ReactionTheory.h"
#include <iostream>


class TestTheory : public ReactionTheory {
 public:
  virtual string getReactionName() const 
  {
     std::cout << " Init test reaction " << std::endl;
     return "test"; 
  };
  virtual int compute(valarray<double> &val, map<string, valarray<double> > &err) 
  {   
    std::cout << " Evaluate prediction " << std::endl;
    return 0;
  };

 protected:

};
