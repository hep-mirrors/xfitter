extern"C"{
  void init_theory_repository_();
}

#include <iostream>
#include <string>

#include "ReactionTheory.h"
#include "TheoryRepository.h"

// Headers of all theory modules:
#include "TestTheory.h"

void init_theory_repository_() {
  std::cout << "Iniitialize theory interfaces ... \n";

// Add test theory interface:
  TestTheory* test = new TestTheory();
  xfitter::TheoryRepository.insert(std::pair<std::string,ReactionTheory*>(test->getReactionName(), test));
} 
