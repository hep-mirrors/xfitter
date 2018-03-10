#pragma once

/*!
 @file actions.h
 @date 10 Mar 2018
 @author SG

  Loadable "actions" to study xFitter information during different stages of execution. 
  Contains python interface
*/

#include <list>
#include "action.h"

using std::list;

namespace XFITTER_ACTION {
  extern list <Action*> gActionsList;
}

// Fortran interface:
extern "C" {
  void actions_initialize_();
  void actions_at_iteration_();
  void actions_finalize_();
}

// Add action to the list
void addAction(Action& a) {
  XFITTER_ACTION::gActionsList.push_back(&a);
}
