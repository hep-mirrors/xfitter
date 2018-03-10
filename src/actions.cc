#include "actions.h"

namespace  XFITTER_ACTION {
  list <Action*> gActionsList;
}

void actions_initialize_() {
  for ( auto action : XFITTER_ACTION::gActionsList) {
    action->Initialize();
  }
}

void actions_at_iteration_() {
  for ( auto action : XFITTER_ACTION::gActionsList) {
    action->AtIteration();
  }
}

void actions_finalize_() {
  for ( auto action : XFITTER_ACTION::gActionsList) {
    action->Finalize();
  }
}
