#pragma once

#include <map>
#include "ReactionTheory.h"

// Map to store pointers to the modules

namespace xfitter {

static std::map<std::string, ReactionTheory*> TheoryRepository;

}
