#pragma once

#include <string>
/**
     @class xfitter_steer

     @brief basic building blocks of xFitter running sequence

     @version 0.1
     @date 2018-07-13
*/

namespace xfitter
{
  /// Load named evolution code.  
  void load_evolution(std::string name="");  
}
