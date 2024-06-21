#pragma once

#include <vector>

namespace xfitter
{
  /**
     @class PartonSens

     @brief Class to calculate sensitivity to PDFs as function of x

     @version 0.1
     @date 2024-06-17
  */

  class PartonSens
  {
  public:
    /// parse yaml, perform analysis
    void doPartonSens();
    /// get partons
    static std::vector<int> getPartons(int i=1) {
      if(i == 1) return _global_partons;
      else return _global_partons1;
    }
    /// get x1
    static double getX1() { return _global_x1; }
    /// get x2
    static double getX2() { return _global_x2; }
  private:
    static std::vector<int> _global_partons;
    static std::vector<int> _global_partons1;
    static double _global_x1;
    static double _global_x2;
  };
  
} //namespace xfitter
