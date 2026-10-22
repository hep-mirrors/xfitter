#pragma once

// Accumulating wall-clock timer for the iteration log

#include <chrono>

namespace xfitter {
namespace gn {

class Timer {
public:
  /// Adds the time from construction to destruction to an accumulator
  class Scope {
  public:
    explicit Scope(double& total) : _total(total), _start(std::chrono::steady_clock::now()) {}
    ~Scope() {
      _total += std::chrono::duration<double>(std::chrono::steady_clock::now() - _start).count();
    }

  private:
    double& _total;
    std::chrono::steady_clock::time_point _start;
  };
};

}
}
