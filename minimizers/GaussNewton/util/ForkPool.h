#pragma once

// Runs independent tasks in forked worker processes; each task returns a fixed number of
// doubles through System V shared memory. With 0 workers the tasks run in this process.

#include <functional>
#include <vector>

namespace xfitter {
namespace gn {

class ForkPool {
public:
  explicit ForkPool(int workers) : _workers(workers) {}

  /// Runs task(i, out) for i = 0..nTasks-1, out pointing to taskSize doubles.
  /// Returns the outputs, task after task (nTasks x taskSize).
  std::vector<double> run(int nTasks, int taskSize,
                          const std::function<void(int task, double* out)>& task) const;

private:
  int _workers;
};

}
}
