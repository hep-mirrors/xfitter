#include "util/ForkPool.h"
#include "xfitter_cpp_base.h"
#include "xfitter_steer.h"

#include <algorithm>
#include <cstdio>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/wait.h>
#include <unistd.h>

namespace xfitter {
namespace gn {

namespace {

/// Shared-memory segment, detached and removed on destruction
class SharedBuffer {
public:
  explicit SharedBuffer(size_t size) {
    _id = shmget(IPC_PRIVATE, std::max<size_t>(size, 1) * sizeof(double), IPC_CREAT | 0600);
    if (_id < 0) hf_errlog(2026102220, "F: GaussNewton: failed to create shared memory");
    void* address = shmat(_id, nullptr, 0);
    if (address == reinterpret_cast<void*>(-1)) hf_errlog(2026102221, "F: GaussNewton: failed to attach shared memory");
    _data = static_cast<double*>(address);
  }
  ~SharedBuffer() {
    shmdt(_data);
    shmctl(_id, IPC_RMID, nullptr);
  }
  SharedBuffer(const SharedBuffer&) = delete;
  SharedBuffer& operator=(const SharedBuffer&) = delete;
  double* data() { return _data; }

private:
  int _id;
  double* _data;
};

}  // namespace

std::vector<double> ForkPool::run(int nTasks, int taskSize,
                                  const std::function<void(int, double*)>& task) const {
  std::vector<double> out(static_cast<size_t>(nTasks) * taskSize);
  const int workers = std::min(xf_ncpu(_workers), nTasks);
  if (workers <= 1) {
    for (int i = 0; i < nTasks; ++i) task(i, out.data() + static_cast<size_t>(i) * taskSize);
    return out;
  }

  SharedBuffer shared(out.size());
  std::fflush(nullptr);  // do not duplicate buffered output in the children
  std::vector<pid_t> children;
  for (int w = 0; w < workers; ++w) {
    const pid_t pid = xf_fork(workers);
    if (pid < 0) hf_errlog(2026102222, "F: GaussNewton: fork failed");
    if (pid == 0) {
      for (int i = w; i < nTasks; i += workers)
        task(i, shared.data() + static_cast<size_t>(i) * taskSize);
      std::fflush(nullptr);
      _exit(0);  // no atexit handlers or destructors of the parent state in the child
    }
    children.push_back(pid);
  }
  bool ok = true;
  for (pid_t pid : children) {
    int status = 0;
    if (waitpid(pid, &status, 0) != pid || !WIFEXITED(status) || WEXITSTATUS(status) != 0) ok = false;
  }
  if (!ok) hf_errlog(2026102223, "F: GaussNewton: a worker process failed");
  std::copy(shared.data(), shared.data() + out.size(), out.begin());
  return out;
}

}
}
