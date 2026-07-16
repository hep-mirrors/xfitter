#include "ForkPool.h"
#include <sys/shm.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>
#include "xfitter_cpp_base.h"
#include "xfitter_steer.h"

ForkPool::SharedMemory::SharedMemory(size_t bytes) {
  shmid_ = shmget(IPC_PRIVATE, bytes, IPC_CREAT | 0600);
  if (shmid_ < 0) {
    hf_errlog(2026061202,"F: ForkPool failed to create shared memory segment");
  }
  ptr_ = shmat(shmid_, nullptr, 0);
  if (ptr_ == reinterpret_cast<void*>(-1)) {
    shmctl(shmid_, IPC_RMID, nullptr);
    hf_errlog(2026061202,"F: ForkPool failed to attach shared memory segment");
  }
}

ForkPool::SharedMemory::~SharedMemory() {
  if (ptr_ && ptr_ != reinterpret_cast<void*>(-1))
    shmdt(ptr_);
  if (shmid_ >= 0)
    shmctl(shmid_, IPC_RMID, nullptr);
}

ForkPool::ForkPool(int ncpu, TaskDistribution distr) : _ncpu(ncpu), _distr(distr) {
  if(ncpu <= 0) {
    hf_errlog(2026062201,"F: ForkPool got invalid ncpu = " + std::to_string(ncpu));
  }
}

void ForkPool::parallel_for(size_t N, const std::function<void(size_t)>& func) const {
  if (N == 0) {
    return;
  }
  int nworkers = std::min<int>(_ncpu, N);
  if (nworkers == 1) {
    for (size_t i = 0; i < N; i++) {
      func(i);
    }
    return;
  }
  std::vector<pid_t> pids;
  pids.reserve(nworkers);
  for (int worker = 0; worker < nworkers; worker++) {
    size_t begin = worker * N / nworkers;
    size_t end   = (worker + 1) * N / nworkers;
    //if(_distr == TaskDistribution::chunky) {
    //  printf("worker = %d [%lu jobs] begin = %lu end = %lu (chunky)\n", worker, end - begin, begin, end);
    //}
    //else if(_distr == TaskDistribution::cyclic) {
    //  printf("worker = %d begin = %d step = %d (cyclic)\n", worker, worker, nworkers);
    //}
    //setbuf(stdout, nullptr);
    fflush(nullptr);
    //pid_t pid = fork();
    pid_t pid = xfitter::xf_fork(std::min(_ncpu, int(N)));
    if (pid < 0) {
      hf_errlog(2026061201,"F: ForkPool failed to fork");
    }
    if (pid == 0) {
      // close all open files (e.g. minuit.out.txt) to avoid multiple buffered output
      //int fdlimit = (int)sysconf(_SC_OPEN_MAX);
      //printf("fdlimit = %d\n", fdlimit);
      //for (int i = STDERR_FILENO + 1; i < fdlimit; i++) {
      //  close(i);
      //}
      if(_distr == TaskDistribution::chunky) {
        for (size_t i = begin; i < end; i++) {
          func(i);
        }
      }
      else if(_distr == TaskDistribution::cyclic) {
        for (size_t i = worker; i < N; i += nworkers) {
          func(i);
        }
      }
      fflush(nullptr);
      _exit(0);
      //exit(0);
    }
    pids.push_back(pid);
  }
  int failed = 0;
  for (pid_t pid : pids) {
    int status;
    if (waitpid(pid, &status, 0) < 0) {
      hf_errlog(2026061204,"F: ForkPool waitpid failed");
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
      failed++;
    }
  }
  if (failed) {
    hf_errlog(2026061205,"F: ForkPool worker failure");
  }
}

ForkPool::TaskDistribution ForkPool::getTaskDistrFromStr(const std::string& str) {
  if(str == "chunky") {
    return TaskDistribution::chunky;
  }
  else if(str == "cyclic") {
    return TaskDistribution::cyclic;
  }
  else {
    auto msg = "F: unknown task distribution " + str;
    hf_errlog(2026061601,msg);
  }
  return TaskDistribution::chunky; // avoid warning
}
