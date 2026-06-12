#pragma once

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

class ForkPool {
  public:
    class SharedMemory {
    public:
      SharedMemory(size_t bytes) {
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
      SharedMemory(const SharedMemory&) = delete;
      SharedMemory& operator=(const SharedMemory&) = delete;
      SharedMemory(SharedMemory&&) = default;
      SharedMemory& operator=(SharedMemory&&) = default;

      ~SharedMemory() {
        if (ptr_ && ptr_ != reinterpret_cast<void*>(-1))
          shmdt(ptr_);
        if (shmid_ >= 0)
          shmctl(shmid_, IPC_RMID, nullptr);
      }

      template<typename T>
      T* data() {return static_cast<T*>(ptr_);}

    private:
      int shmid_ = -1;
      void* ptr_ = nullptr;
    };

    ForkPool(int ncpu) : ncpu_(ncpu) {}

    template<typename T>
    SharedMemory make_shared(size_t n) {
      return SharedMemory(sizeof(T) * n);
    }

    template<typename Func>
    void parallel_for(size_t N, Func func) const {
      if (N == 0)
        return;
      int nworkers = std::min<int>(ncpu_, N);
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
        printf("worker = %d begin = %lu end = %lu\n", worker, begin, end);
        pid_t pid = xfitter::xf_fork(std::min(ncpu_, int(N)));
        if (pid < 0) {
          hf_errlog(2026061201,"F: ForkPool failed to fork");
        }
        if (pid == 0) {
          // close all open files (e.g. minuit.out.txt) to avoid multiple buffered output
          int fdlimit = (int)sysconf(_SC_OPEN_MAX);
          for (int i = STDERR_FILENO + 1; i < fdlimit; i++) {
            close(i);
          }
          for (size_t i = begin; i < end; i++) {
            func(i);
          }
          _exit(0);
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

  private:
    int ncpu_;
};
