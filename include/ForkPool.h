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

  explicit ForkPool(int ncpu)
    : ncpu_(std::max(1, ncpu))
  {}

  class SharedMemory {
  public:

    SharedMemory() = default;

    explicit SharedMemory(size_t bytes) {
      shmid_ = shmget(IPC_PRIVATE, bytes, IPC_CREAT | 0600);

      if (shmid_ < 0)
        throw std::runtime_error(
          "shmget failed: " + std::string(strerror(errno)));

      ptr_ = shmat(shmid_, nullptr, 0);

      if (ptr_ == reinterpret_cast<void*>(-1)) {
        shmctl(shmid_, IPC_RMID, nullptr);

        throw std::runtime_error(
          "shmat failed: " + std::string(strerror(errno)));
      }
    }

    SharedMemory(const SharedMemory&) = delete;
    SharedMemory& operator=(const SharedMemory&) = delete;

    SharedMemory(SharedMemory&& other) noexcept {
      shmid_ = other.shmid_;
      ptr_ = other.ptr_;

      other.shmid_ = -1;
      other.ptr_ = nullptr;
    }

    SharedMemory& operator=(SharedMemory&& other) noexcept {
      if (this != &other) {
        cleanup();

        shmid_ = other.shmid_;
        ptr_ = other.ptr_;

        other.shmid_ = -1;
        other.ptr_ = nullptr;
      }

      return *this;
    }

    ~SharedMemory() {
      cleanup();
    }

    template<typename T>
    T* data() {
      return static_cast<T*>(ptr_);
    }

  private:

    void cleanup() {
      if (ptr_ && ptr_ != reinterpret_cast<void*>(-1))
        shmdt(ptr_);

      if (shmid_ >= 0)
        shmctl(shmid_, IPC_RMID, nullptr);
    }

    int shmid_ = -1;
    void* ptr_ = nullptr;
  };

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
      for (size_t i = 0; i < N; i++)
        func(i);

      return;
    }

    int errpipe[2];

    if (pipe(errpipe) < 0) {
      throw std::runtime_error(
        "pipe failed: " + std::string(strerror(errno)));
    }

    std::vector<pid_t> pids;
    pids.reserve(nworkers);

    for (int worker = 0; worker < nworkers; worker++) {

      size_t begin = worker * N / nworkers;
      size_t end   = (worker + 1) * N / nworkers;

      pid_t pid = fork();

      if (pid < 0) {

        close(errpipe[0]);
        close(errpipe[1]);

        throw std::runtime_error(
          "fork failed: " + std::string(strerror(errno)));
      }

      if (pid == 0) {

        close(errpipe[0]);

        try {

          int fdlimit = static_cast<int>(sysconf(_SC_OPEN_MAX));

          for (int fd = STDERR_FILENO + 1; fd < fdlimit; fd++) {
            if (fd != errpipe[1])
              close(fd);
          }

          for (size_t i = begin; i < end; i++)
            func(i);

          close(errpipe[1]);

          _exit(0);
        }
        catch (const std::exception& e) {

          std::string msg =
            "pid=" + std::to_string(getpid()) +
            ": " + e.what() + "\n";

          write_all(errpipe[1], msg.data(), msg.size());

          close(errpipe[1]);

          _exit(1);
        }
        catch (...) {

          std::string msg =
            "pid=" + std::to_string(getpid()) +
            ": unknown exception\n";

          write_all(errpipe[1], msg.data(), msg.size());

          close(errpipe[1]);

          _exit(1);
        }
      }

      pids.push_back(pid);
    }

    close(errpipe[1]);

    int failed = 0;

    for (pid_t pid : pids) {

      int status;

      if (waitpid(pid, &status, 0) < 0) {

        close(errpipe[0]);

        throw std::runtime_error(
          "waitpid failed: " + std::string(strerror(errno)));
      }

      if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
        failed++;
    }

    std::string errorMessage;

    char buffer[4096];

    ssize_t nread;

    while ((nread = read(errpipe[0], buffer, sizeof(buffer))) > 0)
      errorMessage.append(buffer, nread);

    close(errpipe[0]);

    if (failed) {

      if (!errorMessage.empty()) {
        throw std::runtime_error(
          "ForkPool worker failure:\n" + errorMessage);
      }

      throw std::runtime_error(
        "ForkPool worker terminated unexpectedly");
    }
  }

  private:

  int ncpu_;

  static void write_all(int fd, const char* data, size_t size) {
    while (size > 0) {
      ssize_t n = write(fd, data, size);

      if (n <= 0)
        return;

      data += n;
      size -= n;
    }
  }
};
