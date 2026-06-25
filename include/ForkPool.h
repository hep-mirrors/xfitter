#pragma once

#include <functional>

class ForkPool {
  public:
    enum class TaskDistribution
    {
      chunky,
      cyclic,
    };
    class SharedMemory {
    public:
      SharedMemory(size_t bytes);
      SharedMemory(const SharedMemory&) = delete;
      SharedMemory& operator=(const SharedMemory&) = delete;
      SharedMemory(SharedMemory&&) = default;
      SharedMemory& operator=(SharedMemory&&) = default;
      ~SharedMemory();

      template<typename T>
      T* data() {return static_cast<T*>(ptr_);}

    private:
      int shmid_ = -1;
      void* ptr_ = nullptr;
    };

    ForkPool(int ncpu, TaskDistribution distr = TaskDistribution::chunky);

    void parallel_for(size_t N, const std::function<void(size_t)>& func) const;

    static TaskDistribution getTaskDistrFromStr(const std::string& str);

  private:
    const int _ncpu;
    const TaskDistribution _distr;
};
