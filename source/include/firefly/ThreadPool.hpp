// ====================================================================
// This file is part of FlexibleSUSY.
// Modified work Copyright (C) 2020  Jonas Klappert and Fabian Lange
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#pragma once

#include "firefly/config.hpp"
#include "firefly/Logger.hpp"

#include <condition_variable>
#include <deque>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace firefly {
  /**
   * @class ThreadPool
   * @brief A pool of threads
   *
   * Thread_pool represents a collection of threads. Tasks (callables)
   * can be added to an internal deque. The tasks will be executed as
   * soon as there is an idle thread. The destructor of the Thread_pool
   * will wait until all tasks are finished and the deque is empty.
   */
  class ThreadPool {
  public:
    /**
     * Constructor of the Thread_pool
     *
     * @param pool_size number of threads in the pool
     */
    explicit ThreadPool(std::size_t pool_size = std::thread::hardware_concurrency()) {

      std::unique_lock<std::mutex> lock(mutex);
      for (std::size_t i = 0; i < pool_size; ++i) {
        threads_idle.emplace_back(true);
        threads.emplace_back(
        [this, i]() {
          for (;;) {
            std::function<void(uint32_t thread_number)> task;

            {
              std::unique_lock<std::mutex> lock(mutex);
              threads_idle[i] = true;
              condition_wait.notify_one();
              condition.wait(lock, [this] { return stop || !tasks.empty(); });
              threads_idle[i] = false;

              if (stop && tasks.empty())
                return;

              task = std::move(tasks.front());
              tasks.pop_front();
            }

            task(i);
          }
        });
      }
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool(ThreadPool &&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool& operator=(ThreadPool &&) = delete;

    // waits for all tasks to finish and closes threads
    ~ThreadPool() {
      {
        std::lock_guard<std::mutex> lock(mutex);
        stop = true;
      }

      condition.notify_all();

      try {
        for (auto & t : threads)
          t.join();
      } catch (const std::exception& e) {
        ERROR_MSG(e.what());
      }
    }

    // runs task and returns future
    template <typename Task>
    auto run_packaged_task(Task && task) -> std::future<decltype(task())> {
      using return_t = decltype(task());

      auto ptask = std::make_shared<std::packaged_task<return_t()>>([task](uint32_t thread_id) { return task(thread_id); });

      std::future<return_t> fut = ptask->get_future();

      if (threads.empty()) {
        (*ptask)(0);
      } else {
        {
          std::lock_guard<std::mutex> lock(mutex);
          tasks.emplace_back([ptask](uint32_t thread_id) { (*ptask)(thread_id); });
        }
        condition.notify_one();
      }

      return fut;
    }

    template <typename Task>
    auto run_priority_packaged_task(Task && task) -> std::future<decltype(task())> {
      using return_t = decltype(task());

      auto ptask = std::make_shared<std::packaged_task<return_t()>>([task](uint32_t thread_id) { return task(thread_id); });

      std::future<return_t> fut = ptask->get_future();

      if (threads.empty()) {
        (*ptask)(0);
      } else {
        {
          std::lock_guard<std::mutex> lock(mutex);
          tasks.emplace_front([ptask](uint32_t thread_id) { (*ptask)(thread_id); });
        }
        condition.notify_one();
      }

      return fut;
    }

    // runs task
    template <typename Task>
    void run_task(Task && task) {
      if (threads.empty()) {
        task(0);
      } else {
        {
          std::lock_guard<std::mutex> lock(mutex);
          tasks.emplace_back(std::forward<Task>(task));
        }
        condition.notify_one();
      }
    }

    template <typename Task>
    void run_priority_task(Task && task) {
      if (threads.empty()) {
        task(0);
      } else {
        {
          std::lock_guard<std::mutex> lock(mutex);
          tasks.emplace_front(std::forward<Task>(task));
        }
        condition.notify_one();
      }
    }

    std::size_t size() const { return threads.size(); }

    std::size_t queue_size() { std::lock_guard<std::mutex> lock(mutex); return tasks.size(); }

    // if some threads are working, waits until one finishes and returns true; if all threads are idle, returns false
    bool wait() {
      std::unique_lock<std::mutex> lock(mutex);

      if (!all_threads_idle()) {
        condition_wait.wait(lock);
        return true;
      } else {
        return false;
      }
    }

    void kill_all() {
      std::unique_lock<std::mutex> lock(mutex);

      tasks = std::deque<std::function<void(uint32_t)>>();

      while (!all_threads_idle()) {
        tasks = std::deque<std::function<void(uint32_t)>>();

        condition_wait.wait(lock);
      }

      tasks = std::deque<std::function<void(uint32_t)>>();
    }

  private:
    std::vector<std::thread> threads {};
    std::deque<std::function<void(uint32_t)>> tasks {};
    std::mutex mutex {};
    std::condition_variable condition {};
    bool stop {false};
    std::vector<bool> threads_idle {};
    std::condition_variable condition_wait {};

    // Should only be called when the mutex is locked
    bool all_threads_idle() {
      if (threads.empty()) {
        if (tasks.size() > 0) {
          return false;
        } else {
          return true;
        }
      } else {
        for (auto idle : threads_idle) {
          if (!idle) {
            return false;
          }
        }

        return true;
      }
    }
  };
}
