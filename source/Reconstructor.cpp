#include "Reconstructor.hpp"
#include "Tests.hpp"

namespace firefly {
  Reconstructor::Reconstructor(uint32_t n_, uint32_t thr_n_, uint32_t mode_): n(n_), thr_n(thr_n_), tp(thr_n) {
    mode = mode_;
  }

  void Reconstructor::scan_for_sparsest_shift() {
    scan = true;
  }

  void Reconstructor::reconstruct() {
    std::cout << "New prime: 1\n";
    FFInt::set_new_prime(primes()[prime_it]);
    RatReconst tmp(n);

    start_probe_jobs(std::vector<uint32_t>(n - 1, 1), thr_n);
    started_probes.emplace(std::vector<uint32_t>(n - 1, 1), thr_n);

    FFInt t = 1;
    std::vector<uint32_t> zi_order(n - 1, 1);
    std::vector<FFInt> probe {};

    {
      std::unique_lock<std::mutex> lock(mut);
      find: while (jobs_finished == 0) {
        cond.wait(lock);
      }

      for (auto it = probes.begin(); it != probes.end(); it++) {
        if ((std::get<2>(*it)).wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
          t = std::get<0>(*it);
          zi_order = std::get<1>(*it);
          probe = std::move((std::get<2>(*it)).get());
          it = probes.erase(it);
          --jobs_finished;
          break;
        }
      }

      // sometimes the future is not ready even though the solution job returned already
      // probably there is an additional copy operation involved
      if (probe.empty()) {
        goto find;
      }
    }

    std::cout << "Iteration: 1\n";
    uint32_t items = 0;
    ++fed_ones;

    for (const auto& value : probe) {
      reconst.emplace_back(RatReconst(n));
        // save intermediate results
  //      reconst.set_tag(std::to_string(equation.first) + "_" + std::to_string(term.first));
      reconst.back().feed(t, value, zi_order, prime_it);
      reconst.back().interpolate();
      ++items;
    }

    probe.clear();
    std::cout << "Items: " << items << "\n";

    start_probe_jobs(std::vector<uint32_t>(n - 1, 1), 1);
    ++started_probes.at(std::vector<uint32_t>(n - 1, 1));

    uint32_t iteration = 1;
    uint32_t total_iterations = 0;

    bool done = false;
    bool new_prime = false;

    while(!done) {
      if (new_prime) {
        total_iterations += iteration;
        ++prime_it;
        std::cout << "\nNew prime: " << prime_it + 1 << "\n";
        std::cout << "Iterations for last prime: " << iteration << "\n";
        std::cout << "Iterations in total: " << total_iterations << "\n\n";
        iteration = 0;
        if (prime_it == 4) {
          for (auto& el : reconst) {
            std::cout << el.is_done() << "\n";
          }
          exit(-1);
        }

        tp.kill_all();

        probes.clear();
        jobs_finished = 0;
        started_probes.clear();
        new_prime = false;
        feeding_jobs = 0;

        FFInt::set_new_prime(primes()[prime_it]);

        // if only a small constant is reconstructed it will not ask for new run
        if (probes_for_next_prime == 0) {
          probes_for_next_prime = 1;
        }

        if (!tmp.need_shift()) {
          std::cout << "Disable shift\n";
          tmp.disable_shift();
        }
        tmp.generate_anchor_points();

        // start only coreNumber jobs first, because the reconstruction can be done after the first feed
        if (probes_for_next_prime > thr_n) {
          std::cout << "Starting " << thr_n << " jobs now, the remaining " << probes_for_next_prime - thr_n << " jobs will be started later\n";
          start_probe_jobs(std::vector<uint32_t>(n - 1, 1), thr_n);
          started_probes.emplace(std::vector<uint32_t>(n - 1, 1), thr_n);
        } else {
          std::cout << "Starting " << probes_for_next_prime << " jobs\n";
          start_probe_jobs(std::vector<uint32_t>(n - 1, 1), probes_for_next_prime);
          started_probes.emplace(std::vector<uint32_t>(n - 1, 1), probes_for_next_prime);
        }
        probes_for_next_prime = 0;
      }

      {
        std::unique_lock<std::mutex> lock(mut);
        find_loop: while (jobs_finished == 0) {
          cond.wait(lock);
        }

        for (auto it = probes.begin(); it != probes.end(); it++) {
          if ((std::get<2>(*it)).wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
            t = std::get<0>(*it);
            zi_order = std::get<1>(*it);
            probe = std::move((std::get<2>(*it)).get());
            it = probes.erase(it);
            --jobs_finished;
            break;
          }
        }

        // sometimes the future is not ready even though the solution job returned already
        // probably there is an additional copy operation involved
        if (probe.empty()) {
          goto find_loop;
        }
      }

      ++iteration;

      uint32_t items_done = 0;
      uint32_t items_new_prime = 0;
      {
        std::unique_lock<std::mutex> lock(mut);
        std::cout << "Iteration: " << iteration << "\n";
        if (prime_it == 0 && zi_order == std::vector<uint32_t>(n - 1, 1)) {
          ++fed_ones;
        }
      }

      for (uint32_t i = 0; i != reconst.size(); ++i) {
        if (!reconst[i].is_done()) {
          if (reconst[i].get_prime() == prime_it) {
            {
              std::unique_lock<std::mutex> lock(mut);
              ++feeding_jobs;
            }

            reconst[i].feed(t, probe[i], zi_order, prime_it);

            tp.run_task([this, i](){
              interpolate_job(reconst[i], i);
            });
          } else {
            ++items_new_prime;
          }
        } else {
          ++items_done;
        }
      }

      probe.clear();

      std::unique_lock<std::mutex> lock(mut);
      std::cout << "Done: " << items_done << " / " << items << "\n";
      std::cout << "Want new prime: " << items_new_prime << " / " << items << "\n";

      if (items_done == items) {
        done = true;
      } else if (items_done + items_new_prime == items) {
        new_prime = true;
      } else if (probes.size() == 0) {
        bool cont = false;

        while (feeding_jobs > 0) {
          cond.wait(lock);
          if (probes.size() > 0) {
            cont = true;
            break;
          }
        }

        if (cont) {
          continue;
        }

        // no jobs are running anymore, check if done or new_prime else throw error
        items_done = 0;
        items_new_prime = 0;

        for (uint32_t i = 0; i != reconst.size(); ++i) {
          if (!reconst[i].is_done()) {
            if (reconst[i].get_prime() != prime_it) {
              ++items_new_prime;
            }
          } else {
            ++items_done;
          }
        }

        if (items_done == items) {
          done = true;
        } else if (items_done + items_new_prime == items) {
          new_prime = true;
        } else {
          throw std::runtime_error("No items to feed anymore");
        }
      }
    }

    total_iterations += iteration;
    std::cout << "\nDone\n";
    std::cout << "Iterations for last prime: " << iteration << "\n";
    std::cout << "Primes: " << prime_it + 1 << "\n";
    std::cout << "Iterations: " << total_iterations << "\n\n";
  }

  std::vector<RationalFunction> Reconstructor::get_result_rf() {
    std::vector<RationalFunction> result {};
    for (uint32_t i = 0; i != reconst.size(); ++i) {
      result.emplace_back(reconst[i].get_result());
    }
    return result;
  }

  void Reconstructor::start_probe_jobs(const std::vector<uint32_t> zi_order, const uint32_t start) {
    RatReconst tmp(n);
    std::vector<firefly::FFInt> values(n);
    std::vector<firefly::FFInt> shift = tmp.get_zi_shift_vec();

    for (uint32_t j = 0; j != start; ++j) {
      FFInt t = tmp.get_rand();
      values[0] = t + shift[0];

      std::vector<firefly::FFInt> rand_zi = tmp.get_rand_zi_vec(zi_order);
      for (uint32_t i = 1; i != n; ++i) {
        values[i] = rand_zi[i - 1] * t + shift[i];
      }

      std::unique_lock<std::mutex> lock(mut);

      auto future = tp.run_packaged_task([this, values](){
        std::vector<FFInt> probe {};
        black_box(probe, values);

        std::unique_lock<std::mutex> lock(mut);
        ++jobs_finished;
        cond.notify_one();
        return std::move(probe);
      });

      probes.emplace_back(std::make_tuple(t, zi_order, std::move(future)));
    }
  }

  void Reconstructor::interpolate_job(RatReconst& reconst, uint32_t i) {
    reconst.interpolate();

    // start new jobs if required
    // since interpolate returns immediately if it is already interpolating, too many jobs could be started if the reconstrution is done after the first interpolating finishes
    std::unique_lock<std::mutex> lock(mut);
    if (!reconst.is_done()) {
      if (reconst.get_prime() > prime_it) {
        if (reconst.get_num_eqn() > probes_for_next_prime) {
          probes_for_next_prime = reconst.get_num_eqn();
        }
      } else {
        std::vector<uint32_t> zi_order = reconst.get_zi_order();

        if (prime_it == 0 && zi_order == std::vector<uint32_t>(n - 1, 1)) {
          if (started_probes.at(zi_order) - thr_n <= fed_ones - 1) {
            uint32_t start = fed_ones - started_probes.at(zi_order) + thr_n;
            std::cout << i << " Starting ones: " << start << "\n";
            started_probes.at(zi_order) += start;
            lock.unlock();
            start_probe_jobs(zi_order, start);
            lock.lock();
          }
        } else {
          uint32_t required_probes = reconst.get_num_eqn();
          auto it = started_probes.find(zi_order);

          if (it != started_probes.end()) {
            if (required_probes > started_probes.at(zi_order)) {
              std::cout << i << " Starting ";
              for (auto ele : zi_order) {
                std::cout << ele << " ";
              }
              uint32_t start = required_probes - started_probes.at(zi_order);
              std::cout << "-- " << start << "\n";
              started_probes.at(zi_order) = required_probes;
              lock.unlock();
              start_probe_jobs(zi_order, start);
              lock.lock();
            }
          } else {
            std::cout << i << " Starting ";
            for (auto ele : zi_order) {
              std::cout << ele << " ";
            }
            std::cout << "-- " << required_probes << "\n";
            started_probes.emplace(zi_order, required_probes);
            lock.unlock();
            start_probe_jobs(zi_order, required_probes);
            lock.lock();
          }
        }
      }
    }

    --feeding_jobs;
    cond.notify_one();
  }
}
