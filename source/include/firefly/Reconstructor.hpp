//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2020  Jonas Klappert and Fabian Lange
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//==================================================================================

#pragma once

#include "firefly/config.hpp"
#include "firefly/BlackBoxBase.hpp"
#include "firefly/FFIntVec.hpp"
#include "firefly/gzstream.hpp"
#include "firefly/ParserUtils.hpp"
#include "firefly/RatReconst.hpp"
#include "firefly/ReconstHelper.hpp"
#include "firefly/ThreadPool.hpp"
#include "firefly/tinydir.h"
#include "firefly/utils.hpp"
#include "firefly/version.hpp"

#ifdef WITH_MPI
#include "firefly/MPIWorker.hpp"
#endif

#include <flint/fmpz_poly.h>

#include <chrono>
#include <numeric>
#include <tuple>
#include <algorithm>
#include <sys/stat.h>

namespace firefly {
  typedef std::tuple<uint64_t, std::atomic<int>, RatReconst*> RatReconst_tuple;
  typedef std::list<RatReconst_tuple> RatReconst_list;

  /**
   * @class Reconstructor
   * @brief A class to reconstruct the functions in BlackBoxTemp from its values
   */
  template<typename BlackBoxTemp>
  class Reconstructor {
  public:
    /**
     *  A constructor for the Reconstructor class
     *  @param n_ the number of parameters
     *  @param thr_n_ the number of threads being used during the reconstruction
     *  @param bb_ An instance of a BlackBoxBase class
     *  @param verbosity_ the verbosity level which can be chosen as SILENT (no output), IMPORTANT (only important output), and CHATTY (everything)
     */
    Reconstructor(const uint32_t n_, const uint32_t thr_n_, BlackBoxBase<BlackBoxTemp>& bb_, const int verbosity_ = IMPORTANT);
    /**
     *  A constructor for the Reconstructor class
     *  @param n_ the number of parameters
     *  @param thr_n_ the number of threads being used during the reconstruction
     *  @param bunch_size_ the bunch size
     *  @param bb_ An instance of a BlackBoxBase class
     *  @param verbosity_ the verbosity level which can be chosen as SILENT (no output), IMPORTANT (only important output), and CHATTY (everything)
     */
    Reconstructor(const uint32_t n_, const uint32_t thr_n_, const uint32_t bunch_size_, BlackBoxBase<BlackBoxTemp>& bb_, const int verbosity_ = IMPORTANT);
    /**
     *  A destructor for the Reconstructor class
     */
    ~Reconstructor();
    /**
     *  Default constructor
     */
    Reconstructor() {};
    /**
     *  Starts the reconstruction
     *  @param prime_counter sets how many interpolations have to be performed at most
     */
    void reconstruct(uint32_t prime_counter = 300);
    /**
     *  @param vars the vector of variables as strings
     *  @return the vector of factors for all functions that shall be reconstructed
     */
    std::vector<std::string> get_factors_string(const std::vector<std::string>& vars);
    /**
     *  @return the vector of reconstructed rational functions
     */
    std::vector<RationalFunction> get_result();
    /**
     *  @return the vector of inteprolated rational functions over the last field
     */
    std::vector<RationalFunctionFF> get_result_ff();
    /**
     *  @return a vector of the already reconstructed rational functions and its tag;
     *  they are removed from the internal memory afterwards
     */
    std::vector<std::pair<std::string, RationalFunction>> get_early_results();
    /**
     *  @deprecated
     *  Enables the scan for a sparse shift at the beginning of the reconstruction
     */
    [[deprecated("Deprecated function. Use 'enable_shift_scan' instead.")]]
    void enable_scan();
    /**
     *  Enables the scan for a sparse shift at the beginning of the reconstruction
     */
    void enable_shift_scan();
    /**
     *  Enables the scan for factors
     */
    void enable_factor_scan();
    /**
     *  Activate the safe interpolation mode where the function is completely interpolated in each prime field,
     *  no optimizations are used after the first prime field. Note that this mode cannot handle function changes
     *  which lead to coefficients which will become zero in all but one prime field.
     * < b>This option has to be set before calling resume_from_saved_state().< /b>
     */
    void set_safe_interpolation();
    /**
     *  Sets user defined tags for each reconstruction object and saves intermediate results after each prime field
     *  @param tags_ a vector of user defined tags in an immutable ordering
     */
    void set_tags(const std::vector<std::string>& tags_);
    /**
     *  Sets default tags for each reconstruction object and saves intermediate results after each prime field
     */
    void set_tags();
    /**
     *  @deprecated
     *  Resumes the reconstruction of a set of given functions
     *  @param file_paths_ a vector to the absolute paths to the intermediate results of reconstruction objects
     */
    void resume_from_saved_state(const std::vector<std::string>& file_paths_);
    /**
     *  Resumes the reconstruction of a set of functions which are located in a directory.
     *  The corresponding interpolation objects are created in the same order as they were defined in the prior
     *  run, thus requiring the black box to be probed in the same order.
     */
    void resume_from_saved_state();
    /**
     *  TODO
     */
    void set_anchor_points(const std::vector<std::vector<std::uint64_t>>& anchor_points);
    /**
     *  TODO
     */
    void set_shifts(const std::vector<std::vector<std::uint64_t>>& shifts);
    /**
     *  TODO
     */
    void load_precomputed_probes();
    /**
     *  TODO
     */
    bool reconstruction_done();
    /*
     *  Enables only the scan for factors and returns afterwards
     */
    void stop_after_factor_scan();

    enum verbosity_levels {SILENT, IMPORTANT, CHATTY};
    enum RatReconst_status {RECONSTRUCTING, DONE, DELETE};
  private:
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point prime_start = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point last_print_time = std::chrono::high_resolution_clock::now();
    const uint32_t n;
    const uint32_t thr_n;
    const uint32_t bunch_size = 1;
    uint32_t prime_it = 0;
    uint32_t prime_it_fac = 0;
    uint32_t total_iterations = 0;
    uint32_t iteration = 0;
    uint32_t probes_queued = 0;
    uint32_t balance_of_ones = 0;
    uint32_t probes_for_next_prime = 0;
    uint32_t items = 0;
    uint32_t items_done = 0;
    uint32_t items_new_prime = 0;
    uint32_t feed_jobs = 0;
    uint32_t interpolate_jobs = 0;
    uint64_t ind = 0;
    BlackBoxBase<BlackBoxTemp>& bb;
    int verbosity;
    double average_black_box_time = 0.;
    std::atomic<uint32_t> probes_fed = {0};
    std::atomic<bool> scan = {false};
    std::atomic<bool> factor_scan = {false};
    std::atomic<bool> aborted = {false};
    std::atomic<bool> resumed = {false};
    std::atomic<bool> new_prime = {false};
    std::atomic<bool> done = {false};
    std::atomic<bool> change_var_order = {false};
    static bool printed_logo;
    bool save_states = false;
    bool resume_from_state = false;
    bool safe_mode = false;
    bool one_done = false;
    bool one_new_prime = false;
    bool load_anchor_points = false;
    bool scanned_factors = false;
    bool first_print = true;
    bool stop_after_factors = false;
    bool precomputed_probes = false;
    RatReconst_list reconst;
    std::vector<std::string> tags;
    std::vector<std::string> file_paths;
    std::string curr_var = "";
    std::vector<FFInt> rand_zi_fac {};
    std::ofstream logger;
    std::vector<uint32_t> max_degs {};
    ThreadPool tp;
    // TODO tidy up the mutexes
    std::mutex future_control;
    // average_black_box_time, computed_probes, iteration
    std::mutex job_control;
    // balance_of_ones, started_probes, probes_for_next_prime
    std::mutex feed_control;
    // interpolate_jobs, feed_jobs
    std::mutex print_control;
    // protects nothing, just makes sure that the messages are printed after each other
    std::mutex status_control;
    // items_done, items_new_prime, done, new_prime, one_done, one_new_prime
    std::mutex mutex_probe_queue;
    // continue_communication, proceed, index_map, requested_probes, new_jobs, probes_queued, ind
    std::mutex clean;
    // reconst (only outside access and clean_reconst)
    std::mutex chosen_mutex;
    // chosen_t
    std::condition_variable condition_future;
    std::condition_variable condition_feed;
    std::vector<std::vector<std::string>> factorizations {};
    std::unordered_set<uint32_t> possible_factors_bb_counter {};
    std::unordered_map<uint32_t, std::list<RationalFunction>> factors_rf {};
    std::unordered_map<uint32_t, std::pair<std::list<uint32_t>, std::list<uint32_t>>> factors_degs {};
    std::unordered_map<uint64_t, std::pair<FFInt, std::vector<uint32_t>>> index_map;
    std::unordered_map<std::vector<uint32_t>, uint32_t, UintHasher> started_probes;
    std::deque<std::pair<uint64_t, std::vector<FFInt>>> requested_probes;
    std::queue<std::pair<std::vector<uint64_t>, std::vector<std::vector<FFInt>>>> computed_probes;
    std::unordered_map<std::vector<uint32_t>, std::unordered_set<uint64_t>, UintHasher> chosen_t;
    std::unordered_map<uint32_t, uint32_t> optimal_var_order {}; /**< first is old position, second is new position */
    std::unordered_map<uint32_t, ShuntingYardParser> parsed_factors {};
    std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> max_deg_map_complete {};
    RatReconst tmp_rec;
    std::vector<FFInt> shift;
    std::vector<std::vector<std::uint64_t>> anchor_points {};
    std::vector<std::vector<std::uint64_t>> shifts {};

    /**
    *  Scan the black-box functions for a sparse shift
    */
    void scan_for_shift();
    /**
    *  Scan the black-box functions for a sparse shift
    */
    void scan_for_factors();
    /**
     *  Combines results after factor scan in one prime field
     *  @param poly polynomial in FLINT`s notation
     *  @param combined_ci the map of combined coefficients
     *  @param combined_prime previously used combined prime
     *  @return the combined prime
     */
    fmpzxx combine_primes(const std::unordered_map<uint32_t, uint64_t>& poly,
                          std::unordered_map<uint32_t, fmpzxx>& combined_ci,
                          const fmpzxx& combined_prime);
    /**
     *  Initializes vector of reconstruction objects and starts first probes
     *  @param first bool to indicate whether this functions is called for the first time
     */
    void start_first_runs(bool first = true);
    /**
     *  Queue new probes with zi_order one
     */
    void queue_new_ones();
    /**
     *  Starts new jobs until the reconstruction is done
     *  @param prime_counter sets how many interpolations have to be performed at most
     */
    void run_until_done(uint32_t prime_counter = 300);
    /**
     *  Queues a number of probes for a given zi_order
     *  @param zi_order the order of which a given number of probes should be queued
     *  @param to_start the number of probes which should be queued
     *  @param first bool to indicate whether this functions is called for the first time
     */
#ifndef WITH_MPI
    void queue_probes(const std::vector<uint32_t>& zi_order, const uint32_t to_start);
#else
    void queue_probes(const std::vector<uint32_t>& zi_order, const uint32_t to_start, const bool first = false);
#endif
    /**
     *  Gets a vector of probes and corresponding indices
     *  @param indices vector of indices corresponding to the probes
     *  @param probes vector of probes
     */
    void get_probe(std::vector<uint64_t>& indices, std::vector<std::vector<FFInt>>& probes);
    /**
     *  Feeds the reconstruction objects with the probes
     *  @param indices vector of indices corresponding to the probes
     *  @param probes vector of probes to be fed
     */
    void feed_job(const std::vector<uint64_t>& indices, const std::vector<std::vector<FFInt>>& probes);
    /**
     *  Interpolates a RatReconst and queues new jobs if required
     *  @param rec a reference to the RatReconst
     */
    void interpolate_job(RatReconst_tuple& rec);
    /**
     *  Removes all RatReconst from reconst which are flagged as DELETE
     */
    void clean_reconst();
    /**
     *  Checks whether there are probes requested and starts the calculation
     */
    void get_job(uint32_t thread_id);
    /**
     *  Computes a probe
     *  @param lock_probe_queue the locked mutex of the queue of probes
     */
    void compute_probe(std::unique_lock<std::mutex>& lock_probe_queue, uint32_t thread_id);
    /**
     *  Computes N probes
     *  @param lock_probe_queue the locked mutex of the queue of probes
     */
    template<uint32_t N>
    void compute_probe(std::unique_lock<std::mutex>& lock_probe_queue, uint32_t thread_id);
    /**
     *  Resets all variables when the prime changes
     */
    void reset_new_prime();
    /**
     *  Attempts to continue the calculation if no probes are queued anymore
     */
    void attempt_to_continue();
    /**
     *  TODO
     */
    void load_precomputed_probes_from_file();
    /**
     *  TODO
     */
    void write_requested_probes_to_file();
#ifdef WITH_MPI
    int world_size = -1;
    uint32_t worker_thread_count = 0;
    uint32_t iterations_on_this_node = 0;
    bool proceed = false;
    bool continue_communication = false;
    std::unordered_map<int, uint64_t> nodes;
    std::queue<std::pair<int, uint64_t>> empty_nodes;
    std::condition_variable cond_val;
    std::atomic<bool> new_jobs = {false};
    /**
     *  Performs the setup for MPI
     */
    inline void mpi_setup();
    /**
     *  Send the first probes to the MPIWorker objects
     */
    void send_first_jobs();
    /**
     *  Communicates with the MPIWorker objects
     */
    void mpi_communicate();
#endif
  };

  template<typename BlackBoxTemp>
  bool Reconstructor<BlackBoxTemp>::printed_logo = false;

  template<typename BlackBoxTemp>
  Reconstructor<BlackBoxTemp>::Reconstructor(const uint32_t n_, const uint32_t thr_n_, BlackBoxBase<BlackBoxTemp>& bb_,
#ifndef WITH_MPI
                               const int verbosity_): n(n_), thr_n(thr_n_), bb(bb_), verbosity(verbosity_), tp(thr_n_) {
#else
                               int verbosity_): n(n_), thr_n(thr_n_ - 1), bb(bb_), verbosity(verbosity_), tp(thr_n) {
#endif
    if (n == 0) {
      ERROR_MSG("FireFly does not support functional reconstructions without variables!\n               If you want to continue, set n at least to 1.");
      logger << "FireFly does not support functional reconstructions without variables!\nIf you want to continue, set n at least to 1.\n";
      std::exit(EXIT_FAILURE);
    }

    FFInt::set_new_prime(primes()[prime_it]);
    bb.prime_changed_internal();
    uint64_t seed = static_cast<uint64_t>(std::time(0));
    BaseReconst().set_seed(seed);
    tmp_rec = RatReconst(n);

    logger.open("firefly.log");

    if (verbosity > SILENT) {
      if(!printed_logo) {
        std::cerr << "\nFire\033[1;32mFly\033[0m " << FireFly_VERSION_MAJOR << "."
                  << FireFly_VERSION_MINOR << "." << FireFly_VERSION_RELEASE << "\n\n";
        printed_logo = true;
      }
      INFO_MSG("Launching " << thr_n << " thread(s) with maximum bunch size 1");
      INFO_MSG("Using seed " + std::to_string(seed) + " for random numbers");
      logger << "\nFireFly " << FireFly_VERSION_MAJOR << "."
                << FireFly_VERSION_MINOR << "." << FireFly_VERSION_RELEASE << "\n\n"
                <<"Launching " << thr_n << " thread(s) with maximum bunch size 1\n"
                <<"Using seed " << std::to_string(seed) << " for random numbers\n";
    }
  }

  template<typename BlackBoxTemp>
  Reconstructor<BlackBoxTemp>::Reconstructor(const uint32_t n_, const uint32_t thr_n_, const uint32_t bunch_size_,
#ifndef WITH_MPI
                               BlackBoxBase<BlackBoxTemp>& bb_, const int verbosity_): n(n_), thr_n(thr_n_), bunch_size(bunch_size_), bb(bb_), verbosity(verbosity_), tp(thr_n_) {
#else
                               BlackBoxBase<BlackBoxTemp>& bb_, int verbosity_): n(n_), thr_n(thr_n_ - 1), bunch_size(bunch_size_), bb(bb_), verbosity(verbosity_), tp(thr_n) {
#endif
    if (n == 0) {
      ERROR_MSG("FireFly does not support functional reconstructions without variables!\n               If you want to continue, set n at least to 1.");
      logger << "FireFly does not support functional reconstructions without variables!\nIf you want to continue, set n at least to 1.\n";
      std::exit(EXIT_FAILURE);
    }

    if (bunch_size != 1 && bunch_size != 2 && bunch_size != 4 && bunch_size != 8 && bunch_size != 16 && bunch_size != 32 && bunch_size != 64 && bunch_size != 128) {
      ERROR_MSG("Maximum bunch size " + std::to_string(bunch_size) + " is no supported power of 2!\n               Choose among 1, 2, 4, 8, 16, 32, 64, 128");
      logger << "Maximum bunch size " << std::to_string(bunch_size) << " is no supported power of 2!\nChoose among 1, 2, 4, 8, 16, 32, 64, 128\n";
      std::exit(EXIT_FAILURE);
    }

    FFInt::set_new_prime(primes()[prime_it]);
    bb.prime_changed_internal();
    uint64_t seed = static_cast<uint64_t>(std::time(0));
    BaseReconst().set_seed(seed);
    tmp_rec = RatReconst(n);

    logger.open("firefly.log");

    if (verbosity > SILENT) {
            if(!printed_logo) {
        std::cerr << "\nFire\033[1;32mFly\033[0m " << FireFly_VERSION_MAJOR << "."
                  << FireFly_VERSION_MINOR << "." << FireFly_VERSION_RELEASE << "\n\n";
        printed_logo = true;
      }
      INFO_MSG("Launching " << thr_n << " thread(s) with maximum bunch size " + std::to_string(bunch_size_));
      INFO_MSG("Using seed " + std::to_string(seed) + " for random numbers");
      logger << "\nFireFly " << FireFly_VERSION_MAJOR << "."
                << FireFly_VERSION_MINOR << "." << FireFly_VERSION_RELEASE << "\n\n"
                <<"Launching " << thr_n << " thread(s) with maximum bunch size " << std::to_string(bunch_size_) <<  "\n"
                <<"Using seed " + std::to_string(seed) + " for random numbers\n";
    }
  }

  template<typename BlackBoxTemp>
  Reconstructor<BlackBoxTemp>::~Reconstructor() {
    logger.close();
    tp.kill_all();

    auto it = reconst.begin();

    while (it != reconst.end()) {
      // delete RatReconst
      delete std::get<2>(*it);

      // remove from list
      it = reconst.erase(it);
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::stop_after_factor_scan() {
    stop_after_factors = true;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::enable_scan() {
    enable_shift_scan();
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::enable_shift_scan() {
    if (n == 1) {
      WARNING_MSG("Shift scan disabled for a univariate rational function.");
      logger << "Shift scan disabled for a univariate rational function.\n";
    } else {
      scan = true;
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::enable_factor_scan() {
    factor_scan = true;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::set_tags() {
    save_states = true;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::set_tags(const std::vector<std::string>& tags_) {
    save_states = true;
    tags = tags_;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::resume_from_saved_state() {
    tinydir_dir dir;
    tinydir_open_sorted(&dir, "ff_save/states");

    std::vector<std::string> files;
    std::vector<std::string> paths;

    for (size_t i = 0; i != dir.n_files; ++i) {
      tinydir_file file;
      tinydir_readfile_n(&dir, &file, i);

      if (!file.is_dir) {
        files.emplace_back(file.name);
      }
    }

    tinydir_close(&dir);

    std::sort(files.begin(), files.end(), [](const std::string & l, const std::string & r) {
      return std::stoi(l.substr(0, l.find("_"))) < std::stoi(r.substr(0, r.find("_")));
    });

    for (const auto & file : files) {
      paths.emplace_back("ff_save/states/" + file);
    }

    if (paths.size() != 0) {
      resume_from_saved_state(paths);
    } else {
      save_states = true;
      WARNING_MSG("Directory './ff_save' does not exist or has no content");
      logger << "Directory './ff_save' does not exist or has no content\n";
      INFO_MSG("Starting new reconstruction and saving states");
      logger << "Starting new reconstruction and saving states\n";
      return;
    }

    // Create directories to ensure their existence
    // states and probes should exist at this point, otherwise loading should have failed
    mkdir("ff_save/states", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir("ff_save/tmp", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir("ff_save/probes", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::resume_from_saved_state(const std::vector<std::string>& file_paths_) {
    if (verbosity > SILENT) {
      INFO_MSG("Loading saved states");
      logger << "Loading saved states\n";
    }

    load_anchor_points = false;

    std::ifstream v_file;
    std::string line;
    v_file.open("ff_save/var_order.gz");

    // Variable order
    if (v_file.is_open()) {
      igzstream order_file;
      order_file.open("ff_save/var_order.gz");

      while (std::getline(order_file, line)) {
        auto tmp_vec = parse_vector_32(line, 2);
        optimal_var_order[tmp_vec[0]] = tmp_vec[1];
      }

      order_file.close();

      std::string var_order = "Using optimized variable order: (";

      for (size_t i = 0; i != n; ++i) {
        if (!change_var_order && optimal_var_order[i] != i) {
          change_var_order = true;
        }

        if (i != n - 1) {
          var_order += "x" + std::to_string(optimal_var_order[i] + 1) + ", ";
        } else {
          var_order += "x" + std::to_string(optimal_var_order[i] + 1) + ")";
        }
      }

      if (verbosity > SILENT) {
        INFO_MSG(var_order);
      }

      logger << var_order << "\n";
    } else {
      if (verbosity > SILENT) {
        INFO_MSG("Using default variable order");
      }

      logger << "Using default variable order\n";
    }

    v_file.close();

    // Factors
    parsed_factors = parse_factors(n);
    size_t tmp_size_factors = parsed_factors.size();

    if (verbosity > SILENT) {
      INFO_MSG("Parsed " + std::to_string(tmp_size_factors) + " factor(s)");
    }

    logger << "Parsed " << tmp_size_factors << " factor(s)\n";

    factor_scan = false;

    factors_rf =  parse_factors_rf();

    igzstream validation_file;
    validation_file.open("ff_save/validation.gz");
    v_file.open("ff_save/validation.gz");

    if (v_file.is_open()) {
      if (!precomputed_probes) {
        std::getline(validation_file, line);
        std::vector<FFInt> values = parse_vector_FFInt(line);

        std::vector<FFInt> result = bb.eval(values, 0);
        size_t counter = 0;

        for (const auto & el : parsed_factors) {
          result[el.first] /= el.second.evaluate_pre(values)[0];
        }

        while (std::getline(validation_file, line)) {
          if (std::stoul(line) != result[counter]) {
            ERROR_MSG("Validation failed: Entry " + std::to_string(counter) + " does not match the black-box result!");
            logger << "Validation failed: Entry " + std::to_string(counter) + " does not match the black-box result!\n";
            logger.close();
            std::exit(EXIT_FAILURE);
          }

          ++counter;
        }

        if (counter != result.size()) {
          ERROR_MSG("Validation failed: Number of entries does not match the black box!");
          logger << "Validation failed: Number of entries does not match the black box!\n";
          std::exit(EXIT_FAILURE);
        }
      } else {
        WARNING_MSG("Validating states is disabled with precomputed probes.");
        logger << "Validating states is disabled with precomputed probes.\n";
      }
    } else {
      ERROR_MSG("Validation file not found!");
      logger << "Validation file not found!\n";
      std::exit(EXIT_FAILURE);
    }

    v_file.close();
    validation_file.close();

    save_states = true;
    resume_from_state = true;
    file_paths = file_paths_;
    items = static_cast<uint32_t>(file_paths.size());
    prime_it = 400; // increase so that the minimum is the mininmum of the files

    for (uint32_t i = 0; i != items; ++i) {
      prime_it = std::min(prime_it, parse_prime_number(file_paths[i]));
    }

    FFInt::set_new_prime(primes()[prime_it]);
    bb.prime_changed_internal();

    for (auto& el : parsed_factors) {
      el.second.precompute_tokens();
    }

    tmp_rec.start_from_saved_file(file_paths[0]);

    // Get probe files
    tinydir_dir dir;
    tinydir_open_sorted(&dir, "ff_save/probes");

    std::vector<std::string> files;
    std::vector<std::string> probe_files;

    for (size_t i = 0; i != dir.n_files; ++i) {
      tinydir_file file;
      tinydir_readfile_n(&dir, &file, i);

      if (!file.is_dir) {
        files.emplace_back(file.name);
      }
    }

    tinydir_close(&dir);

    std::sort(files.begin(), files.end(), [](const std::string & l, const std::string & r) {
      return std::stoi(l.substr(0, l.find("_"))) < std::stoi(r.substr(0, r.find("_")));
    });

    for (const auto & file : files) {
      probe_files.emplace_back("ff_save/probes/" + file);
    }

    if (probe_files.size() != items) {
      ERROR_MSG("Mismatch in number of probe files. Found " + std::to_string(probe_files.size()) + " files for " + std::to_string(items) + " functions");
      logger << "Mismatch in number of probe files. Found " << std::to_string(probe_files.size()) << " files for " << std::to_string(items) + " functions\n";
      std::exit(EXIT_FAILURE);
    }

    for (uint32_t i = 0; i != items; ++i) {
      RatReconst* rec = new RatReconst(n);
      /*std::pair<bool, uint32_t> shift_prime = */rec->start_from_saved_file(file_paths[i]);

      // Fill in already used ts from prior run
      auto tmp = rec->read_in_probes(probe_files[i]);

      for (const auto & already_chosen_t : tmp) {
        auto it = chosen_t.find(already_chosen_t.first);

        if (it != chosen_t.end()) {
          // Set the counter for the computed probes to the minimum of all RatReconst, i.e. the common ground
          if (static_cast<uint32_t>(already_chosen_t.second.size()) < started_probes[already_chosen_t.first]) {
            started_probes[already_chosen_t.first] = static_cast<uint32_t>(already_chosen_t.second.size());
          }

          // chosen_t shall contain the union of all t
          for (const auto & t : already_chosen_t.second) {
            if (it->second.find(t) == it->second.end()) {
              it->second.emplace(t);
            }
          }
        } else {
          chosen_t.emplace(std::make_pair(already_chosen_t.first, already_chosen_t.second));
          started_probes.emplace(std::make_pair(already_chosen_t.first, static_cast<uint32_t>(already_chosen_t.second.size())));
        }
      }

      if (safe_mode)
        rec->set_safe_interpolation();

      rec->set_tag(std::to_string(i));

      if (rec->is_done()) {
        ++items_done;

        reconst.emplace_back(std::make_tuple(i, DONE, rec));
      } else {
        if (rec->is_new_prime()) {
          probes_for_next_prime = std::max(probes_for_next_prime, rec->get_num_eqn());
          ++items_new_prime;

          if (rec->get_prime() == prime_it + 1) {
            load_anchor_points = true;
          }
        }

        reconst.emplace_back(std::make_tuple(i, RECONSTRUCTING, rec));
      }

      if (verbosity > SILENT) {
        std::cerr << "\033[1;34mFireFly info:\033[0m " << i << " / " << items << "\r";
      }
    }

    auto it = started_probes.find(std::vector<uint32_t> (n - 1, 1));

    if (it == started_probes.end()) {
      started_probes.emplace(std::vector<uint32_t> (n - 1, 1), 0);
    }

    if (prime_it == 0 && items != items_new_prime + items_done) {
      load_anchor_points = false;
      std::ifstream anchor_point_file;
      anchor_point_file.open("ff_save/anchor_points");

      if (anchor_point_file.is_open()) {
        std::getline(anchor_point_file, line);
        tmp_rec.set_anchor_points(parse_vector_FFInt(line, static_cast<int>(n)));

        const auto tmp_an_vec = tmp_rec.get_anchor_points();

        for (auto & rec : reconst) {
          std::get<2>(rec)->set_anchor_points(tmp_an_vec);
        }
      } else {
        ERROR_MSG("Anchor point file not found!");
        logger << "Anchor point file not found!\n";
        logger.close();
        std::exit(EXIT_FAILURE);
      }

      anchor_point_file.close();

      std::ifstream shift_file;
      shift_file.open("ff_save/shift");

      if (shift_file.is_open()) {
        std::getline(shift_file, line);
        tmp_rec.set_shift(parse_vector_FFInt(line, static_cast<int>(n)));
        shift = tmp_rec.get_zi_shift_vec();
      } else {
        ERROR_MSG("Shift file not found!");
        logger << "Shift file not found!\n";
        logger.close();
        std::exit(EXIT_FAILURE);
      }
    }

    if (safe_mode) {
      scan = false;
    }

    if (scan) {
      if (prime_it == 0 && items_new_prime != items) {
        std::ifstream file;
        file.open("ff_save/scan");

        if (file.is_open()) {
          scan = false;
        } else {
          ERROR_MSG("Cannot resume from saved state because the scan was not completed.");
          ERROR_MSG("Please remove the directory 'ff_save' and start from the beginning.");
          logger << "Cannot resume from saved state because the scan was not completed.\n";
          logger << "Please remove the directory 'ff_save' and start from the beginning.\n";
          logger.close();
          std::exit(EXIT_FAILURE);
        }

        file.close();
      } else {
        scan = false;
      }
    }

    logger << "All files loaded | Done: " << std::to_string(items_done) << " / " << std::to_string(items) <<
      " | " << "Requires new prime field: " << std::to_string(items_new_prime) << " / " << std::to_string(items - items_done) << "\n";

    if (verbosity > SILENT) {
      INFO_MSG("All files loaded | Done: " + std::to_string(items_done) + " / " + std::to_string(items) +
               " | " + "Requires new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items - items_done));
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::set_anchor_points(const std::vector<std::vector<std::uint64_t>>& anchor_points_) {
    anchor_points = anchor_points_;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::set_shifts(const std::vector<std::vector<std::uint64_t>>& shifts_) {
    shifts = shifts_;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::load_precomputed_probes() {
#ifndef WITH_MPI
    precomputed_probes = true;
#else
    WARNING_MSG("Precomputed probes are not supported with MPI enabled!");
#endif
  }

  template<typename BlackBoxTemp>
  bool Reconstructor<BlackBoxTemp>::reconstruction_done() {
    std::lock_guard<std::mutex> lock(status_control);
    return done;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::set_safe_interpolation() {
    safe_mode = true;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::reconstruct(uint32_t prime_counter) {
    start = std::chrono::high_resolution_clock::now();

    if (!aborted || resumed)
      done = false;

    if (resumed)
      resumed = false;

    if (aborted)
      aborted = false;

#ifdef WITH_MPI
    if (precomputed_probes) {
      ERROR_MSG("Precomputed probes are not supported with MPI support enabled!");
      return;
    }

    ThreadPool tp_comm(1);
    tp_comm.run_priority_task([this](uint32_t thread_id) {
      {
        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        cond_val.wait(lock_probe_queue, [this](){return continue_communication;});

        continue_communication = false;
        proceed = true;

        cond_val.notify_one();
      }

      mpi_communicate();
    });
    bool empty_bb = false;
#endif

    if (!resume_from_state) {
      logger << "\n" << "Promote to new prime field: F(" << std::to_string(primes()[prime_it]) << ")\n";

      if (verbosity > SILENT) {
        std::cerr << "\n";
        INFO_MSG("Promote to new prime field: F(" + std::to_string(primes()[prime_it]) + ")");
      }

      if (safe_mode) {
        tmp_rec.set_safe_interpolation();

        if (factor_scan) {
          WARNING_MSG("Disabled factor scan in safe mode!");
          logger << "Disabled factor scan in safe mode!\n";
          factor_scan = false;
        }

        if (scan) {
          WARNING_MSG("Disabled shift scan in safe mode!");
          logger << "Disabled shift scan in safe mode!\n";
          scan = false;
        }
      }

      if (factor_scan) {
        RatReconst::reset();
        scan_for_factors();

        tmp_rec = RatReconst(n);
        scanned_factors = true;
        if (items == 0) {
          scan = false;
          done = true;
#ifdef WITH_MPI
          empty_bb = true;
#endif
        }

	if (stop_after_factors) {
	  return;
	}
      }

      if (scan) {
        scan_for_shift();

        if (items == 0) {
          scan = false;
          done = true;
#ifdef WITH_MPI
          empty_bb = true;
#endif
        } else {
          queue_new_ones();
        }
      } else {
        if (scanned_factors && items == 0) {
          done = true;
#ifdef WITH_MPI
          empty_bb = true;
#endif
        } else {
          start_first_runs(!scanned_factors);
        }

        if (items == 0) {
          done = true;
#ifdef WITH_MPI
          empty_bb = true;
#endif
        }
      }
    } else {
      scan = false;

      if (items_done == items) {
        done = true;
      }
    }

    if (!done) {
      if (save_states && !load_anchor_points) {
        std::string tmp_str = "";
        std::ofstream file;
        file.open("ff_save/shift");
        tmp_str = "";

        for (const auto & el : tmp_rec.get_zi_shift_vec()) {
          tmp_str += std::to_string(el.n) + " ";
        }

        tmp_str += std::string("\n");
        file << tmp_str;
        file.close();
      }

      run_until_done(prime_counter);
    }
#ifdef WITH_MPI
    else {
      if (!empty_bb) {
        mpi_setup();

        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        continue_communication = true;

        cond_val.notify_one();

        cond_val.wait(lock_probe_queue, [this](){return proceed;});

        proceed = false;
      }

      {
        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        requested_probes = std::deque<std::pair<uint64_t, std::vector<FFInt>>>();
        new_jobs = false;

        cond_val.wait(lock_probe_queue, [this](){return proceed;});

        proceed = false;
      }

      tp.kill_all();

      new_jobs = false; // to make sure that no new jobs have been started

      {
        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        continue_communication = true;

        cond_val.notify_one();

        cond_val.wait(lock_probe_queue, [this](){return proceed;});

        proceed = false;
      }
    }
#endif

    if (one_done || one_new_prime) {
      logger << "Probe: " << std::to_string(probes_fed) <<
                 " | Done: " << std::to_string(items_done) << " / " << std::to_string(items) <<
                 " | " << "Requires new prime field: " << std::to_string(items_new_prime) << " / " << std::to_string(items - items_done) << "\n";
    }

    if (done) {
      logger << "Completed reconstruction in "
        << std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count())
        << " s | " << std::to_string(total_iterations) << " probes in total\n"
        << "Required prime fields: " << std::to_string(prime_it) << " + 1\n"
        << "Average time of the black-box probe: " << std::to_string(average_black_box_time) << " s\n";
      logger.close();
    }

    if (verbosity > SILENT) {
      if (one_done || one_new_prime) {
        INFO_MSG("Probe: " + std::to_string(probes_fed) +
                 " | Done: " + std::to_string(items_done) + " / " + std::to_string(items) +
                 " | " + "Requires new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items - items_done) + "\n");
      }

      if (done) {
        INFO_MSG("Completed reconstruction in " +
                 std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count()) + " s | " + std::to_string(total_iterations) + " probes in total");
        INFO_MSG("Required prime fields: " + std::to_string(prime_it) + " + 1");
        INFO_MSG("Average time of the black-box probe: " + std::to_string(average_black_box_time) + " s");
      }
    }
  }

  template<typename BlackBoxTemp>
  std::vector<std::string> Reconstructor<BlackBoxTemp>::get_factors_string(const std::vector<std::string>& vars) {
    std::vector<std::string> result {};

    for (size_t i = 0; i != items; ++ i) {
      std::unordered_map<std::vector<uint32_t>, RationalNumber, UintHasher> rnmap {};
      rnmap.emplace(std::vector<uint32_t> (n, 0), RationalNumber (1,1));
      Polynomial dummy_pol (rnmap);

      RationalFunction rf (dummy_pol, dummy_pol);

      if (change_var_order) {
	rf.set_var_order(optimal_var_order);
      }

      if (factors_rf.find(i) != factors_rf.end()) {
	for (const auto& factor : factors_rf[i]) {
	  rf.add_factor(factor);
	}
      }

      result.emplace_back(rf.to_string(vars));
    }

    return result;
  }

  template<typename BlackBoxTemp>
  std::vector<RationalFunction> Reconstructor<BlackBoxTemp>::get_result() {
    std::vector<RationalFunction> result {};
    uint32_t counter = 0;

    for (auto & rec : reconst) {
      if (std::get<1>(rec) == DONE) {
        result.emplace_back(std::get<2>(rec)->get_result());

        if (change_var_order) {
          result.back().set_var_order(optimal_var_order);
        }

        if (factors_rf.find(counter) != factors_rf.end()) {
          for (const auto& factor : factors_rf[counter]) {
            result.back().add_factor(factor);
          }
        }

        ++counter;
      }
    }

    return result;
  }

  template<typename BlackBoxTemp>
  std::vector<RationalFunctionFF> Reconstructor<BlackBoxTemp>::get_result_ff() {
    std::vector<RationalFunctionFF> result {};

    for (auto & rec : reconst) {
      result.emplace_back(std::get<2>(rec)->get_result_ff());
    }

    return result;
  }

  template<typename BlackBoxTemp>
  std::vector<std::pair<std::string, RationalFunction>> Reconstructor<BlackBoxTemp>::get_early_results() {
    if (factor_scan || scan) {
      return std::vector<std::pair<std::string, RationalFunction>> {};
    }

    std::lock_guard<std::mutex> lock_clean(clean);

    std::vector<std::pair<std::string, RationalFunction>> result;

    for (auto & rec : reconst) {
      if (std::get<1>(rec) == DONE) {
        if (save_states) {
          result.emplace_back(std::make_pair(std::get<2>(rec)->get_tag_name(), std::get<2>(rec)->get_result()));
        } else {
          result.emplace_back(std::make_pair(std::to_string(std::get<0>(rec)), std::get<2>(rec)->get_result()));
        }

        if (change_var_order) {
          result.back().second.set_var_order(optimal_var_order);
        }

        if (factors_rf.find(std::get<0>(rec)) != factors_rf.end()) {
          for (const auto& factor : factors_rf[std::get<0>(rec)]) {
            result.back().second.add_factor(factor);
          }
        }

        // TODO delete or clear memory?
        std::get<1>(rec) = DELETE;
      }
    }

    return result;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::scan_for_shift() {
    logger << "Scanning for a sparse shift\n";

    if (verbosity > SILENT)
      INFO_MSG("Scanning for a sparse shift");

    std::vector<uint32_t> current_shift(n, 0);
    bool first = true;
    bool found_shift = false;
    uint32_t counter = 0;

    tmp_rec.scan_for_sparsest_shift();

    start_first_runs(!scanned_factors);

    if (items == 0) {
      return;
    }

    uint32_t max_deg_num = 0;
    uint32_t max_deg_den = 0;

    // Run this loop until a proper shift is found
    while (!found_shift) {
      if (!first) {
        if (counter != 0) {
          auto shift_pair = generate_next_permutation(current_shift);
          if (shift_pair.first) {
            current_shift = shift_pair.second;
            tmp_rec.set_zi_shift(current_shift);
            shift = tmp_rec.get_zi_shift_vec();
            queue_new_ones();
          } else {
            break;
          }
        } else {
          tmp_rec.set_zi_shift(current_shift);
          shift = tmp_rec.get_zi_shift_vec();
          queue_new_ones();
        }
      }

      run_until_done();

      found_shift = true;

      for (auto & rec : reconst) {
        std::get<1>(rec) = RECONSTRUCTING;

        if (!(std::get<2>(rec)->is_shift_working())) {
          found_shift = false;
        }

        if (first) {
          std::pair<uint32_t, uint32_t> degs = std::get<2>(rec)->get_max_deg();

          max_deg_num = std::max(max_deg_num, degs.first);
          max_deg_den = std::max(max_deg_den, degs.second);
        }
      }

      if (first) {
        found_shift = false;
        first = false;

        logger << "Maximum degree of numerator: " << std::to_string(max_deg_num)
          << " | Maximum degree of denominator: " << std::to_string(max_deg_den) << "\n";
        logger.close();
        logger.open("firefly.log", std::ios_base::app);

        if (verbosity > SILENT) {
          INFO_MSG("Maximum degree of numerator: " + std::to_string(max_deg_num) + " | Maximum degree of denominator: " + std::to_string(max_deg_den));
        }
      } else {
        ++counter;
      }

      items_done = 0;
      done = false;
    }

    if (found_shift) {
      tmp_rec.set_zi_shift(current_shift);
    } else {
      tmp_rec.set_zi_shift(std::vector<uint32_t> (n, 1));
    }

    shift = tmp_rec.get_zi_shift_vec();

    for (auto & rec : reconst) {
      std::get<1>(rec) = RECONSTRUCTING;
      std::get<2>(rec)->accept_shift();
    }

    scan = false;

    if (save_states == true) {
      std::ofstream file;
      file.open("ff_save/scan");
      file.close();
    }

    if (found_shift) {
      std::string msg = "";

      for (const auto & el : current_shift) {
        msg += std::to_string(el) + ", ";
      }

      msg = msg.substr(0, msg.size() - 2);

      logger << "Found a sparse shift after " << std::to_string(counter + 1) << " scans\n"
        << "Shift the variable tuple (" << msg << ")\n";

      if (verbosity > SILENT) {
        INFO_MSG("Found a sparse shift after " + std::to_string(counter + 1) + " scans");
        INFO_MSG("Shift the variable tuple (" + msg + ")");
      }
    } else {
      logger << "Found no sparse shift after " << std::to_string(counter + 1) << " scans\n";

      if (verbosity > SILENT)
        INFO_MSG("Found no sparse shift after " + std::to_string(counter + 1) + " scans");
    }

    logger << "Completed scan in "
      << std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - prime_start).count())
      << " s | " << std::to_string(total_iterations) << " probes in total\n";

    logger << "Average time of the black-box probe: " << std::to_string(average_black_box_time) << " s\n\n";
    logger << "Proceeding with interpolation over prime field F(" << std::to_string(primes()[prime_it]) << ")\n";
    logger.close();
    logger.open("firefly.log", std::ios_base::app);

    if (verbosity > SILENT) {
      INFO_MSG("Completed scan in " + std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - prime_start).count()) +
               " s | " + std::to_string(total_iterations) + " probes in total");
      INFO_MSG("Average time of the black-box probe: " + std::to_string(average_black_box_time) + " s\n");
      INFO_MSG("Proceeding with interpolation over prime field F(" + std::to_string(primes()[prime_it]) + ")");
    }

    prime_start = std::chrono::high_resolution_clock::now();
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::scan_for_factors() {
    auto clock_1 = std::chrono::high_resolution_clock::now();
    factorizations.reserve(n);

    if (save_states) {
      mkdir("ff_save", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir("ff_save/factors", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir("ff_save/factors_rf", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    if (verbosity > SILENT) {
      INFO_MSG("Scanning for factors");
    }

    logger << "Scanning for factors\n";

    bool first = true;
    uint32_t total_number_of_factors = 0;
    uint32_t number_of_factors = 0;
    max_degs = std::vector<uint32_t> (n, 0);
    shift = std::vector<FFInt> (n, 0);
    rand_zi_fac = std::vector<FFInt> (n, 0);
    std::unordered_map<uint32_t, std::pair<std::unordered_set<std::string>, std::unordered_set<std::string>>> factors_str {};
    std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> max_deg_map_complete_tmp {};

    // Scan in all variables
    for (size_t i = 0; i != n; ++i) {
      std::unordered_map<uint32_t, std::unordered_map<uint32_t, fmpzxx>> combined_ni {};
      std::unordered_map<uint32_t, std::unordered_map<uint32_t, fmpzxx>> combined_di {};
      std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> max_deg_map {};
      fmpzxx combined_prime = fmpzxx(FFInt::p);
      bool fac_done = false;
      curr_var = "x" + std::to_string(i + 1);
      possible_factors_bb_counter.clear();
      prime_it_fac = 0;
      factors_degs.clear();

      while (!fac_done) {
        std::unordered_map<uint32_t,std::pair<std::unordered_set<std::string>, std::unordered_set<std::string>>> possible_factors {};
        fmpzxx tmp_combined_prime = combined_prime;

        for (int scan_n = 0; scan_n != 2; ++scan_n) {
          number_of_factors = 0;

          for (size_t j = 0; j != n; ++j) {
            if (j == i) {
              rand_zi_fac[j] = 1;
            } else {
              rand_zi_fac[j] = tmp_rec.get_rand_64();
            }
          }

          start_first_runs(first);
          first = false;

          if (items == 0) {
            scan = false;
            scanned_factors = true;
            return;
          }

          if (prime_it_fac == 0 && scan_n == 0) {
            if (verbosity > SILENT) {
              INFO_MSG("Scanning for factors in " + curr_var);
            }

            logger << "Scanning for factors in " << curr_var << "\n";
          }

          run_until_done();
          prime_it = 0;

          uint32_t counter = 0;

          // Get factors
          for (auto& rec : reconst) {
            if (scan_n == 0) {
              if (possible_factors_bb_counter.find(counter) != possible_factors_bb_counter.end()) {
                auto tmp_fac_num_den = std::get<2>(rec)->get_factors_ff();
                std::unordered_set<std::string> tmp_fac_num {};
                std::unordered_set<std::string> tmp_fac_den {};
                tmp_fac_num.insert(tmp_fac_num_den.first.begin(), tmp_fac_num_den.first.end());
                tmp_fac_den.insert(tmp_fac_num_den.second.begin(), tmp_fac_num_den.second.end());
                possible_factors.emplace(counter, std::make_pair(tmp_fac_num, tmp_fac_den));
              }

              if (prime_it_fac == 0) {
                auto tmp_fac_num_den = std::get<2>(rec)->get_factors_ff();
                std::unordered_set<std::string> tmp_fac_num {};
                std::unordered_set<std::string> tmp_fac_den {};
                tmp_fac_num.insert(tmp_fac_num_den.first.begin(), tmp_fac_num_den.first.end());
                tmp_fac_den.insert(tmp_fac_num_den.second.begin(), tmp_fac_num_den.second.end());
                possible_factors.emplace(counter, std::make_pair(tmp_fac_num, tmp_fac_den));

                // Get maximum degrees
                auto tmp_max_degs = std::get<2>(rec)->get_max_deg();
                max_deg_map.emplace(std::make_pair(counter, tmp_max_degs));
                max_deg_map_complete_tmp[counter].emplace_back(tmp_max_degs);

                uint32_t max_val = std::max(tmp_max_degs.first, tmp_max_degs.second);
                if (max_val > max_degs[i]) {
                  max_degs[i] = max_val;
                }

                uint32_t tmp_n_fac = tmp_fac_num_den.first.size() + tmp_fac_num_den.second.size();
                number_of_factors += tmp_n_fac;

                // Save occurring degrees to perform Thiele just once
                auto tmp_rf = std::get<2>(rec)->get_result_ff();
                std::list<uint32_t> deg_list_num {};
                std::list<uint32_t> deg_list_den {};

                for (const auto& el : tmp_rf.numerator.coefs) {
                  deg_list_num.emplace_back(el.first[0]);
                }

                for (const auto& el : tmp_rf.denominator.coefs) {
                  deg_list_den.emplace_back(el.first[0]);
                }

                factors_degs.emplace(std::make_pair(counter, std::make_pair(deg_list_num, deg_list_den)));

                if (tmp_n_fac != 0) {
                  possible_factors_bb_counter.emplace(counter);
                }
              }
            } else if (scan_n == 1 && possible_factors_bb_counter.find(counter) != possible_factors_bb_counter.end()) {
              // Rewrite and store in result objects. Compare to previous factorizations
              std::unordered_set<uint32_t> fac_nums_c {};
              std::unordered_set<uint32_t> fac_dens_c {};

              auto tmp_factors = std::get<2>(rec)->get_factors_ff();

              // Numerator
              uint32_t fac_counter = 0;
              for (const auto& tmp_factor : tmp_factors.first) {
                if (possible_factors[counter].first.find(tmp_factor) != possible_factors[counter].first.end()) {
                  fac_nums_c.emplace(fac_counter);
                  ++number_of_factors;
                }

                ++fac_counter;
              }

              fac_counter = 0;
              // Denominator
              for (const auto& tmp_factor : tmp_factors.second) {
                if (possible_factors[counter].second.find(tmp_factor) != possible_factors[counter].second.end()) {
                  fac_dens_c.emplace(fac_counter);
                  ++number_of_factors;
                }

                ++fac_counter;
              }

              if (fac_nums_c.empty() && fac_dens_c.empty()) {
                possible_factors_bb_counter.erase(counter);
              } else {
                auto canonical_factors = std::get<2>(rec)->get_canonical_factors(std::make_pair(fac_nums_c, fac_dens_c));

                if (prime_it_fac == 0) {
                  // Init first combinations
                  std::unordered_map<uint32_t, fmpzxx> tmp_combined_ni {};
                  std::unordered_map<uint32_t, fmpzxx> tmp_combined_di {};
                  uint32_t fac_max_deg_num = 0, fac_max_deg_den = 0;

                  for (const auto& mon : canonical_factors.first) {
                    tmp_combined_ni.emplace(std::make_pair(mon.first, fmpzxx(mon.second)));
                    fac_max_deg_num = std::max(fac_max_deg_num, mon.first);
                  }

                  for (const auto& mon : canonical_factors.second) {
                    tmp_combined_di.emplace(std::make_pair(mon.first, fmpzxx(mon.second)));
                    fac_max_deg_den = std::max(fac_max_deg_den, mon.first);
                  }

                  max_deg_map[counter].first -= fac_max_deg_num;
                  max_deg_map[counter].second -= fac_max_deg_den;
                  max_deg_map_complete_tmp[counter].back().first = max_deg_map[counter].first;
                  max_deg_map_complete_tmp[counter].back().second = max_deg_map[counter].second;

                  combined_ni[counter] = tmp_combined_ni;
                  combined_di[counter] = tmp_combined_di;
                } else {
                  // First, check if done, else combine results
                  bool run_test = true;
                  bool combine_results = false;
                  std::unordered_map<std::vector<uint32_t>, RationalNumber, UintHasher> tmp_gni {};
                  std::unordered_map<std::vector<uint32_t>, RationalNumber, UintHasher> tmp_gdi {};

                  // Reconstruct numerator
                  for (const auto& ci : combined_ni[counter]) {
                    fmpzxx a = ci.second;
                    auto res = get_rational_coef(a, combined_prime);

                    if (res.first) {
                      tmp_gni.emplace(std::make_pair(std::vector<uint32_t> (1, ci.first), res.second));
                    } else {
                      run_test = false;
                      break;
                    }
                  }

                  // Reconstruct denominator
                  if (run_test) {
                    for (const auto& ci : combined_di[counter]) {
                      fmpzxx a = ci.second;
                      auto res = get_rational_coef(a, combined_prime);

                      if (res.first) {
                        tmp_gdi.emplace(std::make_pair(std::vector<uint32_t> (1, ci.first), res.second));
                      } else {
                        run_test = false;
                        combine_results = true;
                        break;
                      }
                    }
                  } else {
                    combine_results = true;
                  }

                  if (run_test) {
                    FFInt tmp_rand = tmp_rec.get_rand_64();

                    if (!tmp_gni.empty()) {
                      ff_map gi_ffi;
                      FFInt num = 0;

                      for (const auto& g_i : tmp_gni) {
                        FFInt n(g_i.second.numerator);
                        FFInt d(g_i.second.denominator);
                        gi_ffi.emplace(std::make_pair(g_i.first, n / d));
                      }

                      for (const auto& coeff : canonical_factors.first) {
                        num += coeff.second*tmp_rand.pow(coeff.first);
                      }

                      combine_results = !(PolynomialFF(1, gi_ffi).calc({tmp_rand}) == num);
                    }

                    if (!combine_results && !tmp_gdi.empty()) {
                      ff_map gi_ffi;
                      FFInt num = 0;

                      for (const auto& g_i : tmp_gdi) {
                        FFInt n(g_i.second.numerator);
                        FFInt d(g_i.second.denominator);
                        gi_ffi.emplace(std::make_pair(g_i.first, n / d));
                      }

                      for (const auto& coeff : canonical_factors.second) {
                        num += coeff.second*tmp_rand.pow(coeff.first);
                      }

                      combine_results = !(PolynomialFF(1, gi_ffi).calc({tmp_rand}) == num);
                    }
                  }

                  // combine results
                  if (combine_results) {
                    // numerator
                    if (!canonical_factors.first.empty()) {
                      combined_prime = combine_primes(canonical_factors.first, combined_ni[counter], tmp_combined_prime);
                    }

                    // denominator
                    if (!canonical_factors.second.empty()) {
                      combined_prime = combine_primes(canonical_factors.second, combined_di[counter], tmp_combined_prime);
                    }
                  } else {
                    possible_factors_bb_counter.erase(counter);
                    combined_ni.erase(counter);
                    combined_di.erase(counter);
                    Polynomial tmp_numerator;
                    Polynomial tmp_denominator;

                    // Rewrite to result
                    if (!tmp_gni.empty()) {
                      tmp_numerator = Polynomial(tmp_gni);
                      factors_str[counter].first.emplace(tmp_numerator.to_string({curr_var}));
                      tmp_numerator.set_var_pos(i);
                    }

                    if (!tmp_gdi.empty()) {
                      tmp_denominator = Polynomial(tmp_gdi);
                      factors_str[counter].second.emplace(tmp_denominator.to_string({curr_var}));
                      tmp_denominator.set_var_pos(i);
                    }

                    factors_rf[counter].emplace_back(RationalFunction(tmp_numerator, tmp_denominator));
                  }
                }
              }
            }

            ++counter;

            std::get<1>(rec) = DELETE;
          }

          uint32_t old_max_deg = max_degs[i];
          if (prime_it_fac == 0 && scan_n == 1) {
            max_degs[i] = 0;

            for (const auto& el : max_deg_map) {
              uint32_t max_val = std::max(el.second.first, el.second.second);

              if (max_val > max_degs[i]) {
                max_degs[i] = max_val;
              }
            }
          }

          if (verbosity > SILENT) {
            if (prime_it_fac == 0 && scan_n == 0) {
              INFO_MSG("Maximum degree of x" + std::to_string(i + 1)
                + ": " + std::to_string(max_degs[i]));
              if (max_degs[i] != 0) {
                INFO_MSG("Potential factors in x"
                  + std::to_string(i + 1) + ": " + std::to_string(number_of_factors));
              } else {
                INFO_MSG("No factors in x"
                  + std::to_string(i + 1) + ": " + std::to_string(number_of_factors));
              }
            } else if (prime_it_fac == 0 && scan_n == 1) {
              INFO_MSG("Identified factors in x" + std::to_string(i + 1) + ": "
                + std::to_string(number_of_factors));
              if (max_degs[i] != old_max_deg) {
                INFO_MSG("Maximum degree of x" + std::to_string(i + 1) + " after factoring: "
                  + std::to_string(max_degs[i]));
              }
             INFO_MSG("Starting reconstruction of coefficients");
            }
          }

          if (prime_it_fac == 0 && scan_n == 0) {
            logger << "Maximum degree of x" << std::to_string(i + 1) << ": "
              << std::to_string(max_degs[i]) << "\n";

            if (max_degs[i] != 0) {
              logger << "Potential factors in x"
                << std::to_string(i + 1) << ": " << std::to_string(number_of_factors)
                << "\n";
            } else {
              logger << "No factors in x"
                << std::to_string(i + 1) << ": " << std::to_string(number_of_factors)
                << "\n";
            }
          } else if (prime_it_fac == 0 && scan_n == 1) {
            total_number_of_factors += number_of_factors;
            logger << "Identified factors in x" << std::to_string(i + 1) << ": "
              << std::to_string(number_of_factors) << "\n";

            if (max_degs[i] != old_max_deg) {
              logger << "Maximum degree of x" << std::to_string(i + 1) << " after factoring: "
                << std::to_string(max_degs[i]) << "\n";
            }

            logger << "Starting reconstruction of coefficients\n";
            logger.close();
            logger.open("firefly.log", std::ios_base::app);
          }

          RatReconst::reset(false);
          clean_reconst();
          reconst.clear();

          reset_new_prime();
          items_done = 0;
          done = false;

          if (possible_factors_bb_counter.empty()) {
            fac_done = true;
            break;
          }
        }

        // Promote to new prime field
        if (!fac_done) {
          ++prime_it_fac;
          FFInt::set_new_prime(primes()[prime_it_fac]);
          bb.prime_changed_internal();
        } else {
          if (verbosity > SILENT) {
            INFO_MSG("Completed factor scan in " + curr_var + " | "
              + std::to_string(total_iterations) + " probes in total");
            INFO_MSG("Required prime fields: " + std::to_string(prime_it_fac + 1));
            INFO_MSG("Average time of the black-box probe: " + std::to_string(average_black_box_time) + " s\n");
          }

          logger << "Completed factor scan in " << curr_var << " | "
            << total_iterations << " probes in total\n";
          logger << "Required prime fields: " << prime_it_fac + 1 << "\n"
            << "Average time of the black-box probe: " << std::to_string(average_black_box_time) << " s\n\n";
        }
      }

      // Reset prime
      prime_it_fac = 0;
      FFInt::set_new_prime(primes()[prime_it_fac]);
      bb.prime_changed_internal();
    }

    std::vector<std::string> vars (n);

    // Reorder variables with regards to their maximum degree
    std::vector<uint32_t> indices (n);
    std::iota(indices.begin(), indices.end(), 0);

    std::stable_sort(indices.begin(), indices.end(),
                     [&](uint32_t i1, uint32_t i2) {return max_degs[i1] > max_degs[i2];});

    std::sort(max_degs.begin(), max_degs.end(), std::greater<uint32_t>());

    for (size_t i = 0; i != n; ++i) {
      optimal_var_order.emplace(std::make_pair(i, indices[i]));
      vars[i] = "x" + std::to_string(i + 1);
    }

    for (const auto& el : max_deg_map_complete_tmp) {
      max_deg_map_complete[el.first] = std::vector<std::pair<uint32_t, uint32_t>> (n);

      for (size_t i = 0; i != n; ++i) {
        max_deg_map_complete[el.first][i] = el.second[optimal_var_order[i]];
      }
    }

    std::string var_order = "Using optimized variable order: (";

    for (size_t i = 0; i != n; ++i) {
      if (!change_var_order && optimal_var_order[i] != i) {
        change_var_order = true;
      }

      if (i != n - 1) {
        var_order += "x" + std::to_string(optimal_var_order[i] + 1) + ", ";
      } else {
        var_order += "x" + std::to_string(optimal_var_order[i] + 1) + ")";
      }
    }

    if (save_states && change_var_order) {
      ogzstream file;
      std::string file_name = "ff_save/var_order.gz";
      file.open(file_name.c_str());

      for (size_t i = 0; i != n; ++i) {
        file << i << " " << optimal_var_order[i] << " \n";
      }

      file.close();
    }

    if (save_states) {
      for (const auto& tmp_rfs : factors_rf) {
        ogzstream file;
        std::string file_name = "ff_save/factors_rf/" + std::to_string(tmp_rfs.first) + ".gz";
        file.open(file_name.c_str());
        int tmp_var_pos = -1;

        for (const auto& tmp_rf : tmp_rfs.second) {
          tmp_var_pos = tmp_rf.numerator.get_var_pos();

          if (tmp_var_pos == -1)
            tmp_var_pos = tmp_rf.denominator.get_var_pos();

          file << "var\n";
          file << tmp_var_pos << "\n";
          file << "numerator\n";

          for (const auto& el : tmp_rf.numerator.coefs) {
            file << el.powers[0] << " ";
            file << el.coef.numerator.to_string() << " " << el.coef.denominator.to_string() << "\n";
          }

          file << "denominator\n";

          for (const auto& el : tmp_rf.denominator.coefs) {
            file << el.powers[0] << " ";
            file << el.coef.numerator.to_string() << " " << el.coef.denominator.to_string() << "\n";
          }
        }

        file.close();
      }
    }

#ifdef WITH_MPI
    MPI_Request* requests = new MPI_Request[world_size - 1];

    for (int i = 1; i != world_size; ++i) {
      MPI_Isend(&FFInt::p, 1, MPI_UINT64_T, i, FACTORS, MPI_COMM_WORLD, &requests[i - 1]);
    }

    MPI_Waitall(world_size - 1, requests, MPI_STATUSES_IGNORE);

    delete[] requests;
#endif

    for (const auto& tmp_fac : factors_str) {
      std::string tmp_fac_s = "";

      if (!tmp_fac.second.first.empty()) {
        tmp_fac_s += "(";
        for (const auto& tmp_fac_num : tmp_fac.second.first) {
          tmp_fac_s += "(" + tmp_fac_num + ")*";
        }

        tmp_fac_s.pop_back();
        tmp_fac_s += ")";
      }

      if (!tmp_fac.second.second.empty()) {
        if (!tmp_fac.second.first.empty()) {
          tmp_fac_s += "/(";
        } else {
          tmp_fac_s += "1/(";
        }

        for (const auto& tmp_fac_den : tmp_fac.second.second) {
          tmp_fac_s += "(" + tmp_fac_den + ")*";
        }

        tmp_fac_s.pop_back();
        tmp_fac_s += ")";
      }

      ShuntingYardParser parser = ShuntingYardParser();
      parser.parse_function(tmp_fac_s, vars);
      parser.precompute_tokens();
      parsed_factors.emplace(tmp_fac.first, parser);

#ifdef WITH_MPI
      const int batch_size = 2147483645; // 2147483647 is the largest signed 32-bit integer
      const int split = static_cast<int>(tmp_fac_s.size()) / batch_size;

      if (split > 1) {
        int tmp;
        if (static_cast<int>(tmp_fac_s.size()) % batch_size == 0) {
          tmp = -split;
        } else {
          tmp = - 1 - split;
        }
        MPI_Bcast(&tmp, 1, MPI_INT, master, MPI_COMM_WORLD);
      } else if (split == 1 && static_cast<int>(tmp_fac_s.size()) != batch_size) {
        int tmp = - 1 - split;
        MPI_Bcast(&tmp, 1, MPI_INT, master, MPI_COMM_WORLD);
      }

      for (int i = 0; i != 1 + split; ++i) {
        if (i == split) {
          int amount = static_cast<int>(tmp_fac_s.size()) - i * batch_size;
          if (amount != 0) {
            MPI_Bcast(&amount, 1, MPI_INT, master, MPI_COMM_WORLD);
            MPI_Bcast(&tmp_fac_s[i * batch_size], amount, MPI_CHAR, master, MPI_COMM_WORLD);
          }
        } else {
          int amount = batch_size;
          MPI_Bcast(&amount, 1, MPI_INT, master, MPI_COMM_WORLD);
          MPI_Bcast(&tmp_fac_s[i * batch_size], amount, MPI_CHAR, master, MPI_COMM_WORLD);
        }
      }

      uint32_t function_number = tmp_fac.first;
      MPI_Bcast(&function_number, 1, MPI_UINT32_T, master, MPI_COMM_WORLD);
#endif

      tmp_fac_s += ";\n";

      if (save_states) {
        ogzstream file;
        std::string file_name = "ff_save/factors/" + std::to_string(tmp_fac.first) + ".gz";
        file.open(file_name.c_str());
        file << tmp_fac_s;
        file.close();
      }
    }

#ifdef WITH_MPI
    int end = -1;
    MPI_Bcast(&end, 1, MPI_INT, master, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
#endif

    logger << "Completed factor scan in "
      << std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - clock_1).count())
      << " s | " << std::to_string(total_iterations) << " probes in total\n";

    logger << "Average time of the black-box probe: " << std::to_string(average_black_box_time) << " s\n";

    logger << "Found " << std::to_string(total_number_of_factors) << " factors in total\n";
    logger << var_order << "\n\n";

    if (!scan) {
      logger << "Proceeding with interpolation over prime field F(" << std::to_string(primes()[prime_it]) << ")\n";
    }

    logger.close();
    logger.open("firefly.log", std::ios_base::app);

    if (verbosity > SILENT) {
      INFO_MSG("Completed factor scan in " + std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - clock_1).count()) +
               " s | " + std::to_string(total_iterations) + " probes in total");
      INFO_MSG("Average time of the black-box probe: " + std::to_string(average_black_box_time) + " s");
      INFO_MSG("Found " + std::to_string(total_number_of_factors) + " factors in total");
      INFO_MSG(var_order + "\n");

      if (!scan) {
        INFO_MSG("Proceeding with interpolation over prime field F(" + std::to_string(primes()[prime_it]) + ")");
      }
    }

    factor_scan = false;
    factors_degs.clear();

    prime_start = std::chrono::high_resolution_clock::now();
  }

  template<typename BlackBoxTemp>
  fmpzxx Reconstructor<BlackBoxTemp>::combine_primes(const std::unordered_map<uint32_t, uint64_t>& poly,
                                                        std::unordered_map<uint32_t, fmpzxx>& combined_ci,
                                                        const fmpzxx& combined_prime) {
    std::pair<fmpzxx, fmpzxx> p3;

    for (const auto &it : poly) {
      p3 = run_chinese_remainder(combined_ci[it.first], combined_prime, it.second, FFInt::p, FFInt::p_inv);
      combined_ci[it.first] = p3.first;
    }

    return p3.second;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::start_first_runs(bool first) {
    prime_start = std::chrono::high_resolution_clock::now();

    if (anchor_points.size() > static_cast<std::size_t>(prime_it)) {
      std::vector<FFInt> tmp_an_vec;

      for (const auto value : anchor_points[prime_it]) {
        tmp_an_vec.emplace_back(value);
      }

      tmp_rec.set_anchor_points(tmp_an_vec);
    }

    if (shifts.size() > static_cast<std::size_t>(prime_it)) {
      shift.clear();
      shift.reserve(n);

      for (const auto value : shifts[prime_it]) {
        shift.emplace_back(value);
      }

      tmp_rec.set_shift(shift);
    }

    std::vector<uint32_t> zi_order;

    if (!factor_scan) {
      shift = tmp_rec.get_zi_shift_vec();
      zi_order = std::vector<uint32_t>(n - 1, 1);
    } else {
      zi_order = std::vector<uint32_t>(0, 1);
    }

#ifndef WITH_MPI
    (void) first; // void cast to silence unused parameter warning

    uint32_t to_start = thr_n;//* bunch_size;
    queue_probes(zi_order, to_start);
    started_probes.emplace(zi_order, to_start);

    if (precomputed_probes) {
      load_precomputed_probes_from_file();
    }
#else
    if (first) {
      send_first_jobs();

      uint32_t to_start = thr_n;//* bunch_size;
      queue_probes(zi_order, to_start);
      started_probes[zi_order] += to_start;
    } else {
      uint32_t to_start = buffer * worker_thread_count + thr_n;//* bunch_size;
      queue_probes(zi_order, to_start, true);
      started_probes.emplace(zi_order, to_start);

      cond_val.notify_one();

      {
        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        cond_val.wait(lock_probe_queue, [this](){return proceed;});

        proceed = false;
      }

      for (uint32_t j = 0; j != to_start; ++j) {
        tp.run_task([this](uint32_t thread_id) {
          get_job(thread_id);
        });
      }
    }
#endif

    std::vector<uint64_t> indices;
    std::vector<std::vector<FFInt>> probes;

    get_probe(indices, probes);

    std::vector<FFInt> t_vec;
    t_vec.reserve(indices.size());
    std::vector<std::vector<uint32_t>> zi_order_vec;
    zi_order_vec.reserve(indices.size());

    uint32_t count_ones = 0;

    {
      std::lock_guard<std::mutex> lock_probe_queue(mutex_probe_queue);

      for (const auto& index : indices) {
        auto tmp = std::move(index_map[index]);
        index_map.erase(index);
        t_vec.emplace_back(tmp.first);
        zi_order_vec.emplace_back(std::move(tmp.second));

        if ((prime_it == 0 || safe_mode == true) && (zi_order_vec.back() == std::vector<uint32_t>(n - 1, 1) || (factor_scan && zi_order_vec.back() == std::vector<uint32_t>(0, 1)))) {
          ++count_ones;
        }
      }
    }

    if (count_ones != 0) {
      std::lock_guard<std::mutex> lock_status(job_control);

      balance_of_ones += count_ones;
    }

    {
      std::lock_guard<std::mutex> lock(future_control);

      if (first_print) {
        logger << "Time for the first black-box probe: " << std::to_string(average_black_box_time) << " s\n";

        if (verbosity > SILENT) {
          INFO_MSG("Time for the first black-box probe: " + std::to_string(average_black_box_time) + " s");
        }
      }
    }

    items = static_cast<uint32_t>(probes.size());
    size_t tag_size = tags.size();

#ifdef WITH_MPI
    if (first) {
      std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

      continue_communication = true;

      cond_val.notify_one();

      cond_val.wait(lock_probe_queue, [this](){return proceed;});

      proceed = false;
    }
#endif

    if (tag_size != 0 && tag_size != items) {
      logger << "Number of tags does not match the black box!\n";
      ERROR_MSG("Number of tags does not match the black box!");
      logger.close();
      std::exit(EXIT_FAILURE);
    }

    ogzstream file;

    if (!factor_scan && save_states) {
      mkdir("ff_save", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir("ff_save/states", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir("ff_save/tmp", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir("ff_save/probes", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      std::ofstream anchor_file;
      anchor_file.open("ff_save/anchor_points");

      std::string tmp_str = "";

      for (const auto & el : tmp_rec.get_anchor_points()) {
        tmp_str += std::to_string(el.n) + " ";
      }

      tmp_str += std::string("\n");
      anchor_file << tmp_str;
      anchor_file.close();

      file.open("ff_save/validation.gz");

      std::vector<FFInt> rand_zi= tmp_rec.get_rand_zi_vec(zi_order);

      if (change_var_order) {
        std::vector<uint64_t> tmp_val (n);
        for (const auto & el : optimal_var_order) {
          if (el.first == 0) {
            tmp_val[el.second] = (t_vec[0] + shift[0]).n;
          } else {
            tmp_val[el.second] = (rand_zi[el.first - 1] * t_vec[0] + shift[el.first]).n;
          }
        }

        for (size_t i = 0; i != n; ++i) {
          file << tmp_val[i] << " ";
        }
      } else {
        file << (t_vec[0] + shift[0]).n << " ";

        for (uint32_t i = 1; i != n; ++i) {
          file << (rand_zi[i - 1] * t_vec[0] + shift[i]).n << " ";
        }
      }

      file << "\n";
    }

    for (uint32_t i = 0; i != items; ++i) {
      RatReconst* rec;
      if (!factor_scan) {
        rec = new RatReconst(n);

        if (safe_mode) {
          rec->set_safe_interpolation();
        }

        if (scan) {
          rec->scan_for_sparsest_shift();
        }

        if (!parsed_factors.empty()) {
          rec->set_individual_degree_bounds(max_deg_map_complete[i]);
        }

        if (save_states) {
          rec->set_tag(std::to_string(i));

          if (tag_size > 0) {
            rec->set_tag_name(tags[i]);
          } else {
            rec->set_tag_name(std::to_string(i));
          }

          file << ((probes)[i][0]).n << "\n";
        }
      } else {
        rec = new RatReconst(1);

        if (factors_degs.empty()) {
          rec->calc_factors(curr_var);
        } else {
          rec->calc_factors(curr_var, factors_degs[i]);
        }

        // Remove functions that are irreducible
        if (!possible_factors_bb_counter.empty() &&  possible_factors_bb_counter.find(i) == possible_factors_bb_counter.end()) {
          rec->set_prime_to_max();
        }
      }

      rec->feed(t_vec, probes[i], zi_order_vec, prime_it);
      std::tuple<bool, bool, uint32_t> interpolated_done_prime = rec->interpolate();

      int status = RECONSTRUCTING;

      if (std::get<2>(interpolated_done_prime) > prime_it) {
        ++items_new_prime;
      } else if (std::get<1>(interpolated_done_prime)) {
        status = DONE;
        ++items_done;
        one_done = true;
      }

      reconst.emplace_back(std::make_tuple(i, status, rec));
    }

    if (!factor_scan && save_states) {
      file.close();
      tags.clear();
    }

    if (verbosity > SILENT && first_print) {
      first_print = false;

      if (items == 0) {
        INFO_MSG("Black box has no entries");
        logger << "Black box has no entries\n";
        logger.close();
        logger.open("firefly.log", std::ios_base::app);
        return;
      }

      std::string msg = std::to_string(items) + " function(s) will be interpolated";

      if (factor_scan) {
        msg += "\n";
      }

      logger << msg << "\n";
      logger.close();
      logger.open("firefly.log", std::ios_base::app);

      INFO_MSG(msg);
    }

    queue_probes(zi_order, static_cast<uint32_t>(probes.front().size()));
    started_probes[zi_order] += static_cast<uint32_t>(probes.front().size());
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::queue_new_ones() {
#ifndef WITH_MPI
    uint32_t to_start = thr_n ;//* bunch_size; // TODO
    queue_probes(std::vector<uint32_t> (n - 1, 1), to_start);
#else
    uint32_t to_start = buffer * worker_thread_count + thr_n; // TODO: start even more? * bunch_size
    queue_probes(std::vector<uint32_t> (n - 1, 1), to_start, true);
#endif

    {
      std::lock_guard<std::mutex> lock(job_control);

      started_probes.emplace(std::vector<uint32_t> (n - 1, 1), to_start);
    }

#ifdef WITH_MPI
    cond_val.notify_one();

    {
      std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

      cond_val.wait(lock_probe_queue, [this](){return proceed;});

      proceed = false;
    }

    for (uint32_t j = 0; j != to_start; ++j) {
      tp.run_task([this](uint32_t thread_id) {
        get_job(thread_id);
      });
    }
#endif
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::run_until_done(uint32_t prime_counter) {
    std::vector<uint32_t> zi_order;
    if (!factor_scan) {
      zi_order = std::vector<uint32_t> (n - 1, 1);
    } else {
      zi_order = std::vector<uint32_t> (0, 1);
    }

    new_prime = false;

#ifdef WITH_MPI
    bool mpi_first_send = false;
#endif

    if (resume_from_state) {
      if (prime_it == 0 && items != items_new_prime + items_done) {
        logger << "Resuming in prime field: F(" << std::to_string(primes()[prime_it]) << ")\n";
        INFO_MSG("Resuming in prime field: F(" + std::to_string(primes()[prime_it]) + ")");

        {
          std::lock_guard<std::mutex> lock_status(feed_control);

          interpolate_jobs += items;
        }

        uint32_t counter = 0;

        for (auto & rec : reconst) {
          if (std::get<2>(rec)->get_prime() == 0) {
            ++counter;

            tp.run_priority_task([this, &rec](uint32_t thread_id) {
              (void)thread_id;
              interpolate_job(rec);
            });
          }
        }

        {
          std::lock_guard<std::mutex> lock_status(feed_control);

          interpolate_jobs -= (items - counter);
        }

#ifndef WITH_MPI
        // TODO don't start as much ones
        queue_probes(zi_order, thr_n);//* bunch_size);
        started_probes.at(zi_order) += thr_n;//* bunch_size);
#else
        // TODO optimize
        send_first_jobs();

        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          continue_communication = true;

          cond_val.notify_one();

          cond_val.wait(lock_probe_queue, [this](){return proceed;});

          proceed = false;
        }

        uint32_t to_start = thr_n;//* bunch_size; // TODO

        queue_probes(zi_order, to_start);
        started_probes[zi_order] += to_start;

#endif
      } else {
        new_prime = true;
#ifdef WITH_MPI
        proceed = true;
        mpi_first_send = true;
#endif
      }
    }

    while (!done) {
      if (new_prime) {
        if (factor_scan) {
          break;
        }

#ifdef WITH_MPI
        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          requested_probes = std::deque<std::pair<uint64_t, std::vector<FFInt>>>();
          new_jobs = false;

          cond_val.wait(lock_probe_queue, [this](){return proceed;});

          proceed = false;
        }
#endif

        tp.kill_all();

#ifdef WITH_MPI
        new_jobs = false; // to make sure that no new jobs have been started
#endif

        clean_reconst();

        if (!factor_scan && save_states) {
          for (uint32_t item = 0; item != items; ++item) {
            std::string probes_file_old = "ff_save/probes/" + std::to_string(item) + "_" + std::to_string(prime_it) + ".gz";
            std::string probes_file_new = "ff_save/probes/" + std::to_string(item) + "_" + std::to_string(prime_it + 1) + ".gz";
            std::rename(probes_file_old.c_str(), probes_file_new.c_str());

            if (prime_it) {
              std::string file_name_old = "ff_save/states/" + std::to_string(item) + "_" + std::to_string(prime_it - 1) + ".gz";
              std::string file_name_new = "ff_save/states/" + std::to_string(item) + "_" + std::to_string(prime_it) + ".gz";
              std::rename(file_name_old.c_str(), file_name_new.c_str());
            }
          }
        }

        ++prime_it;
        if(prime_it >= prime_counter || factor_scan) {
          done = true;
        }

#ifdef WITH_MPI
        if (!mpi_first_send) {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          continue_communication = true;
          cond_val.notify_one();

          cond_val.wait(lock_probe_queue, [this](){return proceed;});

          proceed = false;
        }
#endif

        total_iterations += iteration;

        {
          std::lock_guard<std::mutex> lock_print(print_control);

          if ((one_done || one_new_prime) && !factor_scan) {
            logger << "Probe: " << std::to_string(probes_fed)
            << " | Done: " << std::to_string(items_done)
            << " / " + std::to_string(items) << " | " << "Requires new prime field: "
            << std::to_string(items_new_prime) << " / " << std::to_string(items - items_done) << "\n";

            if (verbosity > SILENT) {
            INFO_MSG("Probe: " + std::to_string(probes_fed) +
                     " | Done: " + std::to_string(items_done) + " / " + std::to_string(items) +
                     " | " + "Requires new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items - items_done));
            }
          }

          if (!factor_scan) {
            logger << "Completed current prime field in "
            << std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - prime_start).count())
            << " s | " + std::to_string(total_iterations) << " probes in total\n"
            << "Average time of the black-box probe: " << std::to_string(average_black_box_time) + " s\n\n"
            << "Promote to new prime field: F(" << std::to_string(primes()[prime_it]) << ")\n";
            logger.close();
            logger.open("firefly.log", std::ios_base::app);

            if (verbosity > SILENT) {
              INFO_MSG("Completed current prime field in " +
                     std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - prime_start).count()) +
                     " s | " + std::to_string(total_iterations) + " probes in total");
              INFO_MSG("Average time of the black-box probe: " + std::to_string(average_black_box_time) + " s\n");
              INFO_MSG("Promote to new prime field: F(" + std::to_string(primes()[prime_it]) + ")");
            }
          }
        }

        prime_start = std::chrono::high_resolution_clock::now();

        reset_new_prime();
        FFInt::set_new_prime(primes()[prime_it]);
        bb.prime_changed_internal();

        if (!parsed_factors.empty()) {
          for (auto& el : parsed_factors) {
            el.second.precompute_tokens();
          }
        }

        // if only a small constant is reconstructed it will not ask for new run
        if (probes_for_next_prime == 0) {
#ifndef WITH_MPI
          probes_for_next_prime = thr_n ;//* bunch_size;
#else
          probes_for_next_prime = buffer * worker_thread_count + thr_n; // TODO start even more? * bunch_size
#endif
        }

        if (!safe_mode && (!save_states || (save_states && !load_anchor_points)) && !tmp_rec.need_shift(prime_it)) {
          if (tmp_rec.get_zi_shift_vec() != std::vector<FFInt> (n, 0)) {
            logger << "Disable shift\n";
            logger.close();
            logger.open("firefly.log", std::ios_base::app);

            if (verbosity > SILENT)
              INFO_MSG("Disable shift");

            tmp_rec.disable_shift();
          }
        }

        // Load anchor points and the shift to resume from saved probes
        if (save_states && load_anchor_points) {
          load_anchor_points = false;
          std::string line;
          std::ifstream anchor_point_file;
          anchor_point_file.open("ff_save/anchor_points");

          if (anchor_point_file.is_open()) {
            std::getline(anchor_point_file, line);
            tmp_rec.set_anchor_points(parse_vector_FFInt(line, static_cast<int>(n)));
            const auto tmp_an_vec = tmp_rec.get_anchor_points();

            for (auto & rec : reconst) {
              std::get<2>(rec)->set_anchor_points(tmp_an_vec);
            }
          } else {
            logger << "Anchor point file not found!\n";
            ERROR_MSG("Anchor point file not found!");
            logger.close();
            std::exit(EXIT_FAILURE);
          }

          anchor_point_file.close();

          std::ifstream shift_file;
          shift_file.open("ff_save/shift");

          if (shift_file.is_open()) {
            std::getline(shift_file, line);
            tmp_rec.set_shift(parse_vector_FFInt(line, static_cast<int>(n)));
          } else {
            logger << "Shift file not found!\n";
            logger.close();
            ERROR_MSG("Shift file not found!");
            std::exit(EXIT_FAILURE);
          }

          for (auto & rec : reconst) {
            if (!std::get<2>(rec)->is_done() && std::get<2>(rec)->get_prime() > prime_it) {
              ++items_new_prime;
            }
          }
        } else if (anchor_points.size() > static_cast<std::size_t>(prime_it)) {
          std::vector<FFInt> tmp_an_vec;

          for (const auto value : anchor_points[prime_it]) {
            tmp_an_vec.emplace_back(value);
          }

          tmp_rec.set_anchor_points(tmp_an_vec);

          for (auto & rec : reconst) {
            std::get<2>(rec)->set_anchor_points(tmp_an_vec);
          }
        } else {
          tmp_rec.generate_anchor_points();
          const auto tmp_an_vec = tmp_rec.get_anchor_points();

          for (auto & rec : reconst) {
            std::get<2>(rec)->set_anchor_points(tmp_an_vec);
          }
        }

        if (shifts.size() > static_cast<std::size_t>(prime_it)) {
          shift.clear();
          shift.reserve(n);

          for (const auto value : shifts[prime_it]) {
            shift.emplace_back(value);
          }

          tmp_rec.set_shift(shift);
        }

        shift = tmp_rec.get_zi_shift_vec();

        if (save_states) {
          std::remove("ff_save/anchor_points");
          std::ofstream file;
          file.open("ff_save/anchor_points");
          std::string tmp_str = "";

          for (const auto & el : tmp_rec.get_anchor_points()) {
            tmp_str += std::to_string(el.n) + " ";
          }

          tmp_str += std::string("\n");
          file << tmp_str;
          file.close();

          std::remove("ff_save/shift");
          file.open("ff_save/shift");
          tmp_str = "";

          for (const auto & el : tmp_rec.get_zi_shift_vec()) {
            tmp_str += std::to_string(el.n) + " ";
          }

          tmp_str += std::string("\n");
          file << tmp_str;
          file.close();
        }

        // start only thr_n jobs first, because the reconstruction can be done after the first feed
        uint32_t to_start = 0;

#ifndef WITH_MPI
        if (probes_for_next_prime > thr_n /** bunch_size*/) {
          to_start = thr_n /** bunch_size*/;

          {
#else
        if (mpi_first_send) {
          mpi_first_send = false;
          send_first_jobs(); // TODO send only start

          queue_probes(zi_order, thr_n); // * bunch_size
          started_probes[zi_order] += thr_n; // * bunch_size

          new_jobs = true;
          continue_communication = true;
          cond_val.notify_one();

          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          cond_val.wait(lock_probe_queue, [this](){return proceed;});

          proceed = false;
        } else {
          if (probes_for_next_prime > buffer * worker_thread_count + thr_n) { // * bunch_size
            to_start = buffer * worker_thread_count + thr_n; // * bunch_size
#endif

            if (verbosity == CHATTY) {
              INFO_MSG("Starting " + std::to_string(to_start) + " jobs now, the remaining " + std::to_string(probes_for_next_prime - to_start) + " jobs will be started later");
            }

#ifndef WITH_MPI
            queue_probes(zi_order, to_start);
#else
            queue_probes(zi_order, to_start, true);
#endif
            started_probes.emplace(zi_order, to_start);
#ifndef WITH_MPI
          }
        } else {
          {
#else
          } else {
#endif
            to_start = probes_for_next_prime;

            if (verbosity == CHATTY) {
              INFO_MSG("Starting " + std::to_string(to_start) + " jobs");
            }

#ifndef WITH_MPI
            queue_probes(zi_order, to_start);
#else
            queue_probes(zi_order, to_start, true);
#endif
            started_probes.emplace(zi_order, to_start);
          }

#ifndef WITH_MPI
        }

        if (precomputed_probes) {
          load_precomputed_probes_from_file();
        }
#else
          cond_val.notify_one();

          {
            std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

            cond_val.wait(lock_probe_queue, [this](){return proceed;});

            proceed = false;
          }

          for (uint32_t j = 0; j != to_start; ++j) {
            tp.run_task([this](uint32_t thread_id) {
              get_job(thread_id);
            });
          }
        }
#endif

        probes_for_next_prime = 0;
      }

      std::vector<uint64_t> indices;
      std::vector<std::vector<FFInt>> probes;

      get_probe(indices, probes);

      if (verbosity > SILENT && !factor_scan && !scan && std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - last_print_time).count() > 2.) {
        std::lock_guard<std::mutex> lock_print(print_control);
        last_print_time = std::chrono::high_resolution_clock::now();
        std::cerr << "\033[1;34mFireFly info:\033[0m Probe: " << probes_fed << "\r";
      }

      {
        std::lock_guard<std::mutex> lock_feed(feed_control);

        ++feed_jobs;
      }

      tp.run_priority_task([this, indices = std::move(indices), probes = std::move(probes)](uint32_t thread_id) {
        (void)thread_id;
        feed_job(indices, probes);
      });

      {
        std::lock_guard<std::mutex> lock_status(status_control);

        if (items_done == items) {
          done = true;
          continue;
        } else if (items_done + items_new_prime == items) {
          new_prime = true;
          continue;
        }
      }

      {
        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        if (probes_queued == 0) {
          lock_probe_queue.unlock();

          {
            std::unique_lock<std::mutex> lock_feed(feed_control);

            while (feed_jobs > 0 || interpolate_jobs > 0) {
              condition_feed.wait(lock_feed);

              lock_feed.unlock();
              lock_probe_queue.lock();

              if (probes_queued != 0) {
                lock_probe_queue.unlock();

                break;
              }

              lock_probe_queue.unlock();
              lock_feed.lock();
            }

            lock_probe_queue.lock();

            if (probes_queued != 0) {
              continue;
            }

            lock_probe_queue.unlock();
          }

          // make sure that no threads are running anymore
          while (tp.wait());

          // no jobs are running anymore, check if done or new_prime else throw error
          if (items_done == items) {
            done = true;
            continue;
          } else if (items_done + items_new_prime == items) {
            new_prime = true;
            continue;
          } else if (precomputed_probes) {
            write_requested_probes_to_file();
            return;
          } else {
            std::string msg = "Nothing left to feed: "
                              + std::to_string(items)
                              + " " + std::to_string(items_new_prime)
                              + " " + std::to_string(items_done) + " | "
                              + std::to_string(feed_jobs) + " "
                              + std::to_string(interpolate_jobs) + " | "
                              + std::to_string(probes_fed) + " "
                              + std::to_string(balance_of_ones) + " | "
                              + std::to_string(probes_queued) + " "
                              + std::to_string(computed_probes.size()) + " "
                              + std::to_string(requested_probes.size());
            WARNING_MSG(msg);
            WARNING_MSG("Please report this error");
            WARNING_MSG("Attempting to continue");
            logger << msg << "\nPlease report this error\nAttempting to continue\n";

            attempt_to_continue();
          }
        } else if (precomputed_probes) {
          lock_probe_queue.unlock();
          std::unique_lock<std::mutex> lock_future(future_control);

          if (computed_probes.empty()) {
            lock_future.unlock();

            {
              std::unique_lock<std::mutex> lock_feed(feed_control);

              condition_feed.wait(lock_feed, [this](){return feed_jobs == 0 && interpolate_jobs == 0;});
            }

            {
              std::lock_guard<std::mutex> lock_status(status_control);

              if (items_done == items) {
                done = true;
                continue;
              } else if (items_done + items_new_prime == items) {
                new_prime = true;
                continue;
              }
            }

            write_requested_probes_to_file();
            return;
          }
        }
      }
    }

#ifdef WITH_MPI
    {
      std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

      requested_probes = std::deque<std::pair<uint64_t, std::vector<FFInt>>>();
      new_jobs = false;

      cond_val.wait(lock_probe_queue, [this](){return proceed;});

      proceed = false;
    }
#endif

    tp.kill_all();

#ifdef WITH_MPI
    new_jobs = false; // to make sure that no new jobs have been started
#endif

    if (!factor_scan && !scan && save_states) {
      for (uint32_t item = 0; item != items; ++item) {
        // remove probe files if the interpolation is done
        std::remove(("ff_save/probes/" + std::to_string(item) + "_" + std::to_string(prime_it) + ".gz").c_str());
        ogzstream gzfile;
        std::string file_name = "ff_save/probes/" + std::to_string(item) + "_" + std::to_string(prime_it + 1) + ".gz";
        gzfile.open(file_name.c_str());
        gzfile.close();

        if (prime_it) {
          std::string file_name_old = "ff_save/states/" + std::to_string(item) + "_" + std::to_string(prime_it - 1) + ".gz";
          std::string file_name_new = "ff_save/states/" + std::to_string(item) + "_" + std::to_string(prime_it) + ".gz";
          std::rename(file_name_old.c_str(), file_name_new.c_str());
        }
      }
    }

#ifdef WITH_MPI
    {
      std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

      continue_communication = true;
      cond_val.notify_one();

      cond_val.wait(lock_probe_queue, [this](){return proceed;});

      proceed = false;
    }
#endif

    total_iterations += iteration;
    reset_new_prime();
  }

  // TODO optimize for bunch_size 1?
  template<typename BlackBoxTemp>
#ifndef WITH_MPI
  void Reconstructor<BlackBoxTemp>::queue_probes(const std::vector<uint32_t>& zi_order, const uint32_t to_start) {
#else
  void Reconstructor<BlackBoxTemp>::queue_probes(const std::vector<uint32_t>& zi_order, const uint32_t to_start, const bool first) {
#endif
    bool ones = false;

    if ((prime_it == 0 || safe_mode == true) && (zi_order == std::vector<uint32_t> (n - 1, 1) || (factor_scan && zi_order == std::vector<uint32_t> (0, 1)))) {
      ones = true;
    }

    std::vector<FFInt> rand_zi;
    rand_zi.reserve(zi_order.size());

    if (factor_scan) {
      rand_zi = rand_zi_fac;
    } else {
      rand_zi = tmp_rec.get_rand_zi_vec(zi_order);
    }

    for (uint32_t j = 0; j != to_start; ++j) {
      std::vector<FFInt> values(n);
      FFInt t = tmp_rec.get_rand_64();

      // check if t was already used for this zi_order
      {
        std::lock_guard<std::mutex> chosen_lock(chosen_mutex);

        auto it = chosen_t.find(zi_order);

        if (it != chosen_t.end()) {
          auto itt = it->second.find(t.n);

          if (itt != it->second.end()) {
            --j;
            continue;
          } else {
            it->second.emplace(t.n);
          }
        } else {
          chosen_t.emplace(std::make_pair(zi_order, std::unordered_set<uint64_t>( {t.n})));
        }
      }

      if (!factor_scan) {
        if (change_var_order) {
          for (const auto & el : optimal_var_order) {
            if (el.first == 0) {
              values[el.second] = t + shift[0];
            } else {
              values[el.second] = rand_zi[el.first - 1] * t + shift[el.first];
            }
          }
        } else {
          values[0] = t + shift[0];

          for (uint32_t i = 1; i != n; ++i) {
            values[i] = rand_zi[i - 1] * t + shift[i];
          }
        }
      } else {
        for (uint32_t i = 0; i != n; ++i) {
          if (rand_zi_fac[i] == 1) {
            values[i] = t;
          } else {
            values[i] = FFInt(rand_zi_fac[i]);
          }
        }
      }

      std::lock_guard<std::mutex> lock_probe_queue(mutex_probe_queue);

      if (ones) {
        requested_probes.emplace_front(std::make_pair(ind, std::move(values)));
      } else {
        requested_probes.emplace_back(std::make_pair(ind, std::move(values)));
      }

      index_map.emplace(std::make_pair(ind, std::make_pair(t, zi_order)));
      ++ind;

      ++probes_queued;
    }

#ifdef WITH_MPI
    if (!first) {
#endif
    if (!precomputed_probes) {
      for (uint32_t j = 0; j != to_start; ++j) {
        tp.run_task([this](uint32_t thread_id) {
          get_job(thread_id);
        });
      }
    }
#ifdef WITH_MPI
    }
#endif

#ifdef WITH_MPI
    new_jobs = true;
#endif
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::get_probe(std::vector<uint64_t>& indices, std::vector<std::vector<FFInt>>& probes) {
    {
      std::unique_lock<std::mutex> lock_future(future_control);

      condition_future.wait(lock_future, [this](){return computed_probes.size() != 0;});

      indices = std::move(computed_probes.front().first);
      probes = std::move(computed_probes.front().second);
      computed_probes.pop();
    }

    probes_fed += indices.size();

    std::lock_guard<std::mutex> lock_probe_queue(mutex_probe_queue);

    probes_queued -= static_cast<uint32_t>(indices.size());
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::feed_job(const std::vector<uint64_t>& indices, const std::vector<std::vector<FFInt>>& probes) {
    {
      std::lock_guard<std::mutex> lock(feed_control);

      interpolate_jobs += items;
    }

    std::vector<FFInt> t_vec;
    t_vec.reserve(indices.size());
    std::vector<std::vector<uint32_t>> zi_order_vec;
    zi_order_vec.reserve(indices.size());

    uint32_t count_ones = 0;

    {
      std::lock_guard<std::mutex> lock_probe_queue(mutex_probe_queue);

      for (const auto& index : indices) {
        auto tmp = std::move(index_map[index]);
        index_map.erase(index);
        t_vec.emplace_back(tmp.first);
        zi_order_vec.emplace_back(std::move(tmp.second));

        if ((factor_scan && static_cast<uint32_t>(zi_order_vec.back().size()) != 0) || (!factor_scan && static_cast<uint32_t>(zi_order_vec.back().size()) != n - 1)) {
          logger << "zi_order of probe has wrong length: " << std::to_string(zi_order_vec.back().size()) << "\n";
          ERROR_MSG("zi_order of probe has wrong length: " + std::to_string(zi_order_vec.back().size()));
          logger.close();
          std::exit(EXIT_FAILURE);
        }

        if ((prime_it == 0 || safe_mode == true) && (zi_order_vec.back() == std::vector<uint32_t>(n - 1, 1) || (factor_scan && zi_order_vec.back() == std::vector<uint32_t>(0, 1)))) {
          ++count_ones;
        }
      }
    }

    if (count_ones != 0) {
      std::lock_guard<std::mutex> lock_status(job_control);

      balance_of_ones += count_ones;
    }

    uint32_t counter = 0;

    for (auto & rec : reconst) {
      if (std::get<1>(rec) == RECONSTRUCTING) {
        std::pair<bool, uint32_t> done_prime = std::get<2>(rec)->get_done_and_prime();

        if (!done_prime.first) {
          if (done_prime.second == prime_it) {
            auto interpolate_and_write = std::get<2>(rec)->feed(t_vec, probes[std::get<0>(rec)], zi_order_vec, prime_it);

            if (interpolate_and_write.first) {
              ++counter;

              tp.run_priority_task([this, &rec](uint32_t thread_id) {
                (void)thread_id;
                interpolate_job(rec);
              });
            }

            if (interpolate_and_write.second) {
              tp.run_priority_task([&rec](uint32_t thread_id) {
                (void)thread_id;
                std::get<2>(rec)->write_food_to_file();
              });
            }
          }
        }
      }
    }

    {
      std::lock_guard<std::mutex> lock_status(status_control);

      if (!factor_scan && !scan && (one_done || one_new_prime)) {
        one_done = false;
        one_new_prime = false;

        std::lock_guard<std::mutex> lock_print(print_control);
        logger << "Probe: " + std::to_string(probes_fed)
        << " | Done: " << std::to_string(items_done) << " / " << std::to_string(items)
        << " | " << "Requires new prime field: " << std::to_string(items_new_prime)
        << " / " << std::to_string(items - items_done) << "\n";
        logger.close();
        logger.open("firefly.log", std::ios_base::app);

        if (verbosity > SILENT) {
          INFO_MSG("Probe: " + std::to_string(probes_fed) +
                   " | Done: " + std::to_string(items_done) + " / " + std::to_string(items) +
                   " | " + "Requires new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items - items_done));
        }
      }
    }

    std::lock_guard<std::mutex> lock(feed_control);

    interpolate_jobs -= (items - counter);
    --feed_jobs;
    condition_feed.notify_one();
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::interpolate_job(RatReconst_tuple& rec) {
    if (std::get<1>(rec) == RECONSTRUCTING) {
      std::tuple<bool, bool, uint32_t> interpolated_done_prime = std::get<2>(rec)->interpolate();

      if (std::get<0>(interpolated_done_prime)) { // interpolated
        // start new jobs if required
        if (!std::get<1>(interpolated_done_prime)) { // not done
          if (std::get<2>(interpolated_done_prime) > prime_it) { // compare prime counters
            {
              std::lock_guard<std::mutex> lock_status(status_control);

              one_new_prime = true;
              ++items_new_prime;
            }

            std::lock_guard<std::mutex> lock(job_control);

            if (std::get<2>(rec)->get_num_eqn() > probes_for_next_prime) {
              probes_for_next_prime = std::get<2>(rec)->get_num_eqn();
            }
          } else if (!safe_mode && prime_it != 0) {
            std::vector<std::pair<uint32_t, uint32_t>> all_required_probes = std::get<2>(rec)->get_needed_feed_vec();

            if (!all_required_probes.empty()) {
              for (const auto & some_probes : all_required_probes) {
                if (some_probes.second == 0) {
                  continue;
                }

                uint32_t required_probes = some_probes.second;

                std::vector<uint32_t> zi_order;
                if (!factor_scan) {
                  zi_order = std::vector<uint32_t>(n - 1, some_probes.first);
                } else {
                  zi_order = std::vector<uint32_t>(0, some_probes.first);
                }

                std::unique_lock<std::mutex> lock(job_control);

                auto it = started_probes.find(zi_order);

                if (it != started_probes.end()) {
                  if (required_probes > started_probes[zi_order]) {
                    uint32_t to_start = required_probes - started_probes[zi_order];

                    started_probes[zi_order] = required_probes;

                    lock.unlock();

                    if (verbosity == CHATTY) {
                      std::string msg = "Starting zi_order (";

                      for (const auto & ele : zi_order) {
                        msg += std::to_string(ele) + ", ";
                      }

                      msg = msg.substr(0, msg.length() - 2);
                      msg += ") " + std::to_string(to_start) + " time(s)";

                      std::lock_guard<std::mutex> lock_print(print_control);

                      INFO_MSG(msg);
                    }

                    queue_probes(zi_order, to_start);
                  }
                } else {
                  started_probes.emplace(zi_order, required_probes);

                  lock.unlock();

                  if (verbosity == CHATTY) {
                    std::string msg = "Starting zi_order (";

                    for (const auto & ele : zi_order) {
                      msg += std::to_string(ele) + ", ";
                    }

                    msg = msg.substr(0, msg.length() - 2);
                    msg += ") " + std::to_string(required_probes) + " time(s)";

                    std::lock_guard<std::mutex> lock_print(print_control);

                    INFO_MSG(msg);
                  }

                  queue_probes(zi_order, required_probes);
                }
              }
            }
          } else {
            std::pair<std::vector<std::pair<std::vector<uint32_t>, uint32_t>>, uint32_t> next_orders_pair = std::get<2>(rec)->get_zi_orders();
            auto next_orders = next_orders_pair.first;
            uint32_t tmp_system_size = next_orders_pair.second;

            if ((prime_it == 0 || safe_mode == true) && (factor_scan || (next_orders.size() == 1 && next_orders.front().first == std::vector<uint32_t>(n - 1, 1)))) {
              std::unique_lock<std::mutex> lock(job_control);

              if (balance_of_ones) {
                uint32_t to_start = balance_of_ones;
                balance_of_ones = 0;
                started_probes[next_orders.front().first] += to_start;

                lock.unlock();

                if (verbosity == CHATTY) {
                  std::lock_guard<std::mutex> lock_print(print_control);

                  INFO_MSG("Starting ones: " + std::to_string(to_start));
                }

                queue_probes(next_orders.front().first, to_start);
              }
            } else {
              for (size_t i = 0; i != next_orders.size(); ++i) {
                std::unique_lock<std::mutex> lock(job_control);

                auto it = started_probes.find(next_orders[i].first);

                if (it != started_probes.end()) {
                  if (tmp_system_size > started_probes[next_orders[i].first]) {
                    uint32_t to_start = std::min(tmp_system_size - started_probes[next_orders[i].first], next_orders[i].second);

                    started_probes[next_orders[i].first] += to_start;

                    lock.unlock();

                    if (verbosity == CHATTY) {
                      std::string msg = "Starting zi_order (";

                      for (const auto & ele : next_orders[i].first) {
                        msg += std::to_string(ele) + ", ";
                      }

                      msg = msg.substr(0, msg.length() - 2);
                      msg += ") " + std::to_string(to_start) + " time(s) ";

                      std::lock_guard<std::mutex> lock_print(print_control);

                      INFO_MSG(msg);
                    }

                    queue_probes(next_orders[i].first, to_start);
                  }
                } else {
                  started_probes.emplace(next_orders[i].first, next_orders[i].second);

                  lock.unlock();

                  if (verbosity == CHATTY) {
                    std::string msg = "Starting zi_order (";

                    for (const auto & ele : next_orders[i].first) {
                      msg += std::to_string(ele) + ", ";
                    }

                    msg = msg.substr(0, msg.length() - 2);
                    msg += ") " + std::to_string(next_orders[i].second) + " time(s)";

                    std::lock_guard<std::mutex> lock_print(print_control);

                    INFO_MSG(msg);
                  }

                  queue_probes(next_orders[i].first, next_orders[i].second);
                }
              }
            }
          }
        } else {
          // to be sure that no other thread does the same
          std::lock_guard<std::mutex> lock_status(status_control);

          if (std::get<1>(rec) == RECONSTRUCTING) {
            std::get<1>(rec) = DONE;

            ++items_done;
            one_done = true;
          }
        }
      }
    }

    std::lock_guard<std::mutex> lock(feed_control);

    --interpolate_jobs;
    condition_feed.notify_one();
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::clean_reconst() {
    std::lock_guard<std::mutex> lock_clean(clean);

    auto it = reconst.begin();

    while (it != reconst.end()) {
      if (std::get<1>(*it) == DELETE) {
        // delete RatReconst
        delete std::get<2>(*it);

        // remove from list
        it = reconst.erase(it);
      } else {
        ++it;
      }
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::get_job(uint32_t thread_id) {
    std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

    if (!requested_probes.empty()) {
      switch(compute_bunch_size(static_cast<uint32_t>(requested_probes.size()), thr_n, bunch_size)) {
        case 1:
          compute_probe(lock_probe_queue, thread_id);
          break;
        case 2:
          compute_probe<2>(lock_probe_queue, thread_id);
          break;
        case 4:
          compute_probe<4>(lock_probe_queue, thread_id);
          break;
        case 8:
          compute_probe<8>(lock_probe_queue, thread_id);
          break;
        case 16:
          compute_probe<16>(lock_probe_queue, thread_id);
          break;
        case 32:
          compute_probe<32>(lock_probe_queue, thread_id);
          break;
        case 64:
          compute_probe<64>(lock_probe_queue, thread_id);
          break;
        case 128:
          compute_probe<128>(lock_probe_queue, thread_id);
          break;
/*        case 256:
          compute_probe<256>(lock_probe_queue);
          break;*/
      }

#ifdef WITH_MPI
      if (requested_probes.empty()) {
        new_jobs = false;
      }
#endif
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::compute_probe(std::unique_lock<std::mutex>& lock_probe_queue, uint32_t thread_id) {
    std::vector<uint64_t> indices;
    indices.reserve(1);

    indices.emplace_back(requested_probes.front().first);

    std::vector<FFInt> values_vec(std::make_move_iterator(requested_probes.front().second.begin()), std::make_move_iterator(requested_probes.front().second.end()));

    requested_probes.pop_front();

    lock_probe_queue.unlock();

    auto time0 = std::chrono::high_resolution_clock::now();

    std::vector<FFInt> probe = bb.eval(values_vec, thread_id);//todo if for user-defined probes

    auto time1 = std::chrono::high_resolution_clock::now();

    auto time = std::chrono::duration<double>(time1 - time0).count();

    std::vector<std::vector<FFInt>> tmp;
    tmp.reserve(probe.size());

    for (size_t i = 0; i != probe.size(); ++i) {
      if (!factor_scan && parsed_factors.find(i) != parsed_factors.end()) {
        auto res = parsed_factors[i].evaluate_pre(values_vec);
        probe[i] /= res[0];
      }

      tmp.emplace_back(std::vector<FFInt> (1, probe[i]));
    }

    std::lock_guard<std::mutex> lock(future_control);

    iteration += 1;
#ifndef WITH_MPI
    int tmp_iterations = total_iterations + iteration;
    average_black_box_time = (average_black_box_time * (tmp_iterations - 1) + time) / tmp_iterations;
#else
    iterations_on_this_node += 1;
    int tmp_iterations = total_iterations + iterations_on_this_node;
    average_black_box_time = (average_black_box_time * (tmp_iterations - 1) + time) / tmp_iterations;
#endif

    computed_probes.emplace(std::make_pair(std::move(indices), std::move(tmp)));

    condition_future.notify_one();
  }

  template<typename BlackBoxTemp>
  template<uint32_t N>
  void Reconstructor<BlackBoxTemp>::compute_probe(std::unique_lock<std::mutex>& lock_probe_queue, uint32_t thread_id) {
    std::vector<uint64_t> indices;
    indices.reserve(N);

    std::vector<FFIntVec<N>> values_vec(n);

    for (uint32_t i = 0; i != N; ++i) {
      indices.emplace_back(requested_probes.front().first);

      for (uint32_t j = 0; j != n; ++j) {
        values_vec[j][i] = requested_probes.front().second[j];
      }

      requested_probes.pop_front();
    }

    lock_probe_queue.unlock();

    auto time0 = std::chrono::high_resolution_clock::now();

    std::vector<FFIntVec<N>> probe = bb.eval(values_vec, thread_id);

    auto time1 = std::chrono::high_resolution_clock::now();

    auto time = std::chrono::duration<double>(time1 - time0).count();

    std::vector<std::vector<FFInt>> tmp;
    tmp.reserve(probe.size());

    for (size_t i = 0; i != probe.size(); ++i) {
      // Remove factor from bb result
      if (!factor_scan && parsed_factors.find(i) != parsed_factors.end()) {
        auto res = parsed_factors[i].evaluate_pre(values_vec);

        for (size_t j = 0; j != N; ++j) {
          probe[i][j] /= res[0][j];
        }
      }

      tmp.emplace_back(std::vector<FFInt>(std::make_move_iterator(probe[i].begin()), std::make_move_iterator(probe[i].end())));
    }

    std::lock_guard<std::mutex> lock(future_control);

    iteration += N;
#ifndef WITH_MPI
    int tmp_iterations = total_iterations + iteration;
    average_black_box_time = (average_black_box_time * (tmp_iterations - N) + time) / tmp_iterations;
#else
    iterations_on_this_node += N;
    int tmp_iterations = total_iterations + iterations_on_this_node;
    average_black_box_time = (average_black_box_time * (tmp_iterations - N) + time) / tmp_iterations;
#endif

    computed_probes.emplace(std::make_pair(std::move(indices), std::move(tmp)));

    condition_future.notify_one();
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::attempt_to_continue() {
    uint64_t items_done_tmp = 0;
    uint64_t items_done_tmp2 = 0;
    uint64_t items_new_prime_tmp = 0;
    uint64_t items_new_prime_tmp2 = 0;
    uint64_t interpolated = 0;

    std::unordered_map<std::vector<uint32_t>, uint32_t, UintHasher> started_probes_copy = std::move(started_probes);
    started_probes = std::unordered_map<std::vector<uint32_t>, uint32_t, UintHasher> {};

    for (auto & rec : reconst) {
      if (std::get<1>(rec) == DONE) {
        ++items_done_tmp;
        one_done = true;
      } else if (std::get<1>(rec) == RECONSTRUCTING) {
        std::pair<bool, uint32_t> done_prime = std::get<2>(rec)->get_done_and_prime();

        if (done_prime.first) {
         ++items_done_tmp2;
         std::get<1>(rec) = DONE;
         one_done = true;
        } else {
          if (done_prime.second != prime_it) {
           ++items_new_prime_tmp;
           one_new_prime = true;
          } else {
            // try to use all stored probes
            std::tuple<bool, bool, uint32_t> interpolated_done_prime = std::get<2>(rec)->interpolate();

            if (std::get<0>(interpolated_done_prime)) {
              ++interpolated;
            }

            // force to write probes to the file
            std::get<2>(rec)->write_food_to_file();

            // start new jobs if required
            if (!std::get<1>(interpolated_done_prime)) { // not done
              if (std::get<2>(interpolated_done_prime) > prime_it) { // compare prime counters
                one_new_prime = true;
                ++items_new_prime_tmp2;

                if (std::get<2>(rec)->get_num_eqn() > probes_for_next_prime) {
                  probes_for_next_prime = std::get<2>(rec)->get_num_eqn();
                }
              } else if (!safe_mode && prime_it != 0) {
                std::vector<std::pair<uint32_t, uint32_t>> all_required_probes = std::get<2>(rec)->get_needed_feed_vec();

                if (!all_required_probes.empty()) {
                  for (const auto & some_probes : all_required_probes) {
                    if (some_probes.second == 0) {
                      continue;
                    }

                    uint32_t required_probes = some_probes.second;
                    std::vector<uint32_t> zi_order;
                    if (!factor_scan) {
                      zi_order = std::vector<uint32_t>(n - 1, some_probes.first);
                    } else {
                      zi_order = std::vector<uint32_t>(0, some_probes.first);
                    }

                    auto it = started_probes.find(zi_order);

                    if (it != started_probes.end()) {
                      if (required_probes > started_probes[zi_order]) {
                        uint32_t to_start = required_probes - started_probes[zi_order];
                        started_probes[zi_order] = required_probes;

                        if (verbosity == CHATTY) {
                         std::string msg = "Starting zi_order (";

                         for (const auto & ele : zi_order) {
                           msg += std::to_string(ele) + ", ";
                         }

                         msg = msg.substr(0, msg.length() - 2);
                         msg += ") " + std::to_string(to_start) + " time(s)";

                         INFO_MSG(msg);
                        }

                        queue_probes(zi_order, to_start);
                      }
                    } else {
                      started_probes.emplace(zi_order, required_probes);

                      if (verbosity == CHATTY) {
                        std::string msg = "Starting zi_order (";

                        for (const auto & ele : zi_order) {
                          msg += std::to_string(ele) + ", ";
                        }

                        msg = msg.substr(0, msg.length() - 2);
                        msg += ") " + std::to_string(required_probes) + " time(s)";

                        INFO_MSG(msg);
                      }

                      queue_probes(zi_order, required_probes);
                    }
                  }
                }
              } else {
                std::pair<std::vector<std::pair<std::vector<uint32_t>, uint32_t>>, uint32_t> next_orders_pair = std::get<2>(rec)->get_zi_orders();
                auto next_orders = next_orders_pair.first;
                uint32_t tmp_system_size = next_orders_pair.second;

                if ((prime_it == 0 || safe_mode == true) && (factor_scan || (next_orders.size() == 1 && next_orders.front().first == std::vector<uint32_t>(n - 1, 1)))) {
                  auto it = started_probes.find(next_orders.front().first);

                  if (it == started_probes.end()) {
                    balance_of_ones = 0;
#ifndef WITH_MPI
                    uint32_t to_start = thr_n /** bunch_size*/;
#else
                    uint32_t to_start = static_cast<uint32_t>(buffer) * worker_thread_count + thr_n; //TODO * bunch_size
#endif
                    started_probes[next_orders.front().first] += to_start;

                    if (verbosity == CHATTY) {
                      INFO_MSG("Starting ones: " + std::to_string(to_start));
                    }

                    queue_probes(next_orders.front().first, to_start);
                  }
                } else {
                  for (size_t i = 0; i != next_orders.size(); ++i) {
                    auto it = started_probes.find(next_orders[i].first);

                    if (it != started_probes.end()) {
                      if (tmp_system_size > started_probes[next_orders[i].first]) {
                        uint32_t to_start = std::min(tmp_system_size - started_probes[next_orders[i].first], next_orders[i].second);
                        started_probes[next_orders[i].first] += to_start;

                        if (verbosity == CHATTY) {
                          std::string msg = "Starting zi_order (";

                          for (const auto & ele : next_orders[i].first) {
                            msg += std::to_string(ele) + ", ";
                          }

                          msg = msg.substr(0, msg.length() - 2);
                          msg += ") " + std::to_string(to_start) + " time(s) ";

                          INFO_MSG(msg);
                        }

                        queue_probes(next_orders[i].first, to_start);
                      }
                    } else {
                      started_probes.emplace(next_orders[i].first, next_orders[i].second);

                      if (verbosity == CHATTY) {
                        std::string msg = "Starting zi_order (";

                        for (const auto & ele : next_orders[i].first) {
                          msg += std::to_string(ele) + ", ";
                        }

                        msg = msg.substr(0, msg.length() - 2);
                        msg += ") " + std::to_string(next_orders[i].second) + " time(s)";

                        INFO_MSG(msg);
                      }

                      queue_probes(next_orders[i].first, next_orders[i].second);
                    }
                  }
                }
              }
            } else {
              std::get<1>(rec) = DONE;
              ++items_done_tmp2;
              one_done = true;
            }
          }
        }
      }
    }

    if (interpolated) {
      WARNING_MSG(std::to_string(interpolated) + " interpolations");
      logger << std::to_string(interpolated) + " interpolations\n";
    }

    if (items_done != items_done_tmp || items_done_tmp2 != 0) {
      WARNING_MSG("Some items were not counted as done: " + std::to_string(items_done) + " " + std::to_string(items_done_tmp) + " " + std::to_string(items_done_tmp2));
      logger << "Some items were not counted as done: " + std::to_string(items_done) + " " + std::to_string(items_done_tmp) + " " + std::to_string(items_done_tmp2) + "\n";
      items_done = items_done_tmp + items_done_tmp2;
    }

    if (items_new_prime != items_new_prime_tmp || items_new_prime_tmp2 != 0) {
      WARNING_MSG("Some items were not counted as new prime: " + std::to_string(items_new_prime) + " " + std::to_string(items_new_prime_tmp) + " " + std::to_string(items_new_prime_tmp2));
      logger << "Some items were not counted as new prime: " + std::to_string(items_new_prime) + " " + std::to_string(items_new_prime_tmp) + " " + std::to_string(items_new_prime_tmp2) + "\n";
      items_new_prime = items_new_prime_tmp + items_new_prime_tmp2;
    }

    if (items_done == items) {
      WARNING_MSG("All items done");
      logger << "All items done\n";
      done = true;
    } else if (items_done + items_new_prime == items) {
      WARNING_MSG("All items require next prime");
      logger << "All items require next prime\n";
      new_prime = true;
    } else {
      if (started_probes.size() == 0) {
        WARNING_MSG("Cannot continue");
        logger << "Cannot continue\n";
        throw std::runtime_error("Cannot continue");
      } else {
        std::vector<uint32_t> zi_order;

        if (!factor_scan) {
          zi_order = std::vector<uint32_t>(n - 1, 1);
        } else {
          zi_order = std::vector<uint32_t>(0, 1);
        }

        auto it = started_probes.find(zi_order);

        if (it != started_probes.end() && started_probes.size() == 1) {
          WARNING_MSG("Started only ones");
          logger << "Started only ones\n";
          started_probes_copy.at(it->first) = it->second;
          started_probes = std::move(started_probes_copy);
        } else {
          WARNING_MSG("Started probes other than ones: Deleting probe history");
          logger << "Started probes other than ones: Deleting probe history\n";
          started_probes.emplace(std::make_pair(zi_order, 0));
        }
      }
    }

    WARNING_MSG("Continue normal procedure");
    logger << "Continue normal procedure\n";
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::reset_new_prime() {
    iteration = 0;
    probes_fed = 0;
#ifdef WITH_MPI
    iterations_on_this_node = 0;
    new_jobs = false;
#endif

    balance_of_ones = 0;
    probes_queued = 0;
    started_probes.clear();
    index_map.clear();
    ind = 0;
    feed_jobs = 0;
    interpolate_jobs = 0;
    new_prime = false;
    items_new_prime = 0;
    one_done = false;
    one_new_prime = false;

    // Only reset chosen_t when not resuming from a saved state
    if (!load_anchor_points) {
      chosen_t.clear();
    }

    requested_probes = std::deque<std::pair<uint64_t, std::vector<FFInt>>>();
    computed_probes = std::queue<std::pair<std::vector<uint64_t>, std::vector<std::vector<FFInt>>>>();
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::load_precomputed_probes_from_file() {
    std::string file_name = "probes/" + std::to_string(prime_it) + ".gz";
    std::ifstream ifile(file_name.c_str());
    igzstream gzfile;
    gzfile.open(file_name.c_str());

    if (ifile.is_open()) {
      ifile.close();
      std::string line;

      while (std::getline(gzfile, line)) {
        std::size_t first_delim = line.find(" | ");

        if (first_delim == std::string::npos) {
          WARNING_MSG(file_name + " has a wrong format!");
          logger << file_name + " has a wrong format!\n";
          throw std::runtime_error(file_name + " has a wrong format!");
        }

        std::string zi_order_string = line.substr(0, first_delim + 1);
        std::string remainder = line.substr(first_delim + 3);

        std::size_t second_delim = remainder.find(" | ");

        if (second_delim == std::string::npos) {
          WARNING_MSG(file_name + " has a wrong format!");
          logger << file_name + " has a wrong format!\n";
          throw std::runtime_error(file_name + " has a wrong format!");
        }

        FFInt t(std::stoul(remainder.substr(0, second_delim)));

        std::vector<uint32_t> zi_order = parse_vector_32(zi_order_string, n - 1);

        if (zi_order.size() != n - 1) {
          WARNING_MSG(file_name + " has a wrong format!");
          logger << file_name + " has a wrong format!\n";
          throw std::runtime_error(file_name + " has a wrong format!");
        }

        size_t pos = 0;
        std::string delimiter = " ";
        std::vector<std::vector<FFInt>> probes {};

        std::string probe_string = remainder.substr(second_delim + 3);

        if (probe_string.back() != ' ') {
          probe_string.append(" ");
        }

        while ((pos = probe_string.find(delimiter)) != std::string::npos) {
          probes.emplace_back(std::vector<FFInt> (1, std::stoul(probe_string.substr(0, pos))));
          probe_string.erase(0, pos + 1);
        }

        probes.shrink_to_fit();

        computed_probes.emplace(std::make_pair(std::vector<uint64_t> (1, ind), probes));
        index_map.emplace(std::make_pair(ind, std::make_pair(t, zi_order)));
        ++ind;
        ++probes_queued;
      }
    } else {
      WARNING_MSG("Cannot find " + file_name + "!");
      logger << "Cannot find " + file_name + "!\n";
      std::exit(EXIT_FAILURE);
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::write_requested_probes_to_file() {
    INFO_MSG("Writing requested probes to requested_probes.gz.");
    logger << "Writing requested probes to requested_probes.gz.\n";

    std::string file_name = "requested_probes.gz";
    ogzstream gzfile;
    gzfile.open(file_name.c_str());

    for (const auto & probe : requested_probes) {
      auto t_and_zi_order = index_map.at(probe.first);

      for (const auto & zi : t_and_zi_order.second) {
        gzfile << zi << " ";
      }

      gzfile << "| " << t_and_zi_order.first << " |";

      for (const auto & value : probe.second) {
        gzfile << " " << value;
      }

      gzfile << "\n";
    }

    gzfile.close();
  }

#ifdef WITH_MPI
  template<typename BlackBoxTemp>
  inline void Reconstructor<BlackBoxTemp>::mpi_setup() {
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Bcast(&prime_it, 1, MPI_UINT32_T, master, MPI_COMM_WORLD);

    if (!parsed_factors.empty()) {
      std::vector<std::string> vars;
      vars.reserve(n);

      for (std::uint32_t i = 1; i != n + 1; ++i) {
        vars.emplace_back("x" + std::to_string(i));
      }

      std::vector<std::string> factors = get_factors_string(vars);
      uint32_t counter = 0;

      for (auto & factor : factors) {
        const int batch_size = 2147483645; // 2147483647 is the largest signed 32-bit integer
        const int split = static_cast<int>(factor.size()) / batch_size;

        if (split > 1) {
          int tmp;
          if (static_cast<int>(factor.size()) % batch_size == 0) {
            tmp = -split;
          } else {
            tmp = - 1 - split;
          }
          MPI_Bcast(&tmp, 1, MPI_INT, master, MPI_COMM_WORLD);
        } else if (split == 1 && static_cast<int>(factor.size()) != batch_size) {
          int tmp = - 1 - split;
          MPI_Bcast(&tmp, 1, MPI_INT, master, MPI_COMM_WORLD);
        }

        for (int i = 0; i != 1 + split; ++i) {
          if (i == split) {
            int amount = static_cast<int>(factor.size()) - i * batch_size;
            if (amount != 0) {
              MPI_Bcast(&amount, 1, MPI_INT, master, MPI_COMM_WORLD);
              MPI_Bcast(&factor[i * batch_size], amount, MPI_CHAR, master, MPI_COMM_WORLD);
            }
          } else {
            int amount = batch_size;
            MPI_Bcast(&amount, 1, MPI_INT, master, MPI_COMM_WORLD);
            MPI_Bcast(&factor[i * batch_size], amount, MPI_CHAR, master, MPI_COMM_WORLD);
          }
        }

        uint32_t function_number = counter;
        ++counter;
        MPI_Bcast(&function_number, 1, MPI_UINT32_T, master, MPI_COMM_WORLD);
      }
    }

    int end = -1;
    MPI_Bcast(&end, 1, MPI_INT, master, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::send_first_jobs() {
    mpi_setup();

    std::vector<uint32_t> zi_order;

    if (!factor_scan) {
      zi_order = std::vector<uint32_t>(n - 1, 1);
    } else {
      zi_order = std::vector<uint32_t>(0, 1);
    }

    std::lock_guard<std::mutex> lock_probe_queue(mutex_probe_queue);
    std::lock_guard<std::mutex> lock_job(job_control);

    if (started_probes.find(zi_order) != started_probes.end()) {
      started_probes.emplace(std::make_pair(zi_order, 0));
    }

    std::vector<FFInt> rand_zi;
    if (factor_scan) {
      rand_zi = rand_zi_fac;
    } else {
      rand_zi = tmp_rec.get_rand_zi_vec(zi_order);
    }

    std::lock_guard<std::mutex> chosen_lock(chosen_mutex);

    for (int i = 1; i != world_size; ++i) {
      uint64_t to_start;
      MPI_Recv(&to_start, 1, MPI_UINT64_T, i, RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      nodes.emplace(std::make_pair(i, to_start));
      worker_thread_count += static_cast<uint32_t>(to_start / buffer);

      probes_queued += to_start;
      started_probes[zi_order] += to_start;

      std::vector<uint64_t> values;
      values.reserve(static_cast<uint32_t>(to_start) * (n + 1));

      for (uint32_t ii = 0; ii != static_cast<uint32_t>(to_start); ++ii) {
        //values[ii * (n + 1)] = ind;
        values.emplace_back(ind);

        FFInt t = tmp_rec.get_rand_64();

        // check if t was already used for this zi_order
        auto it = chosen_t.find(zi_order);

        if (it != chosen_t.end()) {
          auto itt = it->second.find(t.n);

          if (itt != it->second.end()) {
            --ii;
            continue;
          } else {
            it->second.emplace(t.n);
          }
        } else {
          chosen_t.emplace(std::make_pair(zi_order, std::unordered_set<uint64_t>( {t.n})));
        }

        if (!factor_scan) {
          if (change_var_order) {
            std::vector<uint64_t> tmp_values (n);

            for (const auto & el : optimal_var_order) {
              if (el.first == 0) {
                tmp_values[el.second] = (t + shift[0]).n;
              } else {
                tmp_values[el.second] = (rand_zi[el.first - 1] * t + shift[el.first]).n;
              }
            }

            values.insert(values.end(), tmp_values.begin(), tmp_values.end());
          } else {
            //values[ii * (n + 1) + 1] = (t + shift[0]).n;
            values.emplace_back((t + shift[0]).n);

            for (uint32_t j = 1; j != n; ++j) {
              //values[ii * (n + 1) + j + 1] = (t * rand_zi[j - 1] + shift[j]).n;
              values.emplace_back((t * rand_zi[j - 1] + shift[j]).n);
            }
          }
        } else {
          for (uint32_t j = 0; j != n; ++j) {
            if (rand_zi[j] == 1) {
              values.emplace_back(t.n);
            } else {
              values.emplace_back(rand_zi_fac[j].n);
            }
          }
        }

        index_map.emplace(std::make_pair(ind, std::make_pair(t, zi_order)));
        ++ind;
      }

      MPI_Send(&values[0], static_cast<int>(static_cast<uint32_t>(to_start) * (n + 1)), MPI_UINT64_T, i, VALUES, MPI_COMM_WORLD);
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::mpi_communicate() {
    while (true) {
      int flag_ext = 0;
      MPI_Status status;

      bool restart_empty_nodes = false;

      while (!flag_ext) {
        MPI_Iprobe(MPI_ANY_SOURCE, RESULT, MPI_COMM_WORLD, &flag_ext, &status);

        if (new_prime || done) {
          break;
        } else if (!empty_nodes.empty() && new_jobs) {
          restart_empty_nodes = true;
          break;
        }
      }

      if (done && !scan) {
        MPI_Request* requests = new MPI_Request[world_size - 1];
        uint64_t tmp;

        for (int i = 1; i != world_size; ++i) {
          MPI_Isend(&tmp, 1, MPI_UINT64_T, i, END, MPI_COMM_WORLD, &requests[i - 1]);
        }

        MPI_Waitall(world_size - 1, requests, MPI_STATUSES_IGNORE);

        delete[] requests;

        MPI_Barrier(MPI_COMM_WORLD);

        {
          std::lock_guard<std::mutex> lock_probe_queue(mutex_probe_queue);

          proceed = true;

          cond_val.notify_one();
        }

        for (int i = 1; i != world_size; ++i) {
          while (true) {
            MPI_Status status_rec;
            MPI_Probe(i, RESULT, MPI_COMM_WORLD, &status_rec);

            int amount;
            MPI_Get_count(&status_rec, MPI_UINT64_T, &amount);

            std::vector<uint64_t> tmp_vec;
            tmp_vec.reserve(amount);

            MPI_Recv(&tmp_vec[0], amount, MPI_UINT64_T, status_rec.MPI_SOURCE, status_rec.MPI_TAG, MPI_COMM_WORLD, &status_rec);

            if (amount == 1 && tmp_vec[0] == 0) {
              break;
            }
          }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        std::vector<double> timings;
        std::vector<double> weights;
        timings.reserve(world_size);
        weights.reserve(world_size);
        double total_weight = 0.;

        for (int i = 1; i != world_size; ++i) {
          MPI_Status status_time;
          double timing[2];
          MPI_Recv(&timing, 2, MPI_DOUBLE, i, TIMING, MPI_COMM_WORLD, &status_time);
          timings.emplace_back(timing[0]);
          weights.emplace_back(timing[1]);
          total_weight += timing[1];
        }

        MPI_Barrier(MPI_COMM_WORLD);

        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          cond_val.wait(lock_probe_queue, [this](){return continue_communication;});

          continue_communication = false;
        }

        iteration = total_weight + iterations_on_this_node;

        timings.emplace_back(average_black_box_time);
        weights.emplace_back(static_cast<double>(total_iterations + iterations_on_this_node));
        total_weight += static_cast<double>(total_iterations + iterations_on_this_node);

        average_black_box_time = 0.;

        for (int i = 0; i != world_size; ++i) {
          average_black_box_time += weights[i] / total_weight * timings[i];
        }

        std::lock_guard<std::mutex> lock_probe_queue(mutex_probe_queue);

        proceed = true;

        cond_val.notify_one();

        break;
      } else if (new_prime || (done && scan)) {
        MPI_Request* requests = new MPI_Request[world_size - 1];

        // send the new-prime signal but not the actual prime
        uint64_t prime_tmp = 1;

        for (int i = 1; i != world_size; ++i) {
          MPI_Isend(&prime_tmp, 1, MPI_UINT64_T, i, NEW_PRIME, MPI_COMM_WORLD, &requests[i - 1]);
        }

        MPI_Waitall(world_size - 1, requests, MPI_STATUSES_IGNORE);

        delete[] requests;

        MPI_Barrier(MPI_COMM_WORLD);

        {
          std::lock_guard<std::mutex> lock_probe_queue(mutex_probe_queue);

          proceed = true;

          cond_val.notify_one();
        }

        for (int i = 1; i != world_size; ++i) {
          while (true) {
            MPI_Status status_rec;
            MPI_Probe(i, RESULT, MPI_COMM_WORLD, &status_rec);

            int amount;
            MPI_Get_count(&status_rec, MPI_UINT64_T, &amount);

            std::vector<uint64_t> tmp_vec;
            tmp_vec.reserve(amount);

            MPI_Recv(&tmp_vec[0], amount, MPI_UINT64_T, status_rec.MPI_SOURCE, status_rec.MPI_TAG, MPI_COMM_WORLD, &status_rec);

            if (amount == 1 && tmp_vec[0] == 0) {
              break;
            }
          }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        std::vector<double> timings;
        std::vector<double> weights;
        timings.reserve(world_size);
        weights.reserve(world_size);
        double total_weight = 0.;

        for (int i = 1; i != world_size; ++i) {
          MPI_Status status_time;
          double timing[2];
          MPI_Recv(&timing, 2, MPI_DOUBLE, i, TIMING, MPI_COMM_WORLD, &status_time);
          timings.emplace_back(timing[0]);
          weights.emplace_back(timing[1]);
          total_weight += timing[1];
        }

        MPI_Barrier(MPI_COMM_WORLD);

        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          cond_val.wait(lock_probe_queue, [this](){return continue_communication;});

          continue_communication = false;
        }

        iteration = total_weight + iterations_on_this_node;

        timings.emplace_back(average_black_box_time);
        weights.emplace_back(static_cast<double>(total_iterations + iterations_on_this_node));
        total_weight += static_cast<double>(total_iterations + iterations_on_this_node);

        average_black_box_time = 0.;

        for (int i = 0; i != world_size; ++i) {
          average_black_box_time += weights[i] / total_weight * timings[i];
        }

        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        proceed = true;

        cond_val.notify_one();

        cond_val.wait(lock_probe_queue, [this](){return new_jobs.load();});

        proceed = true;

        cond_val.notify_one();

        // send the actual new prime
        requests = new MPI_Request[world_size - 1];

        if (!factor_scan) {
          prime_tmp = static_cast<uint64_t>(prime_it);
        } else {
          prime_tmp = prime_it_fac;
        }

        for (int i = 1; i != world_size; ++i) {
          MPI_Isend(&prime_tmp, 1, MPI_UINT64_T, i, NEW_PRIME, MPI_COMM_WORLD, &requests[i - 1]);
        }

        MPI_Waitall(world_size - 1, requests, MPI_STATUSES_IGNORE);

        delete[] requests;

        MPI_Barrier(MPI_COMM_WORLD);

        empty_nodes = std::queue<std::pair<int, uint64_t>>();

        for (int k = 1; k != world_size; ++k) {
          uint64_t free_slots;
          MPI_Recv(&free_slots, 1, MPI_UINT64_T, k, RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          if (requested_probes.size() != 0) {
            uint32_t size = compute_job_number(static_cast<uint32_t>(requested_probes.size()), static_cast<uint32_t>(free_slots), thr_n, bunch_size);

            std::vector<uint64_t> values;
            values.reserve(size * (n + 1));

            for (uint64_t i = 0; i != size; ++i) {
              //values[i * (n + 1)] = requested_probes.front().first;
              values.emplace_back(requested_probes.front().first);

              for (uint64_t j = 0; j != n; ++j) {
                //values[i * (n + 1) + 1 + j] = requested_probes.front().second[j].n;
                values.emplace_back(requested_probes.front().second[j].n);
              }

              requested_probes.pop_front();
            }

            if (requested_probes.empty()) {
              new_jobs = false;
            }

            MPI_Send(&values[0], static_cast<int>(size * (n + 1)), MPI_UINT64_T, k, VALUES, MPI_COMM_WORLD);
          } else if (free_slots == nodes[k]) {
            empty_nodes.emplace(std::make_pair(k, free_slots));
          }
        }
      } else if (restart_empty_nodes) {
        std::lock_guard<std::mutex> lock_probe_queue(mutex_probe_queue);

        while (!empty_nodes.empty() && !requested_probes.empty()) {
          auto node = empty_nodes.front();
          empty_nodes.pop();

          uint32_t size = compute_job_number(static_cast<uint32_t>(requested_probes.size()), static_cast<uint32_t>(node.second), thr_n, bunch_size);

          std::vector<uint64_t> values;
          values.reserve(size * (n + 1));

          for (uint64_t i = 0; i != size; ++i) {
            //values[i * (n + 1)] = requested_probes.front().first;
            values.emplace_back(requested_probes.front().first);

            for (uint64_t j = 0; j != n; ++j) {
              //values[i * (n + 1) + 1 + j] = requested_probes.front().second[j].n;
              values.emplace_back(requested_probes.front().second[j].n);
            }

            requested_probes.pop_front();
          }

          MPI_Send(&values[0], static_cast<int>(size * (n + 1)), MPI_UINT64_T, node.first, VALUES, MPI_COMM_WORLD);
        }

        if (requested_probes.empty()) {
          new_jobs = false;
        }
      } else {
        int amount;
        MPI_Get_count(&status, MPI_UINT64_T, &amount);

        if ((static_cast<uint32_t>(amount) - 1) % (items + 1) != 0) {
          logger << "Corrupted results recieved: " + std::to_string(amount - 1) << "\n";
          ERROR_MSG("Corrupted results recieved: " + std::to_string(amount - 1));
          logger.close();
          std::exit(EXIT_FAILURE);
        }

        uint32_t new_results = (static_cast<uint32_t>(amount) - 1) / (items + 1);

        // TODO optimize the format
        std::vector<uint64_t> results_list;
        results_list.reserve(amount);
        MPI_Recv(&results_list[0], amount, MPI_UINT64_T, status.MPI_SOURCE, RESULT, MPI_COMM_WORLD, &status);

        std::vector<uint64_t> indices;
        indices.reserve(new_results);
        indices.emplace_back(results_list[0]);

        std::vector<std::vector<FFInt>> results;
        results.reserve(items);

        for (uint32_t j = 1; j != items + 1; ++j) {
          std::vector<FFInt> result;
          result.reserve(new_results);
          result.emplace_back(results_list[j]);
          results.emplace_back(result);
        }

        for (uint32_t i = 1; i != new_results; ++i) {
          indices.emplace_back(results_list[i * (items + 1)]);

          for (uint32_t j = 1; j != items + 1; ++j) {
            results[j - 1].emplace_back(results_list[i * (items + 1) + j]);
          }
        }

        {
          std::lock_guard<std::mutex> lock_res(future_control);

          computed_probes.emplace(std::make_pair(std::move(indices), std::move(results)));

          condition_future.notify_one();
        }

        uint64_t free_slots = results_list[amount - 1];

        std::lock_guard<std::mutex> lock_probe_queue(mutex_probe_queue);

        if (!requested_probes.empty()) {
          uint32_t size = compute_job_number(static_cast<uint32_t>(requested_probes.size()), static_cast<uint32_t>(free_slots), thr_n, bunch_size);

          std::vector<uint64_t> values;
          values.reserve(size * (n + 1));

          for (uint64_t i = 0; i != size; ++i) {
            //values[i * (n + 1)] = requested_probes.front().first;
            values.emplace_back(requested_probes.front().first);

            for (uint64_t j = 0; j != n; ++j) {
              //values[i * (n + 1) + 1 + j] = requested_probes.front().second[j].n;
              values.emplace_back(requested_probes.front().second[j].n);
            }

            requested_probes.pop_front();
          }

          if (requested_probes.empty()) {
            new_jobs = false;
          }

          MPI_Send(&values[0], static_cast<int>(size * (n + 1)), MPI_UINT64_T, status.MPI_SOURCE, VALUES, MPI_COMM_WORLD);
        } else if (free_slots == nodes[status.MPI_SOURCE]) {
          empty_nodes.emplace(std::make_pair(status.MPI_SOURCE, free_slots));
        } else {
          MPI_Send(NULL, 0, MPI_UINT64_T, status.MPI_SOURCE, VALUES, MPI_COMM_WORLD);
        }
      }
    }
  }
#endif
}
