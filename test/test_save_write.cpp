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

#include "firefly/DenseSolver.hpp"
#include "firefly/Reconstructor.hpp"
#include "firefly/ShuntingYardParser.hpp"
#include "firefly/tinydir.h"

#ifdef WITH_MPI
#include "firefly/MPIWorker.hpp"
#endif

namespace firefly {
  // Example of how one can use the black-box functor for the automatic interface
  class BlackBoxUser : public BlackBoxBase<BlackBoxUser> {
  public:
    BlackBoxUser(const ShuntingYardParser& par_, int mode_) : par(par_), mode(mode_) {};

    template<typename FFIntTemp>
    std::vector<FFIntTemp> operator()(const std::vector<FFIntTemp>& values, uint32_t thread_id) {
      (void)thread_id;

      //std::vector<FFInt> result;

      // Get results from parsed expressions
      std::vector<FFIntTemp> result = par.evaluate_pre(values);

      result.emplace_back(result[0] / result[3]);

      // Build the matrix mat
      mat_ff<FFIntTemp> mat = {{result[0], result[1]}, {result[2], result[3]}};
      std::vector<int> p {};
      // Compute LU decomposition of mat
      calc_lu_decomposition(mat, p, 2);
      // Compute determinant of mat
      result.emplace_back(calc_determinant_lu(mat, p, 2));

      return result;
    }

    inline void prime_changed() {
      par.precompute_tokens();
      c++;
#ifndef WITH_MPI

      if ((mode == 4 && c == 2) || (mode == 5 && c == 3)) {
        throw std::runtime_error("Abort for save test.");
      }

#else

      if (mode == 4 && c == 1) {
        std::exit(0);
      }

#endif
    }

  private:
    // Internal variables for the black box
    // In this example a ShuntingYardParser
    ShuntingYardParser par;
    int mode = 0;
    int c = 0;
  };
}

void remove_states() {
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

  for (const auto & file : files) {
    paths.emplace_back("ff_save/states/" + file);
  }

  for (const auto & el : paths) {
    std::remove(el.c_str());
  }

  std::remove("ff_save/states");
}

void remove_probes() {
  tinydir_dir dir;
  tinydir_open_sorted(&dir, "ff_save/probes");

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

  for (const auto & file : files) {
    paths.emplace_back("ff_save/probes/" + file);
  }

  for (const auto & el : paths) {
    std::remove(el.c_str());
  }

  std::remove("ff_save/probes");
}

using namespace firefly;
int main() {
  // Clean up if there is a save folder
  remove_states();
  remove_probes();
  std::remove("ff_save");

  std::string root_dir = FIREFLY_ROOT_DIR;
#ifndef WITH_MPI

  try {
    INFO_MSG("Test saving states and starting from them in prime 1");
    ShuntingYardParser p_4(root_dir + "/parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_4(p_4, 4);
    Reconstructor<BlackBoxUser> r_4(4, 4, b_4);
    r_4.enable_shift_scan();
    r_4.set_tags();
    r_4.reconstruct();
  } catch (std::exception& e) {
    RatReconst::reset();
    ShuntingYardParser p_5(root_dir + "/parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_5(p_5, 6);
    Reconstructor<BlackBoxUser> r_5(4, 4, b_5);
    r_5.set_tags();
    r_5.resume_from_saved_state();
    r_5.reconstruct();
    std::remove("ff_save/validation.gz");
    std::remove("ff_save/scan");
    std::remove("ff_save/shift");
    std::remove("ff_save/anchor_points");
    INFO_MSG("Starting from saved states passed");
  }

  // Remove files
  remove_states();
  remove_probes();
  std::remove("ff_save");
  RatReconst::reset();

  try {
    INFO_MSG("Test saving states and starting from them in prime 2");
    ShuntingYardParser p_4(root_dir + "/parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_4(p_4, 5);
    Reconstructor<BlackBoxUser> r_4(4, 4, b_4);
    r_4.enable_shift_scan();
    r_4.set_tags();
    r_4.reconstruct();
  } catch (std::exception& e) {
    RatReconst::reset();
    ShuntingYardParser p_5(root_dir + "/parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_5(p_5, 6);
    Reconstructor<BlackBoxUser> r_5(4, 4, b_5);
    r_5.set_tags();
    r_5.resume_from_saved_state();
    r_5.reconstruct();
    std::remove("ff_save/validation.gz");
    std::remove("ff_save/scan");
    std::remove("ff_save/shift");
    std::remove("ff_save/anchor_points");
  }

  // Remove files
  remove_states();
  remove_probes();
  std::remove("ff_save");
  RatReconst::reset();

  try {
    INFO_MSG("Test saving states and starting from them in prime 1 for 1 variable");
    ShuntingYardParser p_6(root_dir + "/parser_test/s_y_1_v.m", {"x"});
    BlackBoxUser b_6(p_6, 4);
    Reconstructor<BlackBoxUser> r_6(1, 4, b_6);
    r_6.enable_shift_scan();
    r_6.set_tags();
    r_6.reconstruct();
  } catch (std::exception& e) {
    RatReconst::reset();
    ShuntingYardParser p_7(root_dir + "/parser_test/s_y_1_v.m", {"x"});
    BlackBoxUser b_7(p_7, 6);
    Reconstructor<BlackBoxUser> r_7(1, 4, b_7);
    r_7.set_tags();
    r_7.resume_from_saved_state();
    r_7.reconstruct();
    std::remove("ff_save/validation.gz");
    std::remove("ff_save/scan");
    std::remove("ff_save/shift");
    std::remove("ff_save/anchor_points");
    INFO_MSG("Starting from saved states passed");
    std::cerr << "\n";
  }

  // Remove files
  remove_states();
  remove_probes();
  std::remove("ff_save");
  RatReconst::reset();
#else
  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_SERIALIZED, &provided);

  int process;
  MPI_Comm_rank(MPI_COMM_WORLD, &process);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  if (process == master) {
    INFO_MSG("Test saving states and starting from them in prime 1");
    ShuntingYardParser p_4(root_dir + "/parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_4(p_4, 4);
    Reconstructor<BlackBoxUser> r_4(4, 4, b_4);
    r_4.enable_shift_scan();
    r_4.set_tags();
    r_4.reconstruct();
  } else {
    ShuntingYardParser p_1(root_dir + "/parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
    BlackBoxUser b_1(p_1, 4);
    MPIWorker<BlackBoxUser>(4, std::thread::hardware_concurrency(), b_1);
  }

  MPI_Finalize();
#endif

  // Remove log file
  std::remove("firefly.log");

  return 0;
}
