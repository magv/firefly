//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2019  Jonas Klappert and Fabian Lange
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

#include "firefly/Reconstructor.hpp"
#include "firefly/AmplitudeParser.hpp"
#include "firefly/tinydir.h"
#include <sys/stat.h>

bool is_number(const std::string& s){
  return !s.empty() && std::find_if(s.begin(),
    s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}

using namespace firefly;
int main(int argc, char* argv[]) {
  auto time0 = std::chrono::high_resolution_clock::now();
  uint32_t n_threads = 1;
  uint32_t bs = 1;
  bool factor_scan = true;
  bool save_mode = false;
  bool no_interpolation = false;
  bool merge = false;
  bool stop_after_factors = false;
  std::string input_file = "";

  std::ofstream logger;
  logger.open("ff_insert.log");

  if (argc == 1) {
    ERROR_MSG("Please provide an input file as last argument");
    logger << "Please provide an input file as last argument\n";
    logger.close();
    std::exit(EXIT_FAILURE);
  }

  for (int i = 1; i != argc; ++i) {
    std::string arg = argv[i];

    if (arg == "-p" || arg == "--parallel") {
      if (i + 1 != argc) {
        if (is_number(argv[i + 1]))
        n_threads = std::stoi(argv[i + 1]);
        else {
          ERROR_MSG("The argument of -p needs to be a number");
          logger << "The argument of -p needs to be a number\n";
          logger.close();
          std::exit(EXIT_FAILURE);
        }
      } else {
        ERROR_MSG("-p needs an argument");
        logger << "-p needs an argument\n";
        logger.close();
        std::exit(EXIT_FAILURE);
      }

      ++i;
    } else if (arg == "-bs" || arg == "--bunchsize") {
      if (i + 1 != argc) {
        if (is_number(argv[i + 1]))
        bs = std::stoi(argv[i + 1]);
        else {
          ERROR_MSG("The argument of -bs needs to be a number");
          logger << "The argument of -bs needs to be a number\n";
          logger.close();
          std::exit(EXIT_FAILURE);
        }
      }  else {
        ERROR_MSG("-bs needs an argument");
        logger << "-bs needs an argument\n";
        logger.close();
        std::exit(EXIT_FAILURE);
      }

      ++i;
    } else if (arg == "-m" || arg == "--merge") {
      merge = true;
    } else if (arg == "-nfs" || arg == "--nofactorscan") {
      factor_scan = false;
    } else if (arg == "-ni" || arg == "--nointerpolation") {
      no_interpolation = true;
    } else if (arg == "-s" || arg == "--save") {
      save_mode = true;
    } else if (arg == "-fs" || arg == "--factorscan") {
      stop_after_factors = true;
    } else if (arg == "-h" || arg == "--help") {
      std::cerr << "Usage: " << argv[0] << "\n"
                << "Options:\n"
	        << "  -p,--parallel           Sets the number of used threads\n"
                << "  -bs,--bunchsize         Sets the maximum bunch size\n"
                << "  -fs,--factorscan        Stops after the factor scan and write out its results\n"
                << "  -m,--merge              Merges expressions in the given directory to one expression\n"
                << "  -nfs,--nofactorscan     Disables the factor scan\n"
                << "  -ni,--nointerpolation   Disables the interpolation and writes coefficients to files\n"
                << "  -s,--save               Enables the storage of intermediate results\n";
      return 1;
    } else if (i == argc - 1) {
      input_file = arg;
      std::ifstream test_file(input_file);

      if (!test_file.good()) {
        ERROR_MSG("Input file '" + input_file + "' does not exist");
        logger << "Input file '" + input_file + "' does not exist\n";
        logger.close();
        std::exit(EXIT_FAILURE);
      }
    } else {
      ERROR_MSG("Unknown option '" + arg + "'");
      logger << "Unknown option '" + arg + "'\n";
      logger.close();
      std::exit(EXIT_FAILURE);
    }
  }

  if (input_file.size() == 0) {
    ERROR_MSG("Please provide an input file as last argument");
    logger << "Please provide an input file as last argument\n";
    logger.close();
    std::exit(EXIT_FAILURE);
  }

  if (!merge) {
    // Parse integral families
    std::ifstream fam_file_test("config/functions");
    std::ifstream fam_file;

    if (!fam_file_test.good()) {
      ERROR_MSG("Please add a file definining the occurring functions in 'config/functions'");
      logger << "Please add a file definining the occurring functions in 'config/functions'\n";
      logger.close();
      std::exit(EXIT_FAILURE);
    }

    std::vector<std::string> families {};

    fam_file.open("config/functions");

    std::string line;

    while (std::getline(fam_file, line, '\n')) {
      line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

      if (line.size() != 0 && line[0] != '#')
        families.emplace_back(line);
    }

    fam_file.close();

    // Parse variables
    std::ifstream var_file_test("config/vars");
    std::ifstream var_file;

    if (!var_file_test.good()) {
      ERROR_MSG("Please add a file definining the occurring variables 'config/vars'");
      logger << "Please add a file definining the occurring variables 'config/vars'\n";
      logger.close();
      std::exit(EXIT_FAILURE);
    }

    std::vector<std::string> vars {};

    var_file.open("config/vars");

    while (std::getline(var_file, line, '\n')) {
      line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

      if (line.size() != 0 && line[0] != '#')
        vars.emplace_back(line);
    }

    var_file.close();

    if (vars.empty()) {
      ERROR_MSG("FireFly does not support functional reconstructions without variables!\n               If you want to continue, declare a variable.");
      logger << "FireFly does not support functional reconstructions without variables!\nIf you want to continue, declare a variable.\n";
      std::exit(EXIT_FAILURE);
    }

    // Check for a skip file
    std::ifstream skip_file_test("config/skip_functions");
    std::ifstream skip_file;
    std::unordered_set<std::string> skip_functions {};

    if (skip_file_test.good()) {
      skip_file.open("config/skip_functions");

      while (std::getline(skip_file, line, '\n')) {
        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

        if (line.size() != 0 && line[0] != '#')
          skip_functions.emplace(line);
      }

      skip_file.close();

      INFO_MSG("Skipping reconstruction of " + std::to_string(skip_functions.size()) + " basis function(s)\n");
      logger << "Skipping reconstruction of " << std::to_string(skip_functions.size()) << " basis function(s)\n\n";
    }

    logger.close();

    std::vector<std::string> expression_paths;
    std::vector<std::string> expression_files;
    std::unordered_map<std::string, std::string> paths_map;
    tinydir_dir expression_dir;
    tinydir_open_sorted(&expression_dir, input_file.c_str());

    if (expression_dir.n_files != 0) {
      if (input_file[input_file.size() - 1] == '/')
        input_file.erase(input_file.size() - 1 , 1);

      for (size_t i = 0; i != expression_dir.n_files; ++i) {
        tinydir_file file;
        tinydir_readfile_n(&expression_dir, &file, i);

        if (!file.is_dir) {
          expression_paths.emplace_back(input_file + "/" + file.name);
          expression_files.emplace_back(file.name);
        }
      }

      tinydir_close(&expression_dir);
    } else {
      expression_paths.emplace_back(input_file);
      size_t tmp_pos = input_file.rfind('/');

      if (tmp_pos != 0) {
        expression_files.emplace_back(input_file.substr(tmp_pos + 1, input_file.size() - 1 - tmp_pos));
      } else {
        expression_files.emplace_back(input_file);
      }
    }

    tinydir_close(&expression_dir);

    for (size_t j = 0; j != expression_paths.size(); ++j) {
      std::string expression = expression_paths[j];
      std::string fn = expression_files[j];
      // Construct the amplitude parser
      AmplitudeParser ap(vars, families);

      // Parse the amplitude
      ap.parse_amplitude_file(expression);

      // Parse IBP tables
      tinydir_dir dir;
      tinydir_open_sorted(&dir, "replacements");

      std::vector<std::string> files;

      for (size_t i = 0; i != dir.n_files; ++i) {
        tinydir_file file;
        tinydir_readfile_n(&dir, &file, i);

        if (!file.is_dir) {
          files.emplace_back(file.name);
        }
      }

      tinydir_close(&dir);

      for (const auto& file : files) {
        ap.parse_ibp_table_file("replacements/" + file);
      }

      // Build the black boxes and start reconstruction
      size_t masters = ap.check_for_unreplaced_masters();

      std::ofstream basis_f_file;
      basis_f_file.open("basis_functions");
      basis_f_file.close();

      bool skipped = false;

      if (!no_interpolation) {
        std::ofstream file;
        std::string file_name = "out_" + fn;
        file.open(file_name.c_str());
        file << "{\n";
        file.close();

        for (size_t i = 0; i != masters; ++i) {
          basis_f_file.open("basis_functions", std::ios_base::app);
          basis_f_file << ap.get_master(i) << "\n";
          basis_f_file.close();

          if (skip_functions.find(ap.get_master(i)) != skip_functions.end()) {
            INFO_MSG("Skipping basis function: " + ap.get_master(i));
            logger.open("ff_insert.log", std::ios_base::app);
            logger << "Skipping basis function: " + ap.get_master(i) + "\n\n";
            logger.close();
            skipped = true;
          } else {
            if (i == 0) {
              INFO_MSG("Reconstructing coefficient of basis function: " + ap.get_master(i) + "\n");
              logger.open("ff_insert.log", std::ios_base::app);
              logger << "Reconstructing coefficient of basis function: " + ap.get_master(i) + "\n";
              logger.close();
            }

            auto bb = ap.build_black_box(i);

            // Construct the reconstructor
            Reconstructor<FFAmplitudeBlackBox> reconst(bb.n, n_threads, bs, bb);

            // Settings
            if (factor_scan)
              reconst.enable_factor_scan();

	    if (stop_after_factors) {
	      reconst.stop_after_factor_scan();
	      save_mode = false;
	    }

            reconst.enable_shift_scan();

            bool renamed_ff_save = false;

            if (save_mode) {
              std::string tmp = "ff_save_" + ap.get_master(i);
              struct stat buffer;

              if (stat(tmp.c_str(), &buffer) == 0) {
                struct stat buffer_2;

                if (stat("ff_save", &buffer_2) == 0) {
                  std::rename("ff_save", "ff_save_tmp");
                  renamed_ff_save = true;
                }

                std::rename(tmp.c_str(), "ff_save");
              }

              reconst.set_tags({ap.get_master(i)});
              reconst.resume_from_saved_state();
            }

            // Reconstruct
            reconst.reconstruct();

	    if (!stop_after_factors) {
              std::vector<RationalFunction> results = reconst.get_result();
              file.open(file_name.c_str(), std::ios_base::app);

              if (!results.back().zero())
                file <<  "+ " << ap.get_master(i) << "*" + results.back().generate_horner(vars) << "\n";

              file.close();
	    } else {
	      std::string factor = reconst.get_factors_string(vars)[0];
              file.open(file_name.c_str(), std::ios_base::app);
              file <<  "+ " << ap.get_master(i) << "*" + factor << "\n";
              file.close();
	    }

            if (save_mode) {
              std::string tmp = "ff_save_" + ap.get_master(i);
              std::rename("ff_save", tmp.c_str());

              if (renamed_ff_save)
                std::rename("ff_save_tmp", "ff_save");
            }

            RatReconst::reset();

            std::ifstream in("firefly.log");

            if (in.good()) {
              std::string old_log((std::istreambuf_iterator<char>(in)),
                                  (std::istreambuf_iterator<char>()));
              logger.open("ff_insert.log", std::ios_base::app);
              logger << old_log << "\n";
              logger << "-------------------------------------------------------------\n\n";
              logger.close();
            }
          }

          if (i + 1 != masters) {
            std::cerr << "\n";
            INFO_MSG("Coefficients done: " + std::to_string(i + 1) + " / " + std::to_string(masters) + "\n");
            INFO_MSG("Reconstructing coefficient of basis function: " + ap.get_master(i + 1) + "\n");
            logger.open("ff_insert.log", std::ios_base::app);
            logger << "Coefficients done: " + std::to_string(i + 1) + " / " + std::to_string(masters) + "\n\n";
            logger << "-------------------------------------------------------------\n\n";
            logger << "Reconstructing coefficient of basis function: " + ap.get_master(i + 1) + "\n";
            logger.close();
          }

        }

        file.open(file_name.c_str(), std::ios_base::app);
        file << "}\n";
        file.close();

        std::cerr << "\n";

        if (masters == 1 && skipped) {
          std::remove("ff_insert.log");
          std::remove("basis_functions");
          std::remove("firefly.log");
          std::remove(file_name.c_str());
        } else {
          auto time1 = std::chrono::high_resolution_clock::now();
          INFO_MSG("Reconstructed expression in " + std::to_string(std::chrono::duration<double>(time1 - time0).count()) + " s");
          INFO_MSG("Result has been written to '" + file_name + "'");
          logger.open("ff_insert.log", std::ios_base::app);
          logger << "Coefficients done: " + std::to_string(masters) + " / " + std::to_string(masters) + "\n\n";
          logger << "-------------------------------------------------------------\n\n";
          logger << "Reconstructed expression in " + std::to_string(std::chrono::duration<double>(time1 - time0).count()) + " s\n";
          logger << "Result has been written to '" << file_name << "'\n";
          logger.close();
          std::string new_log_name = "ff_insert_" + fn + ".log";
          std::rename("ff_insert.log", new_log_name.c_str());
          std::string new_basis_functions_name = "basis_functions_" + fn;
          std::rename("basis_functions", new_basis_functions_name.c_str());
          std::rename("ff_insert.log", new_log_name.c_str());
          std::remove("firefly.log");
          time0 = std::chrono::high_resolution_clock::now();
        }
      } else {
	std::string dir_name = "coefficients_" + fn;
        mkdir(dir_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        for (size_t i = 0; i != masters; ++i) {
          basis_f_file.open("basis_functions", std::ios_base::app);
          basis_f_file << ap.get_master(i) << "\n";
          basis_f_file.close();

          if (skip_functions.find(ap.get_master(i)) != skip_functions.end()) {
            INFO_MSG("Skipping basis function: " + ap.get_master(i));
            logger.open("ff_insert.log", std::ios_base::app);
            logger << "Skipping basis function: " + ap.get_master(i) + "\n\n";
            logger.close();
            skipped = true;
          } else {
            std::ofstream coef_file;
            std::string file_name = "coefficients_" + fn + "/" + ap.get_master(i) + ".m";
            coef_file.open(file_name.c_str());
            coef_file << "{\n + " << ap.get_master(i) + "*(" << ap.get_unsimplified_coef(i) << ")\n}\n";
            coef_file.close();
          }
        }

        if (masters == 1 && skipped) {
          std::remove("ff_insert.log");
          std::remove("basis_functions");
          std::remove("firefly.log");
        } else {
          INFO_MSG("Unsimplified coefficients have been written to 'coefficients_" + fn + "' directory");
          logger.open("ff_insert.log", std::ios_base::app);
          logger << "Unsimplified coefficients have been written to 'coefficients_" << fn << "' directory\n";
          logger.close();
          std::string new_log_name = "ff_insert_" + fn + ".log";
          std::rename("ff_insert.log", new_log_name.c_str());
          std::string new_basis_functions_name = "basis_functions_" + fn;
          std::rename("basis_functions", new_basis_functions_name.c_str());
        }
      }
    }
  } else {
    if (input_file[input_file.size() - 1] == '/')
      input_file.erase(input_file.size() - 1 , 1);

    std::vector<std::string> files;

    tinydir_dir dir;
    tinydir_open_sorted(&dir, input_file.c_str());

    if (dir.n_files == 0) {
      ERROR_MSG("'" + input_file + "' is not a directory");
      logger << "'" << input_file << "' is not a directory\n";
      logger.close();
      std::exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i != dir.n_files; ++i) {
      tinydir_file file;
      tinydir_readfile_n(&dir, &file, i);

      if (!file.is_dir) {
        files.emplace_back(input_file + "/" + file.name);
      }
    }

    tinydir_close(&dir);

    if (files.size() == 0) {
      ERROR_MSG("Directory '" + input_file + "' has no content");
      logger << "Directory '" << input_file << "' has no content";
      logger.close();
      std::exit(EXIT_FAILURE);
    }

    INFO_MSG("Start merging " + std::to_string(files.size()) + " files");
    logger << "Start merging " + std::to_string(files.size()) + " files\n";
    logger.close();

    std::ofstream ofile;
    std::string merged_file_name = input_file + "_merged.out";
    ofile.open(merged_file_name.c_str());
    ofile << "{\n";

    size_t counter = 0;

    for (const auto& file : files) {
      std::ifstream infile;
      infile.open(file);
      std::string line;

      while (std::getline(infile, line, '\0')) {
        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        line.erase(std::remove(line.begin(), line.end(), ';'), line.end());

        if (line[0] == '{')
          line.erase(0, 1);

        if (line[line.size() - 1] == '}')
          line.erase(line.size() - 1, 1);

        ofile << line << "\n";
      }

      infile.close();
      std::cerr << "\033[1;34mFireFly info:\033[0m " << ++counter << " / " << files.size() << "\r";
    }

    ofile << "}\n";
    ofile.close();

    auto time1 = std::chrono::high_resolution_clock::now();
    INFO_MSG("Merged files in " + std::to_string(std::chrono::duration<double>(time1 - time0).count()) + " s");
    INFO_MSG("Result has been written to '" + merged_file_name + "'");
    logger.open("ff_insert.log", std::ios_base::app);
    logger << "Merged files in " + std::to_string(std::chrono::duration<double>(time1 - time0).count()) + " s\n";
    logger << "Result has been written to '" + merged_file_name + "'\n";
    logger.close();
  }

  return 0;
}
