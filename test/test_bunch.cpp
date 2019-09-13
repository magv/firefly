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

#include "DenseSolver.hpp"
#include "Reconstructor.hpp"
#include "ShuntingYardParser.hpp"
#include "tinydir.h"

namespace firefly {
  class BlackBoxUser : public BlackBoxBase {
  public:
    BlackBoxUser(const ShuntingYardParser& par_, int mode_) : par(par_), mode(mode_) {};

    virtual std::vector<FFInt> operator()(const std::vector<FFInt>& values) {
      //std::vector<FFInt> result;

      // Get results from parsed expressions
      std::vector<FFInt> result = par.evaluate_pre(values);

      result.emplace_back(result[0] / result[3]);

      // Build the matrix mat
      mat_ff mat = {{result[0], result[1]}, {result[2], result[3]}};
      std::vector<int> p {};
      // Compute LU decomposition of mat
      calc_lu_decomposition(mat, p, 2);
      // Compute determinant of mat
      result.emplace_back(calc_determinant_lu(mat, p, 2));

      return result;
    }

    virtual std::vector<std::vector<FFInt>> operator()(const std::vector<std::vector<FFInt>>& values) {
      // Get results from parsed expressions
      std::vector<std::vector<FFInt>> result = par.evaluate_pre(values);

      for (size_t i = 0; i != values.size(); ++i) {
        result[i].emplace_back(result[i][0] / result[i][3]);

        // Build the matrix mat
        mat_ff mat = {{result[i][0], result[i][1]}, {result[i][2], result[i][3]}};
        std::vector<int> p {};
        // Compute LU decomposition of mat
        calc_lu_decomposition(mat, p, 2);
        // Compute determinant of mat
        result[i].emplace_back(calc_determinant_lu(mat, p, 2));
      }

      return result;
    }

    virtual void prime_changed() {
      par.precompute_tokens();
    }

  private:
    ShuntingYardParser par;
    int mode = 0;
  };
}

using namespace firefly;
int main() {
  INFO_MSG("Test bunched evaluation");
  ShuntingYardParser p_3("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
  BlackBoxUser b_3(p_3, 3);
  Reconstructor r_3(4, 4, 4, b_3);
  r_3.enable_scan();
  r_3.reconstruct();
  RatReconst::reset();

  ShuntingYardParser p_3_2("../../parser_test/s_y_safe.m", {"x1", "y", "zZ", "W"});
  BlackBoxUser b_3_2(p_3_2, 3);
  Reconstructor r_3_2(4, 4, 4, b_3_2);
  r_3_2.set_safe_interpolation();
  r_3_2.reconstruct();
  RatReconst::reset();
  INFO_MSG("Bunched evaluation passed");

  return 0;
}
