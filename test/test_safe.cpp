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
    BlackBoxUser(const ShuntingYardParser& par_) : par(par_) {};

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

    virtual void prime_changed() {
      par.precompute_tokens();
    }

  private:
    ShuntingYardParser par;
  };
}

using namespace firefly;
int main() {
  INFO_MSG("Test safe mode");
  ShuntingYardParser p_1("../../parser_test/s_y_safe.m", {"x1", "y", "zZ", "W"});
  BlackBoxUser b_1(p_1);
  Reconstructor r_1(4, 4, b_1);
  r_1.set_safe_interpolation();
  r_1.reconstruct();
  RatReconst::reset();

  ShuntingYardParser p_1_2("../../parser_test/s_y_1_v.m", {"x"});
  BlackBoxUser b_1_2(p_1_2);
  Reconstructor r_1_2(1, 4, b_1_2);
  r_1_2.set_safe_interpolation();
  r_1_2.reconstruct();
  RatReconst::reset();

  ShuntingYardParser p_1_3("../../parser_test/s_y_4_v.m", {"x1", "y", "zZ", "W"});
  BlackBoxUser b_1_3(p_1_3);
  Reconstructor r_1_3(4, 4, b_1_3);
  r_1_3.set_safe_interpolation();
  r_1_3.reconstruct();
  INFO_MSG("Safe mode passed");

  return 0;
}
