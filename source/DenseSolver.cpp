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
#include "Logger.hpp"

namespace firefly {
  void calc_inverse(mat_ff& a, uint32_t n_) {
    int n = static_cast<int>(n_);
    // Augment a with unit matrix
    for (int i = 0; i < n; ++i) {
      std::vector<FFInt> dum(n);
      dum[i] = 1;
      std::vector<FFInt> tmp = a[i];
      tmp.insert(tmp.end(), dum.begin(), dum.end());
      a[i] = tmp;
    }

    for (int i = 0; i < n; ++i) {
      // search for maximum in this column
      FFInt max_el = a[i][i];
      int max_row = i;

      for (int k = i + 1; k < n; ++k) {
        if (a[k][i] > max_el) {
          max_el = a[k][i];
          max_row = k;
        }
      }

      // swap maximum row with current row (column by column)
      for (int k = i; k < 2 * n; ++k) {
        FFInt tmp = a[max_row][k];
        a[max_row][k] = a[i][k];
        a[i][k] = tmp;
      }

      // make all rows below this one 0 in current column
      for (int k = i + 1; k < n; ++k) {
        FFInt c = -a[k][i] / a[i][i];

        for (int j = i; j < 2 * n; ++j) {
          if (i == j)
            a[k][j] = 0;
          else
            a[k][j] += c * a[i][j];
        }
      }
    }

    // solve equation ax=b for an upper triangular matrix a
    mat_ff res(n, std::vector<FFInt>(n));

    for (int i = n - 1; i >= 0; i--) {
      for (int k = n; k < 2 * n; ++k) {
        a[i][k] /= a[i][i];
      }

      for (int row_mod = i - 1; row_mod >= 0; row_mod--) {
        for (int column_mod = n; column_mod < 2 * n; ++column_mod) {
          a[row_mod][column_mod] -= a[i][column_mod] * a[row_mod][i];
        }
      }

      // Remove the unit matrix from the result
      for (int k = n; k < 2 * n; ++k) {
        res[i][k - n] = a[i][k];
      }
    }

    a = res;
  }

  std::vector<FFInt> solve_gauss_system(mat_ff& a, uint32_t n) {
    // Transform the matrix in upper triangular form
    for (uint32_t i = 0; i < n; ++i) {
      // search for maximum in this column
      FFInt max_el = a[i][i];
      uint32_t max_row = i;

      for (uint32_t k = i + 1; k < n; ++k) {
        FFInt tmp = a[k][i];

        if (tmp.n > max_el.n) {
          max_el = tmp;
          max_row = k;
        }
      }

      // swap maximum row with current row (column by column)
      for (uint32_t k = i; k < n + 1; ++k) {
        FFInt tmp = a[max_row][k];
        a[max_row][k] = a[i][k];
        a[i][k] = tmp;
      }

      // Make all rows below this one zero in the current column
      for (uint32_t k = i + 1; k < n; ++k) {
        FFInt c = -a[k][i] / a[i][i];

        for (uint32_t j = i; j < n + 1; ++j) {
          if (i == j) a[k][j] = FFInt(0);
          else a[k][j] += c * a[i][j];
        }
      }
    }

    std::vector<FFInt> results(n);

    if (a[n - 1][n - 1] != 0) {
      // Solve equation A * x = b for an upper triangular matrix
      for (int i = n - 1; i >= 0; i--) {
        results[i] = a[i][n] / a[i][i];

        for (int k = i - 1; k >= 0; k--) {
          a[k][n] -= a[k][i] * results[i];
        }
      }
    } else {
      ERROR_MSG("Singular system of equations!");
      std::exit(-1);
    }

    return results;
  }

  void calc_lu_decomposition(mat_ff& a, std::vector<int>& p, uint32_t n) {
    uint32_t i, k, j, max_row;
    FFInt max_el;
    std::vector<FFInt> tmp;
    p = std::vector<int> (n + 1);

    for (i = 0; i <= n; ++i)
      p[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < n; ++i) {
      max_el = 0;
      max_row = i;

      for (k = i; k < n; ++k)
        if (a[k][i] > max_el) {
          max_el = a[k][i];
          max_row = k;
        }

      if (max_row != i) {
        //pivoting p
        j = p[i];
        p[i] = p[max_row];
        p[max_row] = j;

        //pivoting rows of A
        tmp = a[i];
        a[i] = a[max_row];
        a[max_row] = tmp;

        //counting pivots starting from N (for determinant)
        p[n]++;
      }

      for (j = i + 1; j < n; ++j) {
        a[j][i] /= a[i][i];

        for (k = i + 1; k < n; ++k)
          a[j][k] -= a[j][i] * a[i][k];
      }
    }
  }

  void calc_inverse_lu(const mat_ff& a, mat_ff& ia, std::vector<int>& p, uint32_t n) {
    ia = mat_ff(n, std::vector<FFInt> (n));

    for (uint32_t j = 0; j < n; ++j) {
      for (uint32_t i = 0; i < n; ++i) {
        if (p[i] == static_cast<int>(j))
          ia[i][j] = 1;
        else
          ia[i][j] = 0;

        for (uint32_t k = 0; k < i; ++k)
          ia[i][j] -= a[i][k] * ia[k][j];
      }

      for (int i = n - 1; i >= 0; i--) {
        for (int k = i + 1; k < static_cast<int>(n); ++k)
          ia[i][j] -= a[i][k] * ia[k][j];

        ia[i][j] = ia[i][j] / a[i][i];
      }
    }
  }

  FFInt calc_determinant_lu(const mat_ff& a, std::vector<int>& p, uint32_t n) {
    FFInt det = a[0][0];

    for (uint32_t i = 1; i < n; ++i)
      det *= a[i][i];

    if ((p[n] - n) % 2 == 0)
      return det;
    else
      return -det;
  }

  std::vector<FFInt> solve_lu(mat_ff& a, std::vector<int>& p, const std::vector<FFInt>& b, uint32_t n) {
    std::vector<FFInt> x(n);

    for (uint32_t i = 0; i < n; ++i) {
      x[i] = b[p[i]];

      for (uint32_t k = 0; k < i; ++k)
        x[i] -= a[i][k] * x[k];
    }

    for (int i = n - 1; i >= 0; i--) {
      for (int k = i + 1; k < static_cast<int>(n); ++k)
        x[i] -= a[i][k] * x[k];

      x[i] = x[i] / a[i][i];
    }

    return x;
  }
}

