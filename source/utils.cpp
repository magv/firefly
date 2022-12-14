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

#include "firefly/utils.hpp"
#include "firefly/RatReconst.hpp"
#include "firefly/ReconstHelper.hpp"

#include <algorithm>

namespace firefly {
  std::pair<fmpzxx, fmpzxx> run_chinese_remainder(
    const fmpzxx &a1, const fmpzxx &m1, ulong a2, ulong m2, ulong m2inv)
  {
    fmpzxx m1m2;
    fmpz_mul_ui(m1m2._fmpz(), m1._fmpz(), m2);
    ulong c = n_invmod(fmpz_fdiv_ui(m1._fmpz(), m2), m2);
    if (c == 0) abort();
    fmpzxx p1p2;
    _fmpz_CRT_ui_precomp(p1p2._fmpz(), a1._fmpz(), m1._fmpz(), a2, m2, m2inv, 0, c, 0);
    return std::pair<fmpzxx, fmpzxx>(p1p2, m1m2);
  }

  std::pair<bool, RationalNumber> get_rational_coef(const fmpzxx& a, const fmpzxx& p) {
    RationalNumber rn;
    int ok = _fmpq_reconstruct_fmpz(rn.numerator._fmpz(), rn.denominator._fmpz(), a._fmpz(), p._fmpz());
    return std::make_pair(ok, rn);
  }

  /* Implementation of MQRR from
   * Maximal Quotient Rational Reconstruction: An Almost Optimal Algorithm for Rational Reconstruction
   * by M. Monagan
   */
  std::pair<bool, RationalNumber> get_rational_coef_mqrr(const fmpzxx& u, const fmpzxx& p) {
    RationalNumber rn;
    // set to T so that less than one percent will be false positive results
    fmpzxx T = fmpzxx(1024 * p.sizeinbase(2));
    bool ok = true;

    if (u == 0) {
      if (p > T)
        return std::make_pair(ok, RationalNumber(0, 1));
      else
        return std::make_pair(false, RationalNumber(0, 1));
    }

    fmpzxx n = fmpzxx(0);
    fmpzxx d = fmpzxx(0);
    fmpzxx t0 = fmpzxx(0);
    fmpzxx r0 = p;
    fmpzxx t1 = fmpzxx(1);
    fmpzxx r1 = u;
    fmpzxx tmpr;
    fmpzxx tmpt;
    fmpzxx q;

    while ((!r1.is_zero()) && r0 > T) {
      q = r0 / r1; // fmpzxx automatically rounds off

      if (q > T) {
        n = r1;
        d = t1;
        T = q;
      }

      tmpr = r0;
      r0 = r1;
      r1 = tmpr - q * r1;
      tmpt = t0;
      t0 = t1;
      t1 = tmpt - q * t1;
    }

    if (d.is_zero() || gcd(n, d) != 1) {
      ok = false;
      n = 1;
      d = 1;
    }

    if (d < 0)
      return std::make_pair(ok, RationalNumber((-n).evaluate(), abs(d).evaluate()));

    return std::make_pair(ok, RationalNumber(n, d));
  }

  bool a_grt_b(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    for (int i = a.size() - 1; i != -1; --i) {
      if (a[i] == b[i])
        continue;
      else if (a[i] > b[i])
        return true;
      else
        return false;
    }

    return false;
  }

  bool a_eq_b(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    for (int i = a.size() - 1; i != -1; --i) {
      if (a[i] == b[i])
        continue;
      else
        return false;
    }

    return true;
  }

  bool a_grt_b_s(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    uint32_t deg1 = 0;
    uint32_t deg2 = 0;

    for (const auto & el : a) deg1 += el;

    for (const auto & el : b) deg2 += el;

    if (deg1 < deg2)
      return false;
    else if (deg1 == deg2)
      return a_grt_b(b, a);

    return true;
  }

  std::pair<bool, std::vector<uint32_t>> generate_next_permutation(std::vector<uint32_t>& curr_per) {
    uint32_t* curr_per_arr = &curr_per[0];
    int size = curr_per.size();

    bool got_next_per = std::next_permutation(curr_per_arr, curr_per_arr + size);

    if (got_next_per)
      return std::make_pair(true, std::vector<uint32_t> (curr_per_arr, curr_per_arr + size));
    else {
      int num_of_ones = 0;

      for (int i = 0; i != size; ++i) {
        if (curr_per[i] == 1)
          ++num_of_ones;
      }

      if (num_of_ones == size)
        return std::make_pair(false, std::vector<uint32_t> (size, 0));
      else {
        std::vector<uint32_t> new_per (size, 0);

        for (int i = size - 1; i >= size - 1 - num_of_ones; i--) {
          new_per[i] = 1;
        }

        return std::make_pair(true, new_per);
      }
    }
  }

  PolynomialFF solve_vandermonde_system(std::vector<std::vector<uint32_t>>& degs,
                                        const std::vector<FFInt>& nums,
  const std::vector<FFInt>& val) {
    uint32_t num_eqn = degs.size();
    uint32_t n = val.size();
    ff_map poly;
    poly.reserve(num_eqn);

    // calculate base entries of Vandermonde matrix
    std::vector<FFInt> vis;
    vis.reserve(num_eqn);

    for (const auto & el : degs) {
      FFInt vi = 1;

      // z_1 is always = 1 which does not matter while determining the coefficient
      for (uint32_t i = 0; i < n; ++i) {
        // curr_zi_ord starts at 1, thus we need to subtract 1 entry
        vi *= val[i].pow(el[i + 1]);
      }

      vis.emplace_back(vi);
    }

    // Initialize the coefficient vector of the master polynomial
    std::vector<FFInt> cis(num_eqn);

    // The coefficients of the master polynomial are found by recursion
    // where we have
    // P(Z) = (Z - v_0)*(Z - v_1)*...*(Z - v_{n-1})
    //      =  c_0 + c_1*Z + ... + Z^n
    cis[num_eqn - 1] = -vis[0];

    for (uint32_t i = 1; i < num_eqn; ++i) {
      for (uint32_t j = num_eqn - 1 - i; j < num_eqn - 1; ++j) {
        cis[j] -= vis[i] * cis[j + 1];
      }

      cis[num_eqn - 1] -= vis[i];
    }

    // Each subfactor in turn is synthetically divided,
    // matrix-multiplied by the right hand-side,
    // and supplied with a denominator (since all vi should be different,
    // there is no additional check if a coefficient in synthetical division
    // leads to a vanishing denominator)
    for (uint32_t i = 0; i < num_eqn; ++i) {
      FFInt t = 1;
      FFInt b = 1;
      FFInt s = nums[num_eqn - 1];

      for (int j = num_eqn - 1; j > 0; j--) {
        b = cis[j] + vis[i] * b;
        s += nums[j - 1] * b;
        t = vis[i] * t + b;
      }

      poly.emplace(std::make_pair(degs[i], s / t / vis[i]));
    }

    return PolynomialFF(n + 1, poly);
  }

  uint32_t compute_bunch_size(const uint32_t queue_length, const uint32_t thr_n, const uint32_t max_bunch_size) {
    if (max_bunch_size == 1) {
      return 1;
    } else {
      uint32_t tmp = queue_length / thr_n;

      if (tmp != 0) {
        // Compute the floor power of two
        tmp |= tmp >> 1;
        tmp |= tmp >> 2;
        tmp |= tmp >> 4;
        tmp |= tmp >> 8;
        tmp |= tmp >> 16;

        tmp = (tmp + 1) >> 1;

        if ((tmp << 1) < queue_length && tmp * thr_n != queue_length) {
          tmp <<= 1;
        }

        return std::min(max_bunch_size, tmp);
      } else {
        return 1;
      }
    }
  }

  uint32_t compute_job_number(const uint32_t queue_length, const uint32_t threads, const uint32_t total_threads, const uint32_t max_bunch_size) {
    uint32_t tmp_queue_length = queue_length;

    for (uint32_t i = 0; i != threads; ++i) {
      tmp_queue_length -= compute_bunch_size(tmp_queue_length, total_threads, max_bunch_size);

      if (tmp_queue_length == 0) {
        break;
      }
    }

    return queue_length - tmp_queue_length;
  }
}
