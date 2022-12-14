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
#include "firefly/PolynomialFF.hpp"
#include "firefly/RationalNumber.hpp"

namespace firefly {
  /**
   *    Applies the cinese remainder theorem
   *    @param a1 the first coefficient
   *    @param m1 the first modulus
   *    @param a2 the second coefficient
   *    @param m2 the second modulus
   *    @param m2inv n_preinvert_limb(m2)
   *    @returns the combination of a1 and a2 corresponding to the chinese
   *    remainder theorem
   */
  std::pair<fmpzxx, fmpzxx> run_chinese_remainder(
    const fmpzxx &a1, const fmpzxx &m1, ulong a2, ulong m2, ulong m2inv);
  /**
   *    Applies the rational reconstruction algorithm
   *    @param a a number over a finite field
   *    @param p a prime number defining the finite field
   *    @return a RationalNumber which has been reconstruction using the
   *    rational reconstruction algorithm
   */
  std::pair<bool, RationalNumber> get_rational_coef(const fmpzxx& a, const fmpzxx& p);

   /**
   *    Applies the rational reconstruction algorithm MQRR from
   *    Maximal Quotient Rational Reconstruction: An Almost Optimal Algorithm for Rational Reconstruction
   *    by M. Monagan
   *    @param a a number over a finite field
   *    @param p a prime number defining the finite field
   *    @return a RationalNumber which has been reconstruction using the
   *    rational reconstruction algorithm
   */
  std::pair<bool, RationalNumber> get_rational_coef_mqrr(const fmpzxx& a, const fmpzxx& p);

  /**
  *  Solves the given modified Vandermonde system
  *  @param degs the contributing degrees
  *  @param nums the evaluated numerical values
  *  @param val the anchor points
  *  @return the polynomial
  */
  PolynomialFF solve_vandermonde_system(std::vector<std::vector<uint32_t>>& degs,
                                        const std::vector<FFInt>& nums,
                                        const std::vector<FFInt>& val);
  /**
   *  Compares two vetors colexographically, i.e. (1,0,0) < (0,1,0), and returns
   *  true if the first arguement is greater than the second
   *  @param a first vector which should be probed if it its greater
   *  @param b the reference vector for the comparison
   */
  bool a_grt_b(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b);

  /**
   *  Compares two vetors colexographically, i.e. (1,0,0) < (0,1,0), and returns
   *  true if the first arguement is greater than the second. This function
   *  is particulary written for tuples of 0 and 1
   *  @param a first vector which should be probed if it its greater
   *  @param b the reference vector for the comparison
   */
  bool a_grt_b_s(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b);
  /**
   *  Checks whether two vectors are equal
   *  @param a first vector which should be probed if it its greater
   *  @param b the reference vector for the comparison
   */
  bool a_eq_b(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b);
  /**
   *  Generates the next binary permutation
   *  @param curr_per current permutation
   *  @return a pair where the first bool indicates if the second entry is a new permutation
   */
  std::pair<bool, std::vector<uint32_t>> generate_next_permutation(std::vector<uint32_t>& curr_per);
  /**
   *  Compute the bunch size for the next probes
   *  @param queue_length length of the probes queued
   *  @param thr_n number of threads available
   *  @param max_bunch_size the maximum bunch size allowed
   *  @return bunch_size
   */
  uint32_t compute_bunch_size(const uint32_t queue_length, const uint32_t thr_n, const uint32_t bunch_size);
  /**
   *  Distribute probes for a given number threads
   *  @param queue_length length of the probes queued
   *  @param thr_n number of threads
   *  @param total_threads total number of threads available
   *  @param max_bunch_size the maximum bunch size allowed
   *  @return the number of probes which have been assigned
   */
  uint32_t compute_job_number(const uint32_t queue_length, const uint32_t threads, const uint32_t total_threads, const uint32_t bunch_size);
}
