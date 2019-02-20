#pragma once
#include <gmpxx.h>
#include "RationalNumber.hpp"
#include "FFInt.hpp"

namespace firefly {
  /**
   *    Applies the cinese remainder theorem
   *    @param p1 a pair of a coefficient a and a prime p
   *    @param p2 a pair of a coefficient a and a prime p
   *    @returns the combination of p1 and p2 corresponding to the chinese
   *    remainder theorem
   */
  std::pair<mpz_class, mpz_class> run_chinese_remainder(
    const std::pair<mpz_class, mpz_class>& p1,
    const std::pair<mpz_class, mpz_class>& p2);

  /**
   *    Applies the rational reconstruction algorithm
   *    @param a a number over a finite field
   *    @param p a prime number defining the finite field
   *    @return a RationalNumber which has been reconstruction using the
   *    rational reconstruction algorithm
   */
  std::pair<bool, RationalNumber> get_rational_coef(const mpz_class& a, const mpz_class& p);

  /**
   *    Applies the rational reconstruction algorithm MQRR from
   *    Maximal Quotient Rational Reconstruction: An Almost Optimal Algorithm for Rational Reconstruction
   *    by M. Monagan
   *    @param a a number over a finite field
   *    @param p a prime number defining the finite field
   *    @return a RationalNumber which has been reconstruction using the
   *    rational reconstruction algorithm
   */
  std::pair<bool, RationalNumber> get_rational_coef_mqrr(const mpz_class& a, const mpz_class& p);

  /**
   *
   *
   */
  std::vector<FFInt> solve_gauss_system(uint32_t num_eqn,
                                        std::vector<std::vector<FFInt>>& coef_mat);

  bool a_grt_b(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b);

  bool a_grt_b_s(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b);

  std::vector<std::vector<uint32_t>> generate_possible_shifts(uint32_t r);
}
