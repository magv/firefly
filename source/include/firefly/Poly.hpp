//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2020  Jonas Klappert, Sven Yannick Klein, and Fabian Lange
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
#include "firefly/FFInt.hpp"

namespace firefly {

  /**
  * @class Poly
  * @brief a class for univariate polynomials over finite fields using FireFly integers
  */

  class Poly {

  public:
    /**
    * Default constructor
    */
    Poly();
    /**
    * Constructor for a vector of FFInts
    * @param coeff_vector the vector of the FFInts used as coefficients for the new polynomial, first entry is constant term
    */
    Poly(std::vector<FFInt>& coeff_vector);
    /**
    * Copy constructor
    * @param old_poly the polynomial to be copied
    */
    Poly(const Poly& old_poly);
    /**
    * Destructor of the class
    */
    ~Poly();

    std::vector<FFInt> coeff; // The coefficients stored as a vector of FFInts, the first entry corresponds to the constant term
    /**
    * Returns the degree of the polynomial
    * @return the degree of the polynomial ignoring zeros
    */
    size_t get_deg() const;
    /**
    * shrinks the coeff vector to the smallest size possible, deletes zeros
    */
    void shrink_to_fit();
    /**
    * reverses the polynomial by reversing the coeff vector
    */
    void rev();
    /**
    * finds the roots of the polynomial, only works if the polynomial completly factorizes in linear factors
    * @return a vector of the found roots
    */
    std::vector<FFInt> roots();

    Poly& operator=(const Poly&) = default;
    Poly& operator-=(const Poly&);
    Poly& operator+=(const Poly&);
    Poly& operator*=(const FFInt&);
    Poly& operator/=(const FFInt&);
    Poly& operator*=(const Poly&);
  };

  Poly operator+(const Poly& a, const Poly& b);
  Poly operator-(const Poly& a, const Poly& b);
  Poly operator*(const Poly& a, const FFInt&);
  Poly operator/(const Poly& a, const FFInt&);
  Poly operator*(const Poly& a, const Poly& b);
  Poly operator/(const Poly& a, const Poly& b);
  Poly operator%(const Poly& a, const Poly& b);
  std::ostream& operator<<(std::ostream& out, const Poly& a);
  /**
  * divides the polynomials by fast Euclidean division
  * @param a the numerator polynomial
  * @param z the denominator polynomial
  * returns a pair of polynomials where the first is the quotient and second is the remainder
  */
  std::pair<Poly, Poly> fast_euclidean_division(const Poly& a, const Poly& z);
  /**
  * gives the greatest common divisor of two polynomials
  * @param a one of the polynomials
  * @param b the second polynomial
  * @return the greatest common divisor of the two polynomials
  */
  Poly gcd(const Poly& a, const Poly& b);
}
