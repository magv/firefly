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
#include "firefly/FFInt.hpp"

#include <iostream>
#include <string>

namespace firefly {
  /**
   * @class RationalNumber
   * @brief A container class representing rational numbers
   */
  class RationalNumber {
  public:
    /**
     *  Constructor of a RationalNumber object
     *  @param numerator_ the numerator as a fmpzxx
     *  @param denominator_ the denominator as a fmpzxx
     */
    RationalNumber(int numerator_, int denominator_) : RationalNumber(fmpzxx(numerator_), fmpzxx(denominator_)) {};
    RationalNumber(const fmpzxx& numerator_, const fmpzxx& denominator_);
    RationalNumber();
    RationalNumber operator*(const RationalNumber&);
    RationalNumber& operator+=(const RationalNumber& rn);
    RationalNumber& operator-=(const RationalNumber& rn);
    RationalNumber& operator*=(const RationalNumber& rn);
    bool operator==(const RationalNumber&) const;
    RationalNumber operator-() const;
    std::string string() const;

    fmpzxx numerator; /**< The numerator of the rational number */
    fmpzxx denominator; /**< The denominator of the rational number */
  };

  RationalNumber gcd(const RationalNumber& a, const RationalNumber& b);
  std::ostream& operator<< (std::ostream& out, const RationalNumber&);
}
