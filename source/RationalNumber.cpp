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

#include "firefly/RationalNumber.hpp"

namespace firefly {

  RationalNumber::RationalNumber(const fmpzxx& numerator_, const fmpzxx& denominator_) {
    fmpzxx gcd_(gcd(numerator_, denominator_));
    numerator = numerator_ / gcd_;
    denominator = denominator_ / gcd_;

    if (denominator < 0) {
      numerator = -numerator;
      denominator = -denominator;
    }
  }

  RationalNumber::RationalNumber() {}

  RationalNumber RationalNumber::operator*(const RationalNumber& rn) {
    fmpzxx num(numerator * rn.numerator);
    fmpzxx den(denominator * rn.denominator);

    if (den < 0) {
      num = -num;
      den = -den;
    }

    fmpzxx gcd_(gcd(num, den));
    numerator = num / gcd_;
    denominator = den / gcd_;
    return *this;
  }

  RationalNumber& RationalNumber::operator-=(const RationalNumber& rn) {
    if (rn.denominator != denominator) {
      numerator = numerator * rn.denominator - rn.numerator * denominator;
      denominator = denominator * rn.denominator;
    } else numerator -= rn.numerator;

    if (denominator < 0) {
      numerator = -numerator;
      denominator = -denominator;
    }

    fmpzxx gcd_(gcd(numerator, denominator));
    numerator = numerator / gcd_;
    denominator = denominator / gcd_;
    return *this;
  }

  RationalNumber& RationalNumber::operator+=(const RationalNumber& rn) {
    if (rn.denominator != denominator) {
      numerator = numerator * rn.denominator + rn.numerator * denominator;
      denominator = denominator * rn.denominator;
    } else numerator += rn.numerator;

    fmpzxx gcd_(gcd(numerator, denominator));
    numerator = numerator / gcd_;
    denominator = denominator / gcd_;
    return *this;
  }

  RationalNumber& RationalNumber::operator*=(const RationalNumber& rn) {
    fmpzxx num(numerator * rn.numerator);
    fmpzxx den(denominator * rn.denominator);

    if (den < 0) {
      num = -num;
      den = -den;
    }

    fmpzxx gcd_(gcd(num, den));
    numerator = num / gcd_;
    denominator = den / gcd_;
    return *this;
  }

  RationalNumber RationalNumber::operator-() const {
    return RationalNumber((-numerator).evaluate(), denominator);
  }

  bool RationalNumber::operator==(const RationalNumber& b) const {
    return (numerator == b.numerator && denominator == b.denominator);
  }

  std::string RationalNumber::string() const {
    std::string str;

    if (denominator == 1) {
      if (numerator < 1) {
        str += "(" + numerator.to_string() + ")";
      } else {
        str += numerator.to_string();
      }
    } else {
      if (numerator < 1) {
        str += "(" + numerator.to_string() + "/" + denominator.to_string() + ")";
      } else {
        str += numerator.to_string() + "/" + denominator.to_string();
      }
    }

    return str;
  }

  RationalNumber gcd(const RationalNumber& a, const RationalNumber& b) {
    fmpzxx numerator(gcd(a.numerator * b.denominator, b.numerator * a.denominator));
    fmpzxx denominator(a.denominator * b.denominator);
    return RationalNumber(numerator, denominator);
  }

  std::ostream& operator<< (std::ostream& out, const RationalNumber& a) {
    if (a.denominator == 1) {
      if (a.numerator < 1) {
        out << "(" << a.numerator.to_string() << ")";
      } else {
        out << a.numerator.to_string();
      }
    } else {
      out << "(" << a.numerator.to_string() << "/" << a.denominator.to_string() << ")";
    }

    return out;
  }
}
