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
//====================================================================================

#pragma once

#include "firefly/config.hpp"

#include <flint/ulong_extras.h>
#include <flint/fmpzxx.h>
#include <flint/fmpqxx.h>
#include <iostream>
#include <string>
#include <vector>

namespace firefly {

  typedef flint::fmpzxx fmpzxx;
  typedef flint::fmpqxx fmpqxx;

  /**
  * @class FFInt
  * @brief A class for finite field integers
  */
  class FFInt {
  public:
    /**
     *    A constructor
     *    @param n_ an integer which is a member of the finite field
     */
    template<class T, typename=typename std::enable_if<(std::is_enum<T>::value || std::is_integral<T>::value)>::type>
    FFInt(const T n_) {
      if (n_ >= 0) {
        if (static_cast<uint64_t>(n_) < mod.n) {
          n = static_cast<uint64_t>(n_);
        } else {
          n = static_cast<uint64_t>(n_)  % mod.n;
        }
      } else if (n_ < 0) {
        n = mod.n - static_cast<uint64_t>(-n_) % mod.n;
      }
    }
    /**
     *    A constructor
     *    @param ffint a FFInt object
     */
    FFInt(const FFInt& ffint) : n(ffint.n) {};
    /**
     *  A constructor for a flint::fmpzxx object
     *  @param in the flint::fmpzxx object which should be converted to an FFInt
     */
    FFInt(const fmpzxx &in) : n(fmpz_get_nmod(in._fmpz(), FFInt::mod)) {};
    /**
     *    Default constructor
     */
    FFInt() : n(0) { };
    /**
     *  A direct constructor that bypasses the modulo operation.
     *  @param n an integer, 0 <= n < mod.n
     */
    static inline FFInt from(const uint64_t n) { FFInt x; x.n = n; return x; }
    /**
     *  A static function to define a new prime field
     *  @param prime the defining prime of the field
     */
    static inline void set_new_prime(uint64_t prime) {
      nmod_init(&FFInt::mod, prime);
      FFInt::p = mod.n;
      FFInt::p_inv = mod.ninv;
    }
    /**
     *  Converts the uint64_t to a negative integer (only used for negative exponents)
     */
    int to_neg_int() const { return -static_cast<int>(n); }

    // defining new operators for finite field arithmetic
    FFInt &operator=(const FFInt&) = default;
    FFInt &operator+=(const FFInt &x) { n = _nmod_add(n, x.n, mod); return *this; }
    FFInt &operator-=(const FFInt &x) { n = _nmod_sub(n, x.n, mod); return *this; };
    FFInt &operator*=(const FFInt &x) { n = nmod_mul(n, x.n, mod); return *this; };
    FFInt &operator/=(const FFInt &x) { n = nmod_mul(n, nmod_inv(x.n, mod), mod); return *this; };
    FFInt operator-() const { return FFInt::from(nmod_neg(n, mod)); }
    FFInt operator+() const { return *this; };
    FFInt operator++() { n = _nmod_add(n, 1, mod); return *this; };
    FFInt operator++(int) { FFInt tmp = *this; n = _nmod_add(n, 1, mod); return tmp; };
    FFInt operator--() { n = _nmod_sub(n, 1, mod); return *this; };
    FFInt operator--(int) { FFInt tmp = *this; n = _nmod_sub(n, 1, mod); return tmp; };
    bool operator!() const { return !n; };
    FFInt pow(const FFInt &ffint) const { return FFInt::from(nmod_pow_ui(n, ffint.n, mod)); }
    template<class T, typename=typename std::enable_if<(std::is_enum<T>::value || std::is_integral<T>::value)>::type>
    FFInt pow(const T &power) const {
        return FFInt::from(n_powmod2_preinv(n, power, mod.n, mod.ninv));
      }
    FFInt invert() const { return FFInt::from(nmod_inv(n, mod)); };

    uint64_t n; /**< the integer member of the finite field */
    static uint64_t p; /**< the prime defining the finite field */
    static uint64_t p_inv; /**< the inverse of the defining prime needed for FFTs*/
    static nmod_t mod;
  };

  inline bool operator<(const FFInt &a, const FFInt &b) { return a.n < b.n; };
  inline bool operator<=(const FFInt &a, const FFInt &b) { return a.n <= b.n; };
  inline bool operator>(const FFInt& a, const FFInt &b) { return a.n > b.n; };
  inline bool operator>=(const FFInt &a, const FFInt &b) { return a.n >= b.n; };
  inline bool operator==(const FFInt &a, const FFInt &b) { return a.n == b.n; };
  inline bool operator!=(const FFInt &a, const FFInt &b) { return a.n != b.n; };
  inline FFInt operator/(const FFInt &a, const FFInt &b) { return FFInt::from(nmod_mul(a.n, nmod_inv(b.n, FFInt::mod), FFInt::mod)); };
  inline FFInt operator+(const FFInt &a, const FFInt &b) { return FFInt::from(_nmod_add(a.n, b.n, FFInt::mod)); };
  inline FFInt operator-(const FFInt &a, const FFInt &b) { return FFInt::from(_nmod_sub(a.n, b.n, FFInt::mod)); };
  inline FFInt operator*(const FFInt &a, const FFInt &b) { return FFInt::from(nmod_mul(a.n, b.n, FFInt::mod)); };
  inline FFInt pow(const FFInt& ffint, const FFInt& power) { return ffint.pow(power); };
  template<class T, typename=typename std::enable_if<(std::is_enum<T>::value || std::is_integral<T>::value)>::type>
  inline FFInt pow(const FFInt& ffint, const T& power) { return ffint.pow(power); };
  inline std::ostream &operator<<(std::ostream& out, const FFInt &x) { out << x.n; return out; }

  extern "C" {
    void firefly_exists(void);
  }
}
