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
#include "firefly/PolynomialFF.hpp"
#include "firefly/UintHasher.hpp"

#include <mutex>
#include <queue>

namespace firefly {
  typedef std::unordered_map<std::vector<uint32_t>, fmpzxx, UintHasher> mpz_map;
  typedef std::unordered_map<uint32_t, std::unordered_map<std::vector<uint32_t>, fmpzxx, UintHasher>> mpz_map_map;
  typedef std::unordered_map<std::pair<uint32_t, uint32_t>, FFInt, UintPairHasher> ff_pair_map;
  typedef std::unordered_map<std::vector<uint32_t>, RationalNumber, UintHasher> rn_map;
  typedef std::unordered_map<std::vector<uint32_t>, fmpqxx, UintHasher> fmpq_map;
  typedef std::unordered_map<std::vector<uint32_t>, std::unordered_map<std::vector<uint32_t>, FFInt, UintHasher>, UintHasher> ff_map_map;
  typedef std::unordered_map<std::vector<uint32_t>, std::queue<std::pair<FFInt, FFInt>>, UintHasher> ff_queue_map;
  typedef std::unordered_map<std::vector<uint32_t>, uint32_t, UintHasher> uint32_t_map;
  typedef std::unordered_map<uint32_t, std::vector<PolynomialFF>> polff_vec_map;
  typedef std::unordered_map<uint32_t, PolynomialFF> polff_map;

  /**
   * @class BaseReconst
   * @brief A base class for reconstruction objects
   */
  class BaseReconst {
  public:
    /**
     *  Constructor of the base reconstruction class
     */
    BaseReconst();
    BaseReconst(const BaseReconst& other);
    BaseReconst(BaseReconst && other);
    BaseReconst& operator=(const BaseReconst& other);
    BaseReconst& operator=(BaseReconst && other);
    /**
     *  Resets the states of the PRNG
     */
    static void reset();
    /**
     *  @return a 32-bit (64-bit) random number as an FFInt
     */
    FFInt get_rand_32();
    /**
     *  @return a 64-bit random number as an FFInt
     */
    FFInt get_rand_64();
    /**
     *  Calls the get_rand_64() function
     *  @deprecated since 1.3.2
     *  @return a 64-bit random number as an FFInt
     */
    FFInt get_rand();
    /**
     *  @return the number of equations needed for the current zi_order
     */
    uint32_t get_num_eqn() const;
    /**
     *  @return a bool if the current reconstruction is done
     */
    bool is_done() const;
    /**
     *  @return a bool if the current reconstruction needs a new prime
     */
    bool is_new_prime() const;
    /**
     *  @return the counter of the currently used prime
     */
    uint32_t get_prime() const;
    /**
     *  @returns a pair of a bool if the current reconstruction is done and the
     *  counter of the currently used prime
     */
    std::pair<bool, uint32_t> get_done_and_prime() const;
    /**
     *  @returns the currently used zi_order
     */
    std::vector<uint32_t> get_zi_order() const;
    /**
     *  @return returns the currently interpoleded i of the corresponding zi
     */
    uint32_t get_zi() const;
    /**
     *  Sets a new seed for the random number generator
     */
    void set_seed(uint64_t seed);
  protected:
    std::vector<uint32_t> curr_zi_order {}; /**< A vector which holds the current zi_order. It's length is n - 1 since the first variable is always set to 1. */
    fmpzxx combined_prime; /**< The combination of the used prime numbers with the chinese remained theorem */
    mutable std::mutex mutex_status; /**< A mutex to make all status variables thread safe */
    uint32_t prime_number = 0; /**< An integer which indicates the currently used prime counter */
    uint32_t num_eqn = 0; /**< An integer which indicates the currently needed number of equations to get to the next zi_order  */
    uint32_t n = 0; /**< The number of parameters */
    uint32_t type; /**< A flag which indicates if the object is used for polynomial or rational function reconstruction */
    uint32_t zi = 1; /**< An integer which indicates the current zi which is being interpolated */
    bool use_chinese_remainder = false; /**< A bool which indicates if the Chinese Remainder Theorem is used */
    bool check = false; /**< A bool which indicates if one could check the current interpolation */
    bool done = false; /**< A bool which indicates if the current reconstruction is done */
    bool new_prime = false; /**< A bool which indicates if the reconstruction needs a new prime */
    bool is_interpolating = false; /**< A bool which indicates if the object is currently interpolating the black box */
    /**
     *    Converts the coefficients of a rational function from FFInts to fmpzxx
     *    objects
     *    @param coefs a ff_map over a finite field
     *    @return the coefficients of the given rational function converted to
     *    fmpzxx objects
     */
    mpz_map convert_to_mpz(const ff_map& coefs) const;
    /**
     *    Converts the elements of a vector of RationalNumber objects to FFInts
     *    @param ri the vector of RationalNumber objects
     *    @return elements of ri converted to FFInts
     */
    ff_map convert_to_ffint(const rn_map& ri) const;
    enum reconst_type {POLY, RAT}; /**< An enum to hold the type flags */
  private:
    static std::mutex mutex_state; /**< A mutex to make all state variables thread safe */
    static uint64_t state; /**< An integer which is needed for the random number generator */
    static uint64_t const multiplier; /**< An integer which is needed for the random number generator */
    static uint64_t const increment; /**< An integer which is needed for the random number generator */
    static uint64_t splitmix64_state; /**< The sate of the splitmix64 PRNG */
    static uint64_t s[4]; /**< The state of the Xoshiro256** PRNG */
    /**
     *  Helper function needed by the pcg32 random number generator
     */
    uint32_t rotr32(uint32_t x, uint32_t r);
    /**
     *  @return a 32-bit random number with the pcg32 algorithm
     */
    uint32_t pcg32();
    /**
     *  Initializes the random number generator with a seed
     *  @param seed the seed
     */
    void pc32_init(uint64_t seed);
    /**
     *  Helper for the Xoshiro256**
     */
    uint64_t rol64(uint64_t x, int k);
    /**
     *  @return a 64 bit pseudo-random number
     */
    uint64_t xoshiro256ss();
    /**
     *  Initialize one seed for Xoshiro256**
     *  @return a 64-bit seed
     */
    uint64_t splitmix64();
    /**
     *  Set the seeds for the Xoshiro256** PRNG
     *  @param seed the seed
     */
    void xoshiro256ss_init(uint64_t seed);
  };
}
