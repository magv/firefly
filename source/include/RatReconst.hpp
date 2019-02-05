#pragma once

#include "PolyReconst.hpp"
#include "RationalFunction.hpp"

namespace firefly {

  class RatReconst : public BaseReconst {
  public:
    /**
     *    A constructor
     *    @param n_ the number of parameters
     */
    RatReconst(uint n_);
    RatReconst(const RatReconst& other);
    RatReconst(RatReconst && other);
    RatReconst& operator=(const RatReconst& other);
    RatReconst& operator=(RatReconst && other);
    void feed(const FFInt& new_ti, const FFInt& num, const std::vector<uint>& feed_zi_ord, const uint& fed_prime);
    RationalFunction get_result();
    void interpolate();
    void disable_shift();
    void generate_anchor_points();
    FFInt get_rand_zi(uint zi, uint order);
    std::vector<FFInt> get_rand_zi_vec(std::vector<uint> order);
    FFInt get_zi_shift(uint zi);
    std::vector<FFInt> get_zi_shift_vec();
    bool need_shift();
  private:
    void interpolate(const FFInt& new_ti, const FFInt& num, const std::vector<uint>& feed_zi_ord);
    FFInt comp_ai(int i, int ip, const FFInt& num);
    /**
     *    Normalize the rational function such that the first non-zero coefficient
     *    of the denominator is normalized to 1
     *    @param ratFun the rational function
     *    @param prime a prime number defining the finite field
     *    @return A normalized version of ratFun
     */
    RationalFunction normalize(RationalFunction& rf);
    /**
     *    Constructs the canonical form of the rational function recursivly
     *    @param ai a vector of the coefficients ai as FFints
     *    @param prime a prime number defining the current finite field
     *    @return the rational function in its canonical form
     */
    std::pair<ff_map, ff_map>  construct_canonical();
    /**
     *    Iterates Thiele's interpolation formula to get the canonical form
     *    of the rational function
     *    @param ai a vector of the coefficients ai as FFInts
     *    @param i an integer telling the current degree of the rational function
     *    @param prime a prime number defining the current finite field
     *    @return the recursivly iterated rational function in its canonical form
     */
    std::pair<PolynomialFF, PolynomialFF> iterate_canonical(uint i);
    /**
     *    Calculates f(y_i) using  Thiele's interpolation formula
     *    @param ai a vector of coefficients ai
     *    @param i order of the highest coefficient a_i
     *    @param ip order of sub coefficient a_ip
     *    @param y y_i
     *    @param prime a prime number defining the finite field
     *    @returns f(y_i)
     */
    FFInt comp_fyi(uint i, uint ip, const FFInt& y);
    /**
     *    Test if the guess yields the same answer for the function in the finite
     *    field of prime
     *    @param prime The prime number which defines the finite field
     *    @return true or false
     */
    bool test_guess(const FFInt& num);
    bool rec_rat_coef();
    std::pair<ff_map, ff_map> solve_gauss();
    std::pair<ff_map, ff_map> solve_homogenized_multi_gauss();
    std::tuple<int, uint, std::vector<uint>> feed_poly(int curr_deg,
                                                       uint max_deg, std::unordered_map<uint, PolyReconst>& coef,
                                                       PolyReconst& rec, ff_map_map& saved_num,
                                                       std::unordered_map<uint, PolynomialFF>& sub_save, bool is_num);
    void combine_primes(std::pair<mpz_map, mpz_map>& tmp);
    void build_uni_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, const std::vector<FFInt>& yis);
    void build_homogenized_multi_gauss(const FFInt& tmp_ti, const FFInt& tmp_num, const std::vector<FFInt>& yis);
    bool first_run = true;
    std::list<std::tuple<FFInt, FFInt, std::vector<uint>>> queue;
    std::vector<std::vector<FFInt>> coef_mat {};
    std::unordered_map<uint, std::vector<FFInt>> coef_mat_num {}; //new
    std::unordered_map<uint, std::vector<FFInt>> coef_mat_den {}; //new
    PolynomialFF solved_num; //new
    PolynomialFF solved_den; //new
    uint curr_zi = 2;
    ff_vec_map saved_ti {};
    std::vector<FFInt> ai {};
    std::unordered_map<uint, PolyReconst> coef_n {};
    std::unordered_map<uint, PolyReconst> coef_d {};
    std::unordered_map<uint, PolynomialFF> sub_num {};
    std::unordered_map<uint, PolynomialFF> sub_den {};
    std::unordered_map<uint, std::vector<std::vector<uint>>> non_solved_degs_num {};// a vector entry should be just a pointer to save memory
    std::unordered_map<uint, std::vector<std::vector<uint>>> non_solved_degs_den {};
    ff_map_map saved_num_num {};
    ff_map_map saved_num_den {};
    int max_deg_num = -1;
    int max_deg_den = -1;
    int curr_deg_num = -1;
    int curr_deg_den = -1;
    bool is_singular_system = false;
    std::vector<uint> curr_zi_order_num {};
    std::vector<uint> curr_zi_order_den {};
    uint tmp_solved_coefs_num = 0;
    uint tmp_solved_coefs_den = 0;
    void remove_ni(const std::vector<uint>& deg_vec, RationalNumber& rn);
    void remove_di(const std::vector<uint>& deg_vec, RationalNumber& rn);
    RationalFunction result;
    std::vector<FFInt> ti {}; /**< A vector which holds all arguments t_i */
    rn_map g_ni {}; /**< rational coefficient guesses for the numerator*/
    rn_map g_di {}; /**< rational coefficient guesses for the denominator*/
    mpz_map combined_ni {};  /**< The combination of the coefficients of the numerator over finite field with the chinese remained theorem */
    mpz_map combined_di {};  /**< The combination of the coefficients of the denominator over finite field with the chinese remained theorem */
    static std::mutex mutex_statics;
    void add_non_solved_num(const std::vector<uint>& deg);
    void add_non_solved_den(const std::vector<uint>& deg);
    void check_for_solved_degs(std::vector<uint>& uni_degs, const bool is_num);
    PolynomialFF solve_transposed_vandermonde(std::vector<std::vector<uint>>& degs,
                                              const std::vector<FFInt>& nums);
    void set_new_curr_deg_num_singular(uint key);
    void set_new_curr_deg_den_singular(uint key);
    std::unordered_map<uint, PolynomialFF> solved_degs_num {}; //new
    std::unordered_map<uint, PolynomialFF> solved_degs_den {}; //new
    std::vector<uint> min_deg_den_vec {};//new
    static std::vector<FFInt> shift;
    static ff_pair_map rand_zi;
    static bool need_prime_shift;
  };
}
