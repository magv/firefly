#pragma once

#include "FFInt.hpp"
#include "PolynomialFF.hpp"

namespace firefly {

  class ThieleInterpolator {
  public:
    /**
     *  Constructor
     */
    ThieleInterpolator();
    /**
     *  Adds a point to the Thiele interpolation
     *  @param num the numerical value of the black box
     *  @param yi the numerical value at which the black box has been probed
     *  @return a bool which indicates if the interpolation terminated
     */
    bool add_point(const FFInt& num, const FFInt& yi);
    /**
     *  @return returns a pair of ff_maps which correspond to the numerator and denominator
     */
    std::pair<ff_map, ff_map> get_result();
    ThieleInterpolator& operator=(const ThieleInterpolator&) = default;
  private:
    std::vector<FFInt> ai {}; /**< A vector which holds all coefficients a_i */
    std::vector<FFInt> ti {}; /**< A vector which holds all arguments t_i */
    /**
    *    Computes the coefficient a(i) = ai.at(i) recursively
    *    @param i The order of a(i)
    *    @param ip Recursion order
    *    @param num f(y_i)
    *    @return a(i)
    */
    FFInt comp_ai(int i, int ip, const FFInt& num);
    /**
    *    Constructs the canonical form of the rational function recursivly
    *    @return the rational function in its canonical form
    */
    std::pair<ff_map, ff_map>  construct_canonical();
    /**
     *    Iterates Thiele's interpolation formula to get the canonical form
     *    of the rational function
     *    @param i an integer telling the current degree of the rational function
     *    @return the recursivly iterated rational function in its canonical form
     */
    std::pair<PolynomialFF, PolynomialFF> iterate_canonical(uint32_t i);
    /**
     *    Calculates f(y_i) using  Thiele's interpolation formula
     *    @param i order of the highest coefficient a_i
     *    @param ip order of sub coefficient a_ip
     *    @param y y_i
     *    @returns f(y_i)
     */
    FFInt comp_fyi(uint32_t i, uint32_t ip, const FFInt& y);
  };
}
