#include "Polynomial.hpp"
#include <chrono>

namespace firefly {

  Polynomial::Polynomial(const rn_map& coef) {
    for(const auto& el : coef){
      coefs.emplace_back(Monomial(el.first, el.second));
    }
    n = coefs[0].powers.size();
  }

  Polynomial::Polynomial(const Monomial& coef) {
    coefs.emplace_back(coef);
    n = coef.powers.size();
  }

  Polynomial::Polynomial() {}

  void Polynomial::sort() {
    std::sort(coefs.begin(), coefs.end());
  }

  void Polynomial::clear() {
    coefs.clear();
  }

  std::string Polynomial::string(const std::vector<std::string>& symbols) const {
    std::string str;

    for (const auto & mono : coefs) {
      str += mono.coef.string() + "*";

      for (uint i = 0; i < mono.powers.size(); i++) {
        if (mono.powers[i] > 1) {
          str += symbols[i] + "^" + std::to_string(mono.powers[i]) + "*";
        } else if (mono.powers[i] == 1) {
          str += symbols[i] + "*";
        }
      }

      str.erase(--str.end());
      str += "+";
    }

    str.erase(--str.end());
    return str;
  }

  std::ostream& operator<<(std::ostream& out, const Polynomial& pol) {
    bool first = true;

    for (const auto & mono : pol.coefs) {
      if (first) {
        out <<  mono.coef << "*x^(";

        for (const auto i : mono.powers) {
          out << i << ",";
        }

        out << "\b)";
        first = false;
      } else {
        out << " + " << mono.coef << "*x^(";

        for (const auto i : mono.powers) {
          out << i << ",";
        }

        out << "\b)";
      }

    }

    out << "\n";
    return out;
  }

  PolynomialFF Polynomial::convert_to_PolynomialFF() {
    ff_map coefs_ff;
    uint n = coefs[0].powers.size();

    for (auto & coef : coefs) {
      mpz_class numerator = coef.coef.numerator % FFInt::p;

      if (numerator < 0) numerator += FFInt::p;

      mpz_class denominator = coef.coef.denominator % FFInt::p;
      FFInt coef_ff = FFInt(std::stoull(numerator.get_str())) / FFInt(std::stoull(denominator.get_str()));

      if (coef_ff.n > 0)
        coefs_ff.emplace(std::make_pair(coef.powers, coef_ff));
    }

    return PolynomialFF(n, coefs_ff);
  }

  Polynomial Polynomial::operator*(const RationalNumber& rn) {
    for (auto & mon : coefs) {
      mon.coef *= rn;
    }

    return *this;
  }

}

