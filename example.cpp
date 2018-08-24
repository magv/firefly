#include <iostream>
#include "PolyReconst.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"

int main() {
  uint64_t prime = firefly::primes()[0];
  //firefly::RatReconst rec_rat(1, prime);
  firefly::PolyReconst rec_pol(5, prime);

  try {
    int i = 1;
    std::vector<firefly::FFInt> yis;
    yis.emplace_back(firefly::FFInt(std::rand() % prime, prime));
    yis.emplace_back(firefly::FFInt(std::rand() % prime, prime));
    yis.emplace_back(firefly::FFInt(std::rand() % prime, prime));
    yis.emplace_back(firefly::FFInt(std::rand() % prime, prime));
    yis.emplace_back(firefly::FFInt(std::rand() % prime, prime));
    while (!rec_pol.done) {
      if(rec_pol.new_prime) {
        prime = firefly::primes()[i];
        i++;
      }
      yis[rec_pol.next_zi - 1] = firefly::FFInt(std::rand() % prime, prime);
      mpz_class test;
      test = "1234567891098987998798709805302432022989874343098";
      test = test % prime;
      firefly::FFInt a7(std::stoull(test.get_str()), prime);
      firefly::FFInt a1(1,prime);
      firefly::FFInt a2(18, prime);
      firefly::FFInt a3(25, prime);
      firefly::FFInt a4(300, prime);
      firefly::FFInt a5(2, prime);
      firefly::FFInt a6(7, prime);
      firefly::FFInt num = (a6-a7/a6*yis[0].pow(a5) + a1*yis[2] + a3*yis[0]*yis[1] + a4*yis[3]+ a2*yis[4]);
      rec_pol.feed(prime, yis, num);
    }

    std::cout << rec_pol.get_result();
    //auto rat_fun = rec_rat.reconst();
    //auto pol_fun = rec_pol.reconst();
    //std::cout << rat_fun;
    //std::cout << "f(x) = " << pol_fun;
  } catch (std::exception &e) {
    ERROR_MSG(e.what());
  }

  return 0;
}

