#include <iostream>
#include "PolyReconst.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"

int main() {
  uint64_t prime = firefly::primes()[0];
  firefly::FFInt::p = prime;
  //firefly::RatReconst rec(1);
  firefly::PolyReconst rec(1);

  try {
    int i = 1;
    std::vector<firefly::FFInt> yis;
    yis.emplace_back(firefly::FFInt(std::rand() % prime));
    //yis.emplace_back(firefly::FFInt(std::rand() % prime));
    //yis.emplace_back(firefly::FFInt(std::rand() % prime));
    //yis.emplace_back(firefly::FFInt(std::rand() % prime, prime));
    //yis.emplace_back(firefly::FFInt(std::rand() % prime, prime));
    while (!rec.done) {
      if(rec.new_prime) {
        prime = firefly::primes()[i];
        firefly::FFInt::p = prime;

        /*for(int i = 0; i < yis.size(); i++){
          yis[i] = firefly::FFInt(std::rand() % prime);
        }*/
        i++;
      }
      yis[0] = firefly::FFInt(std::rand() % prime);
      //yis[rec_pol.next_zi - 1] = firefly::FFInt(std::rand() % prime);
      mpz_class test;
      test = "1234567891098987998798709805302432022989874343098";
      test = test % prime;
      firefly::FFInt a7(std::stoull(test.get_str()));
      firefly::FFInt a1(1);
      firefly::FFInt a2(18);
      firefly::FFInt a3(25);
      firefly::FFInt a4(100);
      firefly::FFInt a5(2);
      firefly::FFInt a6(7);
      firefly::FFInt num = (a6-a4/a6*yis[0].pow(a4));//(a1 + a1*yis[0]);//+ a1*yis[2] + a3*yis[0]*yis[1].pow(a6));
      rec.feed(yis, num);
    }

    std::cout << rec.get_result();
    //auto rat_fun = rec_rat.reconst();
    //auto pol_fun = rec_pol.reconst();
    //std::cout << rat_fun;
    //std::cout << "f(x) = " << pol_fun;
  } catch (std::exception &e) {
    ERROR_MSG(e.what());
  }

  return 0;
}