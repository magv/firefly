#include "RatReconst.hpp"
#include "Logger.hpp"
#include "ReconstHelper.hpp"

namespace firefly {

RatReconst::RatReconst(int n_) : n(n_), prime(primes().at(0)) {
   ai.reserve(1000);
   yi.reserve(1000);
}

std::vector<FFInt> RatReconst::reconst(){
   uint maxDegree = 1000;
   const int breakCondition = 3;

   yi.emplace_back(FFInt(std::rand() % prime, prime));
   ai.emplace_back(num(prime, yi.back()));

   for(uint i = 1; i < maxDegree; i++){
      yi.emplace_back(FFInt(std::rand() % prime, prime));
      FFInt fyi;
      bool spuriousPole = true;
      while(spuriousPole){
	 try{
	    fyi = num(prime, yi.back());
	    spuriousPole = false;
	 } catch(const std::exception&){
	    yi.pop_back();
	    yi.emplace_back(FFInt(std::rand() % prime, prime));
	 }
      }

      if(fyi == compFyi(i - 1, yi.back())){
	 bool nonequal = false;
	 for(uint j = 0; j < breakCondition; j++){
	    const FFInt y = FFInt(std::rand() % prime, prime);
	    FFInt fy = num(prime, y);
	    if(fy != compFyi(i - 1, y)){
	       nonequal = true;
	       break;
	    }
	 }
	 if(!nonequal){
	    yi.pop_back();
	    break;
	 }
      }

      spuriousPole = true;
      while(spuriousPole){
	 try{
	    ai.emplace_back(compAi(i, i, fyi));
	    spuriousPole = false;
	 } catch (const std::exception&){
	    yi.pop_back();
	    yi.emplace_back(FFInt(std::rand() % prime, prime));
	 }
      }
      if(i == maxDegree - 1) {
	 maxDegree += 1000;
	 yi.resize(maxDegree);
	 ai.resize(maxDegree);
      }
   }

   return ai;
}

FFInt RatReconst::compAi(int i, int ip, const FFInt& num){
   if(ip == 0){
      return num;
   } else{
      if(compAi(i, ip - 1, num) == ai.at(ip - 1)) throw std::runtime_error("Divide by 0 error!");
      return (yi.at(i) - yi.at(ip - 1))/(compAi(i, ip - 1, num) - ai.at(ip - 1));
   }
}

FFInt RatReconst::compFyi(int i, const FFInt& y){
   return iterateCanonicalNum(i, i, y);
}

FFInt RatReconst::iterateCanonicalNum(uint i, uint ip, const FFInt& y){
   if(ip == 0){
      return ai.at(i);
   } else{
      return ai.at(i - ip) + (FFInt(0, prime) - yi.at(i - ip) + y)/iterateCanonicalNum(i, ip - 1, y);
   }
}

void RatReconst::constrCanonical(){
   if(ai.size() == 0){
      INFO_MSG("Rational function not yet reconstructed or 0.");
   } else if(ai.size() == 1){
      std::vector<FFInt> coefNom{ai.at(0)};
      std::vector<FFInt> coefDen{FFInt(1, prime)};
      Polynomial nom(coefNom);
      Polynomial den(coefDen);
      canonical = std::pair<Polynomial, Polynomial> (nom, den);
   } else{
      std::pair<Polynomial, Polynomial> r = iterateCanonical(1);
      std::vector<FFInt> a0{ai.at(0)};
      std::vector<FFInt> coefZ{FFInt(0, prime) - yi.at(0), FFInt(1, prime)};
      Polynomial constant(a0);
      Polynomial zPol(coefZ);
      canonical = std::pair<Polynomial, Polynomial> (constant*r.first
	 + zPol*r.second, r.first);
   }
}

std::pair<Polynomial, Polynomial> RatReconst::iterateCanonical(uint i){
   if(i < ai.size() - 1){
      std::pair<Polynomial, Polynomial> fnp1 = iterateCanonical(i + 1);
      Polynomial p1 (std::vector<FFInt> {ai.at(i)});
      Polynomial p2 (std::vector<FFInt> {FFInt(0, prime) - yi.at(i), FFInt(1, prime)});
      return std::pair<Polynomial, Polynomial> (fnp1.first*p1 + fnp1.second*p2,
						fnp1.first);
   } else{
      Polynomial p1 (std::vector<FFInt> {ai.at(i)});
      Polynomial p2 (std::vector<FFInt> {FFInt(1, prime)});
      return std::pair<Polynomial, Polynomial> (p1, p2);
   }
}


FFInt RatReconst::num(uint64_t p, const FFInt& y){
   FFInt a0 (2, p);
   FFInt a1 (6, p);
   FFInt a2 (2, p);
   FFInt a3 (7, p);
   FFInt a4 (30, p);
   FFInt a5 (2, p);
   FFInt a6 (7, p);
   FFInt a7 (100, p);
   FFInt a8 (13, p);
   FFInt exp2 (2, p);
   FFInt exp3 (3, p);
   FFInt exp4 (4, p);
   FFInt exp5 (5, p);
   FFInt exp6 (6, p);
   FFInt exp7 (7, p);
   FFInt exp8 (8, p);

   return (a0)/(a2 + a3*y + a4*y.pow(exp2));
}

}