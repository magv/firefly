#include "RatReconst.hpp"
#include "Logger.hpp"
#include "ReconstHelper.hpp"
#include "utils.hpp"
#include <chrono>
#include <algorithm>

namespace firefly {
  std::vector<FFInt> RatReconst::shift {};
  bool RatReconst::shifted = false;
  ff_pair_map RatReconst::rand_zi;
  std::vector<FFInt> RatReconst::anchor_points {};

  RatReconst::RatReconst(uint n_) : n(n_) {
    std::unique_lock<std::mutex> lock_status(mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_feed(mutex_feed, std::defer_lock);
    std::lock(lock_status, lock_feed);

    ti.reserve(300);
    ai.reserve(300);
    combined_prime = FFInt::p;
//     std::srand(std::time(0));

    if (!shifted) {
      shift = std::vector<FFInt> (n);

      if (n > 1) {
        for (auto & el : shift) el = FFInt(std::rand() % 1000000) + FFInt(1);

        curr_zi_order_num = std::vector<uint> (n, 1);
        curr_zi_order_den = std::vector<uint> (n, 1);
        shifted = true;
      }
    }

    if (n > 1) {
      deg_num.emplace_back(-1);
      deg_den.emplace_back(-1);
      curr_zi_order = std::vector<uint> (n, 1);
      curr_zi_order[n - 1] = 0;

      // fill in the rand_vars for zi_order = 1
      if (rand_zi.empty()) {
        for (uint i = 2; i <= n; i++) {
          const FFInt rand = get_rand();
          rand_zi.emplace(std::make_pair(std::make_pair(i, 1), rand));
          anchor_points.emplace_back(rand);
        }
      }
    }
  }

  void RatReconst::feed(const FFInt& new_ti, const FFInt& num, const std::vector<uint>& feed_zi_ord, const uint& fed_prime) {
    std::unique_lock<std::mutex> lock(mutex_feed);
    feed(new_ti, num, feed_zi_ord, fed_prime, lock);
  }

  void RatReconst::feed(const FFInt& new_ti, const FFInt& num, const std::vector<uint>& feed_zi_ord, const uint& fed_prime, std::unique_lock<std::mutex>& lock) {
    if (!done && fed_prime == prime_number) {
      std::vector<uint> tmp_vec;
      std::vector<uint> tmp_vec_rev;
      std::vector<uint> feed_zi_ord_rev;

      if (n > 1) {
        tmp_vec = std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end() - 1);
        tmp_vec_rev = tmp_vec;
        feed_zi_ord_rev = feed_zi_ord;
        std::reverse(feed_zi_ord_rev.begin(), feed_zi_ord_rev.end());
        std::reverse(tmp_vec_rev.begin(), tmp_vec_rev.end());
      }

      // Compare if the food is the expected food; if not, store it for later use
      if (feed_zi_ord == tmp_vec) {
        // first check if we are done. If not start the reconstruction again using
        // the chinese remainder theorem in combining the previous results
        if (new_prime) {
          ti.emplace_back(new_ti);
          sub_num.clear();
          sub_den.clear();

          if (rec_rat_coef()) {
            {
              std::unique_lock<std::mutex> lock(mutex_status);
              done = test_guess(num);
            }

            if (done) {
              std::unique_lock<std::mutex> lock(mutex_status);
              coef_n.clear();
              coef_d.clear();
              combined_di.clear();
              combined_ni.clear();
              degs_n.clear();
              degs_d.clear();
              combined_prime = 0;
              new_prime = false;
              deg_num.clear();
              deg_den.clear();
              curr_zi_order.clear();
              saved_num_num.clear();
              saved_num_den.clear();
              saved_num_den.clear();
              saved_num_num.clear();
              non_solved_coef_num.clear();
              non_solved_coef_den.clear();
              use_chinese_remainder = false;
              return;
            }
          }

          if (!use_chinese_remainder) use_chinese_remainder = true;

          new_prime = false;
          ti.pop_back();
        }

        // basic reconstruction algorithm, check if reconstructed function is equal
        // to numeric input and calculate coefficients a_i, check chinese chinese remainder
        // theorem

        {
          std::unique_lock<std::mutex> lock(mutex_status);
          zi = 1;
        }

        if (max_deg_num == -1) {
          ti.emplace_back(new_ti);
          const uint i = ti.size() - 1;

          if (i == 0) {
            ai.emplace_back(num);
          } else {
            if (num == comp_fyi(i - 1, i - 1, ti.back())) check = true;

            ai.emplace_back(comp_ai(i, i, num));
          }
        } else {
          uint size = coef_mat.size();

          if (size == 0)
            coef_mat.reserve(num_eqn);

          std::vector<std::pair<FFInt, FFInt>> t_food = {std::make_pair(new_ti, num)};

          if (n > 1) {
            try {
              std::vector<uint> tmp_vec = std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end() - 1);
              t_food.insert(t_food.end(), saved_ti.at(tmp_vec).begin(), saved_ti.at(tmp_vec).end());
              saved_ti.erase(tmp_vec);
            } catch (std::out_of_range& e) {
              // do nothing
            }
          }

          for (auto food : t_food) {
            std::vector<FFInt> solved_coef_sub_num {};
            std::vector<FFInt> solved_coef_sub_den {};

            // fill matrix
            std::vector<FFInt> eq;
            eq.reserve(num_eqn + 1);
            FFInt tmp_ti = food.first;
            FFInt tmp_num = food.second;

            std::vector<FFInt> yis;

            if (n > 1) {
              for (uint i = 0; i < curr_zi_order.size() - 1; i ++) {
                yis.emplace_back(rand_zi[std::make_pair(i + 2, curr_zi_order[i])]);
              }
            }

            std::vector<FFInt> ns_yis = yis;
            std::vector<FFInt> ns_yis_wo_t = yis;
            ns_yis.insert(ns_yis.begin(), FFInt(1));

            yis.insert(yis.begin(), tmp_ti);

            for (uint i = 1; i < n; i++) {
              yis[i] = yis[0] * yis[i] + shift[i];
            }

            yis[0] += shift[0];

            for (int r = 0; r <= max_deg_num; r++) {
              if (shift[0] != 0 && coef_n.size() > 0 && coef_n[r].new_prime) {
                FFInt sub;

                if (r < max_deg_num)
                  sub = (coef_n[r].get_result().convert_to_PolynomialFF().calc(ns_yis_wo_t) + sub_num[r].convert_to_PolynomialFF().calc(ns_yis)) * tmp_ti.pow(FFInt(r));
                else
                  sub = (coef_n[r].get_result().convert_to_PolynomialFF().calc(ns_yis_wo_t)) * tmp_ti.pow(FFInt(r));

                solved_coef_sub_num.emplace_back(sub);
              } else if (std::find(non_solved_coef_num.begin(), non_solved_coef_num.end(), r) != non_solved_coef_num.end()) {
                eq.emplace_back(tmp_ti.pow(FFInt(r)));
                solved_coef_sub_num.emplace_back(solved_coefs_num[r].convert_to_PolynomialFF().calc(yis));
              } else {
                FFInt sub = solved_coefs_num[r].convert_to_PolynomialFF().calc(yis);

                if (sub.n > 0)
                  solved_coef_sub_num.emplace_back(sub);
              }
            }

            for (int rp = min_deg_den + 1; rp <= max_deg_den; rp++) {
              if (shift[0] != 0  && coef_d.size() > 0 && coef_d[rp].new_prime) {
                FFInt sub;

                if (rp < max_deg_den)
                  sub = (coef_d[rp].get_result().convert_to_PolynomialFF().calc(ns_yis_wo_t) + sub_den[rp].convert_to_PolynomialFF().calc(ns_yis)) * tmp_ti.pow(FFInt(rp));
                else
                  sub = (coef_d[rp].get_result().convert_to_PolynomialFF().calc(ns_yis_wo_t)) * tmp_ti.pow(FFInt(rp));

                solved_coef_sub_den.emplace_back(sub);
              } else if (std::find(non_solved_coef_den.begin(), non_solved_coef_den.end(), rp) != non_solved_coef_den.end()) {
                eq.emplace_back(-tmp_ti.pow(FFInt(rp)) * tmp_num);
                solved_coef_sub_den.emplace_back(solved_coefs_den[rp - (min_deg_den + 1)].convert_to_PolynomialFF().calc(yis));
              } else {
                FFInt sub = solved_coefs_den[rp - (min_deg_den + 1)].convert_to_PolynomialFF().calc(yis);

                if (sub.n > 0)
                  solved_coef_sub_den.emplace_back(solved_coefs_den[rp - (min_deg_den + 1)].convert_to_PolynomialFF().calc(yis));
              }
            }

            eq.emplace_back(tmp_ti.pow(FFInt(min_deg_den)) * tmp_num);

            for (auto & solved_coef_num : solved_coef_sub_num) {
              eq.back() += -solved_coef_num;
            }

            for (auto & solved_coef_den : solved_coef_sub_den) {
              eq.back() += solved_coef_den * tmp_num;
            }

            coef_mat.emplace_back(std::move(eq));

            if (coef_mat.size() == num_eqn) {
              check = true;
              break;
            }
          }
        }

        if (check) {
          check = false;

          // todo not needed anymore. Only if one wants to check twice
          // if (num == comp_fyi(i - 1, i - 1, ti.back())) {

          std::pair<PolynomialFF, PolynomialFF> canonical;

          if (max_deg_num == -1) {
            if (ai.capacity() != ai.size()) {
              ai.shrink_to_fit();
              ti.shrink_to_fit();
            }

            ti.pop_back();
            ai.pop_back();

            canonical = construct_canonical();
            PolynomialFF denominator = canonical.second;

            //TODO catch new shift
            if (n > 1 && denominator.min_deg()[0] > 0) {
              INFO_MSG("No constant term in denominator! Trying again with new paramter shift...");

              for (uint j = 0; j < n; j++) {
                shift[j] = FFInt(std::rand() % 1000000) + FFInt(1);
              }

              {
                std::unique_lock<std::mutex> lock(mutex_status);
                done = false;
              }
              ai.clear();
              ti.clear();
              return;
            }

            max_deg_num = canonical.first.max_deg()[0];
            max_deg_den = canonical.second.max_deg()[0];
            curr_deg_num = max_deg_num;
            curr_deg_den = max_deg_den;
            min_deg_den = canonical.second.min_deg()[0];
            min_deg_num = canonical.first.min_deg()[0];

            non_solved_coef_num = std::vector<uint> (max_deg_num + 1);
            non_solved_coef_den = std::vector<uint> (max_deg_den + 1);

            FFInt equializer = FFInt(1) / denominator.coef[denominator.min_deg()];

            canonical.first = canonical.first * equializer;
            canonical.second = denominator * equializer;
            std::vector<uint> zero_deg(n);
            Monomial zero_mon(zero_deg, RationalNumber(0, 1));
            solved_coefs_num = std::vector<Polynomial> (max_deg_num + 1, Polynomial(zero_mon));
            solved_coefs_den = std::vector<Polynomial> (max_deg_den - min_deg_den, Polynomial(zero_mon));

            PolynomialFF numerator = canonical.first;
            uint deleted_coefs = 0;
            uint solved_coef_num = 0;
            uint solved_coef_den = 0;

            // check for coefficients which are zero and remove them to save numerical runs
            for (int i = 0; i <= max_deg_num; i++) {
              try {
                std::vector<uint> pow = {(uint) i};
                numerator.coef.at(pow);
                non_solved_coef_num[i - deleted_coefs] = i;
              } catch (std::out_of_range& e) {
                if (shift[0] == 0) {
                  std::vector<uint> pow(n, 0);
                  pow[0] = i;
                  solved_coefs_num[i] = Polynomial(Monomial(pow, RationalNumber(0, 1)));
                  non_solved_coef_num.erase(non_solved_coef_num.begin() + i - deleted_coefs);
                  deleted_coefs ++;
                  solved_coef_num ++;
                }
              }
            }

            deleted_coefs = 0;

            for (int i = min_deg_den + 1; i <= max_deg_den; i++) {
              try {
                std::vector<uint> pow = {(uint) i};
                denominator.coef.at(pow);
                non_solved_coef_den[i - deleted_coefs - 1 - min_deg_den] = i;
              } catch (std::out_of_range& e) {
                std::vector<uint> pow(n, 0);

                if (shift[0] == 0) {
                  pow[0] = i;
                  solved_coefs_den[i - 1 - min_deg_den] = Polynomial(Monomial(pow, RationalNumber(0, 1)));
                  non_solved_coef_den.erase(non_solved_coef_den.begin() + i - 1 - deleted_coefs - min_deg_den);
                  deleted_coefs ++;
                  solved_coef_den ++;
                }
              }

              solved_coefs = solved_coef_num + solved_coef_den;
              num_eqn = max_deg_den + max_deg_num + 1 - min_deg_num - min_deg_den
                        - solved_coefs - tmp_solved_coefs_num - tmp_solved_coefs_den;
            }

            ai.clear();
            ti.clear();
          } else
            canonical = solve_gauss();

          if (n == 1) {
            std::pair<mpz_map, mpz_map> tmp = convert_to_mpz(canonical);
            combine_primes(tmp);
            {
              std::unique_lock<std::mutex> lock(mutex_status);
              prime_number++;
            }
            saved_ti.clear();
            new_prime = true;
            return;
          } else {
            {
              std::unique_lock<std::mutex> lock(mutex_status);
              zi = curr_zi;
            }

            ff_map num_coef = canonical.first.coef;
            ff_map den_coef = canonical.second.coef;

            // save the current results to the map to access them later
            for (uint i = 0; i < non_solved_coef_num.size() - tmp_solved_coefs_num; i ++) {
              const uint deg = non_solved_coef_num[i];

              if (first_run) {
                PolyReconst rec(n - 1, anchor_points, deg);
                coef_n.emplace(std::make_pair(deg, std::move(rec)));
                deg_num.emplace_back(deg);

                if ((int) deg < max_deg_num) {
                  std::vector<uint> zero_deg(n);
                  Monomial zero_mon(zero_deg, RationalNumber(0, 1));
                  sub_num.emplace(std::make_pair(deg, Polynomial(zero_mon)));
                }
              }

              if ((int) deg <= curr_deg_num) {
                // this saves some memory since we only need one numerical value
                // for the constant coefficient
                if (deg == 0 && first_run) {
                  std::vector<uint> key = {deg, zi};
                  saved_num_num[curr_zi_order][key] = num_coef[ {deg}];
                } else {
                  std::vector<uint> key = {deg, zi};
                  saved_num_num[curr_zi_order][key] = num_coef[ {deg}];
                }
              }
            }

            for (uint i = 0; i < non_solved_coef_den.size() - tmp_solved_coefs_den; i ++) {
              const uint deg = non_solved_coef_den[i];

              if (first_run) {
                PolyReconst rec(n - 1, anchor_points, deg);
                coef_d.emplace(std::make_pair(deg, std::move(rec)));
                deg_den.emplace_back(deg);

                if ((int) deg < max_deg_den) {
                  std::vector<uint> zero_deg(n);
                  Monomial zero_mon(zero_deg, RationalNumber(0, 1));
                  sub_den.emplace(std::make_pair(deg, Polynomial(zero_mon)));
                }

              }

              if ((int) deg <= curr_deg_den) {
                // this saves some memory since we only need one numerical value
                // for the constant coefficient
                if (deg == 0 && first_run) {
                  std::vector<uint> key = {deg, zi};
                  saved_num_den[curr_zi_order][key] = den_coef[ {deg}];
                } else {
                  std::vector<uint> key = {deg, zi};
                  saved_num_den[curr_zi_order][key] = den_coef[ {deg}];
                }
              }
            }

            if (first_run) {
              std::sort(deg_num.begin(), deg_num.end());
              std::sort(deg_den.begin(), deg_den.end());
              first_run = false;
            }

            uint zi_num = 0;

            if (curr_deg_num > 0) zi_num = coef_n[curr_deg_num].next_zi + 1;

            uint zi_den = 0;

            if (curr_deg_den > 0) zi_den = coef_d[curr_deg_den].next_zi + 1;

            // reconstruct the numerator
            if (curr_deg_num >= 0) {
              PolyReconst rec = coef_n[curr_deg_num];

              if (rec.curr_zi_order == std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end() - 1)) {
                auto res = feed_poly(curr_deg_num, max_deg_num, coef_n, rec,
                                     saved_num_num, sub_num, true);
                curr_deg_num = std::get<0>(res);
                zi_num = std::get<1>(res);
                curr_zi_order_num = std::get<2>(res);
              }
            }

            // reconstruct the denominator
            if (curr_deg_den >= 0) {
              PolyReconst rec = coef_d[curr_deg_den];

              // if the numerator is done, get the current zi order of the
              // denominator
              if (curr_deg_num == -1) {
                std::unique_lock<std::mutex> lock(mutex_status);
                curr_zi_order = rec.curr_zi_order;
                curr_zi_order.emplace_back(prime_number);
                curr_zi = rec.next_zi + 1;
                zi = curr_zi;
              }

              if (rec.curr_zi_order == std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end() - 1)) {
                auto res = feed_poly(curr_deg_den, max_deg_den, coef_d, rec,
                                     saved_num_den, sub_den, false);
                curr_deg_den = std::get<0>(res);
                zi_den = std::get<1>(res);
                curr_zi_order_den = std::get<2>(res);

                // if the denominator is done, check if the numerator is still undone
                if (curr_deg_den == -1 && curr_deg_num >= 0) {
                  PolyReconst rec_new = coef_n[curr_deg_num];
                  {
                    std::unique_lock<std::mutex> lock(mutex_status);
                    curr_zi_order = rec_new.curr_zi_order;
                    curr_zi_order.emplace_back(prime_number);
                    curr_zi = rec_new.next_zi + 1;
                    zi = curr_zi;
                  }

                  if (rec_new.curr_zi_order == std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end() - 1)) {
                    auto res_new = feed_poly(curr_deg_num, max_deg_num, coef_n, rec_new,
                                             saved_num_num, sub_num, true);
                    curr_deg_num = std::get<0>(res_new);
                    zi_num = std::get<1>(res_new);
                    curr_zi_order_num = std::get<2>(res_new);
                  }
                }
              }
            }

            // check which poly reconst should be feeded next
            // it is promising that feeding a PolyReconst with a higher
            // zi degree will be finished next leading to less numerical runs
            if (curr_deg_den >= 0 && curr_deg_num >= 0) {
              std::vector<uint> zi_order_num_rev = curr_zi_order_num;
              std::reverse(zi_order_num_rev.begin(), zi_order_num_rev.end());
              std::vector<uint> zi_order_den_rev = curr_zi_order_den;
              std::reverse(zi_order_den_rev.begin(), zi_order_den_rev.end());

              if (zi_order_num_rev > zi_order_den_rev) {
                std::unique_lock<std::mutex> lock(mutex_status);
                curr_zi_order = curr_zi_order_num;
                curr_zi = zi_num;
                zi = zi_num;
              } else {
                std::unique_lock<std::mutex> lock(mutex_status);
                curr_zi_order = curr_zi_order_den;
                curr_zi = zi_den;
                zi = zi_den;
              }
            } else if (curr_deg_num >= 0) {
              std::unique_lock<std::mutex> lock(mutex_status);
              curr_zi_order = curr_zi_order_num;
              curr_zi = zi_num;
              zi = zi_num;
            } else if (curr_deg_den >= 0) {
              std::unique_lock<std::mutex> lock(mutex_status);
              curr_zi_order = curr_zi_order_den;
              curr_zi = zi_den;
              zi = zi_den;
            }

            // combine results
            if (curr_deg_den == - 1 && curr_deg_num == -1) {
              saved_num_num.clear();
              saved_num_den.clear();

              first_run = true;

              // Remove normalization due to the shift
              PolynomialFF numerator;
              PolynomialFF denominator;

              for (auto & el : coef_n) {
                std::vector<std::vector<uint>> degs;
                Polynomial result = el.second.get_result();
                for(const auto & mon : result.coefs){
                  degs.emplace_back(mon.powers);
                }
                degs_n[el.first] = degs;
                numerator = numerator + result.homogenize(el.first).convert_to_PolynomialFF();
              }

              for (auto & el : coef_d) {
                std::vector<std::vector<uint>> degs;
                Polynomial result = el.second.get_result();
                for(const auto & mon : result.coefs){
                  degs.emplace_back(mon.powers);
                }
                degs_d[el.first] = degs;
                denominator = denominator + result.homogenize(el.first).convert_to_PolynomialFF();
              }

              coef_n.clear();
              coef_d.clear();

              FFInt first_coef = denominator.coef[denominator.min_deg()];

              // normalize
              FFInt equializer = FFInt(1) / first_coef;

              numerator = numerator * equializer;
              denominator = denominator * equializer;

              std::pair<mpz_map, mpz_map> tmp = convert_to_mpz(std::make_pair(numerator, denominator));

              combine_primes(tmp);

              {
                std::unique_lock<std::mutex> lock(mutex_status);
                prime_number++;
              }
              saved_ti.clear();
              std::unique_lock<std::mutex> lock(mutex_status);
              std::fill(curr_zi_order.begin(), curr_zi_order.end() - 1, 1);
              curr_zi_order[n - 1] = prime_number;
              curr_zi = 2;
              zi = 1;
              new_prime = true;
            } else if (zi > 0) {
              // set new random
              for(uint zi = 2; zi <= n; zi ++){
                auto key = std::make_pair(zi, curr_zi_order[zi - 2]);
                set_new_rand(key);
              }
            }

            /*
             * If not finished, check if we can use some saved runs
             */
            try {
              std::vector<uint> tmp_vec = std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end() - 1);
              std::pair<FFInt, FFInt> key_val = saved_ti.at(tmp_vec).back();
              saved_ti.at(tmp_vec).pop_back();
              feed(key_val.first, key_val.second, tmp_vec, prime_number, lock);
            } catch (std::out_of_range& e) {
              // do nothing
            }

            return;
          }
        }
      } else if (n > 1 && feed_zi_ord != tmp_vec) {
        try {
          saved_ti.at(feed_zi_ord).emplace_back(std::make_pair(new_ti, num));
        } catch (std::out_of_range& e) {
          std::vector<std::pair<FFInt, FFInt>> tmp_ti = {std::make_pair(new_ti, num)};
          saved_ti[feed_zi_ord] = tmp_ti;
        }
      }
    }
  }

  std::tuple<int, uint, std::vector<uint>> RatReconst::feed_poly(int curr_deg,
                                                                 uint max_deg, std::unordered_map<uint, PolyReconst>& coef,
                                                                 PolyReconst& rec, ff_map_map& saved_num,
  std::unordered_map<uint, Polynomial>& sub_save, bool is_num) {
    uint tmp_zi = rec.next_zi + 1;
    std::vector<uint> tmp_zi_ord = curr_zi_order;

    while (!rec.new_prime) {
      try {
        std::vector<uint> key = {(uint) curr_deg, tmp_zi};
        FFInt food = saved_num.at(tmp_zi_ord).at(key);
        // delete unused saved data
        saved_num[tmp_zi_ord].erase(key);
        // set random values for the yis
        std::vector<FFInt> yis {};

        for (uint i = 0; i < tmp_zi_ord.size() - 1; i ++) {
          yis.emplace_back(rand_zi[std::make_pair(i + 2, tmp_zi_ord[i])]);
        }

        // feed to PolyReconst
        // since the constant is just a constant, we do not have to get mutliple
        // numerical values to reconstruct the coefficient
        if (curr_deg == 0) {
          while (!rec.new_prime) {
            if (curr_deg == (int) max_deg)
              rec.feed(yis, food);
            else {
              yis.emplace(yis.begin(), FFInt(1));
              FFInt num_subtraction = sub_save[curr_deg].convert_to_PolynomialFF().calc(yis);
              yis.erase(yis.begin());
              rec.feed(yis, food - num_subtraction);
            }
          }
        } else {
          if (curr_deg == (int) max_deg){
            if(prime_number == 0)
              rec.feed(yis, food);
            else if(is_num)
              rec.feed(degs_n[curr_deg], yis, food);
            else
              rec.feed(degs_d[curr_deg], yis, food);
          }
          else {
            yis.emplace(yis.begin(), FFInt(1));
            FFInt num_subtraction = sub_save[curr_deg].convert_to_PolynomialFF().calc(yis);
            yis.erase(yis.begin());
            if(prime_number == 0)
              rec.feed(yis, food - num_subtraction);
            else if(is_num)
              rec.feed(degs_n[curr_deg], yis, food - num_subtraction);
            else
              rec.feed(degs_d[curr_deg], yis, food - num_subtraction);
          }

          tmp_zi = rec.next_zi + 1;
          tmp_zi_ord = rec.curr_zi_order;
          tmp_zi_ord.emplace_back(prime_number);

        }
      } catch (std::out_of_range& e) {
        coef[curr_deg] = rec;
        return std::make_tuple(curr_deg, tmp_zi, tmp_zi_ord);
      }

      /*
       * Check for new prime & save subtraction term
       */
      if (rec.new_prime) {
        coef[curr_deg] = rec;

        if (curr_deg > 0) {
          Polynomial sub_pol = rec.get_result().homogenize(curr_deg).add_shift(shift);
          sub_pol -= rec.get_result().homogenize(curr_deg);

          for (auto & el : sub_pol.coefs) {
            uint tmp_deg = 0;

            for (auto & n : el.powers) {
              tmp_deg += n;
            }

            sub_save[tmp_deg] += el;
          }
        }

        /*
         * Remove already solved coefficients from Gauss eliminiation
         */
        if (is_num) {
          if (shift[0] != 0) {
            tmp_solved_coefs_num ++;
          }

          deg_num.pop_back();
          curr_deg = deg_num.back();
        } else {
          if (shift[0] != 0 && curr_deg > min_deg_den) {
            tmp_solved_coefs_den ++;
          }

          deg_den.pop_back();
          curr_deg = deg_den.back();
        }

        num_eqn = max_deg_den + max_deg_num + 1 - min_deg_num - min_deg_den
                  - solved_coefs - tmp_solved_coefs_num - tmp_solved_coefs_den;

        if (curr_deg >= 0) {
          rec = coef[curr_deg];
          std::fill(tmp_zi_ord.begin(), tmp_zi_ord.end() - 1, 1);
          tmp_zi_ord[n - 1] = prime_number;
          tmp_zi = rec.next_zi + 1;
        } else break;
      }
    }

    tmp_zi = rec.next_zi + 1;
    return std::make_tuple(curr_deg, tmp_zi, tmp_zi_ord);
  }

  void RatReconst::combine_primes(std::pair<mpz_map, mpz_map>& tmp) {
    if (!use_chinese_remainder) {
      combined_ni = tmp.first;
      combined_di = tmp.second;

      // if the coefficient is not a rational number thus divided by 1,
      // it will not change in the next run and can be omitted to save
      // numerical runs
      if (shift[0].n == 0) {
        mpz_map combined_ni_back = combined_ni;

        for (auto & c_ni : combined_ni_back) {
          try {
            RationalNumber rn = get_rational_coef(c_ni.second, combined_prime);

            if (rn.numerator == c_ni.second && rn.denominator == 1) {
              uint deg = 0;

              for (auto & el : c_ni.first) deg += el;

              remove_ni(deg, c_ni.first, rn);
            }
          } catch (std::exception& e) {
            // do nothing
          }
        }

        mpz_map combined_di_back = combined_di;

        for (auto & c_di : combined_di_back) {
          try {
            RationalNumber rn = get_rational_coef(c_di.second, combined_prime);

            if (rn.numerator == c_di.second && rn.denominator == 1) {
              uint deg = 0;

              for (auto & el : c_di.first) deg += el;

              if ((int) deg != min_deg_den)
                remove_di(deg, c_di.first, rn);
            }
          } catch (std::exception& e) {
            // do nothing
          }
        }

        combined_ni_back.clear();
        combined_di_back.clear();
      }
    } else {
      mpz_map combined_ni_back = combined_ni;
      mpz_map combined_di_back = combined_di;
      mpz_class combined_prime_back = combined_prime;
      std::pair<mpz_class, mpz_class> p1;
      std::pair<mpz_class, mpz_class> p2;
      std::pair<mpz_class, mpz_class> p3;

      //numerator
      for (auto it = combined_ni.begin(); it != combined_ni.end(); ++it) {
        p1 = std::make_pair(it->second, combined_prime);
        p2 = std::make_pair(tmp.first[it->first], FFInt::p);
        p3 = run_chinese_remainder(p1, p2);
        combined_ni[it->first] = p3.first;
      }

      // denominator
      for (auto it = combined_di.begin(); it != combined_di.end(); ++it) {
        p1 = std::make_pair(it->second, combined_prime);
        p2 = std::make_pair(tmp.second[it->first], FFInt::p);
        p3 = run_chinese_remainder(p1, p2);
        combined_di[it->first] = p3.first;
      }

      combined_prime = p3.second;

      // Remove already known coefficients from solve algorithm to save numerical runs
      if (shift[0].n == 0) {
        for (auto & c_ni : combined_ni_back) {
          uint deg = 0;

          for (auto & el : c_ni.first) deg += el;

          try {
            RationalNumber last_rn = get_rational_coef(c_ni.second, combined_prime_back);
            RationalNumber curr_rn = get_rational_coef(combined_ni[c_ni.first], combined_prime);

            if (last_rn == curr_rn)
              remove_ni(deg, c_ni.first, curr_rn);
          } catch (std::exception& e) {
            if (c_ni.second == combined_ni[c_ni.first]) {
              RationalNumber rn = RationalNumber(c_ni.second, 1);
              remove_ni(deg, c_ni.first, rn);
            }
          }
        }

        for (auto & c_di : combined_di_back) {
          uint deg = 0;

          for (auto & el : c_di.first) deg += el;

          if ((int) deg != min_deg_den) {
            try {
              RationalNumber last_rn = get_rational_coef(c_di.second, combined_prime_back);
              RationalNumber curr_rn = get_rational_coef(combined_di[c_di.first], combined_prime);

              if (last_rn == curr_rn)
                remove_di(deg, c_di.first, curr_rn);
            } catch (std::exception& e) {
              if (c_di.second == combined_di[c_di.first]) {
                RationalNumber rn = RationalNumber(c_di.second, 1);
                remove_di(deg, c_di.first, rn);
              }
            }
          }
        }

        combined_ni_back.clear();
        combined_di_back.clear();
        combined_prime_back = 0;
      }
    }

    tmp_solved_coefs_num = 0;
    tmp_solved_coefs_den = 0;
    num_eqn = max_deg_den + max_deg_num + 1 - min_deg_num - min_deg_den
              - solved_coefs - tmp_solved_coefs_num - tmp_solved_coefs_den;

    curr_deg_num = *std::max_element(non_solved_coef_num.begin(), non_solved_coef_num.end());
    curr_deg_den = *std::max_element(non_solved_coef_den.begin(), non_solved_coef_den.end());

    sub_num.clear();
    sub_den.clear();
  }

  RationalFunction RatReconst::get_result() {
    std::unique_lock<std::mutex> lock_status(mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_feed(mutex_feed, std::defer_lock);
    std::lock(lock_status, lock_feed);

    if (done) {
      if (result.numerator.coefs.empty()) {
        Polynomial numerator;
        Polynomial denominator;

        numerator = Polynomial(g_ni);
        denominator = Polynomial(g_di);
        g_ni.clear();
        g_di.clear();

        numerator.sort();
        denominator.sort();
        result = RationalFunction(numerator, denominator);

        RationalNumber first_coef = result.denominator.coefs[0].coef;

        if (first_coef.numerator != 1 || first_coef.denominator != 1) result = normalize(result);
      }

      return result;
    } else {
      throw std::runtime_error("Access to unfinished result");
    }
  }

  bool RatReconst::rec_rat_coef() {
    bool run_test = true;

    for (const auto ci : combined_ni) {
      mpz_class a = ci.second;

      try {
        g_ni[ci.first] = get_rational_coef(a, combined_prime);
      } catch (const std::exception&) {
        run_test = false;
        break;
      }
    }

    for (const auto ci : combined_di) {
      mpz_class a = ci.second;

      try {
        g_di[ci.first] = get_rational_coef(a, combined_prime);
      } catch (const std::exception&) {
        run_test = false;
        break;
      }
    }

    return run_test;
  }

  FFInt RatReconst::comp_ai(int i, int ip, const FFInt& num) {
    if (ip == 0) {
      return num;
    } else {
      FFInt ai_i = comp_ai(i, ip - 1, num);
      return (ti[i] - ti[ip - 1]) / (ai_i - ai[ip - 1]);
    }
  }

  FFInt RatReconst::comp_fyi(uint i, uint ip, const FFInt& y) {
    if (ip == 0) {
      return ai[i];
    } else {
      return ai[i - ip] + (-ti[i - ip] + y) / comp_fyi(i, ip - 1, y);
    }
  }

  std::pair<PolynomialFF, PolynomialFF> RatReconst::solve_gauss() {
    std::vector<FFInt> results = solve_gauss_system(num_eqn, coef_mat);
    coef_mat.clear();

    // Bring result in canonical form
    ff_map numerator;
    ff_map denominator;

    const std::vector<uint> min_power = {(uint) min_deg_den};
    denominator.emplace(std::make_pair(std::move(min_power), FFInt(1)));

    uint non_solved_num_size = non_solved_coef_num.size();

    if (shift[0] != FFInt(0)) non_solved_num_size -= tmp_solved_coefs_num;

    if (non_solved_num_size == 0) {
      const std::vector<uint> min_power_num = {(uint) min_deg_num};
      numerator.emplace(std::make_pair(std::move(min_power_num), FFInt(0)));
    } else {
      for (uint i = 0; i < non_solved_num_size; i ++) {
        std::vector<uint> power = {non_solved_coef_num[i]};
        numerator.emplace(std::make_pair(std::move(power), results[i]));
      }
    }

    for (int i = 0; i < (int)(non_solved_coef_den.size() - 1 - tmp_solved_coefs_den); i ++) {
      uint pow = non_solved_coef_den[i];

      if ((int) pow != min_deg_den) {
        std::vector<uint> power = {pow};
        denominator.emplace(std::make_pair(std::move(power), results[i + non_solved_num_size]));
      }
    }

    return std::make_pair(PolynomialFF(1, numerator), PolynomialFF(1, denominator));
  }

  std::pair<PolynomialFF, PolynomialFF> RatReconst::construct_canonical() {
    if (ai.size() == 1) {
      ff_map numerator_ff;
      std::vector<uint> zero_deg = {0};
      numerator_ff.emplace(std::make_pair(zero_deg, ai[0]));
      ff_map denominator_ff;
      denominator_ff.emplace(std::make_pair(zero_deg, FFInt(1)));
      return std::make_pair(PolynomialFF(1, numerator_ff), PolynomialFF(1, denominator_ff));
    } else {
      std::pair<PolynomialFF, PolynomialFF> r = iterate_canonical(1);
      FFInt mti = -ti[0];
      std::pair<PolynomialFF, PolynomialFF> ratFun(r.first * ai[0] + r.second * mti + r.second.mul(1),
                                                   r.first);
      return ratFun;
    }
  }

  std::pair<PolynomialFF, PolynomialFF> RatReconst::iterate_canonical(uint i) {
    if (i < ai.size() - 1) {
      std::pair<PolynomialFF, PolynomialFF> fnp1 = iterate_canonical(i + 1);
      FFInt mti = -ti[i];
      return std::pair<PolynomialFF, PolynomialFF> (fnp1.first * ai[i] + fnp1.second.mul(1) + fnp1.second * mti,
                                                    fnp1.first);
    } else {
      ff_map numerator_ff;
      std::vector<uint> zero_deg = {0};
      numerator_ff.emplace(std::make_pair(zero_deg, ai[i]));
      ff_map denominator_ff;
      denominator_ff.emplace(std::make_pair(zero_deg, FFInt(1)));
      return std::make_pair(PolynomialFF(1, numerator_ff), PolynomialFF(1, denominator_ff));
    }
  }

  RationalFunction RatReconst::normalize(RationalFunction& rf) {
    RationalNumber equializer = rf.denominator.coefs[0].coef;
    RationalNumber terminator(equializer.denominator, equializer.numerator);

    rf.numerator = rf.numerator * terminator;
    rf.denominator = rf.denominator * terminator;
    return rf;
  }

  std::pair<mpz_map, mpz_map> RatReconst::convert_to_mpz(const std::pair<PolynomialFF, PolynomialFF>& rf) const {
    mpz_map ci_mpz_1;
    mpz_map ci_mpz_2;

    for (const auto coef : rf.first.coef) {
      ci_mpz_1.emplace(coef.first, mpz_class(coef.second.n));
    }

    for (const auto coef : rf.second.coef) {
      ci_mpz_2.emplace(coef.first, mpz_class(coef.second.n));
    }

    return std::make_pair(ci_mpz_1, ci_mpz_2);
  }

  ff_map RatReconst::convert_to_ffint(const rn_map& ri) const {
    ff_map gi_ffi;

    for (const auto & g_i : ri) {
      mpz_class tmp(g_i.second.numerator % FFInt::p);

      if (tmp < 0) tmp = tmp + FFInt::p;

      FFInt n(std::stoull(tmp.get_str()));

      tmp = g_i.second.denominator % FFInt::p;

      FFInt d(std::stoull(tmp.get_str()));

      gi_ffi.emplace(std::make_pair(g_i.first, n / d));
    }

    return gi_ffi;
  }

  bool RatReconst::test_guess(const FFInt& num) {
    ff_map g_ff_ni = convert_to_ffint(g_ni);
    ff_map g_ff_di = convert_to_ffint(g_di);
    PolynomialFF g_ny(n, g_ff_ni);
    PolynomialFF g_dy(n, g_ff_di);
    std::vector<FFInt> yis = std::vector<FFInt> (n, FFInt(1));
    yis[0] = ti[0];

    for (uint i = 1; i < n; i++) {
      yis[i] = yis[0] * rand_zi[std::make_pair(i + 1, curr_zi_order[i - 1])] + shift[i];
    }

    yis[0] += shift[0];

    return (g_ny.calc(yis) / g_dy.calc(yis)) == num;
  }

  void RatReconst::remove_ni(uint deg, const std::vector<uint>& deg_vec, RationalNumber& rn) {
    g_ni[deg_vec] =  rn;
    combined_ni.erase(deg_vec);
    solved_coefs_num[deg] += Monomial(deg_vec, rn);
    bool remove = true;

    for (auto & c_ni_test : combined_ni) {
      uint deg_test = 0;

      for (auto & el : c_ni_test.first) deg_test += el;

      if (deg_test == deg) remove = false;
    }

    if (remove) {
      solved_coefs ++;
      non_solved_coef_num.erase(std::remove(non_solved_coef_num.begin(),
                                            non_solved_coef_num.end(), deg),
                                non_solved_coef_num.end());
    }
  }

  void RatReconst::remove_di(uint deg, const std::vector<uint>& deg_vec, RationalNumber& rn) {
    g_di[deg_vec] =  rn;
    combined_di.erase(deg_vec);
    solved_coefs_den[deg - 1 - min_deg_den] += Monomial(deg_vec, rn);
    bool remove = true;

    for (auto & c_di_test : combined_di) {
      uint deg_test = 0;

      for (auto & el : c_di_test.first) deg_test += el;

      if (deg_test == deg) remove = false;
    }

    if (remove) {
      solved_coefs ++;
      non_solved_coef_den.erase(std::remove(non_solved_coef_den.begin(),
                                            non_solved_coef_den.end(), deg),
                                non_solved_coef_den.end());
    }
  }

  //TODO allow for seed with std::srand(std::time(nullptr));
  FFInt RatReconst::get_rand() {
    return FFInt(std::rand() % (FFInt::p - 1)) + FFInt(1);
  }

  void RatReconst::set_new_rand(std::pair<uint, uint>& key) {
    try {
      rand_zi.at(key);
    } catch (std::out_of_range& e) {
      rand_zi.emplace(std::make_pair(key, get_rand()));
    }
  }

  uint RatReconst::get_num_eqn() {
    std::unique_lock<std::mutex> lock(mutex_status);
    return num_eqn;
  }

  void RatReconst::generate_anchor_points() {
    std::unique_lock<std::mutex> lock_status(mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_feed(mutex_feed, std::defer_lock);
    std::lock(lock_status, lock_feed);

    rand_zi.clear();
    anchor_points.clear();

    for (uint i = 2; i <= n; i++) {
      const FFInt rand = get_rand();
      rand_zi.emplace(std::make_pair(std::make_pair(i, 1), rand));
      anchor_points.emplace_back(rand);
    }
  }

  bool RatReconst::is_done() {
    std::unique_lock<std::mutex> lock(mutex_status);
    return done;
  }

  uint RatReconst::get_prime() {
    std::unique_lock<std::mutex> lock(mutex_status);
    return prime_number;
  }

  std::vector<uint> RatReconst::get_zi_order() {
    std::unique_lock<std::mutex> lock(mutex_status);
    if (n == 1) {
      return std::vector<uint> {};
    } else {
      return std::vector<uint>(curr_zi_order.begin(), curr_zi_order.end() - 1);
    }
  }

  uint RatReconst::get_zi() {
    std::unique_lock<std::mutex> lock(mutex_status);
    return zi;
  }

  RatReconst::RatReconst(const RatReconst& other) {
    std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_my_feed(mutex_feed, std::defer_lock);
    std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_other_feed(other.mutex_feed, std::defer_lock);
    std::lock(lock_my_status, lock_my_feed, lock_other_status, lock_other_feed);

    n = other.n;
    done = other.done;
    zi = other.zi;
    prime_number = other.prime_number;
    curr_zi_order = other.curr_zi_order;
    check = other.check;
    use_chinese_remainder = other.use_chinese_remainder;
    new_prime = other.new_prime;
    first_run = other.first_run;
    coef_mat = other.coef_mat;
    curr_zi = other.curr_zi;
    saved_ti = other.saved_ti;
    ai = other.ai;
    degs_n = other.degs_n;
    degs_d = other.degs_d;
    coef_n = other.coef_n;
    coef_d = other.coef_d;
    deg_num = other.deg_num;
    deg_den = other.deg_den;
    non_solved_coef_num = other.non_solved_coef_num;
    non_solved_coef_den = other.non_solved_coef_den;
    sub_num = other.sub_num;
    sub_den = other.sub_den;
    solved_coefs_num = other.solved_coefs_num;
    solved_coefs_den = other.solved_coefs_den;
    saved_num_num = other.saved_num_num;
    saved_num_den = other.saved_num_den;
    max_deg_num = other.max_deg_num;
    max_deg_den = other.max_deg_den;
    min_deg_den = other.min_deg_den;
    min_deg_num = other.min_deg_num;
    curr_deg_num = other.curr_deg_num;
    curr_deg_den = other.curr_deg_den;
    curr_zi_order_num = other.curr_zi_order_num;
    curr_zi_order_den = other.curr_zi_order_den;
    solved_coefs = other.solved_coefs;
    tmp_solved_coefs_num = other.tmp_solved_coefs_num;
    tmp_solved_coefs_den = other.tmp_solved_coefs_den;
    num_eqn = other.num_eqn;
    result = other.result;
    combined_prime = other.combined_prime;
    ti = other.ti;
    g_ni = other.g_ni;
    g_di = other.g_di;
    combined_ni = other.combined_ni;
    combined_di = other.combined_di;
  }

  RatReconst::RatReconst(RatReconst&& other) {
    std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_my_feed(mutex_feed, std::defer_lock);
    std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_other_feed(other.mutex_feed, std::defer_lock);
    std::lock(lock_my_status, lock_my_feed, lock_other_status, lock_other_feed);

    n = std::move(other.n);
    done = std::move(other.done);
    zi = std::move(other.zi);
    prime_number = std::move(other.prime_number);
    curr_zi_order = std::move(other.curr_zi_order);
    check = std::move(other.check);
    use_chinese_remainder = std::move(other.use_chinese_remainder);
    new_prime = std::move(other.new_prime);
    first_run = std::move(other.first_run);
    coef_mat = std::move(other.coef_mat);
    curr_zi = std::move(other.curr_zi);
    saved_ti = std::move(other.saved_ti);
    ai = std::move(other.ai);
    degs_n = std::move(other.degs_n);
    degs_d = std::move(other.degs_d);
    coef_n = std::move(other.coef_n);
    coef_d = std::move(other.coef_d);
    deg_num = std::move(other.deg_num);
    deg_den = std::move(other.deg_den);
    non_solved_coef_num = std::move(other.non_solved_coef_num);
    non_solved_coef_den = std::move(other.non_solved_coef_den);
    sub_num = std::move(other.sub_num);
    sub_den = std::move(other.sub_den);
    solved_coefs_num = std::move(other.solved_coefs_num);
    solved_coefs_den = std::move(other.solved_coefs_den);
    saved_num_num = std::move(other.saved_num_num);
    saved_num_den = std::move(other.saved_num_den);
    max_deg_num = std::move(other.max_deg_num);
    max_deg_den = std::move(other.max_deg_den);
    min_deg_den = std::move(other.min_deg_den);
    min_deg_num = std::move(other.min_deg_num);
    curr_deg_num = std::move(other.curr_deg_num);
    curr_deg_den = std::move(other.curr_deg_den);
    curr_zi_order_num = std::move(other.curr_zi_order_num);
    curr_zi_order_den = std::move(other.curr_zi_order_den);
    solved_coefs = std::move(other.solved_coefs);
    tmp_solved_coefs_num = std::move(other.tmp_solved_coefs_num);
    tmp_solved_coefs_den = std::move(other.tmp_solved_coefs_den);
    num_eqn = std::move(other.num_eqn);
    result = std::move(other.result);
    combined_prime = std::move(other.combined_prime);
    ti = std::move(other.ti);
    g_ni = std::move(other.g_ni);
    g_di = std::move(other.g_di);
    combined_ni = std::move(other.combined_ni);
    combined_di = std::move(other.combined_di);
  }

  RatReconst& RatReconst::operator=(const RatReconst& other) {
    if (this != &other) {
      std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
      std::unique_lock<std::mutex> lock_my_feed(mutex_feed, std::defer_lock);
      std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
      std::unique_lock<std::mutex> lock_other_feed(other.mutex_feed, std::defer_lock);
      std::lock(lock_my_status, lock_my_feed, lock_other_status, lock_other_feed);

      n = other.n;
      done = other.done;
      zi = other.zi;
      prime_number = other.prime_number;
      curr_zi_order = other.curr_zi_order;
      check = other.check;
      use_chinese_remainder = other.use_chinese_remainder;
      new_prime = other.new_prime;
      first_run = other.first_run;
      coef_mat = other.coef_mat;
      curr_zi = other.curr_zi;
      saved_ti = other.saved_ti;
      ai = other.ai;
      degs_n = other.degs_n;
      degs_d = other.degs_d;
      coef_n = other.coef_n;
      coef_d = other.coef_d;
      deg_num = other.deg_num;
      deg_den = other.deg_den;
      non_solved_coef_num = other.non_solved_coef_num;
      non_solved_coef_den = other.non_solved_coef_den;
      sub_num = other.sub_num;
      sub_den = other.sub_den;
      solved_coefs_num = other.solved_coefs_num;
      solved_coefs_den = other.solved_coefs_den;
      saved_num_num = other.saved_num_num;
      saved_num_den = other.saved_num_den;
      max_deg_num = other.max_deg_num;
      max_deg_den = other.max_deg_den;
      min_deg_den = other.min_deg_den;
      min_deg_num = other.min_deg_num;
      curr_deg_num = other.curr_deg_num;
      curr_deg_den = other.curr_deg_den;
      curr_zi_order_num = other.curr_zi_order_num;
      curr_zi_order_den = other.curr_zi_order_den;
      solved_coefs = other.solved_coefs;
      tmp_solved_coefs_num = other.tmp_solved_coefs_num;
      tmp_solved_coefs_den = other.tmp_solved_coefs_den;
      num_eqn = other.num_eqn;
      result = other.result;
      combined_prime = other.combined_prime;
      ti = other.ti;
      g_ni = other.g_ni;
      g_di = other.g_di;
      combined_ni = other.combined_ni;
      combined_di = other.combined_di;
    }
    return *this;
  }

  RatReconst& RatReconst::operator=(RatReconst&& other) {
    if (this != &other) {
      std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
      std::unique_lock<std::mutex> lock_my_feed(mutex_feed, std::defer_lock);
      std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
      std::unique_lock<std::mutex> lock_other_feed(other.mutex_feed, std::defer_lock);
      std::lock(lock_my_status, lock_my_feed, lock_other_status, lock_other_feed);

      n = std::move(other.n);
      done = std::move(other.done);
      zi = std::move(other.zi);
      prime_number = std::move(other.prime_number);
      curr_zi_order = std::move(other.curr_zi_order);
      check = std::move(other.check);
      use_chinese_remainder = std::move(other.use_chinese_remainder);
      new_prime = std::move(other.new_prime);
      first_run = std::move(other.first_run);
      coef_mat = std::move(other.coef_mat);
      curr_zi = std::move(other.curr_zi);
      saved_ti = std::move(other.saved_ti);
      ai = std::move(other.ai);
      degs_n = std::move(other.degs_n);
      degs_d = std::move(other.degs_d);
      coef_n = std::move(other.coef_n);
      coef_d = std::move(other.coef_d);
      deg_num = std::move(other.deg_num);
      deg_den = std::move(other.deg_den);
      non_solved_coef_num = std::move(other.non_solved_coef_num);
      non_solved_coef_den = std::move(other.non_solved_coef_den);
      sub_num = std::move(other.sub_num);
      sub_den = std::move(other.sub_den);
      solved_coefs_num = std::move(other.solved_coefs_num);
      solved_coefs_den = std::move(other.solved_coefs_den);
      saved_num_num = std::move(other.saved_num_num);
      saved_num_den = std::move(other.saved_num_den);
      max_deg_num = std::move(other.max_deg_num);
      max_deg_den = std::move(other.max_deg_den);
      min_deg_den = std::move(other.min_deg_den);
      min_deg_num = std::move(other.min_deg_num);
      curr_deg_num = std::move(other.curr_deg_num);
      curr_deg_den = std::move(other.curr_deg_den);
      curr_zi_order_num = std::move(other.curr_zi_order_num);
      curr_zi_order_den = std::move(other.curr_zi_order_den);
      solved_coefs = std::move(other.solved_coefs);
      tmp_solved_coefs_num = std::move(other.tmp_solved_coefs_num);
      tmp_solved_coefs_den = std::move(other.tmp_solved_coefs_den);
      num_eqn = std::move(other.num_eqn);
      result = std::move(other.result);
      combined_prime = std::move(other.combined_prime);
      ti = std::move(other.ti);
      g_ni = std::move(other.g_ni);
      g_di = std::move(other.g_di);
      combined_ni = std::move(other.combined_ni);
      combined_di = std::move(other.combined_di);
    }
    return *this;
  }
}
