#include "Hadronizer.h"
#include "PhysicsLists.h"
#include "RecombinationTable.h"
#include "Lorentz.h"
#include "Frag.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <random>

Hadronizer::Hadronizer(const Config& cfg, RNG& rng, const RecombinationTable& table)
: cfg_(cfg), rng_(rng),  table_(table), frag_()
{
  current_max_3body_s_ = cfg_.max_3body_s;
  current_max_3body_p_ = cfg_.max_3body_p;
  if (cfg_.HQid == 4) {
    m_HQ_ = cfg_.m_c;
    omega_HM_ = cfg_.omega_cM;

    // fragmentation settings for charm 
    frag_.setDefaultCharmChemistry();

  } else if (cfg_.HQid == 5) {
    m_HQ_ = cfg_.m_b;
    omega_HM_ = cfg_.omega_bM;

    // No bottom fragmentation chemistry in this HQET version so far.

  } else {
    // 
    m_HQ_ = cfg_.m_c;
    omega_HM_ = cfg_.omega_cM;
  }
}
//////////////////////
// returns true if recombination happened and outHadron filled

void Hadronizer::process(const std::vector<Particle>& input,
                         std::vector<Particle>& output)
{
  output.clear();
  output.reserve(input.size());

  for (const auto& p : input) {
	 
//	  std::cout << "p.wt= " << p.wt << std::endl;//test

    // only act on heavy quarks
    if (!(std::abs(p.pdg) == 4 || std::abs(p.pdg) == 5)) {
      output.push_back(p);
      continue;
    }

    Particle had;
    // -----------------------
    // hadr_flag == 1 : pure fragmentation
    // -----------------------
   // std::cout << "hadr_flag = " << cfg_.hadr_flag << std::endl; //test
    if (cfg_.hadr_flag == 1) {
	    //need to implement fragmentation; now just return original HQ
	    if (frag_.fragment(p, rng_, had)) {
                      output.push_back(had);
//		      std::cout<< " okay" << std::endl; //test
              }else
              {output.push_back(p);}
      //output.push_back(p);
      continue;
    }

    //Particle had;
    bool did = false;

    // ============================================================
    // hadr_flag == 3 : frag + recomb (table-driven)
    // ============================================================
    if (cfg_.hadr_flag == 3) {
      Particle hqReco = p;
      const double originalMass = hqReco.m;
      if (cfg_.rescale_m == 1) {
        hqReco.m = m_HQ_;
        hqReco.E = std::sqrt(std::max(0.0, hqReco.px*hqReco.px +
                                           hqReco.py*hqReco.py +
                                           hqReco.pz*hqReco.pz +
                                           hqReco.m *hqReco.m));
      }
      Lorentz::boost(-hqReco.cvx, -hqReco.cvy, -hqReco.cvz,
                     hqReco.px, hqReco.py, hqReco.pz, hqReco.E);
//	    std::cout << " check " << std::endl;
//           std::cout << "check hqReco " << hqReco.px << " " <<  hqReco.py << " " << hqReco.pz << " " << hqReco.E << std::endl;
      std::vector<double> cdf;
      if (!interpolatedCDF(hqReco, cdf)) {
		      //std::cout << " must frag: outside the pmax-bound" << std::endl; //test
        // outside table => recomb_prob = 0
//	      std::cout << "what!" <<std::endl;
          Lorentz::boost(hqReco.cvx, hqReco.cvy, hqReco.cvz,
                         hqReco.px, hqReco.py, hqReco.pz, hqReco.E);
          if (cfg_.rescale_m == 1) {
            hqReco.m = originalMass;
            hqReco.E = std::sqrt(std::max(0.0, hqReco.px*hqReco.px +
                                               hqReco.py*hqReco.py +
                                               hqReco.pz*hqReco.pz +
                                               hqReco.m *hqReco.m));
          }

		      if (frag_.fragment(hqReco, rng_, had)) {
			      output.push_back(had); 
		      }else 
		      {output.push_back(hqReco);}
		      
		      //  output.push_back(p);
	        continue;
      }
      if (cdf.size() < (size_t)cfg_.nChannel) {
//	      std::cout << "cdf.size= " << cdf.size() << "cfg_.nChannel= " << cfg_.nChannel <<std::endl; //test
        Lorentz::boost(hqReco.cvx, hqReco.cvy, hqReco.cvz,
                       hqReco.px, hqReco.py, hqReco.pz, hqReco.E);
        if (cfg_.rescale_m == 1) {
          hqReco.m = originalMass;
          hqReco.E = std::sqrt(std::max(0.0, hqReco.px*hqReco.px +
                                             hqReco.py*hqReco.py +
                                             hqReco.pz*hqReco.pz +
                                             hqReco.m *hqReco.m));
        }
        if (frag_.fragment(hqReco, rng_, had)) output.push_back(had);
        else output.push_back(hqReco);
        continue;
      }
//      std::cout << " check2 " << std::endl;

      // Fortran: recomb_random = rlu(0)
      const double recomb_random = rng_.uniform();

      // Fortran: recomb_prob = recomb_prob_br(nChannel)
      const double recomb_prob = cdf.back();

      // Fortran: if(recomb_random.lt.recomb_prob) then do recomb else frag
      //std::cout << "recomb_random: " << recomb_random << " " << "recomb_prob: " << recomb_prob << std::endl;
	      if (recomb_random >= recomb_prob) {
		      //std::cout << " must frag: MC determined" << std::endl; //test      
          Lorentz::boost(hqReco.cvx, hqReco.cvy, hqReco.cvz,
                         hqReco.px, hqReco.py, hqReco.pz, hqReco.E);
          if (cfg_.rescale_m == 1) {
            hqReco.m = originalMass;
            hqReco.E = std::sqrt(std::max(0.0, hqReco.px*hqReco.px +
                                               hqReco.py*hqReco.py +
                                               hqReco.pz*hqReco.pz +
                                               hqReco.m *hqReco.m));
          }
		      if (frag_.fragment(hqReco, rng_, had)) {
			      output.push_back(had);
		      }else
		      {output.push_back(hqReco);}

	      
	      //output.push_back(p);
        continue;
      }

      //  apply charge conjugation convention
      auto pdgSigned = [&](int pdgWanted) -> int {
        return (p.pdg < 0) ? -std::abs(pdgWanted) : std::abs(pdgWanted);
      };

	      auto FR2 = [&](int lf, int orb, int pdgWanted, double mWanted) -> bool {
	        return tryRecombine2(hqReco, lf, orb, pdgSigned(pdgWanted), mWanted, had);
	      };
	      auto FR3 = [&](int lf1, int lf2, int orb, int pdgWanted, double mWanted) -> bool {
	        return tryRecombine3(hqReco, lf1, lf2, orb, pdgSigned(pdgWanted), mWanted, had);
	      };
       //std::cout << " must reco" << std::endl; //test 
      //  recomb_prob_br(i) are cumulative thresholds.
      // Here: br(i) = cdf[i-1], i=1..nChannel
      auto br = [&](int i) -> double { return cdf.at((size_t)i - 1); };

      // ------------------------
      // Charm channels and Beauty channels 
      // ------------------------
      if (cfg_.HQid == 4) {
        // ----------------------------
        // charm
        // ----------------------------
        if (recomb_random < br(1)) {
          // (1) D0 or D+
          if (rng_.uniform() < 0.5) did = FR2(1, 0, 421, 1.86); // D0
          else                      did = FR2(1, 0, 411, 1.87); // D+
        }
        else if (recomb_random < br(2)) {
          // (2) D*0 or D*+
          if (rng_.uniform() < 0.5) did = FR2(1, 0, 423, 2.01); // D*0
          else                      did = FR2(1, 0, 413, 2.01); // D*+
        }
        else if (recomb_random < br(3)) {
          // (3) D1(2420)  (p-wave)
          did = FR2(1, 1, 10423, 2.42);
        }
        else if (recomb_random < br(4)) {
          // (4) D0*(2300) (p-wave)
          did = FR2(1, 1, 10411, 2.30);
        }
        else if (recomb_random < br(5)) {
          // (5) D1*(2440) (p-wave)  (PDG here you may refine; placeholder)
          did = FR2(1, 1, 20413, 2.44);
        }
        else if (recomb_random < br(6)) {
          // (6) D2*(~2460) (p-wave)
          did = FR2(1, 1, 415, 2.46);
        }
        else if (recomb_random < br(7)) {
          // (7) Ds(1968)
          did = FR2(3, 0, 431, 1.97);
        }
        else if (recomb_random < br(8)) {
          // (8) Ds*(2122)
          did = FR2(3, 0, 433, 2.12);
        }
        else if (recomb_random < br(9)) {
          // (9) Ds*(2536) (p-wave)
          did = FR2(3, 1, 10433, 2.54);
        }
        else if (recomb_random < br(10)) {
          // (10) Ds0(2317) (p-wave)
          did = FR2(3, 1, 10431, 2.32);
        }
        else if (recomb_random < br(11)) {
          // (11) Ds1(2460) (p-wave)
          did = FR2(3, 1, 20433, 2.46);
        }
        else if (recomb_random < br(12)) {
          // (12) Ds2(2573) (p-wave)
          did = FR2(3, 1, 435, 2.57);
        }
        else if (recomb_random < br(13)) {
          // (13) Lambda_c s-wave 1/2+
          did = FR3(1, 1, 0, 4122, 2.286);
        }
        else if (recomb_random < br(14)) {
          // (14) Lambda_c* s-wave 3/2+  
      	  did = FR3(1, 1, 0, 4124, 2.86);
        }
        else if (recomb_random < br(15)) {
          // (15) Lambda_c p-wave together 
          did = FR3(1, 1, 1, 14122, 2.60);
        }
        else if (recomb_random < br(16)) {
          // (16) Sigma_c s-wave 1/2+ 
          did = FR3(1, 1, 0, 4212, 2.46);
        }
        else if (recomb_random < br(17)) {
          // (17) Sigma_c s-wave 3/2+ 
          did = FR3(1, 1, 0, 4214, 2.52);
        }
        else if (recomb_random < br(18)) {
          // (18) Sigma_c p-wave together 
          did = FR3(1, 1, 1, 14212, 2.80);
        }
        else {
          // If table has more channels but charm logic ends here, fallback to no recomb
          did = false;
        }
      }
      else if (cfg_.HQid == 5) {
        // ----------------------------
        // beauty
        // ----------------------------
        if (recomb_random < br(1)) {
          // (1) B+ / B0 
          if (rng_.uniform() < 0.5) did = FR2(1, 0, 521, 5.28); // B+
          else                      did = FR2(1, 0, 511, 5.28); // B0
        }
        else if (recomb_random < br(2)) {
          // (2) B*  
          did = FR2(1, 0, 513, 5.32);
        }
        else if (recomb_random < br(3)) {
          // (3) B1(5721) (p-wave) 
          did = FR2(1, 1, 10513, 5.72);
        }
        else if (recomb_random < br(4)) {
          // (4) B0*(5721) (p-wave) 
          did = FR2(1, 1, 10511, 5.72);
        }
        else if (recomb_random < br(5)) {
          // (5) B1*(5747) (p-wave) 
          did = FR2(1, 1, 20513, 5.75);
        }
        else if (recomb_random < br(6)) {
          // (6) B2*(5747) (p-wave) 
          did = FR2(1, 1, 515, 5.75);
        }
        else {
          did = false;
        }
      }

	      if (did) output.push_back(had);
	      else {
            Lorentz::boost(hqReco.cvx, hqReco.cvy, hqReco.cvz,
                           hqReco.px, hqReco.py, hqReco.pz, hqReco.E);
            if (cfg_.rescale_m == 1) {
              hqReco.m = originalMass;
              hqReco.E = std::sqrt(std::max(0.0, hqReco.px*hqReco.px +
                                                 hqReco.py*hqReco.py +
                                                 hqReco.pz*hqReco.pz +
                                                 hqReco.m *hqReco.m));
            }
            if (frag_.fragment(hqReco, rng_, had)) output.push_back(had);
            else output.push_back(hqReco);
          }
	      continue;
	    }

    // ============================================================
    // hadr_flag == 2 : recomb-only
    // ============================================================
    if (cfg_.hadr_flag == 2) {
      std::vector<double> cdf;
      if (!interpolatedCDF(p, cdf) || cdf.empty()) {
        output.push_back(p);
        continue;
      }

      auto pdgSigned = [&](int pdgWanted) -> int {
        return (p.pdg < 0) ? -std::abs(pdgWanted) : std::abs(pdgWanted);
      };
      auto FR2 = [&](int lf, int orb, int pdgWanted, double mWanted) -> bool {
        return tryRecombine2(p, lf, orb, pdgSigned(pdgWanted), mWanted, had);
      };
      auto FR3 = [&](int lf1, int lf2, int orb, int pdgWanted, double mWanted) -> bool {
        return tryRecombine3(p, lf1, lf2, orb, pdgSigned(pdgWanted), mWanted, had);
      };
      auto br = [&](int i) -> double { return cdf.at((size_t)i - 1); };
      const double recomb_random = rng_.uniform();

      if (cfg_.HQid == 4) {
        if (recomb_random < br(1)) {
          did = (rng_.uniform() < 0.5) ? FR2(1, 0, 421, 1.86) : FR2(1, 0, 411, 1.87);
        } else if (recomb_random < br(2)) {
          did = (rng_.uniform() < 0.5) ? FR2(1, 0, 423, 2.01) : FR2(1, 0, 413, 2.01);
        } else if (recomb_random < br(3)) {
          did = FR2(1, 1, 10423, 2.42);
        } else if (recomb_random < br(4)) {
          did = FR2(1, 1, 10411, 2.30);
        } else if (recomb_random < br(5)) {
          did = FR2(1, 1, 20413, 2.44);
        } else if (recomb_random < br(6)) {
          did = FR2(1, 1, 415, 2.46);
        } else if (recomb_random < br(7)) {
          did = FR2(3, 0, 431, 1.97);
        } else if (recomb_random < br(8)) {
          did = FR2(3, 0, 433, 2.12);
        } else if (recomb_random < br(9)) {
          did = FR2(3, 1, 10433, 2.54);
        } else if (recomb_random < br(10)) {
          did = FR2(3, 1, 10431, 2.32);
        } else if (recomb_random < br(11)) {
          did = FR2(3, 1, 20433, 2.46);
        } else if (recomb_random < br(12)) {
          did = FR2(3, 1, 435, 2.57);
        } else if (recomb_random < br(13)) {
          did = FR3(1, 1, 0, 4122, 2.286);
        } else if (recomb_random < br(14)) {
          did = FR3(1, 1, 0, 4124, 2.86);
        } else if (recomb_random < br(15)) {
          did = FR3(1, 1, 1, 14122, 2.60);
        } else if (recomb_random < br(16)) {
          did = FR3(1, 1, 0, 4212, 2.46);
        } else if (recomb_random < br(17)) {
          did = FR3(1, 1, 0, 4214, 2.52);
        } else if (recomb_random < br(18)) {
          did = FR3(1, 1, 1, 14212, 2.80);
        }
      } else if (cfg_.HQid == 5) {
        if (recomb_random < br(1)) {
          did = (rng_.uniform() < 0.5) ? FR2(1, 0, 521, 5.28) : FR2(1, 0, 511, 5.28);
        } else if (recomb_random < br(2)) {
          did = FR2(1, 0, 513, 5.32);
        } else if (recomb_random < br(3)) {
          did = FR2(1, 1, 10513, 5.72);
        } else if (recomb_random < br(4)) {
          did = FR2(1, 1, 10511, 5.72);
        } else if (recomb_random < br(5)) {
          did = FR2(1, 1, 20513, 5.75);
        } else if (recomb_random < br(6)) {
          did = FR2(1, 1, 515, 5.75);
        }
      }

      if (did) output.push_back(had);
      else output.push_back(p);
      continue;
    }
    // fallback
    output.push_back(p);
  }
}

double Hadronizer::pLengthInCellLRF(const Particle& HQ) const {
  double px = HQ.px, py = HQ.py, pz = HQ.pz, E = HQ.E;
  // boost to cell local rest frame
  Lorentz::boost(-HQ.cvx, -HQ.cvy, -HQ.cvz, px, py, pz, E);
  return std::sqrt(std::max(0.0, px*px + py*py + pz*pz));
}
double Hadronizer::pLength(const Particle& HQ) const {
  double px = HQ.px, py = HQ.py, pz = HQ.pz, E = HQ.E;
  // boost to cell local rest frame
  // Lorentz::boost(-HQ.cvx, -HQ.cvy, -HQ.cvz, px, py, pz, E);
  return std::sqrt(std::max(0.0, px*px + py*py + pz*pz));
}

bool Hadronizer::interpolatedCDF(const Particle& HQ, std::vector<double>& cdf) const {
  const int nCh = table_.nChannel;
  if (nCh <= 0) return false;
  if (table_.pc_raw.size() < 2) return false;
  if (table_.prob_cum_ch.size() != table_.pc_raw.size()) return false;

  // choose bin width (prob_int_c / prob_int_b)
  const double prob_int = (cfg_.HQid == 4 ? cfg_.prob_int_c : cfg_.prob_int_b); //cfg_.prob_int_c=0.5
  // std::cout << "prob_int = " << prob_int  << std::endl; //test

  // plength in cell LRF
 // const double pL = pLengthInCellLRF(HQ);
   const double pL = pLength(HQ);
  // ind_p=int(pL/prob_int); ind_l=ind_p+1 (1-based)
  const int ind_p = static_cast<int>(pL / prob_int);
  const int i0 = ind_p;       // 0-based lower index
  const int i1 = ind_p + 1;   // 0-based upper index

  // bounds: if outside table -> recomb_prob=0 (Fortran ind_l>=prob_num)
  const int N = static_cast<int>(table_.pc_raw.size());
  //std::cout << "N= " << N<< std::endl; //test

  //const double pL_max = (N-1)*prob_int;  //test
  //std::cout << "recomb table pL_max ~ " << pL_max << "\n";   //test
  //std::cout << "pL=" << pL << " (threshold " << pL_max << ")\n"; //test

  if (i1 >= N) {return false;}
  if (i0 < 0) {return false;}

  //if (i1 >= N) { std::cout << "i1>=N ?! skip reco" << std::endl; return false;}
  //if (i0 < 0) { std::cout << "i0<0 this?" << std::endl; return false;}
  const double p0 = ind_p * prob_int;
  const double frac = (prob_int > 0.0) ? (pL - p0) / prob_int : 0.0;
  const double t = std::clamp(frac, 0.0, 1.0);
  cdf.assign(nCh, 0.0);
  for (int c = 0; c < nCh; ++c) {
    const double y0 = table_.prob_cum_ch[i0][c];
    const double y1 = table_.prob_cum_ch[i1][c];
    cdf[c] = y0 + (y1 - y0) * t;
  }
  // ensure monotonic & clamp to [0,1]
  double prev = 0.0;
  for (int c = 0; c < nCh; ++c) {
    cdf[c] = std::clamp(cdf[c], prev, 1.0);
    prev = cdf[c];
  }
  if (i0 < static_cast<int>(table_.wigner_max.size()) &&
      table_.wigner_max[i0].size() >= 6) {
    current_max_3body_s_ = table_.wigner_max[i0][4];
    current_max_3body_p_ = table_.wigner_max[i0][5];
  } else {
    current_max_3body_s_ = cfg_.max_3body_s;
    current_max_3body_p_ = cfg_.max_3body_p;
  }
  return true;
}

int Hadronizer::drawChannelFromCDF(const std::vector<double>& cdf, double u01) const {
  const int nCh = static_cast<int>(cdf.size());
  if (nCh == 0) return -1;
  // u01 assumed in [0,1)
  for (int c = 0; c < nCh; ++c) {
    if (u01 < cdf[c]) return c;
  }
  return -1;
}





// ===== decay2body  =====
void Hadronizer::decay2body(double& px, double& py, double& pz, double& E,
                            double m_final)
{
  const double m_pi = 0.13;

  // boost to rest frame of initial (px,py,pz,E)
  const double vx = px / E;
  const double vy = py / E;
  const double vz = pz / E;

  Lorentz::boost(-vx, -vy, -vz, px, py, pz, E);

  const double m_initial = E;

  double p_companion = 0.0;
  if (std::abs(m_final - m_initial) > m_pi) {
    const double term = (m_initial*m_initial + m_pi*m_pi - m_final*m_final) / (2.0*m_initial);
    p_companion = std::sqrt(std::max(0.0, term*term - m_pi*m_pi));
  } else {
    p_companion = std::abs(m_initial*m_initial - m_final*m_final) / (2.0*m_initial);
  }

  const double u = rng_.uniform()*2.0 - 1.0;
  const double phi = rng_.uniform() * 2.0 * 3.14159265358979323846;

  const double s = std::sqrt(std::max(0.0, 1.0 - u*u));
  px = -p_companion * s * std::cos(phi);
  py = -p_companion * s * std::sin(phi);
  pz = -p_companion * u;

  E = std::sqrt(px*px + py*py + pz*pz + m_final*m_final);

  // boost back
 Lorentz::boost(vx, vy, vz, px, py, pz, E);
}

/*
// ===== boost by beta vector (bx,by,bz) =====
void Hadronizer::boost(double bx, double by, double bz,
                       double& px, double& py, double& pz, double& E) const
{
  const double b2 = bx*bx + by*by + bz*bz;
  if (b2 < 1e-30) return;

  const double gamma = 1.0 / std::sqrt(std::max(1e-15, 1.0 - b2));
  const double bp = bx*px + by*py + bz*pz;
  const double gamma2 = (gamma - 1.0) / b2;

  const double px_new = px + gamma2*bp*bx + gamma*bx*E;
  const double py_new = py + gamma2*bp*by + gamma*by*E;
  const double pz_new = pz + gamma2*bp*bz + gamma*bz*E;
  const double E_new  = gamma*(E + bp);

  px = px_new; py = py_new; pz = pz_new; E = E_new;
}
*/
