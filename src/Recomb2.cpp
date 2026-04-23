#include "Hadronizer.h"
#include "PhysicsLists.h"
#include "Lorentz.h"
#include <cmath>
#include <algorithm>
#include <iostream>



bool Hadronizer::tryRecombine2(const Particle& HQ, int lightFlavor, int orbital,
                               int pdgWanted, double mWanted,
                               Particle& outHadron)
{
  int todoReco = 1;
  return forceRecomb2(todoReco, HQ, lightFlavor, orbital, pdgWanted, mWanted, outHadron);
}

bool Hadronizer::forceRecomb2(int& todoReco, const Particle& HQ,
                              int lightFlavor, int orbital,
                              int pdgWanted, double mWanted,
                              Particle& outHadron)
{
  
  const double temp = HQ.Thydro;
  if (temp < 0.1) return false;

  // set mlight based on flavor (1=ud, 3=s)
  double mlight = 0.0;
  if (lightFlavor == 1) mlight = cfg_.m_ud;
  else if (lightFlavor == 3) mlight = cfg_.m_s;
  else {
    std::cerr << "forceRecomb2: unsupported lightFlavor=" << lightFlavor << "\n";
    return false;
  }

  // compute kmin2 estimation (used only to set probmax)
  // plmax = 15*temp
  const double plmax = 15.0 * temp;
  const double Elmax = std::sqrt(plmax*plmax + mlight*mlight);

  const double becell = std::sqrt(sqr(HQ.cvx)+sqr(HQ.cvy)+sqr(HQ.cvz));
  const double gacell = 1.0 / std::sqrt(std::max(1e-15, 1.0 - becell*becell));

  const double plmaxcm = gacell*(plmax + becell*Elmax);
  const double Elmaxcm = gacell*(Elmax + becell*plmax);

  const double pHQ = std::sqrt(std::max(0.0, HQ.E*HQ.E - HQ.m*HQ.m));
  const double becm = (pHQ + plmaxcm) / (HQ.E + Elmaxcm);
  const double gacm = 1.0 / std::sqrt(std::max(1e-15, 1.0 - becm*becm));

  const double kmin = gacm*(pHQ - becm*HQ.E);
  double kmin2 = kmin*kmin;
  if (pHQ < plmaxcm) kmin2 = -1.0;

  
  const double reducedM0 = mlight*HQ.m / (mlight + HQ.m);
  const double sigma0 = 1.0 / std::sqrt(reducedM0 * omega_HM_);

  
  int ind_temp = int((temp - 0.160) / 0.002 + 0.5);
  if (ind_temp < 0) return false;
  if (ind_temp > 35) ind_temp = 35;
  int idx = std::clamp(ind_temp * 2, 0, PhysicsLists::kMaxIndex);

  const double mxvalue = (lightFlavor == 1)
    ? PhysicsLists::max_ud[idx]
    : PhysicsLists::max_s[idx];

  
  double directM=0, vx_HM=0, vy_HM=0, vz_HM=0;
  double px_l=0, py_l=0, pz_l=0, E_l=0;

  long long loopcount = 0;
  while (todoReco == 1) {
    judgeHadMech2(orbital, todoReco, HQ, mlight, sigma0, kmin2, mxvalue,
                  directM, vx_HM, vy_HM, vz_HM,
                  px_l, py_l, pz_l, E_l);

    if (++loopcount > 100000000LL) {
      std::cerr << "forceRecomb2: too many loops, no recomb\n";
      return false;
    }
  }

  if (todoReco != 2) {
    // did not recombine
    return false;
  }

  //Fill output hadron 
  outHadron = Particle{};
  outHadron.x = HQ.x;
  outHadron.y = HQ.y;
  outHadron.z = HQ.z;
  outHadron.t = HQ.t;
  outHadron.Thydro = HQ.Thydro;
  outHadron.cvx = HQ.cvx;
  outHadron.cvy = HQ.cvy;
  outHadron.cvz = HQ.cvz;
  outHadron.ipx = HQ.ipx;
  outHadron.ipy = HQ.ipy;
  outHadron.ipz = HQ.ipz;
  outHadron.iE = HQ.iE;
  outHadron.wt  = HQ.wt;

  outHadron.pdg = pdgWanted;
  outHadron.m   = mWanted;

  // recombine in 2-body CM first: p=(0,0,0,E=directM), then boost by vx_HM
  outHadron.px = 0.0;
  outHadron.py = 0.0;
  outHadron.pz = 0.0;
  outHadron.E  = directM;

  Lorentz::boost(vx_HM, vy_HM, vz_HM, outHadron.px, outHadron.py, outHadron.pz, outHadron.E);

  // boost by cell velocity back to global CM
  Lorentz::boost(HQ.cvx, HQ.cvy, HQ.cvz, outHadron.px, outHadron.py, outHadron.pz, outHadron.E);

  // decay2body:adjust mass shell by emitting "pion/photon"
  decay2body(outHadron.px, outHadron.py, outHadron.pz, outHadron.E, outHadron.m);

  return true;
}


/////// judgeHadMech2: 

void Hadronizer::judgeHadMech2(int orbital, int& todoReco,
                               const Particle& HQ,
                               double mlight,
                               double sigma0, double kmin2,
                               double mxvalue,
                               double& directM,
                               double& vx_HM, double& vy_HM, double& vz_HM,
                               double& px_l, double& py_l, double& pz_l, double& E_l)
{
  const double temp = HQ.Thydro;

  // get probmax
  double probmax = 1.0;
  if (orbital == 0) {
    probmax = (kmin2 < 0.0) ? 1.0 : std::exp(-kmin2 * sigma0 * sigma0);
  } else if (orbital == 1) {
    const double x = kmin2 * sigma0 * sigma0;
    probmax = (x < 1.0) ? std::exp(-1.0) : x * std::exp(-x);
  } else {
    std::cerr << "judgeHadMech2: only s/p waves supported\n";
    todoReco = 1;
    return;
  }

  // sample thermal light parton (rejection)
  const double pmax = 15.0 * temp;

  auto distribution = [&](double p) -> double {
    const double Eq = std::sqrt(p*p + mlight*mlight);
    const double Eg = std::sqrt(p*p + cfg_.m_g*cfg_.m_g);
    const double term_q = 6.0 * p*p / (std::exp(Eq/temp) + 1.0);
    const double term_g = (p*p / (std::exp(Eg/temp) - 1.0)) * (16.0/6.0);
    return term_q + term_g;
  };

  double pabs = rng_.uniform() * pmax;
  while (rng_.uniform() * mxvalue > distribution(pabs)) {
    pabs = rng_.uniform() * pmax;
  }

  const double u   = rng_.uniform()*2.0 - 1.0;
  const double phi = rng_.uniform() * 2.0 * 3.14159265358979323846;
  const double s   = std::sqrt(std::max(0.0, 1.0 - u*u));

  px_l = pabs * s * std::cos(phi);
  py_l = pabs * s * std::sin(phi);
  pz_l = pabs * u;
  E_l  = std::sqrt(pabs*pabs + mlight*mlight);

  // pair beta in LAB
  const double cmbx = (HQ.px + px_l) / (HQ.E + E_l);
  const double cmby = (HQ.py + py_l) / (HQ.E + E_l);
  const double cmbz = (HQ.pz + pz_l) / (HQ.E + E_l);

  // boost both to CM frame
  double lpx=px_l, lpy=py_l, lpz=pz_l, lE=E_l;
  double hpx=HQ.px, hpy=HQ.py, hpz=HQ.pz, hE=HQ.E;

  Lorentz::boost(-cmbx, -cmby, -cmbz, lpx, lpy, lpz, lE);
  Lorentz::boost(-cmbx, -cmby, -cmbz, hpx, hpy, hpz, hE);

  const double denom = lE + hE;
  if (denom <= 0.0) { todoReco = 1; return; }

  const double qx = (lE*hpx - hE*lpx) / denom;
  const double qy = (lE*hpy - hE*lpy) / denom;
  const double qz = (lE*hpz - hE*lpz) / denom;
  const double q2 = qx*qx + qy*qy + qz*qz;

  // get probWigner
  double probWigner = 1.0;
  if (orbital == 0) probWigner = std::exp(-q2 * sigma0 * sigma0);
  else              probWigner = std::exp(-q2 * sigma0 * sigma0) * (q2 * sigma0 * sigma0);

  // accept/reject
  const bool pureRecomb = (cfg_.hadr_flag == 2);
  const bool accept = pureRecomb ? true : (rng_.uniform() * probmax < probWigner);

  if (accept) {
    todoReco = 2;
    directM = hE + lE;
    vx_HM = cmbx;
    vy_HM = cmby;
    vz_HM = cmbz;
    return;
  }

  todoReco = 1;
}

