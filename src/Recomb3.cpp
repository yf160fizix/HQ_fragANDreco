#include "Hadronizer.h"
#include "PhysicsLists.h"
#include "Lorentz.h"
#include <cmath>
#include <algorithm>
#include <iostream>


bool Hadronizer::tryRecombine3(const Particle& HQ,
                               int lightFlavor1, int lightFlavor2, int orbital,
                               int pdgWanted, double mWanted,
                               Particle& outHadron)
{
  int todoReco = 1;
  return forceRecomb3(todoReco, HQ,
                      lightFlavor1, lightFlavor2, orbital,
                      pdgWanted, mWanted,
                      outHadron);
}


// Core 3-body recombination 
bool Hadronizer::forceRecomb3(int& todoReco, const Particle& HQ,
                              int lightFlavor1, int lightFlavor2, int orbital,
                              int pdgWanted, double mWanted,
                              Particle& outHadron)
{
  const double temp = HQ.Thydro;
  if (temp < 0.1) return false;

  // light quark masses
  auto massLight = [&](int lf) -> double {
    if (lf == 1) return cfg_.m_ud;
    if (lf == 3) return cfg_.m_s;
    std::cerr << "forceRecomb3: unsupported lightFlavor=" << lf << "\n";
    return -1.0;
  };

  const double m1 = massLight(lightFlavor1);
  const double m2 = massLight(lightFlavor2);
  if (m1 < 0.0 || m2 < 0.0) return false;


  // Wigner sigmas
  const double omega_HB = (std::abs(HQ.pdg) == 4) ? cfg_.omega_cB : cfg_.omega_bB;

  double reducedM1, reducedM2;

  reducedM1 = m1*m2/(m1+m2);
  reducedM2 = (m1+m2)*HQ.m/(m1+m2+HQ.m);
  const double sigma11 = 1.0/std::sqrt(reducedM1*omega_HB);
  const double sigma12 = 1.0/std::sqrt(reducedM2*omega_HB);

  reducedM1 = m1*HQ.m/(m1+HQ.m);
  reducedM2 = (m1+HQ.m)*m2/(m1+m2+HQ.m);
  const double sigma21 = 1.0/std::sqrt(reducedM1*omega_HB);
  const double sigma22 = 1.0/std::sqrt(reducedM2*omega_HB);

  reducedM1 = m2*HQ.m/(m2+HQ.m);
  reducedM2 = (m2+HQ.m)*m1/(m1+m2+HQ.m);
  const double sigma31 = 1.0/std::sqrt(reducedM1*omega_HB);
  const double sigma32 = 1.0/std::sqrt(reducedM2*omega_HB);

 
  // mxvalue from temperature table
  int ind_temp = int((temp - 0.160)/0.002 + 0.5);
  if (ind_temp < 0) return false;
  if (ind_temp > 35) ind_temp = 35;
  const int idx = std::clamp(ind_temp*2, 0, PhysicsLists::kMaxIndex);

  const double mxvalue1 =
    (lightFlavor1 == 1) ? PhysicsLists::max_ud[idx] : PhysicsLists::max_s[idx];
  const double mxvalue2 =
    (lightFlavor2 == 1) ? PhysicsLists::max_ud[idx] : PhysicsLists::max_s[idx];


  double directM = 0.0;
  double vx_HM = 0.0, vy_HM = 0.0, vz_HM = 0.0;

  long long loopcount = 0;
  while (todoReco == 1) {

    judgeHadMech3(orbital, todoReco, HQ,
                  m1, m2,
                  mxvalue1, mxvalue2,
                  sigma11, sigma12,
                  sigma21, sigma22,
                  sigma31, sigma32,
                  directM, vx_HM, vy_HM, vz_HM);

    if (++loopcount > 100000000LL) {
      std::cerr << "forceRecomb3: too many loops, no recombination\n";
      return false;
    }
  }

  if (todoReco != 2) return false;


  // Fill output hadron
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

  // Kinematics: CM --> lab --> decay
  outHadron.px = 0.0;
  outHadron.py = 0.0;
  outHadron.pz = 0.0;
  outHadron.E  = directM;

  Lorentz::boost(vx_HM, vy_HM, vz_HM,
        outHadron.px, outHadron.py, outHadron.pz, outHadron.E);

  Lorentz::boost(HQ.cvx, HQ.cvy, HQ.cvz,
        outHadron.px, outHadron.py, outHadron.pz, outHadron.E);

  decay2body(outHadron.px, outHadron.py, outHadron.pz,
             outHadron.E, outHadron.m);

  return true;
}


// judgeHadMech3 
void Hadronizer::judgeHadMech3(int orbital, int& todoReco,
                               const Particle& HQ,
                               double m1, double m2,
                               double mxvalue1, double mxvalue2,
                               double sigma11, double sigma12,
                               double sigma21, double sigma22,
                               double sigma31, double sigma32,
                               double& directM,
                               double& vx_HM, double& vy_HM, double& vz_HM)
{
  const double temp = HQ.Thydro;
  const double pmax = 15.0 * temp;
  const double PI   = 3.14159265358979323846;

  auto dist = [&](double p, double m) -> double {
    const double Eq = std::sqrt(p*p + m*m);
    const double Eg = std::sqrt(p*p + cfg_.m_g*cfg_.m_g);
    const double term_q = 6.0 * p*p / (std::exp(Eq/temp) + 1.0);
    const double term_g = (p*p / (std::exp(Eg/temp) - 1.0)) * (16.0/6.0);
    return term_q + term_g;
  };

  auto sample_light = [&](double m, double mx,
                          double& px, double& py, double& pz, double& E) {
    double p = rng_.uniform() * pmax;
    while (rng_.uniform() * mx > dist(p, m)) {
      p = rng_.uniform() * pmax;
    }
    const double u   = rng_.uniform()*2.0 - 1.0;
    const double phi = rng_.uniform() * 2.0 * PI;
    const double s   = std::sqrt(std::max(0.0, 1.0 - u*u));
    px = p * s * std::cos(phi);
    py = p * s * std::sin(phi);
    pz = p * u;
    E  = std::sqrt(p*p + m*m);
  };

  // sample two light quarks
  double l1px,l1py,l1pz,l1E;
  double l2px,l2py,l2pz,l2E;
  sample_light(m1, mxvalue1, l1px, l1py, l1pz, l1E);
  sample_light(m2, mxvalue2, l2px, l2py, l2pz, l2E);

  // CM beta in LAB
  const double totE = HQ.E + l1E + l2E;
  const double cmbx = (HQ.px + l1px + l2px) / totE;
  const double cmby = (HQ.py + l1py + l2py) / totE;
  const double cmbz = (HQ.pz + l1pz + l2pz) / totE;

  // boost to CM
  double hpx=HQ.px, hpy=HQ.py, hpz=HQ.pz, hE=HQ.E;
  Lorentz::boost(-cmbx,-cmby,-cmbz, hpx,hpy,hpz,hE);
  Lorentz::boost(-cmbx,-cmby,-cmbz, l1px,l1py,l1pz,l1E);
  Lorentz::boost(-cmbx,-cmby,-cmbz, l2px,l2py,l2pz,l2E);

  // Jacobi momenta (three configs)
  const double q11x = (l2E*l1px - l1E*l2px) / (l1E + l2E);
  const double q11y = (l2E*l1py - l1E*l2py) / (l1E + l2E);
  const double q11z = (l2E*l1pz - l1E*l2pz) / (l1E + l2E);
  const double q11_2 = q11x*q11x + q11y*q11y + q11z*q11z;

  const double q12x = (hE*(l1px+l2px) - (l1E+l2E)*hpx) / (l1E + l2E + hE);
  const double q12y = (hE*(l1py+l2py) - (l1E+l2E)*hpy) / (l1E + l2E + hE);
  const double q12z = (hE*(l1pz+l2pz) - (l1E+l2E)*hpz) / (l1E + l2E + hE);
  const double q12_2 = q12x*q12x + q12y*q12y + q12z*q12z;

  const double q21x = (hE*l1px - l1E*hpx) / (l1E + hE);
  const double q21y = (hE*l1py - l1E*hpy) / (l1E + hE);
  const double q21z = (hE*l1pz - l1E*hpz) / (l1E + hE);
  const double q21_2 = q21x*q21x + q21y*q21y + q21z*q21z;

  const double q22x = (l2E*(l1px+hpx) - (l1E+hE)*l2px) / (l1E + l2E + hE);
  const double q22y = (l2E*(l1py+hpy) - (l1E+hE)*l2py) / (l1E + l2E + hE);
  const double q22z = (l2E*(l1pz+hpz) - (l1E+hE)*l2pz) / (l1E + l2E + hE);
  const double q22_2 = q22x*q22x + q22y*q22y + q22z*q22z;

  const double q31x = (hE*l2px - l2E*hpx) / (l2E + hE);
  const double q31y = (hE*l2py - l2E*hpy) / (l2E + hE);
  const double q31z = (hE*l2pz - l2E*hpz) / (l2E + hE);
  const double q31_2 = q31x*q31x + q31y*q31y + q31z*q31z;

  const double q32x = (l1E*(l2px+hpx) - (l2E+hE)*l1px) / (l1E + l2E + hE);
  const double q32y = (l1E*(l2py+hpy) - (l2E+hE)*l1py) / (l1E + l2E + hE);
  const double q32z = (l1E*(l2pz+hpz) - (l2E+hE)*l1pz) / (l1E + l2E + hE);
  const double q32_2 = q32x*q32x + q32y*q32y + q32z*q32z;

  // Wigner functions
  const double probmax = (orbital == 0) ? current_max_3body_s_ : current_max_3body_p_;

  auto wigner_s = [&](double s1, double s2, double qA2, double qB2) {
    return std::pow(s1*s2, 3.0) * std::exp(-qA2*s1*s1 - qB2*s2*s2);
  };
  auto wigner_p = [&](double s1, double s2, double qA2, double qB2) {
    const double base = std::pow(s1*s2, 3.0) * std::exp(-qA2*s1*s1 - qB2*s2*s2);
    return base * (1.0/3.0) * (qA2*s1*s1 + qB2*s2*s2);
  };

  double W1, W2, W3;
  if (orbital == 0) {
    W1 = wigner_s(sigma11, sigma12, q11_2, q12_2);
    W2 = wigner_s(sigma21, sigma22, q21_2, q22_2);
    W3 = wigner_s(sigma31, sigma32, q31_2, q32_2);
  } else {
    W1 = wigner_p(sigma11, sigma12, q11_2, q12_2);
    W2 = wigner_p(sigma21, sigma22, q21_2, q22_2);
    W3 = wigner_p(sigma31, sigma32, q31_2, q32_2);
  }

  const double common = std::pow(2.0*std::sqrt(PI), 6.0);
  const double probWigner = common * (W1 + W2 + W3) / 3.0;

  const bool pureRecomb = (cfg_.hadr_flag == 2);
  const bool accept = pureRecomb ? true : (rng_.uniform() * probmax < probWigner);
  if (accept) {
    todoReco = 2;
    directM = hE + l1E + l2E;
    vx_HM = cmbx;
    vy_HM = cmby;
    vz_HM = cmbz;
    return;
  }

  todoReco = 1;
}
