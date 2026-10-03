#include "Hadronizer.h"

#include "Lorentz.h"

#include <algorithm>
#include <cmath>
#include <iostream>

bool Hadronizer::tryRecombMeson(const Particle& HQ,
                                LightFlavor light_flavor,
                                OrbitalState orbital,
                                int hadron_pid,
                                double hadron_mass,
                                Particle& hadron)
{
  const double T = HQ.Thydro;
  if (T < 0.1) return false;

  const double m_light = light_flavor == LightFlavor::UpDown
      ? config_.m_q
      : config_.m_s;

  const double p_light_max = 15.0 * T;
  const double E_light_max = std::sqrt(
      p_light_max * p_light_max + m_light * m_light);

  const double cell_beta = std::sqrt(
      HQ.cvx * HQ.cvx + HQ.cvy * HQ.cvy + HQ.cvz * HQ.cvz);
  const double cell_gamma =
      1.0 / std::sqrt(std::max(1e-15, 1.0 - cell_beta * cell_beta));

  const double p_light_max_cm =
      cell_gamma * (p_light_max + cell_beta * E_light_max);
  const double E_light_max_cm =
      cell_gamma * (E_light_max + cell_beta * p_light_max);

  const double p_HQ = std::sqrt(std::max(
      0.0, HQ.E * HQ.E - HQ.m * HQ.m));
  const double beta_cm =
      (p_HQ + p_light_max_cm) / (HQ.E + E_light_max_cm);
  const double gamma_cm =
      1.0 / std::sqrt(std::max(1e-15, 1.0 - beta_cm * beta_cm));

  const double q_rel_min = gamma_cm * (p_HQ - beta_cm * HQ.E);
  double q_rel_sq_min = q_rel_min * q_rel_min;
  if (p_HQ < p_light_max_cm) {
    q_rel_sq_min = -1.0;
  }

  const double mu = m_light * HQ.m / (m_light + HQ.m);
  const double sigma = 1.0 / std::sqrt(mu * omega_M_);

  const auto thermal_max = thermalMax(light_flavor, T);
  if (!thermal_max) return false;

  RecombKinematics kinematics;
  long long attempt_count = 0;
  while (!sampleMesonKinematics(orbital,
                                HQ,
                                m_light,
                                sigma,
                                q_rel_sq_min,
                                *thermal_max,
                                kinematics)) {
    if (++attempt_count > 100000000LL) {
      std::cerr << "Meson recombination exceeded 100000000 attempts.\n";
      return false;
    }
  }

  buildRecombHadron(HQ, hadron_pid, hadron_mass, kinematics, hadron);
  return true;
}

bool Hadronizer::sampleMesonKinematics(
    OrbitalState orbital,
    const Particle& HQ,
    double m_light,
    double sigma,
    double q_rel_sq_min,
    double thermal_max,
    RecombKinematics& result)
{
  const double T = HQ.Thydro;
  double wigner_envelope = 1.0;
  if (orbital == OrbitalState::S) {
    wigner_envelope = q_rel_sq_min < 0.0
        ? 1.0
        : std::exp(-q_rel_sq_min * sigma * sigma);
  } else {
    const double x = q_rel_sq_min * sigma * sigma;
    wigner_envelope = x < 1.0 ? std::exp(-1.0) : x * std::exp(-x);
  }

  FourMomentum p_light = sampleThermalParton(
      T, m_light, thermal_max);

  const double cm_beta_x =
      (HQ.px + p_light.px) / (HQ.E + p_light.E);
  const double cm_beta_y =
      (HQ.py + p_light.py) / (HQ.E + p_light.E);
  const double cm_beta_z =
      (HQ.pz + p_light.pz) / (HQ.E + p_light.E);

  double p_HQ_x = HQ.px;
  double p_HQ_y = HQ.py;
  double p_HQ_z = HQ.pz;
  double E_HQ = HQ.E;
  Lorentz::boost(-cm_beta_x, -cm_beta_y, -cm_beta_z, p_light);
  Lorentz::boost(-cm_beta_x, -cm_beta_y, -cm_beta_z,
                 p_HQ_x, p_HQ_y, p_HQ_z, E_HQ);

  const double E_tot = p_light.E + E_HQ;
  if (E_tot <= 0.0) return false;

  const double qx =
      (p_light.E * p_HQ_x - E_HQ * p_light.px) /
      E_tot;
  const double qy =
      (p_light.E * p_HQ_y - E_HQ * p_light.py) /
      E_tot;
  const double qz =
      (p_light.E * p_HQ_z - E_HQ * p_light.pz) /
      E_tot;
  const double q_rel_sq = qx * qx + qy * qy + qz * qz;

  const double gaussian = std::exp(-q_rel_sq * sigma * sigma);
  const double wigner_probability = orbital == OrbitalState::S
      ? gaussian
      : gaussian * (q_rel_sq * sigma * sigma);
  const bool recomb_only = config_.mode == HadMode::Recomb;
  if (!recomb_only &&
      rng_.uniform() * wigner_envelope >= wigner_probability) {
    return false;
  }

  result.M = E_HQ + p_light.E;
  result.beta_x = cm_beta_x;
  result.beta_y = cm_beta_y;
  result.beta_z = cm_beta_z;
  return true;
}
