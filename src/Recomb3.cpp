#include "Hadronizer.h"

#include "Lorentz.h"

#include <cmath>
#include <iostream>

namespace {

double lightQuarkMass(LightFlavor flavor,
                      const HadConfig& config)
{
  return flavor == LightFlavor::UpDown
      ? config.m_q
      : config.m_s;
}

} // namespace

bool Hadronizer::tryRecombBaryon(const Particle& HQ,
                                 LightFlavor light_flavor1,
                                 LightFlavor light_flavor2,
                                 OrbitalState orbital,
                                 int hadron_pid,
                                 double hadron_mass,
                                 Particle& hadron)
{
  const double T = HQ.Thydro;
  if (T < 0.1) return false;

  const double m_light1 = lightQuarkMass(light_flavor1, config_);
  const double m_light2 = lightQuarkMass(light_flavor2, config_);
  const double omega = std::abs(HQ.pid) == 4
      ? config_.omega_cB
      : config_.omega_bB;

  const BaryonWignerWidths widths = computeBaryonWignerWidths(
      HQ.m, m_light1, m_light2, omega);
  const auto thermal_max_1 = thermalMax(light_flavor1, T);
  const auto thermal_max_2 = thermalMax(light_flavor2, T);
  if (!thermal_max_1 || !thermal_max_2) return false;

  RecombKinematics kinematics;
  long long attempt_count = 0;
  while (!sampleBaryonKinematics(orbital,
                                 HQ,
                                 m_light1,
                                 m_light2,
                                 *thermal_max_1,
                                 *thermal_max_2,
                                 widths,
                                 kinematics)) {
    if (++attempt_count > 100000000LL) {
      std::cerr << "Baryon recombination exceeded 100000000 attempts.\n";
      return false;
    }
  }

  buildRecombHadron(HQ, hadron_pid, hadron_mass, kinematics, hadron);
  return true;
}

bool Hadronizer::sampleBaryonKinematics(
    OrbitalState orbital,
    const Particle& HQ,
    double m_light1,
    double m_light2,
    double thermal_max_1,
    double thermal_max_2,
    const BaryonWignerWidths& widths,
    RecombKinematics& result)
{
  const double T = HQ.Thydro;
  FourMomentum p_light1 = sampleThermalParton(
      T, m_light1, thermal_max_1);
  FourMomentum p_light2 = sampleThermalParton(
      T, m_light2, thermal_max_2);
  FourMomentum p_HQ{HQ.px, HQ.py, HQ.pz, HQ.E};

  const double E_tot =
      HQ.E + p_light1.E + p_light2.E;
  const double cm_beta_x =
      (HQ.px + p_light1.px + p_light2.px) / E_tot;
  const double cm_beta_y =
      (HQ.py + p_light1.py + p_light2.py) / E_tot;
  const double cm_beta_z =
      (HQ.pz + p_light1.pz + p_light2.pz) / E_tot;

  Lorentz::boost(-cm_beta_x, -cm_beta_y, -cm_beta_z, p_HQ);
  Lorentz::boost(-cm_beta_x, -cm_beta_y, -cm_beta_z, p_light1);
  Lorentz::boost(-cm_beta_x, -cm_beta_y, -cm_beta_z, p_light2);

  const BaryonJacobiMomenta jacobi_momenta =
      computeJacobiMomenta(p_HQ, p_light1, p_light2);

  const double wigner_envelope = orbital == OrbitalState::S
      ? baryon_s_wave_envelope_
      : baryon_p_wave_envelope_;
  const double wigner_probability =
      baryonWignerWeight(orbital, widths, jacobi_momenta);
  const bool recomb_only = config_.mode == HadMode::Recomb;
  if (!recomb_only &&
      rng_.uniform() * wigner_envelope >= wigner_probability) {
    return false;
  }

  result.M =
      p_HQ.E + p_light1.E + p_light2.E;
  result.beta_x = cm_beta_x;
  result.beta_y = cm_beta_y;
  result.beta_z = cm_beta_z;
  return true;
}

Hadronizer::BaryonJacobiMomenta Hadronizer::computeJacobiMomenta(
    const FourMomentum& p_HQ,
    const FourMomentum& p_light1,
    const FourMomentum& p_light2)
{
  const double q11x =
      (p_light2.E * p_light1.px - p_light1.E * p_light2.px) /
      (p_light1.E + p_light2.E);
  const double q11y =
      (p_light2.E * p_light1.py - p_light1.E * p_light2.py) /
      (p_light1.E + p_light2.E);
  const double q11z =
      (p_light2.E * p_light1.pz - p_light1.E * p_light2.pz) /
      (p_light1.E + p_light2.E);
  const double q11_sq = q11x * q11x + q11y * q11y + q11z * q11z;

  const double q12x =
      (p_HQ.E * (p_light1.px + p_light2.px) -
       (p_light1.E + p_light2.E) * p_HQ.px) /
      (p_light1.E + p_light2.E + p_HQ.E);
  const double q12y =
      (p_HQ.E * (p_light1.py + p_light2.py) -
       (p_light1.E + p_light2.E) * p_HQ.py) /
      (p_light1.E + p_light2.E + p_HQ.E);
  const double q12z =
      (p_HQ.E * (p_light1.pz + p_light2.pz) -
       (p_light1.E + p_light2.E) * p_HQ.pz) /
      (p_light1.E + p_light2.E + p_HQ.E);
  const double q12_sq = q12x * q12x + q12y * q12y + q12z * q12z;

  const double q21x =
      (p_HQ.E * p_light1.px - p_light1.E * p_HQ.px) /
      (p_light1.E + p_HQ.E);
  const double q21y =
      (p_HQ.E * p_light1.py - p_light1.E * p_HQ.py) /
      (p_light1.E + p_HQ.E);
  const double q21z =
      (p_HQ.E * p_light1.pz - p_light1.E * p_HQ.pz) /
      (p_light1.E + p_HQ.E);
  const double q21_sq = q21x * q21x + q21y * q21y + q21z * q21z;

  const double q22x =
      (p_light2.E * (p_light1.px + p_HQ.px) -
       (p_light1.E + p_HQ.E) * p_light2.px) /
      (p_light1.E + p_light2.E + p_HQ.E);
  const double q22y =
      (p_light2.E * (p_light1.py + p_HQ.py) -
       (p_light1.E + p_HQ.E) * p_light2.py) /
      (p_light1.E + p_light2.E + p_HQ.E);
  const double q22z =
      (p_light2.E * (p_light1.pz + p_HQ.pz) -
       (p_light1.E + p_HQ.E) * p_light2.pz) /
      (p_light1.E + p_light2.E + p_HQ.E);
  const double q22_sq = q22x * q22x + q22y * q22y + q22z * q22z;

  const double q31x =
      (p_HQ.E * p_light2.px - p_light2.E * p_HQ.px) /
      (p_light2.E + p_HQ.E);
  const double q31y =
      (p_HQ.E * p_light2.py - p_light2.E * p_HQ.py) /
      (p_light2.E + p_HQ.E);
  const double q31z =
      (p_HQ.E * p_light2.pz - p_light2.E * p_HQ.pz) /
      (p_light2.E + p_HQ.E);
  const double q31_sq = q31x * q31x + q31y * q31y + q31z * q31z;

  const double q32x =
      (p_light1.E * (p_light2.px + p_HQ.px) -
       (p_light2.E + p_HQ.E) * p_light1.px) /
      (p_light1.E + p_light2.E + p_HQ.E);
  const double q32y =
      (p_light1.E * (p_light2.py + p_HQ.py) -
       (p_light2.E + p_HQ.E) * p_light1.py) /
      (p_light1.E + p_light2.E + p_HQ.E);
  const double q32z =
      (p_light1.E * (p_light2.pz + p_HQ.pz) -
       (p_light2.E + p_HQ.E) * p_light1.pz) /
      (p_light1.E + p_light2.E + p_HQ.E);
  const double q32_sq = q32x * q32x + q32y * q32y + q32z * q32z;

  return {q11_sq, q12_sq, q21_sq, q22_sq, q31_sq, q32_sq};
}

double Hadronizer::baryonWignerWeight(
    OrbitalState orbital,
    const BaryonWignerWidths& widths,
    const BaryonJacobiMomenta& momenta)
{
  const auto s_wave_wigner = [](double sigma_1,
                                double sigma_2,
                                double q1_sq,
                                double q2_sq) {
    return std::pow(sigma_1 * sigma_2, 3.0) *
           std::exp(-q1_sq * sigma_1 * sigma_1 -
                    q2_sq * sigma_2 * sigma_2);
  };
  const auto p_wave_wigner = [](double sigma_1,
                                double sigma_2,
                                double q1_sq,
                                double q2_sq) {
    const double base =
        std::pow(sigma_1 * sigma_2, 3.0) *
        std::exp(-q1_sq * sigma_1 * sigma_1 -
                 q2_sq * sigma_2 * sigma_2);
    return base * (1.0 / 3.0) *
           (q1_sq * sigma_1 * sigma_1 +
            q2_sq * sigma_2 * sigma_2);
  };

  double W1 = 0.0;
  double W2 = 0.0;
  double W3 = 0.0;
  if (orbital == OrbitalState::S) {
    W1 = s_wave_wigner(
        widths.sigma_11, widths.sigma_12,
        momenta.q11_sq, momenta.q12_sq);
    W2 = s_wave_wigner(
        widths.sigma_21, widths.sigma_22,
        momenta.q21_sq, momenta.q22_sq);
    W3 = s_wave_wigner(
        widths.sigma_31, widths.sigma_32,
        momenta.q31_sq, momenta.q32_sq);
  } else {
    W1 = p_wave_wigner(
        widths.sigma_11, widths.sigma_12,
        momenta.q11_sq, momenta.q12_sq);
    W2 = p_wave_wigner(
        widths.sigma_21, widths.sigma_22,
        momenta.q21_sq, momenta.q22_sq);
    W3 = p_wave_wigner(
        widths.sigma_31, widths.sigma_32,
        momenta.q31_sq, momenta.q32_sq);
  }

  constexpr double pi = 3.14159265358979323846;
  const double normalization = std::pow(2.0 * std::sqrt(pi), 6.0);
  return normalization * (W1 + W2 + W3) / 3.0;
}

Hadronizer::BaryonWignerWidths Hadronizer::computeBaryonWignerWidths(
    double m_HQ,
    double m_light1,
    double m_light2,
    double omega)
{
  BaryonWignerWidths widths;
  double mu_1 = m_light1 * m_light2 / (m_light1 + m_light2);
  double mu_2 =
      (m_light1 + m_light2) * m_HQ / (m_light1 + m_light2 + m_HQ);
  widths.sigma_11 = 1.0 / std::sqrt(mu_1 * omega);
  widths.sigma_12 = 1.0 / std::sqrt(mu_2 * omega);

  mu_1 = m_light1 * m_HQ / (m_light1 + m_HQ);
  mu_2 =
      (m_light1 + m_HQ) * m_light2 / (m_light1 + m_light2 + m_HQ);
  widths.sigma_21 = 1.0 / std::sqrt(mu_1 * omega);
  widths.sigma_22 = 1.0 / std::sqrt(mu_2 * omega);

  mu_1 = m_light2 * m_HQ / (m_light2 + m_HQ);
  mu_2 =
      (m_light2 + m_HQ) * m_light1 / (m_light1 + m_light2 + m_HQ);
  widths.sigma_31 = 1.0 / std::sqrt(mu_1 * omega);
  widths.sigma_32 = 1.0 / std::sqrt(mu_2 * omega);
  return widths;
}
