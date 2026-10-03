#include "Hadronizer.h"

#include "Lorentz.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace {

constexpr int kMaxEnvelopeIndex = 70;

const std::array<double, kMaxEnvelopeIndex + 1> max_ud = {{
  0.063, 0.064, 0.065, 0.067, 0.068, 0.069,
  0.071, 0.072, 0.073, 0.075, 0.076, 0.077,
  0.079, 0.080, 0.081, 0.083, 0.084, 0.086,
  0.087, 0.088, 0.090, 0.091, 0.093, 0.094,
  0.096, 0.097, 0.099, 0.100, 0.102, 0.103,
  0.105, 0.107, 0.108, 0.110, 0.111, 0.113,
  0.115, 0.116, 0.118, 0.120, 0.121, 0.123,
  0.125, 0.126, 0.128, 0.130, 0.132, 0.133,
  0.135, 0.137, 0.139, 0.141, 0.142, 0.144,
  0.146, 0.148, 0.150, 0.152, 0.154, 0.155,
  0.157, 0.159, 0.161, 0.163, 0.165, 0.167,
  0.169, 0.171, 0.173, 0.175, 0.177
}};

const std::array<double, kMaxEnvelopeIndex + 1> max_s = max_ud;

} // namespace

std::optional<double> Hadronizer::thermalMax(
    LightFlavor flavor,
    double T) const
{
  int iT =
      static_cast<int>((T - 0.160) / 0.002 + 0.5);
  if (iT < 0) return std::nullopt;
  if (iT > 35) iT = 35;

  const auto envelope_index = static_cast<std::size_t>(
      std::clamp(iT * 2, 0, kMaxEnvelopeIndex));
  return flavor == LightFlavor::UpDown
      ? max_ud[envelope_index]
      : max_s[envelope_index];
}

FourMomentum Hadronizer::sampleThermalParton(
    double T,
    double mass,
    double thermal_max)
{
  const double p_max = 15.0 * T;
  const auto thermal_distribution = [&](double momentum) {
    const double E_light =
        std::sqrt(momentum * momentum + mass * mass);
    const double E_g = std::sqrt(
        momentum * momentum + config_.m_g * config_.m_g);
    const double quark_term =
        6.0 * momentum * momentum /
        (std::exp(E_light / T) + 1.0);
    const double gluon_term =
        (momentum * momentum /
         (std::exp(E_g / T) - 1.0)) *
        (16.0 / 6.0);
    return quark_term + gluon_term;
  };

  double momentum = rng_.uniform() * p_max;
  while (rng_.uniform() * thermal_max >
         thermal_distribution(momentum)) {
    momentum = rng_.uniform() * p_max;
  }

  constexpr double pi = 3.14159265358979323846;
  const double cos_theta = rng_.uniform() * 2.0 - 1.0;
  const double phi = rng_.uniform() * 2.0 * pi;
  const double sin_theta =
      std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));

  return {
      momentum * sin_theta * std::cos(phi),
      momentum * sin_theta * std::sin(phi),
      momentum * cos_theta,
      std::sqrt(momentum * momentum + mass * mass),
  };
}

void Hadronizer::buildRecombHadron(
    const Particle& HQ,
    int hadron_pid,
    double hadron_mass,
    const RecombKinematics& kinematics,
    Particle& hadron)
{
  hadron = makeHadron(HQ, hadron_pid, hadron_mass);
  hadron.E = kinematics.M;
  Lorentz::boost(kinematics.beta_x,
                 kinematics.beta_y,
                 kinematics.beta_z,
                 hadron.px,
                 hadron.py,
                 hadron.pz,
                 hadron.E);
  Lorentz::boost(HQ.cvx,
                 HQ.cvy,
                 HQ.cvz,
                 hadron.px,
                 hadron.py,
                 hadron.pz,
                 hadron.E);
  adjustMassShellByTwoBodyDecay(
      hadron.px, hadron.py, hadron.pz, hadron.E, hadron.m);
}
