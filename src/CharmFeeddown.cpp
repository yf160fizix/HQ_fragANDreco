#include "CharmFeeddown.h"

#include <algorithm>
#include <cmath>

namespace Frag {
namespace {

constexpr double kMD0 = 1.86484;
constexpr double kMDp = 1.86966;
constexpr double kMDs = 1.96834;
constexpr double kMDst0 = 2.00685;
constexpr double kMDstp = 2.01026;
constexpr double kMLc = 2.28646;
constexpr double kMPiC = 0.13957;
constexpr double kMPi0 = 0.13498;
constexpr double kMKp = 0.493677;
constexpr double kMK0 = 0.497611;
constexpr double kMP = 0.938272;
constexpr double kMN = 0.939565;

constexpr double kBrD2ToDpi = 0.62;
constexpr double kBrDs2ToDstK = 0.10;

// The unknown inclusive Lambda_c(2860) -> DN branch is set to 50%; its DN
// component is divided equally between D0 p and D+ n.
constexpr double kBrLc2860ToDN = 1.0;//0.50; //one can also set this to be 0.5

} // namespace

bool CharmFeeddown::decay(Particle& hadron,
                          const std::vector<DecayChannel>& channels,
                          RNG& rng) const
{
  if (channels.empty()) return false;

  double BR_sum = 0.0;
  for (const DecayChannel& channel : channels) {
    BR_sum += channel.BR;
  }
  if (!(BR_sum > 0.0)) return false;

  const double deviate = rng.uniform() * BR_sum;
  double BR_cumulative = 0.0;
  DecayChannel selected_channel = channels.back();
  for (const DecayChannel& channel : channels) {
    BR_cumulative += channel.BR;
    if (deviate < BR_cumulative) {
      selected_channel = channel;
      break;
    }
  }

  const int sign = hadron.pid >= 0 ? +1 : -1;
  twoBodyDecay(hadron,
               selected_channel.m_daughter,
               selected_channel.m_companion,
               rng);
  hadron.pid = sign * std::abs(selected_channel.daughter_pid);
  hadron.m = selected_channel.m_daughter;

  for (int depth = 0; depth < 4; ++depth) {
    if (!decayKnownState(hadron, rng)) break;
  }
  return true;
}

bool CharmFeeddown::decayKnownState(Particle& hadron, RNG& rng) const
{
  const int abs_pid = std::abs(hadron.pid);
  std::vector<DecayChannel> channels;

  if (abs_pid == 423) {
    channels = {{421, kMD0, kMPi0, 1.0}};
  } else if (abs_pid == 413) {
    channels = {
        {421, kMD0, kMPiC, 0.677},
        {411, kMDp, kMPi0, 0.323},
    };
  } else if (abs_pid == 433) {
    channels = {{431, kMDs, kMPi0, 1.0}};
  } else if (abs_pid == 10423 || abs_pid == 20413 ||
             abs_pid == 10413) {
    channels = {
        {423, kMDst0, kMPiC, 0.5},
        {413, kMDstp, kMPiC, 0.5},
    };
  } else if (abs_pid == 10411) {
    channels = {
        {421, kMD0, kMPiC, 0.5},
        {411, kMDp, kMPi0, 0.5},
    };
  } else if (abs_pid == 415) {
    const double BR_Dpi = kBrD2ToDpi;
    const double BR_Dstpi = 1.0 - BR_Dpi;
    channels = {
        {421, kMD0, kMPiC, 0.5 * BR_Dpi},
        {411, kMDp, kMPi0, 0.5 * BR_Dpi},
        {423, kMDst0, kMPiC, 0.5 * BR_Dstpi},
        {413, kMDstp, kMPiC, 0.5 * BR_Dstpi},
    };
  } else if (abs_pid == 10433) {
    constexpr double BR_Dst0K = 36.0 / (36.0 + 31.0);
    constexpr double BR_DstpK = 31.0 / (36.0 + 31.0);
    channels = {
        {423, kMDst0, kMKp, BR_Dst0K},
        {413, kMDstp, kMK0, BR_DstpK},
    };
  } else if (abs_pid == 10431 || abs_pid == 20433) {
    channels = {{431, kMDs, kMPi0, 1.0}};
  } else if (abs_pid == 435) {
    const double BR_DstK = kBrDs2ToDstK;
    const double BR_DK = 1.0 - BR_DstK;
    channels = {
        {421, kMD0, kMKp, 0.5 * BR_DK},
        {411, kMDp, kMK0, 0.5 * BR_DK},
        {423, kMDst0, kMKp, 0.5 * BR_DstK},
        {413, kMDstp, kMK0, 0.5 * BR_DstK},
    };
  } else if (abs_pid == 4124) {
    const double BR_DN = kBrLc2860ToDN;
    channels = {
        {421, kMD0, kMP, 0.5 * BR_DN},
        {411, kMDp, kMN, 0.5 * BR_DN},
        {4122, kMLc, 2.0 * kMPiC, 1.0 - BR_DN},
    };
  } else if (abs_pid == 14122 || abs_pid == 4212 ||
             abs_pid == 4214 || abs_pid == 14212) {
    channels = {{4122, kMLc, kMPiC, 1.0}};
  } else {
    return false;
  }

  return decay(hadron, channels, rng);
}

void CharmFeeddown::twoBodyDecay(Particle& mother,
                                 double m_daughter,
                                 double m_companion,
                                 RNG& rng)
{
  constexpr double pi = 3.14159265358979323846;

  const double initial_px = mother.px;
  const double initial_py = mother.py;
  const double initial_pz = mother.pz;
  const double initial_E = mother.E;
  const double p_initial_sq =
      initial_px * initial_px + initial_py * initial_py + initial_pz * initial_pz;

  if (!(initial_E > 0.0)) {
    mother.m = m_daughter;
    mother.E = std::sqrt(p_initial_sq + m_daughter * m_daughter);
    return;
  }

  const double M_sq = initial_E * initial_E - p_initial_sq;
  if (!(M_sq > 0.0)) {
    mother.m = m_daughter;
    mother.E = std::sqrt(p_initial_sq + m_daughter * m_daughter);
    return;
  }

  const double M = std::sqrt(M_sq);
  if (m_companion <= 0.0 ||
      std::abs(M - m_daughter) < 1e-8) {
    mother.m = m_daughter;
    mother.E = std::sqrt(p_initial_sq + m_daughter * m_daughter);
    return;
  }

  // Kinematically closed effective channels reset the daughter mass shell.
  if (M <= m_daughter + m_companion) {
    mother.m = m_daughter;
    mother.E = std::sqrt(p_initial_sq + m_daughter * m_daughter);
    return;
  }

  const double first_kallen_factor =
      M * M -
      std::pow(m_daughter + m_companion, 2);
  const double second_kallen_factor =
      M * M -
      std::pow(m_daughter - m_companion, 2);
  const double p_star =
      std::sqrt(std::max(0.0, first_kallen_factor * second_kallen_factor)) /
      (2.0 * M);

  const double cos_theta = 2.0 * rng.uniform() - 1.0;
  const double sin_theta =
      std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
  const double phi = 2.0 * pi * rng.uniform();

  double px = p_star * sin_theta * std::cos(phi);
  double py = p_star * sin_theta * std::sin(phi);
  double pz = p_star * cos_theta;
  double E = std::sqrt(
      p_star * p_star + m_daughter * m_daughter);

  const double beta_x = initial_px / initial_E;
  const double beta_y = initial_py / initial_E;
  const double beta_z = initial_pz / initial_E;
  const double beta_squared =
      beta_x * beta_x + beta_y * beta_y + beta_z * beta_z;

  if (beta_squared > 1e-16) {
    const double gamma =
        1.0 / std::sqrt(std::max(1e-16, 1.0 - beta_squared));
    const double beta_dot_momentum =
        beta_x * px + beta_y * py + beta_z * pz;
    const double boost_factor =
        ((gamma - 1.0) * beta_dot_momentum / beta_squared) + gamma * E;

    px += boost_factor * beta_x;
    py += boost_factor * beta_y;
    pz += boost_factor * beta_z;
    E = gamma * (E + beta_dot_momentum);
  }

  mother.px = px;
  mother.py = py;
  mother.pz = pz;
  mother.E = E;
  mother.m = m_daughter;
}

} // namespace Frag
