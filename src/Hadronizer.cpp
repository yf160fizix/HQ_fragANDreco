#include "Hadronizer.h"

#include "Lorentz.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

double momentumMagnitude(const Particle& particle)
{
  return std::sqrt(std::max(0.0, particle.px * particle.px +
                                 particle.py * particle.py +
                                 particle.pz * particle.pz));
}

} // namespace

Hadronizer::Hadronizer(const HadConfig& config,
                       RNG& rng,
                       const RecombTable& table)
    : config_(config),
      rng_(rng),
      recomb_table_(table),
      frag_(config.eps_M, config.eps_B)
{
  baryon_s_wave_envelope_ = config_.default_baryon_s_wave_envelope;
  baryon_p_wave_envelope_ = config_.default_baryon_p_wave_envelope;

  if (config_.HQ_pid == 4) {
    m_HQ_ = config_.m_c;
    omega_M_ = config_.omega_cM;
    const bool needs_frag = config_.mode == HadMode::Frag ||
                            config_.mode == HadMode::FragAndRecomb;
    if (needs_frag && !frag_.loadCharmChemistry(
            config_.charm_meson_table,
            config_.charm_baryon_table,
            config_.Tchem,
            config_.gamma_s,
            config_.gamma_HB)) {
      std::cerr << "Failed to initialize charm fragmentation chemistry.\n";
      ready_ = false;
    }
  } else if (config_.HQ_pid == 5) {
    m_HQ_ = config_.m_b;
    omega_M_ = config_.omega_bM;
  } else {
    m_HQ_ = config_.m_c;
    omega_M_ = config_.omega_cM;
  }
}

void Hadronizer::process(const std::vector<Particle>& input,
                         std::vector<Particle>& output)
{
  output.clear();
  output.reserve(input.size());

  for (const Particle& particle : input) {
    if (std::abs(particle.pid) != 4 && std::abs(particle.pid) != 5) {
      output.push_back(particle);
      continue;
    }

    processHQ(particle, output);
  }
}

void Hadronizer::processHQ(const Particle& HQ,
                           std::vector<Particle>& output)
{
  ++stats_.n_HQ;

  if (config_.mode == HadMode::Frag) {
    fragOrKeep(HQ, output);
    return;
  }

  // Recombination-only mode retains underconstruction.
  if (config_.mode != HadMode::FragAndRecomb) {
    ++stats_.n_unchanged_HQ;
    appendResult(HQ, output);
    return;
  }

  Particle HQ_lrf = HQ;
  const double original_mass = HQ_lrf.m;
  if (config_.rescale_HQ_mass) {
    HQ_lrf.m = m_HQ_;
    HQ_lrf.E = std::sqrt(std::max(
        0.0, HQ_lrf.px * HQ_lrf.px +
             HQ_lrf.py * HQ_lrf.py +
             HQ_lrf.pz * HQ_lrf.pz +
             HQ_lrf.m * HQ_lrf.m));
  }

  Lorentz::boost(-HQ_lrf.cvx,
                 -HQ_lrf.cvy,
                 -HQ_lrf.cvz,
                 HQ_lrf.px,
                 HQ_lrf.py,
                 HQ_lrf.pz,
                 HQ_lrf.E);

  std::vector<double> channel_cdf;
  if (!interpolateChannelCDF(HQ_lrf, channel_cdf) ||
      channel_cdf.size() <
          static_cast<std::size_t>(config_.n_recomb_channels)) {
    restoreLabState(HQ_lrf, original_mass);
    fragOrKeep(HQ_lrf, output);
    return;
  }

  const double recomb_draw = rng_.uniform();
  if (recomb_draw >= channel_cdf.back()) {
    restoreLabState(HQ_lrf, original_mass);
    fragOrKeep(HQ_lrf, output);
    return;
  }

  const auto channel = drawChannel(channel_cdf, recomb_draw);
  Particle hadron;
  const ChannelOutcome outcome = channel
      ? executeChannel(*channel, HQ_lrf, hadron)
      : ChannelOutcome::Failed;

  if (outcome == ChannelOutcome::Dropped) {
    ++stats_.n_dropped;
    return;
  }

  if (outcome == ChannelOutcome::Recomb) {
    hadron.origin = HadOrigin::Recomb;
    if (config_.enable_feeddown) {
      frag_.feedDownByPDG(hadron, rng_);
    }
    ++stats_.n_recomb;
    appendResult(hadron, output);
    return;
  }

  restoreLabState(HQ_lrf, original_mass);
  fragOrKeep(HQ_lrf, output);
}

Hadronizer::ChannelOutcome Hadronizer::executeChannel(
    std::size_t channel_index,
    const Particle& HQ,
    Particle& hadron)
{
  const RecombChannel* channel =
      findRecombChannel(config_.HQ_pid, channel_index);
  if (channel == nullptr) return ChannelOutcome::Failed;
  if (channel->constituents == ConstituentSystem::Unsupported) {
    return ChannelOutcome::Dropped;
  }

  HadronState state = channel->primary;
  if (channel->alternative.pid != 0 && rng_.uniform() >= 0.5) {
    state = channel->alternative;
  }
  const int signed_pid = HQ.pid < 0
      ? -std::abs(state.pid)
      : std::abs(state.pid);

  bool did_recomb = false;
  if (channel->constituents == ConstituentSystem::Meson) {
    did_recomb = tryRecombMeson(HQ,
                                channel->light_flavor1,
                                channel->orbital,
                                signed_pid,
                                state.mass,
                                hadron);
  } else {
    did_recomb = tryRecombBaryon(HQ,
                                 channel->light_flavor1,
                                 channel->light_flavor2,
                                 channel->orbital,
                                 signed_pid,
                                 state.mass,
                                 hadron);
  }
  return did_recomb ? ChannelOutcome::Recomb : ChannelOutcome::Failed;
}

bool Hadronizer::interpolateChannelCDF(const Particle& HQ,
                                       std::vector<double>& cdf)
{
  const std::size_t channel_count = recomb_table_.channelCount();
  if (channel_count == 0 || recomb_table_.momentum_grid.size() < 2 ||
      recomb_table_.channel_cdf.size() != recomb_table_.momentum_grid.size()) {
    return false;
  }

  const double spacing = config_.HQ_pid == 4
      ? config_.dp_c
      : config_.dp_b;
  if (!(spacing > 0.0)) return false;

  const double momentum = momentumMagnitude(HQ);
  const int lower_index = static_cast<int>(momentum / spacing);
  const int upper_index = lower_index + 1;
  const int row_count =
      static_cast<int>(recomb_table_.momentum_grid.size());
  if (lower_index < 0 || upper_index >= row_count) return false;

  const auto lower = static_cast<std::size_t>(lower_index);
  const auto upper = static_cast<std::size_t>(upper_index);
  const double lower_momentum = lower_index * spacing;
  const double fraction = (momentum - lower_momentum) / spacing;
  const double interpolation = std::clamp(fraction, 0.0, 1.0);

  cdf.assign(channel_count, 0.0);
  for (std::size_t channel = 0; channel < channel_count; ++channel) {
    const double lower_probability =
        recomb_table_.channel_cdf[lower][channel];
    const double upper_probability =
        recomb_table_.channel_cdf[upper][channel];
    cdf[channel] = lower_probability +
                   (upper_probability - lower_probability) * interpolation;
  }

  double previous = 0.0;
  for (double& probability : cdf) {
    probability = std::clamp(probability, previous, 1.0);
    previous = probability;
  }

  if (lower < recomb_table_.wigner_maxima.size() &&
      recomb_table_.wigner_maxima[lower].size() >
          RecombTable::baryon_p_wave_column) {
    baryon_s_wave_envelope_ = recomb_table_.wigner_maxima[lower]
        [RecombTable::baryon_s_wave_column];
    baryon_p_wave_envelope_ = recomb_table_.wigner_maxima[lower]
        [RecombTable::baryon_p_wave_column];
  } else {
    baryon_s_wave_envelope_ = config_.default_baryon_s_wave_envelope;
    baryon_p_wave_envelope_ = config_.default_baryon_p_wave_envelope;
  }
  return true;
}

std::optional<std::size_t> Hadronizer::drawChannel(
    const std::vector<double>& cdf,
    double uniform_deviate)
{
  for (std::size_t channel = 0; channel < cdf.size(); ++channel) {
    if (uniform_deviate < cdf[channel]) return channel;
  }
  return std::nullopt;
}

void Hadronizer::fragOrKeep(const Particle& HQ,
                            std::vector<Particle>& output)
{
  Particle hadron;
  if (frag_.fragment(HQ, rng_, hadron)) {
    ++stats_.n_frag;
    appendResult(hadron, output);
  } else {
    ++stats_.n_unchanged_HQ;
    appendResult(HQ, output);
  }
}

void Hadronizer::restoreLabState(Particle& HQ,
                                 double original_mass) const
{
  Lorentz::boost(HQ.cvx,
                 HQ.cvy,
                 HQ.cvz,
                 HQ.px,
                 HQ.py,
                 HQ.pz,
                 HQ.E);
  if (config_.rescale_HQ_mass) {
    HQ.m = original_mass;
    HQ.E = std::sqrt(std::max(
        0.0, HQ.px * HQ.px +
             HQ.py * HQ.py +
             HQ.pz * HQ.pz +
             HQ.m * HQ.m));
  }
}

void Hadronizer::appendResult(const Particle& particle,
                              std::vector<Particle>& output)
{
  ++stats_.final_species[std::abs(particle.pid)];
  output.push_back(particle);
}

Particle Hadronizer::makeHadron(const Particle& HQ,
                                int hadron_pid,
                                double hadron_mass)
{
  Particle hadron;
  hadron.pid = hadron_pid;
  hadron.m = hadron_mass;
  hadron.x = HQ.x;
  hadron.y = HQ.y;
  hadron.z = HQ.z;
  hadron.t = HQ.t;
  hadron.Thydro = HQ.Thydro;
  hadron.cvx = HQ.cvx;
  hadron.cvy = HQ.cvy;
  hadron.cvz = HQ.cvz;
  hadron.ipx = HQ.ipx;
  hadron.ipy = HQ.ipy;
  hadron.ipz = HQ.ipz;
  hadron.iE = HQ.iE;
  hadron.wt = HQ.wt;
  return hadron;
}

void Hadronizer::adjustMassShellByTwoBodyDecay(double& px,
                                                double& py,
                                                double& pz,
                                                double& E,
                                                double final_mass)
{
  constexpr double pion_mass = 0.13;
  constexpr double pi = 3.14159265358979323846;

  const double beta_x = px / E;
  const double beta_y = py / E;
  const double beta_z = pz / E;
  Lorentz::boost(-beta_x, -beta_y, -beta_z, px, py, pz, E);

  const double initial_mass = E;
  double companion_momentum = 0.0;
  if (std::abs(final_mass - initial_mass) > pion_mass) {
    const double E_comp =
        (initial_mass * initial_mass + pion_mass * pion_mass -
         final_mass * final_mass) /
        (2.0 * initial_mass);
    companion_momentum = std::sqrt(std::max(
        0.0, E_comp * E_comp - pion_mass * pion_mass));
  } else {
    companion_momentum =
        std::abs(initial_mass * initial_mass - final_mass * final_mass) /
        (2.0 * initial_mass);
  }

  const double cos_theta = rng_.uniform() * 2.0 - 1.0;
  const double phi = rng_.uniform() * 2.0 * pi;
  const double sin_theta =
      std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
  px = -companion_momentum * sin_theta * std::cos(phi);
  py = -companion_momentum * sin_theta * std::sin(phi);
  pz = -companion_momentum * cos_theta;
  E = std::sqrt(px * px + py * py + pz * pz + final_mass * final_mass);

  Lorentz::boost(beta_x, beta_y, beta_z, px, py, pz, E);
}
