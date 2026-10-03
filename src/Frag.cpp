#include "Frag.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>

namespace Frag {

PetersonFrag::PetersonFrag(double eps_M, double eps_B)
    : eps_M_(eps_M), eps_B_(eps_B)
{
  if (!(eps_M > 0.0) || !std::isfinite(eps_M) ||
      !(eps_B > 0.0) || !std::isfinite(eps_B)) {
    std::cerr << "PetersonFragmentation: both epsilon values must be finite "
                 "and positive.\n";
  }
}

void PetersonFrag::setChemistry(const std::vector<CharmHadronState>& table)
{
  species_.clear();
  species_cdf_.clear();
  species_epsilon_.clear();
  kernel_maxima_.clear();
  ready_ = false;

  double sum = 0.0;

  for (const CharmHadronState& species : table) {
    if (species.pid == 0) continue;
    if (!(species.mass > 0.0)) continue;
    if (!(species.statistical_weight > 0.0)) continue;
    if (species.decay_channels.empty()) continue;

    species_.push_back(species);
    sum += species.statistical_weight;
  }

  if (species_.empty() || !(sum > 0.0)) {
    std::cerr << "PetersonFragmentation::setChemistry: empty or invalid charm-state table.\n";
    return;
  }

  species_cdf_.reserve(species_.size());
  species_epsilon_.reserve(species_.size());
  kernel_maxima_.reserve(species_.size());

  double cumulative_weight = 0.0;

  for (CharmHadronState& species : species_) {
    species.statistical_weight /= sum;
    cumulative_weight += species.statistical_weight;
    species_cdf_.push_back(cumulative_weight);

    const double eps = epsForSpecies(species);
    species_epsilon_.push_back(eps);
    kernel_maxima_.push_back(estimateKernelMax(eps));
  }

  species_cdf_.back() = 1.0;

  ready_ = (species_epsilon_.size() == species_.size())
        && (kernel_maxima_.size() == species_.size())
        && std::all_of(kernel_maxima_.begin(),
                       kernel_maxima_.end(),
                       [](double maximum) { return maximum > 0.0; });

}

bool PetersonFrag::loadCharmChemistry(
    const std::string& meson_path,
    const std::string& baryon_path,
    double Tchem,
    double gamma_s,
    double gamma_HB)
{
  std::vector<CharmHadronState> table;
  try {
    table = loadCharmStates(
        meson_path,
        baryon_path,
        Tchem,
        gamma_s,
        gamma_HB
    );
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return false;
  }

  setChemistry(table);
  return ready_;
}

// Unnormalized Peterson fragmentation kernel,
//   D(z) \propto 1 / { z [ 1 - 1/z - epsilon/(1-z) ]^2 }.
double PetersonFrag::petersonKernel(double z, double epsilon)
{
  if (z <= 0.0 || z >= 1.0) return 0.0;
  if (!(epsilon > 0.0)) return 0.0;

  const double omz = 1.0 - z;
  if (!(omz > 0.0)) return 0.0;

  const double denom = 1.0 - 1.0 / z - epsilon / omz;
  const double val = 1.0 / (z * denom * denom);

  if (!std::isfinite(val) || val <= 0.0) return 0.0;
  return val;
}

double PetersonFrag::estimateKernelMax(
    double epsilon,
    double zmin,
    double zmax)
{
  zmin = std::max(1e-8, zmin);
  zmax = std::min(1.0 - 1e-8, zmax);

  if (!(zmin < zmax)) return 0.0;

  double ymax = 0.0;

  for (int i = 0; i < 2000; ++i) {
    const double t = (i + 0.5) / 2000.0;
    const double z = zmin + (zmax - zmin) * t;
    ymax = std::max(ymax, petersonKernel(z, epsilon));
  }

  // Add a tiny safety margin for the accept-reject envelope.
  return 1.05 * ymax;
}

double PetersonFrag::epsForSpecies(
    const CharmHadronState& species) const
{
  return species.baryon ? eps_B_ : eps_M_;
}

double PetersonFrag::sampleZ(RNG& rng,
                             double epsilon,
                             double kernel_max,
                             double zmin,
                             double zmax) const
{
  zmin = std::max(1e-8, zmin);
  zmax = std::min(1.0 - 1e-8, zmax);

  if (!(zmin < zmax)) return 1.0;
  if (!(epsilon > 0.0)) return 1.0;
  if (!(kernel_max > 0.0)) return 1.0;

  for (int tries = 0; tries < 300000; ++tries) {
    const double z = zmin + (zmax - zmin) * rng.uniform();
    const double y = kernel_max * rng.uniform();

    if (y < petersonKernel(z, epsilon)) return z;
  }

  std::cerr << "PetersonFragmentation::sampleZ: rejection sampler failed; returning z=1.\n";
  return 1.0;
}

int PetersonFrag::sampleSpeciesIndex(RNG& rng) const
{
  const double r = rng.uniform();

  auto it = std::lower_bound(species_cdf_.begin(), species_cdf_.end(), r);

  if (it == species_cdf_.end()) {
    return static_cast<int>(species_cdf_.size()) - 1;
  }

  return static_cast<int>(it - species_cdf_.begin());
}

int PetersonFrag::applyHQSign(int positive_pid, int HQ_pid)
{
  return HQ_pid >= 0
      ? std::abs(positive_pid)
      : -std::abs(positive_pid);
}

bool PetersonFrag::fragment(const Particle& HQ,
                            RNG& rng,
                            Particle& hadron) const
{
  if (!ready_) return false;

  if (std::abs(HQ.pid) != 4) return false;

  const int species_index = sampleSpeciesIndex(rng);
  const CharmHadronState& species =
      species_[static_cast<std::size_t>(species_index)];

  const double eps =
      (species_index >= 0 &&
       species_index < static_cast<int>(species_epsilon_.size()))
      ? species_epsilon_[static_cast<std::size_t>(species_index)]
      : epsForSpecies(species);

  const double kernel_max =
      (species_index >= 0 &&
       species_index < static_cast<int>(kernel_maxima_.size()))
      ? kernel_maxima_[static_cast<std::size_t>(species_index)]
      : estimateKernelMax(eps);

  const double z = sampleZ(rng, eps, kernel_max);

  hadron = HQ;
  hadron.pid =
      applyHQSign(std::abs(species.pid), HQ.pid);
  hadron.m = species.mass;

  hadron.px = z * HQ.px;
  hadron.py = z * HQ.py;
  hadron.pz = z * HQ.pz;

  const double momentum_squared =
      hadron.px * hadron.px + hadron.py * hadron.py + hadron.pz * hadron.pz;

  hadron.E =
      std::sqrt(momentum_squared + species.mass * species.mass);

  hadron.origin = HadOrigin::Frag;
  feeddown_.decay(hadron, species.decay_channels, rng);

  return true;
}

bool PetersonFrag::feedDownByPDG(Particle& hadron, RNG& rng) const
{
  return feeddown_.decayKnownState(hadron, rng);
}

} // namespace Frag
