#include "Frag.h"
#include "RNG.h"
#include <algorithm>
#include <cmath>

namespace Frag {
namespace {
constexpr double kMc  = 1.5;
constexpr double kMD0 = 1.86484;
constexpr double kMLc = 2.28646;
constexpr double kR_D0 = 0.10;
constexpr double kR_Lc = 0.16;
} // namespace

void HQETFrag::setChemistry(const std::vector<Species>& table)
{
  tab_.clear();
  cdf_.clear();
  rvals_.clear();
  ymax_.clear();
  ready_ = false;

  double sum = 0.0;
  for (const auto& s : table) {
    if (s.pdg == 0) continue;
    if (!(s.mass > 0.0)) continue;
    if (!(s.prob > 0.0)) continue;
    tab_.push_back(s);
    sum += s.prob;
  }
  if (tab_.empty() || !(sum > 0.0)) return;

  cdf_.reserve(tab_.size());
  rvals_.reserve(tab_.size());
  ymax_.reserve(tab_.size());
  double acc = 0.0;
  for (auto& s : tab_) {
    s.prob /= sum;
    acc += s.prob;
    cdf_.push_back(acc);

    const double rH = defaultRForSpecies(s.pdg, s.mass);
    rvals_.push_back(rH);
    ymax_.push_back(estimateYmax(rH));
  }
  cdf_.back() = 1.0;
  ready_ = true;
}

void HQETFrag::setDefaultCharmChemistry()
{
  std::vector<Species> t = {
    {421,  1.86484, 0.4175},
    {411,  1.86966, 0.1834},
    {431,  1.96834, 0.1138},
    {4122, 2.28646, 0.2379},
  };
  setChemistry(t);
}

bool HQETFrag::isBaryonPDG(int pdgAbs)
{
  return pdgAbs >= 1000;
}

double HQETFrag::defaultRForSpecies(int pdg, double mass)
{
  const int apdg = std::abs(pdg);
  if (!(mass > kMc)) return isBaryonPDG(apdg) ? kR_Lc : kR_D0;

  if (isBaryonPDG(apdg)) {
    const double ref = (kMLc - kMc) / kMLc;
    const double cur = (mass - kMc) / mass;
    return kR_Lc * (cur / ref);
  }

  const double ref = (kMD0 - kMc) / kMD0;
  const double cur = (mass - kMc) / mass;
  return kR_D0 * (cur / ref);
}

// Unnormalized HQET/Braaten-type fragmentation kernel.
double HQETFrag::hqetKernel(double z, double r)
{
  if (z <= 0.0 || z >= 1.0) return 0.0;
  if (!(r > 0.0) || !(r < 1.0)) return 0.0;

  const double omz = 1.0 - z;
  const double a = 1.0 - (1.0 - r) * z;
  if (!(a > 0.0)) return 0.0;

  const double z2 = z * z;
  const double z3 = z2 * z;
  const double z4 = z2 * z2;
  const double poly =
      6.0
    - 18.0 * (1.0 - 2.0 * r) * z
    + (21.0 - 74.0 * r + 68.0 * r * r) * z2
    - 2.0 * (1.0 - r) * (6.0 - 19.0 * r + 18.0 * r * r) * z3
    + 3.0 * (1.0 - r) * (1.0 - r)
        * (1.0 - 2.0 * r + 2.0 * r * r) * z4;

  const double val = r * z * omz * omz * poly / std::pow(a, 6);
  if (!std::isfinite(val) || val <= 0.0) return 0.0;
  return val;
}


double HQETFrag::estimateYmax(double r, double zmin, double zmax)
{
  zmin = std::max(1e-8, zmin);
  zmax = std::min(1.0 - 1e-8, zmax);
  if (!(zmin < zmax)) return 0.0;

  double ymax = 0.0;
  for (int i = 0; i < 2000; ++i) {
    const double t = (i + 0.5) / 2000.0;
    const double z = zmin + (zmax - zmin) * t;
    ymax = std::max(ymax, hqetKernel(z, r));
  }
  return ymax;
}

double HQETFrag::sampleZ(RNG& rng, double r, double ymax, double zmin, double zmax) const
{
  zmin = std::max(1e-8, zmin);
  zmax = std::min(1.0 - 1e-8, zmax);
  if (!(zmin < zmax)) return 1.0;

  if (!(ymax > 0.0)) return 1.0;

  for (int tries = 0; tries < 300000; ++tries) {
    const double z = zmin + (zmax - zmin) * rng.uniform();
    const double y = ymax * rng.uniform();
    if (y < hqetKernel(z, r)) return z;
  }

  return 1.0;
}

int HQETFrag::sampleSpeciesIndex(RNG& rng) const
{
  const double r = rng.uniform();
  auto it = std::lower_bound(cdf_.begin(), cdf_.end(), r);
  if (it == cdf_.end()) return static_cast<int>(cdf_.size()) - 1;
  return static_cast<int>(it - cdf_.begin());
}

int HQETFrag::applyHQSign(int basePdgPos, int hqPdg)
{
  const int s = signOf(hqPdg);
  return (s >= 0) ? basePdgPos : -basePdgPos;
}

bool HQETFrag::fragment(const Particle& HQ, RNG& rng, Particle& outHad) const
{
  if (!ready_) return false;

  // Charm only for now.
  if (std::abs(HQ.pdg) != 4) return false;

  const int idx = sampleSpeciesIndex(rng);
  const int basePdg = tab_[idx].pdg;
  const double mH = tab_[idx].mass;
  const double rH = rvals_[idx];
  const double ymax = ymax_[idx];
  const int outPdg = applyHQSign(basePdg, HQ.pdg);

  const double z = sampleZ(rng, rH, ymax);

  outHad = HQ;
  outHad.pdg = outPdg;
  outHad.m   = mH;

  outHad.px = z * HQ.px;
  outHad.py = z * HQ.py;
  outHad.pz = z * HQ.pz;

  const double p2 = outHad.px*outHad.px + outHad.py*outHad.py + outHad.pz*outHad.pz;
  outHad.E  = std::sqrt(p2 + mH*mH);

  return true;
}

} // namespace Frag
