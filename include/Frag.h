#pragma once
#include <vector>
#include "Event.h"
#include "RNG.h"

namespace Frag {

struct Species {
  int pdg = 0;        // e.g. 421 (D0), 411 (D+), 431 (Ds+), 4122 (Lc+)
  double mass = 0.0;  // GeV
  double prob = 0.0;  // relative weight (need not sum to 1; we'll normalize)
};

// HQET/Braaten-type charm fragmentation sampler.
// Species are selected from a chemistry table; z is sampled from the
// HQET-inspired fragmentation function used in arXiv:1902.08889.
class HQETFrag {
public:
  HQETFrag() = default;

  void setChemistry(const std::vector<Species>& table);

  // Default chemistry: Table 2, RQM(170), final-state fractions for
  // D0, D+, Ds+, Lambda_c+. The sum is < 1 because Xi_c and Omega_c are
  // omitted here; probabilities are renormalized internally.
  void setDefaultCharmChemistry();

  // Fragment a charm quark/anti-quark Particle HQ into a charmed hadron.
  // - Assumes collinear fragmentation: p_h = z * p_Q (same direction)
  // - Keeps HQ space-time / hydro fields by copying HQ -> outHad
  // - Chooses hadron PDG sign consistent with HQ pdg sign.
  bool fragment(const Particle& HQ, RNG& rng, Particle& outHad) const;

private:
  std::vector<Species> tab_;      // normalized probs
  std::vector<double> cdf_;       // cumulative distribution [0,1]
  std::vector<double> rvals_;     // per-species HQET r parameter
  std::vector<double> ymax_;      // cached rejection envelope maxima
  bool ready_ = false;

  static bool isBaryonPDG(int pdgAbs);
  static double defaultRForSpecies(int pdg, double mass);
  static double hqetKernel(double z, double r);
  static double estimateYmax(double r,
                             double zmin = 1e-4, double zmax = 1.0 - 1e-6);
  static int signOf(int x) { return (x >= 0) ? +1 : -1; }
  static int applyHQSign(int basePdgPos, int hqPdg);

  double sampleZ(RNG& rng, double r, double ymax,
                 double zmin = 1e-4, double zmax = 1.0 - 1e-6) const;
  int sampleSpeciesIndex(RNG& rng) const;
};

} // namespace Frag
