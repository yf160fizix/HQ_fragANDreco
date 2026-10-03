#pragma once

#include <vector>

#include "CharmStateTable.h"
#include "Particle.h"
#include "RNG.h"

namespace Frag {

class CharmFeeddown {
public:
  bool decay(Particle& hadron,
             const std::vector<DecayChannel>& channels,
             RNG& rng) const;

  bool decayKnownState(Particle& hadron, RNG& rng) const;

private:
  static void twoBodyDecay(Particle& mother,
                           double m_daughter,
                           double m_companion,
                           RNG& rng);
};

} // namespace Frag
