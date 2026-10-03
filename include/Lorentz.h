#pragma once

#include "Particle.h"

namespace Lorentz {

// Applies p' = Lambda(beta) p using the (+,-,-,-) energy convention.
void boost(double bx, double by, double bz,
           double& px, double& py, double& pz, double& E);

void boost(double bx, double by, double bz, FourMomentum& momentum);

} // namespace Lorentz
