#include "Lorentz.h"

#include <algorithm>
#include <cmath>

namespace Lorentz {

void boost(double bx, double by, double bz,
           double& px, double& py, double& pz, double& E)
{
  double b2 = bx * bx + by * by + bz * bz;
  if (b2 < 1e-30) return;

  // Clamp invalid input velocities just below the speed of light.
  if (b2 >= 1.0) {
    const double scale = 0.999999 / std::sqrt(b2);
    bx *= scale;
    by *= scale;
    bz *= scale;
    b2 = bx * bx + by * by + bz * bz;
  }

  const double one_minus_b2 = 1.0 - b2;
  const double gamma = 1.0 / std::sqrt(std::max(1e-15, one_minus_b2));

  const double bp = bx * px + by * py + bz * pz;
  const double gamma2 = (gamma - 1.0) / b2;

  const double px_new = px + gamma2 * bp * bx + gamma * bx * E;
  const double py_new = py + gamma2 * bp * by + gamma * by * E;
  const double pz_new = pz + gamma2 * bp * bz + gamma * bz * E;
  const double E_new = gamma * (E + bp);

  px = px_new;
  py = py_new;
  pz = pz_new;
  E = E_new;
}

void boost(double bx, double by, double bz, FourMomentum& momentum)
{
  boost(bx, by, bz,
        momentum.px, momentum.py, momentum.pz, momentum.E);
}

} // namespace Lorentz
