#include "Lorentz.h"
#include <cmath>
#include <algorithm>

namespace Lorentz {

void boost(double bx, double by, double bz,
           double& px, double& py, double& pz, double& E)
{
  double b2 = bx*bx + by*by + bz*bz;
  if (b2 < 1e-30) return;

  // Clamp superluminal beta to avoid invalid gamma / crazy numbers
  if (b2 >= 1.0) {
    const double scale = 0.999999 / std::sqrt(b2);
    bx *= scale; by *= scale; bz *= scale;
    b2 = bx*bx + by*by + bz*bz;
  }

  // gamma = 1/sqrt(1-b^2)
  const double one_minus_b2 = 1.0 - b2;
  const double gamma = 1.0 / std::sqrt(std::max(1e-15, one_minus_b2));

  const double bp = bx*px + by*py + bz*pz;          // beta · p
  const double gamma2 = (gamma - 1.0) / b2;         // (gamma-1)/b^2

  const double px_new = px + gamma2*bp*bx + gamma*bx*E;
  const double py_new = py + gamma2*bp*by + gamma*by*E;
  const double pz_new = pz + gamma2*bp*bz + gamma*bz*E;
  const double E_new  = gamma*(E + bp);

  px = px_new; py = py_new; pz = pz_new; E = E_new;
}

} // namespace Lorentz


/*
FourVec Lorentz::boost( double bx, double by, double bz, const FourVec& v) {
  const double b2 = bx*bx + by*by + bz*bz;
  if (b2 <= 0.0) return v;

  // Protect against invalid beta
  if (b2 >= 1.0) {
    // clamp slightly below 1
    const double scale = 0.999999 / std::sqrt(b2);
    bx *= scale; by *= scale; bz *= scale;
  }

  const double b2c = bx*bx + by*by + bz*bz;
  const double gamma = 1.0 / std::sqrt(1.0 - b2c);

  const double bp = bx*v.px + by*v.py + bz*v.pz;
  const double gamma2 = (b2c > 0.0) ? ( (gamma - 1.0) / b2c ) : 0.0;

  FourVec out;
  out.px = v.px + gamma2*bp*bx + gamma*bx*v.E;
  out.py = v.py + gamma2*bp*by + gamma*by*v.E;
  out.pz = v.pz + gamma2*bp*bz + gamma*bz*v.E;
  out.E  = gamma*(v.E + bp);

  return out;
}
*/

