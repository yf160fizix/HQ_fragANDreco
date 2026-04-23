#pragma once
#include <cmath>

namespace Lorentz {

// Boost by beta vector (bx,by,bz), in-place modifies (px,py,pz,E)
void boost(double bx, double by, double bz,
           double& px, double& py, double& pz, double& E);

} // namespace Lorentz

/*
struct FourVec {
  double px = 0.0, py = 0.0, pz = 0.0, E = 0.0;
};

namespace Lorentz {
  // Boost by velocity (bx,by,bz) with |b|<1.
  // This is the C++ analog of your rotbos usage for 4-momenta.
  FourVec boost( double bx, double by, double bz, const FourVec& v);
}
*/

