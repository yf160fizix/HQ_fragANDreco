#pragma once

struct FourMomentum {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double E = 0.0;
};

enum class HadOrigin {
  Unknown = 0,
  Frag = 1,
  Recomb = 2,
};

struct Particle {
  int pid = 0;

  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double E = 0.0;
  double m = 0.0;

  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double t = 0.0;

  // Local fluid-cell state.
  double Thydro = 0.0;
  double cvx = 0.0;
  double cvy = 0.0;
  double cvz = 0.0;

  // Initial momentum retained in the output record.
  double ipx = 0.0;
  double ipy = 0.0;
  double ipz = 0.0;
  double iE = 0.0;

  double ix = 0.0;
  double iy = 0.0;
  double iz = 0.0;

  double wt = 1.0;
  HadOrigin origin = HadOrigin::Unknown;
};
