#pragma once

struct Particle {
  int pdg = 0;

  // momentum (px py pz p0)
  double px=0, py=0, pz=0, E=0, m=0;

  // position (rx ry rz r0) 
  double x=0, y=0, z=0, t=0;

  // hydro cell info
  double Thydro = 0.0;
  double cvx = 0.0, cvy = 0.0, cvz = 0.0;

  double ipx=0, ipy=0, ipz=0, iE=0;
  double ix=0, iy=0, iz=0;

  // weight
  double wt = 1.0;

};

