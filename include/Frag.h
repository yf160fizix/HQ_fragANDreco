#pragma once

#include <string>
#include <vector>

#include "CharmStateTable.h"
#include "CharmFeeddown.h"
#include "Particle.h"
#include "RNG.h"

namespace Frag {

class PetersonFrag {
public:
  PetersonFrag(double eps_M, double eps_B);

  void setChemistry(const std::vector<CharmHadronState>& table);

  bool loadCharmChemistry(const std::string& meson_path,
                          const std::string& baryon_path,
                          double Tchem,
                          double gamma_s,
                          double gamma_HB);

  bool fragment(const Particle& HQ,
                RNG& rng,
                Particle& hadron) const;

  bool feedDownByPDG(Particle& hadron, RNG& rng) const;

private:
  std::vector<CharmHadronState> species_;
  std::vector<double> species_cdf_;
  std::vector<double> species_epsilon_;
  std::vector<double> kernel_maxima_;

  double eps_M_;
  double eps_B_;
  bool ready_ = false;

  static double petersonKernel(double z, double epsilon);
  static double estimateKernelMax(double epsilon,
                                  double zmin = 1e-4,
                                  double zmax = 1.0 - 1e-6);

  static int applyHQSign(int positive_pid, int HQ_pid);

  double epsForSpecies(const CharmHadronState& state) const;

  double sampleZ(RNG& rng,
                 double epsilon,
                 double kernel_max,
                 double zmin = 1e-4,
                 double zmax = 1.0 - 1e-6) const;

  int sampleSpeciesIndex(RNG& rng) const;

  CharmFeeddown feeddown_;
};

} // namespace Frag
