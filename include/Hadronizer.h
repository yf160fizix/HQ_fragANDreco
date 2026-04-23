#pragma once
#include <cstdint>
#include <vector>
#include <random>
#include "Particle.h"
#include "RNG.h"
#include "Config.h"
#include "RecombinationTable.h"
#include "Frag.h"

class Hadronizer {
public:	

  Hadronizer(const Config& cfg, RNG& rng,const RecombinationTable& table);
 
  // =====================
  // (used by hadr_flag==3, channel decides species)
  // =====================
  bool tryRecombine2(const Particle& HQ, int lightFlavor, int orbital,
                     int pdgWanted, double mWanted, Particle& outHadron);

  bool tryRecombine3(const Particle& HQ,
                     int lightFlavor1, int lightFlavor2, int orbital,
                     int pdgWanted, double mWanted, Particle& outHadron);

  void process(const std::vector<Particle>& input,
               std::vector<Particle>& output);

private:

  // Core implementations for new interfaces
  bool forceRecomb2(int& todoReco, const Particle& HQ,
                    int lightFlavor, int orbital,
                    int pdgWanted, double mWanted,
                    Particle& outHadron);

  bool forceRecomb3(int& todoReco, const Particle& HQ,
                    int lightFlavor1, int lightFlavor2, int orbital,
                    int pdgWanted, double mWanted,
                    Particle& outHadron);

  

  void judgeHadMech2(int orbital, int& todoReco,
                   const Particle& HQ,
                   double mlight,
                   double sigma0, double kmin2,
                   double mxvalue,
                   double& directM, double& vx_HM, double& vy_HM, double& vz_HM,
                   double& px_l, double& py_l, double& pz_l, double& E_l);

  void judgeHadMech3(int orbital, int& todoReco,
                     const Particle& HQ,
                     double mlight1, double mlight2,
                     double mxvalue1, double mxvalue2,
                     double sigma11, double sigma12,
                     double sigma21, double sigma22,
                     double sigma31, double sigma32,
		    // double probmax_s, double probmax_p,
                     double& directM, double& vx_HM, double& vy_HM, double& vz_HM);

  double pLengthInCellLRF(const Particle& HQ) const;
  double pLength(const Particle& HQ) const;
  bool interpolatedCDF(const Particle& HQ, std::vector<double>& cdf) const;
  int drawChannelFromCDF(const std::vector<double>& cdf, double u01) const;


  void decay2body(double& px, double& py, double& pz, double& E,
                  double m_final);

  // helpers
  static double sqr(double x) { return x*x; }

  // Fortran rotbos equivalent (boost by beta vector)
  void boost(double bx, double by, double bz,
             double& px, double& py, double& pz, double& E) const;

  Config cfg_;
  RNG& rng_;
  RecombinationTable table_;

  Frag::HQETFrag frag_;

  // derived params
  double omega_HM_;
  double m_HQ_ = 1.8;
  mutable double current_max_3body_s_ = 0.1;
  mutable double current_max_3body_p_ = 0.2;
};
