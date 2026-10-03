#pragma once

#include <cstddef>
#include <map>
#include <optional>
#include <vector>

#include "HadConfig.h"
#include "Frag.h"
#include "Particle.h"
#include "RNG.h"
#include "RecombChannel.h"
#include "RecombTable.h"

struct HadStats {
  std::size_t n_HQ = 0;
  std::size_t n_frag = 0;
  std::size_t n_recomb = 0;
  std::size_t n_dropped = 0;
  std::size_t n_unchanged_HQ = 0;
  std::map<int, std::size_t> final_species;
};

class Hadronizer {
public:
  Hadronizer(const HadConfig& config,
             RNG& rng,
             const RecombTable& table);

  void process(const std::vector<Particle>& input,
               std::vector<Particle>& output);

  bool ready() const noexcept { return ready_; }
  const HadStats& stats() const noexcept { return stats_; }

private:
  enum class ChannelOutcome { Recomb, Failed, Dropped };

  struct RecombKinematics {
    double M = 0.0;
    double beta_x = 0.0;
    double beta_y = 0.0;
    double beta_z = 0.0;
  };

  struct BaryonWignerWidths {
    double sigma_11 = 0.0;
    double sigma_12 = 0.0;
    double sigma_21 = 0.0;
    double sigma_22 = 0.0;
    double sigma_31 = 0.0;
    double sigma_32 = 0.0;
  };

  struct BaryonJacobiMomenta {
    double q11_sq = 0.0;
    double q12_sq = 0.0;
    double q21_sq = 0.0;
    double q22_sq = 0.0;
    double q31_sq = 0.0;
    double q32_sq = 0.0;
  };

  bool tryRecombMeson(const Particle& HQ,
                      LightFlavor light_flavor,
                      OrbitalState orbital,
                      int hadron_pid,
                      double hadron_mass,
                      Particle& hadron);

  bool tryRecombBaryon(const Particle& HQ,
                       LightFlavor light_flavor1,
                       LightFlavor light_flavor2,
                       OrbitalState orbital,
                       int hadron_pid,
                       double hadron_mass,
                       Particle& hadron);

  bool sampleMesonKinematics(OrbitalState orbital,
                             const Particle& HQ,
                             double m_light,
                             double sigma,
                             double q_rel_sq_min,
                             double thermal_max,
                             RecombKinematics& result);

  bool sampleBaryonKinematics(OrbitalState orbital,
                              const Particle& HQ,
                              double m_light1,
                              double m_light2,
                              double thermal_max_1,
                              double thermal_max_2,
                              const BaryonWignerWidths& widths,
                              RecombKinematics& result);

  std::optional<double> thermalMax(
      LightFlavor flavor, double T) const;
  FourMomentum sampleThermalParton(double T,
                                   double mass,
                                   double thermal_max);
  void buildRecombHadron(const Particle& HQ,
                         int hadron_pid,
                         double hadron_mass,
                         const RecombKinematics& kinematics,
                         Particle& hadron);
  static BaryonWignerWidths computeBaryonWignerWidths(
      double m_HQ,
      double m_light1,
      double m_light2,
      double omega);
  static BaryonJacobiMomenta computeJacobiMomenta(
      const FourMomentum& p_HQ,
      const FourMomentum& p_light1,
      const FourMomentum& p_light2);
  static double baryonWignerWeight(
      OrbitalState orbital,
      const BaryonWignerWidths& widths,
      const BaryonJacobiMomenta& momenta);

  ChannelOutcome executeChannel(std::size_t channel_index,
                                const Particle& HQ,
                                Particle& hadron);
  bool interpolateChannelCDF(const Particle& HQ,
                             std::vector<double>& cdf);
  static std::optional<std::size_t> drawChannel(
      const std::vector<double>& cdf, double uniform_deviate);

  void processHQ(const Particle& HQ,
                 std::vector<Particle>& output);
  void fragOrKeep(const Particle& HQ,
                  std::vector<Particle>& output);
  void restoreLabState(Particle& HQ,
                       double original_mass) const;
  void appendResult(const Particle& particle,
                    std::vector<Particle>& output);
  static Particle makeHadron(const Particle& HQ,
                             int hadron_pid,
                             double hadron_mass);

  void adjustMassShellByTwoBodyDecay(double& px,
                                     double& py,
                                     double& pz,
                                     double& E,
                                     double final_mass);

  HadConfig config_;
  RNG& rng_;
  const RecombTable& recomb_table_;
  Frag::PetersonFrag frag_;
  HadStats stats_;

  double m_HQ_ = 1.8;
  double omega_M_ = 0.20;
  double baryon_s_wave_envelope_ = 0.1;
  double baryon_p_wave_envelope_ = 0.2;
  bool ready_ = true;
};
