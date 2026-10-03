#pragma once

#include <string>

enum class HadMode {
  Frag = 1,
  Recomb = 2,
  FragAndRecomb = 3,
};

struct HadConfig {
  static constexpr double default_omega_cM = 0.20;
  static constexpr double default_omega_cB = 0.267;
  static constexpr double default_omega_bM = 0.14;
  static constexpr double default_omega_bB = 0.14;

  HadMode mode = HadMode::FragAndRecomb;
  int HQ_pid = 4;
  bool rescale_HQ_mass = true;
  bool enable_feeddown = true;

  // Masses, temperatures, oscillator scales, and momentum spacings use GeV.
  double m_c = 1.8;
  double m_b = 5.2;
  double m_q = 0.30;
  double m_s = 0.40;
  double m_g = 0.30;

  double omega_cM = default_omega_cM;
  double omega_cB = default_omega_cB;
  double omega_bM = default_omega_bM;
  double omega_bB = default_omega_bB;

  double Tchem = 0.160;
  double gamma_s = 0.7;
  double gamma_HB = 1.0;
  double eps_M = 0.01;
  double eps_B = 0.03;

  std::string charm_meson_table = "data/charm_meson_states.csv";
  std::string charm_baryon_table = "data/charm_baryon_states.csv";
  std::string charm_recomb_table = "data/recomb_c_raw_M020_B0267.dat";
  std::string charm_wigner_table = "data/max_wigner_c_M020_B0267.dat";

  int n_recomb_channels = 24;
  double dp_c = 0.5;
  double dp_b = 0.5;
  double default_baryon_s_wave_envelope = 0.1;
  double default_baryon_p_wave_envelope = 0.2;
};
