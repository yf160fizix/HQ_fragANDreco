#pragma once

#include <string>
#include <vector>

namespace Frag {

struct DecayChannel {
  int daughter_pid = 0;
  double m_daughter = 0.0;  // GeV
  double m_companion = 0.13957; // GeV
  double BR = 0.0;
};

struct CharmHadronState {
  std::string name;
  int pid = 0;
  double mass = 0.0; // GeV
  int strangeness = 0;
  double spin = 0.0;
  double isospin = 0.0;
  bool baryon = false;
  double statistical_weight = 0.0;
  std::vector<DecayChannel> decay_channels;
};

std::vector<CharmHadronState> loadCharmStates(
    const std::string& meson_path,
    const std::string& baryon_path,
    double Tchem,
    double gamma_s,
    double gamma_HB);

} // namespace Frag
