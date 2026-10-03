#pragma once

#include <cstddef>

enum class LightFlavor { UpDown, Strange };
enum class OrbitalState { S, P };

enum class RecombChannelId {
  D,
  DStar,
  D1_2420,
  D0Star_2300,
  D1_2430,
  D2Star_2460,
  Ds,
  DsStar,
  Ds1_2536,
  Ds0Star_2317,
  Ds1_2460,
  Ds2Star_2573,
  LambdaC,
  LambdaC2860,
  LambdaCPWave,
  SigmaC,
  SigmaCStar,
  SigmaCPWave,
  XiCGroup1,
  XiCGroup2,
  XiCGroup3,
  OmegaCGroup1,
  OmegaCGroup2,
  OmegaCGroup3,
  B,
  BStar,
  B1,
  B0Star,
  B1Star,
  B2Star,
};

enum class ConstituentSystem { Meson, Baryon, Unsupported };

struct HadronState {
  int pid = 0;
  double mass = 0.0;
};

struct RecombChannel {
  RecombChannelId id;
  ConstituentSystem constituents;
  LightFlavor light_flavor1;
  LightFlavor light_flavor2;
  OrbitalState orbital;
  HadronState primary;
  HadronState alternative;
};

const RecombChannel* findRecombChannel(
    int HQ_pid, std::size_t channel_index) noexcept;
