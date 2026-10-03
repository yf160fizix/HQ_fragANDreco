#pragma once

#include <string>
#include <vector>

#include "HadConfig.h"

struct RecombTable {
  static constexpr std::size_t wigner_column_count = 10;
  static constexpr std::size_t baryon_s_wave_column = 4;
  static constexpr std::size_t baryon_p_wave_column = 5;

  std::vector<double> momentum_grid;
  std::vector<std::vector<double>> channel_probabilities;
  std::vector<std::vector<double>> channel_cdf;
  std::vector<std::vector<double>> wigner_maxima;

  std::size_t channelCount() const noexcept;
};

bool loadRecombTable(const std::string& probability_file,
                     const std::string& wigner_file,
                     const HadConfig& config,
                     RecombTable& table);
