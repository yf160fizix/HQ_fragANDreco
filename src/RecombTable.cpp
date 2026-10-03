#include "RecombTable.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {

bool isBlankOrComment(const std::string& line)
{
  const auto first = line.find_first_not_of(" \t\r");
  return first == std::string::npos || line[first] == '#';
}

bool readProbabilityFile(const std::string& path,
                         std::size_t channel_count,
                         std::vector<double>& momentum_grid,
                         std::vector<std::vector<double>>& probabilities)
{
  std::ifstream input(path);
  if (!input) {
    std::cerr << "Cannot open recombination probability table '" << path << "'.\n";
    return false;
  }

  momentum_grid.clear();
  probabilities.clear();

  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    if (isBlankOrComment(line)) continue;

    std::istringstream fields(line);
    double momentum = 0.0;
    double declared_total = 0.0;
    if (!(fields >> momentum >> declared_total) ||
        !std::isfinite(momentum) || !std::isfinite(declared_total)) {
      std::cerr << path << ':' << line_number
                << ": invalid momentum or total probability.\n";
      return false;
    }

    std::vector<double> row(channel_count, 0.0);
    double calculated_total = 0.0;
    for (std::size_t channel = 0; channel < channel_count; ++channel) {
      if (!(fields >> row[channel]) || !std::isfinite(row[channel]) ||
          row[channel] < 0.0) {
        std::cerr << path << ':' << line_number << ": invalid probability for channel "
                  << channel + 1 << ".\n";
        return false;
      }
      calculated_total += row[channel];
    }

    std::string extra;
    if (fields >> extra) {
      std::cerr << path << ':' << line_number << ": unexpected extra column.\n";
      return false;
    }

    constexpr double total_tolerance = 5e-4;
    if (std::abs(calculated_total - declared_total) >
        total_tolerance * std::max(1.0, std::abs(declared_total))) {
      std::cerr << path << ':' << line_number
                << ": declared total probability " << declared_total
                << " differs from channel sum " << calculated_total << ".\n";
      return false;
    }

    momentum_grid.push_back(momentum);
    probabilities.push_back(std::move(row));
  }

  return !momentum_grid.empty();
}

bool readWignerFile(const std::string& path,
                    std::vector<double>& momentum_grid,
                    std::vector<std::vector<double>>& maxima)
{
  std::ifstream input(path);
  if (!input) {
    std::cerr << "Cannot open Wigner-envelope table '" << path << "'.\n";
    return false;
  }

  momentum_grid.clear();
  maxima.clear();

  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    if (isBlankOrComment(line)) continue;

    std::istringstream fields(line);
    double momentum = 0.0;
    if (!(fields >> momentum) || !std::isfinite(momentum)) {
      std::cerr << path << ':' << line_number << ": invalid momentum.\n";
      return false;
    }

    std::vector<double> row(RecombTable::wigner_column_count, 0.0);
    for (std::size_t column = 0; column < row.size(); ++column) {
      if (!(fields >> row[column]) || !std::isfinite(row[column]) ||
          row[column] < 0.0) {
        std::cerr << path << ':' << line_number
                  << ": invalid Wigner envelope in column " << column + 1 << ".\n";
        return false;
      }
    }

    std::string extra;
    if (fields >> extra) {
      std::cerr << path << ':' << line_number << ": unexpected extra column.\n";
      return false;
    }

    momentum_grid.push_back(momentum);
    maxima.push_back(std::move(row));
  }

  return !momentum_grid.empty();
}

void buildCDF(
    const std::vector<std::vector<double>>& probabilities,
    std::vector<std::vector<double>>& cumulative)
{
  const std::size_t row_count = probabilities.size();
  const std::size_t channel_count = row_count == 0 ? 0 : probabilities.front().size();
  cumulative.assign(row_count, std::vector<double>(channel_count, 0.0));

  for (std::size_t row = 0; row < row_count; ++row) {
    double sum = 0.0;
    for (std::size_t channel = 0; channel < channel_count; ++channel) {
      sum += probabilities[row][channel];
      cumulative[row][channel] = sum;
    }

    if (sum > 1.0) {
      for (double& value : cumulative[row]) value /= sum;
    }
  }
}

bool compatibleGrids(const std::vector<double>& first,
                     const std::vector<double>& second,
                     double tolerance = 1e-12)
{
  if (first.size() != second.size()) return false;
  for (std::size_t i = 0; i < first.size(); ++i) {
    if (std::abs(first[i] - second[i]) >
        tolerance * (1.0 + std::abs(first[i]))) {
      return false;
    }
  }
  return true;
}

bool validUniformGrid(const std::vector<double>& grid, double spacing)
{
  if (grid.size() < 2 || !(spacing > 0.0)) return false;
  for (std::size_t i = 1; i < grid.size(); ++i) {
    if (!(grid[i] > grid[i - 1])) return false;
    const double delta = grid[i] - grid[i - 1];
    if (std::abs(delta - spacing) > 1e-10 * (1.0 + spacing)) return false;
  }
  return true;
}

} // namespace

std::size_t RecombTable::channelCount() const noexcept
{
  return channel_probabilities.empty() ? 0 : channel_probabilities.front().size();
}

bool loadRecombTable(const std::string& probability_file,
                     const std::string& wigner_file,
                     const HadConfig& config,
                     RecombTable& table)
{
  if (config.n_recomb_channels <= 0) {
    std::cerr << "Recombination channel count must be positive.\n";
    return false;
  }

  const auto channel_count =
      static_cast<std::size_t>(config.n_recomb_channels);
  if (!readProbabilityFile(probability_file, channel_count,
                           table.momentum_grid, table.channel_probabilities)) {
    return false;
  }

  std::vector<double> wigner_grid;
  if (!readWignerFile(wigner_file, wigner_grid, table.wigner_maxima)) {
    return false;
  }

  if (!compatibleGrids(table.momentum_grid, wigner_grid)) {
    std::cerr << "Momentum grids differ between '" << probability_file
              << "' and '" << wigner_file << "'.\n";
    return false;
  }

  const double expected_spacing = config.HQ_pid == 5
      ? config.dp_b
      : config.dp_c;
  if (!validUniformGrid(table.momentum_grid, expected_spacing)) {
    std::cerr << "Recombination table momentum grid is not strictly increasing with "
              << expected_spacing << " GeV spacing.\n";
    return false;
  }

  buildCDF(table.channel_probabilities, table.channel_cdf);
  return true;
}
