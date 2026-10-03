#pragma once

#include <random>

class RNG {
public:
  explicit RNG(unsigned seed = 12345);

  double uniform(); // [0, 1)

private:
  std::mt19937 engine_;
  std::uniform_real_distribution<double> unit_distribution_;
};
