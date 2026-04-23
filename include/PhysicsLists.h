#pragma once
#include <array>

class PhysicsLists {
public:
  static constexpr int kMaxIndex = 70; // 0..70

  static const std::array<double, kMaxIndex + 1> max_ud;
  static const std::array<double, kMaxIndex + 1> max_s;
};

