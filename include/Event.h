#pragma once
#include <vector>
#include "Particle.h"

struct Event {
  int event_id = 0;
  std::vector<Particle> particles;
};

