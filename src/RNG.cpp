#include "RNG.h"

RNG::RNG(unsigned seed)
    : engine_(seed), unit_distribution_(0.0, 1.0) {}

double RNG::uniform()
{
  return unit_distribution_(engine_);
}
