#include "RNG.h"

RNG::RNG(unsigned seed)
    : gen(seed), dist(0.0, 1.0) {}

double RNG::uniform() {
    return dist(gen);
}

double RNG::uniform(double a, double b) {
    return a + (b - a) * uniform();
}

int RNG::uniformInt(int lo, int hi) {
    std::uniform_int_distribution<int> d(lo, hi);
    return d(gen);
}

