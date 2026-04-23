#pragma once

#include <random>

class RNG {
public:
    explicit RNG(unsigned seed = 12345);

    double uniform();               // [0,1)
    double uniform(double a, double b);
    int    uniformInt(int lo, int hi);

private:
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist;
};

