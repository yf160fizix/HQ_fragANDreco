#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "Config.h"

struct RecombinationTable {

  
  std::vector<double> pc_raw;
  std::vector<double> prob_raw;

  std::vector<double> pc_wigner;
  std::vector<double> prob_wigner;

  
  int nChannel = 0;
  std::vector<std::vector<double>> prob_raw_ch; // [Npc][nChannel]

  // cumulative probabilities (CDF)
  std::vector<std::vector<double>> prob_cum_ch; // [Npc][nChannel]

  // wigner max table 
  static constexpr int kNWigner = 10;
  std::vector<std::vector<double>> wigner_max;  // [Npc][10]

  
  void printFirst(size_t N = 5) const {
//    std::cout << "[RecombTable] Npc=" << pc_raw.size()
//              << " nChannel=" << nChannel << "\n";
//    std::cout << " raw: pc  cum(last)  cum(ch=24)\n";
    for (size_t i = 0; i < std::min(N, pc_raw.size()); ++i) {
      double last = (nChannel > 0 && i < prob_cum_ch.size()) ? prob_cum_ch[i][nChannel-1] : -1.0;
      double ch24 = (nChannel >= 24 && i < prob_cum_ch.size()) ? prob_cum_ch[i][23] : -1.0;
      //std::cout << " " << pc_raw[i] << "  " << last << "  " << ch24 << "\n";
    }
//    std::cout << " wigner: pc  col1\n";
    for (size_t i = 0; i < std::min(N, pc_wigner.size()); ++i) {
      double c1 = (i < wigner_max.size() && !wigner_max[i].empty()) ? wigner_max[i][0] : -1.0;
//      std::cout << " " << pc_wigner[i] << "  " << c1 << "\n";
    }
  }
};

bool readRecombTable(int HQid,
                     const std::string& rawFile,
                     const std::string& wignerFile,
                     RecombinationTable& table,
                     const Config& cfg);

