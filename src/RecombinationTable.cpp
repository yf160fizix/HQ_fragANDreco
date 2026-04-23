#include "RecombinationTable.h"
#include <fstream>
#include <sstream>
#include <cmath>

static bool readRawFileMultiChannel(const std::string& file,
                                    int nChannel,
                                    std::vector<double>& pc,
                                    std::vector<std::vector<double>>& prob_raw_ch)
{
  std::ifstream fin(file);
  
  if (!fin) {
    std::cerr << "[readRecombTable] Cannot open raw file: " << file << "\n";
    return false;
  }
  //std::cout<< "test" << std::endl;
  pc.clear();
  prob_raw_ch.clear();

  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    // allow comment lines
    if (line[0] == '#') continue;

    std::istringstream iss(line);
    double dummy_pc = 0.0, dummy_prob = 0.0;
    if (!(iss >> dummy_pc >> dummy_prob)) continue;

    std::vector<double> ch(nChannel, 0.0);
    for (int c = 0; c < nChannel; ++c) {
      if (!(iss >> ch[c])) {
        std::cerr << "[readRecombTable] raw format error: need "
                  << nChannel << " channel probs per line.\n"
                  << "  line: " << line << "\n";
        return false;
      }
    }

    pc.push_back(dummy_pc);
    prob_raw_ch.push_back(std::move(ch));
  }

  return !pc.empty();
}

static bool readWignerFile10(const std::string& file,
                             std::vector<double>& pc,
                             std::vector<std::vector<double>>& wigner_max)
{
  std::ifstream fin(file);
  if (!fin) {
    std::cerr << "[readRecombTable] Cannot open wigner file: " << file << "\n";
    return false;
  }

  pc.clear();
  wigner_max.clear();

  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream iss(line);
    double dummy_pc = 0.0;
    if (!(iss >> dummy_pc)) continue;

    std::vector<double> w(RecombinationTable::kNWigner, 0.0);
    for (int j = 0; j < RecombinationTable::kNWigner; ++j) {
      if (!(iss >> w[j])) {
        std::cerr << "[readRecombTable] wigner format error: need 10 cols per line.\n"
                  << "  line: " << line << "\n";
        return false;
      }
    }

    pc.push_back(dummy_pc);
    wigner_max.push_back(std::move(w));
  }

  return !pc.empty();
}

static void buildCDFAndNormalize(const std::vector<std::vector<double>>& prob_raw_ch,
                                 std::vector<std::vector<double>>& prob_cum_ch)
{
  const size_t Npc = prob_raw_ch.size();
  const int nChannel = (Npc > 0) ? (int)prob_raw_ch[0].size() : 0;

  prob_cum_ch.assign(Npc, std::vector<double>(nChannel, 0.0));

  for (size_t i = 0; i < Npc; ++i) {
    double acc = 0.0;
    for (int c = 0; c < nChannel; ++c) {
      acc += prob_raw_ch[i][c];
      prob_cum_ch[i][c] = acc;
    }

    // Fortran: if last > 1, normalize by last
    const double last = prob_cum_ch[i][nChannel - 1];
    if (last > 1.0) {
      for (int c = 0; c < nChannel; ++c) prob_cum_ch[i][c] /= last;
    }
  }
}

static bool pcGridsCompatible(const std::vector<double>& a,
                              const std::vector<double>& b,
                              double tol = 1e-12)
{
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i) {
    if (std::fabs(a[i] - b[i]) > tol * (1.0 + std::fabs(a[i]))) return false;
  }
  return true;
}

bool readRecombTable(int HQid,
                     const std::string& rawFile,
                     const std::string& wignerFile,
                     RecombinationTable& table,
                     const Config& cfg)
{
 

  const int nChannel = cfg.nChannel;
  if (nChannel <= 0) {
    std::cerr << "[readRecombTable] cfg.nChannel invalid: " << nChannel << "\n";
    return false;
  }

  table.nChannel = nChannel;
 // std::cout << "nChannel= " << nChannel << std::endl;
  //std::cout << rawFile << std::endl;

  // 1) raw: pc, dummy_prob, prob_raw_ch(1:nChannel)
  if (!readRawFileMultiChannel(rawFile, nChannel, table.pc_raw, table.prob_raw_ch)) {
    std::cerr << "[readRecombTable] raw read failed.\n";
    return false;
  }

  // 2) CDF + normalize (Fortran prob_c)
  buildCDFAndNormalize(table.prob_raw_ch, table.prob_cum_ch);

  // 3) wigner: pc, wigner_max(1:10)
  if (!readWignerFile10(wignerFile, table.pc_wigner, table.wigner_max)) {
    std::cerr << "[readRecombTable] wigner read failed.\n";
    return false;
  }

  // 4) sanity: pc grids match
  if (!pcGridsCompatible(table.pc_raw, table.pc_wigner)) {
    std::cerr << "[readRecombTable] pc grid mismatch between raw and wigner files.\n";
    std::cerr << "  raw size=" << table.pc_raw.size()
              << " wigner size=" << table.pc_wigner.size() << "\n";
    return false;
  }

  // 5) fill legacy vectors (optional, but keeps old code happy)
  // prob_raw: store cum prob at last channel (typically 1 after normalize)
  table.prob_raw.resize(table.pc_raw.size());
  for (size_t i = 0; i < table.pc_raw.size(); ++i) {
    table.prob_raw[i] = table.prob_cum_ch[i][nChannel - 1];
  }

  // prob_wigner: store first wigner column 
  table.prob_wigner.resize(table.pc_wigner.size());
  for (size_t i = 0; i < table.pc_wigner.size(); ++i) {
    table.prob_wigner[i] = table.wigner_max[i][0];
  }

  return true;
}

