#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>

#include "RNG.h"
#include "Config.h"
#include "Hadronizer.h"
#include "IO.h"
#include "Event.h"
#include "RecombinationTable.h"

int main(int argc, char** argv)
{
  // -----------------------------
  // 0) Args
  // -----------------------------
  if (argc < 3) {
    std::cerr
      << "Usage:\n"
      << "  " << argv[0] << " input.txt output.txt [seed] [HQid] [recomb_raw] [max_wigner]\n\n"
      << "Example:\n"
      << "  " << argv[0] << " input.dat output.dat 12345 4 recomb_c_raw.dat max_wigner_c.dat\n\n";
    return 1;
  }

  const std::string inFile  = argv[1];
  const std::string outFile = argv[2];

  const unsigned seed =
    (argc >= 4) ? static_cast<unsigned>(std::stoul(argv[3]))
                : std::random_device{}();

  const int hqid = (argc >= 5) ? std::stoi(argv[4]) : 4;
  const std::string recombRawFile = (argc >= 6) ? argv[5] :
      ((hqid == 5) ? "recomb_b_raw.dat" : "recomb_c_raw.dat");
  const std::string maxWignerFile = (argc >= 7) ? argv[6] :
      ((hqid == 5) ? "max_wigner_b.dat" : "max_wigner_c.dat");

  std::ifstream fin(inFile);
  if (!fin) {
    std::cerr << "Error: cannot open input file: " << inFile << "\n";
    return 2;
  }

  std::ofstream fout(outFile);
  if (!fout) {
    std::cerr << "Error: cannot open output file: " << outFile << "\n";
    return 3;
  }

  // -----------------------------
  // 1) RNG
  // -----------------------------
  RNG rng(seed);

  // -----------------------------
  // 2) Config 
  // -----------------------------
  Config cfg;

  // light parton masses
  cfg.m_ud = 0.30;
  cfg.m_s  = 0.40;
  cfg.m_g  = 0.30;

  // heavy-quark choice 
  cfg.HQid = hqid;      // 4=c, 5=b
  cfg.m_c  = 1.8; //1.3;
  cfg.m_b  = 5.2; //4.8;

  // oscillator params 
  cfg.omega_cM = 0.24;
  cfg.omega_cB = 0.24;
  cfg.omega_bM = 0.14;
  cfg.omega_bB = 0.14;

  // hadronization control
  cfg.hadr_flag  = 3; //3;     // 1: fragmentation, 2: recomb-only, 3: frag+reco
  cfg.prob_int_c = 0.5;   

  // -----------------------------
  // 3) Hadronizer
  // -----------------------------
  RecombinationTable recombTable;
  if (!readRecombTable(cfg.HQid, recombRawFile, maxWignerFile, recombTable, cfg)) {
    std::cerr << "Failed to read recombination data files: "
              << recombRawFile << " and " << maxWignerFile << "\n";
    return 4;
  }
  
 /* 
  //print the first 10 lines as a test
  size_t nPrint = std::min<size_t>(10, recombTable.pc_raw.size());

  for (size_t i = 0; i < nPrint; ++i) {
  std::cout << recombTable.pc_raw[i] << "  ";

  for (int c = 0; c < recombTable.nChannel; ++c) {
    std::cout << recombTable.prob_raw_ch[i][c] << " ";
  }
  std::cout << "\n";
  }*/

/*
  for (size_t i = 0; i < std::min<size_t>(10, recombTable.pc_raw.size()); ++i) {
	  std::cout << recombTable.pc_raw[i] << "  ";
	  for (int c = 0; c < recombTable.nChannel; ++c) {
		  std::cout << recombTable.prob_cum_ch[i][c] << " ";
	  }
	  std::cout << "\n";
  }
*/
  Hadronizer had(cfg, rng, recombTable);

  // -----------------------------
  // 4) Event loop: read -> hadronize -> write
  // -----------------------------
  Event ev;
  ev.event_id = 0;

  long long nEvents = 0;
  long long nInParticles = 0;
  long long nOutParticles = 0;

 // std::cout << "start" << std::endl;
  while (true) {
    const bool read_flag = IO::readEventText(fin, ev);
    if (!read_flag){ 
	    //std::cout << "break" << std::endl;
	    //std::cout << "nEvents= " << nEvents << std::endl;
	    break;
    }
 //  std::cout << "readed"<< std::endl;
    ++nEvents;
    nInParticles += static_cast<long long>(ev.particles.size());

    std::vector<Particle> outParticles;
    had.process(ev.particles, outParticles);
    nOutParticles += static_cast<long long>(outParticles.size());

    Event evOut = ev;
    evOut.particles = std::move(outParticles);

    IO::writeEventText(fout, evOut);
  }

  // -----------------------------
  // 5) If no events were read, run a toy event
  // -----------------------------
  if (nEvents == 0) {
    std::cerr << "Warning: no events read from input. "
              << "Input format may be wrong.\n"
              << "Running a toy event and writing it to output for debugging.\n";

    Event toy = IO::makeToyEvent(0);
    std::vector<Particle> outParticles;
    had.process(toy.particles, outParticles);

    toy.particles = std::move(outParticles);
    IO::writeEventText(fout, toy);

    nEvents = 1;
  }

  // -----------------------------
  // 6) Summary
  // -----------------------------
  std::cerr
    << " Input : " << inFile  << "\n"
    << " Output: " << outFile << "\n"
    << " Seed  : " << seed    << "\n"
    << " HQ id : " << cfg.HQid << "\n"
    << " Hadronization mode : " << cfg.hadr_flag << "\n"
    << " Recomb table : " << recombRawFile << "\n"
    << " Wigner table : " << maxWignerFile << "\n"
    << " Events processed       : " << nEvents << "\n"
    << " Total input particles  : " << nInParticles << "\n"
    << " Total output particles : " << nOutParticles << "\n"
    << " Done." << std::endl ;

  return 0;
}
