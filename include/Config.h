#pragma once
#include <string>

struct Config {
  int hadr_flag =1;          // 1: frag, 2: recomb-only, 3: frag+recomb
  int HQid = 4;               // 4 charm, 5 bottom
  int rescale_m = 1;         
  bool heavyMesonDecay = false;

  //Mass parameters
  double m_c = 1.3;
  double m_b = 5.2;
  double m_ud = 0.30; //Thermal mass of ud quark
  double m_s  = 0.40; //Thermal mass of s quark
  double m_g  = 0.30; //Thermal mass of gluon


  // omega parameters in the recombination formula
  double omega_c = 0.215; // Angular frequancy of the S.H.O. of c-hadron
  double omega_b = 0.102; // Angular frequancy of the S.H.O. of b-hadron
  double omega_cM = 0.24; // Angular frequancy for c-Meson
  double omega_cB = 0.24; // Angular frequancy for c-Baryon
  double omega_bM = 0.14; // Angular frequancy for b-Meson
  double omega_bB = 0.14; // Angular frequancy for b-Baryon


  

  // I/O 
  bool use_stdin = true;      // if false and input_path empty --> use toy event
  std::string input_path;     // optional
  std::string output_path;    // optional
  unsigned long long seed = 12345;

  // parameters of the recombination probability table

  int prob_num = 61; // Total number of points
  double prob_int_c = 0.5; // Bin size for c quark
  double prob_int_b = 0.5; // Bin size for b quark

  // Others

  double max_3body_s=0.1;
  double max_3body_p=0.2;

  int nChannel =24;

  double eps_c =0.05;
  double eps_b =0.05;

};
