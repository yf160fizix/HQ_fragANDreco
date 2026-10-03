#include "HadConfig.h"
#include "Event.h"
#include "Hadronizer.h"
#include "IO.h"
#include "RNG.h"
#include "RecombTable.h"

#include <boost/program_options.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace po = boost::program_options;

namespace {

std::string displayDefault(double value)
{
  std::ostringstream text;
  text << value;
  return text.str();
}

void printUsage(const char* executable,
                const po::options_description& options)
{
  std::cout << "Usage:\n  " << executable
            << " input.oscar output.oscar [options]\n\n"
            << options;
}

bool validConfig(const HadConfig& config)
{
  if (config.HQ_pid != 4 && config.HQ_pid != 5) {
    std::cerr << "HQ_pid must be 4 (charm) or 5 (bottom).\n";
    return false;
  }
  if (!(config.m_c > 0.0) || !(config.m_b > 0.0) ||
      !(config.m_q > 0.0) || !(config.m_s > 0.0) ||
      !(config.m_g > 0.0)) {
    std::cerr << "All constituent masses must be positive.\n";
    return false;
  }
  if (!(config.omega_cM > 0.0) || !(config.omega_cB > 0.0) ||
      !(config.omega_bM > 0.0) || !(config.omega_bB > 0.0)) {
    std::cerr << "All oscillator frequencies must be positive.\n";
    return false;
  }
  if (!std::isfinite(config.Tchem) || !(config.Tchem > 0.0)) {
    std::cerr << "Tchem must be finite and positive.\n";
    return false;
  }
  if (!std::isfinite(config.gamma_s) || !(config.gamma_s >= 0.0) ||
      !std::isfinite(config.gamma_HB) || !(config.gamma_HB > 0.0) ||
      !(config.eps_M > 0.0) ||
      !(config.eps_B > 0.0)) {
    std::cerr << "Invalid charm chemistry or fragmentation parameter.\n";
    return false;
  }
  return true;
}

} // namespace

int main(int argc, char** argv)
{
  HadConfig config;
  std::string input_path;
  std::string output_path;
  std::string recomb_path;
  std::string wigner_path;
  std::string io_format_name = "urqmd";
  int mode = static_cast<int>(config.mode);
  double omega_M = config.omega_cM;
  double omega_B = config.omega_cB;
  unsigned seed_option = 0;

  po::options_description options("Options");
  options.add_options()
      ("help,h", "Show this help")
      ("mode", po::value<int>(&mode)
                   ->value_name("<1|2|3>")
                   ->default_value(mode),
       "Hadronization mode: 1=Frag, 2=Recomb (current pass-through), "
       "3=Frag+Recomb")
      ("io-format", po::value<std::string>(&io_format_name)
                        ->value_name("<standard|urqmd|analysis>")
                        ->default_value(io_format_name),
       "Particle I/O format: standard, urqmd, or analysis")
      ("hq", po::value<int>(&config.HQ_pid)
                 ->default_value(config.HQ_pid),
       "Heavy-quark PDG: 4=charm, 5=bottom")
      ("recomb-table", po::value<std::string>(&recomb_path),
       "Recombination probability table "
       "(charm default: data/recomb_c_raw_M020_B0267.dat)")
      ("wigner-table", po::value<std::string>(&wigner_path),
       "Wigner-envelope table "
       "(charm default: data/max_wigner_c_M020_B0267.dat)")
      ("meson-table", po::value<std::string>(&config.charm_meson_table)
                            ->default_value(config.charm_meson_table),
       "Primary charm-meson state table")
      ("baryon-table", po::value<std::string>(&config.charm_baryon_table)
                             ->default_value(config.charm_baryon_table),
       "Primary charm-baryon state table")
      ("epsilon-meson",
       po::value<double>(&config.eps_M)
           ->default_value(
               config.eps_M,
               displayDefault(config.eps_M)),
       "Charm-meson Peterson epsilon")
      ("epsilon-baryon",
       po::value<double>(&config.eps_B)
           ->default_value(
               config.eps_B,
               displayDefault(config.eps_B)),
       "Charm-baryon Peterson epsilon")
      ("gamma-s", po::value<double>(&config.gamma_s)
                        ->default_value(config.gamma_s,
                                        displayDefault(config.gamma_s)),
       "Strange suppression factor")
      ("gamma-hb", po::value<double>(&config.gamma_HB)
                         ->default_value(config.gamma_HB,
                                         displayDefault(config.gamma_HB)),
       "Primary heavy-baryon fragmentation suppression factor")
      ("tchem", po::value<double>(&config.Tchem)
                    ->value_name("<value>")
                    ->default_value(config.Tchem,
                                    displayDefault(config.Tchem)),
       "Charm chemistry/hadronization temperature in GeV")
      ("omega-m", po::value<double>(&omega_M)
                        ->default_value(
                            omega_M,
                            displayDefault(omega_M)),
       "Active meson oscillator scale in GeV (shown default is charm)")
      ("omega-b", po::value<double>(&omega_B)
                        ->default_value(
                            omega_B,
                            displayDefault(omega_B)),
       "Active baryon oscillator scale in GeV (shown default is charm)")
      ("seed", po::value<unsigned>(&seed_option),
       "Deterministic RNG seed (default: generated)");

  po::options_description positional_options("Positional arguments");
  positional_options.add_options()
      ("input", po::value<std::string>(&input_path)->required(),
       "Input OSCAR file")
      ("output", po::value<std::string>(&output_path)->required(),
       "Output OSCAR file");

  po::options_description all_options;
  all_options.add(options).add(positional_options);
  po::positional_options_description positional;
  positional.add("input", 1).add("output", 1);

  po::variables_map variables;
  try {
    po::store(po::command_line_parser(argc, argv)
                  .options(all_options)
                  .positional(positional)
                  .run(),
              variables);
    if (variables.count("help") != 0U) {
      printUsage(argv[0], options);
      return 0;
    }
    po::notify(variables);
  } catch (const po::error& error) {
    std::cerr << "Command-line error: " << error.what() << "\n\n";
    printUsage(argv[0], options);
    return 1;
  }

  if (mode < static_cast<int>(HadMode::Frag) ||
      mode > static_cast<int>(HadMode::FragAndRecomb)) {
    std::cerr << "Hadronization mode must be 1, 2, or 3.\n";
    return 1;
  }
  config.mode = static_cast<HadMode>(mode);
  const bool uses_frag = config.mode == HadMode::Frag ||
                         config.mode == HadMode::FragAndRecomb;
  const bool uses_recomb = config.mode == HadMode::FragAndRecomb;

  IO::Format io_format = IO::Format::UrQMD;
  if (io_format_name == "standard") {
    io_format = IO::Format::Standard;
  } else if (io_format_name == "analysis") {
    io_format = IO::Format::Analysis;
  } else if (io_format_name != "urqmd") {
    std::cerr << "I/O format must be 'standard', 'urqmd', or 'analysis'.\n";
    return 1;
  }

  if (!variables["omega-m"].defaulted()) {
    if (config.HQ_pid == 5) {
      config.omega_bM = omega_M;
    } else {
      config.omega_cM = omega_M;
    }
  }
  if (!variables["omega-b"].defaulted()) {
    if (config.HQ_pid == 5) {
      config.omega_bB = omega_B;
    } else {
      config.omega_cB = omega_B;
    }
  }

  if (!validConfig(config)) return 1;

  const double default_omega_M = config.HQ_pid == 5
      ? HadConfig::default_omega_bM
      : HadConfig::default_omega_cM;
  const double default_omega_B = config.HQ_pid == 5
      ? HadConfig::default_omega_bB
      : HadConfig::default_omega_cB;
  const double active_omega_M = config.HQ_pid == 5
      ? config.omega_bM
      : config.omega_cM;
  const double active_omega_B = config.HQ_pid == 5
      ? config.omega_bB
      : config.omega_cB;
  const bool explicit_recomb_table =
      variables.count("recomb-table") != 0U;
  const bool explicit_wigner_table =
      variables.count("wigner-table") != 0U;

  if (uses_recomb) {
    const bool non_default_omega =
        active_omega_M != default_omega_M ||
        active_omega_B != default_omega_B;
    if (non_default_omega &&
        (!explicit_recomb_table || !explicit_wigner_table)) {
      std::cerr << "Non-default omega values require explicit "
                   "--recomb-table and --wigner-table.\n";
      return 1;
    }

    if (recomb_path.empty()) {
      recomb_path = config.HQ_pid == 5
          ? "recomb_b_raw.dat"
          : config.charm_recomb_table;
    }
    if (wigner_path.empty()) {
      wigner_path = config.HQ_pid == 5
          ? "max_wigner_b.dat"
          : config.charm_wigner_table;
    }
  }
  const unsigned seed = variables.count("seed") != 0U
      ? seed_option
      : static_cast<unsigned>(std::random_device{}());

  std::ifstream input(input_path);
  if (!input) {
    std::cerr << "Cannot open input file '" << input_path << "'.\n";
    return 2;
  }
  std::ofstream output(output_path);
  if (!output) {
    std::cerr << "Cannot open output file '" << output_path << "'.\n";
    return 3;
  }

  RecombTable recomb_table;
  if (uses_recomb &&
      !loadRecombTable(recomb_path, wigner_path, config, recomb_table)) {
    return 4;
  }

  RNG rng(seed);
  Hadronizer hadronizer(config, rng, recomb_table);
  if (!hadronizer.ready()) return 4;
  Event event;
  std::size_t n_events = 0;
  std::size_t n_input = 0;
  std::size_t n_output = 0;

  while (true) {
    std::string read_error;
    const IO::ReadStatus status =
        IO::readEventText(input, event, read_error, io_format);
    if (status == IO::ReadStatus::EndOfFile) break;
    if (status == IO::ReadStatus::Error) {
      std::cerr << "Failed to read '" << input_path << "' after "
                << n_events << " complete event(s): " << read_error << '\n';
      return 5;
    }

    ++n_events;
    n_input += event.particles.size();
    std::vector<Particle> hadrons;
    hadronizer.process(event.particles, hadrons);
    n_output += hadrons.size();
    event.particles = std::move(hadrons);
    IO::writeEventText(output, event, io_format);
  }

  const HadStats& stats = hadronizer.stats();
  const bool charm = config.HQ_pid == 4;
  const double m_HQ = charm ? config.m_c : config.m_b;
  const double omega_M_active = charm ? config.omega_cM : config.omega_bM;
  const double omega_B_active = charm ? config.omega_cB : config.omega_bB;
  std::cerr
      << "Hadronization summary\n"
      << "  input: " << input_path << '\n'
      << "  output: " << output_path << '\n'
      << "  RNG seed: " << seed << '\n'
      << "  heavy-quark PDG: " << config.HQ_pid << '\n'
      << "  mode: " << static_cast<int>(config.mode) << '\n'
      << "  io-format: " << io_format_name << '\n';
  if (uses_recomb) {
    std::cerr << "  model heavy-quark mass: " << m_HQ << " GeV\n";
  }
  if (uses_frag) {
    std::cerr
        << "  Tchem: " << config.Tchem << " GeV\n"
        << "  gamma_s: " << config.gamma_s << '\n'
        << "  gamma_HB: " << config.gamma_HB << '\n';
  }
  if (uses_recomb) {
    std::cerr
        << "  omega_M: " << omega_M_active << " GeV\n"
        << "  omega_B: " << omega_B_active << " GeV\n";
  }
  if (uses_frag) {
    std::cerr
        << "  epsilon_M: " << config.eps_M << '\n'
        << "  epsilon_B: " << config.eps_B << '\n';
  }
  if (uses_recomb) {
    std::cerr
        << "  recombination table: " << recomb_path << '\n'
        << "  Wigner table: " << wigner_path << '\n';
  }
  if (uses_frag) {
    std::cerr
        << "  meson table: " << config.charm_meson_table << '\n'
        << "  baryon table: " << config.charm_baryon_table << '\n';
  }
  std::cerr
      << "  events: " << n_events << '\n'
      << "  input particles: " << n_input << '\n'
      << "  output particles: " << n_output << '\n'
      << "  heavy quarks: " << stats.n_HQ << '\n'
      << "  recombined: " << stats.n_recomb << '\n'
      << "  fragmented: " << stats.n_frag << '\n'
      << "  dropped by unsupported channel: " << stats.n_dropped << '\n'
      << "  unchanged heavy quarks: " << stats.n_unchanged_HQ << '\n';

  if (!stats.final_species.empty()) {
    std::cerr << "  final heavy-flavor species:";
    for (const auto& [pid, count] : stats.final_species) {
      std::cerr << ' ' << pid << '=' << count;
    }
    std::cerr << '\n';
  }
  return 0;
}
