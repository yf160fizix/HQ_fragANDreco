#include "CharmStateTable.h"

#include <boost/math/special_functions/bessel.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace Frag {
namespace {

constexpr double kMD0     = 1.86484;
constexpr double kMDp     = 1.86966;
constexpr double kMDs     = 1.96834;
constexpr double kMDst0   = 2.00685;
constexpr double kMDstp   = 2.01026;
constexpr double kMLc     = 2.28646;
constexpr double kMXic    = 2.47000;  // effective Xi_c mass
constexpr double kMOmegac = 2.69520;

constexpr double kMPiC    = 0.13957;
constexpr double kMPi0    = 0.13498;
constexpr double kMN      = 0.93827;  // effective nucleon companion

double thermalWeight(const CharmHadronState& state,
                     double Tchem,
                     double gamma_s);

double modifiedBesselK2(double x)
{
#ifdef ISRM_HAS_STD_CYL_BESSEL_K
  return std::cyl_bessel_k(2.0, x);
#else
  return boost::math::cyl_bessel_k(2.0, x);
#endif
}

std::string trim(const std::string& s)
{
  const auto b = std::find_if_not(s.begin(), s.end(),
                                  [](unsigned char c){ return std::isspace(c); });
  const auto e = std::find_if_not(s.rbegin(), s.rend(),
                                  [](unsigned char c){ return std::isspace(c); }).base();
  if (b >= e) return "";
  return std::string(b, e);
}

std::string normalizeKey(std::string s)
{
  std::string out;
  out.reserve(s.size());

  for (char ch : s) {
    const unsigned char c = static_cast<unsigned char>(ch);

    if (ch == '*') {
      out += "star";
    } else if (ch == '+') {
      out += "plus";
    } else if (ch == '/') {
      out += "over";
    } else if (std::isalnum(c)) {
      out += static_cast<char>(std::tolower(c));
    }
  }

  return out;
}

std::vector<std::string> parseCSVLine(const std::string& line)
{
  std::vector<std::string> cells;
  std::string cell;
  bool inQuotes = false;

  for (size_t i = 0; i < line.size(); ++i) {
    const char ch = line[i];

    if (ch == '"') {
      if (inQuotes && i + 1 < line.size() && line[i + 1] == '"') {
        cell.push_back('"');
        ++i;
      } else {
        inQuotes = !inQuotes;
      }
    } else if (ch == ',' && !inQuotes) {
      cells.push_back(trim(cell));
      cell.clear();
    } else {
      cell.push_back(ch);
    }
  }

  cells.push_back(trim(cell));
  return cells;
}

struct Row {
  std::unordered_map<std::string, std::string> values;
  std::string path;
  std::size_t line_number = 0;
};

double parseNumber(const std::string& text,
                   const std::string& path,
                   std::size_t line_number,
                   const std::string& column)
{
  try {
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed == text.size() && std::isfinite(value)) return value;
  } catch (const std::invalid_argument&) {
  } catch (const std::out_of_range&) {
  }

  std::ostringstream error;
  error << path << ':' << line_number
        << ": invalid numeric value in column '" << column
        << "': '" << text << "'.";
  throw std::runtime_error(error.str());
}

std::vector<Row> readCSV(const std::string& path)
{
  std::vector<Row> rows;

  std::ifstream fin(path);
  if (!fin) {
    throw std::runtime_error(
        "CharmStateTable: cannot open CSV file: " + path);
  }

  std::string headerLine;
  if (!std::getline(fin, headerLine)) return rows;

  auto headersRaw = parseCSVLine(headerLine);
  std::vector<std::string> headers;
  headers.reserve(headersRaw.size());

  for (size_t i = 0; i < headersRaw.size(); ++i) {
    std::string h = normalizeKey(headersRaw[i]);
    if (h.empty() && i == 0) h = "family";
    headers.push_back(h);
  }

  std::string line;
  std::size_t line_number = 1;
  while (std::getline(fin, line)) {
    ++line_number;
    if (trim(line).empty()) continue;

    auto cells = parseCSVLine(line);
    if (cells.size() != headers.size()) {
      std::ostringstream error;
      error << path << ':' << line_number << ": expected "
            << headers.size() << " fields, found " << cells.size() << '.';
      throw std::runtime_error(error.str());
    }

    Row row;
    row.path = path;
    row.line_number = line_number;

    for (size_t i = 0; i < cells.size(); ++i) {
      const std::string text = trim(cells[i]);
      if (!text.empty() && headers[i] != "family" &&
          headers[i] != "descriptor") {
        parseNumber(text, path, line_number, headersRaw[i]);
      }
      row.values[headers[i]] = text;
    }

    rows.push_back(std::move(row));
  }

  return rows;
}

bool hasText(const Row& row, const std::vector<std::string>& keys)
{
  for (const auto& key : keys) {
    const auto it = row.values.find(normalizeKey(key));
    if (it != row.values.end() && !trim(it->second).empty()) return true;
  }
  return false;
}

std::string getString(const Row& row,
                      const std::vector<std::string>& keys,
                      const std::string& def = "")
{
  for (const auto& key : keys) {
    const auto it = row.values.find(normalizeKey(key));
    if (it != row.values.end() && !trim(it->second).empty()) {
      return trim(it->second);
    }
  }
  return def;
}

double getDouble(const Row& row,
                 const std::vector<std::string>& keys,
                 double def = 0.0)
{
  for (const auto& key : keys) {
    const auto it = row.values.find(normalizeKey(key));
    if (it == row.values.end()) continue;

    const std::string text = trim(it->second);
    if (text.empty()) continue;

    return parseNumber(text, row.path, row.line_number, key);
  }
  return def;
}

int getInt(const Row& row,
           const std::vector<std::string>& keys,
           int def = 0)
{
  return static_cast<int>(std::lround(getDouble(row, keys, def)));
}

bool containsNoCase(std::string text, std::string pattern)
{
  std::transform(text.begin(), text.end(), text.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  std::transform(pattern.begin(), pattern.end(), pattern.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  return text.find(pattern) != std::string::npos;
}

double companionOrZero(double m_parent,
                       double m_daughter,
                       double m_companion)
{
  // If the "decay" is really identity/direct production, do not soften.
  if (std::abs(m_parent - m_daughter) < 1e-4) return 0.0;

  // If the nominal companion is kinematically impossible,
  // use a zero-mass effective companion instead of rejecting the channel.
  if (m_parent <= m_daughter + m_companion) return 0.0;

  return m_companion;
}

void addMode(std::vector<DecayChannel>& modes,
             int daughter_pid,
             double m_daughter,
             double m_companion,
             double BR)
{
  if (!(BR > 0.0)) return;
  modes.push_back({daughter_pid, m_daughter, m_companion, BR});
}

int fakeMesonPDG(size_t i)
{
  return 800000 + static_cast<int>(i);
}

int fakeBaryonPDG(size_t i)
{
  return 900000 + static_cast<int>(i);
}

int identifyMesonPDG(double mass,
                     double BR_D0,
                     double BR_Dp,
                     double BR_Dst0,
                     double BR_Dstp,
                     double BR_Ds)
{
  if (BR_D0 > 0.99 && std::abs(mass - kMD0) < 0.02) return 421;
  if (BR_Dp > 0.99 && std::abs(mass - kMDp) < 0.02) return 411;
  if (BR_Dst0 > 0.99 && std::abs(mass - kMDst0) < 0.03) return 423;
  if (BR_Dstp > 0.99 && std::abs(mass - kMDstp) < 0.03) return 413;
  if (BR_Ds > 0.99 && std::abs(mass - kMDs) < 0.03) return 431;
  return 0;
}

std::optional<CharmHadronState> makeMesonState(const Row& row,
                                               std::size_t row_index,
                                               double Tchem,
                                               double gamma_s)
{
  if (!hasText(row, {"mass(GeV)", "mass"})) return std::nullopt;

  const double mass = getDouble(row, {"mass(GeV)", "mass"}, 0.0);
  if (!(mass > 0.0)) return std::nullopt;

  const double BR_D0 = getDouble(row, {"directBR_to_D0"}, 0.0);
  const double BR_Dp =
      getDouble(row, {"directBR_to_D+", "directBR_to_Dplus"}, 0.0);
  const double BR_Dst0 =
      getDouble(row, {"directBR_to_D*0", "directBR_to_Dstar0"}, 0.0);
  const double BR_Dstp =
      getDouble(row, {"directBR_to_D*+", "directBR_to_Dstar+"}, 0.0);
  const double BR_Ds =
      getDouble(row, {"directBR_to_Ds+", "directBR_to_Ds"}, 0.0);

  const double m_D_companion = getDouble(
      row, {"decay_D0", "decay_D", "decayD0", "decayD"}, kMPiC);
  const double m_Ds_companion = getDouble(
      row, {"decay_Ds", "decay_Ds+", "decayDs", "decayDsplus"}, kMPi0);
  const double m_Dst_companion = getDouble(
      row, {"decay_Dstar", "decay_D*", "decayDstar"}, kMPiC);

  CharmHadronState state;
  state.name = "meson_" + std::to_string(row_index);
  state.mass = mass;
  state.strangeness = getInt(row, {"strangeness"}, 0);
  state.spin = getDouble(row, {"spin"}, 0.0);
  state.isospin = getDouble(row, {"isospin"}, 0.0);

  const int known_pid = identifyMesonPDG(
      mass, BR_D0, BR_Dp, BR_Dst0, BR_Dstp, BR_Ds);
  state.pid = known_pid != 0 ? known_pid : fakeMesonPDG(row_index);

  addMode(state.decay_channels, 421, kMD0,
          companionOrZero(mass, kMD0, m_D_companion), BR_D0);
  addMode(state.decay_channels, 411, kMDp,
          companionOrZero(mass, kMDp, m_D_companion), BR_Dp);
  addMode(state.decay_channels, 423, kMDst0,
          companionOrZero(mass, kMDst0, m_Dst_companion), BR_Dst0);
  addMode(state.decay_channels, 413, kMDstp,
          companionOrZero(mass, kMDstp, m_Dst_companion), BR_Dstp);
  addMode(state.decay_channels, 431, kMDs,
          companionOrZero(mass, kMDs, m_Ds_companion), BR_Ds);

  if (state.decay_channels.empty()) {
    std::cerr << "CharmStateTable: meson row " << row_index
              << " has no decay modes; skipped.\n";
    return std::nullopt;
  }

  state.statistical_weight =
      thermalWeight(state, Tchem, gamma_s);
  if (!(state.statistical_weight > 0.0)) return std::nullopt;
  return state;
}

std::optional<CharmHadronState> makeBaryonState(
    const Row& row,
    std::size_t row_index,
    double Tchem,
    double gamma_s)
{
  if (!hasText(row, {"mass(GeV)", "mass"})) return std::nullopt;

  const double mass = getDouble(row, {"mass(GeV)", "mass"}, 0.0);
  if (!(mass > 0.0)) return std::nullopt;

  const std::string family = getString(row, {"family"}, "baryon");
  const std::string descriptor =
      getString(row, {"Descriptor", "descriptor"}, "");
  const double BR_Lambda =
      getDouble(row, {"BR_to_Lambda", "BR_to_Lambda_c"}, 0.0);
  const double BR_D =
      getDouble(row, {"BR_to_D+/D0", "BR_to_Dplus/D0"}, 0.0);
  const double BR_Ds =
      getDouble(row, {"BR_to_Ds+", "BR_to_Ds"}, 0.0);
  const double BR_Xi =
      getDouble(row, {"BR_to_Cascade", "BR_to_Xi_c"}, 0.0);
  const double BR_Omega =
      getDouble(row, {"BR_to_Omega_c", "BR_to_Omega"}, 0.0);

  const double m_Lambda_companion = getDouble(
      row, {"Decay_Lambda", "decay_Lambda", "DecayLambda"}, kMPiC);
  const double m_D_companion =
      getDouble(row, {"Decay_D", "decay_D", "DecayD"}, kMN);
  const double m_Ds_companion = getDouble(
      row, {"Decay_Ds", "decay_Ds", "DecayDs"}, m_D_companion);
  const double m_Xi_companion = getDouble(
      row, {"Decay_Cascade", "decay_Cascade", "DecayXi", "Decay_Xi"},
      kMPiC);
  const double m_Omega_companion = getDouble(
      row, {"Decay_Omega", "decay_Omega", "DecayOmega"}, kMPiC);

  CharmHadronState state;
  state.name = family + "_" + descriptor + "_" + std::to_string(row_index);
  state.mass = mass;
  state.strangeness = getInt(row, {"strangeness"}, 0);
  state.spin = getDouble(row, {"spin"}, 0.0);
  state.isospin = getDouble(row, {"isospin"}, 0.0);
  state.baryon = true;
  
  if (containsNoCase(family, "lambda")) {
    state.pid = std::abs(mass - kMLc) < 0.02
        ? 4122
        : fakeBaryonPDG(row_index);
  } else if (containsNoCase(family, "cascade")) {
    state.pid = std::abs(mass - kMXic) < 0.08
        ? 4232
        : fakeBaryonPDG(row_index);
  } else if (containsNoCase(family, "omega")) {
    state.pid = std::abs(mass - kMOmegac) < 0.08
        ? 4332
        : fakeBaryonPDG(row_index);
  } else {
    state.pid = fakeBaryonPDG(row_index);
  }

  addMode(state.decay_channels, 4122, kMLc,
          companionOrZero(mass, kMLc, m_Lambda_companion), BR_Lambda);
  // The aggregate D+/D0 branch is split equally between charge states.
  addMode(state.decay_channels, 421, kMD0,
          companionOrZero(mass, kMD0, m_D_companion), 0.5 * BR_D);
  addMode(state.decay_channels, 411, kMDp,
          companionOrZero(mass, kMDp, m_D_companion), 0.5 * BR_D);
  addMode(state.decay_channels, 431, kMDs,
          companionOrZero(mass, kMDs, m_Ds_companion), BR_Ds);
  addMode(state.decay_channels, 4232, kMXic,
          companionOrZero(mass, kMXic, m_Xi_companion), BR_Xi);
  addMode(state.decay_channels, 4332, kMOmegac,
          companionOrZero(mass, kMOmegac, m_Omega_companion), BR_Omega);

  if (state.decay_channels.empty()) {
    std::cerr << "CharmStateTable: baryon row " << row_index
              << " has no decay modes; skipped.\n";
    return std::nullopt;
  }

  state.statistical_weight =
      thermalWeight(state, Tchem, gamma_s);
  if (!(state.statistical_weight > 0.0)) return std::nullopt;
  return state;
}

std::vector<CharmHadronState> loadMesons(const std::string& path,
                                         double Tchem,
                                         double gamma_s)
{
  const auto rows = readCSV(path);
  std::vector<CharmHadronState> states;
  states.reserve(rows.size());
  for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
    auto state = makeMesonState(
        rows[row_index], row_index, Tchem, gamma_s);
    if (state) states.push_back(std::move(*state));
  }
  return states;
}

std::vector<CharmHadronState> loadBaryons(const std::string& path,
                                          double Tchem,
                                          double gamma_s)
{
  const auto rows = readCSV(path);
  std::vector<CharmHadronState> states;
  states.reserve(rows.size());
  for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
    auto state = makeBaryonState(
        rows[row_index], row_index, Tchem, gamma_s);
    if (state) states.push_back(std::move(*state));
  }
  return states;
}

double thermalWeight(const CharmHadronState& s, double Tchem, double gamma_s)
{
  if (!(s.mass > 0.0) || !(Tchem > 0.0)) return 0.0;

  const double isoFactor = 2.0 * s.isospin + 1.0;

  // RQM-average rows intentionally use the same expression with their
  // tabulated effective spin values.
  const double spinFactor = 2.0 * s.spin + 1.0;

  const double degeneracy = spinFactor * isoFactor;

  if (!(degeneracy > 0.0)) return 0.0;

  const double strangeFactor =
      std::pow(gamma_s, static_cast<double>(std::max(0, s.strangeness)));

  const double x = s.mass / Tchem;

  double k2 = 0.0;
  try {
    k2 = modifiedBesselK2(x);
  } catch (...) {
    k2 = 0.0;
  }

  if (!(k2 > 0.0) || !std::isfinite(k2)) return 0.0;

  return degeneracy * strangeFactor * s.mass * s.mass * Tchem * k2;
}

} // namespace

std::vector<CharmHadronState>
loadCharmStates(const std::string& meson_path,
                const std::string& baryon_path,
                double Tchem,
                double gamma_s,
                double gamma_HB)
{
  std::vector<CharmHadronState> out;

  auto mesons = loadMesons(meson_path, Tchem, gamma_s);
  auto baryons = loadBaryons(baryon_path, Tchem, gamma_s);
  if (mesons.empty()) {
    throw std::runtime_error(
        "CharmStateTable: no valid meson states loaded from " + meson_path);
  }
  if (baryons.empty()) {
    throw std::runtime_error(
        "CharmStateTable: no valid baryon states loaded from " + baryon_path);
  }
  for (CharmHadronState& baryon : baryons) {
    baryon.statistical_weight *= gamma_HB;
  }

  out.insert(out.end(), mesons.begin(), mesons.end());
  out.insert(out.end(), baryons.begin(), baryons.end());

  return out;
}

} // namespace Frag
