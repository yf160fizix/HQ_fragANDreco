#include "IO.h"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace {

class ff {
public:
  ff(double x) : value(x) {}

  const double value;

  friend std::ostream& operator<<(std::ostream& stream,
                                  const ff& field)
  {
    if (field.value == 0.0) {
      stream << "0.000000E+00";
      return stream;
    }

    int exponent = static_cast<int>(
        std::floor(std::log10(std::abs(field.value))));
    double base = field.value / std::pow(10.0, exponent);
    base /= 10.0;
    ++exponent;

    std::stringstream buffer;
    buffer << std::setw(8) << std::setprecision(6) << std::fixed << base;
    if (base >= 0.0) {
      stream << std::setw(8) << buffer.str();
    } else {
      const std::string negative =
          "-" + buffer.str().substr(2, buffer.str().size() - 1);
      stream << std::setw(8) << negative;
    }

    if (exponent >= 0) {
      stream << "E+" << std::setw(2) << std::setfill('0') << exponent;
    } else {
      stream << "E-" << std::setw(2) << std::setfill('0') << std::abs(exponent);
    }
    return stream;
  }
};

bool readNonBlankLine(std::istream& input, std::string& line)
{
  while (std::getline(input, line)) {
    if (line.find_first_not_of(" \t\r") != std::string::npos) return true;
  }
  return false;
}

} // namespace

IO::ReadStatus IO::readEventText(std::istream& input,
                                 Event& event,
                                 std::string& error_message,
                                 Format format)
{
  event.particles.clear();
  error_message.clear();

  std::string line;
  if (!readNonBlankLine(input, line)) {
    return input.eof() ? ReadStatus::EndOfFile : ReadStatus::Error;
  }
  if (!line.empty() && line.back() == '\r') line.pop_back();
  if (line != "OSC1997A") {
    error_message = "expected OSC1997A event header, found '" + line + "'";
    return ReadStatus::Error;
  }

  std::string field_header;
  std::string system_header;
  std::string event_header;
  if (!std::getline(input, field_header) ||
      !std::getline(input, system_header) ||
      !std::getline(input, event_header)) {
    error_message = "truncated OSCAR event header";
    return ReadStatus::Error;
  }

  int event_label = 0;
  int particle_count = 0;
  std::istringstream event_fields(event_header);
  if (!(event_fields >> event_label >> particle_count) || particle_count < 0) {
    error_message = "invalid OSCAR particle count line: '" + event_header + "'";
    return ReadStatus::Error;
  }

  event.particles.reserve(static_cast<std::size_t>(particle_count));
  for (int particle_index = 0; particle_index < particle_count; ++particle_index) {
    if (!std::getline(input, line)) {
      error_message = "unexpected end of file while reading particle " +
                      std::to_string(particle_index);
      event.particles.clear();
      return ReadStatus::Error;
    }

    Particle p;
    int file_index = 0;
    std::istringstream fields(line);
    bool parsed = static_cast<bool>(
        fields >> file_index >> p.pid
               >> p.px >> p.py >> p.pz
               >> p.E >> p.m
               >> p.x >> p.y >> p.z >> p.t
               >> p.Thydro
               >> p.cvx >> p.cvy >> p.cvz
               >> p.ipx >> p.ipy >> p.ipz);
    if (parsed && format == Format::UrQMD) {
      double trailing_weight = 0.0;
      parsed = static_cast<bool>(fields >> p.iE >> p.wt >> trailing_weight);
    } else if (parsed && format == Format::Analysis) {
      int origin = 0;
      parsed = static_cast<bool>(fields >> p.wt);
      if (parsed) {
        if (!(fields >> origin)) {
          error_message = "invalid origin in particle record " +
                          std::to_string(particle_index);
          event.particles.clear();
          return ReadStatus::Error;
        }
        if (origin == 0) {
          p.origin = HadOrigin::Unknown;
        } else if (origin == 1) {
          p.origin = HadOrigin::Frag;
        } else if (origin == 2) {
          p.origin = HadOrigin::Recomb;
        } else {
          error_message = "invalid origin in particle record " +
                          std::to_string(particle_index) + ": " +
                          std::to_string(origin);
          event.particles.clear();
          return ReadStatus::Error;
        }
      }
    } else if (parsed) {
      parsed = static_cast<bool>(fields >> p.wt);
    }
    if (!parsed) {
      error_message = "malformed particle record " +
                      std::to_string(particle_index) + ": '" + line + "'";
      event.particles.clear();
      return ReadStatus::Error;
    }

    std::string extra;
    if (fields >> extra) {
      error_message = "unexpected extra field in particle record " +
                      std::to_string(particle_index);
      event.particles.clear();
      return ReadStatus::Error;
    }
    event.particles.push_back(p);
  }

  ++event.event_id;
  return ReadStatus::Event;
}

void IO::writeEventText(std::ostream& output,
                        const Event& event,
                        Format format)
{
  output << "OSC1997A\n";
  output << "final_id_p_x\n";
  output << "    lbt  1.0alpha   208    82   208    82   aacm  0.1380E+04        1\n";
  output << "        1  " << std::setw(10) << std::setfill(' ')
         << event.particles.size()
         << "    0.001    0.001    1    1       1\n";

  std::size_t index = 0;
  for (const auto& p : event.particles) {
    output << std::setw(10) << std::setfill(' ') << index << "  "
           << std::setw(10) << p.pid << "  "
           << ff(p.px) << "  "
           << ff(p.py) << "  "
           << ff(p.pz) << "  "
           << ff(p.E) << "  "
           << ff(p.m) << "  "
           << ff(p.x) << "  "
           << ff(p.y) << "  "
           << ff(p.z) << "  "
           << ff(p.t) << "  "
           << ff(p.Thydro) << "  "
           << ff(p.cvx) << "  "
           << ff(p.cvy) << "  "
           << ff(p.cvz) << "  "
           << ff(p.ipx) << "  "
           << ff(p.ipy) << "  "
           << ff(p.ipz) << "  ";
    if (format == Format::UrQMD) {
      output << ff(p.iE) << "  "
             << ff(p.wt) << "  "
             << ff(0.0);
    } else if (format == Format::Analysis) {
      output << ff(p.wt) << "  "
             << static_cast<int>(p.origin);
    } else {
      output << ff(p.wt);
    }
    output << '\n';
    ++index;
  }
}
