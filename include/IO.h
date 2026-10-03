#pragma once

#include <istream>
#include <ostream>
#include <string>

#include "Event.h"

namespace IO {

enum class Format {
  Standard,
  UrQMD,
  Analysis,
};

enum class ReadStatus { Event, EndOfFile, Error };

ReadStatus readEventText(std::istream& input,
                         Event& event,
                         std::string& error_message,
                         Format format);
void writeEventText(std::ostream& output,
                    const Event& event,
                    Format format);

} // namespace IO
