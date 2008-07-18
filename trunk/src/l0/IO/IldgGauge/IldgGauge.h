#ifndef GUARD_IO_ILDGGAUGE_H
#define GUARD_IO_ILDGGAUGE_H

#include <string>

namespace IO
{
  namespace Lime
  {
    class Reader;
    class Writer;
  }

  class IldgGauge
  {
    enum Mode
    {
      READ,
      WRITE
    };

    std::string  d_fileName;
    Mode         d_mode;
    LimeReader  *d_limeReader;
    LimeWriter  *d_limeWriter;

    // ILDG info
    std::string  d_version;
    std::string  d_field;
    size_t       d_precision;
    size_t       d_size[4];

    public:
      IldgGauge(std::string fileName, Mode mode);
  };
}

#endif
