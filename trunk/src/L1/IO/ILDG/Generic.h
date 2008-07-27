#ifndef GUARD_IO_ILDG_GENERIC
#define GUARD_IO_ILDG_GENERIC

#include <string>

#include <L0/Core/Core.h>
#include <L0/IO/Lime/Reader.h>

namespace IO
{
  namespace ILDG
  {
    struct Format
    {
      std::string version;
      std::string field;
      size_t      precision;
      size_t      dimensions[4];
    };

    class Generic
    {
      std::string const  d_filename;
      Lime::Reader      *d_reader;

      Format             d_format;
      std::string        d_lfn;

      uint32_t           d_payloadMessage;
      uint32_t           d_formatRecord;
      uint32_t           d_dataRecord;
      uint32_t           d_lfnRecord;

      public:
        Generic(std::string const &filename);

        int readFile(std::string const &filename);
    };
  }
}

#include "Generic/Generic.inlines"

#endif
