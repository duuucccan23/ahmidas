#pragma once

#include <fstream>
#include <string>
#include <L0/Base/ScidacChecksum.h>
#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L0/Tool/IO/Lime/Reader.h>
#include <L0/Tool/IO/Lime/Writer.h>

// Below is necessary to prevent circular including when we make IO functions friends of field.

namespace Tool
{
  namespace IO
  {
    struct ILDGinfo
    {
      std::string version;
      std::string field;
      std::string precision;
      size_t      dims[4];

      ILDGinfo(Lime::Reader &reader);
    };

    struct MILCinfo
    {
      int32_t     magicNumber;
      int32_t     dims[4];
      char        timestamp[64];
      int32_t     order;
      uint64_t    checksum;

      MILCinfo(std::ifstream &reader);
    };


    template< typename Element, size_t L, size_t T >
    Core::Field< Element, L, T > loadILDG(std::string const &filename);

    template< typename Element, size_t L, size_t T >
    Core::Field< Element, L, T > loadScidac(std::string const &filename);

    template< typename Element, size_t L, size_t T >
    Core::Field< Element, L, T > loadMILC(std::string const &filename);

    template< typename Element, size_t L, size_t T >
    void saveILDG(Core::Field< Element, L, T > const &field, std::string const &filename);

    template< typename Element, size_t L, size_t T >
    void saveScidac(Core::Field< Element, L, T > const &field, std::string const &filename);
  }
}
#include "IO/loadILDG.template"
#include "IO/loadScidac.template"
#include "IO/loadMILC.template"
#include "IO/saveILDG.template"
#include "IO/saveScidac.template"