#pragma once

#include <fstream>
#include <string>
#include <iostream>
#include <L0/Base/ScidacChecksum.h>
#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L0/Tool/IO/Lime/Reader.h>
#include <L0/Tool/IO/Lime/Writer.h>
#include <L1/Tool.h>

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

    struct Scidacinfo
    {
      std::string field;
      std::string precision;
      std::string flavours;
      size_t      dims[4];

      Scidacinfo(Lime::Reader &reader);
    };

    template< typename Element >
    Core::Field< Element > loadILDG(std::string const &filename, size_t L, size_t T);

    template< typename Element >
    Core::Field< Element > loadScidac(std::string const &filename, size_t L, size_t T);

    template< typename Element >
    Core::Field< Element > loadMILC(std::string const &filename, size_t L, size_t T);

    template< typename Element >
    void saveILDG(Core::Field< Element > const &field, std::string const &filename);

    template< typename Element >
    void saveScidac(Core::Field< Element > const &field, std::string const &filename);
  }
}
#include "IO/loadILDG.template"
#include "IO/loadScidac.template"
#include "IO/loadMILC.template"
#include "IO/saveILDG.template"
#include "IO/saveScidac.template"
