#pragma once

#include <string>
#include <iostream>
#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L1/Tool.h>
#include <L1/Tool/IO/Lime/Reader.h>
#include <L1/Tool/IO/Lime/Writer.h>
#include <L1/Tool/ScidacChecksum.h>

namespace Tool
{
  namespace IO
  {
    enum filetype
    {
      fileILDG,
      fileSCIDAC,
      fileMILC
    };

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

    void load(Core::Field< QCD::Gauge > *field, std::string const &filename, Tool::IO::filetype);
    void load(Core::Field< QCD::Spinor > *field, std::string const &filename, Tool::IO::filetype);

    template< typename Element >
    void loadILDG(Core::Field< Element > *field, std::string const &filename);

    template< typename Element >
    void loadScidac(Core::Field< Element > *field, std::string const &filename);

    template< typename Element >
    void loadMILC(Core::Field< Element > *field, std::string const &filename);

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