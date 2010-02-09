#pragma once

#include <string>
#include <iostream>
#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L1/Tool.h>
#include <L1/Tool/IO/Lime/Reader.h>
#include <L1/Tool/ScidacChecksum.h>

namespace Tool
{
  namespace IO
  {
    enum filetype
    {
      fileILDG,
      fileSCIDAC
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
    void loadScidac(Core::Field< Element > *field, std::string const &filename);

    template< typename Element >
    inline void loadILDG(Core::Field< Element > *field, std::string const &filename)
    {
      std::cerr << "loadILDG(...) has not been implemented yet for parallel architecture! Aborting..." << std::endl;
      exit(1);
    }

  }
}

#include "IO/loadScidac.template"
