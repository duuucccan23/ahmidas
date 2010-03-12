#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L0/Core/Propagator.h>
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
      fileSCIDAC,
      fileMILC
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
    void load(Core::Field< QCD::Spinor > *field, std::string const &filename, Tool::IO::filetype, size_t const precision);

    void load(Core::Propagator *propagator, std::vector< std::string > const &filenames,
              Tool::IO::filetype type);

    void load(Core::Propagator *propagator, std::vector< std::string > const &filenames,
              Tool::IO::filetype type, size_t const precision);

    void load(Core::StochasticPropagator< 4 > *sPropagator, std::vector< std::string > const &filenames,
              Tool::IO::filetype type);

    template< typename Element >
    void loadScidac(Core::Field< Element > *field, std::string const &filename);

    template< typename Element >
    void loadScidac(Core::Field< Element > *field, std::string const &filename, size_t const precision);

    template< typename Element >
    inline void loadILDG(Core::Field< Element > *field, std::string const &filename)
    {
      std::cerr << "loadILDG(...) has not been implemented yet for parallel architecture! Aborting..." << std::endl;
      exit(1);
    }

    template< typename Element >
    inline void loadMILC(Core::Field< Element > *field, std::string const &filename)
    {
      std::cerr << "loadMILC(...) has not been implemented yet for parallel architecture! Aborting..." << std::endl;
      exit(1);
    }

    template< typename Element >
    inline void saveILDG(Core::Field< Element > const &field, std::string const &filename)
    {
      std::cerr << "saveILDG(...) has not been implemented yet for parallel architecture! Aborting..." << std::endl;
      exit(1);
    }

    template< typename Element >
    inline void saveScidac(Core::Field< Element > const &field, std::string const &filename)
    {
      std::cerr << "saveScidac(...) has not been implemented yet for parallel architecture! Aborting..." << std::endl;
      exit(1);
    }

  }
}


