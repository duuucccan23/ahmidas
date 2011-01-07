#pragma once

#include <string>
#include <cstring>

#include <lemon.h>
#include <mpi.h>

#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L0/Core/Component.h>
#include <L0/Core/Propagator.h>
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

    struct ScidacInfo
    {
      std::string field;
      std::string precision;
      std::string flavours;
      size_t      dims[4];

      ScidacInfo()
      {
	std::fill(dims, dims + 4, 0);
      }
      ScidacInfo(Lime::Reader &reader);
    };

    ScidacInfo parseScidacInfo(char *scidacCStr);

    //* LEMON VERSIONS *//
    void loadPropagator(Core::Field< QCD::Tensor > &field, std::string const &filename);
    void loadPropagator(Core::Field< QCD::Tensor > *field, size_t amount, std::string const &filename);
    void loadScidacInfo(LemonReader *reader, ScidacInfo &info);
    void savePropagator(Core::Propagator const &prop, std::vector< std::string > const &filename);
    void savePropagator(Core::StochasticPropagator< 4 > const &prop, std::vector< std::string > const &filename);
    void savePropagator(Core::StochasticPropagator< 1 > const &prop, std::vector< std::string > const &filename);
    void savePropagatorInfo(LemonWriter *writer, size_t L, size_t T);
    void savePropagatorType(LemonWriter *writer);

    template< typename Element >
    void loadScidacBinary(LemonReader *reader, Core::Field< Element > &field, ScidacInfo const &info);

    template< typename Element >
    void saveScidacBinary(LemonWriter *writer, Core::Field< Element > const &field);

    template< typename Element >
    void saveScidacBinary(LemonWriter *writer, Element *buffer, size_t volume, size_t const localVolume, int *dims);

    //* NON-LEMON VERSIONS *//
    void load(Core::Field< QCD::Gauge >  *field, std::string const &filename, Tool::IO::filetype);
    void load(Core::Field< QCD::Spinor > *field, std::string const &filename, Tool::IO::filetype);
    void load(Core::Field< QCD::Spinor > *field, std::string const &filename, Tool::IO::filetype, size_t const precision);

    void load(Core::Propagator *propagator, std::vector< std::string > const &filenames,
              Tool::IO::filetype type);

    void load(Core::Propagator *propagator, std::vector< std::string > const &filenames,
              Tool::IO::filetype type, size_t const precision);

    void load(Core::StochasticPropagator< 4 > *sPropagator, std::vector< std::string > const &filenames,
              Tool::IO::filetype type);

    void load(Core::StochasticPropagator< 4 > *sPropagator, std::vector< std::string > const &filenames,
              Tool::IO::filetype type, size_t const precision);

    void save(Core::StochasticPropagator< 4 > const *sPropagator, std::vector< std::string > const &filenames,
              Tool::IO::filetype type);

    void load(Core::StochasticPropagator< 1 > *sPropagator, std::vector< std::string > const &filenames,
              Tool::IO::filetype type);

    void load(Core::StochasticPropagator< 1 > *sPropagator, std::vector< std::string > const &filenames,
              Tool::IO::filetype type, size_t const precision);

    void save(Core::StochasticPropagator< 1 > const *sPropagator, std::vector< std::string > const &filenames,
              Tool::IO::filetype type);

    void save(Core::Field< QCD::Gauge >  const *field, std::string const &filename, Tool::IO::filetype);
    void save(Core::Field< QCD::Spinor > const *field, std::string const &filename, Tool::IO::filetype);

    void save(Core::Propagator const *propagator, std::vector< std::string > const &filenames,
              Tool::IO::filetype type);

    template< typename Element >
    void loadScidac(Core::Field< Element > *field, std::string const &filename);

    template< typename Element >
    void loadScidac(Core::Field< Element > *field, std::string const &filename, size_t const precision);

    template< typename Element >
    void loadILDG(Core::Field< Element > *field, std::string const &filename);


    template< typename Element >
    inline void loadMILC(Core::Field< Element > *field, std::string const &filename)
    {
      std::cerr << "loadMILC(...) has not been implemented yet for parallel architecture! Aborting..." << std::endl;
      exit(1);
    }

    template< typename Element >
    void saveILDG(Core::Field< Element > const &field, std::string const &filename);

    template< typename Element >
    void saveScidac(Core::Field< Element > const &field, std::string const &filename);

  }
}

#include "IO/loadScidac_a.template"
#include "IO/loadScidac_b.template"
#include "IO/loadILDG.template"
#include "IO/saveScidac.template"
#include "IO/saveILDG.template"

//* LEMON VARIETIES *//
#include "IO/loadScidacBinary.template"
#include "IO/saveScidacBinary.template"
