#ifndef GUARD_TOOL_IO_H
#define GUARD_TOOL_IO_H

#include <cassert>
#include <iomanip>
#include <sstream>
#include <string>

#include <L0/Base/IO/Lime/Writer.h>
#include <L0/Base/ScidacChecksum.h>
#include <L1/Source/Full.h>
#include <L1/Source/Stochastic.h>

namespace Tool
{
  namespace IO
  {
    template< size_t L, size_t T >
    void saveScidac(Source::Stochastic< L, T > const &source, std::string const &basefilename, size_t timeIdx);

    template< size_t L, size_t T >
    void saveScidac(Source::Full< L, T > &source, std::string const &basefilename);
  }
}

#include "IO/saveScidac.template"
#include "IO/saveScidac_b.template"

#endif
