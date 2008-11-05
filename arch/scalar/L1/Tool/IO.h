#ifndef GUARD_TOOL_IO_H
#define GUARD_TOOL_IO_H

#include <string>
#include <L1/Source/Stochastic.h>
#include <L0/Base/IO/Lime/Writer.h>
#include <sstream>
#include <iomanip>

namespace Tool
{
  namespace IO
  {
    template< size_t L, size_t T >
    void saveScidac(Source::Stochastic< L, T > const &source, std::string const &basefilename, size_t timeslice);

    template< size_t L, size_t T >
    void saveScidac(Source::Full< L, T > &source, std::string const &basefilename);
  }
}

#include "IO/saveScidac.template"
#include "IO/saveScidac_b.template"

#endif
