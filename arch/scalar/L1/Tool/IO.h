#ifndef GUARD_TOOL_IO_H
#define GUARD_TOOL_IO_H

#include <string>
#include <L1/Source/Stochastic.h>
#include <L0/Base/IO/Lime/Writer.h>
#include <sstream>

namespace Tool
{
  namespace IO
  {
    template< size_t L, size_t T >
    void saveScidac(Source::Stochastic< L, T > const &source, std::string const &basefilename, size_t timeslice);
  }
}

#include "IO/saveScidac.template"

#endif
