#ifndef GUARD_LIME_WRITEDATA_H
#define GUARD_LIME_WRITEDATA_H

#include <string>
#include <L0/IO/Lime/c-lime/lime.h>

namespace Lime
{
  class WriteData
  {
      FILE *stream;
      LimeWriter *writer;
      char *headerType;
    
      WriteData(std::string const &filename);
      ~WriteData();
  };
}

#endif
