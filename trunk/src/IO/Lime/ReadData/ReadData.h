#ifndef GUARD_LIME_READDATA_H
#define GUARD_LIME_READDATA_H

#include <string>
#include <IO/Lime/c-lime/lime.h>

namespace Lime
{
  class ReadData
  {
    public:
      FILE *stream;
      LimeReader *reader;
    
      ReadData(std::string const &filename);
      ~ReadData();
  };
}

#endif