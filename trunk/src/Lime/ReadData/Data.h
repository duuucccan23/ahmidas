#ifndef GUARD_LIME_READDATA_H
#define GUARD_LIME_READDATA_H

#include <string>
#include <Lime/c-lime/lime.h>

namespace Lime
{
  class ReadData
  {
    public:
      FILE *stream;
      LimeReader *reader;
      char *headerType;
    
      ReadData(std::string const &filename);
      ~ReadData();
  };
}

#endif
