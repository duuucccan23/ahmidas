#ifndef GUARD_LIME_DATA_H
#define GUARD_LIME_DATA_H

#include <string>

namespace Lime
{
  class Data
  {
    public:
      FILE *stream;
      LimeReader *reader;
      char *headerType;
    
      Data(std::string const &filename);
      ~Data();
  };
}

#include "Data.inlines"

#endif
