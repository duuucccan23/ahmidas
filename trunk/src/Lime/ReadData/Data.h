#ifndef GUARD_LIME_DATA_H
#define GUARD_LIME_DATA_H

#include <string>
#include <Lime/c-lime/lime.h>

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

#endif
