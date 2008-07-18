#ifndef GUARD_LIME_READDATA_H
#define GUARD_LIME_READDATA_H

#include <string>
#include <L0/IO/Lime/c-lime/lime.h>

namespace Lime
{
  class ReadData
  {
    public:
      FILE             *stream;
      LimeReader       *reader;
      LimeRecordHeader *header;
      std::string      limeType;

      ReadData(std::string const &filename);
      ~ReadData();
  };
}

#endif
