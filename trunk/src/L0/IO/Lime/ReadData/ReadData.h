#ifndef GUARD_LIME_READDATA_H
#define GUARD_LIME_READDATA_H

#include <string>
#include <L0/IO/Lime/c-lime/lime.h>

namespace Lime
{
  struct ReadData
  {
    FILE              *stream;
    LimeReader        *reader;
    LimeRecordHeader  *header;

    ReadData(std::string const &filename);
    ~ReadData();
  };
}

#endif
