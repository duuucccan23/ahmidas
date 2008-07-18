#ifndef GUARD_LIME_READER_H
#define GUARD_LIME_READER_H

#include <string>
#include <Core/Core.h>

namespace Lime
{
  struct ReadData;

  class Reader
  {
    ReadData    *d_data;
    int mutable  d_fail;
    uint64_t     d_size;

    public:
      Reader(std::string const &filename, std::string const &filetype);
      ~Reader();

      int fail() const;
      size_t size() const;

      void read(char *buffer, uint64_t elements) const;

      template< typename DataType >
      void read(DataType *buffer, uint64_t elements) const;
  };
}

#include "Reader.inlines"

#endif
