#ifndef GUARD_LIME_READER_H
#define GUARD_LIME_READER_H

#include <string>
#include <vector>

#include <L0/Core/Core.h>

namespace Lime
{
  // We want to avoid the inclusion of C-level libraries through the header
  struct FILE;
  struct LimeReader;
  struct LimeRecordHeader;

  class Reader
  {
    std::string const  d_filename;
    FILE              *d_stream;
    LimeReader        *d_reader;

    std::vector< int32_t >      d_messageIndices;
    std::vector< std::string >  d_limeTypes;
    std::vector< uint64_t >     d_recordSizes;
    std::vector< uint64_t >     d_recordOffsets;

    int32_t                     d_currentRecord;

    int                         d_fail;

    public:
      Reader(std::string const &filename);
      ~Reader();

      void findRecord(int32_t const message, int32_t const record);

      int fail() const;
      size_t size() const;

      template< typename DataType >
      void read(DataType *buffer, uint64_t elements) const;

      void read(char *buffer, uint64_t elements) const;
  };
}

#include "Reader/Reader.inlines"

#endif
